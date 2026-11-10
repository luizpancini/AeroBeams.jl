using AeroBeams, Random

# Root pitch angle [rad]
θ = deg2rad(5)

# Nominal (constant) airspeed [m/s]
Unom = 47

# Time variables [s]
tf = 60
Δt = 2e-4
trackingFrequency = clamp(round(Int,1e-3/Δt), 1, 5)

# Airspeed profile: constant at Unom, with a small random noise to induce flutter
σn = 1e-3 * Unom
rng = Xoshiro(1234)
smoothingWindow = Int(1e4)
noiseRaw = σn .* randn(rng, round(Int, tf/Δt) + 1)
noise = moving_average(noiseRaw, smoothingWindow)
U = t -> Unom + noise[clamp(round(Int, t/Δt) + 1, 1, length(noise))]

# Tip correction function type
tipLossType = "Exponential"

# Mesh (defines number of elements): :nominal or :refined
mesh = :nominal

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Geometric properties
nElem,L,chord,normSparPos,nodalPositions = geometrical_properties_Pazy(mesh=mesh)

# Airspeed range for flutter problem [m/s]
URange = collect(0.5:0.5:90)

# Number of modes for flutter computation
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed in flutter problem
eigenmethod = aeroSolver isa DynamicStallSolver ? :fulldecomp : :krylovkit
nev = 70

# Set system solver options for dynamic problem
systemSolver = aeroSolver isa DynamicStallSolver ?
               create_NewtonRaphson(absTol=1e-7,displayStatus=false) :
               create_NLsolver(method=:newton,linesearch=false,displayStatus=false)

# Strain gauge coordinates (the spar cs is 60 x 2.25 mm)
SGelem = 3          # element number where the strain gauge is located (see Fig. 12 of https://doi.org/10.2514/1.J060621: FBG^1_3 is @ 127 mm spanwise, but wing root is at 32 mm -> actual spanwise position is 127 - 32 = 95 mm which is at the midpoint of element 3)
x2SG = 50e-3/2      # chordwise position: close to LE
x3SG = 2.25e-3/2    # height: on top

# Case ID
idString = string("Theta", round(rad2deg(θ),digits=1), "_U", Unom, "_tf", tf, "_", tipLossType, "_", aeroSolver.name)

## --- Solve flutter problem ---
println("Solving flutter problem")

# Initialize outputs
flutterProblem = Array{EigenProblem}(undef,length(URange))
eigensols = Array{EigenSolution}(undef,length(URange))
tipOOP_f = Array{Float64}(undef,length(URange))
tipAoA_f = Array{Float64}(undef,length(URange))
epsSG_f = Array{Float64}(undef,length(URange))

# Sweep airspeed
for (j,U) in enumerate(URange)
    # Model
    flutterModel = first(create_Pazy(aeroSolver=aeroSolver,mesh=mesh,airfoil=airfoil,upright=true,hasTipCorrection=true,tipLossType=tipLossType,pitchAngle=θ,airspeed=U))
    # Set initial guess solution as previously available one
    x0 = j == 1 ? zeros(0) : flutterProblem[j-1].systemArrays.x
    # Create and solve problem
    flutterProblem[j] = create_EigenProblem(model=flutterModel,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(flutterProblem[j])
    # Eigensolution
    eigensols[j] = flutterProblem[j].eigensolution
    # Kinematics
    states = flutterProblem[j].loadstepHistory.states
    tipOOP_f[j] = -states[end][nElem].nodal.u2[1]
    tipAoA_f[j] = flutterProblem[j].loadstepHistory.aeroVariables[end][nElem].unsteadyKin.αₑ
    ϵ11SG_f = states[end][SGelem].midpoint.γ[1]
    κ2SG_f = states[end][SGelem].midpoint.κ[2]
    κ3SG_f = states[end][SGelem].midpoint.κ[3]
    epsSG_f[j] = (ϵ11SG_f .- κ2SG_f*x3SG .- κ3SG_f*x2SG)
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(3.75), fmax=Hz2rads(50), maxTracks=nModes, dropOutsideBounds=true)

# Flutter onset and offset speeds, frequencies and kinematics
flutterOnsetSpeeds,flutterOnsetFreqs,_,flutterOnsetTipOOPs,flutterOnsetTipAoAs,flutterOnsetQuarterEpss,flutterOffsetSpeeds,flutterOffsetFreqs = find_flutter(URange,freqsOsc,dampsOsc,tipOOP_f,tipAoA_f,epsSG_f)
if !isempty(flutterOnsetSpeeds)
    flutterOnsetSpeed = flutterOnsetSpeeds[1]
    flutterOnsetFreq = flutterOnsetFreqs[1]
    flutterOnsetTipOOP = flutterOnsetTipOOPs[1]
    flutterOnsetTipAoA = flutterOnsetTipAoAs[1]
    flutterOnsetQuarterEps = flutterOnsetQuarterEpss[1]
    println("Flutter onset speed = $(round(flutterOnsetSpeed,digits=1)) m/s, frequency = $(round(rads2Hz(flutterOnsetFreq),digits=1)) Hz, tip OOP = $(round(flutterOnsetTipOOP/L*100,digits=1))% semispan, tip AoA = $(round(rad2deg(flutterOnsetTipAoA),digits=1)) deg, strains = $(round(Int,flutterOnsetQuarterEps*1e6)) microns")
else
    println("Flutter not found")
end
if !isempty(flutterOffsetSpeeds)
    flutterOffsetSpeed = flutterOffsetSpeeds[1]
    flutterOffsetFreq = flutterOffsetFreqs[1]
    println("Flutter offset speed = $(round(flutterOffsetSpeed,digits=1)) m/s, frequency = $(round(rads2Hz(flutterOffsetFreq),digits=1)) Hz")
end

# Check if nominal airspeed lies within the expected flutter range (between onset and offset speeds)
if isnan(flutterOnsetSpeed)
    println("Nominal airspeed = $(Unom) m/s is outside the expected flutter range (no flutter onset found)")
else
    inFlutterRange = Unom >= flutterOnsetSpeed && (isnan(flutterOffsetSpeed) || Unom <= flutterOffsetSpeed)
    println("Nominal airspeed = $(Unom) m/s is $(inFlutterRange ? "WITHIN" : "OUTSIDE") the expected flutter range (onset = $(round(flutterOnsetSpeed,digits=1)) m/s, offset = $(round(flutterOffsetSpeed,digits=1)) m/s)")
end

# Plots
relPath = "/dev/Pazy/S0/InertiaRod/outputs/PazyWingConstantAirspeed"
absPath = string(pwd(),relPath)
mkpath(absPath)
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 14
lw = 2
gr()
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=extrema(URange), ylims=[0,50], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
for mode in 1:nModes
    plot!(URange, rads2Hz.(freqsOsc[mode,:]), c=modeColors[mode], lw=lw, label=false)
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URange), ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
plot!(URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
for mode in 1:nModes
    plot!(URange, dampsOsc[mode,:]./freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(plt_Vgf,string(absPath,"/PazyWingConstantAirspeed_",idString,"_Vgf.pdf"))

## --- Solve dynamic problem ---
println("Solving dynamic problem")

# Model
dynModel,_ = create_Pazy(aeroSolver=aeroSolver,mesh=mesh,airfoil=airfoil,upright=true,pitchAngle=θ,airspeed=U,hasTipCorrection=true,tipLossType=tipLossType,tipLossFunctionIsAirspeedDependent=true)

# Create and solve steady problem for initial airspeed
steadyProblem = create_SteadyProblem(model=dynModel)
solve!(steadyProblem)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=dynModel,finalTime=tf,Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x,trackingFrequency=trackingFrequency)
solve!(dynamicProblem)

# Unpack numerical solution
sol = dynamicProblem.solutionHistory
t = sol.t
states = sol.states
aeroVars = sol.aeroVariables
tipAoA = [aeroVars[i][nElem].unsteadyKin.αₑ for i in eachindex(t)]
tipOOP = -[states[i][nElem].nodal.u2[1] for i in eachindex(t)]
ϵ11SG = [states[i][SGelem].midpoint.γ[1] for i in eachindex(t)]
κ2SG = [states[i][SGelem].midpoint.κ[2] for i in eachindex(t)]
κ3SG = [states[i][SGelem].midpoint.κ[3] for i in eachindex(t)]
eps_SG = ϵ11SG .- κ2SG*x3SG .- κ3SG*x2SG

## Plots
# Set paths
relPath = "/dev/Pazy/S0/InertiaRod/outputs/PazyWingConstantAirspeed"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
# anim = plot_dynamic_deformation(dynamicProblem,follow=:basisA,plotFrequency=round(Int,1e-2/Δt),fps=24,plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/PazyWingConstantAirspeed_",idString,"_deformation.gif"),displayProgress=true)
# display(anim)

# Plot configurations
interactivePlots = false
colors = cgrad(:rainbow, 2, categorical=true)
epsColors = cgrad(:rainbow, 5, categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
fftWindow = 2^11
fftOverlap = 0.5
psdWindow = 2^13
if interactivePlots
    plotlyjs()
else
    gr()
end

# Time interval for the root strains PSD [s] (must be within [0,tf] and span at least psdWindow samples)
tPSDStart = 15
tPSDEnd = isnothing(findfirst(isnan, tipAoA)) ? tf : t[findfirst(isnan, tipAoA)] - 1

# Airspeed
plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
plot!(t, U.(t), lw=lw, c=:black, label=false)
display(plt_U)
savefig(string(absPath,"/PazyWingConstantAirspeed_",idString,"_U.pdf"))

# Tip OOP displacement
plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
plot!(t, tipOOP./L*100, lw=lw, c=colors[1], label=false)
display(plt_tipOOP)
savefig(string(absPath,"/PazyWingConstantAirspeed_",idString,"_tipOOP.pdf"))

# Tip AoA
plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
plot!(t, rad2deg.(tipAoA), lw=lw, c=colors[1], label=false)
display(plt_tipAOA)
savefig(string(absPath,"/PazyWingConstantAirspeed_",idString,"_tipAoA.pdf"))

# Strains at the selected spanwise location (approximately 1/6-span)
plt_eps = plot(xlabel="Time [s]", ylabel="1/6-span strains (\$\\mu\\epsilon\$)", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(t, eps_SG*1e6, lw=lw, c=colors[1], label=false)
display(plt_eps)
savefig(string(absPath,"/PazyWingConstantAirspeed_",idString,"_eps.pdf"))

# Strains spectrogram over time
f_eps, t_eps, S_eps = spectrogram(t, eps_SG*1e6; window=fftWindow, overlap=fftOverlap)
plt_eps_spectrum = heatmap(t_eps, f_eps, S_eps; ylims=[0,100], title="Strains spectrogram", xlabel="Time [s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
if !isnan(flutterOnsetFreq)
    for k in 1:3
        plot!([0,tf], k*rads2Hz(flutterOnsetFreq)*[1,1], ls=:dash, c=:black, lw=lw, label=false)
    end
end
display(plt_eps_spectrum)
savefig(string(absPath,"/PazyWingConstantAirspeed_",idString,"_eps_SGSpectrogram.pdf"))

# Strains PSD (Welch estimate) for the chosen time interval
@assert 0 <= tPSDStart < tPSDEnd <= tf "Invalid PSD time interval"
iPSD = findall(ti -> tPSDStart <= ti <= tPSDEnd, t)
@assert length(iPSD) >= psdWindow "PSD time interval is too short for the chosen psdWindow"
epsPSD = eps_SG[iPSD]*1e6 .- sum(eps_SG[iPSD]*1e6)/length(iPSD)
f_psd, _, S_psd = spectrogram(t[iPSD], epsPSD; window=psdWindow, overlap=fftOverlap, db=false)
psd_eps = 10 .* log10.(vec(sum(S_psd, dims=2)./size(S_psd,2)) .+ eps())
fMask = 0 .< f_psd .<= 100
iDominant = argmax(psd_eps[fMask])
fDominant = f_psd[fMask][iDominant]
println("Dominant root strain frequency in [$(tPSDStart), $(tPSDEnd)] s: $(round(fDominant,digits=1)) Hz")
plt_eps_PSD = plot(xlabel="Frequency [Hz]", ylabel="Strains PSD [dB]", xlims=[0,100], title="Strains PSD ($(tPSDStart)-$(tPSDEnd) s)", tickfont=font(ts), guidefont=font(fs))
plot!(f_psd[fMask], psd_eps[fMask], lw=lw, c=colors[2], label=false)
scatter!([fDominant], [psd_eps[fMask][iDominant]], marker=:circle, ms=6, msw=1, mc=:white, msc=:red, label="Dominant")
if !isnan(flutterOnsetFreq)
    for k in 1:3
        plot!(k*rads2Hz(flutterOnsetFreq)*[1,1], [minimum(psd_eps[fMask]), maximum(psd_eps[fMask])], ls=:dash, c=:black, lw=lw, label=false)
    end
end
display(plt_eps_PSD)
savefig(string(absPath,"/PazyWingConstantAirspeed_",idString,"_eps_SGPSD.pdf"))

println("Finished PazyWingConstantAirspeed.jl")