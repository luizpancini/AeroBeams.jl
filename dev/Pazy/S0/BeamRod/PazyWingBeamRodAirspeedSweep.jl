using AeroBeams, Random

# Test set
test = 2

# Root pitch angle offset (assumed) [rad]
θoffset = deg2rad(0)

# Root pitch angle [rad]
θ = test == 1 ? deg2rad(7) :
    test == 2 ? deg2rad(5) :
    test == 3 ? deg2rad(3) :
    error("Invalid test number")
θ += θoffset

# Time variables
tf = test == 1 ? 1.5*30 :
     test == 2 ? 25 :
     test == 3 ? 30 :
     error("Invalid test number")
Δt = 1e-4
trackingFrequency = clamp(round(Int,1e-3/Δt), 1, 5)

# Profile for airspeed sweep [m/s] (contains a small random noise to induce flutter)
Ui = test == 1 ? 35 :
     test == 2 ? 40 :
     test == 3 ? 40 :
     error("Invalid test number")
Uf = test == 1 ? 50 :
     test == 2 ? 52 :
     test == 3 ? 60 :
     error("Invalid test number")
σn = 1e-4 * Uf
rng = Xoshiro(1234)
smoothingWindow = Int(1e2)
noiseRaw = σn .* randn(rng, round(Int, tf/Δt) + 1)
noise = moving_average(noiseRaw, smoothingWindow)
U = t -> Ui + (Uf-Ui)*t/tf + noise[clamp(round(Int, t/Δt) + 1, 1, length(noise))]

# Time integration method for the dynamic solver
integrationMethod = GeneralizedAlpha(ρ∞=0.8)

# Mesh (defines number of elements): :nominal or :refined
mesh = :nominal

# Number of elements in the wingtip rod
nElementsWingtipRod = 12

# Tip correction function type
tipLossType = "Exponential"

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Geometric properties
nElem,L,chord,normSparPos,_ = geometrical_properties_Pazy()

# Airspeed range for flutter problem [m/s]
URange = collect(0.5:0.5:90)

# Number of modes for flutter computation
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
nev = 50

# Set system solver options for dynamic problem
systemSolver = aeroSolver isa DynamicStallSolver ?
               create_NewtonRaphson(absTol=1e-7,displayStatus=false) :
               create_NLsolver(method=:newton,linesearch=false,displayStatus=false)

# Strain gauge coordinates (the spar cs is 60 x 2.25 mm)
SGelem = 3          # element number where the strain gauge is located (see Fig. 12 of https://doi.org/10.2514/1.J060621: FBG^1_3 is @ 127 mm spanwise, but wing root is at 32 mm -> actual spanwise position is 127 - 32 = 95 mm which is at the midpoint of element 3)
x2SG = 50e-3/2      # chordwise position: close to LE
x3SG = 2.25e-3/2    # height: on top

# Case ID
idString = string("Theta", round(rad2deg(θ),digits=1), "_Ui", Ui, "_Uf", Uf, "_", tipLossType, "_", aeroSolver.name)

# Set paths
relPath = "/dev/Pazy/S0/BeamRod/outputs/PazyWingBeamRodAirspeedSweep"
absPath = string(pwd(),relPath)
mkpath(absPath)

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
    flutterModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=true,hasTipCorrection=true,tipLossType=tipLossType,pitchAngle=θ,airspeed=U,mesh=mesh,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod))
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
flutterOnsetSpeeds,flutterOnsetFreqs,_,flutterOnsetTipOOPs,flutterOnsetTipAoAs,flutterOnsetEpss,flutterOffsetSpeeds,flutterOffsetFreqs,_,flutterOffsetTipOOPs,flutterOffsetTipAoAs,flutterOffsetEpss = find_flutter(URange,freqsOsc,dampsOsc,tipOOP_f,tipAoA_f,epsSG_f)
if !isempty(flutterOnsetSpeeds)
    flutterOnsetSpeed = flutterOnsetSpeeds[1]
    flutterOnsetFreq = flutterOnsetFreqs[1]
    flutterOnsetTipOOP = flutterOnsetTipOOPs[1]
    flutterOnsetTipAoA = flutterOnsetTipAoAs[1]
    flutterOnsetEps = flutterOnsetEpss[1]
    println("Flutter onset speed = $(round(flutterOnsetSpeed,digits=1)) m/s, frequency = $(round(rads2Hz(flutterOnsetFreq),digits=1)) Hz, tip OOP = $(round(flutterOnsetTipOOP/L*100,digits=1))% semispan, tip AoA = $(round(rad2deg(flutterOnsetTipAoA),digits=1)) deg, 1/6-span strains = $(round(Int,flutterOnsetEps*1e6)) microns")
else
    println("Flutter not found")
end
if !isempty(flutterOffsetSpeeds)
    flutterOffsetSpeed = flutterOffsetSpeeds[1]
    flutterOffsetFreq = flutterOffsetFreqs[1]
    flutterOffsetTipOOP = flutterOffsetTipOOPs[1]
    flutterOffsetTipAoA = flutterOffsetTipAoAs[1]
    flutterOffsetEps = flutterOffsetEpss[1]
    println("Flutter offset speed = $(round(flutterOffsetSpeed,digits=1)) m/s, frequency = $(round(rads2Hz(flutterOffsetFreq),digits=1)) Hz, tip OOP = $(round(flutterOffsetTipOOP/L*100,digits=1))% semispan, tip AoA = $(round(rad2deg(flutterOffsetTipAoA),digits=1)) deg, 1/6-span strains = $(round(Int,flutterOffsetEps*1e6)) microns")
end

# Plots
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
savefig(plt_Vgf,string(absPath,"/PazyWingBeamRodAirspeedSweep_",idString,"_Vgf.pdf"))

## --- Solve dynamic problem ---
println("Solving dynamic problem")

# Model
dynModel = first(create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=true,pitchAngle=θ,airspeed=U,hasTipCorrection=true,tipLossType=tipLossType,tipLossFunctionIsAirspeedDependent=true,mesh=mesh,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod))

# Create and solve steady problem for initial airspeed
steadyProblem = create_SteadyProblem(model=dynModel)
solve!(steadyProblem)

# Create and solve dynamic problem
dynamicProblem = create_DynamicProblem(model=dynModel,integrationMethod=integrationMethod,finalTime=tf,Δt=Δt,systemSolver=systemSolver,x0=steadyProblem.systemArrays.x,trackingFrequency=trackingFrequency)
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
epsSG_LE = (ϵ11SG .- κ2SG*x3SG .- κ3SG*x2SG)
epsSG_TE = (ϵ11SG .- κ2SG*x3SG .- κ3SG*-x2SG)

# Flutter indices
jF = [findfirst(ti -> U(ti) >= speed, t) for speed in [flutterOnsetSpeed, flutterOffsetSpeed]]

## Plots
# Animation
anim = plot_dynamic_deformation(dynamicProblem,follow=:basisA,plotFrequency=round(Int,1e-2/Δt),fps=24,plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/PazyWingBeamRodAirspeedSweep_",idString,"_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
showAirspeedInPlots = true
interactivePlots = false
colors = cgrad(:rainbow, 2, categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
ms = 6
msw = 1
windowFFT = 2^11
if interactivePlots
    plotlyjs()
else
    gr()
end

# Airspeed
plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
plot!(t, U.(t), lw=lw, c=:black, label=false)
scatter!(t[jF], U.(t[jF]), marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
display(plt_U)
savefig(string(absPath,"/PazyWingBeamRodAirspeedSweep_",idString,"_U.pdf"))

# Tip OOP displacement
plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(t, tipOOP./L*100, lw=lw, c=colors[1], label=false)
scatter!(t[jF], [flutterOnsetTipOOP, flutterOffsetTipOOP]/L*100, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
if showAirspeedInPlots && !interactivePlots
    plt_tipOOP_twin = twinx(plt_tipOOP)
    plot!(plt_tipOOP_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t))
    plot!(plt_tipOOP_twin, t, U.(t), lw=lw, c=:black, label=false)
end
display(plt_tipOOP)
savefig(string(absPath,"/PazyWingBeamRodAirspeedSweep_",idString,"_tipOOP.pdf"))

# Tip AoA
plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs))
plot!(t, rad2deg.(tipAoA), lw=lw, c=colors[1], label=false)
scatter!(t[jF], rad2deg.([flutterOnsetTipAoA, flutterOffsetTipAoA]), marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
if showAirspeedInPlots && !interactivePlots
    plt_tipAOA_twin = twinx(plt_tipAOA)
    plot!(plt_tipAOA_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t))
    plot!(plt_tipAOA_twin, t, U.(t), lw=lw, c=:black, label=false)
end
display(plt_tipAOA)
savefig(string(absPath,"/PazyWingBeamRodAirspeedSweep_",idString,"_tipAoA.pdf"))

# 1/6-span strains
plt_eps = plot(xlabel="Time [s]", ylabel="1/6-span strains (\$\\mu\$)", xlims=[0,tf], tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(t, epsSG_LE*1e6, lw=lw, c=colors[1], label=false)
scatter!(t[jF], [flutterOnsetEps, flutterOffsetEps]*1e6, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
if showAirspeedInPlots && !interactivePlots
    plt_eps_twin = twinx(plt_eps)
    plot!(plt_eps_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=extrema(t))
    plot!(plt_eps_twin, t, U.(t), lw=lw, c=:black, label=false)
end
display(plt_eps)
savefig(string(absPath,"/PazyWingBeamRodAirspeedSweep_",idString,"_eps.pdf"))

# LE strains spectrogram
f_eps, t_eps, S_eps = spectrogram(t, epsSG_LE*1e6, window=windowFFT)
U_vec = U.(t)
U_eps = Uf > Ui ? 
        U_vec[searchsortedfirst.(Ref(t), t_eps)] :
        reverse(U_vec[searchsortedfirst.(Ref(t), t_eps)])
plt_eps_spectrum = heatmap(U_eps, f_eps, S_eps; ylims=[0,100], title="1/6-span strains spectrogram", xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
plot!([U_vec[1], U_vec[end]], rads2Hz.([flutterOnsetFreq, flutterOnsetFreq]), ls=:dash, c=:black, lw=lw, label=false)
plot!([U_vec[1], U_vec[end]], 2*rads2Hz.([flutterOnsetFreq, flutterOnsetFreq]), ls=:dash, c=:black, lw=lw, label=false)
plot!([U_vec[1], U_vec[end]], 3*rads2Hz.([flutterOnsetFreq, flutterOnsetFreq]), ls=:dash, c=:black, lw=lw, label=false)
display(plt_eps_spectrum)
savefig(string(absPath,"/PazyWingBeamRodAirspeedSweep_",idString,"_epsSpectrogram.pdf"))

println("Finished PazyWingBeamRodAirspeedSweep.jl")