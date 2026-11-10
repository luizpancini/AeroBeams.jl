using AeroBeams, Random

# Select the flutter test to be reproduced (1, 2 or 3; see Fig. 4 of https://doi.org/10.2514/1.J060621)
test = 1

# Time variables [s]
Δt = 1e-4
trackingFrequency = clamp(round(Int,1e-3/Δt), 1, 5)

# Test data: start time, root AoA before and after the change, time when the AoA change starts, and digitized airspeed profile points to be interpolated by cubic splines
if test == 1
    tStart = 45
    AoAs = [7,5]
    tAoAChange = 90
    tPoints = [45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140]
    UPoints = [30.0,34.5,38.2,41.5,44.5,46.8,48.6,50.0,51.0,51.5,51.7,51.8,49.5,44.5,40.5,36.8,33.5,30.2,27.5,25.0]
elseif test == 2
    tStart = 10
    AoAs = [5,3]
    tAoAChange = 75
    tPoints = [10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115]
    UPoints = [26.0,30.5,34.8,38.8,42.2,45.0,47.5,49.8,51.5,52.8,53.7,54.4,54.7,54.2,49.8,45.2,40.5,36.5,32.8,29.2,26.0,23.5]
elseif test == 3
    tStart = 10
    AoAs = [3,1]
    tAoAChange = 70
    tPoints = [10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105]
    UPoints = [28.5,33.2,37.5,41.8,45.5,49.0,52.2,54.8,56.8,58.2,59.3,59.7,59.6,57.5,52.8,47.8,43.0,38.5,34.5,31.0]
else
    error("test must be 1, 2 or 3")
end
tEndTest = last(tPoints)
tf = tEndTest - tStart
tChange = round((tAoAChange - tStart)/Δt)*Δt # time of the AoA step, on the time grid
tPeak = tPoints[argmax(UPoints)] - tStart # time of maximum airspeed (starts the downsweep)

# Root AoAs [rad]: the test is split in two dynamic problems, one for each AoA, the second continuing the first
θs = deg2rad.(AoAs)

# Profile for airspeed [m/s]: digitized from the test, with a small random noise over the maximum value to induce flutter
σn = 1e-4 * maximum(UPoints)
rng = Xoshiro(1234)
smoothingWindow = Int(1e2)
noiseRaw = σn .* randn(rng, round(Int, tf/Δt) + 1)
noise = moving_average(noiseRaw, smoothingWindow)

# Cubic spline through the points (x,y), returning its value and derivative as functions. The ends are natural (zero curvature) unless the slope is given (clamped)
function cubic_spline(x,y; slopeStart=nothing, slopeEnd=nothing)
    n = length(x)
    h = diff(x)
    # Linear system for the second derivatives M at the knots
    A = zeros(n,n)
    rhs = zeros(n)
    for i in 2:n-1
        A[i,i-1:i+1] = [h[i-1], 2*(h[i-1]+h[i]), h[i]]
        rhs[i] = 6*((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
    end
    if isnothing(slopeStart)
        A[1,1] = 1
    else
        A[1,1:2] = [2*h[1], h[1]]
        rhs[1] = 6*((y[2]-y[1])/h[1] - slopeStart)
    end
    if isnothing(slopeEnd)
        A[n,n] = 1
    else
        A[n,n-1:n] = [h[n-1], 2*h[n-1]]
        rhs[n] = 6*(slopeEnd - (y[n]-y[n-1])/h[n-1])
    end
    M = A \ rhs
    # Evaluation of the spline and its derivative
    function locate(t)
        tt = clamp(t, x[1], x[end])
        i = clamp(searchsortedlast(x,tt), 1, n-1)
        return tt, i, (x[i+1]-tt)/h[i], (tt-x[i])/h[i]
    end
    f = function(t)
        tt,i,a,b = locate(t)
        return a*y[i] + b*y[i+1] + ((a^3-a)*M[i] + (b^3-b)*M[i+1])*h[i]^2/6
    end
    df = function(t)
        tt,i,a,b = locate(t)
        return (y[i+1]-y[i])/h[i] - (3*a^2-1)/6*h[i]*M[i] + (3*b^2-1)/6*h[i]*M[i+1]
    end
    return f,df
end

# Nominal airspeed (in test time): one cubic spline for each part of the test (before and after the AoA change), joined at the time of the change with the value and slope of a spline through all points, so that the profile is smooth
Uglobal,dUglobal = cubic_spline(tPoints,UPoints)
UChange,dUChange = Uglobal(tAoAChange),dUglobal(tAoAChange)
before = tPoints .< tAoAChange
after = tPoints .> tAoAChange
Ubefore,_ = cubic_spline([tPoints[before]; tAoAChange], [UPoints[before]; UChange]; slopeEnd=dUChange)
Uafter,_ = cubic_spline([tAoAChange; tPoints[after]], [UChange; UPoints[after]]; slopeStart=dUChange)
Ubase(t) = t+tStart <= tAoAChange ? Ubefore(t+tStart) : Uafter(t+tStart)

# Airspeed [m/s], with noise
U = t -> Ubase(t) + noise[clamp(round(Int, t/Δt) + 1, 1, length(noise))]

# Time integration method for the dynamic solver
integrationMethod = GeneralizedAlpha(ρ∞=0.8)

# Tip correction function type
tipLossType = "Exponential"

# Mesh (defines number of elements): :nominal or :refined
mesh = :nominal

# Number of elements in the wingtip rod
nElementsWingtipRod = 12

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Geometric properties
nElem,L,chord,normSparPos,nodalPositions = geometrical_properties_Pazy(mesh=mesh)

# Airspeed range for flutter problem [m/s]
URange = collect(0.5:0.5:70)

# Number of modes for flutter computation
nModes = 3

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
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
idString = string("Test", test, "_", tipLossType, "_", aeroSolver.name)

# Set paths
relPath = "/dev/Pazy/S0/BeamRod/outputs/PazyWingBeamRodAirspeedAoASweepTest"
absPath = string(pwd(),relPath)
mkpath(absPath)

## --- Solve flutter problem for each AoA ---
# Initialize outputs
flutterProblem = Array{EigenProblem}(undef,length(θs),length(URange))
eigensols = Array{EigenSolution}(undef,length(θs),length(URange))
tipOOP_f = Array{Float64}(undef,length(θs),length(URange))
tipAoA_f = Array{Float64}(undef,length(θs),length(URange))
epsSG_f = Array{Float64}(undef,length(θs),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef,length(θs))
dampsOsc = Array{Matrix{Float64}}(undef,length(θs))
flutterOnsetSpeeds_AoA = Array{Float64}(undef,length(θs))
flutterOnsetFreqs_AoA = Array{Float64}(undef,length(θs))
flutterOnsetOOPs_AoA = Array{Float64}(undef,length(θs))
flutterOnsetAoAs_AoA = Array{Float64}(undef,length(θs))
flutterOnsetEpss_AoA = Array{Float64}(undef,length(θs))
flutterOffsetSpeeds_AoA = Array{Float64}(undef,length(θs))
flutterOffsetFreqs_AoA = Array{Float64}(undef,length(θs))
flutterOffsetOOPs_AoA = Array{Float64}(undef,length(θs))
flutterOffsetAoAs_AoA = Array{Float64}(undef,length(θs))
flutterOffsetEpss_AoA = Array{Float64}(undef,length(θs))

# Loop over AoAs
for (k,θ) in enumerate(θs)
    println("Solving flutter problem for AoA = $(AoAs[k]) deg")
    # Sweep airspeed
    for (j,Uj) in enumerate(URange)
        # Model
        flutterModel = first(create_Pazy(aeroSolver=aeroSolver,mesh=mesh,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod,airfoil=airfoil,upright=true,hasTipCorrection=true,tipLossType=tipLossType,pitchAngle=θ,airspeed=Uj))
        # Set initial guess solution as previously available one
        x0 = j == 1 ? zeros(0) : flutterProblem[k,j-1].systemArrays.x
        # Create and solve problem
        flutterProblem[k,j] = create_EigenProblem(model=flutterModel,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(flutterProblem[k,j])
        # Eigensolution
        eigensols[k,j] = flutterProblem[k,j].eigensolution
        # Kinematics
        states = flutterProblem[k,j].loadstepHistory.states
        tipOOP_f[k,j] = -states[end][nElem].nodal.u2[1]
        tipAoA_f[k,j] = flutterProblem[k,j].loadstepHistory.aeroVariables[end][nElem].unsteadyKin.αₑ
        ϵ11SG_f = states[end][SGelem].midpoint.γ[1]
        κ2SG_f = states[end][SGelem].midpoint.κ[2]
        κ3SG_f = states[end][SGelem].midpoint.κ[3]
        epsSG_f[k,j] = (ϵ11SG_f .- κ2SG_f*x3SG .- κ3SG_f*x2SG)
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols[k,:])
    # Get frequencies and dampings of oscillatory modes
    freqsOsc[k],dampsOsc[k],_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(3.75), fmax=Hz2rads(50), maxTracks=nModes, dropOutsideBounds=true)
    # Flutter onset and offset speeds, frequencies and kinematics
    onsetSpeeds,onsetFreqs,_,onsetTipOOPs,onsetTipAoAs,onsetEpss,offsetSpeeds,offsetFreqs,_,offsetTipOOPs,offsetTipAoAs,offsetEpss = find_flutter(URange,freqsOsc[k],dampsOsc[k],tipOOP_f[k,:],tipAoA_f[k,:],epsSG_f[k,:])
    if !isempty(onsetSpeeds)
        flutterOnsetSpeeds_AoA[k] = onsetSpeeds[1]
        flutterOnsetFreqs_AoA[k] = onsetFreqs[1]
        flutterOnsetOOPs_AoA[k] = onsetTipOOPs[1]
        flutterOnsetAoAs_AoA[k] = onsetTipAoAs[1]
        flutterOnsetEpss_AoA[k] = onsetEpss[1]
        println("Flutter onset speed = $(round(onsetSpeeds[1],digits=1)) m/s, frequency = $(round(rads2Hz(onsetFreqs[1]),digits=1)) Hz, tip OOP = $(round(onsetTipOOPs[1]/L*100,digits=1))% semispan, tip AoA = $(round(rad2deg(onsetTipAoAs[1]),digits=1)) deg, 1/6-span strains = $(round(Int,onsetEpss[1]*1e6)) microns")
    else
        println("Flutter not found")
    end
    if !isempty(offsetSpeeds)
        flutterOffsetSpeeds_AoA[k] = offsetSpeeds[1]
        flutterOffsetFreqs_AoA[k] = offsetFreqs[1]
        flutterOffsetOOPs_AoA[k] = offsetTipOOPs[1]
        flutterOffsetAoAs_AoA[k] = offsetTipAoAs[1]
        flutterOffsetEpss_AoA[k] = offsetEpss[1]
        println("Flutter offset speed = $(round(offsetSpeeds[1],digits=1)) m/s, frequency = $(round(rads2Hz(offsetFreqs[1]),digits=1)) Hz, tip OOP = $(round(offsetTipOOPs[1]/L*100,digits=1))% semispan, tip AoA = $(round(rad2deg(offsetTipAoAs[1]),digits=1)) deg, 1/6-span strains = $(round(Int,offsetEpss[1]*1e6)) microns")
    end
end

# Plots
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 14
lw = 2
gr()
plt_Vf = plot(ylabel="Frequency [Hz]", xlims=extrema(URange), ylims=[0,50], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URange), ylims=[-0.2,0.1], xticks=vcat(0:10:100), tickfont=font(ts), guidefont=font(fs))
plot!(plt_Vg, URange, zeros(length(URange)), c=:black, lw=lw, ls=:dash, label=false)
for (k,θ) in enumerate(θs)
    for mode in 1:nModes
        plot!(plt_Vf, URange, rads2Hz.(freqsOsc[k][mode,:]), c=modeColors[mode], alpha=0.5+0.5*(k-1), lw=lw, label=false)
        plot!(plt_Vg, URange, dampsOsc[k][mode,:]./freqsOsc[k][mode,:], c=modeColors[mode], alpha=0.5+0.5*(k-1), lw=lw, label=false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(plt_Vgf,string(absPath,"/PazyWingBeamRodAirspeedAoASweepTest_",idString,"_Vgf.pdf"))

## --- Solve dynamic problems ---
println("Solving dynamic problems")

# Function to create the dynamic model for a given AoA. The root pitch angle is kept null and the AoA is set through the flight direction. The tip loss is evaluated at the AoA
create_dynamic_model(θ) = first(create_Pazy(aeroSolver=aeroSolver,mesh=mesh,beamWingtipRod=true,nElementsWingtipRod=nElementsWingtipRod,airfoil=airfoil,upright=true,pitchAngle=0,flightDirection=[sin(θ),cos(θ),0],tipLossPitchAngle=θ,airspeed=U,hasTipCorrection=true,tipLossType=tipLossType,tipLossFunctionIsAirspeedDependent=true))

# Create and solve steady problem for initial airspeed and first AoA
steadyProblem = create_SteadyProblem(model=create_dynamic_model(θs[1]))
solve!(steadyProblem)

# Initial and final times of each problem (the second starts where the first ends, with the new AoA)
initialTime = [0.0, tChange]
finalTime = [tChange, tf]

# Solve the dynamic problem of each AoA, the second continuing from the last point of the first
dynamicProblem = Array{DynamicProblem}(undef,2)
for (k,θ) in enumerate(θs)
    println("Solving dynamic problem for AoA = $(AoAs[k]) deg")
    # Model
    dynModel = create_dynamic_model(θ)
    # Initial solution and reference dynamic problem
    x0 = k == 1 ? steadyProblem.systemArrays.x : dynamicProblem[k-1].systemArrays.x
    refDynamicProblem = k == 1 ? nothing : dynamicProblem[k-1]
    # Create and solve dynamic problem
    dynamicProblem[k] = create_DynamicProblem(model=dynModel,integrationMethod=integrationMethod,initialTime=initialTime[k],finalTime=finalTime[k],Δt=Δt,systemSolver=systemSolver,x0=x0,refDynamicProblem=refDynamicProblem,trackingFrequency=trackingFrequency)
    solve!(dynamicProblem[k])
end

# Unpack numerical solution (skipping the initial sample of the continuation problem to avoid duplicated time points)
sols = [problem.solutionHistory for problem in dynamicProblem]
idx = [eachindex(sols[1].t), 2:length(sols[2].t)]
t = vcat([sols[k].t[idx[k]] for k in 1:2]...)
get_all(f) = vcat([[f(sols[k],i) for i in idx[k]] for k in 1:2]...)
tipAoA = get_all((sol,i) -> sol.aeroVariables[i][nElem].unsteadyKin.αₑ)
tipOOP = get_all((sol,i) -> -sol.states[i][nElem].nodal.u2[1])
ϵ11SG = get_all((sol,i) -> sol.states[i][SGelem].midpoint.γ[1])
κ2SG = get_all((sol,i) -> sol.states[i][SGelem].midpoint.κ[2])
κ3SG = get_all((sol,i) -> sol.states[i][SGelem].midpoint.κ[3])
eps_SG = ϵ11SG .- κ2SG*x3SG .- κ3SG*x2SG

# Root AoA over time [rad]
rootAoA = [ti < tChange ? θs[1] : θs[2] for ti in t]

# Time in test time [s]
tTest = t .+ tStart

# Flutter indices along the time history (onset and offset at both AoA)
jFon1 = findfirst(ti -> ti <= tChange && Ubase(ti) >= flutterOnsetSpeeds_AoA[1], t)
jFon2 = findfirst(ti -> ti >= tChange && Ubase(ti) <= flutterOnsetSpeeds_AoA[2], t)
jFoff1 = findfirst(ti -> ti <= tChange && Ubase(ti) >= flutterOffsetSpeeds_AoA[1], t)
jFoff2 = findfirst(ti -> ti >= tChange && Ubase(ti) <= flutterOffsetSpeeds_AoA[2], t)
isnothing(jFoff2) && (jFoff2 = findfirst(ti -> ti >= tChange, t))
flutterOnIdx = filter(!isnothing, [jFon1, jFon2])
flutterOffIdx = filter(!isnothing, [jFoff1, jFoff2])
flutterIdx = sort(union(flutterOnIdx, flutterOffIdx))

## Plots
# Animation
for k in 1:2
    anim = plot_dynamic_deformation(dynamicProblem[k],follow=:basisA,plotFrequency=round(Int,1e-2/Δt),fps=24,plotLimits=([-L/2,L/2],[-L/2,L/2],[0,1.05*L]),save=true,savePath=string(relPath,"/PazyWingBeamRodAirspeedAoASweepTest_",idString,"_deformation_AoA",AoAs[k],".gif"),displayProgress=true)
    display(anim)
end

# Plot configurations
showAirspeedInPlots = true
interactivePlots = false
colors = cgrad(:rainbow, 2, categorical=true)
epsColors = cgrad(:rainbow, 5, categorical=true)
ts = 10
fs = 16
lfs = 8
lw = 2
ms = 6
msw = 1
fftWindow = 2^13
fftOverlap = 0.5
xlimsTest = [tStart,tEndTest]
if interactivePlots
    plotlyjs()
else
    gr()
end

# Airspeed and root AoA
plt_U = plot(xlabel="Time [s]", ylabel="Airspeed [m/s]", xlims=xlimsTest, tickfont=font(ts), guidefont=font(fs))
plot!(tTest, U.(t), lw=lw, c=:black, label=false)
for jF in flutterIdx
    scatter!([tTest[jF]], [U(t[jF])], marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
end
if showAirspeedInPlots && !interactivePlots
    plt_U_twin = twinx(plt_U)
    plot!(plt_U_twin, ylabel="Root AoA [deg]", tickfont=font(ts), guidefont=font(fs), xlims=xlimsTest, ylims=[0,ceil(maximum(AoAs))+1])
    plot!(plt_U_twin, tTest, rad2deg.(rootAoA), lw=lw, c=colors[2], label=false)
end
display(plt_U)
savefig(string(absPath,"/PazyWingBeamRodAirspeedAoASweepTest_",idString,"_U.pdf"))

# Tip OOP displacement
plt_tipOOP = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]", xlims=xlimsTest, tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(tTest, tipOOP./L*100, lw=lw, c=colors[1], label=false)
scatter!(tTest[flutterOnIdx], flutterOnsetOOPs_AoA/L*100, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
scatter!(tTest[flutterOffIdx], flutterOffsetOOPs_AoA/L*100, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
if showAirspeedInPlots && !interactivePlots
    plt_tipOOP_twin = twinx(plt_tipOOP)
    plot!(plt_tipOOP_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=xlimsTest)
    plot!(plt_tipOOP_twin, tTest, U.(t), lw=lw, c=:black, label=false)
end
display(plt_tipOOP)
savefig(string(absPath,"/PazyWingBeamRodAirspeedAoASweepTest_",idString,"_tipOOP.pdf"))

# Tip AoA
plt_tipAOA = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]", xlims=xlimsTest, tickfont=font(ts), guidefont=font(fs))
plot!(tTest, rad2deg.(tipAoA), lw=lw, c=colors[1], label=false)
scatter!(tTest[flutterOnIdx], rad2deg.(flutterOnsetAoAs_AoA), marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
scatter!(tTest[flutterOffIdx], rad2deg.(flutterOffsetAoAs_AoA), marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
if showAirspeedInPlots && !interactivePlots
    plt_tipAOA_twin = twinx(plt_tipAOA)
    plot!(plt_tipAOA_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=xlimsTest)
    plot!(plt_tipAOA_twin, tTest, U.(t), lw=lw, c=:black, label=false)
end
display(plt_tipAOA)
savefig(string(absPath,"/PazyWingBeamRodAirspeedAoASweepTest_",idString,"_tipAoA.pdf"))

# Strains near 1/6-span
plt_eps = plot(xlabel="Time [s]", ylabel="1/6-span strains (\$\\mu\\epsilon\$)", xlims=xlimsTest, tickfont=font(ts), guidefont=font(fs), legend=:best, legendfontsize=lfs)
plot!(tTest, eps_SG*1e6, lw=lw, c=colors[1], label=false)
scatter!(tTest[flutterOnIdx], flutterOnsetEpss_AoA*1e6, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
scatter!(tTest[flutterOffIdx], flutterOffsetEpss_AoA*1e6, marker=:circle, ms=ms, msw=msw, mc=:white, msc=:red, label=false)
if showAirspeedInPlots && !interactivePlots
    plt_eps_twin = twinx(plt_eps)
    plot!(plt_eps_twin, ylabel="Airspeed [m/s]", tickfont=font(ts), guidefont=font(fs), xlims=xlimsTest)
    plot!(plt_eps_twin, tTest, U.(t), lw=lw, c=:black, label=false)
end
display(plt_eps)
savefig(string(absPath,"/PazyWingBeamRodAirspeedAoASweepTest_",idString,"_eps.pdf"))

# 1/6-span strains spectrogram over the whole test, with the flutter frequency and its harmonics at the first AoA (before the change) and at the second AoA (after it)
f_eps, t_eps, S_eps = spectrogram(t, eps_SG*1e6; window=fftWindow, overlap=fftOverlap)
plt_eps_spectrum = heatmap(t_eps .+ tStart, f_eps, S_eps; ylims=[0,100], xlims=xlimsTest, title="1/6-span strains spectrogram", xlabel="Time [s]", ylabel="Frequency [Hz]", c=:rainbow, colorbar_title="Power [dB]", aspect_ratio=:auto, tickfont=font(ts), guidefont=font(fs))
plot!([tAoAChange, tAoAChange], [0,100], ls=:dot, c=:white, lw=lw, label=false)
display(plt_eps_spectrum)
savefig(string(absPath,"/PazyWingBeamRodAirspeedAoASweepTest_",idString,"_epsSpectrogram.pdf"))

println("Finished PazyWingBeamRodAirspeedAoASweepTest.jl")