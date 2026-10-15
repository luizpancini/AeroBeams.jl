using AeroBeams

# Airspeed range [m/s]
URange = vcat(1e-3,0.5:0.5:60)

# Aerodynamic solver
aeroSolver = Indicial()

# Flag for skin (wing covering)
withSkin = false

# Airfoil section
airfoil = deepcopy(FlatPlate)

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "Exponential"

# Flag for upright position
upright = true

# Pitch angle [rad]
θ = deg2rad(7)

# Gravity
g = 0

# Flag for small angles approximation
smallAngles = true

# Fixed geometrical and discretization properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(URange))
tipOOP = Array{Float64}(undef,length(URange))

# Initial guess solution
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Model
    PazyWingFreqsEvolution,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,withSkin=withSkin,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,upright=upright,pitchAngle=θ,airspeed=U,gravityVector=[0,0,-g],smallAngles=smallAngles)
    # Create and solve problem
    problem = create_EigenProblem(model=PazyWingFreqsEvolution,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(problem)
    # Eigensolution
    eigensols[i] = problem.eigensolution
    # Get OOP displacement at midchord
    tipNodeFinalState = problem.loadstepHistory.states[end][nElem].nodal
    tipTwist = twist_angle(tipNodeFinalState.p2_b, outputUnit=:deg)
    tipOOP[i] = -(tipNodeFinalState.u2[1] - chord*(1/2-normSparPos)*sind(tipTwist))
    # Set initial guess solution for next iteration
    global x0 = problem.systemArrays.x
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(eigensols)
rads2Hz!(freqs)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=eps(), maxTracks=nModes)

# Load reference data
freqsVsU_alpha7 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/freqsVsU_alpha7.txt")
freqsVsDisp_alpha7 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/freqsVsDisp_alpha7.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/PazyWingFreqsEvolution"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 16
lfs = 10
lw = 2
ms = 3
msw = 0
gr()

# Frequencies vs. airspeed
pltU = plot(xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]", title="\$\\alpha_r = 7^\\circ\$", xlims=[0,45], ylims=[0,120], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dash, label="Riso & Cesnik (2023)")
for mode in 1:nModes
    plot!(URange, freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
    plot!(freqsVsU_alpha7[1,:], freqsVsU_alpha7[2,:], lw=lw, ls=:dash, c=:black, label=false)
end
display(pltU)
savefig(string(absPath,"/PazyWingFreqsEvolution_U.pdf"))

# Frequencies vs. tip OOP displacement 
pltDisp = plot(xlabel="Tip displacement [% semispan]", ylabel="Frequency [Hz]", title="\$\\alpha_r = 7^\\circ\$", xlims=[0,45], ylims=[0,120], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dash, label="Riso & Cesnik (2023)")
for mode in 1:nModes
    plot!(tipOOP*100/L, freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
    plot!(freqsVsDisp_alpha7[1,:], freqsVsDisp_alpha7[2,:], lw=lw, ls=:dash, c=:black, label=false)
end
display(pltDisp)
savefig(string(absPath,"/PazyWingFreqsEvolution_disp.pdf"))

println("Finished PazyWingFreqsEvolution.jl")