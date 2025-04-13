using AeroBeams, DelimitedFiles

# Hinge configuration
hingeConfiguration = "free"

# Flare angle range [rad]
ΛRange = π/180*[10,15,20]

# Airspeed range [m/s]
URange = collect(0:0.5:40)

# Pitch angle [rad]
θ = 0*π/180

# Gravity
g = 0

# Stiffness of the spring around the hinge for in-plane bending
kIPBendingHinge = 1e12

# Discretization
nElementsInner = 16
nElementsFFWT = 4

# Tip loss options (assumed, since Healy's analysis uses DLM for aerodynamic)
hasTipCorrection = true
tipLossDecayFactor = 12

# Solution method for hinge constraint
solutionMethod = "appliedMoment"

# System solver
σ0 = 1
maxIter = 100
relTol = 1e-8
NR = create_NewtonRaphson(displayStatus=false,initialLoadFactor=σ0,maximumIterations=maxIter,relativeTolerance=relTol)

# Number of modes
nModes = 5

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(ΛRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(ΛRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(ΛRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(ΛRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(ΛRange),length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,length(ΛRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(ΛRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(ΛRange),nModes)
problem = Array{EigenProblem}(undef,length(ΛRange),length(URange))

# Sweep flare angle
for (i,Λ) in enumerate(ΛRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for Λ=$(round(Int,Λ*180/π)) deg, U=$U m/s")
        # Update model
        model = create_HealyBaselineFFWT(solutionMethod=solutionMethod,hingeConfiguration=hingeConfiguration,flareAngle=Λ,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,g=g,kIPBendingHinge=kIPBendingHinge,nElementsInner=nElementsInner,nElementsFFWT=nElementsFFWT)
        # Create and solve problem
        problem[i,j] = create_EigenProblem(model=model,nModes=nModes,systemSolver=NR,frequencyFilterLimits=[1e-2*U,Inf])
        solve!(problem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = problem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(problem[i,j].dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = problem[i,j].eigenvectorsOscillatoryCplx
    end
    # Apply mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies and damping ratios by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
        modeDampingRatios[i,mode] = modeDampings[i,mode]./modeFrequencies[i,mode]
    end
end

# Load reference data
flare10_mode1_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare10_mode1_damp.txt")
flare10_mode2_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare10_mode2_damp.txt")
flare10_mode1_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare10_mode1_freq.txt")
flare10_mode2_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare10_mode2_freq.txt")

flare15_mode1_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare15_mode1_damp.txt")
flare15_mode2_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare15_mode2_damp.txt")
flare15_mode1_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare15_mode1_freq.txt")
flare15_mode2_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare15_mode2_freq.txt")

flare20_mode1_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare20_mode1_damp.txt")
flare20_mode2_damp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare20_mode2_damp.txt")
flare20_mode1_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare20_mode1_freq.txt")
flare20_mode2_freq = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyBaselineFFWTfreeFlutterFlareRangeURange/flare20_mode2_freq.txt")

println("Finished HealyBaselineFFWTfreeFlutterFlareRangeURange.jl")