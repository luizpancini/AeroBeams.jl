using AeroBeams, LinearInterpolations, DelimitedFiles

# Aerodynamic solver and derivatives method
aeroSolver = Indicial()
derivationMethod = AD()

# Altitude
h = 20e3

# Gravity
g = 0

# Pitch angle
θ = 0

# Discretization
nElem = 16

# Initialize model
SMWFlutterTipDispRange,L = create_SMW(aeroSolver=aeroSolver,derivationMethod=derivationMethod,θ=θ,nElem=nElem,altitude=h,g=g)

# Set system solver options
σ0 = 1
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Set tip force and airspeed ranges, and initialize outputs
F3Range = collect(-44:2:44)
URange = collect(13:0.5:35)
untrackedFreqs = Array{Vector{Float64}}(undef,length(F3Range),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(F3Range),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(F3Range),length(URange))
freqs = Array{Vector{Float64}}(undef,length(F3Range),length(URange))
damps = Array{Vector{Float64}}(undef,length(F3Range),length(URange))
tip_u3 = Array{Float64}(undef,length(F3Range),length(URange))
flutterSpeed = Array{Float64}(undef,length(F3Range))
flutterFreq = Array{Float64}(undef,length(F3Range))
flutterMode = Array{Int64}(undef,length(F3Range))
flutterTipDisp = Array{Float64}(undef,length(F3Range))

# Set number of vibration modes
nModes = 5

# Sweep tip force
for (i,F3) in enumerate(F3Range)
    # Update model with tip force
    SMWFlutterTipDispRange,_ = create_SMW(aeroSolver=aeroSolver,derivationMethod=derivationMethod,θ=θ,nElem=nElem,altitude=h,g=g,tipF3=F3)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for F3 = $F3 N, U = $U m/s")
        # Update velocity of basis A (and update model)
        set_motion_basis_A!(model=SMWFlutterTipDispRange,v_A=[0;U;0])
        # Create and solve problem
        problem = create_EigenProblem(model=SMWFlutterTipDispRange,systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
        solve!(problem)
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = problem.frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(problem.dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = problem.eigenvectorsOscillatoryCplx
        # Tip OOP displacement
        tip_u3[i,j] = problem.nodalStatesOverσ[end][nElem].u_n2[3]
    end
    # Frequencies and dampings after mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Flutter speeds, frequencies and tip displacements of modes at current tip force
    dampsCurrentF3 = Array{Vector{Float64}}(undef,nModes)
    freqsCurrentF3 = Array{Vector{Float64}}(undef,nModes)
    flutterSpeedOfMode = Array{Float64}(undef,nModes)
    flutterFreqOfMode = Array{Float64}(undef,nModes)
    flutterTipDispOfMode = Array{Float64}(undef,nModes)
    for mode in 1:nModes
        dampsCurrentF3[mode] = [damps[i,j][mode] for j in eachindex(URange)]
        freqsCurrentF3[mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        indexInstability = findfirst(x->x>0,dampsCurrentF3[mode])
        if isnothing(indexInstability)
            flutterSpeedOfMode[mode] = NaN
            flutterFreqOfMode[mode] = NaN
            flutterTipDispOfMode[mode] = NaN
            continue
        end
        flutterSpeedOfMode[mode] = interpolate(dampsCurrentF3[mode][indexInstability-1:indexInstability],URange[indexInstability-1:indexInstability],0)
        flutterFreqOfMode[mode] = interpolate(dampsCurrentF3[mode][indexInstability-1:indexInstability],freqsCurrentF3[mode][indexInstability-1:indexInstability],0)
        flutterTipDispOfMode[mode] = interpolate(dampsCurrentF3[mode][indexInstability-1:indexInstability],tip_u3[i,indexInstability-1:indexInstability],0)
    end
    # Set flutter speed as the greatest (for compatibility with reference), and get corresponding flutter mode and frequency, and tip displacement
    flutterSpeed[i] = maximum(filter(!isnan,flutterSpeedOfMode))
    flutterMode[i] = findfirst(x->x==flutterSpeed[i],flutterSpeedOfMode)
    flutterFreq[i] = flutterFreqOfMode[flutterMode[i]]
    flutterTipDisp[i] = flutterTipDispOfMode[flutterMode[i]]
end

# Load reference data
flutterSpeedRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "SMW", "flutterSpeedVsTipDisp.txt"))
flutterFreqRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "SMW", "flutterFreqVsTipDisp.txt"))
flutterSpeedVsDispRef = readdlm(joinpath(dirname(@__DIR__), "referenceData", "SMW", "flutterSpeedVsTipDispFull.txt"))

println("Finished SMWFlutterTipDispRange.jl")