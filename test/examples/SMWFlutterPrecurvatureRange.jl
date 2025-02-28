using AeroBeams, LinearAlgebra, LinearInterpolations, DelimitedFiles

# Option for mode tracking
modeTracking = true

# Aerodynamic solver and derivatives method
aeroSolver = Indicial()
derivationMethod = AD()

# Altitude
h = 20e3

# Gravity
g = 0

# Discretization
nElem = 16

# Pitch angle
θ = 0

# Set system solver options (limit initial load factor)
σ0 = 1
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Number of vibration modes
nModes = 5

# Set precurvature, tip force and airspeed ranges, and initialize outputs
kRange = collect(0:0.018:0.018)
F3Range = collect(0:2:40)
URange = vcat(1e-3,collect(20:0.5:35))
untrackedFreqs = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(kRange),length(F3Range),length(URange))
freqs = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),length(URange))
damps = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),length(URange))
tip_u3 = Array{Float64}(undef,length(kRange),length(F3Range),length(URange))
flutterSpeed = Array{Float64}(undef,length(kRange),length(F3Range))
flutterFreq = Array{Float64}(undef,length(kRange),length(F3Range))
flutterMode = Array{Int64}(undef,length(kRange),length(F3Range))
flutterTipDisp = Array{Float64}(undef,length(kRange),length(F3Range))
modeDampings = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),nModes)
modeFrequencies = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),nModes)

SMWFlutterPrecurvatureRange = Array{Model}(undef,length(kRange),length(F3Range),length(URange))

# Sweep wing precurvature
for (ki,k) in enumerate(kRange)
    # Sweep tip force
    for (i,F3) in enumerate(F3Range)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for k=$k, F3 = $F3 N, U = $U m/s")
            # Update model
            SMWFlutterPrecurvatureRange[ki,i,j],_ = create_SMW(aeroSolver=aeroSolver,derivationMethod=derivationMethod,θ=θ*π/180,k2=k,nElem=nElem,altitude=h,g=g,tipF3=F3,airspeed=U)
            # Create and solve problem
            problem = create_EigenProblem(model=SMWFlutterPrecurvatureRange[ki,i,j],systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
            solve!(problem)
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[ki,i,j] = problem.frequenciesOscillatory
            untrackedDamps[ki,i,j] = round_off!(problem.dampingsOscillatory,1e-8)
            untrackedEigenvectors[ki,i,j] = problem.eigenvectorsOscillatoryCplx
            # Tip OOP displacement
            tip_u3[ki,i,j] = problem.nodalStatesOverσ[end][nElem].u_n2[3]
        end
        # Frequencies and dampings after mode tracking
        if modeTracking
            freqs[ki,i,:],damps[ki,i,:],_ = mode_tracking(URange,untrackedFreqs[ki,i,:],untrackedDamps[ki,i,:],untrackedEigenvectors[ki,i,:])
        else
            freqs[ki,i,:],damps[ki,i,:] = untrackedFreqs[ki,i,:],untrackedDamps[ki,i,:]
        end
        # Separate frequencies and dampings by mode
        for mode in 1:nModes
            modeFrequencies[ki,i,mode] = [freqs[ki,i,j][mode] for j in eachindex(URange)]
            modeDampings[ki,i,mode] = [damps[ki,i,j][mode] for j in eachindex(URange)]
        end
        # Flutter speeds, frequencies and tip displacements of modes at current tip force
        dampsCurrentF3 = Array{Vector{Float64}}(undef,nModes)
        freqsCurrentF3 = Array{Vector{Float64}}(undef,nModes)
        flutterSpeedOfMode = Array{Float64}(undef,nModes)
        flutterFreqOfMode = Array{Float64}(undef,nModes)
        flutterTipDispOfMode = Array{Float64}(undef,nModes)
        for mode in 1:nModes
            dampsCurrentF3[mode] = [damps[ki,i,j][mode] for j in eachindex(URange)]
            freqsCurrentF3[mode] = [freqs[ki,i,j][mode] for j in eachindex(URange)]
            indexInstability = findfirst(x->x>0,dampsCurrentF3[mode])
            if isnothing(indexInstability)
                flutterSpeedOfMode[mode] = NaN
                flutterFreqOfMode[mode] = NaN
                flutterTipDispOfMode[mode] = NaN
                continue
            end
            flutterSpeedOfMode[mode] = interpolate(dampsCurrentF3[mode][indexInstability-1:indexInstability],URange[indexInstability-1:indexInstability],0)
            flutterFreqOfMode[mode] = interpolate(dampsCurrentF3[mode][indexInstability-1:indexInstability],freqsCurrentF3[mode][indexInstability-1:indexInstability],0)
            flutterTipDispOfMode[mode] = interpolate(dampsCurrentF3[mode][indexInstability-1:indexInstability],tip_u3[ki,i,indexInstability-1:indexInstability],0)
        end
        # Set flutter speed as the greatest (for compatibility with reference), and get corresponding flutter mode and frequency, and tip displacement
        flutterSpeed[ki,i] = maximum(filter(!isnan,flutterSpeedOfMode))
        flutterMode[ki,i] = findfirst(x->x==flutterSpeed[ki,i],flutterSpeedOfMode)
        flutterFreq[ki,i] = flutterFreqOfMode[flutterMode[ki,i]]
        flutterTipDisp[ki,i] = flutterTipDispOfMode[flutterMode[ki,i]]
    end
end

# Load reference data
flutterSpeedVsTipLoadk0 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "SMW", "flutterSpeedVsTipLoadk0.txt"))
flutterSpeedVsTipLoadk2 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "SMW", "flutterSpeedVsTipLoadk2.txt"))

println("Finished SMWFlutterPrecurvatureRange.jl")