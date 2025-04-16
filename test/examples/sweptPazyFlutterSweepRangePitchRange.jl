using AeroBeams, LinearInterpolations, DelimitedFiles

# --- Configurations ---
# Flag for tip correction
hasTipCorrectionConfig = [false,true,true,true]
# Tip correction function type
tipLossTypeConfig = ["None","Exponential","VLM-undef","VLM-def"]

# Sweep angle
ΛRange = π/180*vcat(0:10:30)

# Root pitch angle range
θRange = π/180*unique(sort([vcat(0:0.25:0.75)...,vcat(1:0.5:7)...,vcat(0.1,2.35,6.75,4.25,5.75)...]))

# Airspeed range
URange = unique(sort([vcat(30:0.5:60)...,vcat(30:1:120)...]))

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = false

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = false

# Gravity
g = 0

# Number of vibration modes
nModes = 3

# Fixed geometrical and discretization properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Set system solver options
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
modeFrequencies = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),nModes)
modeDampings = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),nModes)
tipOOP = Array{Float64}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
flutterOnsetSpeedsOfMode = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),nModes)
flutterOnsetFreqsOfMode = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),nModes)
flutterOnsetDispOfMode = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),nModes)
flutterOffsetSpeedsOfMode = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),nModes)
flutterOffsetFreqsOfMode = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),nModes)
flutterOffsetDispOfMode = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),nModes)

# Sweep configurations
for (c,hasTipCorrection,tipLossType) in zip(1:length(tipLossTypeConfig),hasTipCorrectionConfig,tipLossTypeConfig)
    # Sweep angle of sweep
    for (i,Λ) in enumerate(ΛRange)
        # Sweep root pitch angle
        for (j,θ) in enumerate(θRange)
            # Sweep airspeed
            for (k,U) in enumerate(URange)
                # Display progress
                println("Solving for configuration $c, Λ = $(round(Int,Λ*180/π)) deg, θ = $(θ*180/π) deg, U = $U m/s")
                # Model
                model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,Λ=Λ,θ=θ,airspeed=U,g=g,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepStructuralCorrections=sweepStructuralCorrections)
                # Create and solve problem
                problem = create_EigenProblem(model=model,nModes=nModes,systemSolver=NR,frequencyFilterLimits=[2π,Inf])
                solve!(problem)
                # Frequencies, dampings and eigenvectors
                untrackedFreqs[c,i,j,k] = problem.frequenciesOscillatory
                untrackedDamps[c,i,j,k] = round_off!(problem.dampingsOscillatory,1e-8)
                untrackedEigenvectors[c,i,j,k] = problem.eigenvectorsOscillatoryCplx
                # Get OOP displacement
                tipOOP[c,i,j,k] = problem.nodalStatesOverσ[end][nElem].u_n2_b[3]
            end
            # Apply mode tracking
            freqs[c,i,j,:],damps[c,i,j,:],_ = mode_tracking(URange,untrackedFreqs[c,i,j,:],untrackedDamps[c,i,j,:],untrackedEigenvectors[c,i,j,:])
            # Separate frequencies and damping ratios by mode
            for mode in 1:nModes
                modeFrequencies[c,i,j,mode] = [freqs[c,i,j,k][mode] for k in eachindex(URange)]
                modeDampings[c,i,j,mode] = [damps[c,i,j,k][mode] for k in eachindex(URange)]
                modeDampingRatios[c,i,j,mode] = modeDampings[c,i,j,mode]./modeFrequencies[c,i,j,mode]
            end
            # Loop over modes: compute flutter onset and offset speeds, respective frequencies and OOP displacements
            for mode in 1:nModes
                # Find flutter onset indices
                onsetIndices = findall((modeDampings[c,i,j,mode][2:end] .> 0) .& (modeDampings[c,i,j,mode][1:end-1] .< 0)) .+ 1 
                nIndOn = length(onsetIndices)
                # Loop flutter onset indices
                flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetDisp = Vector{Float64}(undef,nIndOn),Vector{Float64}(undef,nIndOn),Vector{Float64}(undef,nIndOn)
                for (n,k) in enumerate(onsetIndices)
                    flutterOnsetSpeeds[n] = LinearInterpolations.interpolate(modeDampings[c,i,j,mode][k-1:k],URange[k-1:k],0)
                    flutterOnsetFreqs[n] = LinearInterpolations.interpolate(modeDampings[c,i,j,mode][k-1:k],modeFrequencies[c,i,j,mode][k-1:k],0)
                    flutterOnsetDisp[n] = LinearInterpolations.interpolate(modeDampings[c,i,j,mode][k-1:k],tipOOP[c,i,j,k-1:k]/L*100,0)
                end
                if nIndOn == 0
                    flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetDisp = [NaN],[NaN],[NaN]
                end
                # Set flutter onset variables for current mode
                flutterOnsetSpeedsOfMode[c,i,j,mode] = flutterOnsetSpeeds
                flutterOnsetFreqsOfMode[c,i,j,mode] = flutterOnsetFreqs
                flutterOnsetDispOfMode[c,i,j,mode] = flutterOnsetDisp
                # Find flutter offset indices
                offsetIndices = findall((modeDampings[c,i,j,mode][2:end] .< 0) .& (modeDampings[c,i,j,mode][1:end-1] .> 0)) .+ 1
                nIndOff = length(offsetIndices)
                # Find flutter offset variables
                flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetDisp = Vector{Float64}(undef,nIndOff),Vector{Float64}(undef,nIndOff),Vector{Float64}(undef,nIndOff)
                # Loop flutter offset indices
                for (n,k) in enumerate(offsetIndices)
                    flutterOffsetSpeeds[n] = LinearInterpolations.interpolate(-modeDampings[c,i,j,mode][k-1:k],URange[k-1:k],0)
                    flutterOffsetFreqs[n] = LinearInterpolations.interpolate(-modeDampings[c,i,j,mode][k-1:k],modeFrequencies[c,i,j,mode][k-1:k],0)
                    flutterOffsetDisp[n] = LinearInterpolations.interpolate(-modeDampings[c,i,j,mode][k-1:k],tipOOP[c,i,j,k-1:k]/L*100,0)
                end
                if nIndOff == 0
                    flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetDisp = [NaN],[NaN],[NaN]
                end
                # Set flutter offset variables for current mode
                flutterOffsetSpeedsOfMode[c,i,j,mode] = flutterOffsetSpeeds
                flutterOffsetFreqsOfMode[c,i,j,mode] = flutterOffsetFreqs
                flutterOffsetDispOfMode[c,i,j,mode] = flutterOffsetDisp
            end
        end
    end
end

# Sharpy and Feniax data (from AePW4 meetings)
flutterBoundary_dispVsU_Lambda0_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda0_Sharpy.txt")
flutterBoundary_dispVsU_Lambda10_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda10_Sharpy.txt")
flutterBoundary_dispVsU_Lambda0_Feniax = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda0_Feniax.txt")
flutterBoundary_dispVsU_Lambda10_Feniax = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda10_Feniax.txt")
flutterBoundary_dispVsU_Lambda20_Feniax = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda20_Feniax.txt")
flutterBoundary_dispVsU_Lambda30_Feniax = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda30_Feniax.txt")

# Test data by Avin et al.
flutterOnsetVelUp = [49; 43; 38]
flutterOnsetVelDown = [55; 48]
flutterOnsetDispUp = [23; 24.5; 26]
flutterOnsetDispDown = [28.5; 31.5]

println("Finished sweptPazyFlutterSweepRangePitchRange.jl")