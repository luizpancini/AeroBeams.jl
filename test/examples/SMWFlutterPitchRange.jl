using AeroBeams, LinearAlgebra, LinearInterpolations, DelimitedFiles

# Aerodynamic solver and derivatives method
aeroSolver = Indicial()
derivationMethod = AD()

# Altitude
h = 20e3

# Gravity
g = 9.80665

# Discretization
nElem = 16

# Set system solver options (limit initial load factor)
NR1 = create_NewtonRaphson(initialLoadFactor=0.5,maximumLoadFactorStep=0.25)
NR2 = create_NewtonRaphson(initialLoadFactor=0.5,maximumLoadFactorStep=0.5)

# Set number of vibration modes
nModes = 5

# Set root angle and airspeed ranges, and initialize outputs
θRange = collect(0:0.1:5.0)
URange = collect(0:0.2:35)
untrackedFreqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
tip_u3 = Array{Float64}(undef,length(θRange),length(URange))
flutterOnsetSpeed = Array{Float64}(undef,length(θRange),nModes)
flutterOffsetSpeed = Array{Float64}(undef,length(θRange),nModes)
flutterOnsetFreq = Array{Float64}(undef,length(θRange),nModes)
flutterOffsetFreq = Array{Float64}(undef,length(θRange),nModes)
flutterOnsetTipDisp = Array{Float64}(undef,length(θRange),nModes)
flutterOffsetTipDisp = Array{Float64}(undef,length(θRange),nModes)

# Sweep root angle
for (i,θ) in enumerate(θRange)
    # Update model
    SMWFlutterPitchRange,_ = create_SMW(aeroSolver=aeroSolver,derivationMethod=derivationMethod,θ=θ*π/180,nElem=nElem,altitude=h,g=g)
    # Update system solver
    NR = θ <= 1.6 ? NR1 : NR2
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $θ deg, U = $U m/s")
        # Update velocity of basis A (and update model)
        set_motion_basis_A!(model=SMWFlutterPitchRange,v_A=[0;U;0])
        # Create and solve problem
        problem = create_EigenProblem(model=SMWFlutterPitchRange,systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
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
    # Flutter speeds, frequencies and tip displacements of modes at current root angle
    dampsCurrentθ = Array{Vector{Float64}}(undef,nModes)
    freqsCurrentθ = Array{Vector{Float64}}(undef,nModes)
    for mode in 1:nModes
        # Dampings and frequencies of mode
        dampsCurrentθ[mode] = [damps[i,j][mode] for j in eachindex(URange)]
        freqsCurrentθ[mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        # Flutter onset 
        iOnset = findfirst(x->x>0,dampsCurrentθ[mode])
        if isnothing(iOnset) || iOnset == 1
            flutterOnsetSpeed[i,mode] = NaN
            flutterOnsetFreq[i,mode] = NaN
            flutterOnsetTipDisp[i,mode] = NaN
            flutterOffsetSpeed[i,mode] = NaN
            flutterOffsetFreq[i,mode] = NaN
            flutterOffsetTipDisp[i,mode] = NaN
            continue
        end
        flutterOnsetSpeed[i,mode] = interpolate(dampsCurrentθ[mode][iOnset-1:iOnset],URange[iOnset-1:iOnset],0)
        flutterOnsetFreq[i,mode] = interpolate(dampsCurrentθ[mode][iOnset-1:iOnset],freqsCurrentθ[mode][iOnset-1:iOnset],0)
        flutterOnsetTipDisp[i,mode] = interpolate(dampsCurrentθ[mode][iOnset-1:iOnset],tip_u3[i,iOnset-1:iOnset],0)
        # Flutter offset
        iOffset = findnext(x->x<0,dampsCurrentθ[mode],iOnset)
        if isnothing(iOffset)
            flutterOffsetSpeed[i,mode] = NaN
            flutterOffsetFreq[i,mode] = NaN
            flutterOffsetTipDisp[i,mode] = NaN
            continue
        end
        flutterOffsetSpeed[i,mode] = interpolate(-dampsCurrentθ[mode][iOffset-1:iOffset],URange[iOffset-1:iOffset],0)
        flutterOffsetFreq[i,mode] = interpolate(-dampsCurrentθ[mode][iOffset-1:iOffset],freqsCurrentθ[mode][iOffset-1:iOffset],0)
        flutterOffsetTipDisp[i,mode] = interpolate(-dampsCurrentθ[mode][iOffset-1:iOffset],tip_u3[i,iOffset-1:iOffset],0)
    end
end

# Load reference data
flutterSpeedRef = readdlm("test/referenceData/SMW/flutterSpeedVsRootAoA.txt")
flutterFreqRef = readdlm("test/referenceData/SMW/flutterFreqVsRootAoA.txt")
flutterTipDispRef = readdlm("test/referenceData/SMW/flutterTipDispVsRootAoA.txt")

speedVsDispRootAoA0 = readdlm("test/referenceData/SMW/speedVsDispRootAoA0_0.txt")
speedVsDispRootAoA05 = readdlm("test/referenceData/SMW/speedVsDispRootAoA0_5.txt")
speedVsDispRootAoA1 = readdlm("test/referenceData/SMW/speedVsDispRootAoA1_0.txt")
speedVsDispRootAoA2 = readdlm("test/referenceData/SMW/speedVsDispRootAoA2_0.txt")

freqVsSpeedRootAoA2 = readdlm("test/referenceData/SMW/V-f-RootAoA2.txt")
dampVsSpeedRootAoA2 = readdlm("test/referenceData/SMW/V-g-RootAoA2.txt")
freqVsSpeedRootAoA3 = readdlm("test/referenceData/SMW/V-f-RootAoA3.txt")
dampVsSpeedRootAoA3 = readdlm("test/referenceData/SMW/V-g-RootAoA3.txt")
freqVsSpeedRootAoA5 = readdlm("test/referenceData/SMW/V-f-RootAoA5.txt")
dampVsSpeedRootAoA5 = readdlm("test/referenceData/SMW/V-g-RootAoA5.txt")

println("Finished SMWFlutterPitchRange.jl")