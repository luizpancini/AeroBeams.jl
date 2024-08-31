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
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep,displayStatus=false)

# Set number of vibration modes
nModes = 5

# Set precurvature, root angle, and airspeed ranges, and initialize outputs
kRange = collect(0:0.018:0.018)
θRange = collect(0:0.2:5)
URange = collect(0:0.5:33)
untrackedFreqs = Array{Vector{Float64}}(undef,length(kRange),length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(kRange),length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(kRange),length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(kRange),length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(kRange),length(θRange),length(URange))
tip_u3 = Array{Float64}(undef,length(kRange),length(θRange),length(URange))
flutterOnsetSpeed = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]
flutterOffsetSpeed = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]
flutterOnsetFreq = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]
flutterOffsetFreq = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]
flutterOnsetTipDisp = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]
flutterOffsetTipDisp = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]

# Sweep wing precurvature
for (ki,k) in enumerate(kRange)
    # Sweep root angle
    for (i,θ) in enumerate(θRange)
        # Update model
        SMWFlutterPrecurvatureRange2,_ = create_SMW(aeroSolver=aeroSolver,derivationMethod=derivationMethod,θ=θ*π/180,k2=k,nElem=nElem,altitude=h,g=g)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for k=$k, θ = $θ deg, U = $U m/s")
            # Update velocity of basis A (and update model)
            set_motion_basis_A!(model=SMWFlutterPrecurvatureRange2,v_A=[0;U;0])
            # Create and solve problem
            problem = create_EigenProblem(model=SMWFlutterPrecurvatureRange2,systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
            solve!(problem)
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[ki,i,j] = problem.frequenciesOscillatory
            untrackedDamps[ki,i,j] = round_off!(problem.dampingsOscillatory,1e-8)
            untrackedEigenvectors[ki,i,j] = problem.eigenvectorsOscillatoryCplx
            # Tip OOP displacement
            tip_u3[ki,i,j] = problem.nodalStatesOverσ[end][nElem].u_n2[3]
        end
        # Frequencies and dampings after mode tracking
        freqs[ki,i,:],damps[ki,i,:],_ = mode_tracking(URange,untrackedFreqs[ki,i,:],untrackedDamps[ki,i,:],untrackedEigenvectors[ki,i,:])
        # Flutter speeds, frequencies and tip displacements of modes at current combination of curvature and root pitch
        dampsCurrentkθ = Array{Vector{Float64}}(undef,nModes)
        freqsCurrentkθ = Array{Vector{Float64}}(undef,nModes)
        for mode in 1:nModes
            dampsCurrentkθ[mode] = [damps[ki,i,j][mode] for j in eachindex(URange)]
            freqsCurrentkθ[mode] = [freqs[ki,i,j][mode] for j in eachindex(URange)]
            # Flutter onset
            iOnset = 1 .+ findall(i -> dampsCurrentkθ[mode][i] < 0 && dampsCurrentkθ[mode][i+1] > 0, 1:length(dampsCurrentkθ[mode])-1)
            if isempty(iOnset) || isempty(filter!(x->x!=1,iOnset))
                push!(flutterOnsetSpeed[ki][i][mode],NaN)
                push!(flutterOnsetFreq[ki][i][mode],NaN)
                push!(flutterOnsetTipDisp[ki][i][mode],NaN)
                push!(flutterOffsetSpeed[ki][i][mode],NaN)
                push!(flutterOffsetFreq[ki][i][mode],NaN)
                push!(flutterOffsetTipDisp[ki][i][mode],NaN)
                continue
            end
            for j in iOnset
                push!(flutterOnsetSpeed[ki][i][mode],interpolate(dampsCurrentkθ[mode][j-1:j],URange[j-1:j],0))
                push!(flutterOnsetFreq[ki][i][mode],interpolate(dampsCurrentkθ[mode][j-1:j],freqsCurrentkθ[mode][j-1:j],0))
                push!(flutterOnsetTipDisp[ki][i][mode],interpolate(dampsCurrentkθ[mode][j-1:j],tip_u3[ki,i,j-1:j],0))
            end
            # Flutter offset
            iOffset = 1 .+ findall(i -> dampsCurrentkθ[mode][i] > 0 && dampsCurrentkθ[mode][i+1] < 0, 1:length(dampsCurrentkθ[mode])-1)
            if isempty(iOffset)
                push!(flutterOffsetSpeed[ki][i][mode],NaN)
                push!(flutterOffsetFreq[ki][i][mode],NaN)
                push!(flutterOffsetTipDisp[ki][i][mode],NaN)
                continue
            end
            for j in iOffset
                push!(flutterOffsetSpeed[ki][i][mode],interpolate(-dampsCurrentkθ[mode][j-1:j],URange[j-1:j],0))
                push!(flutterOffsetFreq[ki][i][mode],interpolate(-dampsCurrentkθ[mode][j-1:j],freqsCurrentkθ[mode][j-1:j],0))
                push!(flutterOffsetTipDisp[ki][i][mode],interpolate(-dampsCurrentkθ[mode][j-1:j],tip_u3[ki,i,j-1:j],0))
            end
        end
    end
end

println("Finished SMWFlutterPrecurvatureRange2.jl")