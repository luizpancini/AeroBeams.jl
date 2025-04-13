using AeroBeams

# Sweep angle [rad]
ΛRange = [0,20,30]*π/180

# Tip force range
FRange = collect(0:.25:30)

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Flag for upright position
upright = false

# Gravity
g = 0

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Number of modes
nModes = 3

# System solver
σ0 = 0.5
σstep = 0.5
maxIter = 50
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep,maximumIterations=maxIter)

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(ΛRange),length(FRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(ΛRange),length(FRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(ΛRange),length(FRange))
freqs = Array{Vector{Float64}}(undef,length(ΛRange),length(FRange))
modeFrequencies = Array{Vector{Float64}}(undef,length(ΛRange),nModes)
tipOOP = Array{Float64}(undef,length(ΛRange),length(FRange))

# Loop angle of sweep
for (i,Λ) in enumerate(ΛRange)
    # Loop tip force
    for (j,F) in enumerate(FRange)
        # Display progress
        println("Solving Λ = $(round(Λ*180/π)) deg, F = $F N")
        # Set tip force BC on dummy beam
        dummyBeam = create_Beam(length=1,nElements=nElem,S=[isotropic_stiffness_matrix()])
        tipForce = create_BC(name="tipForce",beam=dummyBeam,node=nElem+1,types=["F3b"],values=[F])
        # Model
        model,_ = create_Pazy(upright=upright,Λ=Λ,g=g,sweepStructuralCorrections=sweepStructuralCorrections,additionalBCs=[tipForce])
        # Create and solve problem
        problem = create_EigenProblem(model=model,nModes=nModes,systemSolver=NR)
        solve!(problem)
        # Outputs
        untrackedFreqs[i,j] = problem.frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(problem.dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = problem.eigenvectorsOscillatoryCplx
        tipOOP[i,j] = problem.nodalStatesOverσ[end][nElem].u_n2_b[3]
    end
    # Apply mode tracking
    freqs[i,:],_ = mode_tracking(FRange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(FRange)]
    end
end

println("Finished sweptPazyEigenStructuralSweepCorrectionDispRange.jl")