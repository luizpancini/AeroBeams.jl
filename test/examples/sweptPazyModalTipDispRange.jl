using AeroBeams

# Sweep angle [rad] range
ΛRange = [0,20,30]*π/180

# Tip load
F3TipRange = vcat(0:1:35)

# Flag for ad hoc sectional stiffness corrections with sweep angle
sweepStructuralCorrections = true

# Number of modes
nModes = 3

# Gravity
g = 0

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Initialize outputs
problem = Array{EigenProblem}(undef,length(ΛRange),length(F3TipRange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(ΛRange),length(F3TipRange))
untrackedDamps = Array{Vector{Float64}}(undef,length(ΛRange),length(F3TipRange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(ΛRange),length(F3TipRange))
freqs = Array{Vector{Float64}}(undef,length(ΛRange),length(F3TipRange))
damps = Array{Vector{Float64}}(undef,length(ΛRange),length(F3TipRange))
modeFrequencies = Array{Vector{Float64}}(undef,length(ΛRange),nModes)
tipOOP = Array{Float64}(undef,length(ΛRange),length(F3TipRange))

# System solver
absTol = 1e-7
NR = create_NewtonRaphson(displayStatus=false,absoluteTolerance=absTol)

# Loop sweep angle
for (i,Λ) in enumerate(ΛRange)
    for (j,F3) in enumerate(F3TipRange)
        display("Solving for Λ = $(round(Λ*180/π)) deg, tip force = $F3 N")
        # Dummy beam for the BC
        dummyBeam = create_Beam(length=1,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1)])
        # Tip force BC
        tipForce = create_BC(name="tipForce",beam=dummyBeam,node=nElem+1,types=["F3b"],values=[F3])
        # Model
        sweptPazyModalTipDispRange,_ = create_Pazy(Λ=Λ,upright=false,g=g,sweepStructuralCorrections=sweepStructuralCorrections,additionalBCs=[tipForce])
        # Create and solve problem
        problem[i,j] = create_EigenProblem(model=sweptPazyModalTipDispRange,nModes=nModes,systemSolver=NR)
        solve!(problem[i,j])
        # Get frequencies in Hz, dampings, eigenvectors, and tip OOP displacement
        untrackedFreqs[i,j] = problem[i,j].frequenciesOscillatory/(2π)
        untrackedDamps[i,j] = round_off!(problem[i,j].dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = problem[i,j].eigenvectorsOscillatoryCplx
        tipOOP[i,j] = problem[i,j].nodalStatesOverσ[end][nElem].u_n2_b[3]
    end
    # Apply mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking_hungarian(F3TipRange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(F3TipRange)]
    end
end

println("Finished sweptPazyModalTipDispRange.jl")