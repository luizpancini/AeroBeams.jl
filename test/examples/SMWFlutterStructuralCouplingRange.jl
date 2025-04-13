using AeroBeams, LinearAlgebra, LinearInterpolations, DelimitedFiles

# Torsion / IP bending coupling factor range
ΨRange = [-0.2, 0, 0.2]

# Tip load range
F3Range = collect(0:1:35)

# Airspeed range
URange = collect(20:0.5:35)

# Aerodynamic solver
aeroSolver = Indicial()

# Altitude 
h = 20e3

# Gravity
g = 0

# Discretization
nElem = 16

# Set number of vibration modes
nModes = 5

# Initialize outputs
freqs = Array{Vector{Float64}}(undef,length(ΨRange),length(F3Range),length(URange))
damps = Array{Vector{Float64}}(undef,length(ΨRange),length(F3Range),length(URange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
flutterSpeed = Array{Float64}(undef,length(ΨRange),length(F3Range))
dampsCurrentΨF3 = Array{Vector{Float64}}(undef,nModes)
flutterSpeedOfMode = Array{Float64}(undef,nModes)

# Sweep structural coupling 
for (i,Ψ) in enumerate(ΨRange)
    # Sweep tip force
    for (j,F3) in enumerate(F3Range)
        # Sweep airspeed
        for (k,U) in enumerate(URange)
            # Display progress
            println("Solving for Ψ=$Ψ, F=$F3 N, U = $U m/s")
            # Create model
            model,_ = create_SMW(aeroSolver=aeroSolver,nElem=nElem,altitude=h,airspeed=U,g=g,Ψ=Ψ,tipF3=F3)
            # Create and solve problem
            problem = create_EigenProblem(model=model,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
            solve!(problem)
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[k] = problem.frequenciesOscillatory
            untrackedDamps[k] = round_off!(problem.dampingsOscillatory,1e-8)
            untrackedEigenvectors[k] = problem.eigenvectorsOscillatoryCplx
        end
        # Frequencies and dampings after mode tracking
        freqs[i,j,:],damps[i,j,:],_ = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
        # Flutter speeds of modes at current structural coupling and tip force
        for mode in 1:nModes
            dampsCurrentΨF3[mode] = [damps[i,j,k][mode] for k in eachindex(URange)]
            indexInstability = findfirst(x->x>0,dampsCurrentΨF3[mode])
            flutterSpeedOfMode[mode] = isnothing(indexInstability) || indexInstability == 1 ? NaN : LinearInterpolations.interpolate(dampsCurrentΨF3[mode][indexInstability-1:indexInstability],URange[indexInstability-1:indexInstability],0)
        end
        # Set flutter speed as the greatest (for compatibility with reference)
        flutterSpeed[i,j] = maximum(filter(!isnan,flutterSpeedOfMode))
    end
end

# Load reference data
flutterSpeedVsTipLoadΨ0 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "SMW", "flutterSpeedVsTipLoadPsi0_0.txt"))
flutterSpeedVsTipLoadΨp02 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "SMW", "flutterSpeedVsTipLoadPsi0_2.txt"))
flutterSpeedVsTipLoadΨm02 = readdlm(joinpath(dirname(@__DIR__), "referenceData", "SMW", "flutterSpeedVsTipLoadPsi-0_2.txt"))

println("Finished SMWFlutterStructuralCouplingRange.jl")