using AeroBeams

# Mode tracking option
modeTracking = true

# Stiffness factor
λ = 1e0

# Aerodynamic solver
aeroSolver = Indicial()

# Tip loss option
hasTipCorrection = true

# Option to update airfoil parameters (apply Prandtl-Glauert lift-slope correction)
updateAirfoilParameters = true

# Set NR system solver 
relaxFactor = 0.5
displayStatus = false
maxiter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxiter,displayStatus=displayStatus)

# Set number of vibration modes
nModes = 8

# Set airspeed range and initialize outputs
URange = vcat(30:2:160)
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))

# Attachment springs
μ = 1e-2
ku = μ*[1; 1; 1]
kp = ku
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)
spring2 = create_Spring(elementsIDs=[3],nodesSides=[2],ku=ku,kp=kp)

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Model for trim problem
    BWBtrim = create_BWB(aeroSolver=aeroSolver,stiffnessFactor=λ,hasTipCorrection=hasTipCorrection,updateAirfoilParameters=updateAirfoilParameters,airspeed=U,δElevIsTrimVariable=true,thrustIsTrimVariable=true)
    # Add springs
    add_springs_to_beam!(beam=BWBtrim.beams[2],springs=[spring1])
    add_springs_to_beam!(beam=BWBtrim.beams[3],springs=[spring2])
    # Update model
    BWBtrim.skipValidationMotionBasisA = true
    update_model!(BWBtrim)
    # Set initial guess solution as previous known solution
    x0Trim = i == 1 ? zeros(0) : trimProblem.x
    # Create and trim problem
    global trimProblem = create_TrimProblem(model=BWBtrim,systemSolver=NR,x0=x0Trim)
    solve!(trimProblem)
    # Extract trim variables
    trimAoA[i] = trimProblem.aeroVariablesOverσ[end][BWBtrim.beams[3].elementRange[1]].flowAnglesAndRates.αₑ
    trimThrust[i] = trimProblem.x[end-1]*BWBtrim.forceScaling 
    trimδ[i] = trimProblem.x[end]
    println("Trim AoA = $(trimAoA[i]*180/π), trim thrust = $(trimThrust[i]), trim δ = $(trimδ[i]*180/π)")
    # Model for eigen problem
    BWBeigen = create_BWB(aeroSolver=aeroSolver,stiffnessFactor=λ,hasTipCorrection=hasTipCorrection,updateAirfoilParameters=updateAirfoilParameters,airspeed=U,δElev=trimδ[i],thrust=trimThrust[i])
    # Add springs
    add_springs_to_beam!(beam=BWBeigen.beams[2],springs=[spring1])
    add_springs_to_beam!(beam=BWBeigen.beams[3],springs=[spring2])
    # Update model
    BWBeigen.skipValidationMotionBasisA = true
    update_model!(BWBeigen)
    # Create and solve eigen problem
    global eigenProblem = create_EigenProblem(model=BWBeigen,nModes=nModes,frequencyFilterLimits=[1e0,Inf64],jacobian=trimProblem.jacobian[1:end,1:end-trimProblem.model.nTrimVariables],inertia=trimProblem.inertia)
    solve_eigen!(eigenProblem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = eigenProblem.frequenciesOscillatory
    untrackedDamps[i] = round_off!(eigenProblem.dampingsOscillatory,1e-12)
    untrackedEigenvectors[i] = eigenProblem.eigenvectorsOscillatoryCplx
end

# Frequencies and dampings after mode tracking
if modeTracking
    freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
else
    freqs,damps = untrackedFreqs,untrackedDamps
end

# Separate frequencies and damping ratios by mode
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeFrequencies =  Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
end

println("Finished BWBflutter.jl")