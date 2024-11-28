abstract type Problem end
abstract type SystemSolver end


"""
    InitialVelocitiesUpdateOptions composite type

Defines variables for the update of the initial velocities states

# Fields
- `Δt::Real` = time step
- `maxIter::Int64` = maximum number of iterations
- `relaxFactor::Float64` = relaxation factor
- `tol::Float64` = convergence tolerance
- `displayProgress::Bool` = flag to display progress
"""
@with_kw mutable struct InitialVelocitiesUpdateOptions

    # Fields
    Δt::Real = 0
    maxIter::Int64 = 2
    relaxFactor::Real = 0.5
    tol::Float64 = 1e-6
    displayProgress::Bool = false

    # Constructor
    function InitialVelocitiesUpdateOptions(Δt::Real,maxIter::Int64,relaxFactor::Real,tol::Float64,displayProgress::Bool)

        # Check relaxation factor
        @assert (0 < relaxFactor <= 1)

        return new(Δt,maxIter,relaxFactor,tol,displayProgress)
    end
end
export InitialVelocitiesUpdateOptions


#
# mutable struct ModeShape{T<:Union{Float64,ComplexF64}}

    # ModeShape composite type

#
mutable struct ModeShape{T<:Union{Float64,ComplexF64}}
    
    # Fields
    mode::Int64
    frequency::Float64
    damping::Float64
    elementalStates::Vector{ElementalStates{T}}
    complementaryElementalStates::Vector{ComplementaryElementalStates{T}}
    nodalStates::Vector{NodalStates{T}}

    # Default constructor 
    function ModeShape{T}() where T<:Union{Float64,ComplexF64}

        mode = 0
        frequency = 0.0
        damping = 0.0
        elementalStates = Vector{ElementalStates{T}}()
        complementaryElementalStates = Vector{ComplementaryElementalStates{T}}()
        nodalStates = Vector{NodalStates{T}}()

        return new{T}(mode,frequency,damping,elementalStates,complementaryElementalStates,nodalStates)

    end

    # Alternate constructor 
    function ModeShape{T}(mode::Int64,frequency::Float64,damping::Float64,nElementsTotal::Int64) where T<:Union{Float64,ComplexF64}

        elementalStates = [ElementalStates{T}() for _ in 1:nElementsTotal]
        complementaryElementalStates = [ComplementaryElementalStates{T}() for _ in 1:nElementsTotal]
        nodalStates = [NodalStates{T}() for _ in 1:nElementsTotal]

        return new{T}(mode,frequency,damping,elementalStates,complementaryElementalStates,nodalStates)

    end
end


#

    # SteadyProblem composite type

#
@with_kw mutable struct SteadyProblem <: Problem

    # Primary (inputs to problem creation)
    # ------------------------------------
    # Model
    model::Model
    # System solver
    systemSolver::SystemSolver
    # TF to get linear solution
    getLinearSolution::Bool
    # Secondary (outputs from problem creation)
    # -----------------------------------------
    # States, residual, Jacobian and inertia arrays
    x::Vector{Float64} = zeros(0)
    Δx::Vector{Float64} = zeros(0)
    residual::Vector{Float64} = zeros(0)
    jacobian::SparseMatrixCSC{Float64,Int64} = spzeros(0,0)
    inertia::SparseMatrixCSC{Float64,Int64} = spzeros(0,0)
    # Dummy time
    timeNow::Float64 = 0.0
    # Load factor
    σ::Float64 = 1.0
    # TF to compute only the external forces array at the current nonlinear step (used only for arclength system solver)
    getExternalForcesArray::Bool = false
    # TF for initial states being input
    initialStatesInput::Bool = false
    # TF to compute only residual array and skip Jacobian update at the current nonlinear step (used only in line search's Newton step)
    skipJacobianUpdate::Bool = false
    # Array of partial load steps and respective solutions
    savedσ::Vector{Float64} = Vector{Float64}()
    xOverσ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    elementalStatesOverσ::Vector{Vector{ElementalStates{Float64}}} = Vector{Vector{ElementalStates{Float64}}}()
    nodalStatesOverσ::Vector{Vector{NodalStates{Float64}}} = Vector{Vector{NodalStates{Float64}}}()
    compElementalStatesOverσ::Vector{Vector{ComplementaryElementalStates{Float64}}} = Vector{Vector{ComplementaryElementalStates{Float64}}}()
    aeroVariablesOverσ::Vector{Vector{Union{Nothing,AeroVariables}}} = Vector{Vector{Union{Nothing,AeroVariables}}}()
    # Maximum absolute value of aerodynamic load coefficients over time
    maxAeroForce::Real = 0
    maxAeroMoment::Real = 0

end
export SteadyProblem


"""
    create_SteadyProblem(; kwargs...)

Steady problem constructor

# Keyword arguments
- `model::Model` = model
- `systemSolver::systemSolver` = nonlinear system solver
- `getLinearSolution::Bool` = flag to solve for linear structural solution
- `x0::Vector{Float64}` = initial states
"""
function create_SteadyProblem(; model::Model,systemSolver::SystemSolver=create_NewtonRaphson(),getLinearSolution::Bool=false,x0::Vector{Float64}=zeros(0))

    # Initialize problem
    problem = SteadyProblem(model=model,systemSolver=systemSolver,getLinearSolution=getLinearSolution)

    # Update initial load factor
    problem.σ = systemSolver.initialLoadFactor

    # Set initial elemental and nodal states
    if !isempty(x0)
        @assert length(x0) == model.systemOrder
        problem.initialStatesInput = true
        problem.x = x0*problem.σ
    else
        set_initial_states!(problem)
    end

    # Initialize system arrays with correct size
    initialize_system_arrays!(problem)

    return problem

end
export create_SteadyProblem


#

    # TrimProblem composite type

#
@with_kw mutable struct TrimProblem <: Problem

    
    # Primary (inputs to problem creation)
    # ------------------------------------
    # Model
    model::Model
    # System solver
    systemSolver::SystemSolver
    # TF to get linear solution
    getLinearSolution::Bool
    # TF to get the inertia matrix upon converged solution
    getInertiaMatrix::Bool
    # Secondary (outputs from problem creation)
    # -----------------------------------------
    # States, residual, Jacobian and inertia arrays
    x::Vector{Float64} = zeros(0)
    Δx::Vector{Float64} = zeros(0)
    residual::Vector{Float64} = zeros(0)
    jacobian::SparseMatrixCSC{Float64,Int64} = spzeros(0,0)
    inertia::SparseMatrixCSC{Float64,Int64} = spzeros(0,0)
    # Dummy time
    timeNow::Float64 = 0.0
    # Load factor
    σ::Float64 = 1.0
    # TF to compute only the external forces array at the current nonlinear step (used only for arclength system solver)
    getExternalForcesArray::Bool = false
    # TF for initial states being input
    initialStatesInput::Bool = false
    # TF to compute only residual array and skip Jacobian update at the current nonlinear step (used only in line search's Newton step)
    skipJacobianUpdate::Bool = false
    # Array of partial load steps and respective solutions
    savedσ::Vector{Float64} = Vector{Float64}()
    xOverσ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    elementalStatesOverσ::Vector{Vector{ElementalStates{Float64}}} = Vector{Vector{ElementalStates{Float64}}}()
    nodalStatesOverσ::Vector{Vector{NodalStates{Float64}}} = Vector{Vector{NodalStates{Float64}}}()
    compElementalStatesOverσ::Vector{Vector{ComplementaryElementalStates{Float64}}} = Vector{Vector{ComplementaryElementalStates{Float64}}}()
    aeroVariablesOverσ::Vector{Vector{Union{Nothing,AeroVariables}}} = Vector{Vector{Union{Nothing,AeroVariables}}}()
    # Maximum absolute value of aerodynamic load coefficients over time
    maxAeroForce::Real = 0
    maxAeroMoment::Real = 0

end
export TrimProblem


"""
    create_TrimProblem(; kwargs...)

Trim problem constructor

# Keyword arguments
- `model::Model` = model
- `systemSolver::systemSolver` = nonlinear system solver
- `getLinearSolution::Bool` = flag to solve for linear structural solution
- `getInertiaMatrix::Bool` = flag to compute inertia matrix
- `x0::Vector{Float64}` = initial states
"""
function create_TrimProblem(; model::Model,systemSolver::SystemSolver=create_NewtonRaphson(),getLinearSolution::Bool=false,getInertiaMatrix::Bool=true,x0::Vector{Float64}=zeros(0))

    # Initialize problem
    problem = TrimProblem(model=model,systemSolver=systemSolver,getLinearSolution=getLinearSolution,getInertiaMatrix=getInertiaMatrix)

    # Update initial load factor
    problem.σ = systemSolver.initialLoadFactor

    # Set initial elemental and nodal states
    if !isempty(x0)
        @assert length(x0) == model.systemOrder+model.nTrimVariables
        problem.initialStatesInput = true
        problem.x = x0*problem.σ
    else
        set_initial_states!(problem)
    end

    # Initialize system arrays with correct size
    initialize_system_arrays!(problem)

    return problem

end
export create_TrimProblem


#

    # EigenProblem composite type

#
@with_kw mutable struct EigenProblem <: Problem

    # Primary (inputs to problem creation)
    # ------------------------------------
    # Model
    model::Model
    # System solver
    systemSolver::SystemSolver
    # TF to get linear solution
    getLinearSolution::Bool
    # Real of desired oscillatory modes
    nModes::Int64
    # Frequency filter limits
    frequencyFilterLimits::Vector{Float64}
    # TF to normalize mode shapes
    normalizeModeShapes::Bool
    # Secondary (outputs from problem creation)
    # -----------------------------------------
    # States, residual, Jacobian and inertia arrays
    x::Vector{Float64} = zeros(0)
    Δx::Vector{Float64} = zeros(0)
    residual::Vector{Float64} = zeros(0)
    jacobian::SparseMatrixCSC{Float64,Int64} = spzeros(0,0)
    inertia::SparseMatrixCSC{Float64,Int64} = spzeros(0,0)
    # Dummy time
    timeNow::Float64 = 0.0
    # Load factor
    σ::Float64 = 1.0
    # TF to compute only the external forces array at the current nonlinear step (used only for arclength system solver)
    getExternalForcesArray::Bool = false
    # TF for initial states being input
    initialStatesInput::Bool = false
    # TF to compute only residual array and skip Jacobian update at the current nonlinear step (used only in line search's Newton step)
    skipJacobianUpdate::Bool = false
    # Array of partial load steps and respective solutions
    savedσ::Vector{Float64} = Vector{Float64}()
    xOverσ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    elementalStatesOverσ::Vector{Vector{ElementalStates{Float64}}} = Vector{Vector{ElementalStates{Float64}}}()
    nodalStatesOverσ::Vector{Vector{NodalStates{Float64}}} = Vector{Vector{NodalStates{Float64}}}()
    compElementalStatesOverσ::Vector{Vector{ComplementaryElementalStates{Float64}}} = Vector{Vector{ComplementaryElementalStates{Float64}}}()
    aeroVariablesOverσ::Vector{Vector{Union{Nothing,AeroVariables}}} = Vector{Vector{Union{Nothing,AeroVariables}}}()
    # Frequencies, dampings and eigenvectors
    frequencies::Vector{Float64} = Vector{Float64}()
    dampings::Vector{Float64} = Vector{Float64}()
    eigenvectors::Matrix{ComplexF64} = zeros(ComplexF64, 0, 0)
    frequenciesFiltered::Vector{Float64} = Vector{Float64}()
    dampingsFiltered::Vector{Float64} = Vector{Float64}()
    eigenvectorsFiltered::Matrix{ComplexF64} = zeros(ComplexF64, 0, 0)
    dampingsNonOscillatory::Vector{Float64} = Vector{Float64}()
    frequenciesOscillatory::Vector{Float64} = Vector{Float64}()
    dampingsOscillatory::Vector{Float64} = Vector{Float64}()
    eigenvectorsNonOscillatory::Matrix{ComplexF64} = zeros(ComplexF64, 0, 0)
    eigenvectorsOscillatoryCplx::Matrix{ComplexF64} = zeros(ComplexF64, 0, 0)
    eigenvectorsOscillatoryAbs::Matrix{Float64} = zeros(0, 0)
    # Mode shapes (complex-valued and absolute-valued)
    modeShapesCplx::Vector{ModeShape{ComplexF64}} = Vector{ModeShape{ComplexF64}}()
    modeShapesAbs::Vector{ModeShape{Float64}} = Vector{ModeShape{Float64}}()
    modeShapesAbsNonOsc::Vector{ModeShape{Float64}} = Vector{ModeShape{Float64}}()
    # Maximum absolute value of aerodynamic load coefficients over time
    maxAeroForce::Real = 0
    maxAeroMoment::Real = 0

end
export EigenProblem


"""
    create_EigenProblem(; kwargs...)

Eigen problem constructor

# Keyword arguments
- `model::Model` = model
- `systemSolver::systemSolver` = nonlinear system solver
- `getLinearSolution::Bool` = flag to solve for linear structural solution
- `nModes::Int64=Inf64` = number of modes to be computed
- `frequencyFilterLimits::Vector{Float64}` = limits of the frequency filter
- `normalizeModeShapes::Bool` = flag to normalize mode shapes
- `x0::Vector{Float64}` = initial states
- `jacobian::SparseMatrixCSC{Float64,Int64}` = Jacobian matrix
- `inertia::SparseMatrixCSC{Float64,Int64}` = inertia matrix
"""
function create_EigenProblem(; model::Model,systemSolver::SystemSolver=create_NewtonRaphson(),getLinearSolution::Bool=false,nModes::Int64=Inf64,frequencyFilterLimits::Vector{Float64}=[0,Inf64],normalizeModeShapes::Bool=true,x0::Vector{Float64}=zeros(0),jacobian::SparseMatrixCSC{Float64,Int64}=spzeros(0,0),inertia::SparseMatrixCSC{Float64,Int64}=spzeros(0,0))

    # Initialize problem
    problem = EigenProblem(model=model,systemSolver=systemSolver,getLinearSolution=getLinearSolution,nModes=nModes,frequencyFilterLimits=frequencyFilterLimits,normalizeModeShapes=normalizeModeShapes)

    # Update initial load factor
    problem.σ = systemSolver.initialLoadFactor

    # Set initial elemental and nodal states
    if !isempty(x0)
        @assert length(x0) == model.systemOrder
        problem.initialStatesInput = true
        problem.x = x0*problem.σ
    else
        set_initial_states!(problem)
    end

    # Initialize system arrays
    @assert size(jacobian) == size(inertia) "jacobian and inertia matrices must have the same size"
    if !isempty(jacobian)
        @assert size(jacobian) == (model.systemOrder,model.systemOrder) "size of the input jacobian does not correspond to the number of states of the model"
        problem.jacobian = jacobian
        problem.inertia = inertia
    else
        initialize_system_arrays!(problem)
    end

    return problem

end
export create_EigenProblem


#

    # DynamicProblem composite type

#
@with_kw mutable struct DynamicProblem <: Problem

    # Primary (inputs to problem creation)
    # ------------------------------------
    # Model
    model::Model
    # System solver
    systemSolver::SystemSolver
    # TF to get linear solution
    getLinearSolution::Bool
    # Time variables
    initialTime::Real
    Δt::Union{Nothing,Real}
    finalTime::Union{Nothing,Real}
    timeVector::Union{Nothing,Vector{Float64}}
    adaptableΔt::Bool
    minΔt::Union{Nothing,Real}
    maxΔt::Union{Nothing,Real}
    # Boundary-crossing tolerance
    δb::Float64
    # Initial states update options
    skipInitialStatesUpdate::Bool
    initialVelocitiesUpdateOptions::InitialVelocitiesUpdateOptions
    # TF to track partial solutions at time steps and its frequency
    trackingTimeSteps::Bool
    trackingFrequency::Int64
    # TF to display progress and its frequency
    displayProgress::Bool
    displayFrequency::Int64
    # Secondary (outputs from problem creation)
    # -----------------------------------------
    # States, residual, Jacobian and inertia arrays
    x::Vector{Float64} = zeros(0)
    Δx::Vector{Float64} = zeros(0)
    residual::Vector{Float64} = zeros(0)
    jacobian::SparseMatrixCSC{Float64,Int64} = spzeros(0,0)
    inertia::SparseMatrixCSC{Float64,Int64} = spzeros(0,0)
    # Time variables
    timeNow::Real = 0
    timeBeginTimeStep::Real = 0
    timeEndTimeStep::Real = 0
    indexBeginTimeStep::Int64 = 1
    indexEndTimeStep::Int64 = 2
    sizeOfTime::Int64 = 1
    # Load factor
    σ::Float64 = 1.0
    # TF to compute only the external forces array at the current nonlinear step (used only for arclength system solver)
    getExternalForcesArray::Bool = false
    # TF for initial states being input
    initialStatesInput::Bool = false
    # TF to compute only residual array and skip Jacobian update at the current nonlinear step (used only in line search's Newton step)
    skipJacobianUpdate::Bool = false
    # Arrays of saved time steps, respective solutions and states
    savedTimeVector::Vector{Float64} = Vector{Float64}()
    xOverTime::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    elementalStatesOverTime::Vector{Vector{ElementalStates{Float64}}} = Vector{Vector{ElementalStates{Float64}}}()
    nodalStatesOverTime::Vector{Vector{NodalStates{Float64}}} = Vector{Vector{NodalStates{Float64}}}()
    compElementalStatesOverTime::Vector{Vector{ComplementaryElementalStates{Float64}}} = Vector{Vector{ComplementaryElementalStates{Float64}}}()
    elementalStatesRatesOverTime::Vector{Vector{ElementalStatesRates}} = Vector{Vector{ElementalStatesRates}}()
    compElementalStatesRatesOverTime::Vector{Vector{ComplementaryElementalStatesRates}} = Vector{Vector{ComplementaryElementalStatesRates}}()
    aeroVariablesOverTime::Vector{Vector{Union{Nothing,AeroVariables}}} = Vector{Vector{Union{Nothing,AeroVariables}}}()
    # Maximum absolute value of aerodynamic load coefficients over time
    maxAeroForce::Real = 0
    maxAeroMoment::Real = 0

end
export DynamicProblem


"""
    create_DynamicProblem(; kwargs...)

Dynamic problem constructor

# Keyword arguments
- `model::Model` = model
- `systemSolver::systemSolver` = nonlinear system solver
- `getLinearSolution::Bool` = flag to solve for linear structural solution
- `initialTime::Real` = initial time
- `Δt::Union{Nothing,Real}` = time step
- `finalTime::Union{Nothing,Real}` = final time
- `timeVector::Union{Nothing,Vector{Float64}}` = time vector
- `adaptableΔt::Bool` = flag for adaptable time step
- `minΔt::Union{Nothing,Real}` = minimum time step (when adaptable)
- `maxΔt::Union{Nothing,Real}` = maximum time step (when adaptable)
- `δb::Float64` = discontinuities boundary convergence norm
- `skipInitialStatesUpdate::Bool` = flag to skip update of initial states
- `initialVelocitiesUpdateOptions::InitialVelocitiesUpdateOptions` = options for the initial velocities update
- `trackingTimeSteps::Bool` = flag to track time step solutions
- `trackingFrequency::Int64` = frequency of time steps in which to track solution
- `displayProgress::Bool` = flag to display progress
- `displayFrequency::Int64` = frequency of time steps in which to display progress
- `x0::Vector{Float64}` = initial states
"""
function create_DynamicProblem(; model::Model,systemSolver::SystemSolver=create_NewtonRaphson(),getLinearSolution::Bool=false,initialTime::Real=0.0,Δt::Union{Nothing,Real}=nothing,finalTime::Union{Nothing,Real}=nothing,timeVector::Union{Nothing,Vector{Float64}}=nothing,adaptableΔt::Bool=false,minΔt::Union{Nothing,Real}=nothing,maxΔt::Union{Nothing,Real}=nothing,δb::Float64=-1e-5,skipInitialStatesUpdate::Bool=false,initialVelocitiesUpdateOptions::InitialVelocitiesUpdateOptions=InitialVelocitiesUpdateOptions(),trackingTimeSteps::Bool=true,trackingFrequency::Int64=1,displayProgress::Bool=true,displayFrequency::Int64=0,x0::Vector{Float64}=zeros(0))

    # Initialize problem
    problem = DynamicProblem(model=model,systemSolver=systemSolver,getLinearSolution=getLinearSolution,initialTime=initialTime,Δt=Δt,finalTime=finalTime,timeVector=timeVector,adaptableΔt=adaptableΔt,minΔt=minΔt,maxΔt=maxΔt,δb=δb,skipInitialStatesUpdate=skipInitialStatesUpdate,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions,trackingTimeSteps=trackingTimeSteps,trackingFrequency=trackingFrequency,displayProgress=displayProgress,displayFrequency=displayFrequency)

    # Check time step data
    if !isnothing(Δt)
        @assert Δt > 0
    end
    if !isnothing(minΔt)
        @assert !isnothing(Δt)
        @assert minΔt <= Δt/2
    end
    if !isnothing(maxΔt)
        @assert !isnothing(Δt)
        @assert maxΔt >= Δt
    end

    # Check boundary-crossing tolerance
    @assert δb < 0

    # Update initial load factor
    problem.σ = systemSolver.initialLoadFactor

    # Set initial elemental and nodal states
    if !isempty(x0)
        @assert length(x0) == model.systemOrder
        problem.initialStatesInput = true
        problem.x = x0*problem.σ
        problem.Δt = Inf64
        update_states!(problem)
        problem.Δt = Δt
    else
        set_initial_states!(problem)
        update_initial_aero_states!(problem,preInitialization=true)
    end   

    # Check and initialize time variables
    initialize_time_variables!(problem)

    # Initialize system arrays with correct size
    initialize_system_arrays!(problem)

    # Update display frequency if not input
    if displayProgress == true && displayFrequency == 0
        problem.displayFrequency = ceil(problem.sizeOfTime/100)
    end

    return problem

end
export create_DynamicProblem


# Sets the initial elemental and nodal states into the states array
function set_initial_states!(problem::Problem)

    @unpack x,model,σ,initialStatesInput = problem
    @unpack elements,specialNodes,systemOrder,nTrimVariables,forceScaling,rotationConstraints = model

    # Skip if states were input
    if initialStatesInput
        return
    end

    # Initialize states array
    x = zeros(systemOrder+nTrimVariables)

    # Loop over elements and assign initial states
    for element in elements
        @unpack DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ,DOF_δ = element
        @unpack u,p,F,M,V,Ω,χ = element.states
        x[DOF_u] = u
        x[DOF_p] = p
        x[DOF_F] = F/forceScaling
        x[DOF_M] = M/forceScaling
        x[DOF_V] = V
        x[DOF_Ω] = Ω
        x[DOF_χ] = χ
        x[DOF_δ] = isempty(DOF_δ) ? Vector{Float64}() : element.aero.δNow
    end

    # Loop over special nodes and assign initial states
    for specialNode in specialNodes
        @unpack BCs,uIsPrescribed,pIsPrescribed,DOF_uF,DOF_pM,DOF_trimLoads,u,p,F,M = specialNode
        # Loop directions and set generalized displacements/forces values accordingly
        for i=1:3
            x[DOF_uF[i]] = uIsPrescribed[i] ? F[i]/forceScaling : u[i]
            x[DOF_pM[i]] = pIsPrescribed[i] ? M[i]/forceScaling : p[i]
        end
        # If there are trim loads
        if any(!iszero(DOF_trimLoads))
            # Loop BCs
            for BC in BCs
                @unpack isTrim,initialTrimValue = BC
                # Set generalized trim forces values
                x[DOF_trimLoads[isTrim]] .= initialTrimValue[isTrim]/forceScaling
            end
        end
    end

    # Scale by initial load factor
    x *= σ

    # Set rotation constraints
    for constraint in rotationConstraints
        @unpack masterElemMasterGlobalDOF,slaveElemSlaveGlobalDOF,value = constraint
        x[slaveElemSlaveGlobalDOF] = x[masterElemMasterGlobalDOF] + value
    end

    @pack! problem = x

end


# Initializes the system arrays to the correct size
function initialize_system_arrays!(problem::Problem)

    @unpack model = problem
    @unpack systemOrder,nTrimVariables = model

    Δx = zeros(systemOrder+nTrimVariables)
    residual = zeros(systemOrder)
    jacobian = spzeros(systemOrder,systemOrder+nTrimVariables)
    inertia = spzeros(systemOrder,systemOrder)

    @pack! problem = Δx,residual,jacobian,inertia

end


"""
    solve!(problem::Problem)

Solves a problem 

# Arguments
- problem::Problem
"""
function solve!(problem::Problem)

    # Precompute distributed loads
    precompute_distributed_loads!(problem)

    # Solve according to problem type
    if typeof(problem) in [SteadyProblem,TrimProblem]
        solve_steady!(problem)
    elseif typeof(problem) == EigenProblem
        solve_steady!(problem)
        solve_eigen!(problem)    
    elseif typeof(problem) == DynamicProblem
        solve_dynamic!(problem)
    end

end
export solve!


# Pre-computes the distributed loads over all elements' nodes, over every time step
function precompute_distributed_loads!(problem::Problem)


    @unpack model = problem

    # Unpack time variables for dynamic problem or define dummy ones otherwise
    if typeof(problem) == DynamicProblem
        @unpack timeVector,sizeOfTime = problem
    else
        sizeOfTime = 1
        timeVector = [0.0]
    end

    # Shape function over element domain
    ϕ(n,ζ) = ifelse.(n==1, 1-ζ, ζ) 

    # Loop over beams
    for beam in model.beams
        # Loop over elements
        for element in beam.elements
            # Unpack
            @unpack Δℓ,f_A_of_ζt,m_A_of_ζt,f_b_of_ζt,m_b_of_ζt,ff_A_of_ζt,mf_A_of_ζt,ff_b_of_ζt,mf_b_of_ζt,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb = element
            # Initialize array of distributed loads on element's nodes over time
            f_A = zeros(2,3,sizeOfTime)
            m_A = zeros(2,3,sizeOfTime)
            f_b = zeros(2,3,sizeOfTime)
            m_b = zeros(2,3,sizeOfTime)
            ff_A = zeros(2,3,sizeOfTime)
            mf_A = zeros(2,3,sizeOfTime)
            ff_b = zeros(2,3,sizeOfTime)
            mf_b = zeros(2,3,sizeOfTime)
            # Loop over element's nodes
            for node = 1:2
                # Loop over time
                for (timeIndex,timeNow) in enumerate(timeVector)
                    # Define integrands
                    integrand_f_A(ζ) = Δℓ * f_A_of_ζt(ζ,timeNow) * ϕ(node,ζ)
                    integrand_m_A(ζ) = Δℓ * m_A_of_ζt(ζ,timeNow) * ϕ(node,ζ)
                    integrand_f_b(ζ) = Δℓ * f_b_of_ζt(ζ,timeNow) * ϕ(node,ζ)
                    integrand_m_b(ζ) = Δℓ * m_b_of_ζt(ζ,timeNow) * ϕ(node,ζ)
                    integrand_ff_A(ζ) = Δℓ * ff_A_of_ζt(ζ,timeNow) * ϕ(node,ζ)
                    integrand_mf_A(ζ) = Δℓ * mf_A_of_ζt(ζ,timeNow) * ϕ(node,ζ)
                    integrand_ff_b(ζ) = Δℓ * ff_b_of_ζt(ζ,timeNow) * ϕ(node,ζ)
                    integrand_mf_b(ζ) = Δℓ * mf_b_of_ζt(ζ,timeNow) * ϕ(node,ζ)
                    # Integrate nodal load at current time step, if applicable
                    if beam.hasDistributedDeadForcesBasisA
                        f_A[node,:,timeIndex], = quadgk(integrand_f_A, 0, 1)
                    end
                    if beam.hasDistributedDeadMomentsBasisA
                        m_A[node,:,timeIndex], = quadgk(integrand_m_A, 0, 1)
                    end
                    if beam.hasDistributedDeadForcesBasisb
                        f_b[node,:,timeIndex], = quadgk(integrand_f_b, 0, 1)
                    end
                    if beam.hasDistributedDeadMomentsBasisb
                        m_b[node,:,timeIndex], = quadgk(integrand_m_b, 0, 1)
                    end
                    if beam.hasDistributedFollowerForcesBasisA
                        ff_A[node,:,timeIndex], = quadgk(integrand_ff_A, 0, 1)
                    end
                    if beam.hasDistributedFollowerMomentsBasisA
                        mf_A[node,:,timeIndex], = quadgk(integrand_mf_A, 0, 1)
                    end
                    if beam.hasDistributedFollowerForcesBasisb
                        ff_b[node,:,timeIndex], = quadgk(integrand_ff_b, 0, 1)
                    end
                    if beam.hasDistributedFollowerMomentsBasisb
                        mf_b[node,:,timeIndex], = quadgk(integrand_mf_b, 0, 1)
                    end
                end
            end
            # Update TFs for any nonzero loads on the element over time
            if any(!iszero(f_A))
                hasDistributedDeadForcesBasisA = true
            end
            if any(!iszero(m_A))
                hasDistributedDeadMomentsBasisA = true
            end
            if any(!iszero(f_b))
                hasDistributedDeadForcesBasisb = true
            end
            if any(!iszero(m_b))
                hasDistributedDeadMomentsBasisb = true
            end
            if any(!iszero(ff_A))
                hasDistributedFollowerForcesBasisA = true
            end
            if any(!iszero(mf_A))
                hasDistributedFollowerMomentsBasisA = true
            end
            if any(!iszero(ff_b))
                hasDistributedFollowerForcesBasisb = true
            end
            if any(!iszero(mf_b))
                hasDistributedFollowerMomentsBasisb = true
            end
            # Pack element variables
            @pack! element = f_A,m_A,f_b,m_b,ff_A,mf_A,ff_b,mf_b,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb
        end
    end
end


# Solves a steady problem (includes trim and the steady part of eigen problems) 
function solve_steady!(problem::Problem)

    @unpack getLinearSolution,systemSolver = problem

    if getLinearSolution
        problem.σ = 1.0
        assemble_system_arrays!(problem)
        solve_linear_system!(problem)
        update_states!(problem)
        if problem isa EigenProblem || (problem isa TrimProblem && problem.getInertiaMatrix)
            for element in problem.model.elements
                element_inertia!(problem,problem.model,element)
            end
        end
        save_load_factor_data!(problem,problem.σ,problem.x)
    else
        if typeof(systemSolver) == NewtonRaphson
            solve_NewtonRaphson!(problem)
        end
    end

end


"""
    solve_eigen!(problem::EigenProblem)

Solves an eigenproblem 

# Arguments
- `problem::EigenProblem`
"""
function solve_eigen!(problem::EigenProblem)

    @unpack jacobian,inertia,frequencyFilterLimits,nModes = problem
       
    ## Process eigenvalues and eigenvectors
    #---------------------------------------------------------------------------
    # Solve eigenproblem to get eigenvectors and inverse of eigenvalues
    if problem.model.nTrimVariables > 0 || isapprox(det(jacobian),0)
        inverseEigenvalues, eigenvectors = eigen(-pinv(Matrix(jacobian))*inertia)
    else
        inverseEigenvalues, eigenvectors = eigen(-jacobian\Matrix(inertia))
    end
    # Get eigenvalues
    eigenvalues = 1.0 ./ inverseEigenvalues 
    # Sort in ascending order of frequency
    ind = sortperm(abs.(imag.(eigenvalues)))      
    eigenvalues = eigenvalues[ind]
    eigenvectors = eigenvectors[:,ind]
    # Frequencies [rad/s] are the absolute imaginary values of eigenvalues pairs
    frequencies = abs.(imag.(eigenvalues[1:2:end]))
    # Dampings [rad/s] are the real values of eigenvalues pairs      
    dampings = real.(eigenvalues[1:2:end]) 
    # Keep only one of the pair of eigenvectors
    eigenvectors = eigenvectors[:,2:2:end]
    
    ## Filter by frequency
    #---------------------------------------------------------------------------
    # Get filtered indices
    filteredIndices = findall( x -> x >= frequencyFilterLimits[1] && x <= frequencyFilterLimits[2], frequencies)
    # Filtered frequencies, dampings and eigenvectors 
    frequenciesFiltered = frequencies[filteredIndices]
    dampingsFiltered = dampings[filteredIndices]
    eigenvectorsFiltered = eigenvectors[:,filteredIndices]

    ## Non-oscillatory (aerodynamic, divergence) modes
    #---------------------------------------------------------------------------
    # Get indices of non-oscillatory (zero frequency) modes
    nonOscillatoryIndices = findall( x -> isapprox(x,0), frequencies)
    # Non-oscillatory dampings 
    dampingsNonOscillatory = real.(eigenvalues[nonOscillatoryIndices])
    # Non-oscillatory eigenvectors
    eigenvectorsNonOscillatory = eigenvectors[:,nonOscillatoryIndices]

    ## Get oscillatory (structural, dynamic aeroelastic, flight dynamics) modes
    #---------------------------------------------------------------------------# Starting index of oscillatory modes 
    indexStartOscillatory = findfirst( x -> x > 0, frequenciesFiltered)
    # Frequencies, dampings and eigenvectors of oscillatory modes
    if isnothing(indexStartOscillatory)
        frequenciesOscillatory,dampingsOscillatory,eigenvectorsOscillatoryCplx = Vector{Float64}(),Vector{Float64}(),zeros(0, 0)
    else
        frequenciesOscillatory = frequenciesFiltered[indexStartOscillatory:end]
        dampingsOscillatory = dampingsFiltered[indexStartOscillatory:end]
        eigenvectorsOscillatoryCplx = eigenvectorsFiltered[:,indexStartOscillatory:end]
    end

    ## Filter by number of desired modes
    #---------------------------------------------------------------------------
    # Available number of modes
    availableNModes = length(frequenciesOscillatory)
    # Filter, if possible                 
    if availableNModes >= nModes
        frequenciesOscillatory = frequenciesOscillatory[1:nModes]
        dampingsOscillatory = dampingsOscillatory[1:nModes]
        eigenvectorsOscillatoryCplx = eigenvectorsOscillatoryCplx[:,1:nModes]
    else
        nModes = availableNModes
        display("Only $availableNModes modes could be calculated, due to limited number of elements or frequency filter limits")
    end
    
    ## Post-process oscillatory eigenvectors
    #---------------------------------------------------------------------------
    # Real and imaginary parts of eigenvector
    eigenvectorsRe = real.(eigenvectorsOscillatoryCplx)
    eigenvectorsIm = imag.(eigenvectorsOscillatoryCplx)
    eigenvectorsNonOscillatoryRe = real.(eigenvectorsNonOscillatory)
    eigenvectorsNonOscillatoryIm = imag.(eigenvectorsNonOscillatory)
    # Indices whose real part is larger
    realLargerOsc = @. abs(eigenvectorsRe) > abs(eigenvectorsIm)   
    realLargerNonOsc = @. abs(eigenvectorsNonOscillatoryRe) > abs(eigenvectorsNonOscillatoryIm)
    # Indices whose imaginary part is larger
    imagLargerOsc = .!realLargerOsc
    imagLargerNonOsc = .!realLargerNonOsc                                       
    # Sign of larger part determines true sign of solution
    signsOsc = @. sign(eigenvectorsRe) * realLargerOsc + sign(eigenvectorsIm) * imagLargerOsc
    signsNonOsc = @. sign(eigenvectorsNonOscillatoryRe) * realLargerNonOsc + sign(eigenvectorsNonOscillatoryIm) * imagLargerNonOsc                  
    # Correctly signed absolute value of eigenvectors
    eigenvectorsOscillatoryAbs = @. signsOsc * abs(eigenvectorsOscillatoryCplx)
    eigenvectorsNonOscillatoryAbs = @. signsNonOsc * abs(eigenvectorsNonOscillatory)
    # Get mode shapes of complex-valued and absolute-valued eigenvectors
    modeShapesCplx = get_mode_shapes!(problem,eigenvectorsOscillatoryCplx,frequenciesOscillatory,dampingsOscillatory)
    modeShapesAbs = get_mode_shapes!(problem,eigenvectorsOscillatoryAbs,frequenciesOscillatory,dampingsOscillatory)
    modeShapesAbsNonOsc = get_mode_shapes!(problem,eigenvectorsNonOscillatoryAbs,zeros(length(dampingsNonOscillatory)),dampingsNonOscillatory)

    @pack! problem = nModes,frequencies,dampings,eigenvectors,frequenciesFiltered,dampingsFiltered,eigenvectorsFiltered,dampingsNonOscillatory,frequenciesOscillatory,dampingsOscillatory,eigenvectorsOscillatoryCplx,eigenvectorsOscillatoryAbs,eigenvectorsNonOscillatory,modeShapesCplx,modeShapesAbs,modeShapesAbsNonOsc

end
export solve_eigen!


# Gets the mode shapes given by the eigenvectors
function get_mode_shapes!(problem::Problem,eigenvectors::Matrix{T},frequencies::Vector{Float64},dampings::Vector{Float64}) where T<:Union{Float64,ComplexF64}

    @unpack model,normalizeModeShapes = problem
    @unpack forceScaling,elements,nElementsTotal = model

    # Initialize array of mode shapes
    modeShapes = Vector{ModeShape{T}}()

    # Loop modes
    for (mode,freq,damp) in zip(1:size(eigenvectors,2),frequencies,dampings)
        # Initialize current mode shape
        modeShape = ModeShape{T}(mode,freq,damp,nElementsTotal)
        # Current mode's eigenvector
        eigenvector = eigenvectors[:,mode]
        # Initialize arrays
        uArray,pArray,FArray,MArray,VArray,ΩArray,γArray,κArray,PArray,HArray,θArray = Vector{T}(),Vector{T}(),Vector{T}(),Vector{T}(),Vector{T}(),Vector{T}(),Vector{T}(),Vector{T}(),Vector{T}(),Vector{T}(),Vector{T}()
        # Loop elements
        for (e,element) in enumerate(elements)
            # Get modal states for current element
            u,p,F,M,V,Ω,γ,κ,P,H,u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2 = element_modal_states(element,eigenvector,forceScaling)
            # Manually correct nodal values
            if e > 1
                u_n1,p_n1,u_n1_b,p_n1_b,F_n1,M_n1,θ_n1 = fix_nodal_modal_states!(model,element,modeShape,u_n1,p_n1,u_n1_b,p_n1_b,F_n1,M_n1,θ_n1)
            end
            # Add states to current mode 
            @pack! modeShape.elementalStates[e] = u,p,F,M,V,Ω
            @pack! modeShape.complementaryElementalStates[e] = γ,κ,P,H
            @pack! modeShape.nodalStates[e] = u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2
            # Update arrays
            append!(uArray,u_n1,u,u_n2)
            append!(pArray,p_n1,p,p_n2)
            append!(FArray,F_n1,F,F_n2)
            append!(MArray,M_n1,M,M_n2)
            append!(VArray,V)
            append!(ΩArray,Ω)
            append!(γArray,γ)
            append!(κArray,κ)
            append!(PArray,P)
            append!(HArray,H)
            append!(θArray,θ_n1,θ_n2)
        end
        # Normalize mode shapes by maxima, if applicable
        if normalizeModeShapes
            # Compute norms
            upNorm = maximum(abs.(vcat(uArray,pArray))) > 0.0 ? maximum(abs.(vcat(uArray,pArray))) : 1.0
            FNorm = maximum(abs.(FArray)) > 0.0 ? maximum(abs.(FArray)) : 1.0
            MNorm = maximum(abs.(MArray)) > 0.0 ? maximum(abs.(MArray)) : 1.0
            VNorm = maximum(abs.(VArray)) > 0.0 ? maximum(abs.(VArray)) : 1.0
            ΩNorm = maximum(abs.(ΩArray)) > 0.0 ? maximum(abs.(ΩArray)) : 1.0
            γNorm = maximum(abs.(γArray)) > 0.0 ? maximum(abs.(γArray)) : 1.0
            κNorm = maximum(abs.(κArray)) > 0.0 ? maximum(abs.(κArray)) : 1.0
            PNorm = maximum(abs.(PArray)) > 0.0 ? maximum(abs.(PArray)) : 1.0
            HNorm = maximum(abs.(HArray)) > 0.0 ? maximum(abs.(HArray)) : 1.0
            θNorm = maximum(abs.(θArray)) > 0.0 ? maximum(abs.(θArray)) : 1.0
            # Loop elements 
            for (e,element) in enumerate(elements)
                # Unpack
                @unpack u,p,F,M,V,Ω = modeShape.elementalStates[e]
                @unpack γ,κ,P,H = modeShape.complementaryElementalStates[e]
                @unpack u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2 = modeShape.nodalStates[e] 
                u,u_n1,u_n2,u_n1_b,u_n2_b = divide_inplace!(upNorm, u,u_n1,u_n2,u_n1_b,u_n2_b) 
                p,p_n1,p_n2,p_n1_b,p_n2_b = divide_inplace!(upNorm, p,p_n1,p_n2,p_n1_b,p_n2_b) 
                F,F_n1,F_n2 = divide_inplace!(FNorm, F,F_n1,F_n2)
                M,M_n1,M_n2 = divide_inplace!(MNorm, M,M_n1,M_n2)
                V ./= VNorm
                Ω ./= ΩNorm
                γ ./= γNorm
                κ ./= κNorm
                P ./= PNorm
                H ./= HNorm
                θ_n1,θ_n2 = divide_inplace!(θNorm, θ_n1,θ_n2)
                @pack! modeShape.elementalStates[e] = u,p,F,M,V,Ω
                @pack! modeShape.complementaryElementalStates[e] = γ,κ,P,H
                @pack! modeShape.nodalStates[e] = u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2
            end
        end
        # Add to mode shapes array
        push!(modeShapes,modeShape)
    end
    
    return modeShapes
end


# Corrects the nodal modal states which may not coincide depending on which element is considered
function fix_nodal_modal_states!(model,element,modeShape,u_n1,p_n1,u_n1_b,p_n1_b,F_n1,M_n1,θ_n1)

    @unpack elementNodes = model

    # Find global ID of the first element that contains the local first node
    firstElem = findfirst(row -> element.nodesGlobalID[1] in row, elementNodes)

    # Skip if this is itself the first element
    if firstElem == element.globalID
        return u_n1,p_n1,u_n1_b,p_n1_b,F_n1,M_n1,θ_n1
    end

    # Set value of local first node as that of the second local node of that element
    u_n1 = modeShape.nodalStates[firstElem].u_n2
    p_n1 = modeShape.nodalStates[firstElem].p_n2
    u_n1_b = modeShape.nodalStates[firstElem].u_n2_b
    p_n1_b = modeShape.nodalStates[firstElem].p_n2_b
    F_n1 = modeShape.nodalStates[firstElem].F_n2
    M_n1 = modeShape.nodalStates[firstElem].M_n2
    θ_n1 = modeShape.nodalStates[firstElem].θ_n2

    return u_n1,p_n1,u_n1_b,p_n1_b,F_n1,M_n1,θ_n1
end


# Solves a dynamic problem 
function solve_dynamic!(problem::Problem)

    # Solve the initial time step to get consistent initial states, if applicable
    if !problem.skipInitialStatesUpdate
        solve_initial_dynamic!(problem)
    end
    # Save initial states
    save_time_step_data!(problem,problem.timeNow)
    # Time march
    if problem.adaptableΔt
        adaptable_time_march!(problem)
    else
        time_march!(problem)
    end

end


# Validates and initializes the time variables
function initialize_time_variables!(problem::Problem)

    @unpack initialTime,finalTime,Δt,timeVector,initialVelocitiesUpdateOptions = problem

    # Validate and assign time vector
    if isnothing(timeVector)
        # Final time and time step must be inputs
        @assert !isnothing(initialTime)
        @assert !isnothing(Δt)
        # Check final and initial time inputs
        @assert finalTime >= initialTime + Δt
        # Set time vector
        timeVector = collect(initialTime:Δt:finalTime)
        @assert length(timeVector) > 1
    else
        # Check that time vector is monotonically increasing
        @assert timeVector == sort!(copy(timeVector))
        @assert length(timeVector) > 1
        # Set initial and final times, and initial time step
        initialTime = timeVector[1]
        Δt = timeVector[2] - timeVector[1]
        finalTime = timeVector[end]
    end
    
    # Assign other time variables
    timeNow = initialTime
    timeBeginTimeStep = initialTime
    timeEndTimeStep = initialTime+Δt
    sizeOfTime = length(timeVector)

    # Update time step for initial velocities update, if not input
    if initialVelocitiesUpdateOptions.Δt == 0
        initialVelocitiesUpdateOptions.Δt = 1e-2*Δt
    end

    @pack! problem = initialTime,Δt,finalTime,timeVector,timeNow,timeBeginTimeStep,timeEndTimeStep,sizeOfTime,initialVelocitiesUpdateOptions

end


# Gets consistent initial conditions and solves the first time step
function solve_initial_dynamic!(problem::Problem)

    problemCopy = deepcopy(problem)
    modelCopy = problemCopy.model
    @unpack BCs,BCedNodes = modelCopy
    
    # Initialize generalized displacements BCs 
    initialDisplacementsBCs = BCs

    # Initialize flag for node having being assigned BCs
    assigned = falses(modelCopy.nNodesTotal)
    assigned[BCedNodes] .= true

    # List of initial generalized displacements types (defined in basis b)
    displacementTypes = ["u1b","u2b","u3b","p1b","p2b","p3b"]

    # Loop elements
    for (e, element) in enumerate(modelCopy.elements)
        @unpack u_n1,u_n2,p_n1,p_n2 = element.nodalStates
        # Loop element's nodes
        for (n, localNode) in enumerate(element.nodesLocalID[1]:element.nodesLocalID[2])
            globalNode = element.nodesGlobalID[n]
            if !assigned[globalNode]
                # Update flag
                assigned[globalNode] = true
                # Initial nodal generalized displacements
                initialNodalDisplacements = n == 1 ? vcat(u_n1,p_n1) : vcat(u_n2,p_n2)
                # Indices of nonzero initial generalized displacements
                nonzeroIndices = findall(!iszero,initialNodalDisplacements)
                # Nonzero initial generalized displacements values and types
                nonzeroValues = initialNodalDisplacements[nonzeroIndices]
                nonzeroTypes = displacementTypes[nonzeroIndices]
                # Create BCs 
                bc = create_BC(beam=element.parent,node=localNode,types=nonzeroTypes,values=nonzeroValues)
                # Add BCs to list
                push!(initialDisplacementsBCs,bc)
            end
        end  
    end

    # Update model as if all nodes were BC'ed to the specified initial displacements
    modelCopy.BCs = initialDisplacementsBCs
    set_BCed_nodes!(modelCopy)
    modelCopy.skipValidationMotionBasisA = true
    update_model!(modelCopy)
    problemCopy.model = modelCopy
 
    # Set initial states
    problemCopy.initialStatesInput = false
    set_initial_states!(problemCopy)

    # Initialize system arrays with correct size
    initialize_system_arrays!(problemCopy)

    # Set infinite Δt
    problemCopy.Δt = Inf64

    # Get equivalent initial states rates
    get_equivalent_states_rates!(problemCopy)

    # Solve to find unspecified initial displacements and/or rotations
    solve_steady!(problemCopy)

    # Copy initial states to the original problem
    copy_initial_states!(problem,problemCopy)

    # Update initial velocities states
    update_initial_velocities!(problem)

    # Update initial aerodynamic states
    update_initial_aero_states!(problem)
    
end


# Copies the initial states to the original problem
function copy_initial_states!(problem::Problem,problemCopy::Problem)

    # Unpack problem copy
    xCopy = problemCopy.x
    modelCopy = problemCopy.model

    # Unpack original problem
    @unpack x,model = problem
    @unpack forceScaling = model

    # Loop over elements
    for (element,elementCopy) in zip(model.elements,modelCopy.elements)
        # Unpack DOFs
        @unpack states,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω = element
        DOF_uCopy,DOF_pCopy,DOF_FCopy,DOF_MCopy,DOF_VCopy,DOF_ΩCopy = elementCopy.DOF_u,elementCopy.DOF_p,elementCopy.DOF_F,elementCopy.DOF_M,elementCopy.DOF_V,elementCopy.DOF_Ω
        # Set elemental states
        u = x[DOF_u] = xCopy[DOF_uCopy]
        p = x[DOF_p] = xCopy[DOF_pCopy]
        F = x[DOF_F] = xCopy[DOF_FCopy]
        M = x[DOF_M] = xCopy[DOF_MCopy]
        V = x[DOF_V] = xCopy[DOF_VCopy]
        Ω = x[DOF_Ω] = xCopy[DOF_ΩCopy]
        F,M = forceScaling*F,forceScaling*M
        # Pack element states
        @pack! element.states = u,p,F,M,V,Ω
    end

    # Loop over special nodes
    for specialNode in model.specialNodes
        # Get equivalent copied special node
        matchingSpecialNode(s) = s.globalID == specialNode.globalID
        specialNodeCopyIndex = findfirst(matchingSpecialNode,modelCopy.specialNodes)
        specialNodeCopy = modelCopy.specialNodes[specialNodeCopyIndex]
        # Unpack DOFs
        @unpack DOF_uF,DOF_pM = specialNode
        DOF_uFCopy,DOF_pMCopy = specialNodeCopy.DOF_uF,specialNodeCopy.DOF_pM
        # Set nodal states
        x[DOF_uF] = xCopy[DOF_uFCopy]
        x[DOF_pM] = xCopy[DOF_pMCopy]
    end

    @pack! problem = x,model
end


# Updates the velocity states by running a very small time step
function update_initial_velocities!(problem::Problem)

    @unpack timeVector,Δt,timeNow,indexBeginTimeStep,indexEndTimeStep,timeBeginTimeStep,timeEndTimeStep,initialVelocitiesUpdateOptions,model = problem
    @unpack maxIter,relaxFactor,tol,displayProgress = initialVelocitiesUpdateOptions
    @unpack elements = model

    # Get original values of time variables and the system solver
    ΔtOriginal = Δt
    timeVectorOriginal = timeVector
    systemSolverOriginal = problem.systemSolver

    # Set very small time step and update time variables
    Δt = initialVelocitiesUpdateOptions.Δt
    timeVector = [timeVector[1],timeVector[1]+Δt]
    timeNow = timeVector[1] + Δt
    indexBeginTimeStep = 1
    indexEndTimeStep = 2
    timeBeginTimeStep = timeNow - Δt 
    timeEndTimeStep = timeNow 
    @pack! problem = timeVector,Δt,timeNow,indexBeginTimeStep,indexEndTimeStep,timeBeginTimeStep,timeEndTimeStep
    # Update basis A orientation
    update_basis_A_orientation!(problem,update_R_A_array=false)
    # Update boundary conditions
    for BC in model.BCs
        update_BC_data!(BC,timeNow)
    end   
    # Set default system solver (with reduced relative tolerance and large number of maximum iterations)
    problem.systemSolver = create_NewtonRaphson(maximumIterations=200,relativeTolerance=1e-10)
    # Initial displacements, sectional velocities and displacement's rates
    disp0 = [[e.states.u; e.states.p] for e in elements]
    vel0 = [[e.states.V; e.states.Ω] for e in elements]
    dispRates0 = [[e.statesRates.udot; e.statesRates.pdot] for e in elements]
    # Set velocity indices to update
    DoFToUpdate = [e.parent.velDoFToUpdate for e in elements]
    # Initialize convergence variables
    ϵ = 1
    iter = 0
    # Loop until converged velocities are found
    while ϵ > tol
        # Update iteration counter
        iter += 1
        # Loop over elements: Set relaxation factor on velocity states in need of update to stabilize
        if iter > 1
            for (e,element) in enumerate(elements)
                # Unpack
                @unpack x = problem
                @unpack DOF_V,DOF_Ω = element
                @unpack V,Ω = element.states
                @unpack udot,pdot = element.statesRates
                # Linear and angular velocities' indices
                Vind = 1:3
                Ωind = 4:6
                # DOF to update for current element
                VIndToUpdate = DoFToUpdate[e][Vind]
                ΩIndToUpdate = DoFToUpdate[e][Ωind]
                # Reset velocities not to be updated
                udot[.!VIndToUpdate] .= dispRates0[e][Vind[.!VIndToUpdate]]
                pdot[.!ΩIndToUpdate] .= dispRates0[e][Ωind[.!ΩIndToUpdate]]
                x[DOF_V[.!VIndToUpdate]] = V[.!VIndToUpdate] = vel0[e][Vind[.!VIndToUpdate]]
                x[DOF_Ω[.!ΩIndToUpdate]] = Ω[.!ΩIndToUpdate] = vel0[e][Ωind[.!ΩIndToUpdate]]
                # Update velocities with relaxation factor, for applicable DOF
                udot[VIndToUpdate] = udot[VIndToUpdate]*relaxFactor + dispRates0[e][Vind[VIndToUpdate]]*(1-relaxFactor)
                pdot[ΩIndToUpdate] = pdot[ΩIndToUpdate]*relaxFactor + dispRates0[e][Ωind[ΩIndToUpdate]]*(1-relaxFactor)
                x[DOF_V[VIndToUpdate]] = V[VIndToUpdate] = V[VIndToUpdate]*relaxFactor + vel0[e][Vind[VIndToUpdate]]*(1-relaxFactor)
                x[DOF_Ω[ΩIndToUpdate]] = Ω[ΩIndToUpdate] = Ω[ΩIndToUpdate]*relaxFactor + vel0[e][Ωind[ΩIndToUpdate]]*(1-relaxFactor)
                # Pack
                @pack! problem = x
                @pack! element.states = V,Ω
                @pack! element.statesRates = udot,pdot
            end
        end
        # Sectional velocities and accelerations arrays at begin of time step 
        velBegin = [[e.states.V; e.states.Ω] for e in elements]
        accBegin = [[e.statesRates.Vdot; e.statesRates.Ωdot] for e in elements]
        # Reset accelerations to zero 
        for element in elements
            Vdot,Ωdot = zeros(3),zeros(3)
            @pack! element.statesRates = Vdot,Ωdot
        end
        # Get equivalent states' rates at the begin of the time step
        get_equivalent_states_rates!(problem)
        # Solve the system at the current time step
        solve_time_step!(problem)
        # Reset displacements to values at the begin of time step
        for (e,element) in enumerate(elements)
            u,p = disp0[e][1:3],disp0[e][4:6]
            @pack! element.states = u,p
        end
        # Sectional velocities and accelerations arrays at end of time step 
        velEnd = [[e.states.V; e.states.Ω] for e in elements]
        accEnd = [[e.statesRates.Vdot; e.statesRates.Ωdot] for e in elements]
        # Update current difference between initial and final velocities 
        ϵ = maximum(maximum([abs.(velBegin[e].-velEnd[e]) .* DoFToUpdate[e] for e in eachindex(elements)]))
        # Print iteration results
        if displayProgress
            println("iter: $iter, ϵ = $ϵ")
        end
        # Check number of iterations
        if iter == maxIter
            println("Unconverged initial velocities - solution of generalized velocities and accelerations may be unstable: iter = $iter, ϵ = $ϵ")
            break
        end
    end

    # Reset original time variables and system solver
    Δt = ΔtOriginal
    timeVector = timeVectorOriginal
    timeNow = timeVector[1] 
    indexBeginTimeStep = 1
    indexEndTimeStep = 2
    timeBeginTimeStep = timeNow
    timeEndTimeStep = timeNow + Δt
    problem.systemSolver = systemSolverOriginal

    @pack! problem = timeVector,Δt,timeNow,indexBeginTimeStep,indexEndTimeStep,timeBeginTimeStep,timeEndTimeStep,model
end


# Marches the dynamic problem in time
function time_march!(problem::Problem)

    @unpack model,sizeOfTime,trackingTimeSteps,trackingFrequency,displayProgress,displayFrequency = problem

    # Advance time
    for timeIndex = 2:sizeOfTime    
        # Update time variables 
        update_time_variables!(problem,timeIndex)
        @unpack timeNow = problem
        # Update basis A orientation
        update_basis_A_orientation!(problem)
        # Update boundary conditions
        for BC in model.BCs
            update_BC_data!(BC,timeNow)
        end   
        # Get equivalent states' rates at the begin of the time step
        get_equivalent_states_rates!(problem)
        # Update DS model complementary variables of previous time step, if applicable
        update_BL_complementary_variables!(problem)
        # Solve the system at the current time step
        solve_time_step!(problem)
        # Stop if unconverged
        if !problem.systemSolver.convergedFinalSolution
            println("Unconverged solution, stopping...")
            break
        end
        # Save time step data, if applicable    
        if trackingTimeSteps && rem(timeIndex,trackingFrequency) == 0
            save_time_step_data!(problem,timeNow)           
        end
        # Display progress, if applicable
        if displayProgress && (rem(timeIndex,displayFrequency) == 0 || timeIndex == sizeOfTime)
            progress = round(timeIndex/sizeOfTime*100,digits=1)
            println("Simulation progress: $progress %")
        end
    end
end


# Marches the dynamic problem in time using an adaptable time step
function adaptable_time_march!(problem::Problem)

    @unpack model,trackingTimeSteps,trackingFrequency,displayProgress,displayFrequency,timeVector,timeNow,Δt,minΔt,maxΔt,δb = problem

    # Time at begin of step
    timeInitStep = deepcopy(timeNow)

    # Initialize time step index, and flag for reduced time step
    timeIndex = 1
    reducedΔtstep = false

    # Advance time
    while timeInitStep < timeVector[end]
        # Store problem state at begin of time step
        stepInitState = copy_state(problem)
        # Limite time step
        Δt = max(minΔt, min(Δt, timeVector[end] - timeInitStep))
        # Update time variables 
        timeNow += Δt
        @pack! problem = timeNow,Δt
        # Update basis A orientation
        update_basis_A_orientation!(problem)
        # Update boundary conditions
        for BC in model.BCs
            update_BC_data!(BC,timeNow)
        end   
        # Get equivalent states' rates at the begin of the time step
        get_equivalent_states_rates!(problem)
        # Update DS model complementary variables of previous time step, if applicable
        if !reducedΔtstep
            update_BL_complementary_variables!(problem)
        end
        # Solve the system at the current time step
        solve_time_step!(problem)
        # Get maximum stall boundary for DS models
        boundaryMax = BL_stall_boundary(problem)
        # Save time step data, if applicable    
        if trackingTimeSteps && rem(timeIndex,trackingFrequency) == 0
            save_time_step_data!(problem,timeNow)           
        end
        # If unconverged or in boundary crossing, reduce time step (if possible), restore problem state at begin of time step and update flag
        if !problem.systemSolver.convergedFinalSolution || boundaryMax < δb 
            println(boundaryMax)
            Δt *= 1/2
            if Δt < minΔt
                println("Minimum time step reached, stopping...")
                return
            end
            restore_state!(problem,stepInitState)
            reducedΔtstep = true
            continue
        else
            reducedΔtstep = false
        end
        # Increase time step, if possible
        Δt = min(2*Δt, maxΔt)
        # Display progress, if applicable
        if displayProgress && rem(timeIndex,displayFrequency) == 0
            progress = round(timeNow/timeVector[end]*100,digits=1)
            println("Simulation progress: $progress %")
        end
        # Update time at begin of time step and time index
        timeInitStep = deepcopy(timeNow)
        timeIndex += 1
        # Update time vector, if applicable
        if reducedΔtstep
            insert!(timeVector,timeIndex,timeNow)
        end
    end

    @pack! problem = timeVector
end


# Copies the current state of the problem
function copy_state(problem::Problem)

    elementStates = [deepcopy(element.states) for element in problem.model.elements]
    elementStatesRates = [deepcopy(element.statesRates) for element in problem.model.elements]
    σ = deepcopy(problem.σ)
    x = deepcopy(problem.x)
    elementBLiCompVars = Vector{BLiComplementaryVariables}(undef,problem.model.nElementsTotal)
    for (e,element) in enumerate(problem.model.elements)
        if isnothing(element.aero) || !(typeof(element.aero.solver) in [BLi])
            continue
        end
        elementBLiCompVars[e] = deepcopy(element.aero.BLiCompVars)
    end
    
    return (elementStates,elementStatesRates,σ,x,elementBLiCompVars)
end


# Restores the state of the problem
function restore_state!(problem::Problem,stepInitState)

    for (e,element) in enumerate(problem.model.elements)
        element.states = stepInitState[1][e]
        element.statesRates = stepInitState[2][e]
        if isnothing(element.aero) || !(typeof(element.aero.solver) in [BLi])
            continue
        end
        element.aero.BLiCompVars = stepInitState[5][e]
    end
    problem.σ = stepInitState[3]
    problem.x = stepInitState[4]
    
end


# Updates the time variables (time, time step, time indices)
function update_time_variables!(problem::Problem,timeIndex::Int64)

    @unpack timeVector = problem

    timeNow = timeVector[timeIndex]
    Δt = timeVector[timeIndex]-timeVector[timeIndex-1]
    indexBeginTimeStep = timeIndex-1
    indexEndTimeStep = timeIndex 
    timeBeginTimeStep = timeVector[indexBeginTimeStep]
    timeEndTimeStep = timeVector[indexEndTimeStep]

    @pack! problem = timeNow,Δt,indexBeginTimeStep,indexEndTimeStep,timeBeginTimeStep,timeEndTimeStep

end


# Updates the orientation of the basis A for the next time step
function update_basis_A_orientation!(problem::Problem;update_R_A_array::Bool=true)

    @unpack model,timeNow,Δt = problem
    @unpack R_A,ω_A = model

    # Skew-symmetric operator of the angular velocity vector
    ω_A_tilde = tilde(ω_A(timeNow))

    # Rotation tensor from basis I to basis A, at the current time, and its transpose
    R_A = inv((2/Δt*I3-ω_A_tilde))*(2/Δt*I3+ω_A_tilde)*R_A
    round_off!(R_A)
    R_AT = Matrix(R_A')

    # Add to rotation tensors array over time, if applicable
    if update_R_A_array
        push!(problem.model.R_A_ofTime, R_A)
    end

    @pack! model = R_A,R_AT

end


# Gets the equivalent states' rates at the begin of the current time step 
function get_equivalent_states_rates!(problem::Problem)

    @unpack model,Δt = problem
    @unpack elements = model

    # Loop over elements
    for element in elements
        # Unpack element data (element states and rates known at the begin of time step)
        @unpack u,p,V,Ω,χ = element.states
        @unpack udot,pdot,Vdot,Ωdot,χdot = element.statesRates
        # Equivalent states' rates at the begin of time step 
        udotEquiv = udot + 2/Δt*u
        pdotEquiv = pdot + 2/Δt*p
        VdotEquiv = Vdot + 2/Δt*V
        ΩdotEquiv = Ωdot + 2/Δt*Ω
        χdotEquiv = χdot + 2/Δt*χ
        # Pack element data
        @pack! element = udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,χdotEquiv
    end
end


# Updates the complementary variables of previous time step for the dynamic stall models
function update_BL_complementary_variables!(problem::Problem)

    # Loop elements
    for element in problem.model.elements
        # Skip elements without an aerodynamic surface with dynamic stall solver
        if isnothing(element.aero) || !(typeof(element.aero.solver) in [BLi,BLo])
            continue
        end
        # Unpack variables for BLi model
        @unpack αlag = element.aero.BLiStates
        @unpack qR = element.aero.BLiKin
        @unpack upstroke,P,stallOnsetRatio = element.aero.BLiFlow
        # Update and pack variables for BLi model
        stallOnsetRatioPrev,αlagPrev,qRPrev,PPrev,upstrokePrev = stallOnsetRatio,αlag,qR,P,upstroke
        @pack! element.aero.BLiCompVars = stallOnsetRatioPrev,αlagPrev,qRPrev,PPrev,upstrokePrev
        # Unpack variables for BLo model
        @unpack stallOnsetRatio = element.aero.BLoFlow
        # Update and pack variables for BLo model
        stallOnsetRatioPrev = stallOnsetRatio
        @pack! element.aero.BLoCompVars = stallOnsetRatioPrev
    end

end


# Computes the maximum stall boundary value over all beam elements in the Beddoes-Leishman model
function BL_stall_boundary(problem::Problem)

    # Initialize
    boundary = Vector{Float64}(undef,problem.model.nElementsTotal)
    hasDSelements = false

    # Loop elements
    for (e,element) in enumerate(problem.model.elements)
        # Skip elements without an aerodynamic surface with dynamic stall solver
        if isnothing(element.aero) || !(typeof(element.aero.solver) in [BLi])
            continue
        end
        # Set flag
        hasDSelements = true
        # Unpack variables
        @unpack stallOnsetRatioPrev = element.aero.BLiCompVars
        @unpack stallOnsetRatio = element.aero.BLiFlow
        # Compute boundary
        boundary[e] = (abs(stallOnsetRatioPrev)-1) * (abs(stallOnsetRatio)-1)
    end

    # Get maximum
    boundaryMax = hasDSelements ? maximum(boundary) : 0.0

    return boundaryMax
end


# Solves the current time step
function solve_time_step!(problem::Problem)

    if problem.getLinearSolution
        solve_linear_system!(problem)
        problem.systemSolver.convergedFinalSolution = true
    else
        if typeof(problem.systemSolver) == NewtonRaphson
            solve_NewtonRaphson!(problem)
        end
    end

end


# Saves the solution at the current time step
function save_time_step_data!(problem::Problem,timeNow::Real)

    @unpack x,savedTimeVector,xOverTime,elementalStatesOverTime,nodalStatesOverTime,compElementalStatesOverTime,elementalStatesRatesOverTime,compElementalStatesRatesOverTime,aeroVariablesOverTime,model = problem
    @unpack elements = model

    # Add current time
    push!(savedTimeVector,timeNow)

    # Add curent system states 
    push!(xOverTime,x)

    # Add current elemental states 
    currentElementalStates = Vector{ElementalStates}()
    for element in elements
        push!(currentElementalStates,deepcopy(element.states))
    end
    push!(elementalStatesOverTime,currentElementalStates)

    # Add current nodal states 
    currentNodalStates = Vector{NodalStates}()
    for element in elements
        push!(currentNodalStates,deepcopy(element.nodalStates))
    end
    push!(nodalStatesOverTime,currentNodalStates)

    # Add current complementary elemental states 
    currentComplementaryElementalStates = Vector{ComplementaryElementalStates}()
    for element in elements
        push!(currentComplementaryElementalStates,deepcopy(element.compStates))
    end
    push!(compElementalStatesOverTime,currentComplementaryElementalStates)

    # Add current states' rates
    currentStatesRates = Vector{ElementalStatesRates}()
    for element in elements
        push!(currentStatesRates,deepcopy(element.statesRates))
    end
    push!(elementalStatesRatesOverTime,currentStatesRates)

    # Add current complementary elemental states' rates
    currentComplementaryElementalStatesRates = Vector{ComplementaryElementalStatesRates}()
    for element in elements
        push!(currentComplementaryElementalStatesRates,deepcopy(element.compStatesRates))
    end
    push!(compElementalStatesRatesOverTime,currentComplementaryElementalStatesRates)

    # Add current aerodynamic variables
    currentAeroVariables = Vector{Union{Nothing,AeroVariables}}()
    for element in elements
        # Skip elements without aero
        if isnothing(element.aero)
            push!(currentAeroVariables,nothing)
            continue
        end
        push!(currentAeroVariables,AeroVariables(deepcopy(element.aero.flowParameters),deepcopy(element.aero.flowAnglesAndRates),deepcopy(element.aero.flowVelocitiesAndRates),deepcopy(element.aero.aeroCoefficients),deepcopy(element.aero.BLiKin),deepcopy(element.aero.BLiFlow),deepcopy(element.aero.BLoFlow)))
    end
    push!(aeroVariablesOverTime,currentAeroVariables)

    @pack! problem = savedTimeVector,xOverTime,elementalStatesOverTime,nodalStatesOverTime,compElementalStatesOverTime,elementalStatesRatesOverTime,compElementalStatesRatesOverTime,aeroVariablesOverTime

end