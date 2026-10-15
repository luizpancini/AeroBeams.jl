## PROBLEM COMPOSITE TYPES
# ------------------------------------------------------------------------------
# TrimWorkspace composite type
# Fields
# - `J`: dense Jacobian matrix
# - `JtJ`: J'*J matrix
# - `Jtr`: J'*residual vector
# - `J_aug`: augmented Jacobian matrix
# - `r_aug`: augmented residual vector
#
struct TrimWorkspace

    J::Matrix{Float64}
    JtJ::Matrix{Float64}
    Jtr::Vector{Float64}
    J_aug::Matrix{Float64}
    r_aug::Vector{Float64}
    
end
TrimWorkspace(nEqs,nStates) = TrimWorkspace(
    Matrix{Float64}(undef,nEqs,nStates),
    Matrix{Float64}(undef,nStates,nStates),
    Vector{Float64}(undef,nStates),
    zeros(Float64,nEqs+nStates,nStates),
    zeros(Float64,nEqs+nStates),
)


# TrimSolution composite type
# Fields
# - `trimLoads`: vector of trim loads
# - `trimδ`: vector of trim control surface deflections
#
struct TrimSolution

    trimLoads::Vector{Float64}
    trimδ::Vector{Float64}

end


# ShiftAndInvert composite type for generalized eigenvalue problem
# Fields
# - `J`: factorization of the Jacobian matrix
# - `I`: inertia matrix
# - `temp`: temporary array for eigenvalue solver
#
struct ShiftAndInvert{TJ,TI,TT}

    J::TJ
    I::TI
    temp::TT

end


# ModeShape composite type
# Fields
# - `frequency`: mode frequency
# - `damping`: mode damping
# - `state`: vector of modal states (with length equal to the number of elements in the model)
#
struct ModeShape
    
    frequency::Float64
    damping::Float64
    state::Vector{ModalState}

end
export ModeShape


# EigenSolution composite type
# Fields
# - `λ`: eigenvalues
# - `V`: eigenvectors
# - `modeShapes`: mode shapes
#
struct EigenSolution

    λ::Vector{ComplexF64}
    V::Matrix{ComplexF64}
    
    modeShapes::Vector{ModeShape}

end
export EigenSolution
EigenSolution() = EigenSolution(
    Vector{ComplexF64}(), Matrix{ComplexF64}(undef, 0, 0), Vector{ModeShape}()
)


# SystemArrays composite type
# Fields
# - `x`: states array
# - `Δx`: states increment array
# - `residual`: residual array
# - `jacobian`: sparse Jacobian matrix
# - `inertia`: (sparse) inertia matrix
# - `luCache`: LU factorization cache for the Jacobian matrix
#
mutable struct SystemArrays

    x::Vector{Float64}
    Δx::Vector{Float64}
    residual::Vector{Float64}
    jacobian::SparseMatrixCSC{Float64,Int}
    inertia::SparseMatrixCSC{Float64,Int}
    luCache::Union{Nothing,SparseArrays.UMFPACK.UmfpackLU{Float64,Int}}
    
end
Base.copy(x::SystemArrays) = copy_fields(x)


# LoadStepHistory composite type
# Fields
# - `σ`: array of load factors at load steps
# - `states`: array of element states over load steps
# - `hacStates`: array of hinge axis constraint states over load steps
# - `aeroVariables`: array of attached-flow aerodynamic variables over load steps
# - `bliAeroVariables`: array of BLi model aerodynamic variables over load steps
# - `bloAeroVariables`: array of BLo model aerodynamic variables over load steps
#
mutable struct LoadStepHistory

    σ::Vector{Float64}
    states::Vector{Vector{ElementState}}
    hacStates::Vector{Vector{Union{Nothing,HACStateUnion}}}
    aeroVariables::Vector{Vector{Union{Nothing,AttachedAeroVariables{Float64}}}}
    bliAeroVariables::Vector{Vector{Union{Nothing,BLiAeroVariables{Float64}}}}
    bloAeroVariables::Vector{Vector{Union{Nothing,BLoAeroVariables{Float64}}}}

end
LoadStepHistory() = LoadStepHistory(
    Vector{Float64}(), Vector{Vector{ElementState}}(), Vector{Vector{Union{Nothing,HACStateUnion}}}(), Vector{Vector{Union{Nothing,AttachedAeroVariables{Float64}}}}(), Vector{Vector{Union{Nothing,BLiAeroVariables{Float64}}}}(), Vector{Vector{Union{Nothing,BLoAeroVariables{Float64}}}}(),
)


# SolutionHistory composite type
# Fields
# - `index`: current saved time step index
# - `t`: array of saved time instances
# - `states`: array of element states over time
# - `hacStates`: array of hinge axis constraint states over time
# - `aeroVariables`: array of attached-flow aerodynamic variables over time
# - `bliAeroVariables`: array of BLi model aerodynamic variables over time
# - `bloAeroVariables`: array of BLo model aerodynamic variables over time
#
mutable struct SolutionHistory

    index::Int

    t::Vector{Float64}
    states::Vector{Vector{ElementState}}
    hacStates::Vector{Vector{Union{Nothing,HACStateUnion}}}
    aeroVariables::Vector{Vector{Union{Nothing,AttachedAeroVariables{Float64}}}}
    bliAeroVariables::Vector{Vector{Union{Nothing,BLiAeroVariables{Float64}}}}
    bloAeroVariables::Vector{Vector{Union{Nothing,BLoAeroVariables{Float64}}}}

end
SolutionHistory() = SolutionHistory(
    1,
    Vector{Float64}(), Vector{Vector{ElementState}}(), Vector{Vector{Union{Nothing,HACStateUnion}}}(), Vector{Vector{Union{Nothing,AttachedAeroVariables{Float64}}}}(), Vector{Vector{Union{Nothing,BLiAeroVariables{Float64}}}}(), Vector{Vector{Union{Nothing,BLoAeroVariables{Float64}}}}(),
)


# TimeState composite type
# Fields
# - `timeNow`: current time
# - `timeIndex`: current time index
# - `Δt`: current time step
#
mutable struct TimeState

    timeNow::Float64
    timeIndex::Int
    Δt::Float64

end
TimeState(timeVector) = TimeState(
    timeVector[1], 1, timeVector[2]-timeVector[1]
)


# Time composite type
# Fields
# - `initialTime`: initial time
# - `finalTime`: final time
# - `timeVector`: time vector
# - `sizeOfTime`: size of time vector
#
struct Time

    initialTime::Float64
    finalTime::Float64
    timeVector::Vector{Float64}
    sizeOfTime::Int
    
    state::TimeState

end
Time(timeVector) = Time(
    timeVector[1], timeVector[end], timeVector, length(timeVector),
    TimeState(timeVector)
)


# SteadyProblem composite type
# Fields
# - `model`: model
# - `systemSolver`: system solver
# - `systemArrays`: system arrays
# - `loadstepHistory`: loadstep solution history
#
@with_kw mutable struct SteadyProblem

    # Definition fields (inputs to SteadyProblem creation)
    model::Model
    systemSolver::SystemSolver

    # System and history fields
    systemArrays::SystemArrays
    loadstepHistory::LoadStepHistory

end
export SteadyProblem


# TrimProblem composite type
# Fields
# - `model`: model
# - `systemSolver`: system solver
# - `getInertiaMatrix`: flag to compute inertia matrix
# - `systemArrays`: system arrays
# - `trimWorkspace`: workspace with arrays for trim problem
# - `trimSolution`: trim solution
# - `loadstepHistory`: loadstep solution history
#
@with_kw mutable struct TrimProblem

    # Definition fields (inputs to TrimProblem creation)
    model::Model
    systemSolver::NonlinearSolver
    getInertiaMatrix::Bool

    # System and history fields
    systemArrays::SystemArrays
    trimWorkspace::TrimWorkspace
    trimSolution::TrimSolution
    loadstepHistory::LoadStepHistory

end
export TrimProblem


# EigenProblem composite type
# Fields
# - `model`: model
# - `systemSolver`: system solver
# - `eigenmethod`: eigenvalue computation method
# - `nev`: minimum number of eigenvalues to be computed
# - `refTrimProblem`: reference trim problem, containing the steady solution
# - `systemArrays`: system arrays
# - `loadstepHistory`: loadstep solution history
# - `eigensolution`: eigensolution
#
@with_kw mutable struct EigenProblem

    # Definition fields (inputs to EigenProblem creation)
    model::Model
    systemSolver::SystemSolver
    eigenmethod::Symbol
    nev::Int
    refTrimProblem::Union{Nothing,TrimProblem}

    # System and history fields
    systemArrays::SystemArrays
    loadstepHistory::LoadStepHistory
    eigensolution::EigenSolution

end
export EigenProblem


# DynamicProblem composite type
# Fields
# - `model`: model
# - `systemSolver`: system solver
# - `skipInitialStatesUpdate`: flag to skip update of initial states
# - `trackingTimeSteps`: flag to track time step solutions
# - `trackingFrequency`: frequency of time steps in which to track solution
# - `displayProgress`: flag to display progress
# - `displayFrequency`: frequency of time steps in which to display progress
# - `time`: time variables
# - `systemArrays`: system arrays
# - `solutionHistory`: solution history over time
#
@with_kw mutable struct DynamicProblem

    # Definition fields (inputs to DynamicProblem creation)
    model::Model
    systemSolver::SystemSolver
    skipInitialStatesUpdate::Bool
    trackingTimeSteps::Bool
    trackingFrequency::Int
    displayProgress::Bool
    displayFrequency::Int

    # Time variables field
    time::Time

    # System and history fields
    systemArrays::SystemArrays
    solutionHistory::SolutionHistory

end
export DynamicProblem


const Problem = Union{SteadyProblem,TrimProblem,EigenProblem,DynamicProblem}

## PUBLIC API
# -----------------------------------------------------------------------------
"""
    create_SteadyProblem(; kwargs...)

Steady problem constructor

# Keyword arguments
- `model::Model`: model
- `systemSolver::SystemSolver`: system solver
- `x0::AbstractVector{<:Real}`: initial states' guess
"""
function create_SteadyProblem(; 
    model::Model,
    systemSolver::SystemSolver=create_NLsolver(),
    x0::AbstractVector{<:Real}=zeros(0)
    )

    # Validate inputs
    validate_steady_problem_inputs(systemSolver,x0)

    # Reset model state, if applicable
    if isempty(x0)
        reset_model_state!(model)
    end

    # Reset system solver state
    reset_solver_state!(systemSolver)

    # Initialize system arrays with correct size and set initial states
    systemArrays = initialize_system_arrays!(model,systemSolver,x0)

    # Initialize loadstep solution history
    loadstepHistory = LoadStepHistory()

    return SteadyProblem(
        model,
        systemSolver,
        systemArrays,
        loadstepHistory
    )

end
export create_SteadyProblem


"""
    create_TrimProblem(; kwargs...)

Trim problem constructor

# Keyword arguments
- `model::Model`: model
- `systemSolver::NonlinearSolver`: system solver
- `getInertiaMatrix::Bool`: flag to compute inertia matrix
- `x0::AbstractVector{<:Real}`: initial states' guess
"""
function create_TrimProblem(; 
    model::Model,
    systemSolver::NonlinearSolver=create_NLsolver(),
    getInertiaMatrix::Bool=true,
    x0::AbstractVector{<:Real}=zeros(0)
    )

    # Validate inputs
    validate_trim_problem_inputs(systemSolver)

    # Reset model state, if applicable
    if isempty(x0)
        reset_model_state!(model)
    end

    # Reset system solver state
    reset_solver_state!(systemSolver)

    # Initialize system arrays with correct size and set initial states
    systemArrays = initialize_system_arrays!(model,systemSolver,x0,isTrimProblem=true)

    # Initialize trim workspace with correct size
    trimWorkspace = TrimWorkspace(
        systemArrays.residual |> length,
        systemArrays.x |> length
    )

    # Initialize trim solution
    trimSolution = TrimSolution(
        zeros(length(model.assemblyData.nTrimLoads)),
        zeros(length(model.assemblyData.nTrimδ))
    )

    # Initialize loadstep solution history
    loadstepHistory = LoadStepHistory()

    return TrimProblem(
        model,
        systemSolver,
        getInertiaMatrix,
        systemArrays,
        trimWorkspace,
        trimSolution,
        loadstepHistory
    )

end
export create_TrimProblem


"""
    create_EigenProblem(; kwargs...)

Eigen problem constructor

# Keyword arguments
- `model::Model`: model
- `eigenmethod::Symbol`: eigenvalue computation method (default: `:krylovkit`)
- `systemSolver::SystemSolver`: system solver
- `nev::Int`: number of eigenvalues to be computed
- `x0::AbstractVector{<:Real}`: initial states' guess
- `refTrimProblem::Union{Nothing,TrimProblem}`: reference trim problem containing steady solution
"""
function create_EigenProblem(; 
    model::Model,
    eigenmethod::Symbol=:krylovkit,
    nev::Int=model.assemblyData.systemOrder,
    x0::AbstractVector{<:Real}=zeros(0),
    refTrimProblem::Union{Nothing,TrimProblem}=nothing,
    systemSolver::SystemSolver = refTrimProblem isa TrimProblem ? refTrimProblem.systemSolver : create_NLsolver(),
    )

    # Validate inputs
    validate_eigenproblem_inputs(systemSolver,eigenmethod,nev,x0,refTrimProblem)

    # Reset model state, if applicable
    if isempty(x0)
        reset_model_state!(model)
    end

    # Reset system solver state, if applicable
    if isnothing(refTrimProblem)
        reset_solver_state!(systemSolver)
    end
    
    # Initialize system arrays with correct size and set initial states
    systemArrays = initialize_system_arrays!(model,systemSolver,x0,isEigenProblem=true,refTrimProblem=refTrimProblem)

    # Initialize loadstep solution history
    loadstepHistory = isnothing(refTrimProblem) ?
                      LoadStepHistory() :                 
                      deepcopy(refTrimProblem.loadstepHistory)

    # Initialize eigensolution
    eigensolution = EigenSolution()                      

    return EigenProblem(
        model,
        systemSolver,
        eigenmethod,
        nev,
        refTrimProblem,
        systemArrays,
        loadstepHistory,
        eigensolution
    )

end
export create_EigenProblem


"""
    create_DynamicProblem(; kwargs...)

Dynamic problem constructor. Unless `skipInitialStatesUpdate` is set, the initial states and their rates (velocities and accelerations) are made consistent with the equations of motion, the motion of basis A, the BCs and, when applicable, the initial conditions of the beams (see [`add_initial_conditions_to_beam!`](@ref add_initial_conditions_to_beam!)), without any user-defined tuning.

# Keyword arguments
- `model::Model`: model
- `systemSolver::SystemSolver`: system solver
- `skipInitialStatesUpdate::Bool`: flag to skip update of initial states
- `trackingTimeSteps::Bool`: flag to track time step solutions
- `trackingFrequency::Int`: frequency of time steps in which to track solution
- `displayProgress::Bool`: flag to display progress
- `displayFrequency::Union{Nothing,Int}`: frequency of time steps in which to display progress
- `initialTime::Real`: initial time
- `Δt::Union{Nothing,Real}`: time step
- `finalTime::Union{Nothing,Real}`: final time
- `timeVector::Union{Nothing,AbstractVector{<:Real}}`: time vector
- `x0::AbstractVector{<:Real}`: initial states' guess
- `refDynamicProblem::Union{DynamicProblem,Nothing}`: reference problem for continuation of dynamic solution
"""
function create_DynamicProblem(; 
    model::Model,
    systemSolver::SystemSolver=create_NLsolver(),
    trackingTimeSteps::Bool=true,
    trackingFrequency::Int=1,
    displayProgress::Bool=true,
    displayFrequency::Union{Nothing,Int}=nothing,
    initialTime::Real=0.0,
    Δt::Union{Nothing,Real}=nothing,
    finalTime::Union{Nothing,Real}=nothing,
    timeVector::Union{Nothing,AbstractVector{<:Real}}=nothing,
    x0::AbstractVector{<:Real}=zeros(0),
    refDynamicProblem::Union{DynamicProblem,Nothing}=nothing,
    skipInitialStatesUpdate::Bool=x0!=zeros(0),
    )

    # Validate inputs
    validate_dynamic_problem_inputs(systemSolver,trackingFrequency,displayFrequency,displayProgress,x0,skipInitialStatesUpdate)

    # Reset model state, if applicable
    if isempty(x0)
        reset_model_state!(model)
    end

    # Reset system solver state
    reset_solver_state!(systemSolver)

    # Validate and initialize time variables
    time = initialize_time_variables(initialTime,finalTime,Δt,timeVector)

    # Initialize system arrays
    systemArrays = initialize_system_arrays!(model,systemSolver,x0,refDynamicProblem=refDynamicProblem)

    # Update initial states from x0, if provided
    if !isempty(x0)
        update_initial_states_from_x0!(model,time,systemSolver,x0)
    end

    # Initialize solution history
    solutionHistory = initialize_solution_history(time,model.assemblyData.elements,model.hingeAxisConstraints,trackingTimeSteps,trackingFrequency)

    # Set display frequency if not input
    if isnothing(displayFrequency)
        displayFrequency = ceil(time.sizeOfTime/100)
    end

    return DynamicProblem(
        model,
        systemSolver,
        skipInitialStatesUpdate,
        trackingTimeSteps,
        trackingFrequency,
        displayProgress,
        displayFrequency,
        time,
        systemArrays,
        solutionHistory
    )

end
export create_DynamicProblem



"""
    solve!(problem)

Solves a problem.

# Arguments
- problem::Problem: the problem instance to be solved.
"""
function solve!(problem::Problem)

    # Precompute prescribed motion of basis A
    precompute_motion_basis_A!(problem)

    # Precompute prescribed BCs
    precompute_BCs_history!(problem)

    # Precompute prescribed distributed loads
    precompute_distributed_loads!(problem)

    # Precompute prescribed flap deflection values
    precompute_flap_deflection!(problem)

    # Update motion of basis A
    update_motion_basis_A!(problem,1)

    # Update boundary conditions
    update_BCs_state!(problem,1)

    # Solve according to problem type
    if problem isa SteadyProblem || problem isa TrimProblem
        solve_steady!(problem)
    elseif problem isa EigenProblem
        solve_steady!(problem)
        solve_eigen!(problem)    
    elseif problem isa DynamicProblem
        solve_dynamic!(problem)
    end

end
export solve!


"""
    solve_eigen!(problem)

Solves an eigenvalue problem whose steady-state solution has already been obtained.

# Arguments
# problem::EigenProblem: the eigenvalue problem instance to be solved.
"""
function solve_eigen!(problem::EigenProblem)

    @unpack eigenmethod,nev = problem
    @unpack luCache,inertia,jacobian = problem.systemArrays
    @unpack converged = problem.systemSolver.state

    # Update state at the steady solution, if the eigenproblem was initialized from a reference trim problem
    if !isnothing(problem.refTrimProblem)
        f!(problem)
        foreach_element(problem.model.assemblyData.elementGroups) do element
            element_rotation_derivatives!(element,0.0,false)
        end
    end
       
    # Check number of available eigenvalues and update nev if necessary
    nevAvailable = size(inertia,1)
    if nevAvailable < nev
        @info("Only $nevAvailable (instead of $nev requested) eigenvalues can be calculated, due to limited number of elements in assembly")
        nev = nevAvailable
    end

    # Small shift for singular Jacobians
    σ = converged ?
        0.0 : 
        1e-6*opnorm(jacobian,Inf)

    # Set LU factorization of the shifted Jacobian, if applicable
    luCache = (converged && !isnothing(luCache)) ?
              luCache : 
              lu(jacobian - σ*I)

    # Solve eigenproblem to get eigenvectors and corresponding eigenvalues in ascending order of frequency
    λ, V = eigs_generalized(luCache,inertia; method=eigenmethod,nev=nev,which=:LM,σ=σ,sortby=x->abs(imag(x)),restarts=nevAvailable)

    # Re-solve with a small shift if the unshifted operator is ill-scaled (see `needs_shift`)
    if converged && eigenmethod != :fulldecomp
        σ2 = needs_shift(λ,σ)
        if σ2 > 0
            try
                luShifted = lu(jacobian + σ2*inertia)
                λ, V = eigs_generalized(luShifted,inertia; method=eigenmethod,nev=nev,which=:LM,σ=σ2,sortby=x->abs(imag(x)),restarts=nevAvailable)
            catch err
                @warn "Could not re-solve the eigenproblem with a shift of $σ2 ($(sprint(showerror, err))): keeping the unshifted solution, whose damping may be contaminated by round-off errors"
            end
        end
    end

    # Get mode shapes
    modeShapes = mode_shapes(problem,λ,V)

    # Set eigensolution
    setfield!(problem, :eigensolution, EigenSolution(λ, V, modeShapes))

end
export solve_eigen!



## PRIVATE API
# -----------------------------------------------------------------------------
# Validate inputs for steady problem
function validate_steady_problem_inputs(systemSolver,x0)

    if systemSolver isa LinearSolver
        @assert isempty(x0) "initial guess `x0` must be empty for linear solvers"
    end

end

# Validate inputs for trim problem
function validate_trim_problem_inputs(systemSolver)

    if systemSolver isa NLsolver
        @assert systemSolver.props.method == :newton "trim problem with NLsolver requires the :newton method"
        # @assert systemSolver.props.linesearch == false "trim problem with NLsolver requires linesearch == false"
    end

end

# Validate inputs for eigenproblem
function validate_eigenproblem_inputs(systemSolver,eigenmethod,nev,x0,refTrimProblem)

    if systemSolver isa LinearSolver
        @assert isempty(x0) "initial guess `x0` must be empty for linear solvers"
    end

    @assert eigenmethod in [:arnoldi, :arpack, :krylovkit, :fulldecomp] "unsupported eigenvalue computation `eigenmethod`, select one of :arnoldi, :arpack, :krylovkit, :fulldecomp"

    @assert nev > 0 "number of eigenvalues `nev` must be positive"

    @assert isempty(x0) || isnothing(refTrimProblem) "cannot input `x0` and `refTrimProblem` at the same time"

end


# Validate inputs for dynamic problem
function validate_dynamic_problem_inputs(systemSolver,trackingFrequency,displayFrequency,displayProgress,x0,skipInitialStatesUpdate)

    if systemSolver isa LinearSolver
        @assert isempty(x0) "initial guess `x0` must be empty for linear solvers"
    end

    @assert trackingFrequency > 0 "`trackingFrequency` must be positive"

    if !isnothing(displayFrequency) && displayProgress
        @assert displayFrequency > 0 "`displayFrequency` must be positive"
    end

    @assert !(!isempty(x0) && !skipInitialStatesUpdate) "cannot input `x0` and set `skipInitialStatesUpdate=false` at the same time"

end


# Initializes the system arrays to the correct size and sets the initial elemental and nodal states, or sets the initial states' guess if input
function initialize_system_arrays!(model,systemSolver,x0; isTrimProblem=false,isEigenProblem=false,refTrimProblem=nothing,refDynamicProblem=nothing)

    # Unpack model data
    @unpack elements,specialNodes,forceScaling,systemOrder,nTrimVariables = model.assemblyData
    @unpack hingeAxisConstraints = model

    # Set number of equations and states of the system
    nEqs = systemOrder
    nStates = isTrimProblem ?
              systemOrder + nTrimVariables :
              systemOrder

    # Set flags
    inertiaWillBeUsed = isTrimProblem || isEigenProblem
    resetMatrices = isnothing(refTrimProblem) && isnothing(refDynamicProblem)

    # Set initial system arrays
    systemArrays = if !isempty(x0)
        initialize_system_arrays_from_x0(x0,nEqs,nStates,elements,specialNodes,hingeAxisConstraints,inertiaWillBeUsed)
    elseif !isnothing(refTrimProblem)
        initialize_system_arrays_from_refTrimProblem(refTrimProblem,nEqs,nStates)
    elseif !isnothing(refDynamicProblem)
        initialize_system_arrays_from_refDynamicProblem(refDynamicProblem,nEqs,nStates)
    else
        initialize_default_system_arrays(systemSolver,forceScaling,nEqs,nStates,elements,specialNodes,hingeAxisConstraints,inertiaWillBeUsed)
    end

    # Set Jacobian indices' maps for each element and special node
    set_jacobian_indices_maps!(systemArrays.jacobian,elements,specialNodes,hingeAxisConstraints,resetMatrices)

    # Set inertia matrix indices' maps for each element, if applicable
    if inertiaWillBeUsed
        set_inertia_indices_maps!(systemArrays.inertia,elements,resetMatrices)
    end

    return systemArrays

end


# Initializes the system arrays from the input initial states' guess `x0`
function initialize_system_arrays_from_x0(x0,nEqs,nStates,elements,specialNodes,hingeAxisConstraints,inertiaWillBeUsed)

    # Validate
    @assert length(x0) == nStates "x0 must have length equal to the number of states of the model (length(x0) = $(length(x0)), nStates = $nStates)"

    # Initialize
    x = x0
    Δx = zeros(nStates)
    residual = zeros(nEqs)
    jacobian = initialize_jacobian_sparsity_pattern(nEqs,nStates,elements,specialNodes,hingeAxisConstraints)
    inertia = inertiaWillBeUsed ? 
              initialize_inertia_sparsity_pattern(nEqs,elements) : 
              spzeros(nEqs,nEqs)
    luCache = nothing

    return SystemArrays(x,Δx,residual,jacobian,inertia,luCache)

end


# Initializes the system arrays from the reference trim problem
function initialize_system_arrays_from_refTrimProblem(refTrimProblem,nEqs,nStates)

    # Sizes from reference trim problem
    nEqsTrimProblem,nStatesTrimProblem = size(refTrimProblem.systemArrays.jacobian)
    nTrimVariables = refTrimProblem.model.assemblyData.nTrimVariables

    # Validate
    @assert nEqs == nEqsTrimProblem "number of equations of the current model is different from that of the reference trim problem"
    @assert nStates == (nStatesTrimProblem-nTrimVariables) "number of states of the current model is different from that of the reference trim problem"

    # Copy system arrays from reference trim problem
    x = copy(refTrimProblem.systemArrays.x[1:nEqs])
    Δx = zeros(nEqs)
    residual = copy(refTrimProblem.systemArrays.residual)
    jacobian = copy(refTrimProblem.systemArrays.jacobian[1:nEqs,1:nEqs])
    inertia = copy(refTrimProblem.systemArrays.inertia)
    local luCache
    try
        luCache = lu(jacobian)
    catch
        luCache = lu(jacobian - 1e-6*opnorm(jacobian,Inf)*I)
    end

    return SystemArrays(x,Δx,residual,jacobian,inertia,luCache)

end


# Initializes the system arrays from the reference dynamic problem
function initialize_system_arrays_from_refDynamicProblem(refDynamicProblem,nEqs,nStates)

    # Sizes from reference dynamic problem
    nEqsDynamicProblem,nStatesDynamicProblem = size(refDynamicProblem.jacobian)

    # Validate
    @assert nEqs == nEqsDynamicProblem "number of equations of the current model is different from that of the reference dynamic problem"
    @assert nStates == nStatesDynamicProblem "number of states of the current model is different from that of the reference dynamic problem"

    # Copy system arrays from reference dynamic problem
    return copy(refDynamicProblem.systemArrays)

end


# Initializes the system arrays to the correct size and sets the initial elemental and nodal states
function initialize_default_system_arrays(systemSolver,forceScaling,nEqs,nStates,elements,specialNodes,hingeAxisConstraints,inertiaWillBeUsed)

    # Initialize
    x = zeros(nStates)
    Δx = zeros(nStates)
    residual = zeros(nEqs)
    jacobian = initialize_jacobian_sparsity_pattern(nEqs,nStates,elements,specialNodes,hingeAxisConstraints)
    inertia = inertiaWillBeUsed ? 
              initialize_inertia_sparsity_pattern(nEqs,elements) : 
              spzeros(nEqs,nEqs)
    luCache = nothing

    # Loop over elements and assign initial states
    for element in elements
        @unpack DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
        @unpack u,p,F,M,V,Ω,χ = element.state.midpoint
        x[DOF_u] .= u
        x[DOF_p] .= p
        x[DOF_F] .= F/forceScaling
        x[DOF_M] .= M/forceScaling
        x[DOF_V] .= V
        x[DOF_Ω] .= Ω
        x[DOF_χ] .= χ
    end

    # Loop over special nodes and assign initial states
    for specialNode in specialNodes
        @unpack BCs = specialNode
        @unpack uIsPrescribed,pIsPrescribed = specialNode.flags
        @unpack u,p,F,M = specialNode.state
        @unpack DOF_uF,DOF_pM,DOF_T = specialNode.assemblyData
        for i=1:3
            x[DOF_uF[i]] = uIsPrescribed[i] ?
                           F[i]/forceScaling :
                           u[i]
            x[DOF_pM[i]] = pIsPrescribed[i] ?
                           M[i]/forceScaling :
                           p[i]
        end
    end

    # Scale states by initial load factor
    x *= systemSolver.loadStepState.σ

    return SystemArrays(x,Δx,residual,jacobian,inertia,luCache)

end


# Initializes the Jacobian sparsity pattern
function initialize_jacobian_sparsity_pattern(nEqs,nStates,elements,specialNodes,hingeAxisConstraints)

    # Initialize row and column indices for the sparse Jacobian matrix
    I = Int[]
    J = Int[]

    # Give size hint to reduce repeated memory allocations with push!()
    sizehint!(I, 37 * 9 * length(elements) + 8 * 9 * length(specialNodes))
    sizehint!(J, 37 * 9 * length(elements) + 8 * 9 * length(specialNodes))

    # Helper function
    add_block!(rows, cols) = for j in cols, i in rows
        push!(I, i); push!(J, j)
    end

    # Loop over elements
    for element in elements
        # Unpack assembly indices
        @unpack eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,eqs_Fχ,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
        @unpack aero = element
        DOF_δ = !isnothing(aero.surface) && aero.geometry.isFlapped ? 
                aero.surface.assemblyData.DOF_δ : 
                Int[]
        # Set indices of node 1
        add_block!(eqs_Fu1, DOF_p); add_block!(eqs_Fu1, DOF_F); add_block!(eqs_Fu1, DOF_V); add_block!(eqs_Fu1, DOF_Ω); add_block!(eqs_Fu1, DOF_χ); add_block!(eqs_Fu1, DOF_δ)
        add_block!(eqs_Fp1, DOF_p); add_block!(eqs_Fp1, DOF_F); add_block!(eqs_Fp1, DOF_M); add_block!(eqs_Fp1, DOF_V); add_block!(eqs_Fp1, DOF_Ω); add_block!(eqs_Fp1, DOF_χ); add_block!(eqs_Fp1, DOF_δ)
        add_block!(eqs_FF1, DOF_u); add_block!(eqs_FF1, DOF_p); add_block!(eqs_FF1, DOF_F); add_block!(eqs_FF1, DOF_M)
        add_block!(eqs_FM1, DOF_p); add_block!(eqs_FM1, DOF_F); add_block!(eqs_FM1, DOF_M)
        # Set indices of node 2
        add_block!(eqs_Fu2, DOF_p); add_block!(eqs_Fu2, DOF_F); add_block!(eqs_Fu2, DOF_V); add_block!(eqs_Fu2, DOF_Ω); add_block!(eqs_Fu2, DOF_χ); add_block!(eqs_Fu2, DOF_δ)
        add_block!(eqs_Fp2, DOF_p); add_block!(eqs_Fp2, DOF_F); add_block!(eqs_Fp2, DOF_M); add_block!(eqs_Fp2, DOF_V); add_block!(eqs_Fp2, DOF_Ω); add_block!(eqs_Fp2, DOF_χ); add_block!(eqs_Fp2, DOF_δ)
        add_block!(eqs_FF2, DOF_u); add_block!(eqs_FF2, DOF_p); add_block!(eqs_FF2, DOF_F); add_block!(eqs_FF2, DOF_M)
        add_block!(eqs_FM2, DOF_p); add_block!(eqs_FM2, DOF_F); add_block!(eqs_FM2, DOF_M)
        # Set indices of the element midpoint
        add_block!(eqs_FV, DOF_u); add_block!(eqs_FV, DOF_p); add_block!(eqs_FV, DOF_V)
        add_block!(eqs_FΩ, DOF_p); add_block!(eqs_FΩ, DOF_Ω)
        add_block!(eqs_Fχ, DOF_p); add_block!(eqs_Fχ, DOF_V); add_block!(eqs_Fχ, DOF_Ω); add_block!(eqs_Fχ, DOF_χ)
    end

    # Loop over special nodes
    for specialNode in specialNodes
        # Unpack special node assembly data
        @unpack connectedElementsGlobalIDs,doubleAttachmentSpringsData,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,DOF_uF,DOF_pM,DOF_T = specialNode.assemblyData
        # Set indices of equilibrium equations
        add_block!(eqs_Fu, DOF_uF); add_block!(eqs_Fu, DOF_pM);
        add_block!(eqs_Fp, DOF_uF); add_block!(eqs_Fp, DOF_pM);
        # Set indices of trim forces equations, if there are trim forces
        DOF_TF = DOF_T[1:3]
        if any(!iszero(DOF_TF))
            nzT = findall(!iszero,DOF_TF)
            add_block!(eqs_Fu, DOF_TF[nzT])
        end
        # Set indices of trim moment equations, if there are trim moments
        DOF_TM = DOF_T[4:6]
        if any(!iszero(DOF_TM))
            nzT = findall(!iszero,DOF_TM)
            add_block!(eqs_Fp, DOF_TM[nzT])
        end
        # Set indices of compatibility equations
        for e in eachindex(connectedElementsGlobalIDs)
            add_block!(eqs_FF[e], DOF_uF); add_block!(eqs_FF[e], DOF_pM)
            add_block!(eqs_FM[e], DOF_uF); add_block!(eqs_FM[e], DOF_pM)
        end
        # Initialize counter for double attachment springs
        c = 0
        # Unpack and loop over special node's springs
        @unpack springs = specialNode
        for spring in springs
            @unpack hasDoubleAttachment = spring.assemblyData
            # Double attachment springs
            if hasDoubleAttachment
                # Update counter
                c += 1
                # Unpack double attachment spring data for this spring
                @unpack otherNodeDOF_uF,otherNodeDOF_pM = doubleAttachmentSpringsData[c]
                # Set indices of equilibrium equations
                add_block!(eqs_Fu, otherNodeDOF_uF); add_block!(eqs_Fu, otherNodeDOF_pM);
                add_block!(eqs_Fp, otherNodeDOF_uF); add_block!(eqs_Fp, otherNodeDOF_pM);
            end
        end
    end

    # Loop hinge axis constraints
    for hac in hingeAxisConstraints
        # Unpack hinge axis constraint assembly data
        @unpack eqs_Fp,eqs_λ,DOF_pM,DOF_pS,DOF_λ = hac.assemblyData
        # Set indices of equilibrium equations
        add_block!(eqs_Fp, DOF_pM); add_block!(eqs_Fp, DOF_pS); add_block!(eqs_Fp, DOF_λ);
        # Set indices of the Lagrange multiplier equations
        add_block!(eqs_λ, DOF_pM); add_block!(eqs_λ, DOF_pS);
    end

    # Set dummy values for the Jacobian
    V = eps()*ones(Float64, length(I))

    # Initialize sparse Jacobian matrix
    jacobian = sparse(I, J, V, nEqs, nStates, +)

    return jacobian

end


# Initializes the inertia matrix sparsity pattern
function initialize_inertia_sparsity_pattern(nEqs,elements)

    # Initialize row and column indices for the sparse inertia matrix
    I = Int[]
    J = Int[]

    # Give size hint to reduce repeated memory allocations with push!()
    sizehint!(I, 14 * 9 * length(elements))
    sizehint!(J, 14 * 9 * length(elements))

    # Helper function
    add_block!(rows, cols) = for j in cols, i in rows
        push!(I, i); push!(J, j)
    end

    # Loop over elements
    for element in elements
        # Unpack assembly indices
        @unpack eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FV,eqs_FΩ,eqs_Fχ,DOF_u,DOF_p,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
        # Set indices of node 1
        add_block!(eqs_Fu1, DOF_p); add_block!(eqs_Fu1, DOF_V); add_block!(eqs_Fu1, DOF_Ω);
        add_block!(eqs_Fp1, DOF_p); add_block!(eqs_Fp1, DOF_V); add_block!(eqs_Fp1, DOF_Ω);
        # Set indices of node 2
        add_block!(eqs_Fu2, DOF_p); add_block!(eqs_Fu2, DOF_V); add_block!(eqs_Fu2, DOF_Ω);
        add_block!(eqs_Fp2, DOF_p); add_block!(eqs_Fp2, DOF_V); add_block!(eqs_Fp2, DOF_Ω);
        # Set indices of the element midpoint
        add_block!(eqs_FV, DOF_u);
        add_block!(eqs_FΩ, DOF_p);
        add_block!(eqs_Fχ, DOF_V); add_block!(eqs_Fχ, DOF_Ω); add_block!(eqs_Fχ, DOF_χ)
    end

    # Set dummy values for the inertia
    V = eps()*ones(Float64, length(I))

    # Initialize sparse inertia matrix
    inertia = sparse(I, J, V, nEqs, nEqs, +)

    return inertia

end


# Updates the initial states of the elements, special nodes, and hinge axis constraints when the state vector x0 is given
function update_initial_states_from_x0!(model,time,systemSolver,x0)

    @unpack σ = systemSolver.loadStepState
    @unpack timeIndex,timeNow = time.state

    update_states!(model,x0,σ,timeIndex,timeNow,Inf,true)

end


# Validates and initializes the time variables
function initialize_time_variables(initialTime,finalTime,Δt,timeVector)

    # Validate inputs and set correct types
    @assert !isnothing(Δt) && !isnothing(finalTime) || !isnothing(timeVector) "either Δt and finalTime, or timeVector must be input"
    initialTime = Float64(initialTime)
    if !isnothing(Δt)
        @assert Δt > 0 "Δt must be positive"
        Δt = Float64(Δt)
    end
    if !isnothing(finalTime)
        @assert finalTime >= initialTime + Δt "finalTime must be greater than or equal to initialTime + Δt"
        finalTime = Float64(finalTime)
    end
    if !isnothing(timeVector)
        @assert length(timeVector) > 1 "timeVector must have length greater than 1"
        @assert timeVector == sort!(copy(timeVector)) "timeVector must be monotonically increasing"
        timeVector = Float64.(timeVector)
        initialTime = timeVector[1]
        Δt = timeVector[2] - timeVector[1]
        finalTime = timeVector[end]
    else
        timeVector = collect(initialTime:Δt:finalTime)
    end

    # Initialize Time object
    time = Time(timeVector)

    return time

end


# Initializes the solution history for a dynamic problem
function initialize_solution_history(time,elements,hingeAxisConstraints,trackingTimeSteps,trackingFrequency)

    # Get total number of elements
    nElementsTotal = length(elements)

    # Initialize solution history
    solutionHistory = SolutionHistory()

    # If tracking time steps, preallocate arrays for saved time and states over time
    if trackingTimeSteps
        # Indices of time steps that will be saved and their size
        savedTimeIndices = collect(1:trackingFrequency:time.sizeOfTime)
        nSavedSteps = length(savedTimeIndices)
        # Set array of saved time instances
        setfield!(solutionHistory, :t, time.timeVector[savedTimeIndices])
        # Preallocate arrays for states and aerodynamic variables over time
        states = Vector{Vector{ElementState}}(undef, nSavedSteps)
        hacStates = Vector{Vector{Union{Nothing,HACStateUnion}}}(undef, nSavedSteps)
        aeroVariables = Vector{Vector{Union{Nothing,AttachedAeroVariables{Float64}}}}(undef, nSavedSteps)
        bliAeroVariables = Vector{Vector{Union{Nothing,BLiAeroVariables{Float64}}}}(undef, nSavedSteps)
        bloAeroVariables = Vector{Vector{Union{Nothing,BLoAeroVariables{Float64}}}}(undef, nSavedSteps)
        for k in 1:nSavedSteps
            states[k] = Vector{ElementState}(undef, nElementsTotal)
            aeroVariables[k] = Vector{Union{Nothing,AttachedAeroVariables{Float64}}}(nothing, nElementsTotal)
            hacStates[k] = Vector{Union{Nothing,HACStateUnion}}(nothing, length(hingeAxisConstraints))
            bliAeroVariables[k] = Vector{Union{Nothing,BLiAeroVariables{Float64}}}(nothing, nElementsTotal)
            bloAeroVariables[k] = Vector{Union{Nothing,BLoAeroVariables{Float64}}}(nothing, nElementsTotal)
            for (e, element) in enumerate(elements)
                states[k][e] = ElementNaNState(element,NaN)
                if !isnothing(element.aero.surface)
                    aeroVariables[k][e] = AttachedNaNAeroVariables()
                    bliAeroVariables[k][e] = copy(element.aero.stateCache.bliVars)
                    bloAeroVariables[k][e] = copy(element.aero.stateCache.bloVars)
                end
            end
            for (h, hac) in enumerate(hingeAxisConstraints)
                hacStates[k][h] = copy(hac.state)
            end
        end
        setfield!(solutionHistory, :states, states)
        setfield!(solutionHistory, :hacStates, hacStates)
        setfield!(solutionHistory, :aeroVariables, aeroVariables)
        setfield!(solutionHistory, :bliAeroVariables, bliAeroVariables)
        setfield!(solutionHistory, :bloAeroVariables, bloAeroVariables)
    end

    return solutionHistory

end


# Sets the Jacobian indices' maps for each element, special node and hinge axis constraint
function set_jacobian_indices_maps!(jacobian,elements,specialNodes,hingeAxisConstraints,resetJacobian)

    # Helper function
    bim(rows, cols) = build_index_map(jacobian, rows, cols)

    # Loop over elements
    for element in elements
        # Unpack
        @unpack eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_Fχ,eqs_FΩ,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
        @unpack aero = element
        DOF_δ = !isnothing(aero.surface) && aero.geometry.isFlapped ? 
                aero.surface.assemblyData.DOF_δ : 
                Int[]
        # Get number of aerodynamic and trim flap states
        N = length(DOF_χ)
        Nδ = length(DOF_δ)
        # Set Jacobian maps
        setfield!(element, :jacobianMaps, 
        JacobianIndicesMaps{3*N,3*Nδ,N*N}(
            bim(eqs_Fu1,DOF_p), bim(eqs_Fu1,DOF_F), bim(eqs_Fu1,DOF_V), bim(eqs_Fu1,DOF_Ω), bim(eqs_Fu1,DOF_χ), bim(eqs_Fu1,DOF_δ),
            bim(eqs_Fp1,DOF_p), bim(eqs_Fp1,DOF_F), bim(eqs_Fp1,DOF_M), bim(eqs_Fp1,DOF_V), bim(eqs_Fp1,DOF_Ω), bim(eqs_Fp1,DOF_χ), bim(eqs_Fp1,DOF_δ),
            bim(eqs_FF1,DOF_u), bim(eqs_FF1,DOF_p), bim(eqs_FF1,DOF_F), bim(eqs_FF1,DOF_M),
            bim(eqs_FM1,DOF_p), bim(eqs_FM1,DOF_F), bim(eqs_FM1,DOF_M),
            bim(eqs_Fu2,DOF_p), bim(eqs_Fu2,DOF_F), bim(eqs_Fu2,DOF_V), bim(eqs_Fu2,DOF_Ω), bim(eqs_Fu2,DOF_χ), bim(eqs_Fu2,DOF_δ),
            bim(eqs_Fp2,DOF_p), bim(eqs_Fp2,DOF_F), bim(eqs_Fp2,DOF_M), bim(eqs_Fp2,DOF_V), bim(eqs_Fp2,DOF_Ω), bim(eqs_Fp2,DOF_χ), bim(eqs_Fp2,DOF_δ),
            bim(eqs_FF2,DOF_u), bim(eqs_FF2,DOF_p), bim(eqs_FF2,DOF_F), bim(eqs_FF2,DOF_M),
            bim(eqs_FM2,DOF_p), bim(eqs_FM2,DOF_F), bim(eqs_FM2,DOF_M),
            bim(eqs_FV,DOF_u), bim(eqs_FV,DOF_p), bim(eqs_FV,DOF_V),
            bim(eqs_FΩ,DOF_p), bim(eqs_FΩ,DOF_Ω),
            bim(eqs_Fχ,DOF_p), bim(eqs_Fχ,DOF_V), bim(eqs_Fχ,DOF_Ω), bim(eqs_Fχ,DOF_χ)
            )
        )
    end

    # Loop over special nodes
    for specialNode in specialNodes
        # Unpack special node data
        @unpack springs = specialNode
        @unpack connectedElementsGlobalIDs,doubleAttachmentSpringsData,eqs_Fu,eqs_Fp,DOF_uF,DOF_pM = specialNode.assemblyData
        # Initialize Jacobian maps for double attachment springs
        mapsDoubleSpring_Fu = Vector{SVector{9,Int}}()
        mapsDoubleSpring_Fp = Vector{SVector{9,Int}}()
        # Initialize counter for double attachment springs
        c = 0
        # Loop over special node's springs
        for spring in springs
            # Double attachment springs
            if spring.assemblyData.hasDoubleAttachment
                # Update counter
                c += 1
                # Unpack double attachment spring data for this spring
                @unpack otherNodeDOF_uF,otherNodeDOF_pM = doubleAttachmentSpringsData[c]
                # Set Jacobian indices maps for this double attachment spring
                push!(mapsDoubleSpring_Fu, bim(eqs_Fu,otherNodeDOF_uF))
                push!(mapsDoubleSpring_Fp, bim(eqs_Fp,otherNodeDOF_pM))
            end
        end
        # Set the Jacobian maps for the special node
        setfield!(specialNode, :jacobianMaps, SpecialNodeJacobianIndicesMaps(bim(eqs_Fu,DOF_uF), bim(eqs_Fp,DOF_pM), mapsDoubleSpring_Fu, mapsDoubleSpring_Fp))
    end

    # Loop over hinge axis constraints
    for hac in hingeAxisConstraints
        # Unpack hinge axis constraint assembly data
        @unpack eqs_Fp,eqs_λ,DOF_pMS,DOF_λ = hac.assemblyData
        N = length(DOF_λ)
        # Set the Jacobian maps for the hinge axis constraint
        setfield!(hac, :jacobianMaps, HingeAxisConstraintJacobianIndicesMaps{N,6*N}(bim(eqs_Fp,DOF_pMS), bim(eqs_Fp,DOF_λ), bim(eqs_λ,DOF_pMS)))
    end

    # Reset Jacobian values to zero once the sparsity pattern has been set, if applicable
    if resetJacobian
        zero_arrays!(jacobian)
    end

end


# Sets the inertia indices' maps for each element
function set_inertia_indices_maps!(inertia,elements,resetInertia)

    # Helper function
    bim(rows, cols) = build_index_map(inertia, rows, cols)

    # Loop over elements and set inertia matrix maps
    for element in elements
        # Unpack
        @unpack eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FV,eqs_Fχ,eqs_FΩ,DOF_u,DOF_p,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
        # Get number of aerodynamic states
        N = length(DOF_χ)
        # Set
        setfield!(element, :inertiaMaps, 
            InertiaIndicesMaps{3*N,N*N}(
                bim(eqs_Fu1,DOF_p), bim(eqs_Fu1,DOF_V), bim(eqs_Fu1,DOF_Ω),
                bim(eqs_Fp1,DOF_p), bim(eqs_Fp1,DOF_V), bim(eqs_Fp1,DOF_Ω),
                bim(eqs_Fu2,DOF_p), bim(eqs_Fu2,DOF_V), bim(eqs_Fu2,DOF_Ω),
                bim(eqs_Fp2,DOF_p), bim(eqs_Fp2,DOF_V), bim(eqs_Fp2,DOF_Ω),
                bim(eqs_FV,DOF_u),
                bim(eqs_FΩ,DOF_p),
                bim(eqs_Fχ,DOF_V), bim(eqs_Fχ,DOF_Ω), bim(eqs_Fχ,DOF_χ)
            )
        )
    end

    # Reset inertia matrix values to zero once the sparsity pattern has been set, if applicable
    if resetInertia
        zero_arrays!(inertia)
    end

end


# Precomputes the motion of basis A over time
function precompute_motion_basis_A!(problem)

    # Skip if not a dynamic problem
    problem isa DynamicProblem || return

    @unpack timeVector,sizeOfTime = problem.time
    @unpack u_A,v_A,ω_A,vdot_A,ωdot_A,u_AVec,v_AVec,ω_AVec,vdot_AVec,ωdot_AVec,R_AVec = problem.model.kinematicsBasisA
    @unpack R_A = problem.model.kinematicsBasisA.state

    # Save initial rotation tensor
    R_A0 = copy(R_A)

    # Resize rotation tensor history vector
    resize!(R_AVec, sizeOfTime)

    # Loop over time steps -- set rotation tensor array
    for i in 2:sizeOfTime
        # Current time
        t = timeVector[i]
        # Current time step and time factor
        Δt = timeVector[i] - timeVector[i-1]
        α = 2 / Δt
        # Skew-symmetric operator of the angular velocity vector
        Ω = tilde(eval_motion(ω_A, t))
        # Update rotation tensor from basis I to basis A, at the current time
        R_A .= (α*I3 - Ω) \ ((α*I3 + Ω) * R_A)
        # Set in rotation tensor history vector
        R_AVec[i] = SMatrix{3,3}(R_A)
    end

    # Resize motion vectors
    resize!(u_AVec, sizeOfTime)
    resize!(v_AVec, sizeOfTime)
    resize!(ω_AVec, sizeOfTime)
    resize!(vdot_AVec, sizeOfTime)
    resize!(ωdot_AVec, sizeOfTime)

    # Set vectors over time
    u_AVec .= eval_motion.(u_A,timeVector)
    v_AVec .= eval_motion.(v_A,timeVector)
    ω_AVec .= eval_motion.(ω_A,timeVector)
    vdot_AVec .= eval_motion.(vdot_A,timeVector)
    ωdot_AVec .= eval_motion.(ωdot_A,timeVector)

    # Reset rotation tensor to initial value
    R_A .= R_A0

end


# Precomputes the history of each prescribed BC over time, including the maximum absolute values of forces and moments
function precompute_BCs_history!(problem)

    # Unpack time variables for dynamic problem or define dummy one otherwise
    timeVector = time_vector(problem)
    Nt = length(timeVector)

    # Loop over BCs
    for bc in problem.model.BCs
        precompute_BC_history!(bc,timeVector,Nt)
    end

end


# Functional unit of precompute_BCs_history!
function precompute_BC_history!(bc,timeVector,Nt)

    # Unpack
    @unpack R0_n,isTrim = bc.deducedData
    @unpack displacementsOnAofTime,deadLoadsOnAofTime,followerLoadsOnAofTime,displacementsOnbofTime,deadLoadsOnbofTime,followerLoadsOnbofTime = bc.loadFuns

    # Initialize history arrays
    dA = [zeros(MVector{6,Float64}) for _ in 1:Nt]
    dl_A = [zeros(MVector{6,Float64}) for _ in 1:Nt]
    fl_A = [zeros(MVector{6,Float64}) for _ in 1:Nt]
    db = [zeros(MVector{6,Float64}) for _ in 1:Nt]
    dl_b = [zeros(MVector{6,Float64}) for _ in 1:Nt]
    fl_b = [zeros(MVector{6,Float64}) for _ in 1:Nt]
    cv = [zeros(MVector{6,Float64}) for _ in 1:Nt]

    # Initialize maxima
    Fmax = 0.0
    Mmax = 0.0

    # Loop over time
    @inbounds for (k,t) in enumerate(timeVector)

        # Evaluate prescribed functions at current time
        @inbounds for i in 1:6
            dA[k][i] = displacementsOnAofTime[i](t)
            dl_A[k][i] = deadLoadsOnAofTime[i](t)
            fl_A[k][i] = followerLoadsOnAofTime[i](t)
            db[k][i] = displacementsOnbofTime[i](t)
            dl_b[k][i] = deadLoadsOnbofTime[i](t)
            fl_b[k][i] = followerLoadsOnbofTime[i](t)
        end

        # Compute current aggregate values
        valuesOnA = dA[k] + dl_A[k] + fl_A[k]
        valuesOnb = db[k] + dl_b[k] + fl_b[k]
        cv[k] .= (valuesOnA +
                    vcat(
                        R0_n*valuesOnb[1:3],
                        R0_n*valuesOnb[4:6]
                        )
                    ) .* .!isTrim

        # Update maxima
        Fmax = max(Fmax, 
                    maximum(abs, @view dl_A[k][1:3]), 
                    maximum(abs, @view dl_b[k][1:3]), 
                    maximum(abs, @view fl_A[k][1:3]), 
                    maximum(abs, @view fl_b[k][1:3])
                    )

        Mmax = max(Mmax, 
                    maximum(abs, @view dl_A[k][4:6]), 
                    maximum(abs, @view dl_b[k][4:6]), 
                    maximum(abs, @view fl_A[k][4:6]), 
                    maximum(abs, @view fl_b[k][4:6])
                    )
    end

    # Set history field
    setfield!(bc, :history, BCHistory(dA,dl_A,fl_A,db,dl_b,fl_b,cv,Fmax,Mmax))

end


# Pre-computes the distributed loads over all elements' nodes, over every time step
function precompute_distributed_loads!(problem)

    # Unpack time variables for dynamic problem or define dummy one otherwise
    timeVector = time_vector(problem)
    Nt = length(timeVector)

    # Loop over elements and precompute distributed loads over time
    for beam in problem.model.beams, element in beam.elements
        precompute_element_distributed_loads!(element,timeVector,Nt)
    end

end


# Functional unit of precompute_distributed_loads!
function precompute_element_distributed_loads!(element,timeVector,Nt)

    # Unpack element data
    @unpack Δℓ = element.geometry
    @unpack f_A_of_ζt,m_A_of_ζt,f_b_of_ζt,m_b_of_ζt,ff_A_of_ζt,mf_A_of_ζt,ff_b_of_ζt,mf_b_of_ζt,has_f_A,has_m_A,has_f_b,has_m_b,has_ff_A,has_mf_A,has_ff_b,has_mf_b = element.distLoads

    # Initialize arrays of distributed loads on element's nodes over time, loop over time and integrate the distributed loads at each time step, then save the precomputed distributed loads in the element's distLoads field
    if has_f_A
        f_A = zeros(2,3,Nt)
        @inbounds for (i,t) in enumerate(timeVector)
            gauss_legendre3_shaped_fun!(@view(f_A[:,:,i]), ζ -> Δℓ*f_A_of_ζt(ζ, t))
        end
        setfield!(element.distLoads, :f_A, f_A)
    end
    if has_m_A
        m_A = zeros(2,3,Nt)
        @inbounds for (i,t) in enumerate(timeVector)
            gauss_legendre3_shaped_fun!(@view(m_A[:,:,i]), ζ -> Δℓ*m_A_of_ζt(ζ, t))
        end
        setfield!(element.distLoads, :m_A, m_A)
    end
    if has_f_b
        f_b = zeros(2,3,Nt)
        @inbounds for (i,t) in enumerate(timeVector)
            gauss_legendre3_shaped_fun!(@view(f_b[:,:,i]), ζ -> Δℓ*f_b_of_ζt(ζ, t))
        end
        setfield!(element.distLoads, :f_b, f_b)
    end
    if has_m_b
        m_b = zeros(2,3,Nt)
        @inbounds for (i,t) in enumerate(timeVector)
            gauss_legendre3_shaped_fun!(@view(m_b[:,:,i]), ζ -> Δℓ*m_b_of_ζt(ζ, t))
        end
        setfield!(element.distLoads, :m_b, m_b)
    end
    if has_ff_A
        ff_A = zeros(2,3,Nt)
        @inbounds for (i,t) in enumerate(timeVector)
            gauss_legendre3_shaped_fun!(@view(ff_A[:,:,i]), ζ -> Δℓ*ff_A_of_ζt(ζ, t))
        end
        setfield!(element.distLoads, :ff_A, ff_A)
    end
    if has_mf_A
        mf_A = zeros(2,3,Nt)
        @inbounds for (i,t) in enumerate(timeVector)
            gauss_legendre3_shaped_fun!(@view(mf_A[:,:,i]), ζ -> Δℓ*mf_A_of_ζt(ζ, t))
        end
        setfield!(element.distLoads, :mf_A, mf_A)
    end
    if has_ff_b
        ff_b = zeros(2,3,Nt)
        @inbounds for (i,t) in enumerate(timeVector)
            gauss_legendre3_shaped_fun!(@view(ff_b[:,:,i]), ζ -> Δℓ*ff_b_of_ζt(ζ, t))
        end
        setfield!(element.distLoads, :ff_b, ff_b)
    end
    if has_mf_b
        mf_b = zeros(2,3,Nt)
        @inbounds for (i,t) in enumerate(timeVector)
            gauss_legendre3_shaped_fun!(@view(mf_b[:,:,i]), ζ -> Δℓ*mf_b_of_ζt(ζ, t))
        end
        setfield!(element.distLoads, :mf_b, mf_b)
    end

end


# Precomputes the flap deflection and its rates over time
function precompute_flap_deflection!(problem)

    # Skip if not a dynamic problem
    problem isa DynamicProblem || return

    # Time vector
    timeVector = time_vector(problem)

    # Set surfaces array
    surfaces = AeroSurface[beam.attachments.surface for beam in problem.model.beams if !isnothing(beam.attachments.surface)]

    # Loop over surfaces
    for surface in surfaces
        # Unpack functions of time
        @unpack δFun,δdotFun,δddotFun,δMultiplier = surface.assemblyData
        # Compute set flap deflection vector and its rates over time and update surface assembly data
        δVec = map(δFun, timeVector) * δMultiplier
        δdotVec = map(δdotFun, timeVector) * δMultiplier
        δddotVec = map(δddotFun, timeVector) * δMultiplier
        setfield!(surface, :assemblyData, AeroSurfaceAssemblyData(surface, δVec, δdotVec, δddotVec))
        # Set current (initial) flap deflection and its rates, and update surface state
        δ = δVec[1]
        δdot = δdotVec[1]
        δddot = δddotVec[1]
        setfield!(surface, :state, AeroSurfaceState(δ,δdot,δddot))
    end

end


# Solves a steady problem (includes trim and the steady part of eigen problems) 
function solve_steady!(problem)

    # Skip the steady solve and go straight to Jacobian and inertia matrices assembly in problems singular by construction (no BCs)
    if problem isa EigenProblem && isempty(problem.model.BCs)
        j!(problem)
        i!(problem)
        problem.systemSolver.state.converged = false
        @info("Jacobian is singular by construction (no BCs): det(J)=$(det(problem.systemArrays.jacobian)) - skipping steady solve")
        return
    end

    # Check if linear or nonlinear solution is required
    if problem.systemSolver isa LinearSolver
        solve_linear!(problem)
    else #if problem.systemSolver isa NonlinearSolver
        solve_nonlinear_system!(problem)
    end

    # Set trim solution for trim problems
    if problem isa TrimProblem
        set_trim_solution!(problem)
    end

    # Assemble inertia matrix in eigen and trim problems, if applicable
    if problem isa EigenProblem || (problem isa TrimProblem && problem.getInertiaMatrix)
        i!(problem)
    end

end


# Sets the trim solution for a trim problem based on the solution array
function set_trim_solution!(problem)
    
    @unpack x = problem.systemArrays
    @unpack rangeTrimLoads,rangeTrimδ = problem.model.assemblyData
    @unpack forceScaling = problem.model.assemblyData

    setfield!(problem, :trimSolution, TrimSolution(x[rangeTrimLoads]*forceScaling, x[rangeTrimδ]))
    
end


# Eigenvalue solver using Arnoldi method
function eigensolve_arnoldi(map; nev, which, tol, restarts, kwargs...)
    decomp, history = ArnoldiMethod.partialschur(map; nev, which, tol, restarts, kwargs...)
    # history.converged || @info("ArnoldiMethod converged only $(history.nconverged) of $nev eigenvalues")
    return ArnoldiMethod.partialeigen(decomp)
end


# Eigenvalue solver using Krylov subspace method
function eigensolve_krylovkit(map; nev, which, tol, restarts, krylovdim=max(2nev, 30), kwargs...)
    n = size(map, 1)
    krylovdim = min(krylovdim, n)
    vals, vecs, info = KrylovKit.eigsolve(x -> map * x, rand(n), nev, which;
        tol=tol, krylovdim=krylovdim, maxiter=restarts, issymmetric=false, verbosity=KrylovKit.SILENT_LEVEL, kwargs...)
    nconv = info.converged
    nconv == 0 && error("KrylovKit did not converge any eigenvalues")
    # nconv < nev && @info("KrylovKit converged only $nconv of requested $nev eigenvalues")
    return vals[1:nconv], reduce(hcat, vecs[1:nconv])
end


# Eigenvalue solver using ARPACK
function eigensolve_arpack(map; nev, which, tol, restarts, kwargs...)
    n = size(map, 1)
    nevA = min(nev, n - 2)
    vals, vecs, nconv = Arpack.eigs(map; nev=nevA, which=which, tol=tol, maxiter=restarts, kwargs...)
    nconv == 0 && error("ARPACK did not converge any eigenvalues")
    # nconv < nev && @info("ARPACK converged only $nconv of requested $nev eigenvalues")
    return vals, vecs
end


# Full eigenvalue decomposition solver for dense matrices
function eigensolve_fulldecomp(map; kwargs...)
    A = Matrix(map)
    all(isfinite, A) || error("Eigenmatrix contains non-finite values")
    F = eigen(A)
    idx = sortperm(abs.(F.values); rev=true)
    return F.values[idx], F.vectors[:, idx]
end


# Dictionary of available eigensolvers for the generalized eigenvalue problem
const EIGS_SOLVERS = (
    arnoldi    = eigensolve_arnoldi,
    arpack     = eigensolve_arpack,
    fulldecomp = eigensolve_fulldecomp,
    krylovkit  = eigensolve_krylovkit,
)

# Callable instance of the ShiftAndInvert operator for the generalized eigenvalue problem
function (A::ShiftAndInvert)(y, x)
    mul!(A.temp, A.I, x)    # temp = I*x
    ldiv!(y, A.J, A.temp)   # y = J\(I*x)
    rmul!(y, -1)            # y = -J\(I*x)
    return y
end

# Solves the generalized eigenvalue problem for the Jacobian and inertia matrices
function eigs_generalized(luCache, inertia;
    nev=10, which=:LM, σ=0.0, tol=sqrt(eps()), restarts=100,
    method::Symbol=:krylovkit, fallback::Union{Nothing,Symbol}=:fulldecomp,
    solver_kwargs=(;), sortby=abs, rev=false)

    # Ensure the specified method is available in the eigensolvers dictionary
    haskey(EIGS_SOLVERS, method) || error("method must be one of $(keys(EIGS_SOLVERS))")

    # Size of the problem
    n = size(inertia, 1)

    # Shift-and-invert operator for the eigensolver
    op = ShiftAndInvert(luCache, inertia, Vector{Float64}(undef, n))

    # Linear map for the eigensolver
    map = LinearMaps.LinearMap{Float64}(op, n; ismutating=true)

    # Arguments for the eigensolver
    args = (; nev, which, tol, restarts, solver_kwargs...)

    # Solve eigenvalue problem using the specified method, or fall back to the alternative if it fails
    local θ, V
    try
        θ, V = EIGS_SOLVERS[method](map; args...)
    catch err
        (fallback === nothing || fallback === method) && rethrow()
        @warn "Method :$method failed ($(sprint(showerror, err)))"
        @info "Using :$fallback"
        θ, V = EIGS_SOLVERS[fallback](map; args...)
    end

    # Eigenvalues of original problem are the inverses of the eigenvalues of the shifted problem
    λ = σ .+ 1 ./ θ

    # Remove non-finite eigenvalues (and corresponding eigenvectors)
    keep = findall(x -> isfinite(real(x)) && isfinite(imag(x)), λ)
    λ = λ[keep]
    V = V[:, keep]

    # Sort eigenvalues and eigenvectors by the specified criterion
    p = sortperm(sortby.(λ); rev=rev)
    λ = λ[p]
    V = V[:, p]

    return λ, V

end

# Returns the shift to be used to re-solve an eigenproblem solved without shift (σ = 0), or zero if no shift is needed. With σ = 0 the shift-and-invert operator has eigenvalues θ = 1/λ, so eigenvalues very close to zero (e.g. the aerodynamic states in problems without airflow, or rigid body modes) have huge |θ| and set the round-off level of the operator to ~eps*|θ|max. The error in an eigenvalue of frequency ω is then ~eps*|θ|max*ω², which shows up as spurious damping in the oscillatory modes. If the estimated error is not negligible, a small shift (relative to the largest frequency found), which moves those eigenvalues away from the singularity while hardly changing which modes are computed, is returned
function needs_shift(λ,σ,rtol=1e-8)

    ω = maximum(abs ∘ imag, λ; init=0.0)
    θ = maximum(x -> abs(1/(x-σ)), λ; init=0.0)

    # Estimated absolute error of the eigenvalues of largest frequency, compared with a relative accuracy of rtol
    errorEstimate = eps() * θ * ω^2

    return errorEstimate > rtol*ω ? rtol*ω : 0.0

end

# Real snapshot of a complex eigenvector, rotated in the complex plane by a phase common to all entries (so that the linear relations among the states are preserved, which is not the case if each entry is treated independently). The phase maximizes the norm of the real part of the entries in `idx`.
function real_eigenvector(v, idx)

    # Eigenvectors of non-oscillatory modes are already real
    if maximum(abs ∘ imag, v; init=0.0) <= 1e-12 * maximum(abs, v; init=0.0)
        return real.(v)
    end

    # Phase that maximizes the norm of the real part of the displacement-like states
    ϕ = -angle(sum(x -> x^2, view(v, idx))) / 2

    vr = real.(v .* cis(ϕ))

    # Remove the remaining sign ambiguity by making the displacement-like entry of largest magnitude positive
    return vr[idx[argmax(abs.(view(vr, idx)))]] < 0 ? -vr : vr

end


# Computes mode shapes from the given eigenvectors
function mode_shapes(problem,λ,V)

    @unpack forceScaling,elements,nElementsTotal,elementNodes = problem.model.assemblyData

    # Indices of the displacement-like states (displacements, rotation parameters, forces and moments) in the eigenvectors
    idxPos = reduce(vcat, vcat(e.assemblyIndices.DOF_u, e.assemblyIndices.DOF_p, e.assemblyIndices.DOF_F, e.assemblyIndices.DOF_M) for e in elements)

    # Real snapshots of the (complex) eigenvectors, preserving the relative phases among their entries
    Vreal = reduce(hcat, real_eigenvector(@view(V[:,mode]), idxPos) for mode in axes(V,2))

    # Number of modes
    nModes = size(Vreal,2)

    # Frequencies and dampings
    frequencies = imag.(λ)
    dampings = real.(λ)

    # Initialize array of mode shapes
    modeShapes = Vector{ModeShape}(undef, nModes)

    # Loop modes
    for (mode,freq,damp) in zip(1:nModes,frequencies,dampings)
        # Initialize current mode shape
        modeShapes[mode] = ModeShape(freq,damp,[ModalState() for _ in 1:nElementsTotal])
        # Current mode's eigenvector
        @views eigenvector = Vreal[:,mode]
        # Initialize maxima
        upNorm = 0.0
        FNorm = 0.0
        MNorm = 0.0
        VNorm = 0.0
        ΩNorm = 0.0
        γNorm = 0.0
        κNorm = 0.0
        PNorm = 0.0
        HNorm = 0.0
        θNorm = 0.0
        # Loop elements
        for (e,element) in enumerate(elements)
            # Get modal states for current element
            element_modal_states!(element,modeShapes[mode].state[e],eigenvector,forceScaling)
            # Update maxima of absolute values of each state
            @unpack u,p,F,M,V,Ω,γ,κ,P,H = modeShapes[mode].state[e].midpoint
            @unpack u1,u2,p1,p2,F1,F2,M1,M2,θ = modeShapes[mode].state[e].nodal
            upNorm = max(upNorm, maximum(abs, vcat(u,u1,u2,p,p1,p2)))
            FNorm = max(FNorm, maximum(abs, vcat(F,F1,F2)))
            MNorm = max(MNorm, maximum(abs, vcat(M,M1,M2)))
            VNorm = max(VNorm, maximum(abs, V))
            ΩNorm = max(ΩNorm, maximum(abs, Ω))
            γNorm = max(γNorm, maximum(abs, γ))
            κNorm = max(κNorm, maximum(abs, κ))
            PNorm = max(PNorm, maximum(abs, P))
            HNorm = max(HNorm, maximum(abs, H))
            θNorm = max(θNorm, maximum(abs, θ))
        end
        # Normalize mode shapes by maxima
        for (e,element) in enumerate(elements)
            # Unpack
            @unpack u,p,F,M,V,Ω,γ,κ,P,H = modeShapes[mode].state[e].midpoint
            @unpack u1,u2,p1,p2,u1_b,u2_b,p1_b,p2_b,F1,F2,M1,M2,θ = modeShapes[mode].state[e].nodal
            # Normalize states: displacements and rotation parameters share the same normalization factor, so that their relative magnitudes are preserved
            divide_inplace!(upNorm, u,u1,u2,u1_b,u2_b,p,p1,p2,p1_b,p2_b) 
            divide_inplace!(FNorm, F,F1,F2)
            divide_inplace!(MNorm, M,M1,M2)
            V ./= VNorm
            Ω ./= ΩNorm
            γ ./= γNorm
            κ ./= κNorm
            P ./= PNorm
            H ./= HNorm
            θ ./= θNorm
        end
    end
    
    return modeShapes

end

# Solves a dynamic problem 
function solve_dynamic!(problem)

    # Reset tracked history before running, if applicable
    if problem.solutionHistory.index > 1
        setfield!(problem, :solutionHistory, initialize_solution_history(problem.time,problem.model.assemblyData.elements,problem.model.hingeAxisConstraints,problem.trackingTimeSteps,problem.trackingFrequency))
    end

    # Solve the initial time step to get consistent initial states, if applicable
    if !problem.skipInitialStatesUpdate
        solve_initial_dynamic!(problem)
    end

    # Save initial time step data, if applicable
    if problem.trackingTimeSteps
        save_time_step_data!(problem)
    end

    # Time march
    time_march!(problem)

end


# Gets consistent initial states and rates
function solve_initial_dynamic!(problem)

    # Problems with initial conditions: get the initial states from the specified initial generalized displacements and velocities, and then the rates
    if has_initial_conditions(problem.model)
        solve_initial_conditions!(problem)
    # Otherwise, the initial displacements and rotations are those of the model: get the velocities and rates
    else
        initial_rates_startup!(problem)
    end

end


# Checks if any of the model's beams has initial conditions
has_initial_conditions(model) = any(beam.initialConditions.hasInitialConditions for beam in model.beams)


# Gets consistent initial states and rates in problems with initial conditions. The initial generalized displacements and velocities are specified at the beams' nodes in any combination, and the components that are not specified (null) follow from the specified ones through the governing equations. The procedure has three stages:
# 1. Initial displacements, rotations and internal loads: solution of the static problem in which the specified nodal displacements and rotations are prescribed
# 2. Initial velocities: rates of the above solution when the prescribed nodal displacements and rotations move at their specified rates (quasi-static response of the unspecified ones)
# 3. Initial accelerations (and consistent states and rates): the solution of a few steps from the above states and rates, extrapolated back to the initial time (see initial_rates_startup!)
function solve_initial_conditions!(problem)

    # Initial time and first time step
    t0 = problem.time.initialTime
    Δt = problem.time.timeVector[2] - problem.time.timeVector[1]

    # Nodal initial generalized displacements and velocities (resolved in basis b), and flags for them being specified
    positions = initial_nodal_data(problem.model,:position)
    velocities = initial_nodal_data(problem.model,:velocity)
    noRates = fill(zeros(SVector{6,Float64}),size(positions))
    isSpecified(data) = map(d -> d .!= 0,data)

    # Stage 1: solve the static problem with the specified displacements and rotations prescribed
    staticProblem = initial_static_problem(problem)
    solve_initial_static!(staticProblem,get_BCs_for_initial_conditions(staticProblem.model,positions,noRates,isSpecified(positions),t0))

    # Stage 2: get the velocities when the nodal displacements and rotations move at the specified rates
    if any(any(!iszero,velocity) for velocity in velocities)
        # Static problem in which the solved nodal displacements and rotations are prescribed where the velocities are specified, moving at those rates
        velocityProblem = initial_static_problem(problem)
        solve_initial_static!(velocityProblem,get_BCs_for_initial_conditions(velocityProblem.model,solved_nodal_data(staticProblem),velocities,isSpecified(velocities),t0))
        # Quasi-static rates of the elements' displacements and rotation parameters
        rates = quasi_static_rates(velocityProblem,t0,Δt)
        # Set the rates and update the velocities states of the static problem accordingly
        if isnothing(rates)
            @warn "The initial velocities could not be obtained from the specified rates (singular system) - the specified rates are used instead"
        else
            for (element,udot,pdot) in zip(staticProblem.model.assemblyData.elements,rates...)
                element.state.rates.udot .= udot
                element.state.rates.pdot .= pdot
            end
            update_equivalent_states_rates!(staticProblem)
            solve_initial_static_system!(staticProblem)
        end
    end

    # Copy the states to the original problem
    copy_initial_states!(problem,staticProblem)

    # Stage 3: get the accelerations
    initial_rates_startup!(problem)

end


# Gets the nodal generalized displacements (type = :position) or generalized velocities (type = :velocity) specified as initial conditions, resolved in basis b, at both nodes of each element
function initial_nodal_data(model,type)

    @unpack elements = model.assemblyData

    # Initialize nodal data
    data = Matrix{SVector{6,Float64}}(undef,length(elements),2)

    # Loop elements
    for (e,element) in enumerate(elements)
        # Unpack
        @unpack x1_n1,x1_n2 = element.geometry
        @unpack u0_of_x1,p0_of_x1,udot0_of_x1,pdot0_of_x1 = model.beams[element.assemblyData.parentBeamID].initialConditions
        # Get the functions of the translational and rotational components
        uFun,pFun = type == :position ?
                    (u0_of_x1,p0_of_x1) :
                    (udot0_of_x1,pdot0_of_x1)
        # Evaluate at the element's nodes
        data[e,1] = SVector{6,Float64}(vcat(uFun(x1_n1),pFun(x1_n1)))
        data[e,2] = SVector{6,Float64}(vcat(uFun(x1_n2),pFun(x1_n2)))
    end

    return data

end


# Gets the nodal generalized displacements of a solved problem, resolved in basis b, at both nodes of each element
function solved_nodal_data(problem)

    @unpack elements = problem.model.assemblyData

    # Initialize nodal data
    data = Matrix{SVector{6,Float64}}(undef,length(elements),2)

    # Loop elements
    for (e,element) in enumerate(elements)
        @unpack u1_b,u2_b,p1_b,p2_b = element.state.nodal
        data[e,1] = SVector{6,Float64}(vcat(u1_b,p1_b))
        data[e,2] = SVector{6,Float64}(vcat(u2_b,p2_b))
    end

    return data

end


# Builds a lightweight temporary static problem (infinite time step) to solve for the initial states, without affecting the original problem
function initial_static_problem(problem)

    # Copy the model and the system solver
    modelCopy = deepcopy(problem.model)
    systemSolverCopy = deepcopy(problem.systemSolver)
    reset_solver_state!(systemSolverCopy)

    # Time variables with infinite time step
    timeCopy = Time(problem.time.timeVector)
    timeCopy.state.Δt = Inf

    return DynamicProblem(
        model=modelCopy,
        systemSolver=systemSolverCopy,
        skipInitialStatesUpdate=true,
        trackingTimeSteps=false,
        trackingFrequency=1,
        displayProgress=false,
        displayFrequency=1,
        time=timeCopy,
        systemArrays=initialize_system_arrays!(modelCopy,systemSolverCopy,zeros(0)),
        solutionHistory=SolutionHistory()
    )

end


# Gets the BCs of a static problem for the initial states: the model's BCs and, at the nodes, BCs that prescribe the components of the nodal data (resolved in basis b) flagged as `hard`. The prescribed values are `values` at the initial time `t0`, and evolve linearly in time at the given `rates`. At nodes that already have BCs, only the components that are not touched by those BCs are prescribed
function get_BCs_for_initial_conditions(model,values,rates,hard,t0)

    @unpack elements,BCedNodes,nNodesTotal = model.assemblyData

    # Initialize BCs with those of the model
    initialBCs = copy(model.BCs)

    # Types of the BCs for each component of the nodal data (defined in basis b)
    types = (:u1b,:u2b,:u3b,:p1b,:p2b,:p3b)

    # Initialize flag for node having being assigned BCs
    assigned = falses(nNodesTotal)

    # Loop elements and their nodes
    for (e,element) in enumerate(elements)
        @unpack nodesLocalIDs,nodesGlobalIDs,parentBeamID = element.assemblyData
        for n in 1:2
            # Get global node ID and skip if the node has already been assigned
            globalNode = nodesGlobalIDs[n]
            assigned[globalNode] && continue
            assigned[globalNode] = true
            # Components to be prescribed
            prescribe = hard[e,n]
            # In BC'ed nodes, skip the components that are touched by the BCs of the model
            if globalNode in BCedNodes
                R0_n = n == 1 ? element.geometry.R0_n1 : element.geometry.R0_n2
                nodeBCs = filter(BC -> BC.definition.beam.assemblyData.nodesGlobalIDs[BC.definition.node] == globalNode, model.BCs)
                prescribe = prescribe .& .!components_blocked_by_BCs(nodeBCs,R0_n)
            end
            # Skip nodes without components to be prescribed
            any(prescribe) || continue
            # Set the conditions as functions of time that start at the specified values
            conditions = Dict{Symbol,Any}()
            for j in findall(prescribe)
                conditions[types[j]] = linear_function_of_time(values[e,n][j],rates[e,n][j],t0)
            end
            # Create BC with the state at the initial time and add to list
            BC = create_BC(beam=model.beams[parentBeamID],node=nodesLocalIDs[n],conditions=conditions)
            update_BC_state!(BC,t0)
            push!(initialBCs,BC)
        end
    end

    return initialBCs

end


# Linear function of time with the given value at time t0 and rate
linear_function_of_time(value,rate,t0) = t -> value + rate * (t - t0)


# Gets the flags for the components of the nodal generalized displacements (resolved in basis b) that cannot be prescribed because of BCs on the node, which touch components of the generalized displacements and loads (resolved in basis A, in the order u1,u2,u3,p1,p2,p3)
function components_blocked_by_BCs(BCs,R0_n)

    # Get the components touched by the BCs (in basis A)
    touched = falses(6)
    for BC in BCs
        for type in BC.definition.types
            name = String(type)
            # Direction of the type and offset for the rotational ones (rotation parameters and moments)
            direction = parse(Int,name[end-1:end-1])
            offset = (startswith(name,"p") || startswith(name,"M")) ? 3 : 0
            # Types defined on basis b touch all components of basis A that are connected to the direction by the nodal rotation tensor
            if name[end] == 'A'
                touched[offset+direction] = true
            else
                touched[offset .+ findall(!iszero,BC.deducedData.R0_n[:,direction])] .= true
            end
        end
    end

    # Components of basis b that are connected to the touched components of basis A
    blocked = falses(6)
    for j in 1:6
        offset = j <= 3 ? 0 : 3
        blocked[j] = any(touched[offset .+ findall(!iszero,R0_n[:,j-offset])])
    end

    return blocked

end


# Solves the static problem of the initial states with the given BCs
function solve_initial_static!(problem,BCs)

    # Update model with the BCs
    update_model_BCs!(problem.model; BCs=BCs)

    # Null accelerations, as the static problem has no inertial loads (they could be left from a previous solution of the problem)
    for element in problem.model.assemblyData.elements
        zero_arrays!(element.state.rates.Vdot,element.state.rates.Ωdot,element.state.rates.Pdot,element.state.rates.Hdot,element.state.rates.χdot)
    end

    # Update system arrays with the new BCs
    setfield!(problem, :systemArrays, initialize_system_arrays!(problem.model,problem.systemSolver,zeros(0)))

    # Update equivalent states' rates with infinite Δt
    update_equivalent_states_rates!(problem)

    # Solve to find the unspecified displacements, rotations and internal loads
    solve_initial_static_system!(problem)

end


# Solves the system of equations of the static problem of the initial states, and warns if it is unconverged
function solve_initial_static_system!(problem)

    solve_steady!(problem)

    if !problem.systemSolver.loadStepState.convergedFinalSolution
        @warn "The static solution used to get the initial states did not converge - the initial states may be inconsistent"
    end

end


# Gets the rates of the midpoint displacements and rotation parameters of the elements of a solved static problem when the displacements prescribed by the BCs move in time at their rates (quasi-static response). The static solution is linear in the BCs' values, so its rates (time derivatives) are exactly obtained by a Newton step with the converged Jacobian and the BCs' values shifted by their rates. Returns nothing if the Jacobian is singular
function quasi_static_rates(problem,t0,Δt)

    @unpack BCs = problem.model
    @unpack elements = problem.model.assemblyData
    @unpack systemSolver,systemArrays = problem

    # Make sure that the Jacobian is fully updated
    if systemSolver isa NewtonRaphson
        systemSolver.state.skipJacobianUpdate = false
        systemSolver.state.skipAeroJacobianUpdate = false
    end

    # Residual and Jacobian at the converged solution
    f!(problem)
    residual0 = copy(systemArrays.residual)
    j!(problem)
    jacobian = copy(systemArrays.jacobian)

    # Shift the displacements prescribed by the BCs by their rates (obtained by a forward finite difference of the functions of time)
    δt = 1e-4 * Δt
    for BC in BCs
        @unpack currentValue = BC.state
        value = copy(currentValue)
        update_BC_state!(BC,t0+δt)
        rate = (currentValue - value) / δt
        update_BC_state!(BC,t0)
        currentValue .+= rate .* .!BC.deducedData.isLoad
    end

    # Residual with the shifted BCs
    f!(problem)
    residual1 = copy(systemArrays.residual)

    # Restore the BCs
    for BC in BCs
        update_BC_state!(BC,t0)
    end

    # Get the rates by the Newton step
    luCache = lu(jacobian; check=false)
    issuccess(luCache) || return nothing
    Δx = luCache \ (residual0 - residual1)

    # Rates of the elements' midpoint displacements and rotation parameters
    udot = [SVector{3,Float64}(Δx[element.assemblyIndices.DOF_u]) for element in elements]
    pdot = [SVector{3,Float64}(Δx[element.assemblyIndices.DOF_p]) for element in elements]

    return udot,pdot

end


# Gets the initial states and rates, which are consistent with the equations of motion and with the motion of basis A and the displacements and rotations prescribed by the BCs. A few steps are solved from the current states: backward-Euler ones followed by BDF2 ones, both of which damp the high-frequency content that the time step cannot resolve, with the motion of basis A and the prescribed displacements and rotations evaluated at the end of each step (all loads are kept at their values at the initial time). The solution of the last steps lies on the slow manifold of the system (algebraic variables and velocities slaved to others are consistent), and is extrapolated back to the initial time by a cubic polynomial, whose derivative gives the rates. Taking the states and their rates from the same polynomial makes them consistent with each other, so that the trapezoidal rule that follows does not carry the memory of an inconsistent initial rate (which would otherwise persist as an oscillation of the rates at every time step)
function initial_rates_startup!(problem)

    @unpack model,systemArrays = problem
    @unpack elements = model.assemblyData
    @unpack x = systemArrays

    # Initial time and size of the steps (half of the first time step)
    t0 = problem.time.initialTime
    h = (problem.time.timeVector[2] - problem.time.timeVector[1]) / 2

    # Number of steps, of which the first ones are backward-Euler (the others are BDF2)
    nSteps = 6
    nBackwardEuler = 2

    # Save initial states and rates
    xStart = copy(x)
    ratesStart = elements_rates_vector(elements)

    # Initialize history of the states at the begin (index 1) and at the end of each step
    xHistory = Vector{Vector{Float64}}(undef,nSteps+1)
    xHistory[1] = copy(x)

    # Solve the steps
    converged = true
    timeOriginal = problem.time
    try
        for step in 1:nSteps
            # Time variables of the step: the trapezoidal rule's rates, 2/Δt*states - equivalent rates, are those of backward Euler (Δt = 2h) or BDF2 (Δt = 4h/3) with the appropriate equivalent rates
            isBDF2 = step > nBackwardEuler
            setfield!(problem, :time, Time([t0,t0+(isBDF2 ? 4h/3 : 2h)]))
            # Update the motion of basis A and the displacements and rotations prescribed by the BCs, at the end of the step
            update_motion_basis_A!(problem; timeNow=t0+step*h, Δt=h)
            update_BCs_displacements!(model.BCs,t0,t0+step*h)
            # Update equivalent states' rates at the begin of the step
            set_startup_rates!(elements,xHistory,step,isBDF2,h)
            update_equivalent_states_rates!(problem)
            # Solve the step
            solve_time_step!(problem)
            # Stop if unconverged
            if !problem.systemSolver.loadStepState.convergedFinalSolution
                converged = false
                break
            end
            # Save states at the end of the step
            xHistory[step+1] = copy(x)
        end
    finally
        # Reset original time variables, motion of basis A and BCs
        setfield!(problem, :time, timeOriginal)
        update_motion_basis_A!(problem,1)
        update_BCs_state!(problem,1)
    end

    # Load factor
    σ = problem.systemSolver.loadStepState.σ

    # Keep the initial states and rates if the steps did not converge (at full load factor, as the failed solution may have left it reduced)
    if !converged
        @warn "The first steps used to get the initial rates did not converge - the initial rates may be inconsistent"
        x .= xStart
        set_elements_rates!(elements,ratesStart)
        update_states!(model,x,1.0,1,t0,Inf,true)
        return
    end

    # Cubic polynomial through the states at the end of steps 3 to 6 (the first steps are discarded because they have the transient): extrapolate to the initial time, and differentiate there, to get the states and their rates
    x .= 20*xHistory[4] - 45*xHistory[5] + 36*xHistory[6] - 10*xHistory[7]
    xdot = (-37/3*xHistory[4] + 63/2*xHistory[5] - 27*xHistory[6] + 47/6*xHistory[7]) / h
    set_elements_rates_from_derivative!(elements,xdot)

    # Update the states of the elements, hinge axis constraints and special nodes, at the initial time
    update_states!(model,x,σ,1,t0,Inf,true)

end


# Updates the state of the BCs at time t, with the loads kept at their values at time t0 (only the prescribed displacements and rotations advance)
function update_BCs_displacements!(BCs,t0,t)

    for BC in BCs
        @unpack currentValue = BC.state
        @unpack isLoad = BC.deducedData
        update_BC_state!(BC,t0)
        loads = copy(currentValue)
        update_BC_state!(BC,t)
        for i in 1:6
            isLoad[i] && (currentValue[i] = loads[i])
        end
    end

end


# Sets the rates at the begin of a startup step, so that the equivalent rates (rates plus 2/Δt times the states) are those of backward Euler (null rates) or BDF2 (rates (x₁-x₂)/(2h), with x₁ and x₂ the states at the previous two steps)
function set_startup_rates!(elements,xHistory,step,isBDF2,h)

    for element in elements
        @unpack udot,pdot,Vdot,Ωdot,χdot = element.state.rates
        if isBDF2
            @unpack DOF_u,DOF_p,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
            x₁ = xHistory[step]
            x₂ = xHistory[step-1]
            udot .= (x₁[DOF_u] - x₂[DOF_u]) / 2h
            pdot .= (x₁[DOF_p] - x₂[DOF_p]) / 2h
            Vdot .= (x₁[DOF_V] - x₂[DOF_V]) / 2h
            Ωdot .= (x₁[DOF_Ω] - x₂[DOF_Ω]) / 2h
            χdot .= (x₁[DOF_χ] - x₂[DOF_χ]) / 2h
        else
            zero_arrays!(udot,pdot,Vdot,Ωdot,χdot)
        end
    end

end


# Gets the rates of the elements' states as a single vector
function elements_rates_vector(elements)

    rates = Float64[]

    for element in elements
        @unpack udot,pdot,Vdot,Ωdot,χdot = element.state.rates
        append!(rates,udot,pdot,Vdot,Ωdot,χdot)
    end

    return rates

end


# Sets the rates of the elements' states from a single vector (in the format of elements_rates_vector)
function set_elements_rates!(elements,rates)

    i = 0

    for element in elements
        @unpack udot,pdot,Vdot,Ωdot,χdot = element.state.rates
        # Set the rates of the states
        for rate in (udot,pdot,Vdot,Ωdot,χdot)
            rate .= rates[i+1:i+length(rate)]
            i += length(rate)
        end
        # Set the rates of the momenta
        set_momenta_rates!(element)
    end

end


# Sets the rates of the elements' states from the vector of the time derivative of the system's states
function set_elements_rates_from_derivative!(elements,xdot)

    for element in elements
        @unpack udot,pdot,Vdot,Ωdot,χdot = element.state.rates
        @unpack DOF_u,DOF_p,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
        # Set the rates of the states
        udot .= xdot[DOF_u]
        pdot .= xdot[DOF_p]
        Vdot .= xdot[DOF_V]
        Ωdot .= xdot[DOF_Ω]
        χdot .= xdot[DOF_χ]
        # Set the rates of the momenta
        set_momenta_rates!(element)
    end

end


# Sets the rates of the momenta of the element from the rates of the velocities
function set_momenta_rates!(element)

    @unpack Vdot,Ωdot,Pdot,Hdot = element.state.rates
    @unpack I11,I12,I21,I22 = element.sectional

    Pdot .= I11 * Vdot + I12 * Ωdot
    Hdot .= I21 * Vdot + I22 * Ωdot

end


# Updates the equivalent states' rates at the begin of the current time step 
function update_equivalent_states_rates!(problem)

    @unpack Δt = problem.time.state
    @unpack elementGroups = problem.model.assemblyData

    # Loop over elements and update equivalent states' rates
    foreach_element(elementGroups) do element
        update_element_equivalent_states_rates!(element,Δt)
    end

end


# Functional unit of update_equivalent_states_rates!
function update_element_equivalent_states_rates!(element,Δt)

    @unpack u,p,V,Ω,χ = element.state.midpoint
    @unpack udot,pdot,Vdot,Ωdot,χdot = element.state.rates
    @unpack udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,χdotEquiv = element.stateCache.equivStatesRates

    # Set equivalent rates
    udotEquiv .= udot + 2/Δt * u
    pdotEquiv .= pdot + 2/Δt * p
    VdotEquiv .= Vdot + 2/Δt * V
    ΩdotEquiv .= Ωdot + 2/Δt * Ω
    χdotEquiv .= χdot + 2/Δt * χ

end


# Copies the initial states and the rates of the displacements and rotation parameters to the original problem
function copy_initial_states!(problem,problemCopy)

    # Unpack problem copy
    xCopy = problemCopy.systemArrays.x
    modelCopy = problemCopy.model

    # Unpack original problem
    @unpack x = problem.systemArrays
    @unpack forceScaling = problem.model.assemblyData

    # Loop over elements
    for (element,elementCopy) in zip(problem.model.assemblyData.elements,modelCopy.assemblyData.elements)
        # Unpack
        @unpack u,p,F,M,V,Ω,χ = element.state.midpoint
        @unpack udot,pdot = element.state.rates
        @unpack DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
        DOF_uCopy = elementCopy.assemblyIndices.DOF_u
        DOF_pCopy = elementCopy.assemblyIndices.DOF_p
        DOF_FCopy = elementCopy.assemblyIndices.DOF_F
        DOF_MCopy = elementCopy.assemblyIndices.DOF_M
        DOF_VCopy = elementCopy.assemblyIndices.DOF_V
        DOF_ΩCopy = elementCopy.assemblyIndices.DOF_Ω
        DOF_χCopy = elementCopy.assemblyIndices.DOF_χ
        # Set elemental states
        u .= x[DOF_u] .= xCopy[DOF_uCopy]
        p .= x[DOF_p] .= xCopy[DOF_pCopy]
        F .= x[DOF_F] .= xCopy[DOF_FCopy]
        M .= x[DOF_M] .= xCopy[DOF_MCopy]
        V .= x[DOF_V] .= xCopy[DOF_VCopy]
        Ω .= x[DOF_Ω] .= xCopy[DOF_ΩCopy]
        χ .= x[DOF_χ] .= xCopy[DOF_χCopy]
        F .*= forceScaling
        M .*= forceScaling
        # Set elemental rates of the displacements and rotation parameters
        udot .= elementCopy.state.rates.udot
        pdot .= elementCopy.state.rates.pdot
    end

    # Loop over special nodes
    for specialNode in problem.model.assemblyData.specialNodes
        # Get equivalent copied special node
        matchingSpecialNode(s) = s.assemblyData.globalID == specialNode.assemblyData.globalID
        specialNodeCopyIndex = findfirst(matchingSpecialNode,modelCopy.assemblyData.specialNodes)
        specialNodeCopy = modelCopy.assemblyData.specialNodes[specialNodeCopyIndex]
        # Unpack
        @unpack uIsPrescribed,pIsPrescribed = specialNode.flags
        @unpack DOF_uF,DOF_pM = specialNode.assemblyData
        @unpack u,p,F,M = specialNodeCopy.state
        # Set nodal states (generalized displacements, or forces where displacements are prescribed)
        for i=1:3
            x[DOF_uF[i]] = uIsPrescribed[i] ?
                           F[i]/forceScaling :
                           u[i]
            x[DOF_pM[i]] = pIsPrescribed[i] ?
                           M[i]/forceScaling :
                           p[i]
        end
    end

end


# # Updates the velocity states by running a very small time step (legacy version)
# function update_initial_velocities!(problem)

#     @unpack initialVelocitiesUpdateOptions,model = problem
#     @unpack Δt,maxIter,relaxFactor,tol,displayProgress = initialVelocitiesUpdateOptions
#     @unpack beams = model
#     @unpack elements = model.assemblyData

#     # Get original time object
#     timeOriginal = problem.time

#     # Set new time object
#     timeVector0 = [timeOriginal.initialTime, timeOriginal.initialTime + Δt]
#     timeNow = timeVector0[2]
#     setfield!(problem, :time, Time(timeVector0))

#     # Update motion of basis A
#     update_motion_basis_A!(problem; timeNow=timeNow,Δt=Δt)

#     # Update boundary conditions
#     update_BCs_state_at_current_time!(problem,timeNow)

#     # Initial sectional velocities and displacement's rates
#     V0 = deepcopy([e.state.midpoint.V for e in elements])
#     Ω0 = deepcopy([e.state.midpoint.Ω for e in elements])
#     udot0 = deepcopy([e.state.rates.udot for e in elements])
#     pdot0 = deepcopy([e.state.rates.pdot for e in elements])

#     # Initialize convergence variables
#     ϵ = 1
#     iter = 0

#     # Loop until converged velocities are found
#     while ϵ > tol
#         # Update counter
#         iter += 1
#         # Loop over elements: set relaxation factor on velocity states in need of update to stabilize
#         if iter > 1
#             for (e,element) in enumerate(elements)
#                 # Unpack
#                 @unpack x = problem.systemArrays
#                 @unpack DOF_V,DOF_Ω = element.assemblyIndices
#                 @unpack V,Ω = element.state.midpoint
#                 @unpack udot,pdot = element.state.rates
#                 # Update velocities with relaxation factor
#                 udot .= udot * relaxFactor + udot0[e] * (1-relaxFactor)
#                 pdot .= pdot * relaxFactor + pdot0[e] * (1-relaxFactor)
#                 x[DOF_V] .= V .= V * relaxFactor + V0[e] * (1-relaxFactor)
#                 x[DOF_Ω] .= Ω .= Ω * relaxFactor + Ω0[e] * (1-relaxFactor)
#             end
#         end
#         # Sectional velocities at begin of time step 
#         velBegin = [[e.state.midpoint.V; e.state.midpoint.Ω] for e in elements]
#         # Reset sectional accelerations to zero 
#         for element in elements
#             zero_arrays!(element.state.rates.Vdot, element.state.rates.Ωdot)
#         end
#         # Update equivalent states' rates at the begin of the time step
#         update_equivalent_states_rates!(problem)
#         # Solve the system at the current time step
#         solve_time_step!(problem)
#         # Sectional velocities at end of time step 
#         velEnd = [[e.state.midpoint.V; e.state.midpoint.Ω] for e in elements]
#         # Update current difference between initial and final velocities 
#         ϵ = maximum(maximum([abs.(velBegin[e].-velEnd[e]) for e in eachindex(elements)]))
#         # Print iteration results
#         displayProgress && println("iter: $iter, ϵ = $ϵ")
#         # Check number of iterations
#         if iter == maxIter
#             println("Unconverged initial velocities - solution of generalized velocities and accelerations may be unstable: iter = $iter, ϵ = $ϵ")
#             break
#         end
#     end

#     # Reset original time variables
#     setfield!(problem, :time, timeOriginal)

# end


# Marches the dynamic problem in time
function time_march!(problem)

    @unpack trackingTimeSteps,trackingFrequency,displayProgress,displayFrequency = problem
    @unpack sizeOfTime = problem.time

    # Advance time
    for timeIndex = 2:sizeOfTime
        # Update time variables 
        update_time_state!(problem,timeIndex)
        # Update motion of basis A
        update_motion_basis_A!(problem,timeIndex)
        # Update boundary conditions
        update_BCs_state!(problem,timeIndex)
        # Update equivalent states' rates at the begin of the time step
        update_equivalent_states_rates!(problem)
        # Update DS models complementary variables of previous time step, if applicable
        update_BL_complementary_variables!(problem)
        # Solve the system at the current time step
        solve_time_step!(problem)
        # Stop if unconverged
        if !problem.systemSolver.loadStepState.convergedFinalSolution
            println("Unconverged solution, stopping...")
            break
        end
        # Save time step data, if applicable    
        if trackingTimeSteps && iszero(rem(timeIndex-1,trackingFrequency))
            save_time_step_data!(problem)           
        end
        # Display progress, if applicable
        if displayProgress && (iszero(rem(timeIndex-1,displayFrequency)) || timeIndex == sizeOfTime)
            progress = round(timeIndex/sizeOfTime*100,digits=1)
            println("Simulation progress: $progress %")
        end
    end

end


# Updates the current time state
function update_time_state!(problem,timeIndex)

    @unpack timeVector = problem.time

    setfield!(problem.time.state, :timeNow, timeVector[timeIndex])
    setfield!(problem.time.state, :timeIndex, timeIndex)
    setfield!(problem.time.state, :Δt, timeVector[timeIndex] - timeVector[timeIndex-1])

end


# Updates the motion of basis A at the specified time instance and time step
function update_motion_basis_A!(problem; timeNow=problem.time.state.timeNow, Δt=problem.time.state.Δt)

    @unpack u_A,v_A,ω_A,vdot_A,ωdot_A = problem.model.kinematicsBasisA
    @unpack R_A,R_AT = problem.model.kinematicsBasisA.state

    # Time step factor
    α = 2 / Δt

    # Skew-symmetric operator of the angular velocity vector
    Ω = tilde(ω_A(timeNow))

    # Rotation tensor from basis I to basis A, at the current time, and its transpose
    R_A .= (α*I3 - Ω) \ ((α*I3 + Ω) * R_A)
    R_AT .= R_A'

    # Update kinematics of basis A
    problem.model.kinematicsBasisA.state.u_A .= MVector{3,Float64}(u_A(timeNow))
    problem.model.kinematicsBasisA.state.v_A .= MVector{3,Float64}(v_A(timeNow))
    problem.model.kinematicsBasisA.state.ω_A .= MVector{3,Float64}(ω_A(timeNow))
    problem.model.kinematicsBasisA.state.vdot_A .= MVector{3,Float64}(vdot_A(timeNow))
    problem.model.kinematicsBasisA.state.ωdot_A .= MVector{3,Float64}(ωdot_A(timeNow))

end


# Updates the motion of basis A at the current time index, using pre-computed arrays
function update_motion_basis_A!(problem,timeIndex)

    @unpack u_AVec,v_AVec,ω_AVec,vdot_AVec,ωdot_AVec,R_AVec = problem.model.kinematicsBasisA
    @unpack u_A,v_A,ω_A,vdot_A,ωdot_A,R_A,R_AT = problem.model.kinematicsBasisA.state

    # Set current values
    u_A .= u_AVec[timeIndex]
    v_A .= v_AVec[timeIndex]
    ω_A .= ω_AVec[timeIndex]
    vdot_A .= vdot_AVec[timeIndex]
    ωdot_A .= ωdot_AVec[timeIndex]
    R_A .= R_AVec[timeIndex]
    R_AT .= R_A'

end


# Updates the boundary conditions at the current time index
function update_BCs_state!(problem,timeIndex)

    for BC in problem.model.BCs
        update_BC_state!(BC,timeIndex)
    end

end


# Updates the complementary variables of previous time step for the dynamic stall models
function update_BL_complementary_variables!(problem)

    # Loop elements
    for element in problem.model.assemblyData.elements
        # Unpack element data
        @unpack surface,stateCache = element.aero
        # Skip elements without aerodynamic surface or whose solver is not a DynamicStallSolver
        (isnothing(surface) || !isa(surface.solvers.ppSolver, DynamicStallSolver)) && continue
        # Update according to solver
        update_element_complementary_variables!(stateCache,surface.solvers.ppSolver)
    end

end

# Functional units of update_BL_complementary_variables! for the BLi and BLo solvers
function update_element_complementary_variables!(s,solver::BLi)
    setfield!(s.bliVars.comp, :stallOnsetRatioPrev, s.bliVars.flow.stallOnsetRatio)
    setfield!(s.bliVars.comp, :αlagPrev, s.bliVars.states.αlag)
    setfield!(s.bliVars.comp, :qRPrev, s.bliVars.flow.qR)
    setfield!(s.bliVars.comp, :PPrev, s.bliVars.flow.P)
    setfield!(s.bliVars.comp, :upstrokePrev, s.bliVars.flow.upstroke)
end
function update_element_complementary_variables!(s,solver::BLo)
    setfield!(s.bloVars.comp, :stallOnsetRatioPrev, s.bloVars.flow.stallOnsetRatio)
end


# Solves the current time step
function solve_time_step!(problem)

    # Check if linear or nonlinear solution is required
    if problem.systemSolver isa LinearSolver
        solve_linear!(problem)
    else #if problem.systemSolver isa NonlinearSolver
        solve_nonlinear_system!(problem)
    end

end


# Saves the solution at the current load step
function save_loadstep_data!(problem)

    @unpack σ,states,hacStates,aeroVariables,bliAeroVariables,bloAeroVariables = problem.loadstepHistory
    @unpack hingeAxisConstraints = problem.model
    @unpack elements,nElementsTotal = problem.model.assemblyData
    @unpack systemSolver = problem

    # Pre-allocate space for arrays
    if length(σ) == 1 && systemSolver isa NonlinearSolver
        nStepsMax = round(Int, 1/systemSolver.loadStepProps.minLoadFac)
        sizehint!(σ, nStepsMax)
        sizehint!(states, nStepsMax)
        sizehint!(hacStates, nStepsMax)
        sizehint!(aeroVariables, nStepsMax)
        sizehint!(bliAeroVariables, nStepsMax)
        sizehint!(bloAeroVariables, nStepsMax)
    end

    # Add to saved load factors array
    push!(σ, systemSolver.loadStepState.σ)

    # Pre-allocate this step's data
    stepStates = Vector{ElementState}(undef, nElementsTotal)
    stepHACStates = Vector{Union{Nothing,HACStateUnion}}(nothing, length(hingeAxisConstraints))
    stepAeroVars = Vector{Union{Nothing,AttachedAeroVariables{Float64}}}(nothing, nElementsTotal)
    stepBLiAeroVars = Vector{Union{Nothing,BLiAeroVariables{Float64}}}(nothing, nElementsTotal)
    stepBLoAeroVars = Vector{Union{Nothing,BLoAeroVariables{Float64}}}(nothing, nElementsTotal)

    # Loop over elements
    for (e, element) in enumerate(elements)
        save_element_loadstep_data!(e,element,stepStates,stepAeroVars,stepBLiAeroVars,stepBLoAeroVars)
    end

    # Loop over hinge axis constraints
    for (h,hac) in enumerate(hingeAxisConstraints)
        save_hac_loadstep_data!(h,hac,stepHACStates)
    end

    # Push this step's vectors into the history
    push!(states, stepStates)
    push!(hacStates, stepHACStates)
    push!(aeroVariables, stepAeroVars)
    push!(bliAeroVariables, stepBLiAeroVars)
    push!(bloAeroVariables, stepBLoAeroVars)
    
end


# Functional unit of save_loadstep_data! for elements
function save_element_loadstep_data!(e,element,states,aeroVariables,bliAeroVariables,bloAeroVariables)

    # Alias
    s = element.aero.surface

    # Set states
    states[e] = copy(element.state)

    # Set attached-flow aerodynamic quantities
    aeroVariables[e] = isnothing(s) ? 
                       nothing : 
                       copy(element.aero.stateCache.attachedVars)

    # Set BLi aerodynamic quantities
    bliAeroVariables[e] = isnothing(s) || !isa(s.solvers.ppSolver, BLi) ? 
                          nothing : 
                          copy(element.aero.stateCache.bliVars)

    # Set BLo aerodynamic quantities
    bloAeroVariables[e] = isnothing(s) || !isa(s.solvers.ppSolver, BLo) ? 
                          nothing : 
                          copy(element.aero.stateCache.bloVars)

end


# Functional unit of save_loadstep_data! for hinge axis constraints
function save_hac_loadstep_data!(h,hac,stepHACStates)
    stepHACStates[h] = copy(hac.state)
end


# Saves the solution at the current time step
function save_time_step_data!(problem)

    @unpack states,hacStates,aeroVariables,bliAeroVariables,bloAeroVariables,index = problem.solutionHistory

    # Loop over elements
    for (e, element) in enumerate(problem.model.assemblyData.elements)
        save_element_time_step_data!(e,element,states,aeroVariables,bliAeroVariables,bloAeroVariables,index)
    end

    # Loop over hinge axis constraints
    for (h,hac) in enumerate(problem.model.hingeAxisConstraints)
        save_hac_time_step_data!(h,hac,hacStates,index)
    end

    # Update solution history index
    setfield!(problem.solutionHistory, :index, index + 1)

end


# Functional unit of save_time_step_data!
function save_element_time_step_data!(e,element,states,aeroVariables,bliAeroVariables,bloAeroVariables,index)

    # Alias
    s = element.aero.surface

    # Set states
    dest = states[index][e]
    copy_fields!(dest.midpoint, element.state.midpoint)
    copy_fields!(dest.nodal,    element.state.nodal)
    copy_fields!(dest.rates,    element.state.rates)

    # Set attached-flow aerodynamic quantities
    aeroVariables[index][e] = isnothing(s) ? 
                              nothing : 
                              copy(element.aero.stateCache.attachedVars)

    # Set BLi aerodynamic quantities
    bliAeroVariables[index][e] = isnothing(s) || !isa(s.solvers.ppSolver, BLi) ? 
                                 nothing : 
                                 copy(element.aero.stateCache.bliVars)

    # Set BLo aerodynamic quantities
    bloAeroVariables[index][e] = isnothing(s) || !isa(s.solvers.ppSolver, BLo) ? 
                                 nothing : 
                                 copy(element.aero.stateCache.bloVars)

end

# Functional unit of save_time_step_data! for hinge axis constraints
function save_hac_time_step_data!(h,hac,hacStates,index)
    hacStates[index][h] = copy(hac.state)
end

# Returns the current time index of the problem
@inline time_index(problem::DynamicProblem) = problem.time.state.timeIndex
@inline time_index(problem) = 1

# Returns the current time of the problem
@inline time_now(problem::DynamicProblem) = problem.time.state.timeNow
@inline time_now(problem) = 0.0

# Returns the time vector of the problem
@inline time_vector(problem::DynamicProblem) = problem.time.timeVector
@inline time_vector(problem) = [0.0]

# Returns the time step of the problem
@inline time_step(problem::DynamicProblem) = problem.time.state.Δt
@inline time_step(problem) = Inf