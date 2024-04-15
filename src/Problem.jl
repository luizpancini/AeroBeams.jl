
"""
abstract type Problem

Defines a general Problem type

"""
abstract type Problem end


"""
abstract type SystemSolver

Defines a general problem for the system of equations

"""
abstract type SystemSolver end


"""
@with_kw mutable struct InitialVelocitiesUpdateOptions

InitialVelocitiesUpdateOptions composite type

Defines variables for the update of the initial velocities states

# Fields
- Δt::Number 
- maxIter::Int64 
- relaxFactor::Float64 
- tol::Float64
"""
@with_kw mutable struct InitialVelocitiesUpdateOptions

    # Fields
    Δt::Number = 0
    maxIter::Int64 = 2
    relaxFactor::Number = 0.5
    tol::Float64 = 1e-6
    displayProgress::Bool = false

    # Constructor
    function InitialVelocitiesUpdateOptions(Δt::Number,maxIter::Int64,relaxFactor::Number,tol::Float64,displayProgress::Bool)

        # Check relaxation factor
        @assert (0 < relaxFactor <= 1)

        return new(Δt,maxIter,relaxFactor,tol,displayProgress)
    end
end
export InitialVelocitiesUpdateOptions


"""
mutable struct ModeShape{T<:Union{Float64,ComplexF64}}

ModeShape composite type

# Fields
- mode::Int64
- frequency::Float64
- damping::Float64
- elementalStates::Vector{ElementalStates{T}}
- complementaryElementalStates::Vector{ComplementaryElementalStates{T}}
- nodalStates::Vector{NodalStates{T}}
"""
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


"""
@with_kw mutable struct SteadyProblem <: Problem

SteadyProblem composite type

Defines the problem of steady type

# Fields
-
"""
@with_kw mutable struct SteadyProblem <: Problem

    # Model
    model::Model = Model()
    # States, residual, Jacobian and inertia arrays
    x::Vector{Float64} = zeros(0)
    Δx::Vector{Float64} = zeros(0)
    residual::Vector{Float64} = zeros(0)
    jacobian::Matrix{Float64} = zeros(0,0)
    inertia::Matrix{Float64} = zeros(0,0)
    jacobianDeterminant::Float64 = 0.0
    # Dummy time
    timeNow::Float64 = 0.0
    # TF to get linear solution
    getLinearSolution::Bool = false
    # System solver
    systemSolver::SystemSolver = NewtonRaphson()
    # Load factor
    σ::Float64 = 1.0
    # TF to compute only the external forces array at the current nonlinear step (used only for arclength system solver)
    getExternalForcesArray::Bool = false
    # Array of partial load steps and respective solutions
    savedσ::Vector{Float64} = Vector{Float64}()
    xOverσ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    elementalStatesOverσ::Vector{Vector{ElementalStates{Float64}}} = Vector{Vector{ElementalStates{Float64}}}()
    nodalStatesOverσ::Vector{Vector{NodalStates{Float64}}} = Vector{Vector{NodalStates{Float64}}}()
    compElementalStatesOverσ::Vector{Vector{ComplementaryElementalStates{Float64}}} = Vector{Vector{ComplementaryElementalStates{Float64}}}()


    # Constructor
    function SteadyProblem(model::Model,x::Vector{Float64},Δx::Vector{Float64},residual::Vector{Float64},jacobian::Matrix{Float64},inertia::Matrix{Float64},jacobianDeterminant::Float64,timeNow::Float64,getLinearSolution::Bool,systemSolver::SystemSolver,σ::Float64,getExternalForcesArray::Bool,savedσ::Vector{Float64},xOverσ::Vector{Vector{Float64}},elementalStatesOverσ::Vector{Vector{ElementalStates{Float64}}},nodalStatesOverσ::Vector{Vector{NodalStates{Float64}}},compElementalStatesOverσ::Vector{Vector{ComplementaryElementalStates{Float64}}})

        # Initialize problem
        problem = new(model,x,Δx,residual,jacobian,inertia,jacobianDeterminant,timeNow,getLinearSolution,systemSolver,σ,getExternalForcesArray,savedσ,xOverσ,elementalStatesOverσ,nodalStatesOverσ,compElementalStatesOverσ)

        # Set initial elemental and nodal states
        set_initial_states!(problem)

        # Initialize system arrays with correct size
        initialize_system_arrays!(problem)

        # Update initial load factor
        problem.σ = systemSolver.initialLoadFactor

        return problem

    end
end
export SteadyProblem


"""
@with_kw mutable struct TrimProblem <: Problem

TrimProblem composite type

Defines the problem of trim type

# Fields
-
"""
@with_kw mutable struct TrimProblem <: Problem

    # Model
    model::Model = Model()
    # States, residual, Jacobian and inertia arrays
    x::Vector{Float64} = zeros(0)
    Δx::Vector{Float64} = zeros(0)
    residual::Vector{Float64} = zeros(0)
    jacobian::Matrix{Float64} = zeros(0,0)
    inertia::Matrix{Float64} = zeros(0,0)
    jacobianDeterminant::Float64 = 0.0
    # Dummy time
    timeNow::Float64 = 0.0
    # TF to get linear solution
    getLinearSolution::Bool = false
    # System solver
    systemSolver::SystemSolver = NewtonRaphson()
    # Load factor
    σ::Float64 = 1.0
    # TF to compute only the external forces array at the current nonlinear step (used only for arclength system solver)
    getExternalForcesArray::Bool = false
    # Array of partial load steps and respective solutions
    savedσ::Vector{Float64} = Vector{Float64}()
    xOverσ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    elementalStatesOverσ::Vector{Vector{ElementalStates{Float64}}} = Vector{Vector{ElementalStates{Float64}}}()
    nodalStatesOverσ::Vector{Vector{NodalStates{Float64}}} = Vector{Vector{NodalStates{Float64}}}()
    compElementalStatesOverσ::Vector{Vector{ComplementaryElementalStates{Float64}}} = Vector{Vector{ComplementaryElementalStates{Float64}}}()

    
    # Constructor
    function TrimProblem(model::Model,x::Vector{Float64},Δx::Vector{Float64},residual::Vector{Float64},jacobian::Matrix{Float64},inertia::Matrix{Float64},jacobianDeterminant::Float64,timeNow::Float64,getLinearSolution::Bool,systemSolver::SystemSolver,σ::Float64,getExternalForcesArray::Bool,savedσ::Vector{Float64},xOverσ::Vector{Vector{Float64}},elementalStatesOverσ::Vector{Vector{ElementalStates{Float64}}},nodalStatesOverσ::Vector{Vector{NodalStates{Float64}}},compElementalStatesOverσ::Vector{Vector{ComplementaryElementalStates{Float64}}})

        # Initialize problem
        problem = new(model,x,Δx,residual,jacobian,inertia,jacobianDeterminant,timeNow,getLinearSolution,systemSolver,σ,getExternalForcesArray,savedσ,xOverσ,elementalStatesOverσ,nodalStatesOverσ,compElementalStatesOverσ)

        # Set initial elemental and nodal states
        set_initial_states!(problem)

        # Initialize system arrays with correct size
        initialize_system_arrays!(problem)

        # Update initial load factor
        problem.σ = systemSolver.initialLoadFactor

        return problem

    end
end
export TrimProblem


"""
@with_kw mutable struct EigenProblem <: Problem

EigenProblem composite type

Defines the problem of eigen type

# Fields
-
"""
@with_kw mutable struct EigenProblem <: Problem

    # Model
    model::Model = Model()
    # States, residual, Jacobian and inertia arrays
    x::Vector{Float64} = zeros(0)
    Δx::Vector{Float64} = zeros(0)
    residual::Vector{Float64} = zeros(0)
    jacobian::Matrix{Float64} = zeros(0,0)
    inertia::Matrix{Float64} = zeros(0,0)
    jacobianDeterminant::Float64 = 0.0
    # Dummy time
    timeNow::Float64 = 0.0
    # TF to get linear solution
    getLinearSolution::Bool = false
    # System solver
    systemSolver::SystemSolver = NewtonRaphson()
    # Load factor
    σ::Float64 = 1.0
    # TF to compute only the external forces array at the current nonlinear step (used only for arclength system solver)
    getExternalForcesArray::Bool = false
    # Array of partial load steps and respective solutions
    savedσ::Vector{Float64} = Vector{Float64}()
    xOverσ::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    elementalStatesOverσ::Vector{Vector{ElementalStates{Float64}}} = Vector{Vector{ElementalStates{Float64}}}()
    nodalStatesOverσ::Vector{Vector{NodalStates{Float64}}} = Vector{Vector{NodalStates{Float64}}}()
    compElementalStatesOverσ::Vector{Vector{ComplementaryElementalStates{Float64}}} = Vector{Vector{ComplementaryElementalStates{Float64}}}()
    # Frequencies, dampings and eigenvectors
    frequencies::Vector{Float64} = Vector{Float64}()
    dampings::Vector{Float64} = Vector{Float64}()
    eigenvectors::Matrix{ComplexF64} = zeros(ComplexF64, 0, 0)
    frequenciesFiltered::Vector{Float64} = Vector{Float64}()
    dampingsFiltered::Vector{Float64} = Vector{Float64}()
    eigenvectorsFiltered::Matrix{ComplexF64} = zeros(ComplexF64, 0, 0)
    frequenciesOscillatory::Vector{Float64} = Vector{Float64}()
    dampingsOscillatory::Vector{Float64} = Vector{Float64}()
    eigenvectorsOscillatoryCplx::Matrix{ComplexF64} = zeros(ComplexF64, 0, 0)
    eigenvectorsOscillatoryAbs::Matrix{Float64} = zeros(0, 0)
    # Number of desired oscillatory modes
    nModes::Int64 = Inf64
    # Frequency filter limits
    frequencyFilterLimits::Vector{Float64} = [0,Inf64]
    # Mode shapes (complex-valued and absolute-valued)
    modeShapesCplx::Vector{ModeShape{ComplexF64}} = Vector{ModeShape{ComplexF64}}()
    modeShapesAbs::Vector{ModeShape{Float64}} = Vector{ModeShape{Float64}}()
    # TF to normalize mode shapes
    normalizeModeShapes::Bool = false

    # Constructor
    function EigenProblem(model::Model,x::Vector{Float64},Δx::Vector{Float64},residual::Vector{Float64},jacobian::Matrix{Float64},inertia::Matrix{Float64},jacobianDeterminant::Float64,timeNow::Float64,getLinearSolution::Bool,systemSolver::SystemSolver,σ::Float64,getExternalForcesArray::Bool,savedσ::Vector{Float64},xOverσ::Vector{Vector{Float64}},elementalStatesOverσ::Vector{Vector{ElementalStates{Float64}}},nodalStatesOverσ::Vector{Vector{NodalStates{Float64}}},compElementalStatesOverσ::Vector{Vector{ComplementaryElementalStates{Float64}}},frequencies::Vector{Float64},dampings::Vector{Float64},eigenvectors::Matrix{ComplexF64},frequenciesFiltered::Vector{Float64},dampingsFiltered::Vector{Float64},eigenvectorsFiltered::Matrix{ComplexF64},frequenciesOscillatory::Vector{Float64},dampingsOscillatory::Vector{Float64},eigenvectorsOscillatoryCplx::Matrix{ComplexF64},eigenvectorsOscillatoryAbs::Matrix{Float64},nModes::Int64,frequencyFilterLimits::Vector{Float64},modeShapesCplx::Vector{ModeShape{ComplexF64}},modeShapesAbs::Vector{ModeShape{Float64}},normalizeModeShapes::Bool)

        # Initialize problem
        problem = new(model,x,Δx,residual,jacobian,inertia,jacobianDeterminant,timeNow,getLinearSolution,systemSolver,σ,getExternalForcesArray,savedσ,xOverσ,elementalStatesOverσ,nodalStatesOverσ,compElementalStatesOverσ,frequencies,dampings,eigenvectors,frequenciesFiltered,dampingsFiltered,eigenvectorsFiltered,frequenciesOscillatory,dampingsOscillatory,eigenvectorsOscillatoryCplx,eigenvectorsOscillatoryAbs,nModes,frequencyFilterLimits,modeShapesCplx,modeShapesAbs,normalizeModeShapes)

        # Set initial elemental and nodal states
        set_initial_states!(problem)

        # Initialize system arrays with correct size
        initialize_system_arrays!(problem)

        # Update initial load factor
        problem.σ = systemSolver.initialLoadFactor

        return problem

    end
end
export EigenProblem


"""
@with_kw mutable struct DynamicProblem <: Problem

DynamicProblem composite type

Defines the problem of dynamic type

# Fields
-
"""
@with_kw mutable struct DynamicProblem <: Problem

    # Model
    model::Model = Model()
    # States, residual, Jacobian and inertia arrays
    x::Vector{Float64} = zeros(0)
    Δx::Vector{Float64} = zeros(0)
    residual::Vector{Float64} = zeros(0)
    jacobian::Matrix{Float64} = zeros(0,0)
    inertia::Matrix{Float64} = zeros(0,0)
    jacobianDeterminant::Float64 = 0.0
    # Time variables
    initialTime::Number = 0.0
    Δt::Number 
    finalTime::Number
    timeVector::Union{Nothing,Vector{Float64}} = nothing
    timeNow::Number = initialTime
    timeBeginTimeStep::Number = timeNow
    timeEndTimeStep::Number = Δt
    indexBeginTimeStep::Int64 = 1
    indexEndTimeStep::Int64 = 2
    sizeOfTime::Int64 = 0
    # TF to get linear solution
    getLinearSolution::Bool = false
    # System solver
    systemSolver::SystemSolver = NewtonRaphson()
    # Load factor
    σ::Float64 = 1.0
    # TF to compute only the external forces array at the current nonlinear step (used only for arclength system solver)
    getExternalForcesArray::Bool = false
    # Initial velocities update options
    initialVelocitiesUpdateOptions::InitialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions()
    # TF to track partial solutions at time steps and its frequency
    trackingTimeSteps::Bool = true
    trackingFrequency::Int64 = 1
    # TF to display progress and its frequency
    displayProgress::Bool = true
    displayFrequency::Int64 = 0
    # Arrays of saved time steps, respective solutions and states
    savedTimeVector::Vector{Float64} = Vector{Float64}()
    xOverTime::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    elementalStatesOverTime::Vector{Vector{ElementalStates{Float64}}} = Vector{Vector{ElementalStates{Float64}}}()
    nodalStatesOverTime::Vector{Vector{NodalStates{Float64}}} = Vector{Vector{NodalStates{Float64}}}()
    compElementalStatesOverTime::Vector{Vector{ComplementaryElementalStates{Float64}}} = Vector{Vector{ComplementaryElementalStates{Float64}}}()
    elementalStatesRatesOverTime::Vector{Vector{ElementalStatesRates}} = Vector{Vector{ElementalStatesRates}}()
    compElementalStatesRatesOverTime::Vector{Vector{ComplementaryElementalStatesRates}} = Vector{Vector{ComplementaryElementalStatesRates}}()

    # Constructor
    function DynamicProblem(model::Model,x::Vector{Float64},Δx::Vector{Float64},residual::Vector{Float64},jacobian::Matrix{Float64},inertia::Matrix{Float64},jacobianDeterminant::Float64,initialTime::Number,Δt::Number,finalTime::Number,timeVector::Union{Nothing,Vector{Float64}},timeNow::Number,timeBeginTimeStep::Number,timeEndTimeStep::Number,indexBeginTimeStep::Int64,indexEndTimeStep::Int64,sizeOfTime::Int64,getLinearSolution::Bool,systemSolver::SystemSolver,σ::Float64,getExternalForcesArray::Bool,initialVelocitiesUpdateOptions::InitialVelocitiesUpdateOptions,trackingTimeSteps::Bool,trackingFrequency::Int64,displayProgress::Bool,displayFrequency::Int64,savedTimeVector::Vector{Float64},xOverTime::Vector{Vector{Float64}},elementalStatesOverTime::Vector{Vector{ElementalStates{Float64}}},nodalStatesOverTime::Vector{Vector{NodalStates{Float64}}},compElementalStatesOverTime::Vector{Vector{ComplementaryElementalStates{Float64}}},elementalStatesRatesOverTime::Vector{Vector{ElementalStatesRates}},compElementalStatesRatesOverTime::Vector{Vector{ComplementaryElementalStatesRates}})

        # Initialize problem
        problem = new(model,x,Δx,residual,jacobian,inertia,jacobianDeterminant,initialTime,Δt,finalTime,timeVector,timeNow,timeBeginTimeStep,timeEndTimeStep,indexBeginTimeStep,indexEndTimeStep,sizeOfTime,getLinearSolution,systemSolver,σ,getExternalForcesArray,initialVelocitiesUpdateOptions,trackingTimeSteps,trackingFrequency,displayProgress,displayFrequency,savedTimeVector,xOverTime,elementalStatesOverTime,nodalStatesOverTime,compElementalStatesOverTime,elementalStatesRatesOverTime,compElementalStatesRatesOverTime)

        # Set initial elemental and nodal states
        set_initial_states!(problem)   

        # Check and initialize time variables
        initialize_time_variables!(problem)

        # Initialize system arrays with correct size
        initialize_system_arrays!(problem)

        # Update initial load factor
        problem.σ = systemSolver.initialLoadFactor

        # Update display frequency if not input
        if displayProgress == true && displayFrequency == 0
            problem.displayFrequency = ceil(problem.sizeOfTime/100)
        end

        return problem

    end
end
export DynamicProblem


"""
set_initial_states!(problem::Problem,skipSizeAssertion::Bool=false)

Sets the initial elemental and nodal states into the states array

# Arguments
- problem::Problem
- skipSizeAssertion::Bool
"""
function set_initial_states!(problem::Problem,skipSizeAssertion::Bool=false)

    @unpack x,model = problem
    @unpack elements,specialNodes,systemOrder,nTrimVariables,forceScaling = model

    # Check if states array was input with correct size
    if !skipSizeAssertion && !isempty(x) 
        @assert length(x) == systemOrder+nTrimVariables
        return
    end

    # Initialize states array
    x = zeros(systemOrder+nTrimVariables)

    # Loop over elements and assign initial states
    for element in elements
        @unpack DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω = element
        @unpack u,p,F,M,V,Ω = element.states
        x[DOF_u] = u
        x[DOF_p] = p
        x[DOF_F] = F/forceScaling
        x[DOF_M] = M/forceScaling
        x[DOF_V] = V
        x[DOF_Ω] = Ω
    end

    # Loop over special nodes and assign initial states
    for specialNode in specialNodes
        @unpack connectedElements,ζonElements,BCs,DOF_uF,DOF_pM,DOF_trimLoads = specialNode
        # Get first connected element and the node's side on it
        firstConnectedElement = connectedElements[1]
        sideOnElement = ζonElements[1] == -1 ? 1 : 2
        # Set generalized displacements/forces values accordingly
        if sideOnElement == 1
            x[DOF_uF] = firstConnectedElement.nodalStates.u_n1
            x[DOF_pM] = firstConnectedElement.nodalStates.p_n1
        else
            x[DOF_uF] = firstConnectedElement.nodalStates.u_n2
            x[DOF_pM] = firstConnectedElement.nodalStates.p_n2
        end
        # If there are trim loads
        if any(x->x>0,DOF_trimLoads)
            # Loop BCs
            for BC in BCs
                @unpack isTrim,currentValue = BC
                # Set generalized trim forces values
                x[DOF_trimLoads[isTrim]] .= currentValue[isTrim]/forceScaling
            end
        end
    end

    @pack! problem = x

end


"""
initialize_system_arrays!(problem::Problem)

Initializes the system arrays to the correct size

# Arguments
- problem::Problem
"""
function initialize_system_arrays!(problem::Problem)

    @unpack model = problem
    @unpack systemOrder,nTrimVariables = model

    Δx = zeros(systemOrder+nTrimVariables)
    residual = zeros(systemOrder)
    jacobian = zeros(systemOrder,systemOrder+nTrimVariables)
    inertia = zeros(systemOrder,systemOrder)

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


"""
precompute_distributed_loads!(problem::Problem)

Pre-computes the distributed loads over all elements' nodes, over every time step

# Arguments
- problem::Problem
"""
function precompute_distributed_loads!(problem::Problem)


    @unpack model = problem

    # Unpack time variables for dynamic problem or define dummy ones otherwise
    if typeof(problem) == DynamicProblem
        @unpack timeVector,sizeOfTime = problem
    else
        sizeOfTime = 1
        timeVector = [0.0]
    end

    # Define shape function over element domain
    ϕ(node,ζ) = (1 - ζ) * (node == 1) + ζ * (node == 2)

    # Loop over beams
    for beam in model.beams
        # Loop over elements
        for element in beam.elements
            # Unpack
            @unpack Δℓ,f_A_of_ζt,m_A_of_ζt,f_b_of_ζt,m_b_of_ζt,ff_A_of_ζt,mf_A_of_ζt,ff_b_of_ζt,mf_b_of_ζt,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb = element
            # Initialize array of distributed loads on element's nodes over time
            f_A = zeros(Float64,2,3,sizeOfTime)
            m_A = zeros(Float64,2,3,sizeOfTime)
            f_b = zeros(Float64,2,3,sizeOfTime)
            m_b = zeros(Float64,2,3,sizeOfTime)
            ff_A = zeros(Float64,2,3,sizeOfTime)
            mf_A = zeros(Float64,2,3,sizeOfTime)
            ff_b = zeros(Float64,2,3,sizeOfTime)
            mf_b = zeros(Float64,2,3,sizeOfTime)
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


"""
solve_steady!(problem::Problem)

Solves a steady problem (includes trim and the steady part of eigen problems) 

# Arguments
- problem::Problem
"""
function solve_steady!(problem::Problem)

    @unpack getLinearSolution,systemSolver = problem

    if getLinearSolution
        problem.σ = 1.0
        assemble_system_arrays!(problem)
        solve_linear_system!(problem)
        update_states!(problem)
        save_load_factor_data!(problem,problem.σ,problem.x)
    else
        if typeof(systemSolver) == NewtonRaphson
            solve_NewtonRaphson!(problem)
        end
    end

end


"""
solve_eigen!(problem::Problem)

Solves an eigenproblem 

# Arguments
- problem::Problem
"""
function solve_eigen!(problem::Problem)

    @unpack jacobianDeterminant,jacobian,inertia,frequencyFilterLimits,nModes = problem
       
    ## Process eigenvalues and eigenvectors
    #---------------------------------------------------------------------------
    # Solve eigenproblem to get eigenvectors and inverse of eigenvalues
    if isapprox(jacobianDeterminant,0)
        inverseEigenvalues, eigenvectors = eigen(-pinv(jacobian)*inertia)
    else
        inverseEigenvalues, eigenvectors = eigen(-jacobian\inertia)
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

    ## Get oscillatory (structural, flight dynamics) modes
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
    if availableNModes > nModes
        frequenciesOscillatory = frequenciesOscillatory[1:nModes]
        dampingsOscillatory = dampingsOscillatory[1:nModes]
        eigenvectorsOscillatoryCplx = eigenvectorsOscillatoryCplx[:,1:nModes]
    else
        display("Only $availableNModes modes could be calculated, due to limited number of elements")
    end
    
    ## Post-process oscillatory eigenvectors
    #---------------------------------------------------------------------------
    # Real and imaginary parts of eigenvector
    eigenvectorsOscillatoryRe = real.(eigenvectorsOscillatoryCplx)
    eigenvectorsOscillatoryIm = imag.(eigenvectorsOscillatoryCplx)
    # Indices whose real part is larger
    realLarger = @. abs(eigenvectorsOscillatoryRe) > abs(eigenvectorsOscillatoryIm)   
    # Indices whose imaginary part is larger
    imagLarger = @. !realLarger                                       
    # Sign of larger part determines true sign of solution
    signs = @. sign(eigenvectorsOscillatoryRe) * realLarger + sign(eigenvectorsOscillatoryIm) * imagLarger                    
    # Correctly signed absolute value of eigenvectors
    eigenvectorsOscillatoryAbs = @. signs * abs(eigenvectorsOscillatoryCplx)   
    # Get mode shapes of complex-valued and absolute-valued eigenvectors
    get_mode_shapes!(problem,eigenvectorsOscillatoryCplx,frequenciesOscillatory,dampingsOscillatory)
    get_mode_shapes!(problem,eigenvectorsOscillatoryAbs,frequenciesOscillatory,dampingsOscillatory)

    @pack! problem = frequencies,dampings,eigenvectors,frequenciesFiltered,dampingsFiltered,eigenvectorsFiltered,frequenciesOscillatory,dampingsOscillatory,eigenvectorsOscillatoryCplx,eigenvectorsOscillatoryAbs

end


"""
get_mode_shapes!(problem::Problem,eigenvectorsOscillatory::Matrix{T},frequenciesOscillatory::Vector{Float64},dampingsOscillatory::Vector{Float64})

Gets the mode shapes given by the eigenvectors

# Arguments
- problem::Problem
- eigenvectorsOscillatory::Matrix{T}
- frequenciesOscillatory::Vector{Float64}
- dampingsOscillatory::Vector{Float64}
"""
function get_mode_shapes!(problem::Problem,eigenvectorsOscillatory::Matrix{T},frequenciesOscillatory::Vector{Float64},dampingsOscillatory::Vector{Float64}) where T<:Union{Float64,ComplexF64}

    @unpack model,normalizeModeShapes = problem
    @unpack forceScaling,elements,nElementsTotal = model

    if T === ComplexF64
        @unpack modeShapesCplx = problem
    else
        @unpack modeShapesAbs = problem
    end

    # Loop modes
    for (mode,freq,damp) in zip(1:size(eigenvectorsOscillatory,2),frequenciesOscillatory,dampingsOscillatory)
        # Initialize current mode shape
        modeShape = ModeShape{T}(mode,freq,damp,nElementsTotal)
        # Current mode's eigenvector
        eigenvector = eigenvectorsOscillatory[:,mode]
        # Initialize maxima
        displacementMax,rotationMax,forceMax,momentMax,velocityMax,strainMax,momentumMax,angleMax = 0,0,0,0,0,0,0,0
        # Loop elements
        for (e,element) in enumerate(elements)
            # Get modal states for current element
            u,p,F,M,V,Ω,γ,κ,P,H,u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2 = element_modal_states(element,eigenvector,forceScaling)
            # Update maxima
            displacementMax = maximum(abs.(vcat([displacementMax,u_n1,u_n2]...)))
            rotationMax = maximum(abs.(vcat([rotationMax,p_n1,p_n2]...)))
            forceMax = maximum(abs.(vcat([forceMax,F_n1,F_n2]...)))
            momentMax = maximum(abs.(vcat([momentMax,M_n1,M_n2]...)))
            velocityMax = maximum(abs.(vcat([velocityMax,V,Ω]...)))
            strainMax = maximum(abs.(vcat([strainMax,γ,κ]...)))
            momentumMax = maximum(abs.(vcat([momentumMax,P,H]...)))
            angleMax = maximum(abs.(vcat([angleMax,θ_n1,θ_n2]...)))
            # Add states to current mode 
            @pack! modeShape.elementalStates[e] = u,p,F,M,V,Ω
            @pack! modeShape.complementaryElementalStates[e] = γ,κ,P,H
            @pack! modeShape.nodalStates[e] = u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2
        end
        # Normalize mode shapes by maxima, if applicable
        if normalizeModeShapes
            # Loop elements 
            for (e,element) in enumerate(elements)
                @unpack u,p,F,M,V,Ω = modeShape.elementalStates[e]
                @unpack γ,κ,P,H = modeShape.complementaryElementalStates[e]
                @unpack u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2 = modeShape.nodalStates[e] 
                u,u_n1,u_n2,u_n1_b,u_n2_b = divide_inplace(displacementMax, u,u_n1,u_n2,u_n1_b,u_n2_b) 
                p,p_n1,p_n2,p_n1_b,p_n2_b = divide_inplace(rotationMax, p,p_n1,p_n2,p_n1_b,p_n2_b) 
                F,F_n1,F_n2 = divide_inplace(forceMax, F,F_n1,F_n2)
                M,M_n1,M_n2 = divide_inplace(momentMax, M,M_n1,M_n2)
                V,Ω = divide_inplace(velocityMax, V,Ω)
                γ,κ = divide_inplace(strainMax, γ,κ)
                P,H = divide_inplace(momentumMax, P,H)
                θ_n1,θ_n2 = divide_inplace(angleMax, θ_n1,θ_n2)
                @pack! modeShape.elementalStates[e] = u,p,F,M,V,Ω
                @pack! modeShape.complementaryElementalStates[e] = γ,κ,P,H
                @pack! modeShape.nodalStates[e] = u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,F_n1,F_n2,M_n1,M_n2,θ_n1,θ_n2
            end
        end
        # Add to mode shapes array
        if T === ComplexF64
            push!(modeShapesCplx,modeShape)
        else
            push!(modeShapesAbs,modeShape)
        end
    end

    if T === ComplexF64
        @pack! problem = modeShapesCplx
    else
        @pack! problem = modeShapesAbs
    end
end


"""
solve_dynamic!(problem::Problem)

Solves a dynamic problem 

# Arguments
- problem::Problem
"""
function solve_dynamic!(problem::Problem)

    # Solve the initial time step, to get consistent initial states
    solve_initial_dynamic!(problem)
    # Time march
    time_march!(problem)

end


"""
initialize_time_variables!(problem::Problem)

Validates and initializes the time variables

# Arguments
- problem::Problem
"""
function initialize_time_variables!(problem::Problem)

    @unpack initialTime,finalTime,Δt,timeVector,initialVelocitiesUpdateOptions = problem

    # Check final and initial time inputs
    @assert finalTime >= initialTime + Δt

    # Assign and validate time vector
    if isnothing(timeVector)
        timeVector = collect(initialTime:Δt:finalTime)
    end
    @assert timeVector == sort!(copy(timeVector))

    # Assign other time variables
    timeNow = timeVector[1]
    timeBeginTimeStep = timeVector[1]
    timeEndTimeStep = timeVector[2]
    sizeOfTime = length(timeVector)

    # Update time step for initial velocities update, if not input
    if initialVelocitiesUpdateOptions.Δt == 0
        initialVelocitiesUpdateOptions.Δt = 0.1*Δt
    end

    @pack! problem = timeVector,timeNow,timeBeginTimeStep,timeEndTimeStep,sizeOfTime,initialVelocitiesUpdateOptions

end


"""
solve_initial_dynamic!(problem::Problem)

Gets consistent initial conditions and solves the first time step

# Arguments
- problem::Problem
"""
function solve_initial_dynamic!(problem::Problem)

    problemCopy = deepcopy(problem)
    @unpack model,trackingTimeSteps = problemCopy
    @unpack BCs,BCedNodes = model
    
    # Initialize generalized displacements BCs 
    initialDisplacementsBCs = BCs

    # Initialize flag for node having being assigned BCs
    assigned = falses(model.nNodesTotal)
    assigned[BCedNodes] .= true

    # List of initial generalized displacements types (defined in basis b)
    displacementTypes = ["u1b","u2b","u3b","p1b","p2b","p3b"]

    # Loop elements
    for (e, element) in enumerate(model.elements)
        @unpack u_n1,u_n2,p_n1,p_n2 = element.nodalStates
        # Loop element's nodes
        for (n, localNode) in enumerate(element.nodesLocalID[1]:element.nodesLocalID[2])
            globalNode = element.nodesGlobalID[n]
            if !assigned[globalNode]
                # Update flag
                assigned[globalNode] = true
                # Initial nodal generalized displacements
                if n == 1 
                    initialNodalDisplacements = vcat(u_n1,p_n1)
                else
                    initialNodalDisplacements = vcat(u_n2,p_n2)
                end
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
    set_BCs!(model,initialDisplacementsBCs)
    model.skipValidationMotionBasisA = true
    update_model!(model)
    @pack! problemCopy = model
 
    # Set initial states
    skipSizeAssertion = true
    set_initial_states!(problemCopy,skipSizeAssertion)

    # Initialize system arrays with correct size
    initialize_system_arrays!(problemCopy)

    # Update Δt
    problemCopy.Δt = Inf64

    # Get equivalent initial states rates
    get_equivalent_states_rates!(problemCopy)

    # Solve to find unspecified initial states
    solve_steady!(problemCopy)

    # Copy initial states to the original problem
    copy_initial_states!(problem,problemCopy)

    # Update initial velocities states
    update_initial_velocities!(problem)

    # Save time step data, if applicable    
    if trackingTimeSteps 
        save_time_step_data!(problem,problem.timeNow)           
    end
    
end


"""
copy_initial_states!(problem::Problem,problemCopy::Problem)

Copies the initial states to the original problem

# Arguments
- problem::Problem
- problemCopy::Problem
"""
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
        x[DOF_F] = xCopy[DOF_FCopy]
        x[DOF_M] = xCopy[DOF_MCopy]
        V = x[DOF_V] = xCopy[DOF_VCopy]
        Ω = x[DOF_Ω] = xCopy[DOF_ΩCopy]
        F,M = forceScaling*x[DOF_F],forceScaling*x[DOF_M]
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


"""
update_initial_velocities!(problem::Problem)

Updates the velocity states by running a very small time step

# Arguments
- problem::Problem
"""
function update_initial_velocities!(problem::Problem)

    @unpack timeVector,Δt,timeNow,indexBeginTimeStep,indexEndTimeStep,timeBeginTimeStep,timeEndTimeStep,initialVelocitiesUpdateOptions,model = problem
    @unpack maxIter,relaxFactor,tol,displayProgress = initialVelocitiesUpdateOptions
    @unpack elements,nElementsTotal = model

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
    update_basis_A_orientation!(problem)
    # Update boundary conditions
    for BC in model.BCs
        update_BC_data!(BC,timeNow)
    end   
    # Set default system solver (with default convergence tolerances and large # of maximum iterations)
    problem.systemSolver = NewtonRaphson(maximumIterations=100)
    # Initialize velocity arrays at begin and end of time step
    velInitial = zeros(6*nElementsTotal)
    velFinal = zeros(6*nElementsTotal)
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
                @unpack u,V,Ω = element.states
                @unpack udot,pdot = element.statesRates
                # Indices to update
                VToUpdate = DoFToUpdate[e][1:3]
                ΩToUpdate = DoFToUpdate[e][4:6]
                # Multiply by relaxation factor, if applicable
                udot[VToUpdate] = udot[VToUpdate]*relaxFactor
                pdot[ΩToUpdate] = pdot[ΩToUpdate]*relaxFactor
                V[VToUpdate] = V[VToUpdate]*relaxFactor
                Ω[ΩToUpdate] = Ω[ΩToUpdate]*relaxFactor
                # Pack
                @pack! element.states = V,Ω
                @pack! element.statesRates = udot,pdot
            end
        end
        # Get velocities array at begin of time step 
        velInitial = vcat([[e.states.V; e.states.Ω] for e in elements]...)
        # Reset accelerations to zero
        for element in elements
            Vdot,Ωdot,uddot,pddot = zeros(3),zeros(3),zeros(3),zeros(3)
            @pack! element.statesRates = Vdot,Ωdot,uddot,pddot
        end
        # Get equivalent states' rates at the begin of the time step
        get_equivalent_states_rates!(problem)
        # Solve the system at the current time step
        solve_time_step!(problem)
        # Get velocities array at end of time step 
        velFinal = vcat([[e.states.V; e.states.Ω] for e in elements]...)
        # Update current difference between initial and final velocities 
        ϵ = maximum(abs.(velInitial-velFinal) .* vcat(DoFToUpdate...))
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
    timeNow = timeVector[1] + Δt
    indexBeginTimeStep = 1
    indexEndTimeStep = 2
    timeBeginTimeStep = timeNow - Δt
    timeEndTimeStep = timeNow 
    problem.systemSolver = systemSolverOriginal

    @pack! problem = timeVector,Δt,timeNow,indexBeginTimeStep,indexEndTimeStep,timeBeginTimeStep,timeEndTimeStep,model
end


"""
time_march!(problem::Problem)

Marches the dynamic problem in time

# Arguments
- problem::Problem
"""
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
        # Solve the system at the current time step
        solve_time_step!(problem)
        # Stop if unconverged
        if !problem.systemSolver.convergedFinalSolution
            break
        end
        # Save time step data, if applicable    
        if trackingTimeSteps && rem(timeIndex,trackingFrequency) == 0
            save_time_step_data!(problem,timeNow)           
        end
        # Display progress, if applicable
        if displayProgress && rem(timeIndex,displayFrequency) == 0
            progress = round(timeIndex/sizeOfTime*100,digits=1)
            println("Simulation progress: $progress %")
        end
    end
end


"""
update_time_variables!(problem::Problem,timeIndex::Int64)

Updates the time variables (time, time step, time indices)

# Arguments
- problem::Problem
- timeIndex::Int64
"""
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


"""
update_basis_A_orientation!(problem::Problem)

Updates the orientation of the basis A for the next time step

# Arguments
- problem::Problem
"""
function update_basis_A_orientation!(problem::Problem)

    @unpack model,timeNow,Δt = problem
    @unpack R_A,ω_A = model

    # Skew-symmetric operator of the angular velocity vector
    ω_A_tilde = tilde(ω_A(timeNow))

    # Rotation tensor from basis I to basis A, at the current time, and its transpose
    R_A = inv((2/Δt*I3-ω_A_tilde))*(2/Δt*I3+ω_A_tilde)*R_A
    round_off!(R_A)
    R_AT = Matrix(R_A')

    @pack! model = R_A,R_AT

end


"""
get_equivalent_states_rates!(problem::Problem)

Gets the equivalent states' rates at the begin of the current time step 

# Arguments
- problem::Problem
"""
function get_equivalent_states_rates!(problem::Problem)

    @unpack model,Δt = problem
    @unpack elements = model

    # Loop over elements
    for element in elements
        # Unpack element data (element states and rates known at the begin of time step)
        @unpack u,p,V,Ω = element.states
        @unpack udot,pdot,Vdot,Ωdot,uddot,pddot = element.statesRates
        # Equivalent states' rates at the begin of time step 
        udotEquiv = udot + 2/Δt*u
        pdotEquiv = pdot + 2/Δt*p
        VdotEquiv = Vdot + 2/Δt*V
        ΩdotEquiv = Ωdot + 2/Δt*Ω
        uddotEquiv = uddot + 2/Δt*udot
        pddotEquiv = pddot + 2/Δt*pdot
        # Pack element data
        @pack! element = udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,uddotEquiv,pddotEquiv
    end
end


"""
solve_time_step!(problem::Problem)

Solves the current time step  

# Arguments
- problem::Problem
"""
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


"""
save_time_step_data!(problem::Problem,timeNow::Number)

Saves the solution at the current time step

# Arguments
- problem::Problem
"""
function save_time_step_data!(problem::Problem,timeNow::Number)

    @unpack x,savedTimeVector,xOverTime,elementalStatesOverTime,nodalStatesOverTime,compElementalStatesOverTime,elementalStatesRatesOverTime,compElementalStatesRatesOverTime,model = problem

    # Add current time
    push!(savedTimeVector,timeNow)

    # Add curent system states 
    push!(xOverTime,x)

    # Add current elemental states 
    currentElementalStates = Vector{ElementalStates}()
    for element in model.elements
        push!(currentElementalStates,deepcopy(element.states))
    end
    push!(elementalStatesOverTime,currentElementalStates)

    # Add current nodal states 
    currentNodalStates = Vector{NodalStates}()
    for element in model.elements
        push!(currentNodalStates,deepcopy(element.nodalStates))
    end
    push!(nodalStatesOverTime,currentNodalStates)

    # Add current complementary elemental states 
    currentComplementaryElementalStates = Vector{ComplementaryElementalStates}()
    for element in model.elements
        push!(currentComplementaryElementalStates,deepcopy(element.compStates))
    end
    push!(compElementalStatesOverTime,currentComplementaryElementalStates)

    # Add current states' rates
    currentStatesRates = Vector{ElementalStatesRates}()
    for element in model.elements
        push!(currentStatesRates,deepcopy(element.statesRates))
    end
    push!(elementalStatesRatesOverTime,currentStatesRates)

    # Add current complementary elemental states' rates
    currentComplementaryElementalStatesRates = Vector{ComplementaryElementalStatesRates}()
    for element in model.elements
        push!(currentComplementaryElementalStatesRates,deepcopy(element.compStatesRates))
    end
    push!(compElementalStatesRatesOverTime,currentComplementaryElementalStatesRates)


    @pack! problem = savedTimeVector,xOverTime,elementalStatesOverTime,nodalStatesOverTime,compElementalStatesOverTime,elementalStatesRatesOverTime,compElementalStatesRatesOverTime

end