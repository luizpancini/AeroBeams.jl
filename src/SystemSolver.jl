"""
@with_kw mutable struct NewtonRaphson <: SystemSolver

    Newton-Raphson system solver composite type

# Fields
- 
"""
@with_kw mutable struct NewtonRaphson <: SystemSolver

    # User input variables
    absoluteTolerance::Float64 = 1e-8
    relativeTolerance::Float64 = 1e-8
    maximumIterations::Int64 = 20
    desiredIterations::Int64 = 5
    maximumRelativeError::Float64 = 1e6
    initialLoadFactor::Number = 1.0
    minimumLoadFactor::Float64 = 0.01
    maximumLoadFactorStep::Float64 = 0.5
    minimumLoadFactorStep::Float64 = 0.01
    trackingLoadSteps::Bool = true
    displayStatus::Bool = false

    # Algorithm variables
    loadFactorStep::Float64 = max(min(maximumLoadFactorStep,0.5),minimumLoadFactorStep)
    convergedPartialSolution::Bool = false
    divergedPartialSolution::Bool = false
    convergedFinalSolution::Bool = false

end

# Constructor
function create_NewtonRaphson(;absoluteTolerance::Float64=1e-8,relativeTolerance::Float64=1e-8,maximumIterations::Int64=20,desiredIterations::Int64=5,maximumRelativeError::Float64=1e6,initialLoadFactor::Number=1.0,minimumLoadFactor::Float64=0.01,maximumLoadFactorStep::Float64=0.5,minimumLoadFactorStep::Float64=0.01,trackingLoadSteps::Bool=true,displayStatus::Bool=false)

    return NewtonRaphson(absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance,maximumIterations=maximumIterations,desiredIterations=desiredIterations,maximumRelativeError=maximumRelativeError,initialLoadFactor=initialLoadFactor,minimumLoadFactor=minimumLoadFactor,maximumLoadFactorStep=maximumLoadFactorStep,minimumLoadFactorStep=minimumLoadFactorStep,trackingLoadSteps=trackingLoadSteps,displayStatus=displayStatus)

end
export create_NewtonRaphson


"""
solve_NewtonRaphson!(problem::Problem)

Solves the nonlinear system of equations at current time step 

# Arguments
- problem::Problem
"""
function solve_NewtonRaphson!(problem::Problem)

    # Unpack system solver
    @unpack absoluteTolerance,relativeTolerance,maximumIterations,desiredIterations,maximumRelativeError,minimumLoadFactor,loadFactorStep,maximumLoadFactorStep,minimumLoadFactorStep,trackingLoadSteps,displayStatus,convergedFinalSolution,convergedPartialSolution,divergedPartialSolution = problem.systemSolver

    # Reset converged solution flag and load factor step
    convergedFinalSolution = false
    loadFactorStep = max(min(maximumLoadFactorStep,0.5),minimumLoadFactorStep)

    ## Load steps
    #---------------------------------------------------------------------------
    loadstep = 0
    while !convergedFinalSolution   
        @unpack x,residual,σ = problem
        # Update loadstep
        loadstep += 1 
        # Initialize iteration count and this loadstep's known states, residual and load factor
        iter = 0
        x_known = x
        residual_known = residual
        loadFactor_known = loadstep == 1 ? 0.0 : (convergedPartialSolution ? σ-loadFactorStep : σ)
        # Reset partial convergence and divergence flags
        convergedPartialSolution,divergedPartialSolution = false,false
        # Update load factor 
        σ = loadstep == 1 ? σ : min(1.0,σ+loadFactorStep)
        @pack! problem = σ
        ## Nonlinear solution procedure
        # ----------------------------------------------------------------------
        while !convergedPartialSolution       
            # Update iteration count
            iter += 1     
            # Assemble system matrices
            assemble_system_arrays!(problem)  
            # Stop if Jacobian becomes singular
            @unpack jacobianDeterminant = problem
            if isapprox(jacobianDeterminant,0)
                divergedPartialSolution = true
                println("Singular Jacobian, stopping...") 
                return
            end 
            # Solve the linear system (update states array)
            solve_linear_system!(problem)
            @unpack x,Δx,residual = problem
            # Calculate residual (absolute) and solution (relative) convergence norms
            ϵ_abs = norm(residual)
            ϵ_rel = norm(Δx)/norm(x)
            # Check convergence norms
            if ϵ_abs < absoluteTolerance || ϵ_rel < relativeTolerance
                convergedPartialSolution = true
            elseif ϵ_rel > maximumRelativeError || isnan(ϵ_abs)
                divergedPartialSolution = true 
                println("Residual too large...")
                break
            end            
            # Check iterations count
            if iter == maximumIterations
                break
            end
        end   
        ## Update load factor step
        #-----------------------------------------------------------------------
        # In case of convergence trouble
        if !convergedPartialSolution                            
            # Case of unconverged iterations not dealt yet from trim analyses
            if typeof(problem) == TrimProblem
                error("Trim analysis unconverged")
            end
            # Reduce load factor (limit to minimum) 
            σ = max(σ/2,minimumLoadFactor)
            # Reduce load factor step (limit to minimum)
            if divergedPartialSolution
                loadFactorStep = loadFactorStep/2
            else
                loadFactorStep = loadFactorStep*sqrt(desiredIterations/iter)
            end
            loadFactorStep = max(min(1.0-loadFactor_known,loadFactorStep),minimumLoadFactorStep)
            # Check if algorithm is stuck
            if loadFactorStep == minimumLoadFactorStep 
                println("NR algorithm stuck, stopping...")
                return
            end
            # Reset to previously converged solution
            x = x_known 
            residual = residual_known
            σ = loadFactor_known
            println("Unconverged, reducing load factor and trying again...")
            @pack! problem = x,residual,σ 
        # In case of converged partial solution    
        else
            # Update states on element level
            update_states!(problem)
            # Adjust load factor step (and impose bounds)
            loadFactorStep = max(minimumLoadFactorStep,min(maximumLoadFactorStep, loadFactorStep*sqrt(desiredIterations/iter)))
            # Save solution at current load factor, if tracking
            if trackingLoadSteps && typeof(problem) in [SteadyProblem,TrimProblem,EigenProblem]
                save_load_factor_data!(problem,σ,x)
            end
        end
        # Display status
        if displayStatus
            loadFactorPercent = σ*100
            println("Load factor at $loadFactorPercent %, after $iter iterations")
        end
        # Check for full load and convergence reached
        if abs(σ-1) < 1e-6 && convergedPartialSolution 
            convergedFinalSolution = true
        end
    end    

    @pack! problem.systemSolver = convergedPartialSolution,divergedPartialSolution,convergedFinalSolution,loadFactorStep
end


"""
assemble_system_arrays!(problem::Problem)

Assembles the residual vector, Jacobian and inertia matrices of the system of equations

# Arguments
- problem::Problem
"""
function assemble_system_arrays!(problem::Problem)

    @unpack model = problem
    @unpack elements,specialNodes,systemOrder,nTrimVariables = model

    # Reset Jacobian and inertia matrices (only for debugging purposes)
    problem.jacobian = zeros(systemOrder,systemOrder+nTrimVariables)
    problem.inertia = zeros(systemOrder,systemOrder)

    # Get contributions from the elements
    for element in elements
        element_arrays!(problem,model,element)
    end

    # Get contributions from the special nodes
    for specialNode in specialNodes
        special_node_arrays!(problem,model,specialNode)
    end

    # Determinant of Jacobian
    jacobianDeterminant = typeof(problem) != TrimProblem ? det(problem.jacobian) : NaN

    @pack! problem = jacobianDeterminant 

end


"""
solve_linear_system!(problem::Problem)

Solves the linear system of equations at current time step and load factor

# Arguments
- problem::Problem
"""
function solve_linear_system!(problem::Problem)

    @unpack x,Δx,residual,jacobian,jacobianDeterminant = problem

    # Check Jacobian's determinant
    if isapprox(jacobianDeterminant,0)
        println("Singular Jacobian: Δx not updated")
        return 
    end

    # Solve the linear system according to problem type
    if typeof(problem) != TrimProblem
        try
            Δx = -sparse(jacobian)\residual
        catch
            Δx = -jacobian\residual
        end
    else
        Δx = -pinv(jacobian)*residual
    end

    # Update states array
    x += Δx

    @pack! problem = x,Δx

end


"""
save_load_factor_data!(problem::Problem,σ::Float64,x::Vector{Float64})

Saves the solution at the current load factor

# Arguments
- problem::Problem
- σ::Float64
- x::Vector{Float64}
"""
function save_load_factor_data!(problem::Problem,σ::Float64,x::Vector{Float64})


    @unpack savedσ,xOverσ,elementalStatesOverσ,nodalStatesOverσ,compElementalStatesOverσ,model = problem    

    # Add current load factor
    push!(savedσ,σ)

    # Add curent system states 
    push!(xOverσ,x)

    # Add current elemental states 
    currentElementalStates = Vector{ElementalStates}()
    for element in model.elements
        push!(currentElementalStates,deepcopy(element.states))
    end
    push!(elementalStatesOverσ,currentElementalStates)

    # Add current nodal states 
    currentNodalStates = Vector{NodalStates}()
    for element in model.elements
        push!(currentNodalStates,deepcopy(element.nodalStates))
    end
    push!(nodalStatesOverσ,currentNodalStates)

    # Add current complementary elemental states 
    currentComplementaryElementalStates = Vector{ComplementaryElementalStates}()
    for element in model.elements
        push!(currentComplementaryElementalStates,deepcopy(element.compStates))
    end
    push!(compElementalStatesOverσ,currentComplementaryElementalStates)

    @pack! problem = savedσ,xOverσ,elementalStatesOverσ,nodalStatesOverσ,compElementalStatesOverσ

end