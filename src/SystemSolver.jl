#
# @with_kw mutable struct NewtonRaphson <: SystemSolver

#     Newton-Raphson composite type

#
@with_kw mutable struct NewtonRaphson <: SystemSolver

    # User input variables
    absoluteTolerance::Float64
    relativeTolerance::Float64
    maximumIterations::Int64
    desiredIterations::Int64
    maximumAbsoluteError::Real
    maximumRelativeError::Real
    initialLoadFactor::Real
    minimumLoadFactor::Float64
    maximumLoadFactorStep::Float64
    minimumLoadFactorStep::Float64
    ρ::Real
    trackingLoadSteps::Bool
    displayStatus::Bool
    minConvRateAeroJacUpdate::Real
    minConvRateJacUpdate::Real
    alwaysUpdateJacobian::Bool
    ΔλRelaxFactor::Real

    # Algorithm variables
    loadFactorStep::Float64 = max(min(maximumLoadFactorStep,0.5),minimumLoadFactorStep)
    convergedPartialSolution::Bool = false
    convergedFinalSolution::Bool = false
    convRate::Real = 0

end


"""
    create_NewtonRaphson(; kwargs...)

Newton-Raphson nonlinear system solver constructor

# Keyword arguments
- `absoluteTolerance::Float64` = absolute convergence tolerance
- `relativeTolerance::Float64` = relative convergence tolerance
- `maximumIterations::Int64` = maximum number of iterations
- `desiredIterations::Int64` = desired number of iterations
- `maximumAbsoluteError::Real` = maximum absolute error for divergence detection
- `maximumRelativeError::Real` = maximum relative error for divergence detection
- `initialLoadFactor::Real` = initial load factor
- `minimumLoadFactor::Float64` = minimum load factor
- `maximumLoadFactorStep::Float64` = maximum load factor step
- `minimumLoadFactorStep::Float64` = minimum load factor step
- `ρ::Real` = relaxation factor for trim variables
- `trackingLoadSteps::Bool` = flag to track partial load steps solutions
- `displayStatus::Bool` = flag to display status
- `minConvRateAeroJacUpdate::Real` = minimum convergence rate to skip computation of aerodynamic Jacobians
- `minConvRateJacUpdate::Real` = minimum convergence rate to skip computation of structural Jacobians
- `alwaysUpdateJacobian::Bool` = flag to update Jacobians on every iteration
- `ΔλRelaxFactor` = relaxation factor for update of Lagrange multipliers
"""
function create_NewtonRaphson(; absoluteTolerance::Float64=1e-8,relativeTolerance::Float64=1e-8,maximumIterations::Int64=20,desiredIterations::Int64=5,maximumAbsoluteError::Real=1e6,maximumRelativeError::Real=1e6,initialLoadFactor::Real=1.0,minimumLoadFactor::Float64=0.01,maximumLoadFactorStep::Float64=0.5,minimumLoadFactorStep::Float64=0.01,ρ::Real=1.0,trackingLoadSteps::Bool=true,displayStatus::Bool=false,minConvRateAeroJacUpdate::Real=2.0,minConvRateJacUpdate::Real=2.0,alwaysUpdateJacobian::Bool=true,ΔλRelaxFactor::Real=1)

    @assert maximumIterations > 1 "maximumIterations must be greater than 1"
    @assert desiredIterations > 1 "desiredIterations must be greater than 1"
    @assert 0 <= initialLoadFactor <= 1 "initialLoadFactor must be between 0 and 1"
    @assert 0 < minimumLoadFactor <= 0.5 "minimumLoadFactor must be between 0 and 0.5"
    @assert 0 < maximumLoadFactorStep <= 1 "maximumLoadFactorStep must be between 0 and 1"
    @assert 0 < minimumLoadFactorStep <= 1 "minimumLoadFactorStep must be between 0 and 1"
    @assert 0.5 <= ρ <= 1 "relaxation factor for trim variables (ρ) must be between 0.5 and 1.0"
    @assert 0 < ΔλRelaxFactor <= 1 "relaxation factor (ΔλRelaxFactor) for Lagrange multipliers update must be between 0 and 1"
    @assert minConvRateAeroJacUpdate > 1 "minConvRateAeroJacUpdate to skip calculation of aerodynamic derivatives must be greater than 1"
    @assert minConvRateJacUpdate > 1 "minConvRateJacUpdate to skip update of Jacobian matrix must be greater than 1"

    return NewtonRaphson(absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance,maximumIterations=maximumIterations,desiredIterations=desiredIterations,maximumAbsoluteError=maximumAbsoluteError,maximumRelativeError=maximumRelativeError,initialLoadFactor=initialLoadFactor,minimumLoadFactor=minimumLoadFactor,maximumLoadFactorStep=maximumLoadFactorStep,minimumLoadFactorStep=minimumLoadFactorStep,ρ=ρ,trackingLoadSteps=trackingLoadSteps,displayStatus=displayStatus,minConvRateAeroJacUpdate=minConvRateAeroJacUpdate,minConvRateJacUpdate=minConvRateJacUpdate,alwaysUpdateJacobian=alwaysUpdateJacobian,ΔλRelaxFactor=ΔλRelaxFactor)

end
export create_NewtonRaphson


# Solves the nonlinear system of equations at current time step
function solve_NewtonRaphson!(problem::Problem)

    # Unpack system solver
    @unpack absoluteTolerance,relativeTolerance,maximumIterations,desiredIterations,maximumAbsoluteError,maximumRelativeError,minimumLoadFactor,loadFactorStep,maximumLoadFactorStep,minimumLoadFactorStep,trackingLoadSteps,displayStatus,convergedFinalSolution,convergedPartialSolution,alwaysUpdateJacobian,minConvRateJacUpdate = problem.systemSolver

    # Reset converged solution flag and load factor step
    convergedPartialSolution = convergedFinalSolution = false
    loadFactorStep = max(min(maximumLoadFactorStep,0.5),minimumLoadFactorStep)
    @pack! problem.systemSolver = convergedPartialSolution,convergedFinalSolution,loadFactorStep

    ## Load steps
    #---------------------------------------------------------------------------
    loadstep = 0
    while !convergedFinalSolution   
        @unpack x,residual,σ = problem
        # Update loadstep
        loadstep += 1 
        # Initialize iteration count and this loadstep's known states, residual and load factor
        iter = 0
        xKnown = deepcopy(x)
        residualKnown = deepcopy(residual)
        BCsKnown = deepcopy(problem.model.BCs)
        specialNodesKnown = deepcopy(problem.model.specialNodes)
        loadFactorKnown = loadstep == 1 ? 0.0 : deepcopy(σ)
        ϵ_abs_previous = 1
        # Reset partial convergence flag
        convergedPartialSolution = false
        # Reset TF to skip Jacobian update
        problem.skipJacobianUpdate = false
        # Update load factor 
        σ = loadstep == 1 ? σ : max(loadFactorKnown, min(1.0,σ+loadFactorStep))
        @pack! problem = σ
        ## Nonlinear solution procedure
        # ----------------------------------------------------------------------
        while !convergedPartialSolution       
            # Update iteration count
            iter += 1     
            # Assemble system matrices
            assemble_system_arrays!(problem)  
            # Solve the linear system (and get convergence flag)
            converged = solve_linear_system!(problem)
            if !converged
                break
            end
            # Calculate residual (absolute) and solution (relative) convergence norms
            @unpack x,Δx,residual = problem
            ϵ_abs = norm(residual)
            ϵ_rel = norm(Δx)/norm(x)
            # Display status
            if displayStatus
                println("i: $iter, σ: $σ, ϵ_abs: $ϵ_abs, ϵ_rel: $ϵ_rel")
            end
            # Update converge rate
            convRate = ϵ_abs_previous/ϵ_abs
            @pack! problem.systemSolver = convRate
            # Update previous value of ϵ_abs
            ϵ_abs_previous = deepcopy(ϵ_abs)
            # Update TF to skip Jacobian update, if applicable
            problem.skipJacobianUpdate = (!alwaysUpdateJacobian && convRate > minConvRateJacUpdate) ? true : false
            # Check convergence norms
            if ϵ_abs < absoluteTolerance || ϵ_rel < relativeTolerance
                convergedPartialSolution = true
            elseif ϵ_abs > maximumAbsoluteError || ϵ_rel > maximumRelativeError || isnan(ϵ_abs)
                if displayStatus
                    println("Diverging solution...")
                end
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
            # Reduce load factor step (limit to minimum)
            loadFactorStep = max(min(1.0-loadFactorKnown,loadFactorStep/2),minimumLoadFactorStep)
            # Check if algorithm is stuck
            if loadFactorStep == minimumLoadFactorStep 
                println("NR algorithm stuck, stopping...")
                update_states!(problem)
                return
            end
            # Reset to previously known solution
            x = xKnown 
            residual = residualKnown
            σ = loadFactorKnown
            BCs = BCsKnown
            specialNodes = specialNodesKnown
            if displayStatus
                println("Unconverged, reducing load factor and trying again...")
            end
            @pack! problem = x,residual,σ 
            # @pack! problem.model = BCs,specialNodes
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
        # Check for full load and convergence reached
        if abs(σ-1) < 1e-6 && convergedPartialSolution
            # Update flag
            convergedFinalSolution = true
            # Update Jacobian matrix, if applicable
            if !alwaysUpdateJacobian 
                problem.skipJacobianUpdate = false
                for element in problem.model.elements
                    distributed_loads_derivatives_rotation_parameters!(element)
                    aero_derivatives!(problem,problem.model,element)
                    element_jacobian!(problem,problem.model,element)
                    element.aero = reset_dual_numbers(element.aero)
                end
            end
            # Get inertia matrix in eigen and trim problems, if applicable
            if problem isa EigenProblem || (problem isa TrimProblem && problem.getInertiaMatrix)
                for element in problem.model.elements
                    element_inertia!(problem,problem.model,element)
                end
            end
            # Update loads extrema
            update_maximum_aero_loads!(problem)
        end
    end    

    @pack! problem.systemSolver = convergedPartialSolution,convergedFinalSolution,loadFactorStep
end


# Assembles the residual vector and Jacobian matrix of the system of equations
function assemble_system_arrays!(problem::Problem,x::Vector{Float64}=problem.x)

    @unpack model = problem
    @unpack elements,specialNodes = model
    @pack! problem = x

    # Reset Jacobian matrix, if applicable
    if problem isa TrimProblem
        problem.jacobian .= 0
    end

    # Update states of the elements first (for better convergence with relative rotation constraints)
    for element in elements
        element_states!(problem,model,element)
    end

    # Get contributions from the elements
    for element in elements
        element_arrays!(problem,model,element)
    end

    # Update states of the special nodes first (for better convergence with doubly-attached springs)
    for specialNode in specialNodes
        special_node_states!(problem,model,specialNode)
    end

    # Get contributions from the special nodes
    for specialNode in specialNodes
        special_node_arrays!(problem,model,specialNode)
    end

    return problem.residual

end


# Solves the linear system of equations at current time step and load factor
function solve_linear_system!(problem::Problem)

    @unpack x,Δx,residual,jacobian,getLinearSolution,σ = problem
    @unpack systemOrder,nTrimVariables,hasRotationConstraints,hasHingeAxisConstraints,rotationConstraints,hingeAxisConstraints,rotationConstraintBalanceLoadBCid,hingeAxisConstraintBalanceLoadBCid,forceScaling = problem.model
    @unpack ρ,ΔλRelaxFactor = problem.systemSolver

    # Issue warning for singular Jacobian in linear problems
    if getLinearSolution && nTrimVariables == 0 && det(jacobian) ≈ 0
        println("Singular Jacobian")
    end

    # Solve the linear system according to problem type
    #---------------------------------------------------------------------------
    # Trim problem
    if problem isa TrimProblem
        Δx .= -pinv(Matrix(jacobian))*residual
    # Problem with rotation and/or hinge axis constraints    
    elseif hasRotationConstraints || hasHingeAxisConstraints
        Δx,Δλ = linear_solver_with_constraints(x,jacobian,residual,rotationConstraints,hingeAxisConstraints,problem)
    # Regular problem    
    else
        try
            Δx .= -jacobian\residual
        catch
            try
                Δx .= line_search(x,residual,jacobian)
            catch
                return false
            end
        end
    end

    # Update states array (apply relaxation factor on trim variables)
    x[1:systemOrder] .+= Δx[1:systemOrder]
    x[systemOrder+1:systemOrder+nTrimVariables] .+= ρ*Δx[systemOrder+1:systemOrder+nTrimVariables]

    # Update balance loads of rotation constraints with Lagrange multipliers increment
    for (i,constraint) in enumerate(rotationConstraints)
        # Set value
        problem.model.BCs[rotationConstraintBalanceLoadBCid[i]].values[1] = constraint.balanceMoment -= Δλ[i]*ΔλRelaxFactor/σ
        # Update BC
        update_BC_data!(problem.model.BCs[rotationConstraintBalanceLoadBCid[i]],problem.timeNow)
    end

    # Update balance loads of hinge axis constraints with Lagrange multipliers increment and rotation data across the hinge
    for (i,constraint) in enumerate(hingeAxisConstraints)
        # Set values
        nλ = isnothing(constraint.ΔpValue) ? 2 : 3
        problem.model.BCs[hingeAxisConstraintBalanceLoadBCid[i]].values[1:nλ] = constraint.balanceMoment -= Δλ[length(rotationConstraints)+i:length(rotationConstraints)+i+nλ-1]*ΔλRelaxFactor/σ
        # Update BC
        update_BC_data!(problem.model.BCs[hingeAxisConstraintBalanceLoadBCid[i]],problem.timeNow)
        # Update rotation parameters vector and angle of rotation across the hinge
        constraint.Δp = x[constraint.slaveElementGlobalDOFs] - x[constraint.masterElementGlobalDOFs]
        constraint.Δϕ = rotation_angle_limited(constraint.Δp)
    end

    @pack! problem = x,Δx,residual

    return true

end


# Performs a line search to solve the current Newton step
function line_search(x,residual,jacobian,λ=1e-8)

    # Regularize the inverse of the Jacobian
    jacobianInverse = inv(jacobian + λ * I(size(jacobian,1)))

    # Solve for the Newton step
    p = -jacobianInverse*residual

    # Find step size
    α = line_search_step_size(x,p,residual,jacobian)

    # Update the solution increment
    Δx = α * p

    return Δx

end


# Updates the step size (α) of the line search
function line_search_step_size(x,p,residual,jacobian,c1=1e-4,ρ=0.5)

    # Copy problem and update TF to skip Jacobian update
    problemCopy = deepcopy(problem)
    problemCopy.skipJacobianUpdate = true

    # Inline function to get the residual
    res = x -> assemble_system_arrays!(problemCopy,x) 

    # Initialize step size and fixed values
    α = 1.0
    β = dot(p, jacobian*p)
    γ = norm(residual)^2

    # Update step size
    while norm(res(x + α * p))^2 > γ + c1 * α * β
        α *= ρ
    end

    return α
end


# Solves a linear, overdetermined system with constraints
function linear_solver_with_constraints(x,jacobian,residual,rotationConstraints,hingeAxisConstraints,problem)

    # Initialize arrays of overdetermined system
    A = copy(jacobian)
    b = copy(-residual)

    # Size of nominal system (without constraints)
    N = length(b)

    # Loop rotation constraints
    for constraint in rotationConstraints
        @unpack masterGlobalDOF,slaveGlobalDOF = constraint
        # Initialize current value of Lagrange multiplier coefficients array
        c = zeros(size(A,1))
        # Adjust arrays by adding the equation: Δx[slaveGlobalDOF] - Δx[masterGlobalDOF] = 0
        c[slaveGlobalDOF] = 1
        c[masterGlobalDOF] = -1
        A = [A c; c' 0]
        b = [b; 0]
    end

    # Loop hinge axis constraints
    for constraint in hingeAxisConstraints
        @unpack masterElementGlobalDOFs,slaveElementGlobalDOFs,masterElementGlobalMasterDOF,masterElementGlobalSlaveDOFs,slaveElementGlobalMasterDOF,slaveElementGlobalSlaveDOFs,initialHingeAxis,updateHingeAxis,masterDOF,slaveDOFs,masterElementGlobalID,ΔpValue = constraint
        # Get current rotation parameters and rotation tensor (from basis b to basis B, resolved in basis A) of master element of the hinge constraint
        @unpack p = problem.model.elements[masterElementGlobalID].states
        @unpack R = problem.model.elements[masterElementGlobalID]
        # Update current hinge axis vector, resolved in basis A
        currentHingeAxis = R*initialHingeAxis
        # Compute derivatives of normalized hinge axis components w.r.t rotation parameters (of master element)
        ∂hingeAxis∂pMasterElem = updateHingeAxis ? ForwardDiff.jacobian(x -> hinge_axis_components(x,initialHingeAxis,masterDOF), p) : zeros(3,3)
        # Hinge axis for computations
        hingeAxis = updateHingeAxis ? currentHingeAxis : initialHingeAxis
        # Loop slave DOFs and set hinge axis constraints
        for j=1:2
            # Initialize current value of Lagrange multiplier coefficients array
            c = zeros(size(A,1))            
            # Adjust arrays by adding the equation: Δx[slaveElementGlobalSlaveDOFs[j]] - Δx[masterElementGlobalSlaveDOFs[j]] - ( Δx[slaveElementGlobalMasterDOF] - Δx[masterElementGlobalMasterDOF] ) * hingeAxis[slaveDOFs[j]]/hingeAxis[masterDOF] - ( x[slaveElementGlobalMasterDOF] - x[masterElementGlobalMasterDOF] ) * ∂hingeAxis∂pMasterElem[j,:] * Δx[masterElementGlobalDOFs] = 0
            c[slaveElementGlobalSlaveDOFs[j]] = 1
            c[masterElementGlobalSlaveDOFs[j]] = -1
            c[slaveElementGlobalMasterDOF] = -hingeAxis[slaveDOFs[j]]/hingeAxis[masterDOF]
            c[masterElementGlobalMasterDOF] = hingeAxis[slaveDOFs[j]]/hingeAxis[masterDOF]
            c[masterElementGlobalDOFs] += -(x[slaveElementGlobalMasterDOF] - x[masterElementGlobalMasterDOF]) * ∂hingeAxis∂pMasterElem[slaveDOFs[j],:]
            A = [A c; c' 0]
            b = [b; 0]
        end
        # Set rotation norm constraint, if applicable
        if !isnothing(ΔpValue) && !iszero(ΔpValue)
            # Initialize current value of Lagrange multiplier coefficients array
            c = zeros(size(A,1))
            # Adjust arrays by adding the equation: (x[slaveElementGlobalDOFs] - x[masterElementGlobalDOFs])' * (Δx[slaveElementGlobalDOFs] - Δx[masterElementGlobalDOFs]) = 0
            # Note: sign(ΔpValue) is used to yield the correct sign on Δλ
            c[slaveElementGlobalDOFs] = sign(ΔpValue) * (x[slaveElementGlobalDOFs] - x[masterElementGlobalDOFs])
            c[masterElementGlobalDOFs] = -sign(ΔpValue) * (x[slaveElementGlobalDOFs] - x[masterElementGlobalDOFs])
            A = [A c; c' 0]
            b = [b; 0]
        elseif !isnothing(ΔpValue) && iszero(ΔpValue)
            # Initialize current value of Lagrange multiplier coefficients array
            c = zeros(size(A,1))
            # Adjust arrays by adding the equation: Δx[slaveElementGlobalDOFs] - Δx[masterElementGlobalDOFs] = 0
            c[slaveElementGlobalDOFs] .= 1
            c[masterElementGlobalDOFs] .= -1
            A = [A c; c' 0]
            b = [b; 0]
        end
        # Pack data
        @pack! constraint = currentHingeAxis
    end

    # Solve constrained linear system for solution and Lagrange multipliers increments
    sol = A\b
    Δx = sol[1:N]
    Δλ = sol[N+1:end]

    return Δx,Δλ
end

# Computes the current hinge axis normalized by its reference direction (masterDOF) value
function hinge_axis_components(p,initialHingeAxis,masterDOF)

    R,_ = rotation_tensor_WM(p)

    currentHingeAxis = R*initialHingeAxis
        
    f = currentHingeAxis/currentHingeAxis[masterDOF]

    return f
end

# Saves the solution at the current load factor
function save_load_factor_data!(problem::Problem,σ::Float64,x::Vector{Float64})


    @unpack savedσ,xOverσ,elementalStatesOverσ,nodalStatesOverσ,compElementalStatesOverσ,aeroVariablesOverσ = problem 
    @unpack elements = problem.model   

    # Add current load factor
    push!(savedσ,σ)

    # Add curent system states 
    push!(xOverσ,x)

    # Add current elemental states 
    currentElementalStates = Vector{ElementalStates}()
    for element in elements
        push!(currentElementalStates,deepcopy(element.states))
    end
    push!(elementalStatesOverσ,currentElementalStates)

    # Add current nodal states 
    currentNodalStates = Vector{NodalStates}()
    for element in elements
        push!(currentNodalStates,deepcopy(element.nodalStates))
    end
    push!(nodalStatesOverσ,currentNodalStates)

    # Add current complementary elemental states 
    currentComplementaryElementalStates = Vector{ComplementaryElementalStates}()
    for element in elements
        push!(currentComplementaryElementalStates,deepcopy(element.compStates))
    end
    push!(compElementalStatesOverσ,currentComplementaryElementalStates)

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
    push!(aeroVariablesOverσ,currentAeroVariables)

    @pack! problem = savedσ,xOverσ,elementalStatesOverσ,nodalStatesOverσ,compElementalStatesOverσ,aeroVariablesOverσ

end


# Updates the maximum absolute values of aerodynamic loads
function update_maximum_aero_loads!(problem::Problem)

    @unpack maxAeroForce,maxAeroMoment = problem

    # Initialize arrays
    aeroForces = Vector{Float64}()
    aeroMoment = Vector{Float64}()

    # Loop elements
    for element in problem.model.elements
        # Skip elements without aero
        if isnothing(element.aero)
            continue
        end
        push!(aeroForces,element.aero.aeroCoefficients.ct,element.aero.aeroCoefficients.cn)
        push!(aeroMoment,element.aero.aeroCoefficients.cm)
    end

    # Update, if applicable
    if !isempty(aeroForces)
        maxAeroForce = max(maxAeroForce,maximum(abs.(aeroForces)))
        maxAeroMoment = max(maxAeroMoment,maximum(abs.(aeroMoment)))
    end

    @pack! problem = maxAeroForce,maxAeroMoment
end