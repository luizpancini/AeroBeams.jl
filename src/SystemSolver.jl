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
    allowAdvanceThroughUnconvergedAeroStates::Bool
    aeroStatesResidualRatioThreshold::Float64
    aeroStatesRelativeErrorThreshold::Float64

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
- `absoluteTolerance::Float64`: absolute convergence tolerance
- `relativeTolerance::Float64`: relative convergence tolerance
- `maximumIterations::Int64`: maximum number of iterations
- `desiredIterations::Int64`: desired number of iterations
- `maximumAbsoluteError::Real`: maximum absolute error for divergence detection
- `maximumRelativeError::Real`: maximum relative error for divergence detection
- `initialLoadFactor::Real`: initial load factor
- `minimumLoadFactor::Float64`: minimum load factor
- `maximumLoadFactorStep::Float64`: maximum load factor step
- `minimumLoadFactorStep::Float64`: minimum load factor step
- `ρ::Real`: relaxation factor for trim variables
- `trackingLoadSteps::Bool`: flag to track partial load steps solutions
- `displayStatus::Bool`: flag to display status
- `minConvRateAeroJacUpdate::Real`: minimum convergence rate to skip computation of aerodynamic Jacobians
- `minConvRateJacUpdate::Real`: minimum convergence rate to skip computation of structural Jacobians
- `alwaysUpdateJacobian::Bool`: flag to update Jacobians on every iteration
- `allowAdvanceThroughUnconvergedAeroStates::Bool`: flag to allow the advancement of the algorithm through unconverged aerodynamic states
- `aeroStatesResidualRatioThreshold::Float64`: threshold ratio of aerodynamic/total residuals in order to allow the above flag
- `aeroStatesRelativeErrorThreshold::Float64`: threshold relative error in order to allow the above flag
"""
function create_NewtonRaphson(; absoluteTolerance::Float64=1e-8,relativeTolerance::Float64=1e-8,maximumIterations::Int64=20,desiredIterations::Int64=5,maximumAbsoluteError::Real=1e6,maximumRelativeError::Real=1e6,initialLoadFactor::Real=1.0,minimumLoadFactor::Float64=0.01,maximumLoadFactorStep::Float64=0.5,minimumLoadFactorStep::Float64=0.01,ρ::Real=1.0,trackingLoadSteps::Bool=true,displayStatus::Bool=false,minConvRateAeroJacUpdate::Real=2.0,minConvRateJacUpdate::Real=2.0,alwaysUpdateJacobian::Bool=true,allowAdvanceThroughUnconvergedAeroStates::Bool=false,aeroStatesResidualRatioThreshold::Float64=0.95,aeroStatesRelativeErrorThreshold::Float64=1e-2)

    @assert maximumIterations > 1 "maximumIterations must be greater than 1"
    @assert desiredIterations > 1 "desiredIterations must be greater than 1"
    @assert 0 <= initialLoadFactor <= 1 "initialLoadFactor must be between 0 and 1"
    @assert 0 < minimumLoadFactor <= 0.5 "minimumLoadFactor must be between 0 and 0.5"
    @assert 0 < maximumLoadFactorStep <= 1 "maximumLoadFactorStep must be between 0 and 1"
    @assert 0 < minimumLoadFactorStep <= 1 "minimumLoadFactorStep must be between 0 and 1"
    @assert 0.5 <= ρ <= 1 "relaxation factor for trim variables (ρ) must be between 0.5 and 1.0"
    @assert minConvRateAeroJacUpdate > 1 "minConvRateAeroJacUpdate to skip calculation of aerodynamic derivatives must be greater than 1"
    @assert minConvRateJacUpdate > 1 "minConvRateJacUpdate to skip update of Jacobian matrix must be greater than 1"
    @assert 0 < aeroStatesResidualRatioThreshold < 1 "aeroStatesResidualRatioThreshold must be between 0 and 1"

    return NewtonRaphson(absoluteTolerance=absoluteTolerance,relativeTolerance=relativeTolerance,maximumIterations=maximumIterations,desiredIterations=desiredIterations,maximumAbsoluteError=maximumAbsoluteError,maximumRelativeError=maximumRelativeError,initialLoadFactor=initialLoadFactor,minimumLoadFactor=minimumLoadFactor,maximumLoadFactorStep=maximumLoadFactorStep,minimumLoadFactorStep=minimumLoadFactorStep,ρ=ρ,trackingLoadSteps=trackingLoadSteps,displayStatus=displayStatus,minConvRateAeroJacUpdate=minConvRateAeroJacUpdate,minConvRateJacUpdate=minConvRateJacUpdate,alwaysUpdateJacobian=alwaysUpdateJacobian,allowAdvanceThroughUnconvergedAeroStates=allowAdvanceThroughUnconvergedAeroStates,aeroStatesResidualRatioThreshold=aeroStatesResidualRatioThreshold,aeroStatesRelativeErrorThreshold=aeroStatesRelativeErrorThreshold)

end
export create_NewtonRaphson


# Solves the nonlinear system of equations at current time step
function solve_NewtonRaphson!(problem::Problem)

    # Unpack system solver
    @unpack absoluteTolerance,relativeTolerance,maximumIterations,desiredIterations,maximumAbsoluteError,maximumRelativeError,minimumLoadFactor,loadFactorStep,maximumLoadFactorStep,minimumLoadFactorStep,trackingLoadSteps,displayStatus,convergedFinalSolution,convergedPartialSolution,alwaysUpdateJacobian,minConvRateJacUpdate,allowAdvanceThroughUnconvergedAeroStates,aeroStatesResidualRatioThreshold,aeroStatesRelativeErrorThreshold = problem.systemSolver

    # Reset converged solution flag and load factor step
    convergedPartialSolution = convergedFinalSolution = false
    loadFactorStep = max(min(maximumLoadFactorStep,problem.σ),minimumLoadFactorStep)
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
        modelKnown = deepcopy(problem.model)
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
                # Check special case
                if allowAdvanceThroughUnconvergedAeroStates && norm(residual[problem.model.DOF_χ_all])/ϵ_abs > aeroStatesResidualRatioThreshold && ϵ_rel < aeroStatesRelativeErrorThreshold
                    convergedPartialSolution = true
                    println("Advancing through unconverged aero states...")
                end
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
            model = modelKnown
            if displayStatus
                println("Unconverged, reducing load factor and trying again...")
            end
            @pack! problem = x,residual,σ,model
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
            # Assemble inertia matrix in eigen and trim problems, if applicable
            assemble_inertia_matrix!(problem)
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

    # Reset spring loads' Jacobians
    reset_spring_loads_jacobians!(problem,specialNodes)

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


# Assembles the inertia matrix
function assemble_inertia_matrix!(problem::Problem)

    # Assemble only for eigenproblems and trim problems, when required
    if problem isa EigenProblem || (problem isa TrimProblem && problem.getInertiaMatrix)
        # Get contributions from each element
        for element in problem.model.elements
            element_inertia!(problem,problem.model,element)
        end
        # Set contributions from hinge constraint equations
        if problem.model.hasHingeAxisConstraints
            augmented_inertia_matrix!(problem)
        end
    end

    return nothing
end


# Solves the linear system of equations at current time step and load factor
function solve_linear_system!(problem::Problem)

    @unpack x,Δx,residual,jacobian,getLinearSolution,σ,timeNow = problem
    @unpack systemOrder,nTrimVariables,hasHingeAxisConstraints,hingeAxisConstraints,forceScaling = problem.model
    @unpack ρ = problem.systemSolver

    # Issue warning for singular Jacobian in linear problems
    if getLinearSolution && nTrimVariables == 0 && det(jacobian) ≈ 0
        println("Singular Jacobian")
    end

    # Solve the linear system according to problem type
    #---------------------------------------------------------------------------
    # Trim problem
    if problem isa TrimProblem
        Δx .= -pinv(Matrix(jacobian))*residual
    # Problem with hinge axis constraints    
    elseif hasHingeAxisConstraints
        Δx,Δλ = linear_solver_with_constraints(problem,x,jacobian,residual,hingeAxisConstraints)
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

    # Update hinge constraint data
    for (i,constraint) in enumerate(hingeAxisConstraints)
        @unpack solutionMethod,rotationIsFixed,initialHingeAxis,masterElementGlobalDOFs,slaveElementGlobalDOFs,slaveDOFs,balanceMomentBCID,balanceMoment,λ,Jc = constraint
        # Update Lagrange multipliers
        λ += Δλ[i:i+length(λ)-1]
        # Update balance moment
        balanceMoment = -Jc[:,slaveElementGlobalDOFs[slaveDOFs]]'*λ/σ
        # Update BCs, if applicable
        if solutionMethod == "appliedMoment"
            problem.model.BCs[balanceMomentBCID].values = balanceMoment
            update_BC_data!(problem.model.BCs[balanceMomentBCID],timeNow)
        end
        # Current rotation parameters of master and slave elements
        pM = x[masterElementGlobalDOFs]
        pS = x[slaveElementGlobalDOFs]
        # Update current hinge axis, hinge rotation parameter vector, signed magnitude (value) and hinge angle
        currentHingeAxis = current_hinge_axis(pM,initialHingeAxis)
        pHValue = hinge_rotation_value(pM,pS,initialHingeAxis)
        pH = hinge_rotation_parameters(pM,pS)
        ϕ = rotation_angle_limited(pH)
        # Pack data
        @pack! constraint = currentHingeAxis,pH,ϕ,λ,balanceMoment
        if !rotationIsFixed
            @pack! constraint = pHValue
        end
    end
    
    @pack! problem = x,Δx

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
function linear_solver_with_constraints(problem,x,jacobian,residual,hingeAxisConstraints)

    # Size of nominal system (without constraints)
    N = length(residual)

    A = Matrix(jacobian)
    if any(isnan, A) || any(isinf, A)
        println("jacobian contains NaNs")
    end

    # Initialize augmented arrays
    augmentedJacobian = copy(jacobian)
    augmentedResidual = copy(residual)

    # Loop hinge axis constraints
    for constraint in hingeAxisConstraints
        @unpack solutionMethod,rotationIsFixed,initialHingeAxis,pHValue,masterElementGlobalDOFs,slaveElementGlobalDOFs,slaveDOFs,λ = constraint
        # Current rotation parameters of master and slave elements
        pM = x[masterElementGlobalDOFs]
        pS = x[slaveElementGlobalDOFs]
        # If hinge rotation is fixed (known)
        if rotationIsFixed
            # Residual and Jacobian terms
            resC = C(pM,pS,initialHingeAxis,pHValue=pHValue)
            Jc = zeros(3,length(augmentedResidual))
            Jb = spzeros(length(augmentedResidual),length(augmentedResidual))
            Jc[:,masterElementGlobalDOFs] .= ∂C_∂pM(pM,pS,initialHingeAxis,pHValue=pHValue)
            Jc[:,slaveElementGlobalDOFs] .= ∂C_∂pS(pM,pS,initialHingeAxis,pHValue=pHValue)
            Jb[masterElementGlobalDOFs,masterElementGlobalDOFs] .= ∂2CTλ_∂pM2(pM,pS,initialHingeAxis,λ,pHValue=pHValue)
            Jb[slaveElementGlobalDOFs,slaveElementGlobalDOFs] .= ∂2CTλ_∂pS2(pM,pS,initialHingeAxis,λ,pHValue=pHValue)
            Jb[masterElementGlobalDOFs,slaveElementGlobalDOFs] .= ∂2CTλ_∂pMpS(pM,pS,initialHingeAxis,λ,pHValue=pHValue)
            Jb[slaveElementGlobalDOFs,masterElementGlobalDOFs] .= ∂2CTλ_∂pSpM(pM,pS,initialHingeAxis,λ,pHValue=pHValue)
        # If the hinge rotation is not fixed (is unknown)
        else
            # Scaled rotation parameters of slave element
            pSscaled = scaled_rotation_parameters(pS)
            # Residual and Jacobian terms of the constraint
            resC = C(pM,pSscaled,initialHingeAxis,slaveDOFs=slaveDOFs)
            Jc = zeros(2,length(augmentedResidual))
            Jb = spzeros(length(augmentedResidual),length(augmentedResidual))
            Jc[:,masterElementGlobalDOFs] .= ∂C_∂pM(pM,pSscaled,initialHingeAxis,slaveDOFs=slaveDOFs)
            Jc[:,slaveElementGlobalDOFs] .= ∂C_∂pS(pM,pSscaled,initialHingeAxis,slaveDOFs=slaveDOFs)
            Jb[masterElementGlobalDOFs,masterElementGlobalDOFs] .= ∂2CTλ_∂pM2(pM,pSscaled,initialHingeAxis,λ,slaveDOFs=slaveDOFs)
            Jb[slaveElementGlobalDOFs,slaveElementGlobalDOFs] .= ∂2CTλ_∂pS2(pM,pSscaled,initialHingeAxis,λ,slaveDOFs=slaveDOFs)
            Jb[masterElementGlobalDOFs,slaveElementGlobalDOFs] .= ∂2CTλ_∂pMpS(pM,pS,initialHingeAxis,λ,slaveDOFs=slaveDOFs)
            Jb[slaveElementGlobalDOFs,masterElementGlobalDOFs] .= ∂2CTλ_∂pSpM(pM,pS,initialHingeAxis,λ,slaveDOFs=slaveDOFs)
        end
        # Update residual of nominal system, if applicable
        if solutionMethod == "addedResidual"
            residual .+= Jc'*λ
        end
        # Adjust augmented residual and Jacobian arrays by adding contributions from the constraint equations
        augmentedJacobian = [jacobian.+Jb Jc'; Jc zeros(size(Jc,1),size(Jc,1))]
        augmentedResidual = [residual; resC]
        # Pack data
        @pack! constraint = Jc
    end

    A = Matrix(augmentedJacobian)
    if any(isnan, A) || any(isinf, A)
        println("augmentedJacobian contains NaNs")
    end
    if any(isnan, augmentedResidual) || any(isinf, augmentedResidual)
        println("augmentedResidual contains NaNs")
    end

    # Solve constrained linear system for solution and Lagrange multipliers increments
    sol = -augmentedJacobian\augmentedResidual
    Δx = sol[1:N]
    Δλ = sol[N+1:end]

    # Pack data (updated residual and augmented Jacobian)
    @pack! problem = residual
    if problem isa EigenProblem
        @pack! problem = augmentedJacobian
    end

    return Δx,Δλ
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