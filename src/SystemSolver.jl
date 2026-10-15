## SYSTEM SOLVER TYPES
# ------------------------------------------------------------------------------
# Load stepping properties composite type
# Fields
# - `fullLoadFacTol::Float64`: tolerance for the full load factor
# - `initLoadFac::Float64`: initial load factor
# - `minLoadFac::Float64`: minimum load factor
# - `maxLoadFacStep::Float64`: maximum load factor step
# - `minLoadFacStep::Float64`: minimum load factor step
# - `trackLoadSteps::Bool`: flag to track partial load steps solutions
#
struct LoadSteppingProperties

    fullLoadFacTol::Float64

    initLoadFac::Float64
    minLoadFac::Float64
    maxLoadFacStep::Float64
    minLoadFacStep::Float64

    trackLoadSteps::Bool

end


# Load stepping state composite type
# Fields
# - `loadstep::Int`: current load step counter
# - `σ::Float64`: current load factor
# - `σStep::Float64`: current load factor step
# - `convergedLoadStep::Bool`: flag indicating if the current load step has converged
# - `convergedFinalSolution::Bool`: flag indicating if the final solution has converged
# - `xKnown::Vector{Float64}`: vector of known solution components
#
mutable struct LoadSteppingState

    loadstep::Int

    σ::Float64
    σStep::Float64

    convergedLoadStep::Bool
    convergedFinalSolution::Bool

    xKnown::Vector{Float64}

end
LoadSteppingState(minLoadFacStep,maxLoadFacStep,initLoadFac) = 
    LoadSteppingState(
        0,
        initLoadFac, clamp(initLoadFac,minLoadFacStep,maxLoadFacStep), 
        false, false,
        Float64[]
    )
LoadSteppingState(loadStepProps::LoadSteppingProperties) = 
    LoadSteppingState(loadStepProps.minLoadFacStep, loadStepProps.maxLoadFacStep, loadStepProps.initLoadFac)
LoadSteppingState() = LoadSteppingState(0,1.0,1.0)    


# TrimSolver composite type
# Fields
# - `method::Symbol`: method for linear system solution
# - `damping::Union{Nothing,Float64}`: damping term, for damped least squares method
@with_kw struct TrimSolver

    method::Symbol
    damping::Union{Nothing,Float64}

end
export TrimSolver


# LinearSolverState composite type
# Fields
# - `converged::Bool`: flag indicating if the linear system solution has converged
#
mutable struct LinearSolverState

    converged::Bool

end
LinearSolverState() = LinearSolverState(false)


# LinearSolver composite type
# Fields
# - `loadStepState::LoadSteppingState`: state of the load stepping procedure
# - `state::LinearSolverState`: state of the linear solver
#
@with_kw mutable struct LinearSolver

    loadStepState::LoadSteppingState = LoadSteppingState()

    state::LinearSolverState = LinearSolverState()

end
export LinearSolver


# NRProperties composite type
# Fields
# - `method::Symbol`: method for nonlinear solution (:Newton for full Newton iterations or :NewtonChord to use chord method with Jacobian update skipping)
# - `absTol::Float64`: absolute convergence tolerance
# - `relTol::Float64`: relative convergence tolerance
# - `maxIter::Int`: maximum number of iterations
# - `desIter::Int`: desired number of iterations
# - `maxAbsErr::Float64`: maximum absolute error for divergence detection
# - `maxRelErr::Float64`: maximum relative error for divergence detection
# - `minConvRateJacUpdate::Float64`: minimum convergence rate to skip computation of structural Jacobians
# - `minConvRateAeroJacUpdate::Float64`: minimum convergence rate to skip computation of aerodynamic Jacobians
# - `allowAdvanceThroughUnconvergedAeroStates::Bool`: flag to allow the advancement of the algorithm through unconverged aerodynamic states
# - `aeroResRatioThreshold::Float64`: threshold ratio of aerodynamic/total residuals in order to allow the above flag
# - `aeroΔxErrorThreshold::Float64`: threshold relative error in order to allow the above flag
#
struct NRProperties

    method::Symbol

    absTol::Float64
    relTol::Float64

    maxIter::Int
    desIter::Int

    maxAbsErr::Float64
    maxRelErr::Float64
    
    displayStatus::Bool

    minConvRateJacUpdate::Float64
    minConvRateAeroJacUpdate::Float64

    allowAdvanceThroughUnconvergedAeroStates::Bool
    aeroResRatioThreshold::Float64
    aeroΔxErrorThreshold::Float64
    
end


# NRState composite type
# Fields
# - `iter::Int`: current iteration counter
# - `ϵ_abs`: current absolute (residual) convergence norm
# - `ϵ_rel`: current relative (solution) convergence norm
# - `ϵ_abs_previous`: previous absolute (residual) convergence norm
# - `convRate::Float64`: current convergence rate of residual norm
# - `skipJacobianUpdate`: flag to compute only residual array and skip Jacobian update at the current iteration
# - `skipAeroJacobianUpdate::Bool`: flag to skip aerodynamic Jacobian updates at the current iteration
# - `converged::Bool`: flag indicating if the current linear solution has converged
#
mutable struct NRState

    iter::Int

    ϵ_abs::Float64
    ϵ_rel::Float64

    ϵ_abs_previous::Float64
    convRate::Float64

    skipJacobianUpdate::Bool
    skipAeroJacobianUpdate::Bool

    converged::Bool

end
NRState() = NRState(0, 0.0, 0.0, 0.0, 0.0, false, false, false)


# Newton-Raphson solver composite type
# Fields
# - `props::NRProperties`: properties of the Newton-Raphson solver
# - `linearTrimSolver::TrimSolver`: linear trim solver
# - `loadStepProps::LoadSteppingProperties`: properties of the load stepping procedure
# - `loadStepState::LoadSteppingState`: state of the load stepping procedure
# - `state::NRState`: state of the Newton-Raphson solver
#
@with_kw mutable struct NewtonRaphson

    props::NRProperties

    linearTrimSolver::TrimSolver

    loadStepProps::LoadSteppingProperties
    loadStepState::LoadSteppingState  = LoadSteppingState(loadStepProps)

    state::NRState = NRState()

end


# NLsolver properties composite type
# Fields
# - `method::Symbol`: nonlinear solve method
# - `linesearch::Bool`: flag to enable linesearch method for the :newton solver
# - `jcalc::Symbol`: Jacobian calculation method
# - `checkJac::Bool`: flag to compare analytical and numerical Jacobians at each iteration
# - `ftol::Float64`: absolute convergence tolerance
# - `maxIter::Int`: maximum number of iterations
# - `desIter::Int`: desired number of iterations
# - `displayStatus::Bool`: flag to display status
#
struct NLsolverProperties

    method::Symbol
    linesearch::Bool
    jcalc::Symbol
    checkJac::Bool

    ftol::Float64

    maxIter::Int
    desIter::Int

    displayStatus::Bool

end


# NLsolver state composite type
# Fields
# - `iter::Int`: current iteration counter
# - `ϵ_abs::Float64`: current absolute (residual) convergence norm
# - `ϵ_rel::Float64`: current relative (solution) convergence norm
# - `converged::Bool`: flag indicating if the nonlinear solver has converged
#
mutable struct NLsolverState

    iter::Int

    ϵ_abs::Float64
    ϵ_rel::Float64

    converged::Bool

end
NLsolverState() = NLsolverState(0, 0.0, 0.0, false)


# NLsolver composite type
# Fields
# - `props::NLsolverProperties`: properties of the nonlinear solver
# - `loadStepProps::LoadSteppingProperties`: properties of the load stepping procedure
# - `loadStepState::LoadSteppingState`: state of the load stepping procedure
# - `state::NLsolverState`: state of the nonlinear solver
#
@with_kw mutable struct NLsolver

    props::NLsolverProperties

    linearTrimSolver::TrimSolver

    loadStepProps::LoadSteppingProperties
    loadStepState::LoadSteppingState = LoadSteppingState(loadStepProps)

    state::NLsolverState = NLsolverState()

end
export NLsolver


# ArclengthSolver properties composite type
# Fields
# - `absTol::Float64`: absolute convergence tolerance for the residual
# - `maxIter::Int`: maximum number of corrector iterations per arclength step
# - `desIter::Int`: desired number of corrector iterations per arclength step
# - `initArcLength::Float64`: initial arclength step
# - `minArcLength::Float64`: minimum arclength step
# - `maxArcLength::Float64`: maximum arclength step
# - `maxLoadFac::Float64`: load factor at which the continuation stops (first converged point with σ ≥ maxLoadFac)
# - `maxSteps::Int`: maximum number of arclength steps
# - `displayStatus::Bool`: flag to display status
#
struct ArclengthSolverProperties

    absTol::Float64

    maxIter::Int
    desIter::Int

    initArcLength::Float64
    minArcLength::Float64
    maxArcLength::Float64

    maxLoadFac::Float64
    maxSteps::Int

    displayStatus::Bool

end


# ArclengthSolver state composite type
# Fields
# - `iter::Int`: corrector iteration counter of the current arclength step
# - `ϵ_abs::Float64`: current absolute (residual) convergence norm
# - `arcLength::Float64`: current arclength step
# - `converged::Bool`: flag indicating if the current arclength step has converged
#
mutable struct ArclengthSolverState

    iter::Int

    ϵ_abs::Float64
    arcLength::Float64

    converged::Bool

end
ArclengthSolverState(arcLength=0.0) = ArclengthSolverState(0, 0.0, arcLength, false)


# ArclengthSolver composite type (spherical arclength continuation, for problems with limit points such as snap-through)
# Fields
# - `props::ArclengthSolverProperties`: properties of the arclength solver
# - `loadStepProps::LoadSteppingProperties`: properties of the load stepping procedure (only initLoadFac, minLoadFac and trackLoadSteps are used)
# - `loadStepState::LoadSteppingState`: state of the load stepping procedure (σ is the load factor, which is an unknown of the problem)
# - `state::ArclengthSolverState`: state of the arclength solver
#
@with_kw mutable struct ArclengthSolver

    props::ArclengthSolverProperties

    loadStepProps::LoadSteppingProperties
    loadStepState::LoadSteppingState = LoadSteppingState(loadStepProps)

    state::ArclengthSolverState = ArclengthSolverState(props.initArcLength)

end
export ArclengthSolver

const NonlinearSolver = Union{NewtonRaphson, NLsolver, ArclengthSolver}
const SystemSolver = Union{LinearSolver, NewtonRaphson, NLsolver, ArclengthSolver}

## PUBLIC API
# ----------------------------------------------------------------------------
"""
    create_TrimSolver(; kwargs...)

Trim solver constructor

# Keyword arguments
- `method::Symbol`: method for trim (least-squares) system solution: either :LQ (LQ decomposition), :LM (Levenberg-Marquardt with fixed relative damping), or :SVD (Singular Value Decomposition) or :pinv (Moore-Penrose pseudo-inverse, which yields minimum norm and works with rank-deficient Jacobians)
- `damping::Union{Nothing,Real}`: damping term for method = :LM
"""
function create_TrimSolver(; 
    method::Symbol=:LM, 
    damping::Union{Nothing,Float64}=1e-14
    )
            
    # Validate inputs
    @assert method in [:LQ, :LM, :SVD, :pinv] "set method as :LQ, :LM, :SVD, or :pinv"
    if method == :LM
        @assert !isnothing(damping) "damping must be provided for :LM method"
        @assert damping > 0 "damping must be positive for :LM method"
    end

    # Set correct types
    damping = isnothing(damping) ? 
              nothing :
              Float64(damping)

    return TrimSolver(method,damping)

end
export create_TrimSolver


"""
    create_LinearSolver()

Linear system solver constructor

"""
function create_LinearSolver()

    return LinearSolver()

end
export create_LinearSolver


"""
    create_NewtonRaphson(; kwargs...)

Newton-Raphson nonlinear system solver constructor

# Keyword arguments
- `method::Symbol`: method for nonlinear solution
- `absTol::Float64`: absolute convergence tolerance
- `relTol::Float64`: relative convergence tolerance
- `maxIter::Int`: maximum number of iterations
- `desIter::Int`: desired number of iterations
- `maxAbsErr::Real`: maximum absolute error for divergence detection
- `maxRelErr::Real`: maximum relative error for divergence detection
- `fullLoadFacTol::Float64`: load step tolerance for the final solution convergence
- `initLoadFac::Real`: initial load factor
- `minLoadFac::Float64`: minimum load factor
- `maxLoadFacStep::Float64`: maximum load factor step
- `minLoadFacStep::Float64`: minimum load factor step
- `trackLoadSteps::Bool`: flag to track partial load steps solutions
- `displayStatus::Bool`: flag to display status
- `minConvRateJacUpdate::Real`: minimum convergence rate to skip computation of structural Jacobians when using :NewtonChord method
- `minConvRateAeroJacUpdate::Real`: minimum convergence rate to skip computation of aerodynamic Jacobians when using :NewtonChord method
- `allowAdvanceThroughUnconvergedAeroStates::Bool`: flag to allow the advancement of the algorithm through unconverged aerodynamic states
- `aeroResRatioThreshold::Float64`: threshold ratio of aerodynamic/total residuals in order to allow the above flag
- `aeroΔxErrorThreshold::Float64`: threshold ratio of aerodynamic/total Δx in order to allow the above flag
- `trimSolverMethod::Symbol`: method for linear system solution in trim problems
- `dampingTrim::Float64`: damping term for method = :LM
"""
function create_NewtonRaphson(;
    method::Symbol=:Newton, 
    absTol::Float64=1e-8,
    relTol::Float64=1e-8,
    maxIter::Int=50,
    desIter::Int=5,
    maxAbsErr::Real=1e6,
    maxRelErr::Real=1e300,
    fullLoadFacTol::Float64=1e-6,
    initLoadFac::Real=1.0,
    minLoadFac::Float64=0.01,
    maxLoadFacStep::Float64=1.0,
    minLoadFacStep::Float64=0.01,
    trackLoadSteps::Bool=true,
    displayStatus::Bool=false,
    minConvRateJacUpdate::Real=2.0,
    minConvRateAeroJacUpdate::Real=2.0,
    allowAdvanceThroughUnconvergedAeroStates::Bool=true,
    aeroResRatioThreshold::Float64=0.0,
    aeroΔxErrorThreshold::Float64=0.0,
    trimSolverMethod::Symbol=:LM,
    dampingTrim::Float64=1e-14
    )

    # Validate inputs
    @assert method in [:Newton, :NewtonChord] "set method as one of :Newton or :NewtonChord"
    @assert maxIter > 1 "maxIter must be greater than 1"
    @assert desIter > 1 "desIter must be greater than 1"
    @assert 0 <= initLoadFac <= 1 "initLoadFac must be between 0 and 1"
    @assert 0 < minLoadFac <= 0.5 "minLoadFac must be between 0 and 0.5"
    @assert 0 < maxLoadFacStep <= 1 "maxLoadFacStep must be between 0 and 1"
    @assert 0 < minLoadFacStep <= 1 "minLoadFacStep must be between 0 and 1"
    @assert minConvRateJacUpdate > 1 "minConvRateJacUpdate to skip update of Jacobian matrix must be greater than 1"
    @assert minConvRateAeroJacUpdate > 1 "minConvRateAeroJacUpdate to skip calculation of aerodynamic derivatives must be greater than 1"
    @assert 0 <= aeroResRatioThreshold < 1 "aeroResRatioThreshold must be between 0 and 1"
    @assert 0 <= aeroΔxErrorThreshold < 1 "aeroΔxErrorThreshold must be between 0 and 1"
    @assert trimSolverMethod in [:LQ, :LM, :SVD, :pinv] "set trimSolverMethod as :LQ, :LM, :SVD, or :pinv"
    @assert dampingTrim > 0 "dampingTrim must be positive"

    # Set correct types for real-valued inputs
    maxAbsErr = Float64(maxAbsErr)
    maxRelErr = Float64(maxRelErr)
    initLoadFac = Float64(initLoadFac)
    minConvRateJacUpdate = Float64(minConvRateJacUpdate)
    minConvRateAeroJacUpdate = Float64(minConvRateAeroJacUpdate) 

    # Create Newton-Raphson properties
    props = NRProperties(
        method,
        absTol,
        relTol,
        maxIter,
        desIter,
        maxAbsErr,
        maxRelErr,
        displayStatus,
        minConvRateJacUpdate,
        minConvRateAeroJacUpdate,
        allowAdvanceThroughUnconvergedAeroStates,
        aeroResRatioThreshold,
        aeroΔxErrorThreshold
    )

    # Create linear trim solver
    linearTrimSolver = TrimSolver(
        trimSolverMethod,
        dampingTrim
        )

    # Create load stepping properties
    loadStepProps = LoadSteppingProperties(
        fullLoadFacTol,
        initLoadFac,
        minLoadFac,
        maxLoadFacStep,
        minLoadFacStep,
        trackLoadSteps,
    )

    return NewtonRaphson(
        props=props,
        linearTrimSolver=linearTrimSolver,
        loadStepProps=loadStepProps,
        )

end
export create_NewtonRaphson


"""
    create_NLsolver(; kwargs...)

NLSolve.jl nonlinear system solver constructor.

# Keyword arguments
- `method::Symbol`: nonlinear solve method
- `linesearch::Bool`: flag to enable linesearch method for the :newton method
- `jcalc::Symbol`: Jacobian calculation method, either :analytical or :numerical
- `checkJac::Bool`: flag to compare analytical and numerical Jacobians at each iteration
- `ftol::Float64`: absolute convergence tolerance for residual
- `maxIter::Int`: maximum number of iterations
- `desIter::Int`: desired number of iterations for the nonlinear solver
- `fullLoadFacTol::Float64`: load step tolerance for the final solution convergence
- `initLoadFac::Real`: initial load factor for the load stepping procedure
- `minLoadFac::Float64`: minimum load factor for the load stepping procedure
- `maxLoadFacStep::Float64`: maximum load factor step for the load stepping procedure
- `minLoadFacStep::Float64`: minimum load factor step for the load stepping procedure
- `trackLoadSteps::Bool`: flag to track load steps during the nonlinear solve
- `displayStatus::Bool`: flag to display convergence status
- `trimSolverMethod::Symbol`: method for linear system solution in trim problems
- `dampingTrim::Float64`: damping term for method = :LM
"""
function create_NLsolver(;
    method::Symbol=:newton,
    linesearch::Bool=false,
    jcalc::Symbol=:analytical,
    checkJac::Bool=false,
    ftol::Float64=1e-8,
    maxIter::Int=50,
    desIter::Int=5,
    fullLoadFacTol::Float64=1e-6,
    initLoadFac::Real=1.0,
    minLoadFac::Float64=0.01,
    maxLoadFacStep::Float64=1.0,
    minLoadFacStep::Float64=0.01,
    trackLoadSteps::Bool=true,
    displayStatus::Bool=false,
    trimSolverMethod::Symbol=:LM,
    dampingTrim::Float64=1e-14
    )

    # Validate inputs
    @assert method in [:newton,:trust_region] "select method as one of :newton, :trust_region"
    @assert jcalc in [:analytical, :numerical] "select jcalc as one of :analytical, :numerical"
    @assert ftol > 0 "ftol must be positive"
    @assert 0 <= initLoadFac <= 1 "initLoadFac must be between 0 and 1"
    @assert 0 < minLoadFac <= 0.5 "minLoadFac must be between 0 and 0.5"
    @assert 0 < maxLoadFacStep <= 1 "maxLoadFacStep must be between 0 and 1"
    @assert 0 < minLoadFacStep <= 1 "minLoadFacStep must be between 0 and 1"
    @assert trimSolverMethod in [:LQ, :LM, :SVD, :pinv] "set trimSolverMethod as :LQ, :LM, :SVD, or :pinv"
    @assert dampingTrim > 0 "dampingTrim must be positive"

    # Create nonlinear solver properties
    props = NLsolverProperties(
        method,
        linesearch,
        jcalc,
        checkJac,
        ftol,
        maxIter,
        desIter,
        displayStatus,
    )

    # Create lineartrim solver
    linearTrimSolver = TrimSolver(
        trimSolverMethod,
        dampingTrim,
    )

    # Create load stepping properties
    loadStepProps = LoadSteppingProperties(
        fullLoadFacTol,
        initLoadFac,
        minLoadFac,
        maxLoadFacStep,
        minLoadFacStep,
        trackLoadSteps,
    )

    return NLsolver(
        props=props,
        linearTrimSolver=linearTrimSolver,
        loadStepProps=loadStepProps,
    )

end
export create_NLsolver


"""
    create_ArclengthSolver(; kwargs...)

Arclength (Riks/Crisfield spherical) continuation solver constructor. Traces the equilibrium path in the (x,σ) space, so it can go through limit points (snap-through and snap-back), where load-controlled Newton methods fail. Meant for static (non-trim) problems, with the load factor σ scaling the applied loads.

# Keyword arguments
- `absTol::Float64`: absolute convergence tolerance for the residual
- `maxIter::Int`: maximum number of corrector iterations per arclength step
- `desIter::Int`: desired number of corrector iterations per arclength step, used to adapt the arclength step
- `initArcLength::Real`: initial arclength step, in the (x,σ) space
- `minArcLength::Float64`: minimum arclength step
- `maxArcLength::Real`: maximum arclength step
- `maxLoadFac::Real`: load factor at which the continuation stops (it may be greater than 1)
- `maxSteps::Int`: maximum number of arclength steps
- `initLoadFac::Real`: initial load factor
- `minLoadFac::Float64`: minimum load factor (only used to size the load step history)
- `trackLoadSteps::Bool`: flag to track the solution at each step
- `displayStatus::Bool`: flag to display status
"""
function create_ArclengthSolver(;
    absTol::Float64=1e-8,
    maxIter::Int=20,
    desIter::Int=5,
    minArcLength::Float64=1e-8,
    maxArcLength::Real=Inf,
    initArcLength::Real=min(1,maxArcLength),
    maxLoadFac::Real=1,
    maxSteps::Int=typemax(Int),
    initLoadFac::Real=0.0,
    minLoadFac::Float64=0.01,
    trackLoadSteps::Bool=true,
    displayStatus::Bool=false
    )

    # Validate inputs
    @assert absTol > 0 "absTol must be positive"
    @assert maxIter > 1 "maxIter must be greater than 1"
    @assert desIter > 1 "desIter must be greater than 1"
    @assert 0 < minArcLength <= initArcLength <= maxArcLength "arclengths must satisfy 0 < minArcLength <= initArcLength <= maxArcLength"
    @assert maxLoadFac > initLoadFac >= 0 "maxLoadFac must be greater than initLoadFac, which must be non-negative"
    @assert maxSteps > 0 "maxSteps must be positive"
    @assert 0 < minLoadFac <= 0.5 "minLoadFac must be between 0 and 0.5"

    # Create arclength solver properties
    props = ArclengthSolverProperties(
        absTol,
        maxIter,
        desIter,
        Float64(initArcLength),
        minArcLength,
        Float64(maxArcLength),
        Float64(maxLoadFac),
        maxSteps,
        displayStatus,
    )

    # Create load stepping properties (load factor limits are not used by the arclength procedure)
    loadStepProps = LoadSteppingProperties(
        1e-6,
        Float64(initLoadFac),
        minLoadFac,
        1.0,
        minLoadFac,
        trackLoadSteps,
    )

    return ArclengthSolver(
        props=props,
        loadStepProps=loadStepProps,
    )

end
export create_ArclengthSolver


## PRIVATE API
# ----------------------------------------------------------------------------
# Solves a linear problem
function solve_linear!(problem)

    # Assemble system arrays (residual and Jacobian)
    fj!(problem)

    # Solve linear system
    solve_linear_system!(problem)

    # Update states
    update_states!(problem)

    # Save loadstep data
    save_loadstep_data!(problem)

end

# Solves a linear system of equations
function solve_linear_system!(problem)

    # Solve according to problem type and set convergence flag
    problem.systemSolver.state.converged = problem isa TrimProblem ? 
                                           solve_trim_system!(problem) : 
                                           solve_regular_system!(problem)

end

# Factorizes the Jacobian into the problem's LU cache, reusing the symbolic analysis whenever possible
function factorize_jacobian!(problem, J; check::Bool=false)
    luCache = problem.systemArrays.luCache
    if isnothing(luCache)
        luCache = lu(J; check=check)
    else
        try
            lu!(luCache, J; check=check, reuse_symbolic=true)
        catch e
            e isa ArgumentError || rethrow()
            # Symbolic analysis failed (jacobian structure changed), recompute LU factorization from scratch
            luCache = lu(J; check=check)
        end
    end
    problem.systemArrays.luCache = luCache
    return luCache
end

# Solves the linear system of equations for regular (square Jacobian) problems
function solve_regular_system!(problem)

    @unpack x,Δx,residual,jacobian,luCache = problem.systemArrays

    # @info "norm(residual): $(norm(residual))"
    # @info "norm(jacobian): $(norm(jacobian))"
    # @info "det(jacobian): $(det(jacobian))"

    # LU factorize Jacobian
    luCache = factorize_jacobian!(problem, jacobian)

    # Return false if the Jacobian is singular
    issuccess(luCache) || begin
        @warn "Jacobian is singular"
        return false
    end

    # Solve the linear system for the state increment Δx
    ldiv!(Δx, luCache, residual)
    rmul!(Δx, -1)

    # Update states array
    x .+= Δx

    # @info "norm(Δx): $(norm(Δx))"
    # @info "Δx: $(round.(Δx, sigdigits=3))"

    return true

end

 # Solves the linear system of equations for trim (rectangular Jacobian) problems
function solve_trim_system!(problem)

    @unpack x,Δx = problem.systemArrays

    # Compute the Newton step for the trim problem
    trim_newton_step!(problem)

    # Update states array
    x .+= Δx

    return true

end

# Computes the Newton step for the trim problem
function trim_newton_step!(problem)

    @unpack method,damping = problem.systemSolver.linearTrimSolver
    @unpack Δx,residual,jacobian = problem.systemArrays
    @unpack J,JtJ,Jtr,J_aug,r_aug = problem.trimWorkspace

    # @info "norm(residual): $(norm(residual))"
    # @info "norm(jacobian): $(norm(jacobian))"

    # Dense Jacobian
    copyto!(J, jacobian)

    # Solve underdetermined system using selected method
    if method == :LQ
        solve_trim_LQ!(Δx,J,residual)
    elseif method == :LM
        solve_trim_LM!(Δx,J,residual,JtJ,Jtr,J_aug,r_aug,damping)
    elseif method == :SVD
        solve_trim_SVD!(Δx,J,residual)
    elseif method == :pinv
        solve_trim_pinv!(Δx,J,residual)
    end

    # @info "Δx: $(round.(Δx, sigdigits=3))"

    return Δx

end

# Solve for trim problem solution increment with LQ decomposition method
function solve_trim_LQ!(Δx,J,residual)
    F = lq!(J)
    Δx .= -(F \ residual)
end

# Solve for trim problem solution increment with simplified (non-adaptive damping) Levenberg-Marquardt method
function solve_trim_LM!(Δx,J,residual,JtJ,Jtr,J_aug,r_aug,damping)

    # Compute the damping term based on the norm of the Jacobian
    ϵ = damping * norm(J)

    # JtJ = J' * J + ϵ * I
    LinearAlgebra.BLAS.syrk!('U', 'T', 1.0, J, 0.0, JtJ)
    @inbounds for i in axes(JtJ,1)
        JtJ[i,i] += ϵ
    end

    # Jtr = J' * residual
    mul!(Jtr, J', residual)

    # Try to solve the normal equations with Cholesky factorization, if it fails, use LQ decomposition of the augmented system
    try
        F = cholesky!(Symmetric(JtJ, :U))
        Δx .= -(F \ Jtr)
    catch e
        e isa PosDefException || rethrow(e)
        m, n = size(J)
        @views J_aug[1:m, :] .= J
        for i in 1:n
            J_aug[m+i, i] = sqrt(ϵ)
        end
        @views r_aug[1:m] .= residual
        Δx .= -(J_aug \ r_aug)
    end

end

# Solve for trim problem solution increment with Singular Value Decomposition method
function solve_trim_SVD!(Δx,J,residual)
    F = svd!(J)
    Δx .= -(F \ residual)
end

# Solve for trim problem solution increment with Pseudo-Inverse (Moore-Penrose) method
function solve_trim_pinv!(Δx,J,residual)
    Δx .= - pinv(J) * residual
end

# Solves a nonlinear system of equations
function solve_nonlinear_system!(problem)

    # Arclength continuation has its own stepping procedure
    problem.systemSolver isa ArclengthSolver && return solve_arclength!(problem)

    # Unpack solution vector
    @unpack x = problem.systemArrays

    # Unpack system solver
    @unpack systemSolver = problem
    @unpack desIter,displayStatus = systemSolver.props

    # Unpack load stepping properties
    @unpack fullLoadFacTol,initLoadFac,minLoadFac,maxLoadFacStep,minLoadFacStep,trackLoadSteps = systemSolver.loadStepProps

    # Alias for load stepping state
    lss = systemSolver.loadStepState

    # Reset state
    lss.loadstep = 0
    lss.σ = initLoadFac
    lss.σStep = maxLoadFacStep
    lss.convergedFinalSolution = false

    # Pre-allocate known solution array
    if isempty(lss.xKnown)
        lss.xKnown = similar(x)
    end

    ## Load steps
    #---------------------------------------------------------------------------
    while !lss.convergedFinalSolution
        ## Initialization
        # ----------------------------------------------------------------------
        # Update loadstep counter
        lss.loadstep += 1
        # Current loadstep's known solution array and load factor
        copyto!(lss.xKnown, x)
        σKnown = lss.loadstep == 1 ? 
                 0.0 : 
                 lss.σ
        # Reset convergence flags
        lss.convergedLoadStep = false
        # Update current load factor
        lss.σ = lss.loadstep == 1 ? 
                lss.σ :
                clamp(lss.σ+lss.σStep, σKnown, 1.0)
        ## Solve current load step
        # ----------------------------------------------------------------------
        solve_load_step!(problem,lss,systemSolver)
        ## Load factor adjustment
        #-----------------------------------------------------------------------
        # In case of convergence trouble
        if !lss.convergedLoadStep
            # Reduce load factor step (impose bounds)
            lss.σStep = clamp(
                lss.σStep/2,
                minLoadFacStep,
                1-(σKnown+minLoadFac/lss.loadstep)
            )
            # Check if algorithm is stuck
            if isapprox(lss.σStep, minLoadFacStep; atol=1e-9)
                println("Load stepping algorithm stuck, stopping...")
                return
            end
            # Reset to previously known state
            copyto!(x,lss.xKnown)
            lss.σ = σKnown
            # Display status
            displayStatus && println("Unconverged, reducing load factor to $(clamp(lss.σ+lss.σStep, σKnown, 1.0)) and trying again...")
        # In case of converged partial solution    
        else
            # Update states of elements, hinge axis constraints and special nodes
            update_states!(problem)
            # Adjust load factor step (imposing bounds)
            lss.σStep = clamp(
                lss.σStep * sqrt(desIter/max(systemSolver.state.iter, 1)),
                minLoadFacStep,
                maxLoadFacStep
            )
            # Save solution at current load step, if tracking
            if trackLoadSteps && !(problem isa DynamicProblem)
                save_loadstep_data!(problem)
            end
            # Update flag when convergence with full load is reached within the specified tolerance
            if isapprox(lss.σ, 1.0; atol=fullLoadFacTol)
                lss.convergedFinalSolution = true
            end
        end
    end

end

# Solves a load step using NLSolve.jl
function solve_load_step!(problem,lss,systemSolver::NLsolver)

    @unpack x,residual,jacobian = problem.systemArrays
    @unpack method,linesearch,jcalc,checkJac,ftol,maxIter,displayStatus = systemSolver.props

    # Alias for solver state
    s = systemSolver.state

    # Initialize relative convergence norm
    s.ϵ_rel = 0.0

    # Jacobian check
    labels = checkJac ? 
             build_jacobian_labels(problem) : 
             nothing

    # Residual, Jacobian and combined functions for NLSolve.jl
    function r_nls!(F, x)
        copyto!(problem.systemArrays.x, x)
        R = f!(problem)
        F === R || copyto!(F, R)
        return F
    end

    function j_nls!(J, x)
        copyto!(problem.systemArrays.x, x)
        checkJac && f!(problem)
        j!(problem)
        checkJac && compare_jacobians(problem, J, labels; tag="(loadstep $(lss.loadstep), σ=$(round(lss.σ,digits=4))")
        return J
    end

    function fj_nls!(F, J, x)
        copyto!(problem.systemArrays.x, x)
        fj!(problem)
        checkJac && compare_jacobians(problem, J, labels; tag="(loadstep $(lss.loadstep), σ=$(round(lss.σ,digits=4))")
        return nothing
    end

    # Linear solver function for NLSolve.jl (also estimates relative convergence norm)
    function linsolve_nls!(Δx, J, F)
        if problem isa TrimProblem
            step = trim_newton_step!(problem)
            Δx .= .-step  # NLsolve uses inverted sign convention
        else
            luCache = factorize_jacobian!(problem, J; check=false)
            issuccess(luCache) || @warn "Jacobian is singular"
            ldiv!(Δx, luCache, F)
            # rmul!(Δx, -1) # Not necessary: NLsolve uses inverted sign convention
        end
        s.ϵ_rel = norm(Δx) / max(norm(problem.systemArrays.x), eps())
    end

    # Keyword arguments for NLSolve.jl based on the chosen method
    kw = if method == :newton && linesearch
        (; linesearch=LineSearches.BackTracking() )
    else
        (;)
    end

    # Solve based on the chosen Jacobian calculation method
    local sol,df
    try
        if jcalc == :analytical
            # Create OnceDifferentiable object to exploit Jacobian sparsity
            df = NLsolve.OnceDifferentiable(r_nls!, j_nls!, fj_nls!, x, residual, jacobian)
            # OnceDifferentiable stores copies of the residual and Jacobian: alias them to the problem's arrays so they are assembled in place
            df.F = residual
            problem.systemArrays.jacobian = df.DF
            # Solve using NLSolve.jl
            sol = NLsolve.nlsolve(df, x; method=method, ftol=ftol, iterations=maxIter, show_trace=displayStatus, linsolve=linsolve_nls!, kw...)
            # Assemble the Jacobian at the solution and factorize it (not needed in dynamic problems, where it is rebuilt at every time step)
            if !(problem isa DynamicProblem)
                NLsolve.jacobian!(df, sol.zero)
                factorize_jacobian!(problem, problem.systemArrays.jacobian; check=!(problem isa TrimProblem))
            end
        else#if jcalc == :numerical
            # Create OnceDifferentiable object with numerical (FD) Jacobian
            df = NLsolve.OnceDifferentiable(r_nls!, x, residual; autodiff=ADTypes.AutoFiniteDiff(fdtype=Val(:central)))
            # Solve using NLSolve.jl
            sol = NLsolve.nlsolve(df, x; method=method, ftol=ftol, iterations=maxIter, show_trace=displayStatus, kw...)
            # Assemble the Jacobian
            problem.systemArrays.jacobian = NLsolve.jacobian!(df, copy(sol.zero))
        end
    catch e
        # Diverging solution leads to InexactError
        e isa InexactError || rethrow()
        sol = nothing
    end

    # Diverged: flag as unconverged so the load stepping procedure reduces the step
    if isnothing(sol)
        s.converged = false
        displayStatus && println("NLSolve.jl diverged with σ: $(lss.σ)")
        return
    end

    # Update convergence flag, number of iterations and residual convergence norm
    s.converged = sol.f_converged
    s.iter = sol.iterations
    s.ϵ_abs = sol.residual_norm

    # Update solution and load step convergence flag, if applicable
    if s.converged
        copyto!(x, sol.zero)
        lss.convergedLoadStep = true
        displayStatus && println("Time step $(time_index(problem)), load step $(lss.loadstep), σ: $(lss.σ), i: $(s.iter), ϵ_abs: $(round(s.ϵ_abs,sigdigits=6)), ϵ_rel ≈ $(round(s.ϵ_rel,sigdigits=6))")
    else
        displayStatus && println("NLSolve.jl did not converge with σ: $(lss.σ)")
    end

end

# Solves a load step using custom Newton-Raphson method
function solve_load_step!(problem,lss,systemSolver::NewtonRaphson)

    @unpack x,Δx,residual = problem.systemArrays
    @unpack rangeAeroStates = problem.model.assemblyData
    @unpack method,absTol,relTol,maxIter,maxAbsErr,maxRelErr,displayStatus,minConvRateJacUpdate,minConvRateAeroJacUpdate,allowAdvanceThroughUnconvergedAeroStates,aeroResRatioThreshold,aeroΔxErrorThreshold = systemSolver.props

    # Alias for solver state
    s = systemSolver.state

    # Reset Jacobian update flags, iteration counter, and previous absolute convergence norm
    s.skipJacobianUpdate = false
    s.skipAeroJacobianUpdate = false
    s.iter = 0
    s.ϵ_abs_previous = 1.0
    
    # Newton-Raphson iteration loop
    while !lss.convergedLoadStep && s.iter < maxIter

        # Update iteration count
        s.iter += 1

        # Assemble system arrays (residual and Jacobian) at current time step and load factor
        fj!(problem)

        # Solve the linear system (and check linear solution convergence flag)
        solve_linear_system!(problem)
        s.converged || return

        # Residual (absolute) and solution (relative) convergence norms
        s.ϵ_abs = norm(residual, Inf)
        s.ϵ_rel = norm(Δx)/max(norm(x), eps())

        # Convergence rate
        s.convRate = s.ϵ_abs_previous / s.ϵ_abs

        # Update previous value of absolute convergence norm
        s.ϵ_abs_previous = s.ϵ_abs

        # Flags to skip Jacobian update
        if method == :NewtonChord && s.iter > 1
            s.skipJacobianUpdate = s.convRate > minConvRateJacUpdate
            s.skipAeroJacobianUpdate = s.convRate > minConvRateAeroJacUpdate
        end

        # Display status, if applicable
        displayStatus && println("Time step $(time_index(problem)), load step $(lss.loadstep), σ: $(lss.σ), i: $(s.iter), ϵ_abs: $(round(s.ϵ_abs,sigdigits=6)), ϵ_rel: $(round(s.ϵ_rel,sigdigits=6))")

        # Check convergence
        if s.ϵ_abs < absTol && s.ϵ_rel < relTol
            lss.convergedLoadStep = true
        elseif s.ϵ_abs > maxAbsErr || s.ϵ_rel > maxRelErr || !isfinite(s.ϵ_abs)
            displayStatus && println("Diverging solution...")
            return     
        elseif !isempty(rangeAeroStates) && allowAdvanceThroughUnconvergedAeroStates && s.iter == maxIter 
            # Compute error ratios for aerodynamic states
            aeroResErrorRatio = norm(residual[rangeAeroStates], Inf) / s.ϵ_abs
            aeroΔxErrorRatio = norm(Δx[rangeAeroStates]) / norm(Δx)
            if aeroResErrorRatio >= aeroResRatioThreshold && aeroΔxErrorRatio >= aeroΔxErrorThreshold
                # Admit converged partial solution
                lss.convergedLoadStep = true
                @info("Advancing through unconverged aero states (ϵ_abs: $(round(s.ϵ_abs,sigdigits=6)), ϵ_rel: $(round(s.ϵ_rel,sigdigits=6)), ϵ_absAero/ϵ_abs=$(round(aeroResErrorRatio, sigdigits=4)), ϵ_relAero/ϵ_rel=$(round(aeroΔxErrorRatio, sigdigits=4)))...")
            end
        end

    end

end

# Solves a nonlinear problem by arclength continuation (tangent predictor and spherical-constraint Newton corrector)
function solve_arclength!(problem)

    @unpack x,residual,jacobian = problem.systemArrays
    @unpack systemSolver = problem
    @unpack absTol,maxIter,desIter,minArcLength,maxArcLength,maxLoadFac,maxSteps,displayStatus = systemSolver.props
    @unpack initLoadFac,trackLoadSteps = systemSolver.loadStepProps
    lss = systemSolver.loadStepState
    s = systemSolver.state

    # Reset state
    lss.loadstep = 0
    lss.σ = initLoadFac
    lss.convergedFinalSolution = false
    s.arcLength = systemSolver.props.initArcLength
    s.converged = false
    isempty(lss.xKnown) && (lss.xKnown = similar(x))

    # Solution at the initial load factor (load control)
    newton_fixed_load!(problem)
    if !s.converged
        println("Arclength solver: initial load factor solution did not converge, stopping...")
        return
    end
    update_states!(problem)
    lss.loadstep = 1
    trackLoadSteps && save_loadstep_data!(problem)

    # Continuation
    xPrev = zeros(length(x))
    σPrev = 0.0
    v = similar(x)
    while lss.loadstep < maxSteps && !lss.convergedFinalSolution

        lss.loadstep += 1
        copyto!(lss.xKnown, x)
        σKnown = lss.σ

        # Unit tangent to the equilibrium path at the known point
        fj!(problem)
        R = copy(residual)
        q = residual_load_derivative(problem, R)
        luCache = factorize_jacobian!(problem, jacobian)
        issuccess(luCache) || (@warn "Jacobian is singular at load factor $σKnown, stopping..."; return)
        v .= .-(luCache \ q)
        nrm = sqrt(dot(v,v)+1)
        tx = v ./ nrm
        tσ = 1/nrm
        # Keep marching in the same direction as the previous step
        if lss.loadstep > 2 && dot(tx,xPrev) + tσ*σPrev < 0
            tx .*= -1
            tσ = -tσ
        end

        # Try the current arclength, halving it upon failure
        s.converged = false
        while !s.converged
            x .= lss.xKnown .+ s.arcLength .* tx
            lss.σ = σKnown + s.arcLength*tσ
            s.iter = 0
            while s.iter < maxIter
                s.iter += 1
                fj!(problem)
                R = copy(residual)
                s.ϵ_abs = norm(R, Inf)
                displayStatus && println("Arclength step $(lss.loadstep), σ: $(lss.σ), i: $(s.iter), ϵ_abs: $(round(s.ϵ_abs,sigdigits=6))")
                if s.ϵ_abs < absTol
                    s.converged = true
                    break
                elseif !isfinite(s.ϵ_abs)
                    break
                end
                q = residual_load_derivative(problem, R)
                luCache = factorize_jacobian!(problem, jacobian)
                issuccess(luCache) || break
                dx1 = .-(luCache \ R)
                dx2 = .-(luCache \ q)
                # Linearized spherical constraint: |Δx|² + Δσ² = Δs²
                Δx = x .- lss.xKnown
                Δσ = lss.σ - σKnown
                g = dot(Δx,Δx) + Δσ^2 - s.arcLength^2
                λ = (-g/2 - dot(Δx,dx1)) / (dot(Δx,dx2) + Δσ)
                x .+= dx1 .+ λ .* dx2
                lss.σ += λ
            end
            if !s.converged
                s.arcLength /= 2
                if s.arcLength < minArcLength
                    println("Arclength step too small, stopping...")
                    copyto!(x, lss.xKnown)
                    lss.σ = σKnown
                    return
                end
            end
        end

        # If the step overshot the maximum load factor, return to it by load-controlled Newton iterations (retry with smaller step if it fails)
        if lss.σ > maxLoadFac
            t = (maxLoadFac-σKnown)/(lss.σ-σKnown)
            x .= lss.xKnown .+ t .* (x .- lss.xKnown)
            lss.σ = maxLoadFac
            if !newton_fixed_load!(problem)
                s.arcLength /= 2
                if s.arcLength < minArcLength
                    println("Arclength step too small, stopping...")
                    copyto!(x, lss.xKnown)
                    lss.σ = σKnown
                    return
                end
                copyto!(x, lss.xKnown)
                lss.σ = σKnown
                lss.loadstep -= 1
                continue
            end
        end

        # Accept the step and adapt the arclength
        σPrev = lss.σ - σKnown
        xPrev .= x .- lss.xKnown
        update_states!(problem)
        trackLoadSteps && save_loadstep_data!(problem)
        s.arcLength = clamp(s.arcLength*sqrt(desIter/s.iter), minArcLength, maxArcLength)
        lss.convergedFinalSolution = lss.σ >= maxLoadFac

    end

    lss.convergedFinalSolution || println("Arclength solver stopped before reaching maxLoadFac")

end

# Newton iterations at the current (fixed) load factor, starting from the current state. Returns the convergence flag
function newton_fixed_load!(problem)

    @unpack x,residual,jacobian = problem.systemArrays
    @unpack absTol,maxIter = problem.systemSolver.props
    s = problem.systemSolver.state

    s.iter = 0
    s.converged = false
    while s.iter <= maxIter
        fj!(problem)
        s.ϵ_abs = norm(residual, Inf)
        if s.ϵ_abs < absTol
            s.converged = true
            break
        end
        s.iter += 1
        luCache = factorize_jacobian!(problem, jacobian)
        issuccess(luCache) || (@warn "Jacobian is singular"; break)
        x .-= luCache \ residual
    end

    return s.converged

end

# Derivative of the residual with respect to the load factor (forward finite differences), given the residual R at the current state
function residual_load_derivative(problem, R)

    lss = problem.systemSolver.loadStepState
    h = 1e-6
    σ0 = lss.σ
    lss.σ = σ0 + h
    q = (f!(problem) .- R) ./ h
    lss.σ = σ0
    return q

end

# Resets to default solver state
function reset_solver_state!(systemSolver)

    if systemSolver isa LinearSolver
        setfield!(systemSolver, :state, LinearSolverState())
    elseif systemSolver isa NewtonRaphson
        setfield!(systemSolver, :state, NRState())
        setfield!(systemSolver, :loadStepState, LoadSteppingState(systemSolver.loadStepProps))
    elseif systemSolver isa NLsolver
        setfield!(systemSolver, :state, NLsolverState())
        setfield!(systemSolver, :loadStepState, LoadSteppingState(systemSolver.loadStepProps))
    elseif systemSolver isa ArclengthSolver
        setfield!(systemSolver, :state, ArclengthSolverState(systemSolver.props.initArcLength))
        setfield!(systemSolver, :loadStepState, LoadSteppingState(systemSolver.loadStepProps))
    end

end

# Accessor functions for Jacobian update flags
get_skipJacobianUpdate(systemSolver::NewtonRaphson) = systemSolver.state.skipJacobianUpdate
get_skipAeroJacobianUpdate(systemSolver::NewtonRaphson) = systemSolver.state.skipAeroJacobianUpdate
get_skipJacobianUpdate(systemSolver) = false
get_skipAeroJacobianUpdate(systemSolver) = false

function build_jacobian_labels(problem)
    n = length(problem.systemArrays.x)
    @unpack elements,specialNodes = problem.model.assemblyData
    rowKind = fill("?", n); rowOwner = zeros(Int, n)
    colKind = fill("?", n); colOwner = zeros(Int, n)

    function tag!(kinds, owners, idx, kind, id)
        for i in idx
            if kinds[i] == "?"
                kinds[i] = kind; owners[i] = id
            end
        end
    end

    for e in elements
        id = e.assemblyData.globalID
        a = e.assemblyIndices
        for (k, eq) in (("Fu",a.eqs_Fu1),("Fu",a.eqs_Fu2),("Fp",a.eqs_Fp1),("Fp",a.eqs_Fp2),
                        ("FF",a.eqs_FF1),("FF",a.eqs_FF2),("FM",a.eqs_FM1),("FM",a.eqs_FM2),
                        ("FV",a.eqs_FV),("FΩ",a.eqs_FΩ),("Fχ",a.eqs_Fχ))
            tag!(rowKind, rowOwner, eq, k, id)
        end
        for (k, dof) in (("u",a.DOF_u),("p",a.DOF_p),("F",a.DOF_F),("M",a.DOF_M),
                         ("V",a.DOF_V),("Ω",a.DOF_Ω),("χ",a.DOF_χ))
            tag!(colKind, colOwner, dof, k, id)
        end
    end

    for sn in specialNodes
        id = sn.assemblyData.globalID
        d = sn.assemblyData
        tag!(rowKind, rowOwner, d.eqs_Fu, "Fu_node", id)
        tag!(rowKind, rowOwner, d.eqs_Fp, "Fp_node", id)
        for e in eachindex(d.eqs_FF)
            tag!(rowKind, rowOwner, d.eqs_FF[e], "FF_node", id)
            tag!(rowKind, rowOwner, d.eqs_FM[e], "FM_node", id)
        end
        tag!(colKind, colOwner, d.DOF_uF, "uF_node", id)
        tag!(colKind, colOwner, d.DOF_pM, "pM_node", id)
        # tag!(colKind, colOwner, d.DOF_T,  "T_node",  id)
    end

    for c in problem.model.hingeAxisConstraints
        tag!(rowKind, rowOwner, c.assemblyData.eqs_λ, "C", 0)
        tag!(colKind, colOwner, c.assemblyData.DOF_λ,  "λ", 0)
    end

    return (; rowKind, rowOwner, colKind, colOwner)
end

function compare_jacobians(problem, Jan, labels; h0=1e-8, rtol=1e-4, atol=1e-4, ntop=10, tag="")
    xs = problem.systemArrays.x
    x0 = copy(xs)
    n  = length(x0)
    bad  = Tuple{Float64,Int,Int,Float64,Float64}[]           # (erro, i, k, Jan, Jfd)
    blocks = Dict{Tuple{String,String},Float64}()             # erro máx por bloco

    for k in 1:n
        h = h0 * max(1.0, abs(x0[k]))
        xs .= x0; xs[k] = x0[k] + h; fp = copy(f!(problem))
        xs .= x0; xs[k] = x0[k] - h; fm = copy(f!(problem))
        colFD = (fp .- fm) ./ (2h)
        colAN = Vector(Jan[:, k])
        for i in 1:n
            err = abs(colAN[i] - colFD[i])
            if err > atol + rtol * max(abs(colAN[i]), abs(colFD[i]))
                push!(bad, (err, i, k, colAN[i], colFD[i]))
                key = (labels.rowKind[i], labels.colKind[k])
                blocks[key] = max(get(blocks, key, 0.0), err)
            end
        end
    end

    xs .= x0
    f!(problem)                     # restaura caches em x0

    println("── Checagem da Jacobiana $tag: $(length(bad)) entradas fora da tolerância")
    for (key, err) in sort(collect(blocks); by=last, rev=true)
        println("   bloco F_", key[1], "_", key[2], "   erro máx = ", round(err, sigdigits=3))
    end
    sort!(bad; by=first, rev=true)
    for (err, i, k, a, d) in first(bad, ntop)
        println("   row=$i ($(labels.rowKind[i]), id=$(labels.rowOwner[i]))  ",
                "col=$k ($(labels.colKind[k]), id=$(labels.colOwner[k]))  ",
                "analítica=$(round(a,sigdigits=5))  FD=$(round(d,sigdigits=5))")
    end
    return bad
end