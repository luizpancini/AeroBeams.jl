## INTEGRATION METHOD COMPOSITE TYPES
# ------------------------------------------------------------------------------
"""

    Newmark IntegrationMethod composite type

Newmark method (generalized trapezoidal rule) for the time integration of dynamic problems. The equations are solved at the end of the time step, where the rates of the states are given by the increment of the states and by the rates at the begin of the time step, `ẋₙ₊₁ = (xₙ₊₁ - xₙ)/(γΔt) - (1-γ)/γ ẋₙ`. With `γ = 1/2`, this is the trapezoidal rule (Newmark's average acceleration method), which is second-order accurate and has no numerical dissipation: the modes whose frequencies are not resolved by the time step (such as those of very stiff parts of the structure) are not damped, and may grow in nonlinear problems. Values of `γ` greater than 1/2 damp them, at the cost of first-order accuracy (`γ = 1` is the backward Euler method)

# Keyword arguments
- `γ::Real`: Newmark's γ parameter, between 1/2 and 1 (default: 1/2)
"""
struct Newmark

    name::Symbol
    γ::Float64

    function Newmark(; γ::Real=1/2)

        # Validate
        @assert 1/2 <= γ <= 1 "γ must be between 1/2 and 1"

        return new(:newmark, Float64(γ))

    end

end
export Newmark


# BDF2State composite type
# Fields
# - `xPrevious`: states at the begin of the previous time step
# - `ΔtPrevious`: previous time step (NaN if there is no previous time step)
#
mutable struct BDF2State

    xPrevious::Vector{Float64}
    ΔtPrevious::Float64

end
BDF2State() = BDF2State(Vector{Float64}(), NaN)


"""

    BDF2 IntegrationMethod composite type

Second-order backward differentiation formula (BDF2) for the time integration of dynamic problems. The equations are solved at the end of the time step, where the rates of the states are given by the states at the end of the time step and at the begin of the current and previous time steps, `ẋₙ₊₁ = [(1+2ω)/(1+ω) xₙ₊₁ - (1+ω) xₙ + ω²/(1+ω) xₙ₋₁]/Δt`, in which `ω = Δt/Δtₙ₋₁` is the ratio of the current and previous time steps (with a constant time step, `ẋₙ₊₁ = (3xₙ₊₁ - 4xₙ + xₙ₋₁)/(2Δt)`). The method is second-order accurate and L-stable: the modes whose frequencies are not resolved by the time step are strongly damped. The modes that are resolved are also damped, with an algorithmic damping ratio of about (ωΔt)³/4 (e.g., 2.4e-4 for ωΔt = 0.1, but 0.07 for ωΔt = 1), which may be significant for the modes that are not well resolved. At the first time step, the states at the begin of the previous time step are extrapolated from the initial states and rates, `xₙ₋₁ = xₙ - Δt ẋₙ`
"""
struct BDF2

    name::Symbol
    state::BDF2State

    function BDF2()

        return new(:bdf2, BDF2State())

    end

end
export BDF2


# GeneralizedAlphaState composite type
# Fields
# - `xBegin`: states at the begin of the time step
# - `xdotBegin`: rates of the states at the begin of the time step (null for the states without rates)
# - `dynamicStallMemoryBegin`: memory (complementary) variables of the dynamic stall models of the elements at the begin of the time step (nothing for the elements without them)
#
struct GeneralizedAlphaState

    xBegin::Vector{Float64}
    xdotBegin::Vector{Float64}
    dynamicStallMemoryBegin::Vector{Union{Nothing,BLiComplementaryVariables{Float64},BLoComplementaryVariables{Float64}}}

end
GeneralizedAlphaState() = GeneralizedAlphaState(Vector{Float64}(), Vector{Float64}(), Vector{Union{Nothing,BLiComplementaryVariables{Float64},BLoComplementaryVariables{Float64}}}())


"""

    GeneralizedAlpha IntegrationMethod composite type

Generalized-α method for first-order systems (Jansen, Whiting and Hulbert, 2000, https://doi.org/10.1016/S0045-7825(00)00203-6) for the time integration of dynamic problems. The equations are solved within the time step, at the time `tₙ₊αf = tₙ + αf Δt`, with the states `xₙ₊αf = xₙ + αf (xₙ₊₁ - xₙ)` and the rates `ẋₙ₊αm = ẋₙ + αm (ẋₙ₊₁ - ẋₙ)`, and the rates at the end of the time step are given by Newmark's relation, `ẋₙ₊₁ = (xₙ₊₁ - xₙ)/(γΔt) - (1-γ)/γ ẋₙ`. The parameters follow from the spectral radius at infinite frequency `ρ∞`: `αm = (3-ρ∞)/(2(1+ρ∞))`, `αf = 1/(1+ρ∞)` and `γ = 1/2 + αm - αf`, for which the method is second-order accurate and unconditionally stable (in linear problems), and the modes whose frequencies are not resolved by the time step are damped by the factor `ρ∞` at each time step: `ρ∞ = 1` gives no numerical dissipation (midpoint rule), and `ρ∞ = 0` annihilates those modes in a single time step (it is equivalent to BDF2 in linear problems). The numerical dissipation of the modes resolved by the time step is much smaller than with BDF2 (for `ρ∞ = 0.8`, an algorithmic damping ratio of about 3e-7 for ωΔt = 0.1 and 3e-4 for ωΔt = 1). The time-dependent inputs (motion of basis A, BCs, distributed loads and flap deflections) are linearly interpolated at `tₙ₊αf`. The memory variables of the dynamic stall models (such as the stall onset and vortex shedding events) are only updated with the states at the end of the time step

# Keyword arguments
- `ρ∞::Real`: spectral radius at infinite frequency, between 0 and 1 (default: 0.8)
"""
struct GeneralizedAlpha

    name::Symbol
    ρ∞::Float64
    αm::Float64
    αf::Float64
    γ::Float64
    state::GeneralizedAlphaState

    function GeneralizedAlpha(; ρ∞::Real=0.8)

        # Validate
        @assert 0 <= ρ∞ <= 1 "ρ∞ must be between 0 and 1"

        # Parameters for second-order accuracy and spectral radius at infinite frequency ρ∞
        αm = (3-ρ∞)/(2*(1+ρ∞))
        αf = 1/(1+ρ∞)
        γ = 1/2 + αm - αf

        return new(:gen_alpha, Float64(ρ∞), αm, αf, γ, GeneralizedAlphaState())

    end

end
export GeneralizedAlpha


const IntegrationMethod = Union{Newmark,BDF2,GeneralizedAlpha}


## PRIVATE API
# ------------------------------------------------------------------------------
# Gets the integration method from its name (`:newmark`, `:bdf2` or `:gen_alpha`, with the default parameters), or returns it if it is already an integration method
integration_method(method::IntegrationMethod) = method
function integration_method(name::Symbol)

    @assert name in (:newmark,:bdf2,:gen_alpha) "integrationMethod must be :newmark, :bdf2 or :gen_alpha (or one of the integration methods Newmark, BDF2 or GeneralizedAlpha)"

    return name == :newmark ? Newmark() :
           name == :bdf2 ? BDF2() :
           GeneralizedAlpha()

end


# Fraction of the time step at which the equations are solved
evaluation_step_fraction(method::GeneralizedAlpha) = method.αf
evaluation_step_fraction(method) = 1.0


# Resets the history of the integration method, at the begin of the time march of a problem with nStates states
reset_integration_method!(method::Newmark,nStates) = nothing
function reset_integration_method!(method::BDF2,nStates)

    resize!(method.state.xPrevious, nStates)
    method.state.ΔtPrevious = NaN

end
function reset_integration_method!(method::GeneralizedAlpha,nStates)

    resize!(method.state.xBegin, nStates)
    resize!(method.state.xdotBegin, nStates)
    fill!(method.state.xdotBegin, 0)
    empty!(method.state.dynamicStallMemoryBegin)

end


# Begins the current time step, after the update of its equivalent states' rates and of the complementary variables of the dynamic stall models: saves the data at the begin of the time step that are needed to complete it. Only the generalized-α method, whose equations are solved within the time step, needs them
begin_time_step!(problem,method) = nothing
function begin_time_step!(problem,method::GeneralizedAlpha)

    @unpack xBegin,xdotBegin,dynamicStallMemoryBegin = method.state
    @unpack x = problem.systemArrays
    @unpack elements = problem.model.assemblyData

    # States and their rates at the begin of the time step
    xBegin .= x
    set_derivative_from_elements_rates!(xdotBegin,elements)

    # Memory variables of the dynamic stall models at the begin of the time step
    resize!(dynamicStallMemoryBegin, length(elements))
    for (e,element) in enumerate(elements)
        dynamicStallMemoryBegin[e] = dynamic_stall_memory(element)
    end

end


# Gets a copy of the memory (complementary) variables of the dynamic stall model of the element, or nothing if the element has none
function dynamic_stall_memory(element)

    @unpack surface,stateCache = element.aero

    isnothing(surface) && return nothing
    surface.solvers.ppSolver isa BLi && return copy(stateCache.bliVars.comp)
    surface.solvers.ppSolver isa BLo && return copy(stateCache.bloVars.comp)

    return nothing

end


# Restores the memory (complementary) variables of the dynamic stall model of the element
restore_dynamic_stall_memory!(element,memory::Nothing) = nothing
restore_dynamic_stall_memory!(element,memory::BLiComplementaryVariables) = copy_mutable_fields!(element.aero.stateCache.bliVars.comp,memory)
restore_dynamic_stall_memory!(element,memory::BLoComplementaryVariables) = copy_mutable_fields!(element.aero.stateCache.bloVars.comp,memory)


# Copies the values of the fields of a mutable struct into another one of the same type
function copy_mutable_fields!(dest::T,src::T) where T

    for field in fieldnames(T)
        setfield!(dest, field, getfield(src, field))
    end

    return dest

end


# Updates the equivalent states' rates of the elements at the begin of the current time step, such that the rates of their states during the time step are ẋ = c x - ẋEquiv, according to the integration method
function update_equivalent_states_rates!(problem)

    update_equivalent_states_rates!(problem,problem.integrationMethod)

end

# Newmark method: ẋₙ₊₁ = c xₙ₊₁ - (c xₙ + (1-γ)/γ ẋₙ), with c = 1/(γΔt)
function update_equivalent_states_rates!(problem,method::Newmark)

    @unpack Δt = problem.time.state
    @unpack elementGroups = problem.model.assemblyData
    @unpack γ = method

    # Coefficients of the states and of their rates at the begin of the time step
    c = 1/(γ*Δt)
    b = (1-γ)/γ

    # Loop over elements and update equivalent states' rates
    foreach_element(elementGroups) do element
        set_element_equivalent_states_rates!(element,c,c,b)
    end

end

# BDF2 method: ẋₙ₊₁ = c xₙ₊₁ - (a xₙ - d xₙ₋₁), with c = (1+2ω)/((1+ω)Δt), a = (1+ω)/Δt and d = ω²/((1+ω)Δt), where ω = Δt/Δtₙ₋₁
function update_equivalent_states_rates!(problem,method::BDF2)

    @unpack xPrevious,ΔtPrevious = method.state
    @unpack x = problem.systemArrays
    @unpack Δt = problem.time.state
    @unpack elementGroups = problem.model.assemblyData

    # First time step: states at the begin of the previous time step extrapolated from the current ones, xₙ₋₁ = xₙ - Δt ẋₙ, so that ẋₙ₊₁ = c xₙ₊₁ - (c xₙ + 1/2 ẋₙ), with c = 3/(2Δt)
    if isnan(ΔtPrevious)
        c = 3/(2Δt)
        foreach_element(elementGroups) do element
            set_element_equivalent_states_rates!(element,c,c,1/2)
        end
    # Next time steps: states at the begin of the previous time step
    else
        ω = Δt/ΔtPrevious
        c = (1+2ω)/((1+ω)*Δt)
        a = (1+ω)/Δt
        d = ω^2/((1+ω)*Δt)
        foreach_element(elementGroups) do element
            set_element_BDF2_equivalent_states_rates!(element,c,a,d,xPrevious)
        end
    end

    # Save the states at the begin of the time step and the time step, for the next time step
    xPrevious .= x
    method.state.ΔtPrevious = Δt

end

# Generalized-α method: ẋₙ₊αm = c xₙ₊αf - (c xₙ - (1-αm/γ) ẋₙ), with c = αm/(γ αf Δt)
function update_equivalent_states_rates!(problem,method::GeneralizedAlpha)

    @unpack αm,αf,γ = method
    @unpack Δt = problem.time.state
    @unpack elementGroups = problem.model.assemblyData

    # Coefficients of the states and of their rates at the begin of the time step
    c = αm/(γ*αf*Δt)
    b = -(1-αm/γ)

    # Loop over elements and update equivalent states' rates
    foreach_element(elementGroups) do element
        set_element_equivalent_states_rates!(element,c,c,b)
    end

end


# Sets the coefficient c of the states in their rates and the equivalent states' rates of the element, ẋEquiv = a xₙ + b ẋₙ, from its states (xₙ) and rates (ẋₙ) at the begin of the time step, so that the rates of the states during the time step are ẋ = c x - ẋEquiv
function set_element_equivalent_states_rates!(element,c,a,b)

    @unpack u,p,V,Ω,χ = element.state.midpoint
    @unpack udot,pdot,Vdot,Ωdot,χdot = element.state.rates
    @unpack udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,χdotEquiv = element.stateCache.equivStatesRates

    # Set coefficient of the states
    setfield!(element.stateCache.equivStatesRates, :c, c)

    # Set equivalent rates
    udotEquiv .= b * udot + a * u
    pdotEquiv .= b * pdot + a * p
    VdotEquiv .= b * Vdot + a * V
    ΩdotEquiv .= b * Ωdot + a * Ω
    χdotEquiv .= b * χdot + a * χ

end

# Same as above, with the states and rates of the element at the begin of the time step taken from the vectors of states (x) and of their time derivatives (xdot) of the system
function set_element_equivalent_states_rates!(element,c,a,b,x,xdot)

    @unpack DOF_u,DOF_p,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
    @unpack udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,χdotEquiv = element.stateCache.equivStatesRates

    # Set coefficient of the states
    setfield!(element.stateCache.equivStatesRates, :c, c)

    # Set equivalent rates
    udotEquiv .= b * xdot[DOF_u] + a * x[DOF_u]
    pdotEquiv .= b * xdot[DOF_p] + a * x[DOF_p]
    VdotEquiv .= b * xdot[DOF_V] + a * x[DOF_V]
    ΩdotEquiv .= b * xdot[DOF_Ω] + a * x[DOF_Ω]
    χdotEquiv .= b * xdot[DOF_χ] + a * x[DOF_χ]

end

# Sets the coefficient c of the states in their rates and the equivalent states' rates of the element for the BDF2 method, ẋEquiv = a xₙ - d xₙ₋₁, from its states at the begin of the time step (xₙ) and at the begin of the previous time step (xₙ₋₁, taken from the vector of states of the system xPrevious)
function set_element_BDF2_equivalent_states_rates!(element,c,a,d,xPrevious)

    @unpack u,p,V,Ω,χ = element.state.midpoint
    @unpack DOF_u,DOF_p,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
    @unpack udotEquiv,pdotEquiv,VdotEquiv,ΩdotEquiv,χdotEquiv = element.stateCache.equivStatesRates

    # Set coefficient of the states
    setfield!(element.stateCache.equivStatesRates, :c, c)

    # Set equivalent rates
    udotEquiv .= a * u - d * xPrevious[DOF_u]
    pdotEquiv .= a * p - d * xPrevious[DOF_p]
    VdotEquiv .= a * V - d * xPrevious[DOF_V]
    ΩdotEquiv .= a * Ω - d * xPrevious[DOF_Ω]
    χdotEquiv .= a * χ - d * xPrevious[DOF_χ]

end


# Sets the vector of time derivatives of the states of the system (xdot) with the rates of the states of the elements (the null entries are those of the states without rates)
function set_derivative_from_elements_rates!(xdot,elements)

    for element in elements
        @unpack udot,pdot,Vdot,Ωdot,χdot = element.state.rates
        @unpack DOF_u,DOF_p,DOF_V,DOF_Ω,DOF_χ = element.assemblyIndices
        xdot[DOF_u] .= udot
        xdot[DOF_p] .= pdot
        xdot[DOF_V] .= Vdot
        xdot[DOF_Ω] .= Ωdot
        xdot[DOF_χ] .= χdot
    end

end


# Completes the current time step, after the solution of the equations: sets the states and their rates at the end of the time step. With the Newmark and BDF2 methods, the equations are solved at the end of the time step, so there is nothing to do
complete_time_step!(problem,method) = nothing

# Generalized-α method: the equations are solved at tₙ₊αf, so the states at the end of the time step are extrapolated from those at its begin and at tₙ₊αf, xₙ₊₁ = xₙ + (xₙ₊αf - xₙ)/αf, the rates are given by Newmark's relation, and the states of the elements, hinge axis constraints and special nodes are updated with the time-dependent inputs at the end of the time step. The memory variables of the dynamic stall models are restored to those at the begin of the time step before this update, so that they evolve from the begin to the end of the time step
function complete_time_step!(problem,method::GeneralizedAlpha)

    @unpack αf,γ = method
    @unpack xBegin,xdotBegin,dynamicStallMemoryBegin = method.state
    @unpack x = problem.systemArrays
    @unpack timeVector = problem.time
    @unpack timeIndex,Δt = problem.time.state
    @unpack beams = problem.model
    @unpack elements,elementGroups = problem.model.assemblyData

    # States at the end of the time step
    x .= xBegin .+ (x .- xBegin) ./ αf

    # Time variables at the end of the time step
    setfield!(problem.time.state, :timeNow, timeVector[timeIndex])
    setfield!(problem.time.state, :stepFraction, 1.0)

    # Motion of basis A, BCs and flap deflections at the end of the time step
    update_motion_basis_A!(problem,timeIndex)
    update_BCs_state!(problem,timeIndex)
    for beam in beams
        flap_deflection!(beam,x,timeIndex,1.0)
    end

    # Rates at the end of the time step (Newmark's relation): ẋₙ₊₁ = c xₙ₊₁ - (c xₙ + (1-γ)/γ ẋₙ), with c = 1/(γΔt)
    c = 1/(γ*Δt)
    b = (1-γ)/γ
    foreach_element(elementGroups) do element
        set_element_equivalent_states_rates!(element,c,c,b,xBegin,xdotBegin)
    end

    # Memory variables of the dynamic stall models at the begin of the time step
    for (element,memory) in zip(elements,dynamicStallMemoryBegin)
        restore_dynamic_stall_memory!(element,memory)
    end

    # Update the states of the elements, hinge axis constraints and special nodes at the end of the time step
    update_states!(problem)

end
