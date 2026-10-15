# Utility function to extract the primal value from ForwardDiff.Dual or Real numbers
@inline _val(x::Real) = x
@inline _val(x::ForwardDiff.Dual) = ForwardDiff.value(x)

# Computes the nodal resutants from aerodynamic loads and aerodynamic state matrices
function aero_core!(f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,A,B,model,element,aeroStateCache,p,V,Ω,χ,δ,Vdot,Ωdot,timeNow,Δt,isDynamicProblem)

    # Unpack and set flags
    @unpack atmosphere = model

    # Steady aerodynamic kinematics
    α = aero_steady_kinematics!(element,aeroStateCache,V,Ω)

    # Return zero arrays if the angle of attack is undefined
    if isnan(α)
        zero_arrays!(f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,A,B)
        return f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,B
    end
    
    # Unsteady aerodynamic kinematics
    aero_unsteady_kinematics!(element,aeroStateCache,Vdot,Ωdot)

    # Flow parameters
    flow_parameters!(atmosphere,element,aeroStateCache)

    # Update airfoil parameters
    update_airfoil_parameters!(element,aeroStateCache)

    # Update inertial parameters
    update_inertial_parameters!(element,aeroStateCache)

    # Update tip loss factor, if applicable
    isDynamicProblem && update_tip_loss_factor!(model,element)

    # Get local gust velocity, if applicable
    isDynamicProblem && local_gust_velocity!(timeNow,model,element,aeroStateCache,p)

    # Effective (unsteady) angle of attack
    effective_angle_of_attack!(element,aeroStateCache,χ,δ)

    # Aerodynamic coefficients 
    aero_coefficients!(isDynamicProblem,timeNow,Δt,element,aeroStateCache,χ,δ)

    # Aerodynamic state matrices 
    aero_state_matrices!(Ac,Ai,An,Af,Ag,A,B,element,aeroStateCache,δ)
    
    # Aerodynamic loads' nodal resultants
    aero_loads_resultants!(f1χ,f2χ,m1χ,m2χ,atmosphere,element,aeroStateCache,p)

    return f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,B

end

# Computes the steady aerodynamic kinematic variables
function aero_steady_kinematics!(element,aeroStateCache,V::AbstractVector{T},Ω::AbstractVector{T}) where T<:Number

    @unpack RwT,c,normSparPos = element.aero.geometry
    @unpack U = aeroStateCache.attachedVars.steadyKin

    # Relative wind velocity vector at spar position, resolved in basis W
    U .= - RwT * V

    # Spanwise, tangential and normal components of relative wind velocity
    Uₛ,Uₜ,Uₙ = U

    # In-plane relative wind speed
    Uᵢ = sqrt(Uₜ^2 + Uₙ^2)

    # Spanwise angular velocity component, resolved in basis W
    Ωₐ = dot(RwT[1,:], Ω)

    # Relative normal wind velocity component at the airfoil's 1/2- and 3/4-chord
    UₙMid = Uₙ - Ωₐ * c * (normSparPos - 1/2)
    UₙTQC = Uₙ - Ωₐ * c * (normSparPos - 3/4)

    # Angle of attack  
    α = !iszero(Uₜ) ? 
        atan(-Uₙ/Uₜ) : 
        T(NaN)

    # Sideslip angle 
    β = !iszero(Uᵢ) ? 
        atan(Uₛ/Uᵢ) :
        T(NaN)

    @pack! aeroStateCache.attachedVars.steadyKin = Uₛ,Uₜ,Uₙ,Uᵢ,UₙMid,UₙTQC,Ωₐ,α,β
    
    return α

end


# Computes the unsteady aerodynamic kinematic variables
function aero_unsteady_kinematics!(element,aeroStateCache,Vdot,Ωdot)

    @unpack RwT,c,normSparPos = element.aero.geometry
    @unpack Uₜ,Uₙ,Uᵢ = aeroStateCache.attachedVars.steadyKin
    @unpack Udot = aeroStateCache.attachedVars.unsteadyKin

    # Time derivative of relative wind velocity vector, resolved in basis W
    Udot .= - RwT * Vdot

    # Time derivatives of tangential and normal components of relative wind velocity
    Uₜdot,Uₙdot = Udot[2],Udot[3]

    # Time derivative of in-plane relative wind speed
    Uᵢdot = (Uₜ * Uₜdot + Uₙ * Uₙdot) / Uᵢ

    # Spanwise angular acceleration component, resolved in basis W
    Ωₐdot = dot(RwT[1,:], Ωdot)

    # Relative normal wind acceleration component at the airfoil's 1/2- and 3/4-chord
    UₙdotMid = Uₙdot - Ωₐdot * c * (normSparPos - 1/2)
    UₙdotTQC = Uₙdot - Ωₐdot * c * (normSparPos - 3/4)

    # Time derivative of angle of attack 
    αdot = (-Uₜ * Uₙdot + Uₙ * Uₜdot) / Uᵢ^2

    @pack! aeroStateCache.attachedVars.unsteadyKin = Uᵢdot,UₙdotMid,UₙdotTQC,Ωₐdot,αdot

end


# Computes flow parameters
function flow_parameters!(atmosphere,element,aeroStateCache)

    @unpack ρ,μ,a = atmosphere
    @unpack b,c = element.aero.geometry
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin

    # Reynolds
    Re = ρ*Uᵢ*c/μ

    # Mach
    Ma = Uᵢ/a

    # Prandtl-Glauert compressibility factor and its square
    βₚ² = Ma < 1 ? 
          1-Ma^2 :
          one(typeof(Ma))
    βₚ = sqrt(βₚ²)

    # Flow inverse time scale
    Θ = Uᵢ/b*βₚ²

    # Chord-to-sound-speed time scale
    Tᵢ = c/a

    @pack! aeroStateCache.attachedVars.flowParameters = Re,Ma,βₚ,βₚ²,Θ,Tᵢ

end


# Updates the aerodynamic parameters of the airfoil according to current nondimensional flow parameters 
function update_airfoil_parameters!(element,aeroStateCache)

    @unpack ppSolver = element.aero.surface.solvers
    @unpack updateAirfoilParameters,airfoil = element.aero.surface.geometry
    @unpack Ma = aeroStateCache.attachedVars.flowParameters

    # Skip if parameters are not to be updated
    updateAirfoilParameters || return

    # Update parameters based on the solver type
    if ppSolver isa AttachedSolver
        update_attached_flow_parameters!(airfoil,_val(Ma))
    elseif ppSolver isa BLi
        update_BLi_parameters!(airfoil,_val(Ma))
    elseif ppSolver isa BLo
        update_BLo_parameters!(airfoil,_val(Ma))
    end
    
end


# Updates the inertial indicial time scales
function update_inertial_parameters!(element,aeroStateCache)
    
    @unpack surface = element.aero
    @unpack ppSolver = surface.solvers

    # Skip if the solver is not unsteady attached or if incompressible inertial loads are enabled
    if !(ppSolver isa UnsteadySolver) || ppSolver.incompressibleInertialLoads
        return
    end

    @unpack AC,bC,AI,bI = ppSolver
    @unpack Ma,βₚ,Tᵢ = aeroStateCache.attachedVars.flowParameters

    # Set correction factor for circulatory aerodynamic states
    γbC = ppSolver isa DynamicStallSolver ?
          surface.geometry.γbC :
          ones(length(bC))

    # Repeatedly used variable
    ΣγACbC = sum(AC .* bC .* γbC)

    # Inertial indicial time scale parameters
    Tnα = Tᵢ / (1 - Ma + π * βₚ * Ma^2 * ΣγACbC)
    TnM = Tᵢ / (1 - Ma + π / βₚ * Ma^2 * ΣγACbC)
    Tnθ̇ = Tᵢ / (1 - Ma + 2π * βₚ * Ma^2 * ΣγACbC)
    Tmα = Tᵢ * (AI[1]*bI[2] + AI[2]*bI[1]) / (bI[1] * bI[2] * (1-Ma))
    TmM = Tmα
    Tmθ̇ = 7*Tᵢ / (15*(1-Ma) + 3π * βₚ * Ma^2 * AI[3]*bI[3])

    # Apply correction factor (see Jose et al. - Unsteady Aerodynamic Modeling with Time-Varying Free-Stream Mach Numbers - 2006)
    λ = 0.75
    Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = multiply_scalars!(λ,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇)
    
    @pack! aeroStateCache.attachedVars.flowParameters = Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇

end


# Updates the tip loss factor over the element with the current reference airspeed
function update_tip_loss_factor!(model,element)

    # Skip if the tip loss function is not airspeed dependent
    !(element.aero.surface.tipLoss.tipLossFunctionIsAirspeedDependent) && return

    @unpack v_A = model.kinematicsBasisA.state
    @unpack ϖ = element.aero.loadsDistribution

    # Reference airspeed
    Uref = norm(v_A)

    # Tip loss factor at element's midpoint value
    ϖMid = ϖ(Uref,0.5)

    @pack! element.aero.loadsDistribution = Uref,ϖMid

end


# Computes the current gust velocity
function local_gust_velocity!(timeNow,model,element,aeroStateCache,p)

    @unpack gust = model

    # Skip if gust is not defined
    isnothing(gust) && return

    @unpack RwR0 = element.aero.geometry
    @unpack U,UGust = aeroStateCache.attachedVars.steadyKin

    # Rotation tensor
    R = rotation_tensor_WM(p)

    # Reset gust velocity vector
    UGust .= zero(U)

    # Gust defined over time
    if gust.isDefinedOverTime && (gust.initialTime <= timeNow <= gust.finalTime)
        @unpack initialTime,finalTime = gust
        @unpack UGustInertial = gust
        @unpack R_AT = model.kinematicsBasisA.state
        # Transform gust velocity vector from inertial basis to local, deformed aerodynamic basis W
        UGust = (R*RwR0)' * R_AT * UGustInertial(timeNow)
    # Gust defined over space and time
    elseif !gust.isDefinedOverTime
        @unpack UGustInertial = gust
        @unpack u_A,R_A,R_AT = model.kinematicsBasisA.state
        @unpack r = element.assemblyData
        @unpack u = element.state.midpoint
        # Get current position of element's midpoint, resolved in the inertial basis
        r⁰ = u_A + R_A * (r + u)
        # Transform gust velocity vector from inertial basis to local, deformed aerodynamic basis W
        UGust = (R*RwR0)' * R_AT * UGustInertial(r⁰,timeNow)
    end

    # Set tangential and normal gust velocity components
    UₜGust,UₙGust = UGust[2],UGust[3]

    @pack! aeroStateCache.attachedVars.steadyKin = UₜGust,UₙGust

end


# Computes the effective (unsteady) angle of attack
function effective_angle_of_attack!(element,aeroStateCache,χ,δ)

    @unpack ppSolver = element.aero.surface.solvers
    @unpack airfoil = element.aero.surface.geometry
    @unpack Uₜ,UₜGust = aeroStateCache.attachedVars.steadyKin
    ϵₙ = get_airfoil_state_field(airfoil, ppSolver, :ϵₙ)

    # Quasi-steady flap-induced normalwash
    wf = flap_normalwash(element,aeroStateCache,δ)

    # Effective normalwash: pitch-plunge-induced, flap-induced, gust-induced, and total
    wₑp = pitch_plunge_effective_normalwash(element,aeroStateCache,χ,ppSolver)
    wₑf = flap_effective_normalwash(element,χ,wf)
    wₑg = gust_effective_normalwash(element,aeroStateCache,χ)
    wₑ = wₑp + ϵₙ*wₑf + wₑg

    # Effective circulatory AoA
    αₑ = atan( - wₑ / (Uₜ + UₜGust) )

    @pack! aeroStateCache.attachedVars.steadyKin = wf
    @pack! aeroStateCache.attachedVars.unsteadyKin = αₑ

end


# Computes the aerodynamic coefficients
function aero_coefficients!(isDynamicProblem,timeNow,Δt,element,aeroStateCache,χ,δ)

    @unpack ppSolver = element.aero.surface.solvers

    # Compute coefficients according to the solver type
    if ppSolver isa AttachedSolver
        attached_flow_aero_coefficients!(element,aeroStateCache,χ,δ)
    elseif ppSolver isa BLi
        BLi_aero_coefficients!(isDynamicProblem,timeNow,Δt,element,aeroStateCache,χ,δ)
    elseif ppSolver isa BLo
        BLo_aero_coefficients!(isDynamicProblem,timeNow,Δt,element,aeroStateCache,χ,δ)
    end

end


# Computes the aerodynamic state matrices
function aero_state_matrices!(Ac,Ai,An,Af,Ag,A,B,element,aeroStateCache,δ)

    @unpack ppSolver = element.aero.surface.solvers

    # Pitch-plunge state matrices
    if ppSolver isa Indicial
        indicial_state_matrices!(Ac,Ai,A,B,element,aeroStateCache)
    elseif ppSolver isa Inflow
        inflow_state_matrices!(Ac,Ai,A,B,element,aeroStateCache)
    elseif ppSolver isa BLi
        BLi_state_matrices!(Ac,Ai,An,A,B,element,aeroStateCache)
    elseif ppSolver isa BLo
        BLo_state_matrices!(Ac,Ai,An,A,B,element,aeroStateCache)
    end
    
    # Flap state matrices
    flap_state_matrices!(Af,A,B,element,aeroStateCache,δ)

    # Gust state matrices
    gust_state_matrices!(Ag,A,B,element,aeroStateCache)

end


# Computes the pitch-plunge-induced aerodynamic state matrices for the indicial ppSolver
function indicial_state_matrices!(Ac,Ai,A,B,element,aeroStateCache)

    @unpack circulatoryPitchPlungeStatesRange,inertialPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack incompressibleInertialLoads,AC,bC,bI = element.aero.surface.solvers.ppSolver
    @unpack Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack UₙdotTQC = aeroStateCache.attachedVars.unsteadyKin
    @unpack Ma,Θ,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = aeroStateCache.attachedVars.flowParameters

    # Rate of cnα * UₙTQC
    cnαUₙTQCdot = cnαUₙTQC_rate(element,aeroStateCache)

    # Circulatory state matrix diagonal
    Ac .= - Θ * bC

    # Inertial state matrix diagonal
    Ai .= incompressibleInertialLoads ? 
          z0 : 
          - SVector(1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM))

    # State matrix
    A[circulatoryPitchPlungeStatesRange,circulatoryPitchPlungeStatesRange] .= Diagonal(Ac)
    A[inertialPitchPlungeStatesRange,inertialPitchPlungeStatesRange] .= Diagonal(Ai)

    # Input matrix
    B[circulatoryPitchPlungeStatesRange] .= cnαUₙTQCdot * AC
    B[inertialPitchPlungeStatesRange] .= incompressibleInertialLoads ? 
                                         z0 : 
                                         SVector(α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)

end


# Computes the pitch-plunge-induced aerodynamic state matrices for the inflow ppSolver
function inflow_state_matrices!(Ac,Ai,A,B,element,aeroStateCache)

    @unpack circulatoryPitchPlungeStatesRange,inertialPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack incompressibleInertialLoads,AₚInv,AₚInvVec,AₚInvcₚ,bI,nInflowStates = element.aero.surface.solvers.ppSolver
    @unpack Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack UₙdotTQC = aeroStateCache.attachedVars.unsteadyKin
    @unpack Ma,Θ,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = aeroStateCache.attachedVars.flowParameters
    
    # Circulatory state matrix (flattened)
    Ac .= - Θ * AₚInvVec

    # Inertial state matrix diagonal
    Ai .= incompressibleInertialLoads ? 
          z0 : 
          - SVector(1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM))

    # State matrix
    A[circulatoryPitchPlungeStatesRange,circulatoryPitchPlungeStatesRange] .= -Θ * AₚInv
    A[inertialPitchPlungeStatesRange,inertialPitchPlungeStatesRange] .= Diagonal(Ai)

    # Input matrix
    B[circulatoryPitchPlungeStatesRange] .= UₙdotTQC * AₚInvcₚ
    B[inertialPitchPlungeStatesRange] .= incompressibleInertialLoads ? 
                                         z0 : 
                                         SVector(α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)

end


# Computes the aerodynamic state matrices for the modified Beddoes-Leishman model
function BLi_state_matrices!(Ac,Ai,An,A,B,element,aeroStateCache)

    @unpack circulatoryPitchPlungeStatesRange,inertialPitchPlungeStatesRange,nonlinearStatesRange = element.aero.surface.stateRanges
    @unpack γbC = element.aero.surface.geometry
    @unpack incompressibleInertialLoads,AC,bC,bI = element.aero.surface.solvers.ppSolver
    @unpack Ta,Tg = element.aero.surface.geometry.airfoil.state.bli
    @unpack isFlapped,b = element.aero.geometry
    @unpack Uₙ,Uᵢ,UₙGust,Ωₐ,α,wf = aeroStateCache.attachedVars.steadyKin
    @unpack UₙdotTQC = aeroStateCache.attachedVars.unsteadyKin
    @unpack βₚ²,Θ,Ma,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = aeroStateCache.attachedVars.flowParameters
    @unpack R,Ta_SO,TfN,TfM,TfT,fPrimeN,fPrimeM,fPrimeT = aeroStateCache.bliVars.flow
        
    # Rate of cnα * UₙTQC
    cnαUₙTQCdot = cnαUₙTQC_rate(element,aeroStateCache)

    # Circulatory state matrix diagonal
    Ac .= - Θ * bC .* γbC

    # Inertial state matrix diagonal
    Ai .= incompressibleInertialLoads ? 
          z0 :
          - SVector(1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM))

    # Nonlinear state matrix diagonal
    An .= - Uᵢ/b * SVector(1/Ta, 1/TfN, 1/TfM, 1/TfT, 1/(3Ta), 1/Ta_SO, 1/Tg)          

    # State matrix
    A[circulatoryPitchPlungeStatesRange,circulatoryPitchPlungeStatesRange] .= Diagonal(Ac)
    A[inertialPitchPlungeStatesRange,inertialPitchPlungeStatesRange] .= Diagonal(Ai)
    A[nonlinearStatesRange,nonlinearStatesRange] .= Diagonal(An)

    # Input matrix
    B[circulatoryPitchPlungeStatesRange] .= cnαUₙTQCdot * AC
    B[inertialPitchPlungeStatesRange] = incompressibleInertialLoads ? 
                                        z0 :
                                        SVector(α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)
    B[nonlinearStatesRange] .= Uᵢ/b * SVector((Uₙ+wf)/Ta, fPrimeN/TfN, fPrimeM/TfM, fPrimeT/TfT, R/(3Ta), R/Ta_SO, UₙGust/Tg)

end


# Computes the aerodynamic state matrices for the original Beddoes-Leishman model
function BLo_state_matrices!(Ac,Ai,An,A,B,element,aeroStateCache)

    @unpack b = element.aero.geometry
    @unpack circulatoryPitchPlungeStatesRange,inertialPitchPlungeStatesRange,nonlinearStatesRange = element.aero.surface.stateRanges
    @unpack γbC = element.aero.surface.geometry
    @unpack incompressibleInertialLoads,AC,bC,bI = element.aero.surface.solvers.ppSolver
    @unpack Tf₀,Tp = element.aero.surface.geometry.airfoil.state.blo
    @unpack Uᵢ,UₙTQC,Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack Uᵢdot = aeroStateCache.attachedVars.unsteadyKin
    @unpack βₚ²,Θ,Ma,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = aeroStateCache.attachedVars.flowParameters
    @unpack cnP = aeroStateCache.attachedVars.aeroCoefficients
    @unpack Tf,Tv,f,fPrime,cvdotP,cvdotN = aeroStateCache.bloVars.flow

    # Rate of cnα * UₙTQC
    cnαUₙTQCdot = cnαUₙTQC_rate(element,aeroStateCache)

    # Circulatory state matrix diagonal
    Ac .= - Θ * bC .* γbC

    # Inertial state matrix diagonal
    Ai .= incompressibleInertialLoads ? 
          z0 :
          - SVector(1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM))

    # Nonlinear state matrix diagonal
    An .= - Uᵢ/b * SVector(1/Tp, 1/Tf, 2/Tf₀, 1/Tv, 1/Tv)

    # State matrix
    A[circulatoryPitchPlungeStatesRange,circulatoryPitchPlungeStatesRange] .= Diagonal(Ac)
    A[inertialPitchPlungeStatesRange,inertialPitchPlungeStatesRange] .= Diagonal(Ai)
    A[nonlinearStatesRange, nonlinearStatesRange] .= Diagonal(An)

    # Input matrix
    B[circulatoryPitchPlungeStatesRange] .= cnαUₙTQCdot * AC
    B[inertialPitchPlungeStatesRange] = incompressibleInertialLoads ? 
                                        z0 :
                                        SVector(α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)
    B[nonlinearStatesRange] .= SVector(cnP/(Tp*b/Uᵢ), fPrime/(Tf*b/Uᵢ), f/(Tf₀*b/Uᵢ/2), cvdotP, cvdotN)

end


# Computes the aerodynamic state matrices for the flap states
function flap_state_matrices!(Af,A,B,element,aeroStateCache,δ)

    @unpack flapStatesRange = element.aero.stateRanges

    # Skip if there are no flap-induced flow states
    isempty(flapStatesRange) && return nothing

    # Unpack
    @unpack ACf,bCf = element.aero.surface.solvers.flapSolver
    @unpack Θ = aeroStateCache.attachedVars.flowParameters

    # Rate of cnα * wf
    cnαwfdot = cnαwf_rate(element,aeroStateCache,δ)

    # State matrix diagonal
    Af .= - Θ * bCf

    # State matrix
    A[flapStatesRange,flapStatesRange] .= Diagonal(Af)

    # Input matrix
    B[flapStatesRange] .= cnαwfdot * ACf

    return nothing

end


# Computes the aerodynamic state matrices for the gust states
function gust_state_matrices!(Ag,A,B,element,aeroStateCache)

    @unpack gustStatesRange = element.aero.stateRanges

    # Skip if there are no gust-induced flow states
    isempty(gustStatesRange) && return nothing

    # Unpack
    @unpack bG = element.aero.surface.solvers.gustSolver
    @unpack UₙGust = aeroStateCache.attachedVars.steadyKin
    @unpack Θ = aeroStateCache.attachedVars.flowParameters

    # State matrix diagonal
    Ag .= - Θ * bG

    # State matrix
    A[gustStatesRange,gustStatesRange] .= Diagonal(Ag)

    # Input matrix
    B[gustStatesRange] .= UₙGust

    return nothing

end


# Computes the aerodynamic coefficients under attached flow conditions
function attached_flow_aero_coefficients!(element,aeroStateCache,χ,δ)

    # Normal force coefficient
    attached_flow_cn!(element,aeroStateCache,χ,δ)

    # Pitching moment coefficient about the spar position
    attached_flow_cm!(element,aeroStateCache,χ,δ)

    # Tangential force coefficient
    attached_flow_ct!(element,aeroStateCache,δ)

    # Scale aerodynamic coefficients with the loads distribution correction factor
    scale_attached_flow_aero_coefficients!(element,aeroStateCache)

end


# Computes the normal force aerodynamic coefficient for attached flow
function attached_flow_cn!(element,aeroStateCache,χ,δ)

    @unpack ppSolver,flapSolver = element.aero.surface.solvers
    @unpack smallAngles = element.aero.surface.geometry
    @unpack inertialPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack incompressibleInertialLoads = ppSolver
    @unpack ϵₙ,cnα = element.aero.surface.geometry.airfoil.state.attached
    @unpack cnδ = element.aero.surface.geometry.flapParameters
    @unpack isFlapped,Th,b,aₕ = element.aero.geometry
    @unpack Uᵢ,Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack αₑ,UₙdotMid = aeroStateCache.attachedVars.unsteadyKin
    @unpack Ma,Tnα,TnM,Tnθ̇ = aeroStateCache.attachedVars.flowParameters
    if ppSolver isa UnsteadySolver && !incompressibleInertialLoads
        χ1,χ2,χ3,χ4,χ5,χ6,χ7,χ8 = view(χ,inertialPitchPlungeStatesRange)
    end

    # Circulatory component
    cnC = smallAngles ? 
          cnα * αₑ :
          cnα * sin(αₑ) * cos(αₑ)
    if isFlapped && flapSolver isa TableLookup
        cnC += cnδ * δ
    end

    # Inertial component
    cnI = incompressibleInertialLoads ? 
          π*b*UₙdotMid/Uᵢ^2 : 
          4/Ma * (α - χ1/Tnα) + -4*aₕ*b/(Uᵢ*Ma) * (Ωₐ - χ2/Tnθ̇) + 4*α/Ma^2 * (Ma - χ6/TnM)
    if isFlapped && flapSolver isa ThinAirfoilTheory
        @unpack δdot,δddot = element.aero.surface.state
        cnI -= ϵₙ*b/Uᵢ^2 * (Uᵢ*Th[4]*δdot + b*Th[1]*δddot)
    end

    # Non-circulatory
    cnNC = cnI

    # Total
    cn = cnC + cnNC

    @pack! aeroStateCache.attachedVars.aeroCoefficients = cn,cnC,cnI,cnNC

end


# Computes the pitching moment aerodynamic coefficient at the attachment point (i.e., the beam reference line) for attached flow
function attached_flow_cm!(element,aeroStateCache,χ,δ)

    @unpack ppSolver,flapSolver = element.aero.surface.solvers
    @unpack inertialPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack incompressibleInertialLoads = ppSolver
    @unpack ϵₘ,cm₀,cnα,x_ac = element.aero.surface.geometry.airfoil.state.attached
    @unpack cmδ = element.aero.surface.geometry.flapParameters
    @unpack isFlapped,Th,b,aₕ,normSparPos = element.aero.geometry
    @unpack Uᵢ,Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack UₙdotMid,Ωₐdot = aeroStateCache.attachedVars.unsteadyKin
    @unpack Ma,βₚ,Tmα,TmM,Tmθ̇ = aeroStateCache.attachedVars.flowParameters
    @unpack cnC = aeroStateCache.attachedVars.aeroCoefficients
    if ppSolver isa UnsteadySolver && !incompressibleInertialLoads
        @unpack AI,bI = ppSolver
        χ1,χ2,χ3,χ4,χ5,χ6,χ7,χ8 = view(χ,inertialPitchPlungeStatesRange)
    end

    # Circulatory component
    cmC = cnC * (normSparPos - x_ac)
    if isFlapped && flapSolver isa TableLookup
        cmC += cmδ * δ
    end

    # Inertial component
    cmI = incompressibleInertialLoads ? 
          -π*b/(2*Uᵢ^2) * (-aₕ * UₙdotMid + b * Ωₐdot / 8 + Uᵢ * Ωₐ / 2) : 
          2*aₕ/Ma * (α - χ3*AI[1]/(bI[1]*Tmα) - χ4*AI[2]/(bI[2]*Tmα)) + -2*b*(aₕ^2+1/3)/(Uᵢ*Ma) * (Ωₐ - χ5/Tmθ̇) + 2*aₕ*α/Ma^2 * (Ma - χ7*AI[1]/(bI[1]*TmM) - χ8*AI[2]/(bI[2]*TmM)) - π*b*Ωₐ/(4*Uᵢ*βₚ)
    if isFlapped && flapSolver isa ThinAirfoilTheory
        @unpack δdot,δddot = element.aero.surface.state
        cmI -= ϵₘ/(2*Uᵢ^2) * (Uᵢ^2*Th[14]*δ + Uᵢ*b*Th[15]*δdot + b^2*Th[16]*δddot)
    end

    # Non-circulatory
    cmNC = cmI + cm₀

    # Total at attachment point
    cm = cmC + cmNC

    @pack! aeroStateCache.attachedVars.aeroCoefficients = cm,cmC,cmI,cmNC

end


# Computes the tangential force aerodynamic coefficient for attached flow
function attached_flow_ct!(element,aeroStateCache,δ)

    @unpack flapSolver = element.aero.surface.solvers
    @unpack smallAngles,hasInducedDrag = element.aero.surface.geometry
    @unpack η,cd₀ = element.aero.surface.geometry.airfoil.state.attached
    @unpack cdδ = element.aero.surface.geometry.flapParameters
    @unpack isFlapped,b = element.aero.geometry
    @unpack Uᵢ,Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack αₑ = aeroStateCache.attachedVars.unsteadyKin
    @unpack cnC = aeroStateCache.attachedVars.aeroCoefficients

    # Parasite component
    ct₀ = smallAngles ? 
          -cd₀*one(typeof(α)) : 
          -cd₀/cos(α)

    # Circulatory component
    ctC = smallAngles ? 
          η * cnC * (αₑ+b*Ωₐ/(2*Uᵢ)) : 
          η * cnC * tan(αₑ+b*Ωₐ/(2*Uᵢ))
    if isFlapped && flapSolver isa TableLookup
        ctC -= cdδ * abs(δ) / cos(α)
    end

    # Induced (drag) component
    if hasInducedDrag
        @unpack AR = element.aero.surface.assemblyData
        ctC -= cnC^2 / (π*AR) / cos(α)
    end

    # Non-circulatory
    ctNC = ct₀

    # Total
    ct = ctC + ctNC

    @pack! aeroStateCache.attachedVars.aeroCoefficients = ct,ctC,ctNC

end


# Scales the circulatory components of aerodynamic coefficients and recompute the totals for the attached flow solvers
function scale_attached_flow_aero_coefficients!(element,aeroStateCache)

    # Unpack
    @unpack ϖMid = element.aero.loadsDistribution
    @unpack cnC,cmC,ctC,cnNC,cmNC,ctNC = aeroStateCache.attachedVars.aeroCoefficients

    # Scale circulatory components
    cnC,cmC,ctC = multiply_scalars!(ϖMid, cnC,cmC,ctC)

    # Recompute totals
    cn = cnC + cnNC
    cm = cmC + cmNC
    ct = ctC + ctNC

    # Set the scaled aerodynamic coefficients
    @pack! aeroStateCache.attachedVars.aeroCoefficients = cn,cm,ct,cnC,cmC,ctC

end


# Computes the aerodynamic coefficients according to the modified incompressible flow Beddoes-Leishman model
function BLi_aero_coefficients!(isDynamicProblem,timeNow,Δt,element,aeroStateCache,χ,δ)

    # Name nonlinear states
    BLi_nonlinear_states!(element,aeroStateCache,χ)

    # Additional kinematics
    BLi_kinematics!(element,aeroStateCache)

    # Motion qualifiers
    BLi_motion_qualifiers!(timeNow,element,aeroStateCache)

    # Breakpoint angles
    BLi_breakpoint_angles!(element,aeroStateCache)

    # Time delay constants
    BLi_time_delays!(element,aeroStateCache)

    # Flow separation points
    BLi_separation_points!(element,aeroStateCache)

    # Variables at stall onset and DSV loads
    if isDynamicProblem
        BLi_stall_time!(timeNow,Δt,element,aeroStateCache)
        BLi_DSV_loads!(element,aeroStateCache)
    else
        cnV = cmV = ctV = 0.0
        @pack! aeroStateCache.attachedVars.aeroCoefficients = cnV,cmV,ctV
    end

    # Normal force coefficient
    BLi_cn!(element,aeroStateCache,χ,δ)

    # Pitching moment coefficient about the spar position
    BLi_cm!(element,aeroStateCache,χ,δ)

    # Tangential flow coefficient
    BLi_ct!(element,aeroStateCache,δ)

    # Scale circulatory components of aerodynamic coefficients and recompute totals
    scale_BLi_aero_coefficients!(element,aeroStateCache)

end


# Sets the nonlinear states of the Beddoes-Leishman model
function BLi_nonlinear_states!(element,aeroStateCache,χ::AbstractVector{T}) where T<:Number

    @unpack nonlinearStatesRange = element.aero.surface.stateRanges
    @unpack f₀N,f₀M,f₀T = element.aero.surface.geometry.airfoil.state.bli
    @unpack Uₜ,UₜGust = aeroStateCache.attachedVars.steadyKin

    # Zero of current type
    z = zero(T)

    # Aliases for nonlinear states
    χnl1,χnl2,χnl3,χnl4,χnl5,χnl6,χnl7 = view(χ, nonlinearStatesRange)

    # Set nonlinear states
    αlag = !iszero(Uₜ + UₜGust) ? 
           atan( -(χnl1 + χnl7) / (Uₜ + UₜGust) ) :
           z
    f2primeN = clamp(χnl2, f₀N, 1)
    f2primeM = clamp(χnl3, f₀M, 1)
    f2primeT = clamp(χnl4, f₀T, 1)
    RD = max(z, χnl5)
    RD_stallOnsetRatio = max(z, χnl6)

    @pack! aeroStateCache.bliVars.states = αlag,f2primeN,f2primeM,f2primeT,RD,RD_stallOnsetRatio

end


# Computes additional kinematics associated with the Beddoes-Leishman model
function BLi_kinematics!(element,aeroStateCache)

    @unpack c = element.aero.geometry
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin
    @unpack αdot = aeroStateCache.attachedVars.unsteadyKin
    @unpack r₀ = element.aero.surface.geometry.airfoil.state.bli

    # Non-dimensional pitch rate (round off to suppress noise)
    q = round_off(αdot*c/Uᵢ, 1e-6)

    # Reduced non-dimensional pitch rate
    r = q/2

    # Unsigned ratio of reduced pitch rate to critical pitch rate
    qR = abs(r)/r₀

    # Unsigned capped reduced pitch rate ratio 
    R = min(1, qR)

    @pack! aeroStateCache.bliVars.flow = q,r,qR,R

end


# Computes motion qualifiers of the modified Beddoes-Leishman model
function BLi_motion_qualifiers!(timeNow,element,aeroStateCache)

    @unpack b = element.aero.geometry
    @unpack αds₀,αₛₛ,γLS,Tf = element.aero.surface.geometry.airfoil.state.bli
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin
    @unpack αlag,RD,RD_stallOnsetRatio = aeroStateCache.bliVars.states
    @unpack r,qR = aeroStateCache.bliVars.flow
    @unpack αlagPrev,stallOnsetRatioPrev,PPrev,qRPrev,upstrokePrev,stallOnsetRatioMax,stallOnsetRatioMin,qRmax,tv0P,tv0N,RD_tv0P,RD_tv0N = aeroStateCache.bliVars.comp

    # Critical angle for dynamic stall
    αcrit = αₛₛ + (αds₀ - αₛₛ) * RD_stallOnsetRatio

    # Dynamic stall onset ratio
    stallOnsetRatio = αlag / αcrit

    # Flag for upstroke phase
    upstroke = stallOnsetRatio * r >= 0
    
    # Flag for upstroke phase beginning at current time step
    upstrokeBeginning = upstroke && !upstrokePrev
    
    # Flag for stalled condition on current cycle
    stalled = abs(stallOnsetRatioMax) > 1

    # Light stall factor
    P = stalled ? 
        exp(-γLS*(abs(stallOnsetRatioMax)^4-1)) :
        PPrev 

    # Steady condition factor
    SS = exp(-50*RD)

    # Update maximum stall onset ratio and corresponding delayed capped reduced pitch rate (update only for increments in the lagged AoA)
    if abs(stallOnsetRatio) > abs(stallOnsetRatioPrev) && abs(αlag) > abs(αlagPrev) 
        stallOnsetRatioMax = stallOnsetRatio
    end

    # Update minimum stall onset ratio (approximate as that at begin of upstroke)
    if upstrokeBeginning
        stallOnsetRatioMin = stallOnsetRatio
    end

    # Update maximum reduced pitch rate ratio
    if qR > qRPrev
        qRmax = qR
    end

    # Time decay function for tangential force
    tv0 = max(tv0P, tv0N)
    lastRD_tv0 = tv0P > tv0N ? 
                 RD_tv0P : 
                 RD_tv0N
    Ts = isfinite(tv0) ? 
         exp( -5*(timeNow-tv0) * max(0.1,lastRD_tv0)^2 / (Tf*b/Uᵢ)) : 
         0.0      

    @pack! aeroStateCache.bliVars.flow = stallOnsetRatio,upstroke,stalled,P,SS
    @pack! aeroStateCache.bliVars.comp = Ts,lastRD_tv0,stallOnsetRatioMax,stallOnsetRatioMin,qRmax

end


# Computes the breakpoint angles of the modified Beddoes-Leishman model
function BLi_breakpoint_angles!(element,aeroStateCache)

    @unpack αds₀,αₛₛ,α1₀N,α1₀M,α1₀T,δα₀,δα₁,dt,dm,zm,ztd = element.aero.surface.geometry.airfoil.state.bli
    @unpack RD = aeroStateCache.bliVars.states
    @unpack qR,R,upstroke,stalled,P = aeroStateCache.bliVars.flow
    @unpack qRmax = aeroStateCache.bliVars.comp

    # Breakpoint angle offsets
    sqrtOfR = !iszero(R) ? 
              sqrt(R) : 
              zero(R)
    if upstroke
        fR = !iszero(qRmax) ? 
             R / min(1, qRmax) * sqrtOfR :
             zero(R)
        δα1N = (αds₀ - αₛₛ) * RD        
        δα1M = δα1N + dm*fR                                  
        δα1T = δα1N + dt*fR
    else 
        δα1N = stalled ? 
               -min(α1₀N/2, δα₀ + δα₁ * qR * (1+P)) :
               zero(R)
        δα1M = δα1N * (1 - zm * sqrtOfR * (1-P/2))
        δα1T = δα1N * (1 - ztd * sqrtOfR * (1-P/4))
    end

    # Unsteady breakpoint of separation angles
    α1N = α1₀N + δα1N 
    α1M = α1₀M + δα1M
    α1T = α1₀T + δα1T

    @pack! aeroStateCache.bliVars.flow = α1N,α1M,α1T

end


# Computes time delay variables for the modified Beddoes-Leishman model
function BLi_time_delays!(element,aeroStateCache)

    @unpack Ta,Tf,λ₁,λ₂ = element.aero.surface.geometry.airfoil.state.bli
    @unpack RD = aeroStateCache.bliVars.states
    @unpack qR,R,stallOnsetRatio,upstroke,P = aeroStateCache.bliVars.flow
    @unpack qRmax = aeroStateCache.bliVars.comp

    # Stall suppression factor (tested on the primal value only: iszero is false for a ForwardDiff.Dual with nonzero partials, and the derivative of RD^λ₂ with λ₂ < 1 is unbounded at RD = 0)
    σ = iszero(_val(RD)) ? 
        zero(RD) : 
        RD^λ₂

    # Stall supression time delay
    Ta_SO = Ta * (1/4 + 
                  λ₁ * σ * (1-R^(1/4)) * exp( -(abs(stallOnsetRatio) - 1)^2 / 0.01 )
                  )

    TfN = Tf * (1 +
                1/2 * P^6 * (!upstroke) + 
                RD^3 * (4 * (1-qR/qRmax) * (!upstroke) + 
                        1/2 * (qR/qRmax)^4 * upstroke
                        )
                )
    TfM = Tf * (1 - R/2 * (1-P) * (!upstroke))
    TfT = TfN

    @pack! aeroStateCache.bliVars.flow = Ta_SO,TfN,TfM,TfT

end


# Computes the flow separation points of the modified Beddoes-Leishman model
function BLi_separation_points!(element,aeroStateCache)

    @unpack α1₀N,α1₀M,α1₀T,βσ1N,βσ1T,βσ2N,βS2Nlpr,βS2Tlpr,βS1Nu,βS1Mu,βS1Tu,βS1Nd,βS1Md,βS1Td,βS2Nu,βS2Mu,βS2Tu,βS2Nd,βS2Md,βS2Td,ξ,f₀N,f₀M,f₀T,fbN,fbM,fbT,S1N,S1M,S1T,S2N,S2M,S2T = element.aero.surface.geometry.airfoil.state.bli
    @unpack α = aeroStateCache.attachedVars.steadyKin
    @unpack αlag,RD = aeroStateCache.bliVars.states
    @unpack R,stallOnsetRatio,upstroke,stalled,P,SS,α1N,α1M,α1T = aeroStateCache.bliVars.flow
    @unpack stallOnsetRatioMin,lastRD_tv0 = aeroStateCache.bliVars.comp

    # Factor defining the flattening of fPrime for light stall
    ψ = P^3 * (1 - SS)

    # Factor defining how close to stalling angle has been at beginning of upstroke
    σ = stallOnsetRatioMin*stallOnsetRatio > 0 ? 
        min(0.5, abs(stallOnsetRatioMin)^2.5) : 
        0.0 

    # In-and-out of stall factor
    ξRD_tv0 = ξ * ((1 - cos(2π*lastRD_tv0))/2)

    # Quasi-steady separation points: fN, fM, fT
    absα = abs(α)
    fN = absα <= α1₀N ? 
            1 - (1-fbN) * exp((absα-α1₀N)/S1N) : 
            f₀N + (fbN-f₀N) * exp((α1₀N-absα)/S2N)
    fM = absα <= α1₀M ? 
            1 - (1-fbM) * exp((absα-α1₀M)/S1M) : 
            f₀M + (fbM-f₀M) * exp((α1₀M-absα)/S2M)
    fT = absα <= α1₀T ? 
            1 - (1-fbT) * exp((absα-α1₀T)/S1T) : 
            f₀T + (fbT-f₀T) * exp((α1₀T-absα)/S2T)

    # Normal force unsteady lagged separation point based on lagged AoA: fPrimeN
    absαlag = abs(αlag)
    sqrtOfR = sqrt(R)
    if absαlag <= α1N
        αDiffN = absαlag - α1N
        S1PrimeN = upstroke ?
                   S1N * (1 + βS1Nu * RD) : 
                   S1N * (1 + βS1Nd * RD)
        fPrimeN = 1 - (1-fbN) * exp(αDiffN/S1PrimeN)
    else
        αDiffN = α1N - absαlag
        if upstroke
            f₀Nprime = f₀N
            S2PrimeN = S2N * (1 + βS2Nu*RD + βS2Nlpr*(1-RD)^2*(1-SS) + βσ2N*σ*RD*(1-RD))
        else                
            f₀Nprime = min(0.25, (1-sqrtOfR)*((stalled ? f₀N : fbN-0.2)+ξRD_tv0*σ))
            S2PrimeN = S2N * (1 + βS2Nd*RD + βS2Nlpr*(1-RD)^2*(1-SS) + βσ2N*σ*RD*(1-RD) + βσ1N*ψ*R)
        end
        fPrimeN = f₀Nprime + (fbN-f₀Nprime) * exp(αDiffN/S2PrimeN)
    end

    # Pitching moment unsteady lagged separation point based on lagged AoA: fPrimeM
    if absαlag <= α1M
        αDiffM = absαlag - α1M
        S1PrimeM = upstroke ? 
                   S1M * (1 + βS1Mu * RD) : 
                   S1M * (1 + βS1Md * RD)
        fPrimeM = 1-(1-fbM)*exp(αDiffM/S1PrimeM)
    else
        αDiffM = α1M - absαlag
        if upstroke
            f₀Mprime = f₀M
            S2PrimeM = S2M * (1 + βS2Mu * RD)
        else                
            f₀Mprime = min(0.25, (1-sqrtOfR)*((stalled ? f₀M : fbM-0.2)+ξRD_tv0*σ))
            S2PrimeM = S2M * (1 + βS2Md * RD)
        end
        fPrimeM = f₀Mprime + (fbM-f₀Mprime) * exp(αDiffM/S2PrimeM)
    end

    # Tangential force unsteady lagged separation point based on lagged AoA: fPrimeT
    if absαlag <= α1T
        αDiffT = absαlag - α1T
        S1PrimeT = upstroke ? 
                   S1T * (1 + βS1Tu * RD) : 
                   S1T * (1 + βS1Td * RD)
        fPrimeT = 1-(1-fbT)*exp(αDiffT/S1PrimeT)
    else
        αDiffT = α1T - absαlag
        if upstroke
            f₀Tprime = f₀T
            S2PrimeT = S2T * (1 + βS2Tu*RD + βS2Tlpr*(1-RD)^2*(1-SS))
        else
            f₀Tprime = min(0.25, (1-sqrtOfR)*((stalled ? f₀T : fbT-0.2)+ξRD_tv0*σ))
            S2PrimeT = S2T * (1 + βS2Td*RD + βS2Tlpr*(1-RD)^2*(1-SS) + βσ1T*ψ*R)
        end
        fPrimeT = f₀Tprime + (fbT-f₀Tprime) * exp(αDiffT/S2PrimeT)
    end

    @pack! aeroStateCache.bliVars.flow = fN,fM,fT,fPrimeN,fPrimeM,fPrimeT

end


# Computes the quasi-steady flow separation points of the modified Beddoes-Leishman model for an element
function BLi_quasi_steady_separation_points(element,Uₜ,Uₙ)

    # Check tangential velocity
    if isapprox(Uₜ, 0)
        fN,fM,fT = 1.0,1.0,1.0
        return fN,fM,fT
    end

    @unpack α1₀N,α1₀M,α1₀T,f₀N,f₀M,f₀T,fbN,fbM,fbT,S1N,S1M,S1T,S2N,S2M,S2T = element.aero.surface.geometry.airfoil.state.bli

    # Angle of attack  
    α = atan(-Uₙ/Uₜ)
    absα = abs(α)

    # Quasi-steady separation points
    fN = absα <= α1₀N ? 
         1 - (1-fbN) * exp((absα-α1₀N)/S1N) : 
         f₀N + (fbN-f₀N) * exp((α1₀N-absα)/S2N)

    fM = absα <= α1₀M ? 
         1 - (1-fbM) * exp((absα-α1₀M)/S1M) : 
         f₀M + (fbM-f₀M) * exp((α1₀M-absα)/S2M)

    fT = absα <= α1₀T ? 
         1 - (1-fbT) * exp((absα-α1₀T)/S1T) : 
         f₀T + (fbT-f₀T) * exp((α1₀T-absα)/S2T)

    return fN,fM,fT

end


# Computes variables at the time of stall onset for the modified Beddoes-Leishman model
function BLi_stall_time!(timeNow,Δt,element,aeroStateCache)

    @unpack b = element.aero.geometry
    @unpack gᵥ,Tv = element.aero.surface.geometry.airfoil.state.bli
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin
    @unpack f2primeN,RD = aeroStateCache.bliVars.states
    @unpack qR,R,stallOnsetRatio,upstroke,stalled,fN = aeroStateCache.bliVars.flow
    @unpack stallOnsetRatioPrev,fDiff_tv0P,qR_tv0P,R_tv0P,RD_tv0P,upstroke_tv0P,secondPVortexYetToCome,fDiff_tv0P2,RD_tv0P2,upstroke_tv0P2,fDiff_tv0N,qR_tv0N,R_tv0N,RD_tv0N,upstroke_tv0N,secondNVortexYetToCome,fDiff_tv0N2,RD_tv0N2,upstroke_tv0N2,tv0P,tv0N = aeroStateCache.bliVars.comp

    # Time instances
    t1 = timeNow - Δt
    t2 = timeNow

    # Positive AoA
    #---------------------------------------------------------------------------
    # Primary vortex
    if stallOnsetRatio >= 1 && stallOnsetRatioPrev < 1
        # Linear interpolation for time of vortex shedding 
        tv0P = t1 + (1 - stallOnsetRatioPrev) / (stallOnsetRatio - stallOnsetRatioPrev) * (t2 - t1)
        # Variables at time of vortex shedding
        fDiff_tv0P = max(0, f2primeN-fN)
        qR_tv0P = qR
        R_tv0P = R
        RD_tv0P = RD
        upstroke_tv0P = upstroke
        # Update flag
        secondPVortexYetToCome = true
    end
    # Nondimensional time since primary vortex shedding 
    τvP = (timeNow - tv0P) * Uᵢ / b
    # Secondary vortices
    if τvP >= (2+gᵥ)*Tv && secondPVortexYetToCome
        # Update flag
        secondPVortexYetToCome = false
        # Set variables at the time of second vortex
        fDiff_tv0P2 = max(0, f2primeN-fN)
        RD_tv0P2 = RD
        upstroke_tv0P2 = upstroke
    end

    # Negative AoA
    #---------------------------------------------------------------------------
    # Primary vortex
    if stallOnsetRatio <= -1 && stallOnsetRatioPrev > -1 
        # Linear interpolation for time of vortex shedding 
        tv0N = t1 + (1 - (-stallOnsetRatioPrev)) / ((-stallOnsetRatio) - (-stallOnsetRatioPrev)) * (t2 - t1)
        # Variables at time of vortex shedding
        fDiff_tv0N = max(0, f2primeN-fN)
        qR_tv0N = qR
        R_tv0N = R
        RD_tv0N = RD
        upstroke_tv0N = upstroke
        # Update flag
        secondNVortexYetToCome = true
    end
    # Nondimensional time since primary vortex shedding 
    τvN = (timeNow - tv0N) * Uᵢ / b
    # Secondary vortices
    if τvN >= (2+gᵥ)*Tv && secondNVortexYetToCome
        # Update flag
        secondNVortexYetToCome = false
        # Set variables at the time of second vortex
        fDiff_tv0N2 = max(0, f2primeN-fN)
        RD_tv0N2 = RD
        upstroke_tv0N2 = upstroke
    end

    @pack! aeroStateCache.bliVars.comp = tv0P,fDiff_tv0P,qR_tv0P,R_tv0P,RD_tv0P,upstroke_tv0P,secondPVortexYetToCome,τvP,fDiff_tv0P2,RD_tv0P2,upstroke_tv0P2,tv0N,fDiff_tv0N,qR_tv0N,R_tv0N,RD_tv0N,upstroke_tv0N,secondNVortexYetToCome,τvN,fDiff_tv0N2,RD_tv0N2,upstroke_tv0N2

end


# Computes DSV loads for the modified Beddoes-Leishman model
function BLi_DSV_loads!(element,aeroStateCache)

    @unpack b = element.aero.geometry
    @unpack μv₂,ν₁,ν₂,ν₃,ν₄,ν₅,χu,χd,gᵥ,gᵥ₂,Tv,Tv₂,Vn₁,Vn₂,Vn₃,Vm,Vt = element.aero.surface.geometry.airfoil.state.bli
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin
    @unpack qR,R,stallOnsetRatio,upstroke,P = aeroStateCache.bliVars.flow
    @unpack fDiff_tv0P,qR_tv0P,R_tv0P,RD_tv0P,upstroke_tv0P,τvP,fDiff_tv0P2,RD_tv0P2,upstroke_tv0P2,fDiff_tv0N,qR_tv0N,R_tv0N,RD_tv0N,upstroke_tv0N,τvN,fDiff_tv0N2,RD_tv0N2,upstroke_tv0N2 = aeroStateCache.bliVars.comp

    # Set current type
    T = typeof(Uᵢ)

    # Positive AoA
    #---------------------------------------------------------------------------

    # First vortex: normal force and pitching moment
    #---------------------------------------------------------------------------
    if τvP <= 2*Tv 
        # Vortex strength 
        Vs₁ = Vn₁ * RD_tv0P^ν₁ * fDiff_tv0P^ν₂ * max(1, min(3, qR_tv0P^ν₃)) 
        # Vortex stretch factor
        χᵥ = upstroke_tv0P ? 
             1 + χu * R_tv0P : 
             1 + χd
        # Vortex shape function
        Vₓ = τvP/Tv <= 1/χᵥ ? 
             sin(π/2*τvP/Tv*χᵥ)^(ν₄/χᵥ+RD_tv0P) : 
             sin(π/2*(1+(χᵥ/(2*χᵥ-1))*(τvP/Tv-1/χᵥ)))^(ν₄*χᵥ+RD_tv0P)
        # Airloads coefficients
        cnV1P = Vs₁ * Vₓ
        cmV1P = -Vm * cnV1P
    else
        cnV1P = zero(T)
        cmV1P = zero(T)
    end

    # First vortex: tangential force 
    #---------------------------------------------------------------------------
    # Nondimensional time at which vortex begins to affect tangential force
    τvT = Tv
    # Increment in vortex convection time at high pitch rates
    δTvc = min(1/2, qR_tv0P-R_tv0P)
    # Airload coefficient
    ctV1P = (τvP > Tv && τvP <= (2+δTvc)*Tv && upstroke_tv0P) ? 
            -Vt * Vn₁ * RD_tv0P^ν₁ * fDiff_tv0P^ν₂ * min(5, qR_tv0P^ν₅) * ((1-cos(2π*(τvP-τvT)/((1+δTvc)*Tv)))/2) : 
            zero(T)

    # Second vortex 
    #---------------------------------------------------------------------------
    # Nondimensional time at which second vortex begins
    τv2 = (2+gᵥ)*Tv
    if τvP > τv2 && τvP <= τv2+2*Tv₂ && fDiff_tv0P2 > 0.05
        # Vortex strength
        Vs₂ = Vn₂ * RD_tv0P^ν₁ * fDiff_tv0P^ν₂ * max(1, min(3, qR_tv0P^ν₃)) * (1 + μv₂ * R_tv0P^2 * upstroke_tv0P2)
        # Vortex shape function
        Vx₂ = (1-cos(2π*(τvP-τv2)/(2*Tv₂)))/2 
        # Airloads coefficients
        cnV2P = Vs₂ * Vx₂
        cmV2P = -Vm * cnV2P * RD_tv0P2^3
    else
        cnV2P = zero(T)
        cmV2P = zero(T)
    end

    # Third vortex 
    #---------------------------------------------------------------------------
    # Nondimensional time at which third vortex begins
    τv3 = τv2+(2+gᵥ₂)*Tv₂
    if τvP > τv3 && τvP <= τv3+2*Tv₂ && fDiff_tv0P2 > 0.1 && upstroke_tv0P2
        # Vortex strength
        Vs3 = Vn₃ * RD_tv0P2 * RD_tv0P^ν₁ * fDiff_tv0P^ν₂ * (1 + μv₂ * R_tv0P^2)
        # Vortex shape function 
        Vx₃ = (1-cos(2π*(τvP-τv3)/(2*Tv₂)))/2
        # Airloads coefficients
        cnV3P = Vs3 * Vx₃ 
        cmV3P = -Vm * cnV3P * RD_tv0P2^3
    else
        cnV3P = zero(T)
        cmV3P = zero(T)
    end 

    # Negative AoA
    #---------------------------------------------------------------------------

    # First vortex: normal force and pitching moment
    #---------------------------------------------------------------------------
    if τvN <= 2*Tv 
        # Vortex strength
        Vs₁ = Vn₁ * RD_tv0N^ν₁ * fDiff_tv0N^ν₂ * max(1, min(3, qR_tv0N^ν₃))
        # Vortex stretch factor
        χᵥ = upstroke_tv0N ? 
             1 + χu * R_tv0N : 
             1 + χd
        # Vortex shape function
        Vₓ = τvN/Tv <= 1/χᵥ ? 
             sin(π/2*τvN/Tv*χᵥ)^(ν₄/χᵥ+RD_tv0N) : 
             sin(π/2*(1+(χᵥ/(2*χᵥ-1))*(τvN/Tv-1/χᵥ)))^(ν₄*χᵥ+RD_tv0N)
        # Airloads coefficients
        cnV1N = -Vs₁ * Vₓ
        cmV1N = -Vm * cnV1N
    else
        cnV1N = zero(T)
        cmV1N = zero(T)
    end

    # First vortex: tangential force 
    #---------------------------------------------------------------------------
    # Time at which vortex begins to affect tangential force
    τvT = Tv
    # Increment in vortex convection time at high pitch rates
    δTvc = min(1/2, qR_tv0N-R_tv0N)
    # Airload coefficient
    ctV1N = (τvN > Tv && τvN <= (2+δTvc)*Tv && upstroke_tv0N) ? 
            -Vt * Vn₁ * RD_tv0N^ν₁ * fDiff_tv0N^ν₂ * min(5, qR_tv0N^ν₅) * ((1-cos(2π*(τvN-τvT)/((1+δTvc)*Tv)))/2) : 
            zero(T)

    # Second vortex 
    #---------------------------------------------------------------------------
    # Nondimensional time at which second vortex begins 
    τv2 = (2+gᵥ)*Tv
    if τvN > τv2 && τvN <= τv2+2*Tv₂ && fDiff_tv0N2 > 0.05
        # Vortex strength
        Vs₂ = Vn₂ * RD_tv0N^ν₁ * fDiff_tv0N^ν₂ * max(1, min(3, qR_tv0N^ν₃)) * (1 + μv₂ * R_tv0N^2 * upstroke_tv0N2)
        # Vortex shape function 
        Vx₂ = (1-cos(2π*(τvN-τv2)/(2*Tv₂)))/2
        # Airloads coefficients
        cnV2N = -Vs₂ * Vx₂
        cmV2N = -Vm * cnV2N * RD_tv0N2^3
    else
        cnV2N = zero(T)
        cmV2N = zero(T)
    end

    # Third vortex 
    #---------------------------------------------------------------------------
    # Nondimensional time at which third vortex begins
    τv3 = τv2+(2+gᵥ₂)*Tv₂
    if τvN > τv3 && τvN <= τv3+2*Tv₂ && fDiff_tv0N2 > 0.1 && upstroke_tv0N2
        # Vortex strength
        Vs3 = Vn₃ * RD_tv0N2 * RD_tv0N^ν₁ * fDiff_tv0N^ν₂ * (1 + μv₂*R_tv0N^2)
        # Vortex shape function 
        Vx₃ = (1-cos(2π*(τvN-τv3)/(2*Tv₂)))/2 
        # Airloads coefficients
        cnV3N = -Vs3 * Vx₃ 
        cmV3N = -Vm * cnV3N * RD_tv0N2^3
    else
        cnV3N = zero(T)
        cmV3N = zero(T)
    end

    # Total DSV loads
    #---------------------------------------------------------------------------
    cnV = cnV1P + cnV2P + cnV3P + cnV1N + cnV2N + cnV3N
    cmV = cmV1P + cmV2P + cmV3P + cmV1N + cmV2N + cmV3N
    ctV = ctV1P + ctV1N

    @pack! aeroStateCache.attachedVars.aeroCoefficients = cnV,cmV,ctV

end


# Computes normal force coefficient for the modified incompressible Beddoes-Leishman model
function BLi_cn!(element,aeroStateCache,χ,δ)

    @unpack flapSolver = element.aero.surface.solvers
    @unpack inertialPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack incompressibleInertialLoads = element.aero.surface.solvers.ppSolver
    @unpack ϵₙ,cnα = element.aero.surface.geometry.airfoil.state.bli
    @unpack cnδ = element.aero.surface.geometry.flapParameters
    @unpack isFlapped,Th,b,c,aₕ = element.aero.geometry
    @unpack Uᵢ,Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack αₑ,UₙdotMid = aeroStateCache.attachedVars.unsteadyKin
    @unpack Ma,Tnα,TnM,Tnθ̇ = aeroStateCache.attachedVars.flowParameters
    @unpack cnV = aeroStateCache.attachedVars.aeroCoefficients
    @unpack f2primeN = aeroStateCache.bliVars.states
    if !incompressibleInertialLoads
        χ1,χ2,χ3,χ4,χ5,χ6,χ7,χ8 = view(χ,inertialPitchPlungeStatesRange)
    end

    # Circulatory component - attached flow
    cnC = cnα * sin(αₑ) * cos(αₑ)
    if isFlapped && flapSolver isa TableLookup
        cnC += cnδ * δ
    end

    # Circulatory component - separated flow
    cnF = cnC * ((1+sqrt(f2primeN))/2)^2

    # Inertial component
    cnI = incompressibleInertialLoads ? 
          π*b*UₙdotMid/Uᵢ^2 : 
          4/Ma * (α - χ1/Tnα) + -4*aₕ*b/(Uᵢ*Ma) * (Ωₐ - χ2/Tnθ̇) + 4*α/Ma^2 * (Ma - χ6/TnM)
    if isFlapped && flapSolver isa ThinAirfoilTheory
        @unpack δdot,δddot = element.aero.surface.state
        cnI -= ϵₙ*b/Uᵢ^2 * (Uᵢ*Th[4]*δdot + b*Th[1]*δddot)
    end

    # Non-circulatory
    cnNC = cnI

    # Total 
    cn = cnF + cnNC + cnV

    @pack! aeroStateCache.attachedVars.aeroCoefficients = cn,cnC,cnNC,cnI,cnF,cnV

end


# Computes the pitching moment aerodynamic coefficient at the attachment point (i.e., the beam reference line) for the modified incompressible Beddoes-Leishman model
function BLi_cm!(element,aeroStateCache,χ,δ)

    @unpack flapSolver = element.aero.surface.solvers
    @unpack inertialPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack incompressibleInertialLoads,AI,bI = element.aero.surface.solvers.ppSolver
    @unpack ϵₘ,κ₀,κ₁,κ₂,κ₃,cm₀,cnα,K₀,K₁,K₂ = element.aero.surface.geometry.airfoil.state.bli
    @unpack cmδ = element.aero.surface.geometry.flapParameters
    @unpack isFlapped,Th,b,normSparPos,aₕ = element.aero.geometry
    @unpack Uᵢ,Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack UₙdotMid,Ωₐdot = aeroStateCache.attachedVars.unsteadyKin
    @unpack Ma,βₚ,Tmα,TmM,Tmθ̇ = aeroStateCache.attachedVars.flowParameters
    @unpack cnF,cmV = aeroStateCache.attachedVars.aeroCoefficients
    @unpack f2primeM,RD = aeroStateCache.bliVars.states
    @unpack R,stallOnsetRatio,upstroke,stalled,P = aeroStateCache.bliVars.flow
    if !incompressibleInertialLoads
        χ1,χ2,χ3,χ4,χ5,χ6,χ7,χ8 = view(χ,inertialPitchPlungeStatesRange)
    end

    # Center of pressure dynamic variables
    K1Prime = K₁ * (1 - κ₁ * RD * (1-abs(stallOnsetRatio))) - 
              κ₂ * R * abs(stallOnsetRatio) * upstroke
    K2Prime = K₂ * (1 + κ₃ * stalled * R^2 * (1-P) *(!upstroke))

    # Center of pressure offset from quarter-chord
    δCP = K₀ + K1Prime * (1-f2primeM) + K2Prime * sin(π*f2primeM^κ₀)

    # Circulatory component - separated flow
    cmF = cnF * (normSparPos - (1/4-δCP))
    if isFlapped && flapSolver isa TableLookup
        cmF += cmδ * δ
    end

    # Inertial component
    cmI = incompressibleInertialLoads ? 
          -π*b/(2*Uᵢ^2)*(-aₕ*UₙdotMid+b*Ωₐdot/8+Uᵢ*Ωₐ/2) : 
          2*aₕ/Ma * (α - χ3*AI[1]/(bI[1]*Tmα) - χ4*AI[2]/(bI[2]*Tmα)) + -2*b*(aₕ^2+1/3)/(Uᵢ*Ma) * (Ωₐ - χ5/Tmθ̇) + 2*aₕ*α/Ma^2 * (Ma - χ7*AI[1]/(bI[1]*TmM) - χ8*AI[2]/(bI[2]*TmM)) - π*b*Ωₐ/(4*Uᵢ*βₚ)
    if isFlapped && flapSolver isa ThinAirfoilTheory
        @unpack δdot,δddot = element.aero.surface.state
        cmI -= ϵₘ/(2*Uᵢ^2) * (Uᵢ^2*Th[14]*δ + Uᵢ*b*Th[15]*δdot + b^2*Th[16]*δddot)
    end

    # Non-circulatory
    cmNC = cmI + cm₀

    # Total at attachment point
    cm = cmF + cmNC + cmV

    @pack! aeroStateCache.attachedVars.aeroCoefficients = cm,cmF,cmI,cmNC

end


# Computes the tangential force aerodynamic coefficient for the modified incompressible Beddoes-Leishman model
function BLi_ct!(element,aeroStateCache,δ)

    @unpack flapSolver = element.aero.surface.solvers
    @unpack hasInducedDrag = element.aero.surface.geometry
    @unpack α1₀T,η,cd₀,cnα,E₀,E₁ = element.aero.surface.geometry.airfoil.state.bli
    @unpack cdδ = element.aero.surface.geometry.flapParameters
    @unpack isFlapped,b = element.aero.geometry
    @unpack Uᵢ,Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack αₑ = aeroStateCache.attachedVars.unsteadyKin
    @unpack cnF,ctV = aeroStateCache.attachedVars.aeroCoefficients
    @unpack αlag,f2primeT,RD = aeroStateCache.bliVars.states
    @unpack R,stallOnsetRatio,upstroke,stalled = aeroStateCache.bliVars.flow
    @unpack Ts = aeroStateCache.bliVars.comp

    # Ratio of lagged AoA to tangential steady breakpoint angle
    stallOnsetRatioT = abs(αlag) / α1₀T

    # Parasite component
    ct₀ = -cd₀ / cos(α)

    # Circulatory component - separated flow
    ctF = (1-η*RD^2) * cnF * tan(αₑ+b*Ωₐ/(2*Uᵢ)) * 
          f2primeT^(1/2 + 
                    stallOnsetRatioT + 
                    E₀ * RD * stalled * (1-Ts) * R * abs(stallOnsetRatio)^(1/2) * (!upstroke)
                    )
    if stallOnsetRatioT > 1
        ctF -= E₁ * (1-RD^3) * min(1, stallOnsetRatioT^3-1)
    end
    if isFlapped && flapSolver isa TableLookup
        ctF -= cdδ * abs(δ) / cos(α)
    end

    # Induced (drag) component
    if hasInducedDrag
        @unpack AR = element.aero.surface.assemblyData
        ctF -= cnF^2 / (π*AR) / cos(α)
    end

    # Non-circulatory
    ctNC = ct₀

    # Total
    ct = ctF + ctNC + ctV

    @pack! aeroStateCache.attachedVars.aeroCoefficients = ct,ctF,ctV,ctNC

end


# Scales the circulatory components of aerodynamic coefficients and recompute the totals for the BLi ppSolver
function scale_BLi_aero_coefficients!(element,aeroStateCache)

    # Unpack
    @unpack ϖMid = element.aero.loadsDistribution
    @unpack cnC,cmC,ctC,cnNC,cmNC,ctNC,cnF,cmF,ctF,cnV,cmV,ctV = aeroStateCache.attachedVars.aeroCoefficients

    # Scale circulatory components
    cnC,cmC,ctC,cnF,cmF,ctF,cnV,cmV,ctV = multiply_scalars!(ϖMid, cnC,cmC,ctC,cnF,cmF,ctF,cnV,cmV,ctV)

    # Recompute totals
    cn = cnF + cnNC + cnV
    cm = cmF + cmNC + cmV
    ct = ctF + ctNC + ctV

    # Set the scaled aerodynamic coefficients
    @pack! aeroStateCache.attachedVars.aeroCoefficients = cn,cm,ct,cnC,cmC,ctC,cnF,cmF,ctF,cnV,cmV,ctV

end


# Computes the aerodynamic coefficients according to the original Beddoes-Leishman model
function BLo_aero_coefficients!(isDynamicProblem,timeNow,Δt,element,aeroStateCache,χ,δ)

    # Name nonlinear states
    BLo_nonlinear_states!(element,aeroStateCache,χ)

    # Motion qualifiers
    BLo_motion_qualifiers!(element,aeroStateCache)

    # Breakpoint angle
    BLo_breakpoint_angle!(element,aeroStateCache)

    # Flow separation points
    BLo_separation_points!(element,aeroStateCache)

    # Variables at stall onset
    isDynamicProblem && BLo_stall_time!(timeNow,Δt,element,aeroStateCache)

    # Time delay constants
    BLo_time_delays!(element,aeroStateCache)

    # Normal force coefficient
    BLo_cn!(element,aeroStateCache,χ,δ)

    # DSV accumulation rate
    BLo_vortex_accumulation_rate!(element,aeroStateCache,χ)

    # Pitching moment coefficient about the spar position
    BLo_cm!(element,aeroStateCache,χ,δ)

    # Tangential flow coefficient
    BLo_ct!(element,aeroStateCache,δ)

    # Scale 
    scale_BLo_aero_coefficients!(element,aeroStateCache)

end


# Sets the nonlinear states of the original Beddoes-Leishman model
function BLo_nonlinear_states!(element,aeroStateCache,χ)

    @unpack nonlinearStatesRange = element.aero.surface.stateRanges
    @unpack f₀ = element.aero.surface.geometry.airfoil.state.blo

    # Aliases
    χnl1,χnl2,χnl3,χnl4,χnl5 = view(χ,nonlinearStatesRange)

    # Nonlinear states
    cnPprime = χnl1
    f2Prime = clamp(χnl2, f₀, 1.0)
    fPrimeM = clamp(χnl3, f₀, 1.0)
    cnVP = χnl4
    cnVN = χnl5

    @pack! aeroStateCache.bloVars.states = cnPprime,f2Prime,fPrimeM,cnVP,cnVN

end


# Computes motion qualifiers of the original Beddoes-Leishman model
function BLo_motion_qualifiers!(element,aeroStateCache)

    @unpack c = element.aero.geometry
    @unpack cn₁,cnα = element.aero.surface.geometry.airfoil.state.blo
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin
    @unpack αdot = aeroStateCache.attachedVars.unsteadyKin
    @unpack cnPprime = aeroStateCache.bloVars.states

    # Non-dimensional pitch rate (round off to supress noise)
    q = round_off(αdot*c/Uᵢ, 1e-6)

    # Lagged angle of attack
    αlag = cnPprime / cnα

    # Dynamic stall onset ratio
    stallOnsetRatio = cnPprime / cn₁

    # Flag for upstroke phase
    upstroke = stallOnsetRatio * q >= 0

    @pack! aeroStateCache.bloVars.flow = αlag,q,stallOnsetRatio,upstroke

end


# Computes the breakpoint angle of the original Beddoes-Leishman model
function BLo_breakpoint_angle!(element,aeroStateCache)

    @unpack α1₀,δα = element.aero.surface.geometry.airfoil.state.blo
    @unpack upstroke = aeroStateCache.bloVars.flow
    @unpack f2Prime = aeroStateCache.bloVars.states

    # Breakpoint angle offset
    δα1 = upstroke ? 
          zero(f2Prime) : 
          -δα * (1-f2Prime)^(1/4)

    # Unsteady breakpoint of separation angle
    α1 = α1₀ + δα1

    @pack! aeroStateCache.bloVars.flow = α1

end


# Computes the flow separation point of the original Beddoes-Leishman model
function BLo_separation_points!(element,aeroStateCache)

    @unpack α1₀,cnα,f₀,fb,S1,S2 = element.aero.surface.geometry.airfoil.state.blo
    @unpack α = aeroStateCache.attachedVars.steadyKin
    @unpack α1,αlag = aeroStateCache.bloVars.flow
    @unpack cnPprime = aeroStateCache.bloVars.states

    # Absolute values of angles
    absα = abs(α)
    absαlag = abs(αlag)

    # Quasi-steady separation point
    f = absα <= α1₀ ?
        1 - (1-fb) * exp((absα-α1₀)/S1) :
        f₀ + (fb-f₀) * exp((α1₀-absα)/S2)

    # Unsteady lagged separation point based on lagged AoA
    fPrime = absαlag <= α1 ?
             1 - (1-fb) * exp((absαlag-α1)/S1) :
             f₀ + (fb-f₀) * exp((α1-absαlag)/S2)

    @pack! aeroStateCache.bloVars.flow = f,fPrime

end


# Computes variables at the time of stall onset for the original Beddoes-Leishman model
function BLo_stall_time!(timeNow,Δt,element,aeroStateCache)

    @unpack b = element.aero.geometry
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin
    @unpack stallOnsetRatio = aeroStateCache.bloVars.flow
    @unpack stallOnsetRatioPrev,tv0P,tv0N = aeroStateCache.bloVars.comp

    # Time instances
    t1 = timeNow - Δt
    t2 = timeNow

    # Positive AoA
    #---------------------------------------------------------------------------
    # Linear interpolation for time of vortex shedding 
    if stallOnsetRatio >= 1 && stallOnsetRatioPrev < 1 
        tv0P = t1 + (1 - stallOnsetRatioPrev) / (stallOnsetRatio - stallOnsetRatioPrev) * (t2 - t1)
    end
    # Nondimensional time since primary vortex shedding 
    τvP = (timeNow - tv0P) * Uᵢ / b

    # Negative AoA
    #---------------------------------------------------------------------------
    # Linear interpolation for time of vortex shedding
    if stallOnsetRatio <= -1 && stallOnsetRatioPrev > -1 
        tv0N = t1 + (1 - (-stallOnsetRatioPrev)) / ((-stallOnsetRatio) - (-stallOnsetRatioPrev)) * (t2 - t1)
    end
    # Nondimensional time since primary vortex shedding 
    τvN = (timeNow - tv0N) * Uᵢ / b

    @pack! aeroStateCache.bloVars.comp = tv0P,τvP,tv0N,τvN

end


# Computes time delay variables for the original Beddoes-Leishman model
function BLo_time_delays!(element,aeroStateCache)

    @unpack Tf₀,Tv₀,TvL,fb = element.aero.surface.geometry.airfoil.state.blo
    @unpack stallOnsetRatio,upstroke = aeroStateCache.bloVars.flow
    @unpack f2Prime = aeroStateCache.bloVars.states
    @unpack τvP,τvN = aeroStateCache.bloVars.comp

    τv = min(τvP,τvN)
    
    # Vortex shedding phase
    if abs(stallOnsetRatio) >= 1
        # Primary vortex
        if τv <= TvL && upstroke
            # Accelerate the rate of boundary layer dettachment during the vortex convection
            Tf = Tf₀*3/4
            # Nominal rate of vortex accumulation
            Tv = Tv₀
        elseif TvL < τv <= 2*TvL && upstroke
            # Accelerate more the rate of dettachment after the vortex reaches the trailing edge 
            Tf = Tf₀/2
            # Increase the rate of decay of the vortex after it reaches the trailing edge 
            Tv = Tv₀/2
        elseif τv > 2*TvL && upstroke
            # Maintain the rate of boundary layer dettachment after the vortex is totally shed
            Tf = Tf₀/2
            # Maintain a high rate of vortex lift decay after the vortex is totally shed
            Tv = Tv₀/2
        elseif τv <= 2*TvL && !upstroke
            # Maintain the rate of dettachment if the rate of change of AoA changes during the vortex shedding 
            Tf = Tf₀/2
            # Increase the rate of decay of the vortex lift if the rate of change of AoA changes during the vortex shedding 
            Tv = Tv₀/2
        elseif τv > 2*TvL && !upstroke
            # Delay the reattachment of the boundary layer after the vortex is totally shed and the AoA is decreasing
            Tf = 4*Tf₀
            # Maintain a high rate of vortex lift decay after the vortex is totally shed
            Tv = Tv₀/2
        end
    # Reattachment phase    
    else
        # Maintain high rate of vortex lift decay after the vortex is totally shed
        Tv = Tv₀/2
        # Delay the reattachment of the boundary layer
        Tf = 4*Tf₀
        # Dimitriadis' suggestion
        if upstroke && f2Prime >= fb
            # Set to nominal conditions if the rate of change of AoA is increasing and the flow is lightly separated
            Tf = Tf₀
        elseif upstroke && f2Prime < fb
            # Accelerate boundary layer reattachment if the rate of change of AoA is already increasing and the flow is still massively separated
            Tf = Tf₀/2
        end
    end

    @pack! aeroStateCache.bloVars.flow = Tf,Tv

end


# Compute the dynamic stall vortex accumulation rate for the original Beddoes-Leishman model
function BLo_vortex_accumulation_rate!(element,aeroStateCache,χ::AbstractVector{T}) where T<:Number

    @unpack b = element.aero.geometry
    @unpack circulatoryPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack γbC = element.aero.surface.geometry
    @unpack AC,bC = element.aero.surface.solvers.ppSolver
    @unpack cnα,TvL = element.aero.surface.geometry.airfoil.state.blo
    @unpack ϖMid = element.aero.loadsDistribution
    @unpack Uᵢ,Uₜ = aeroStateCache.attachedVars.steadyKin
    @unpack UₙdotTQC = aeroStateCache.attachedVars.unsteadyKin
    @unpack Θ = aeroStateCache.attachedVars.flowParameters
    @unpack cnC = aeroStateCache.attachedVars.aeroCoefficients
    @unpack f2Prime = aeroStateCache.bloVars.states
    @unpack stallOnsetRatio,Kf,fPrime,Tf = aeroStateCache.bloVars.flow
    @unpack τvP,τvN = aeroStateCache.bloVars.comp

    # Precomupte quantities, if vortex is over the chord
    if abs(stallOnsetRatio) >= 1 && (τvP <= TvL || τvN <= TvL)
        # Kirchhoff-Helmholtz factor time rate
        Kfdot = 1/4 * (1+1/sqrt(f2Prime)) * (fPrime-f2Prime) / (Tf*b/Uᵢ)
        # Circulatory unsteady normal force coefficient time rate (assume cnC≈cnα*αₑ, αₑ≈wₑ/Uₜ, and neglect time rate of U and cnα)
        cnCdot = (-Θ*sum(bC .* γbC .* view(χ,circulatoryPitchPlungeStatesRange)) + 
        cnα * UₙdotTQC * (1-sum(AC))
                  ) / Uₜ               
        # Vorticity coefficient rate (scale with loads distribution correction factor)
        cvdot = (cnCdot * ϖMid) * (1-Kf) - (cnC * ϖMid) * Kfdot
    end

    # Zero of current type
    z = zero(T)

    # Positive vorticity coefficient rate (if positive vortex is over the chord)
    cvdotP = (stallOnsetRatio >= 1 && τvP <= TvL) ? 
             max(z,cvdot) :
             z

    # Negative vorticity coefficient rate (if negative vortex is over the chord)
    cvdotN = (stallOnsetRatio <= -1 && τvN <= TvL) ?
             min(z,cvdot) :
             z

    @pack! aeroStateCache.bloVars.flow = cvdotP,cvdotN

end


# Computes normal force coefficient for the original Beddoes-Leishman model
function BLo_cn!(element,aeroStateCache,χ,δ)

    @unpack flapSolver = element.aero.surface.solvers
    @unpack inertialPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack incompressibleInertialLoads = element.aero.surface.solvers.ppSolver
    @unpack ϵₙ,cnα = element.aero.surface.geometry.airfoil.state.blo
    @unpack cnδ = element.aero.surface.geometry.flapParameters
    @unpack isFlapped,Th,b,aₕ = element.aero.geometry
    @unpack Uᵢ,α,Ωₐ = aeroStateCache.attachedVars.steadyKin
    @unpack UₙdotMid,αₑ = aeroStateCache.attachedVars.unsteadyKin
    @unpack Ma,Tnα,TnM,Tnθ̇ = aeroStateCache.attachedVars.flowParameters
    @unpack f2Prime,cnVP,cnVN = aeroStateCache.bloVars.states
    if !incompressibleInertialLoads
        χ1,χ2,χ3,χ4,χ5,χ6,χ7,χ8 = view(χ,inertialPitchPlungeStatesRange)
    end

    # Circulatory component - attached flow
    cnC = cnα * sin(αₑ) * cos(αₑ)
    if isFlapped && flapSolver isa TableLookup
        cnC += cnδ * δ
    end

    # Kirchhoff-Helmholtz factor
    Kf = ((1+sqrt(f2Prime))/2)^2

    # Circulatory component - separated flow
    cnF = cnC * Kf

    # Inertial component
    cnI = incompressibleInertialLoads ?
          π*b*UₙdotMid/Uᵢ^2 : 
          4/Ma * (α - χ1/Tnα) + -4*aₕ*b/(Uᵢ*Ma) * (Ωₐ - χ2/Tnθ̇) + 4*α/Ma^2 * (Ma - χ6/TnM)
    if isFlapped && flapSolver isa ThinAirfoilTheory
        @unpack δdot,δddot = element.aero.surface.state
        cnI -= ϵₙ*b/Uᵢ^2 * (Uᵢ*Th[4]*δdot + b*Th[1]*δddot)
    end

    # Vortex
    cnV = cnVP + cnVN

    # Potential flow component
    cnP = cnC + cnI

    # Non-circulatory
    cnNC = cnI

    # Total 
    cn = cnF + cnNC + cnV

    @pack! aeroStateCache.attachedVars.aeroCoefficients = cn,cnC,cnI,cnF,cnP,cnV,cnNC
    @pack! aeroStateCache.bloVars.flow = Kf

end


# Computes the pitching moment aerodynamic coefficient at the attachment point (i.e., the beam reference line) for the original Beddoes-Leishman model
function BLo_cm!(element,aeroStateCache,χ,δ)

    @unpack flapSolver = element.aero.surface.solvers
    @unpack inertialPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack incompressibleInertialLoads,AI,bI = element.aero.surface.solvers.ppSolver
    @unpack ϵₘ,cm₀,K₀,K₁,K₂,TvL = element.aero.surface.geometry.airfoil.state.blo
    @unpack cmδ = element.aero.surface.geometry.flapParameters
    @unpack isFlapped,normSparPos,b,Th,aₕ = element.aero.geometry
    @unpack Uᵢ,Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack UₙdotMid,Ωₐdot = aeroStateCache.attachedVars.unsteadyKin
    @unpack Ma,βₚ,Tmα,TmM,Tmθ̇ = aeroStateCache.attachedVars.flowParameters
    @unpack cnC,cnF = aeroStateCache.attachedVars.aeroCoefficients
    @unpack f2Prime,fPrimeM,cnVP,cnVN = aeroStateCache.bloVars.states
    @unpack τvP,τvN = aeroStateCache.bloVars.comp
    if !incompressibleInertialLoads
        χ1,χ2,χ3,χ4,χ5,χ6,χ7,χ8 = view(χ,inertialPitchPlungeStatesRange)
    end

    # Circulatory unsteady - attached flow
    cmC = cnC * (normSparPos - (1/4-K₀))
    if isFlapped && flapSolver isa TableLookup
        cmC += cmδ * δ
    end

    # Lagged separation point and center of pressure offset from quarter-chord
    f2PrimeM = max(f2Prime, fPrimeM)
    δCP = K₀ + K₁ * (1-f2PrimeM) + K₂ * sin(π*f2PrimeM^2)

    # Circulatory unsteady - separated flow
    cmF = cnF * (normSparPos - (1/4-δCP))
    if isFlapped && flapSolver isa TableLookup
        cmF += cmδ * δ
    end

    # Inertial
    cmI = incompressibleInertialLoads ? 
          -π*b/(2*Uᵢ^2)*(-aₕ*UₙdotMid+b*Ωₐdot/8+Uᵢ*Ωₐ/2) : 
          2*aₕ/Ma * (α - χ3*AI[1]/(bI[1]*Tmα) - χ4*AI[2]/(bI[2]*Tmα)) + -2*b*(aₕ^2+1/3)/(Uᵢ*Ma) * (Ωₐ - χ5/Tmθ̇) + 2*aₕ*α/Ma^2 * (Ma - χ7*AI[1]/(bI[1]*TmM) - χ8*AI[2]/(bI[2]*TmM)) - π*b*Ωₐ/(4*Uᵢ*βₚ)
    if isFlapped && flapSolver isa ThinAirfoilTheory
        @unpack δdot,δddot = element.aero.surface.state
        cmI -= ϵₘ/(2*Uᵢ^2) * (Uᵢ^2*Th[14]*δ + Uᵢ*b*Th[15]*δdot + b^2*Th[16]*δddot)
    end

    # Vortex
    CPvP = τvP <= 2*TvL ? 
           0.25*(1-cos(π*τvP/TvL)) : 
           0.0
    CPvN = τvN <= 2*TvL ? 
           0.25*(1-cos(π*τvN/TvL)) : 
           0.0
    cmV = - (CPvP * cnVP + CPvN * cnVN)
    
    # Non-circulatory
    cmNC = cmI + cm₀

    # Total at attachment point
    cm = cmF + cmNC + cmV

    @pack! aeroStateCache.attachedVars.aeroCoefficients = cm,cmC,cmF,cmI,cmV,cmNC

end


# Computes the tangential force aerodynamic coefficient for the original Beddoes-Leishman model
function BLo_ct!(element,aeroStateCache,δ)

    @unpack flapSolver = element.aero.surface.solvers
    @unpack hasInducedDrag = element.aero.surface.geometry
    @unpack η,cd₀,cnα,cn₁,Df,E₀ = element.aero.surface.geometry.airfoil.state.blo
    @unpack cdδ = element.aero.surface.geometry.flapParameters
    @unpack isFlapped,b = element.aero.geometry
    @unpack Uᵢ,Ωₐ,α = aeroStateCache.attachedVars.steadyKin
    @unpack αₑ = aeroStateCache.attachedVars.unsteadyKin
    @unpack cnF = aeroStateCache.attachedVars.aeroCoefficients
    @unpack cnPprime,f2Prime = aeroStateCache.bloVars.states
    @unpack stallOnsetRatio = aeroStateCache.bloVars.flow

    # Stall factor
    ϕ = abs(stallOnsetRatio) < 1 ? 
        zero(α) : 
        Df * (abs(cnPprime) - cn₁)

    # Parasite component
    ct₀ = -cd₀ / cos(α)

    # Circulatory component - separated flow
    ctF = η * cnF * tan(αₑ+b*Ωₐ/(2*Uᵢ)) * (sqrt(f2Prime)*f2Prime^ϕ-E₀)
    if isFlapped && flapSolver isa TableLookup
        ctF -= cdδ * abs(δ) / cos(α)
    end

    # Induced (drag) component
    if hasInducedDrag
        @unpack AR = element.aero.surface.assemblyData
        ctF -= cnF^2 / (π*AR) / cos(α)
    end

    # Non-circulatory
    ctNC = ct₀
    
    # Total
    ct = ctF + ctNC

    @pack! aeroStateCache.attachedVars.aeroCoefficients = ct,ctF,ctNC

end


# Scales the circulatory components of aerodynamic coefficients and recompute the totals for the BLo ppSolver
function scale_BLo_aero_coefficients!(element,aeroStateCache)

    # Unpack
    @unpack ϖMid = element.aero.loadsDistribution
    @unpack cnC,cmC,ctC,cnNC,cmNC,ctNC,cnF,cmF,ctF,cnV,cmV,cnI = aeroStateCache.attachedVars.aeroCoefficients

    # Scale circulatory components
    cnC,cmC,ctC,cnF,cmF,ctF,cnV,cmV = multiply_scalars!(ϖMid, cnC,cmC,ctC,cnF,cmF,ctF,cnV,cmV)

    # Recompute totals
    cn = cnF + cnNC + cnV
    cm = cmF + cmNC + cmV
    ct = ctF + ctNC
    cnP = cnC + cnI

    # Set the scaled aerodynamic coefficients
    @pack! aeroStateCache.attachedVars.aeroCoefficients = cn,cm,ct,cnC,cmC,ctC,cnF,cmF,ctF,cnV,cmV,cnP

end


# Computes the aerodynamic nodal loads resultants
function aero_loads_resultants!(f1χ,f2χ,m1χ,m2χ,atmosphere,element,aeroStateCache,p)

    @unpack ρ = atmosphere
    @unpack Δℓ = element.geometry
    @unpack ppSolver = element.aero.surface.solvers
    @unpack hasTipCorrection = element.aero.surface.tipLoss
    @unpack csecΛ,cosΛ²,RwR0 = element.aero.geometry
    @unpack ϖ,ϖMid,Uref = element.aero.loadsDistribution
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin
    @unpack Uᵢdot = aeroStateCache.attachedVars.unsteadyKin
    @unpack ct,cn,cm,ctNC,cnNC,cmNC = aeroStateCache.attachedVars.aeroCoefficients
    @unpack Fa = aeroStateCache.loadsResultants

    # Rotation tensor from basis A to basis B
    R = rotation_tensor_WM(p)

    # Scaled dynamic pressure factor
    fac = 1/2 * ρ * Uᵢ^2 * csecΛ * Δℓ

    # Loads resultants array, resolved in basis W
    if hasTipCorrection
        # Aerodynamic loads at quadrature points
        f1v = aero_load_at_ζ(ζ1, fac, ϖ, Uref, ϖMid, cosΛ², csecΛ, ct, cn, cm, ctNC, cnNC, cmNC)
        f2v = aero_load_at_ζ(ζ2, fac, ϖ, Uref, ϖMid, cosΛ², csecΛ, ct, cn, cm, ctNC, cnNC, cmNC)
        f3v = aero_load_at_ζ(ζ3, fac, ϖ, Uref, ϖMid, cosΛ², csecΛ, ct, cn, cm, ctNC, cnNC, cmNC)
        # Aerodynamic loads resultants
        Fa[2:4]  .= w1*f1v*(1-ζ1) + w2*f2v*(1-ζ2) + w3*f3v*(1-ζ3)
        Fa[8:10] .= w1*f1v*ζ1     + w2*f2v*ζ2     + w3*f3v*ζ3
    else
        # Split midpoint value equally between the two nodes
        fMid = aero_load_at_ζ(0.5, fac, ϖ, Uref, ϖMid, cosΛ², csecΛ, ct, cn, cm, ctNC, cnNC, cmNC)
        Fa[2:4]  .= fMid / 2
        Fa[8:10] .= fMid / 2
    end

    # Aerodynamic contributions to nodal resultants, resolved in basis A
    f1χ .= R * RwR0 * view(Fa, 1:3)
    m1χ .= R * RwR0 * view(Fa, 4:6)
    f2χ .= R * RwR0 * view(Fa, 7:9)
    m2χ .= R * RwR0 * view(Fa, 10:12)

end


# Aerodynamic loads per unit length, as a function of the local element coordinate (ζ)
@inline function aero_load_at_ζ(ζ, fac, ϖ, Uref, ϖMid, cosΛ², csecΛ, ct, cn, cm, ctNC, cnNC, cmNC)
    k = ϖ(Uref, ζ) / ϖMid * cosΛ²
    return fac * SVector(
        k*(ct-ctNC) + ctNC,
        k*(cn-cnNC) + cnNC,
        csecΛ*(k*(cm-cmNC) + cmNC)
    )
end


# Computes the effective (unsteady) pitch-plunge-induced normalwash for the QuasiSteady ppSolver
function pitch_plunge_effective_normalwash(element,aeroStateCache,χ,ppSolver::QuasiSteady)
    return aeroStateCache.attachedVars.steadyKin.UₙTQC
end


# Computes the effective (unsteady) pitch-plunge-induced normalwash for the Indicial solvers
function pitch_plunge_effective_normalwash(element,aeroStateCache,χ,ppSolver::IndicialSolver)

    @unpack airfoil = element.aero.surface.geometry
    @unpack circulatoryPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack UₙTQC = aeroStateCache.attachedVars.steadyKin

    cnα = get_airfoil_state_field(airfoil, ppSolver, :cnα)

    return UₙTQC - sum(view(χ,circulatoryPitchPlungeStatesRange)) / cnα

end


# Computes the effective (unsteady) pitch-plunge-induced normalwash for the Inflow ppSolver
function pitch_plunge_effective_normalwash(element,aeroStateCache,χ,ppSolver::Inflow)

    @unpack circulatoryPitchPlungeStatesRange = element.aero.surface.stateRanges
    @unpack UₙTQC = aeroStateCache.attachedVars.steadyKin

    return UₙTQC - 1/2*dot(ppSolver.bₚ, view(χ,circulatoryPitchPlungeStatesRange))

end


# Computes the instantaneous (quasi-steady) flap-induced normalwash
function flap_normalwash(element,aeroStateCache,δ)

    @unpack δdot = element.aero.surface.state
    @unpack b,Th = element.aero.geometry
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin

    wf = isnothing(Th) ? 
         0.0 : 
         Th[17]*Uᵢ*δ + b*Th[18]*δdot

    return wf

end


# Computes the time rate of the instantaneous flap-induced normalwash
function flap_normalwash_rate(element,aeroStateCache,δ)

    @unpack δdot,δddot = element.aero.surface.state
    @unpack b,Th = element.aero.geometry
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin
    @unpack Uᵢdot = aeroStateCache.attachedVars.unsteadyKin

    wfdot = Th[17]*(Uᵢ*δdot + Uᵢdot*δ) + b*Th[18]*δddot

    return wfdot

end


# Computes the effective (unsteady) flap-induced normalwash
function flap_effective_normalwash(element,χ,wf)

    @unpack ppSolver = element.aero.surface.solvers
    @unpack airfoil = element.aero.surface.geometry
    @unpack flapStatesRange = element.aero.stateRanges
    cnα = get_airfoil_state_field(airfoil, ppSolver, :cnα)

    wₑf = wf - sum(view(χ, flapStatesRange)) / cnα

    return wₑf

end


# Computes the effective (unsteady) gust-induced normalwash
function gust_effective_normalwash(element,aeroStateCache,χ)

    @unpack gustStatesRange = element.aero.stateRanges
    @unpack AGbG = element.aero.surface.solvers.gustSolver
    @unpack UₙGust = aeroStateCache.attachedVars.steadyKin
    @unpack Θ = aeroStateCache.attachedVars.flowParameters

    wₑg = isempty(gustStatesRange) ? 
          UₙGust : 
          Θ * dot(AGbG, view(χ, gustStatesRange))

    return wₑg

end


# Computes the time rate of the normal force coefficient slope, cnα
function cnα_rate(element,aeroStateCache)

    @unpack ppSolver = element.aero.surface.solvers
    @unpack airfoil = element.aero.surface.geometry
    @unpack Ma,βₚ,βₚ² = aeroStateCache.attachedVars.flowParameters
    @unpack Uᵢ = aeroStateCache.attachedVars.steadyKin
    @unpack Uᵢdot = aeroStateCache.attachedVars.unsteadyKin
    cnα = get_airfoil_state_field(airfoil, ppSolver, :cnα)

    # Time derivative of βₚ
    βₚdot = - Uᵢ * Uᵢdot / ((Uᵢ/Ma)^2 * βₚ)

    # Time derivative of cnα (assuming it scales with 1/βₚ)
    cnαdot = cnα * βₚdot / βₚ²

    return cnαdot

end


# Computes the time rate of the product of cnα by UₙTQC
function cnαUₙTQC_rate(element,aeroStateCache)

    @unpack ppSolver = element.aero.surface.solvers
    @unpack airfoil = element.aero.surface.geometry
    @unpack UₙTQC = aeroStateCache.attachedVars.steadyKin
    @unpack UₙdotTQC = aeroStateCache.attachedVars.unsteadyKin
    cnα = get_airfoil_state_field(airfoil, ppSolver, :cnα)

    # Time derivative of cnα
    cnαdot = cnα_rate(element,aeroStateCache)

    # Time derivative of the product cnα * UₙTQC
    cnαUₙTQCdot = cnα * UₙdotTQC + cnαdot * UₙTQC

    return cnαUₙTQCdot

end


# Computes the time rate of the product of cnα by wf
function cnαwf_rate(element,aeroStateCache,δ)

    @unpack ppSolver = element.aero.surface.solvers
    @unpack airfoil = element.aero.surface.geometry
    @unpack wf = aeroStateCache.attachedVars.steadyKin
    cnα = get_airfoil_state_field(airfoil, ppSolver, :cnα)

    # Time derivative of wf
    wfdot = flap_normalwash_rate(element,aeroStateCache,δ)

    # Time derivative of cnα
    cnαdot = cnα_rate(element,aeroStateCache)

    # Time derivative of the product cnα * wf
    cnαwfdot = cnα * wfdot + cnαdot * wf

    return cnαwfdot

end


# Updates the initial aerodynamic states assuming their rates are zero
function initial_aero_states!(atmosphere,element,χ,V,Ω,Vdot,Ωdot)

    # Skip if there are no aerodynamic surface or states
    (isnothing(element.aero.surface) || isempty(χ)) && return χ

    # Unpack
    @unpack stateCache = element.aero
    @unpack Ac,Ai,An,Af,Ag,A,B = stateCache.stateMatrices
    @unpack δ = element.aero.surface.state

    # Initialize time variables and set dynamic problem flag
    timeNow = Δt = 0.0
    isDynamicProblem = false

    # Run twice to update nonlinear states (necessary only for nonlinear solvers)
    for i=1:2

        # Steady aerodynamic kinematics
        α = aero_steady_kinematics!(element,stateCache,V,Ω)

        # Skip if the angle of attack is undefined
        isnan(α) && return χ
        
        # Unsteady aerodynamic kinematics
        aero_unsteady_kinematics!(element,stateCache,Vdot,Ωdot)

        # Flow parameters
        flow_parameters!(atmosphere,element,stateCache)

        # Update airfoil parameters
        update_airfoil_parameters!(element,stateCache)

        # Update inertial parameters
        update_inertial_parameters!(element,stateCache)

        # Effective (unsteady) angle of attack
        effective_angle_of_attack!(element,stateCache,χ,δ)

        # Aerodynamic coefficients
        aero_coefficients!(isDynamicProblem,timeNow,Δt,element,stateCache,χ,δ)

        # State matrices
        aero_state_matrices!(Ac,Ai,An,Af,Ag,A,B,element,stateCache,δ)

        # Set initial aerodynamic states
        χ .= - A \ B

    end

    return χ

end