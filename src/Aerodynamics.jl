"""
aero_steady_kinematics!(element::Element,V,Ω)

Computes the steady aerodynamic kinematic variables

# Arguments
- element::Element
- V
- Ω
"""
function aero_steady_kinematics!(element::Element,V,Ω)

    @unpack RwT,c,normSparPos = element.aero

    # Relative wind velocity vector at spar position, resolved in basis W
    U = -RwT*V
    
    # Relative wind speed 
    U∞ = norm(U)

    # Spanwise, tangential and normal components of relative wind velocity
    Uₛ,Uₜ,Uₙ = U[1],U[2],U[3]

    # In-plane relative wind speed 
    Uᵢ = sqrt(Uₜ^2+Uₙ^2)

    # Angle of attack  
    α = abs(Uₜ) > 0 ? atan(-Uₙ/Uₜ) : NaN64

    # Sideslip angle 
    β = abs(Uᵢ) > 0 ? atan(Uₛ/Uᵢ) : NaN64

    # Spanwise angular velocity component, resolved in basis W
    Ωₐ = dot(RwT[1,:],Ω)

    # Relative normal wind velocity component at the airfoil's 1/2- and 3/4-chord
    UₙMid = Uₙ - Ωₐ*c*(normSparPos-1/2)
    UₙTQC = Uₙ - Ωₐ*c*(normSparPos-3/4)

    @pack! element.aero.flowAnglesAndRates = α,β
    @pack! element.aero.flowVelocitiesAndRates = U,U∞,Uₛ,Uₜ,Uₙ,Uᵢ,Ωₐ,UₙMid,UₙTQC

end


"""
aero_unsteady_kinematics!(element::Element,Vdot,Ωdot)

Computes the unsteady aerodynamic kinematic variables

# Arguments
- element::Element
- Vdot
- Ωdot
"""
function aero_unsteady_kinematics!(element::Element,Vdot,Ωdot)

    @unpack RwT,c,normSparPos = element.aero
    @unpack Uₜ,Uₙ,Uᵢ = element.aero.flowVelocitiesAndRates

    # Time derivative of relative wind velocity vector, resolved in basis W
    Udot = -RwT*Vdot

    # Time derivatives of tangential and normal components of relative wind velocity
    Uₜdot,Uₙdot = Udot[2],Udot[3]

    # Time derivative of in-plane relative wind speed
    Uᵢdot = (Uₜ*Uₜdot+Uₙ*Uₙdot)/Uᵢ

    # Time derivative of angle of attack 
    αdot = (-Uₜ*Uₙdot+Uₙ*Uₜdot)/Uᵢ^2

    # Spanwise angular acceleration component, resolved in basis W
    Ωₐdot = dot(RwT[1,:],Ωdot)

    # Relative normal wind acceleration component at the airfoil's 1/2- and 3/4-chord
    UₙdotMid = Uₙdot - Ωₐdot*c*(normSparPos-1/2)
    UₙdotTQC = Uₙdot - Ωₐdot*c*(normSparPos-3/4)

    @pack! element.aero.flowAnglesAndRates = αdot
    @pack! element.aero.flowVelocitiesAndRates = Udot,Uₜdot,Uₙdot,Uᵢdot,Ωₐdot,UₙdotMid,UₙdotTQC

end


"""
nondimensional_flow_parameters!(model::Model,element::Element)

Computes the nondimensional flow parameters

# Arguments
- model::Model
- element::Element
"""
function nondimensional_flow_parameters!(model::Model,element::Element)

    @unpack c = element.aero
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates
    @unpack ρ,μ,a = model.atmosphere

    # Reynolds
    Re = ρ*Uᵢ*c/μ

    # Mach
    Ma = Uᵢ/a  

    # Prandtl-Glauert compressibility factor and its square
    βₚ² = Ma < 1 ? 1-Ma^2 : 1.0
    βₚ = sqrt(βₚ²)

    @pack! element.aero.flowParameters = Re,Ma,βₚ,βₚ²

end


"""
local_gust_velocity!(problem::Problem,model::Model,element::Element)

Computes the gust velocity in the local, deformed aerodynamic basis (basis W)

# Arguments
- problem::Problem
- model::Model
- element::Element
"""
function local_gust_velocity!(problem::Problem,model::Model,element::Element)

    @unpack timeNow = problem
    @unpack R_AT,gust = model
    @unpack R = element
    @unpack RwR0 = element.aero
    @unpack UₜGust,UₙGust = element.aero.flowVelocitiesAndRates

    # Reset gust velocity components
    UₜGust,UₙGust = 0,0

    # Gust defined over time
    if gust.isDefinedOverTime && (gust.initialTime <= timeNow <= gust.finalTime)
        @unpack UGustInertial = gust
        # Transform gust velocity vector from inertial basis to local deformed aerodynamic basis W
        UGust = (R*RwR0)'*R_AT * UGustInertial(timeNow)
        UₜGust,UₙGust = UGust[2],UGust[3]
    # Gust defined on space    
    elseif !gust.isDefinedOverTime
    end

    @pack! element.aero.flowVelocitiesAndRates = UₜGust,UₙGust

end


"""
flap_deflection_rates!(problem,element::Element)

Computes the current values of flap deflection rates

# Arguments
- problem
- element::Element
"""
function flap_deflection_rates!(problem,element::Element)

    @unpack timeNow = problem
    @unpack δdot,δddot = element.aero
    δdotNow,δddotNow = δdot(timeNow),δddot(timeNow)

    @pack! element.aero = δdotNow,δddotNow

end


"""
aero_coefficients!(problem::Problem,element::Element,χ,δNow)

Computes the aerodynamic coefficients

# Arguments
- problem::Problem
- element::Element
- χ
- δNow
"""
function aero_coefficients!(problem::Problem,element::Element,χ,δNow)

    @unpack solver = element.aero

    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        attached_flow_aero_coefficients!(element,δNow)
    elseif typeof(solver) == BLi
        BLi_aero_coefficients!(problem,element,χ,δNow)
    elseif typeof(solver) == BLc
        BLc_aero_coefficients!(element,χ,δNow)    
    end   

end


"""
aero_state_matrices!(element::Element,δNow)

Computes the aerodynamic state matrices

# Arguments
- element::Element
- δNow
"""
function aero_state_matrices!(element::Element,δNow)

    @unpack solver = element.aero

    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        A,B = attached_flow_state_matrices!(element,δNow)
    elseif typeof(solver) == BLi
        A,B = BLi_state_matrices!(element,δNow)
    elseif typeof(solver) == BLc
        A,B = BLc_state_matrices!(element,δNow)    
    end   
    
    return A,B

end


"""
aero_loads_resultants!(model::Model,element::Element)

Computes the aerodynamic nodal loads resultants

# Arguments
- model::Model
- element::Element
"""
function aero_loads_resultants!(model::Model,element::Element)

    @unpack ρ = model.atmosphere
    @unpack Δℓ,R = element
    @unpack c,RwR0,Λ,ϖ,ϖMid,hasTipCorrection = element.aero
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates
    @unpack cn,cm,ct = element.aero.aeroCoefficients

    # Aerodynamic loads per unit length, as a function of the local element coordinate (ζ)
    f(ζ) = 1/2*ρ*Uᵢ^2*c*cos(Λ)^2 * [ct; cn; c*cm] * ϖ(ζ)/ϖMid

    # Shape function over element domain
    ϕ(n,ζ) = ifelse.(n==1, 1-ζ, ζ) 

    # Loop over nodes - Set aerodynamic loads array, resolved in basis W
    F = zeros(typeof(cn),12)
    for node=1:2  
        # Nodal loads array: perform integration only if tip correction is present, otherwise split equally among nodes
        if hasTipCorrection
            integrand(ζ) = Δℓ * f(ζ) .* ϕ(node,ζ) 
            F[6*node-4:6*node-2], = quadgk(integrand, 0, 1)
        else
            F[6*node-4:6*node-2] = f(1/2)*Δℓ/2
        end
    end

    # Rotation tensor from basis A to the deformed aerodynamic basis W, resolved in basis A
    RRwR0 = R*RwR0

    # Aerodynamic contributions to nodal resultants loads, resolved in basis A
    f1_χ = RRwR0 * F[1:3]
    m1_χ = RRwR0 * F[4:6]
    f2_χ = RRwR0 * F[7:9]
    m2_χ = RRwR0 * F[10:12]

    @pack! element.aero = F

    return f1_χ,f2_χ,m1_χ,m2_χ

end


"""
attached_flow_aero_coefficients!(element::Element,δNow)

Computes the aerodynamic coefficients under attached flow conditions

# Arguments
- element::Element
- δNow
"""
function attached_flow_aero_coefficients!(element::Element,δNow)

    # Normal force coefficient
    attached_flow_cn!(element,δNow)

    # Pitching moment coefficient about the spar position
    attached_flow_cm!(element,δNow)

    # Tangential flow coefficient
    attached_flow_ct!(element,δNow)

end


"""
update_airfoil_parameters!(problem::Problem,element::Element)

Updates the aerodynamic parameters of the airfoil according to current nondimensional flow parameters 

# Arguments
- problem::Problem
- element::Element
"""
function update_airfoil_parameters!(problem::Problem,element::Element)

    @unpack timeNow = problem
    @unpack solver,airfoil,updateAirfoilParameters,flapSiteID,b = element.aero
    @unpack name = airfoil
    @unpack Re,Ma = element.aero.flowParameters
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates

    if problem isa DynamicProblem
        @unpack Δt,initialTime = problem
    end

    # Skip if parameters are not to be updated
    if (!(problem isa DynamicProblem) && !updateAirfoilParameters) || (problem isa DynamicProblem && !updateAirfoilParameters && timeNow > initialTime + Δt)
        return
    end
    
    # Get airfoil parameters
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        attachedFlowParameters = AttachedFlowParameters(name,Re=Re,Ma=Ma,flapSiteID=flapSiteID)
        @pack! airfoil = attachedFlowParameters
    else
        separatedFlowParameters = SeparatedFlowParameters(name,Re=Re,Ma=Ma,flapSiteID=flapSiteID,U=Uᵢ,b=b)
        @pack! airfoil = separatedFlowParameters
    end

    @pack! element.aero = airfoil

end


"""
effective_angle_of_attack!(element::Element,χ,δNow)

Computes the effective (unsteady) angle of attack

# Arguments
- element::Element
- χ
- δNow
"""
function effective_angle_of_attack!(element::Element,χ,δNow)

    @unpack solver,airfoil = element.aero
    @unpack Uₜ,UₜGust = element.aero.flowVelocitiesAndRates
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        @unpack ϵₙ = airfoil.attachedFlowParameters
    else
        @unpack ϵₙ = airfoil.separatedFlowParameters
    end

    # Effective normalwash: pitch-plunge-induced, flap-induced, gust-induced, and total
    wₑp = pitch_plunge_effective_normalwash(element,χ)
    wₑf = flap_effective_normalwash(element,χ,δNow)
    wₑg = gust_effective_normalwash(element,χ)
    wₑ = wₑp + ϵₙ*wₑf + wₑg

    # Effective circulatory AoA 
    αₑ = atan(-wₑ/(Uₜ+UₜGust))

    @pack! element.aero.flowAnglesAndRates = αₑ

end


"""
attached_flow_cn!(element::Element,δNow)

Computes the normal force aerodynamic coefficient for attached flow

# Arguments
- element::Element
- δNow
"""
function attached_flow_cn!(element::Element,δNow)

    @unpack flapLoadsSolver,flapped,b,δdotNow,δddotNow,ϖMid = element.aero
    @unpack αₑ = element.aero.flowAnglesAndRates
    @unpack Uᵢ,UₙdotMid = element.aero.flowVelocitiesAndRates
    @unpack ϵₙ,cnα,cnδ = element.aero.airfoil.attachedFlowParameters

    # Circulatory component 
    cnC = cnα * sin(αₑ) * cos(αₑ)
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cnC += cnδ*δNow
    end

    # Inertial component
    cnI = π*b*UₙdotMid/Uᵢ^2
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cnI -= ϵₙ*b/Uᵢ^2*(Uᵢ*Th[4]*δdotNow+b*Th[1]*δddotNow)
    end

    # Total 
    cn = cnC + cnI

    # Scale by tip loss correction factor at element's midpoint
    cn,cnC,cnI = multiply_inplace!(ϖMid,cn,cnC,cnI)

    @pack! element.aero.aeroCoefficients = cn,cnC,cnI
end


"""
attached_flow_cm!(element::Element,δNow)

Computes the pitching moment aerodynamic coefficient at the attachment point (i.e., the beam reference line) for attached flow

# Arguments
- element::Element
- δNow
"""
function attached_flow_cm!(element::Element,δNow)

    @unpack flapLoadsSolver,flapped,b,normSparPos,normFlapPos,δdotNow,δddotNow,ϖMid = element.aero
    @unpack αₑ = element.aero.flowAnglesAndRates
    @unpack cnC = element.aero.aeroCoefficients
    @unpack Uᵢ,UₙdotMid,Ωₐ,Ωₐdot = element.aero.flowVelocitiesAndRates
    @unpack ϵₘ,cm₀,cmα,cmδ,cnα = element.aero.airfoil.attachedFlowParameters

    # Circulatory component
    cmC = cm₀ + cmα * αₑ + cnC/ϖMid*(normSparPos-1/4)
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cmC += cmδ*δNow
    end

    # Inertial component
    cmI = -π*b/(2*Uᵢ^2)*((1-2*normSparPos)*UₙdotMid+b*Ωₐdot/8)
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cmI -= ϵₘ/(2*Uᵢ^2)*(Uᵢ^2*Th[14]*δNow+Uᵢ*b*Th[15]*δdotNow+b^2*Th[16]*δddotNow)
    end

    # Rotation-induced component (this is the increment due to the fact that UₙdotMid acts at midchord but the normal force is at the 3/4-chord)
    cmRot = -π/4*b*Ωₐ/Uᵢ

    # Total at attachment point
    cm = cmC + cmI + cmRot

    # Scale by tip loss correction factor at element's midpoint
    cm,cmC,cmI,cmRot = multiply_inplace!(ϖMid,cm,cmC,cmI,cmRot)

    @pack! element.aero.aeroCoefficients = cm,cmC,cmI,cmRot

end


"""
attached_flow_ct!(element::Element,δNow)

Computes the tangential force aerodynamic coefficient for attached flow

# Arguments
- element::Element
- δNow
"""
function attached_flow_ct!(element::Element,δNow)

    @unpack flapped,flapLoadsSolver,ϖMid = element.aero
    @unpack αₑ = element.aero.flowAnglesAndRates
    @unpack cnC = element.aero.aeroCoefficients
    @unpack cd₀,cdδ,cnα = element.aero.airfoil.attachedFlowParameters

    # Circulatory component
    ct = -cd₀/cos(αₑ) + cnα * sin(αₑ)^2
    if flapped && typeof(flapLoadsSolver) == TableLookup
        ct -= cdδ*abs(δNow)/cos(αₑ)
    end

    # Scale by tip loss correction factor at element's midpoint
    ct *= ϖMid

    @pack! element.aero.aeroCoefficients = ct

end


"""
pitch_plunge_effective_normalwash(element::Element,χ)

Computes the effective (unsteady) pitch-plunge-induced normalwash

# Arguments
- element::Element
- χ
"""
function pitch_plunge_effective_normalwash(element::Element,χ)

    @unpack solver,linearPitchPlungeStatesRange = element.aero
    @unpack UₙTQC = element.aero.flowVelocitiesAndRates
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        @unpack cnα = element.aero.airfoil.attachedFlowParameters
    else
        @unpack cnα = element.aero.airfoil.separatedFlowParameters
    end

    # Effective pitch-plunge-induced normalwash
    if typeof(solver) == QuasiSteady
        wₑp = UₙTQC
    elseif typeof(solver) == Indicial
        wₑp = UₙTQC-sum(χ[linearPitchPlungeStatesRange])/cnα
    elseif typeof(solver) == Inflow
        @unpack bₚ = solver
        wₑp = UₙTQC-1/2*dot(bₚ,χ[linearPitchPlungeStatesRange])
    elseif typeof(solver) == BLi
        wₑp = UₙTQC-sum(χ[linearPitchPlungeStatesRange])/cnα    
    end

    return wₑp
end


"""
flap_effective_normalwash(element::Element,χ,δNow)

Computes the effective (unsteady) flap-induced normalwash

# Arguments
- element::Element
- δNow
"""
function flap_effective_normalwash(element::Element,χ,δNow)

    @unpack flapStatesRange = element.aero

    # Skip if there are no flap states
    if isnothing(flapStatesRange)
        return 0
    end
    
    # Quasi-steady flap-induced normalwash
    wFlap = flap_normalwash(element,δNow)

    # Effective flap-induced normalwash
    wₑf = wFlap-sum(χ[flapStatesRange])

    return wₑf
end


"""
flap_normalwash(element::Element,δNow)

Computes the instantaneous (quasi-steady) flap-induced normalwash

# Arguments
- element::Element
- δNow
"""
function flap_normalwash(element::Element,δNow)

    @unpack b,δdotNow = element.aero
    @unpack Th = element.aero.flapLoadsSolver
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates

    wFlap = Th[10]/π*Uᵢ*δNow+b*Th[17]*δdotNow
  
    return wFlap
end


"""
gust_effective_normalwash(element::Element,χ)

Computes the effective (unsteady) gust-induced normalwash

# Arguments
- element::Element
- χ
"""
function gust_effective_normalwash(element::Element,χ)

    @unpack linearGustStatesRange = element.aero

    # Skip if there are not gust states
    if isnothing(linearGustStatesRange)
        return 0
    end

    @unpack solver,b = element.aero
    @unpack βₚ² = element.aero.flowParameters
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates
    @unpack AGbG = element.aero.gustLoadsSolver

    # Effective gust-induced normalwash
    wₑg = Uᵢ/b*βₚ²*dot(AGbG,χ[linearGustStatesRange])

    return wₑg
end


"""
flap_normalwash_rate(element::Element,δNow)

Computes the time rate of the instantaneous flap-induced normalwash

# Arguments
- element::Element
- δNow
"""
function flap_normalwash_rate(element::Element,δNow)

    @unpack b,δdotNow,δddotNow = element.aero
    @unpack Th = element.aero.flapLoadsSolver
    @unpack Uᵢ,Uᵢdot = element.aero.flowVelocitiesAndRates

    wdotFlap = Th[18]*(Uᵢ*δdotNow+Uᵢdot*δNow)+b*Th[17]*δddotNow

    return wdotFlap
end


"""
cnαUₙTQC_rate(element::Element)

Computes the time rate of the product of cnα by UₙTQC

# Arguments
- element::Element
"""
function cnαUₙTQC_rate(element::Element)

    @unpack Ma,βₚ,βₚ² = element.aero.flowParameters
    @unpack Uᵢ,Uᵢdot,UₙTQC,UₙdotTQC = element.aero.flowVelocitiesAndRates
    @unpack cnα = element.aero.airfoil.attachedFlowParameters

    # Time derivative of cnα (assuming it scales with 1/βₚ)
    βₚdot = -Uᵢ*Uᵢdot/((Uᵢ/Ma)^2*βₚ)
    cnαdot = cnα*βₚdot/βₚ^2

    # Time derivative of the product cnα * UₙTQC
    return cnα*UₙdotTQC+cnαdot*UₙTQC
end


"""
attached_flow_state_matrices!(element::Element,δNow)

Computes the aerodynamic state matrices for the indicial and inflow solvers

# Arguments
- element::Element
- δNow
"""
function attached_flow_state_matrices!(element::Element,δNow)

    @unpack solver,nTotalAeroStates,pitchPlungeStatesRange,flapStatesRange,gustStatesRange,b = element.aero
    @unpack βₚ² = element.aero.flowParameters
    @unpack Uᵢ,UₙdotTQC = element.aero.flowVelocitiesAndRates
    @unpack cn = element.aero.aeroCoefficients

    # Initialize state matrices with appropriate types
    A = zeros(typeof(cn),nTotalAeroStates,nTotalAeroStates)
    B = zeros(typeof(cn),nTotalAeroStates)

    # Skip for quasi-steady solver
    if typeof(solver) == QuasiSteady
        return A,B
    end

    tmp = -Uᵢ/b*βₚ²

    # Pitch-plunge-induced flow states
    if typeof(solver) == Indicial
        @unpack AW,bWMat = element.aero.solver
        # Get the rate of cnαUₙTQC
        cnαUₙTQCdot = cnαUₙTQC_rate(element)
        # Set state matrices
        A[pitchPlungeStatesRange,pitchPlungeStatesRange] .= tmp*bWMat
        B[pitchPlungeStatesRange] .= cnαUₙTQCdot*AW
    elseif typeof(solver) == Inflow
        @unpack AₚInv,AₚInvcₚ = element.aero.solver
        # Set state matrices
        A[pitchPlungeStatesRange,pitchPlungeStatesRange] .= tmp*AₚInv
        B[pitchPlungeStatesRange] .= UₙdotTQC*AₚInvcₚ
    end

    # Flap-induced flow states
    if !isnothing(flapStatesRange)
        @unpack AWf,bWfMat = element.aero.flapLoadsSolver
        # Get flap normalwash rate
        wdotFlap = flap_normalwash_rate(element,δNow)
        # Set state matrices
        A[flapStatesRange,flapStatesRange] .= tmp*bWfMat
        B[flapStatesRange] .= wdotFlap*AWf
    end

    # Gust-induced flow states
    if !isnothing(gustStatesRange)
        @unpack bGMat = element.aero.gustLoadsSolver
        @unpack UₙGust = element.aero.flowVelocitiesAndRates
        # Set state matrices
        A[gustStatesRange,gustStatesRange] .= tmp*bGMat
        B[gustStatesRange] .= UₙGust*ones(length(gustStatesRange))
    end

    @pack! element.aero = A,B

    return A,B
end


"""
BLi_aero_coefficients!(problem::Problem,element::Element,χ,δNow)

Computes the aerodynamic coefficients according to the modified incompressible flow Beddoes-Leishman model

# Arguments
- problem::Problem
- element::Element
- χ
- δNow
"""
function BLi_aero_coefficients!(problem::Problem,element::Element,χ,δNow)

    # Name nonlinear states
    BL_nonlinear_states!(element,χ)

    # Additional kinematics
    BL_kinematics!(element)

    # Motion qualifiers
    BL_motion_qualifiers!(problem,element)

    # Breakpoint angles
    BL_breakpoint_angles!(element)

    # Time delay constants
    BL_time_delays!(element)

    # Flow separation points
    BL_separation_points!(element)

    # For dynamic problems
    if problem isa DynamicProblem

        # Variables at stall onset
        BL_stall_time!(problem,element)

        # DSV loads
        BL_DSV_loads!(element)
    end

    # Normal force coefficient
    BLi_cn!(element,δNow)

    # Pitching moment coefficient about the spar position
    BLi_cm!(element,δNow)

    # Tangential flow coefficient
    BLi_ct!(element,δNow)

end


"""
BL_kinematics!(element::Element)

Computes additional kinematics associated with the Beddoes-Leishman model

# Arguments
- element::Element
"""
function BL_kinematics!(element::Element)

    @unpack c = element.aero
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates
    @unpack αdot = element.aero.flowAnglesAndRates
    @unpack r₀ = element.aero.airfoil.separatedFlowParameters

    # Reduced non-dimensional pitch rate (round off to supress noise)
    r = round_off!(αdot*c/2/Uᵢ,1e-12)

    # Unsigned ratio of reduced pitch rate to critical pitch rate
    qR = abs(r)/r₀

    # Unsigned capped reduced pitch rate ratio 
    R = min(1, qR)

    @pack! element.aero.BLkin = r,qR,R

end


"""
BL_nonlinear_states!(element::Element,χ)

Sets the nonlinear states of the Beddoes-Leishman model

# Arguments
- element::Element
- χ
"""
function BL_nonlinear_states!(element::Element,χ)

    @unpack nonlinearPitchPlungeStatesRange,nonlinearGustStatesRange = element.aero
    @unpack Uₜ,UₜGust = element.aero.flowVelocitiesAndRates

    if !isnothing(nonlinearGustStatesRange)
        αlag = abs(Uₜ+UₜGust) > 0 ? atan(-(χ[nonlinearPitchPlungeStatesRange[1]]+χ[nonlinearGustStatesRange[1]])/(Uₜ+UₜGust)) : 0.0
    else
        αlag = abs(Uₜ) > 0 ? atan(-χ[nonlinearPitchPlungeStatesRange[1]]/Uₜ) : 0.0
    end
    f2primeN = χ[nonlinearPitchPlungeStatesRange[2]]
    f2primeM = χ[nonlinearPitchPlungeStatesRange[3]]
    f2primeT = χ[nonlinearPitchPlungeStatesRange[4]]
    RD = χ[nonlinearPitchPlungeStatesRange[5]]
    RD_stallOnsetRatio = χ[nonlinearPitchPlungeStatesRange[6]]

    @pack! element.aero.BLstates = αlag,f2primeN,f2primeM,f2primeT,RD,RD_stallOnsetRatio

end


"""
BL_motion_qualifiers!(problem::Problem,element::Element)

Computes motion qualifiers of the modified Beddoes-Leishman model

# Arguments
- problem::Problem
- element::Element
"""
function BL_motion_qualifiers!(problem::Problem,element::Element)

    @unpack timeNow = problem
    @unpack αds₀,αₛₛ,γLS,Tf = element.aero.airfoil.separatedFlowParameters
    @unpack r,qR = element.aero.BLkin
    @unpack αlag,RD,RD_stallOnsetRatio = element.aero.BLstates
    @unpack αlagPrev,stallOnsetRatioPrev,PPrev,qRPrev,upstrokePrev,maxStallOnsetRatio,minStallOnsetRatio,qRmax,tv0P,tv0N,RD_tv0P,RD_tv0N = element.aero.BLcompVars

    # Critical angle for dynamic stall
    αcrit = αₛₛ + (αds₀-αₛₛ)*RD_stallOnsetRatio

    # Dynamic stall onset ratio
    stallOnsetRatio = αlag/αcrit

    # TF for upstroke phase
    upstroke = stallOnsetRatio * r >= 0
    
    # TF for upstroke phase beginning at current time step
    upstrokeBeginning = upstroke && !upstrokePrev
    
    # TF for stalled condition on current cycle
    S = abs(maxStallOnsetRatio) > 1

    # Light stall factor
    P = S == true ? exp(-γLS*(abs(maxStallOnsetRatio)^4-1)) : PPrev 

    # Steady condition factor
    T = exp(-50*RD)

    # Update maximum stall onset ratio and corresponding delayed capped reduced pitch rate (update only for increments in the lagged AoA)
    if abs(stallOnsetRatio) > abs(stallOnsetRatioPrev) && abs(αlag) > abs(αlagPrev) 
        maxStallOnsetRatio = stallOnsetRatio
        @pack! element.aero.BLcompVars = maxStallOnsetRatio
    end

    # Update minimum stall onset ratio (approximate as that at begin of upstroke)
    if upstrokeBeginning
        minStallOnsetRatio = stallOnsetRatio
        @pack! element.aero.BLcompVars = minStallOnsetRatio
    end

    # Update maximum reduced pitch rate ratio
    if qR > qRPrev
        qRmax = qR
        @pack! element.aero.BLcompVars = qRmax
    end

    # Time decay function for tangential force
    tv0 = max(tv0P, tv0N)
    lastRD_tv0 = tv0P > tv0N ? RD_tv0P : RD_tv0N
    Ts = exp(-5*(timeNow-tv0)*max(0.1,lastRD_tv0)^2/Tf)

    @pack! element.aero.BLflow = stallOnsetRatio,upstroke,S,P,T
    @pack! element.aero.BLcompVars = Ts,lastRD_tv0

end


"""
BL_breakpoint_angles!(element::Element)

Computes the breakpoint angles of the modified Beddoes-Leishman model

# Arguments
- element::Element
"""
function BL_breakpoint_angles!(element::Element)

    @unpack αds₀,αₛₛ,α1₀N,α1₀M,α1₀T,δα₀,δα₁,dt,dm,zm,ztd = element.aero.airfoil.separatedFlowParameters
    @unpack RD = element.aero.BLstates
    @unpack qR,R = element.aero.BLkin
    @unpack upstroke,S,P = element.aero.BLflow
    @unpack qRmax = element.aero.BLcompVars

    # Breakpoint angle offsets
    sqrtOfR = sqrt(R)
    if upstroke 
        fR = R/min(1, qRmax)*sqrtOfR
        δα1N = (αds₀-αₛₛ)*RD        
        δα1M = δα1N+dm*fR                                  
        δα1T = δα1N*+dt*fR                                   
    else 
        δα1N = S == true ? -min(α1₀N/2, δα₀+δα₁*qR*(1+P)) : 0
        δα1M = δα1N*(1-zm*sqrtOfR*(1-P/2))                               
        δα1T = δα1N*(1-ztd*sqrtOfR*(1-P/4))    
    end

    # Unsteady breakpoint of separation angles
    α1N = α1₀N + δα1N 
    α1M = α1₀M + δα1M
    α1T = α1₀T + δα1T

    @pack! element.aero.BLflow = α1N,α1M,α1T

end


"""
BL_time_delays!(element::Element)

Computes time delay variables for the modified Beddoes-Leishman model

# Arguments
- element::Element
"""
function BL_time_delays!(element::Element)

    @unpack Ta,Tf,λ₁,λ₂ = element.aero.airfoil.separatedFlowParameters
    @unpack RD = element.aero.BLstates
    @unpack qR,R = element.aero.BLkin
    @unpack stallOnsetRatio,upstroke,P = element.aero.BLflow
    @unpack qRmax = element.aero.BLcompVars

    # Stall supression time delay
    Ta_SO = Ta*(1/4+λ₁*RD^λ₂*(1-R^(1/4))*exp(-(abs(stallOnsetRatio)-1)^2/0.01))

    # Separation points time delays         
    TfN = Tf*(1+1/2*P^6*(!upstroke)+RD^3*(4*(1-qR/qRmax)*(!upstroke)+1/2*(qR/qRmax)^4*upstroke))
    TfM = Tf*(1-1/2*R*(1-P)*(!upstroke))
    TfT = TfN

    @pack! element.aero.BLflow = Ta_SO,TfN,TfM,TfT

end


"""
BL_separation_points!(element::Element)

Computes the flow separation points of the modified Beddoes-Leishman model

# Arguments
- element::Element
"""
function BL_separation_points!(element::Element)

    @unpack α1₀N,α1₀M,α1₀T,βσ1N,βσ1T,βσ2N,βS2Nlpr,βS2Tlpr,βS1Nu,βS1Mu,βS1Tu,βS1Nd,βS1Md,βS1Td,βS2Nu,βS2Mu,βS2Tu,βS2Nd,βS2Md,βS2Td,ξ,f₀N,f₀M,f₀T,fbN,fbM,fbT,S1N,S1M,S1T,S2N,S2M,S2T = element.aero.airfoil.separatedFlowParameters
    @unpack α = element.aero.flowAnglesAndRates
    @unpack αlag,RD = element.aero.BLstates
    @unpack R = element.aero.BLkin
    @unpack stallOnsetRatio,upstroke,S,P,T,α1N,α1M,α1T = element.aero.BLflow
    @unpack minStallOnsetRatio,lastRD_tv0 = element.aero.BLcompVars

    # Factor defining the flattening of fPrime for light stall
    ψ = P^3*(1-T)

    # Factor defining how close to stalling angle has been at beginning of upstroke
    σ = minStallOnsetRatio*stallOnsetRatio > 0 ? min(0.5, abs(minStallOnsetRatio)^2.5) : 0.0 

    # In-and-out of stall factor
    ξRD_tv0 = ξ*((1-cos(2π*lastRD_tv0))/2)

    # Quasi-steady separation points: fN, fM, fT
    absα = abs(α)
    fN = absα <= α1₀N ? 1-(1-fbN)*exp((absα-α1₀N)/S1N) : f₀N+(fbN-f₀N)*exp((α1₀N-absα)/S2N)
    fM = absα <= α1₀M ? 1-(1-fbM)*exp((absα-α1₀M)/S1M) : f₀M+(fbM-f₀M)*exp((α1₀M-absα)/S2M)
    fT = absα <= α1₀T ? 1-(1-fbT)*exp((absα-α1₀T)/S1T) : f₀T+(fbT-f₀T)*exp((α1₀T-absα)/S2T)

    # Normal force unsteady lagged separation point based on lagged AoA: fPrimeN
    absαlag = abs(αlag)
    sqrtOfR = sqrt(R)
    if absαlag <= α1N
        αDiffN = absαlag - α1N
        S1PrimeN = upstroke ? S1N*(1+βS1Nu*RD) : S1N*(1+βS1Nd*RD)
        fPrimeN = 1-(1-fbN)*exp(αDiffN/S1PrimeN)
    else
        αDiffN = α1N-absαlag
        if upstroke  
            S2PrimeN = S2N*(1+βS2Nu*RD+βS2Nlpr*(1-RD)^2*(1-T)+βσ2N*σ*RD*(1-RD))
            f₀N2 = f₀N
        else                
            f₀N2 = !S ? fbN-0.2 : f₀N
            f₀N2 = min(0.25, (1-sqrtOfR)*(f₀N2+ξRD_tv0*σ))
            S2PrimeN = S2N*(1+βS2Nd*RD+βS2Nlpr*(1-RD)^2*(1-T)+βσ2N*σ*RD*(1-RD)+βσ1N*ψ*R)      
        end
        fPrimeN = f₀N2+(fbN-f₀N2)*exp(αDiffN/S2PrimeN)
    end

    # Pitching moment unsteady lagged separation point based on lagged AoA: fPrimeM
    if absαlag <= α1M
        αDiffM = absαlag - α1M
        S1PrimeM = upstroke ? S1M*(1+βS1Mu*RD) : S1M*(1+βS1Md*RD)
        fPrimeM = 1-(1-fbM)*exp(αDiffM/S1PrimeM)
    else
        αDiffM = α1M - absαlag
        if upstroke  
            S2PrimeM = S2M*(1+βS2Mu*RD)
            f₀M2 = f₀M
        else                
            f₀M2 = !S ? fbM-0.2 : f₀M
            f₀M2 = min(0.25, (1-sqrtOfR)*(f₀M2+ξRD_tv0*σ))
            S2PrimeM = S2M*(1+βS2Md*RD)
        end
        fPrimeM = f₀M2+(fbM-f₀M2)*exp(αDiffM/S2PrimeM)
    end

    # Tangential force unsteady lagged separation point based on lagged AoA: fPrimeT
    if absαlag <= α1T
        αDiffT = absαlag - α1T
        S1PrimeT = upstroke ? S1T*(1+βS1Tu*RD) : S1T*(1+βS1Td*RD)
        fPrimeT = 1-(1-fbT)*exp(αDiffT/S1PrimeT)
    else
        αDiffT = α1T - absαlag
        if upstroke  
            S2PrimeT = S2T*(1+βS2Tu*RD+βS2Tlpr*(1-RD)^2*(1-T))
            f₀T2 = f₀T
        else
            f₀T2 = !S ? fbT-0.2 : f₀T
            f₀T2 = min(0.25, (1-sqrtOfR)*(f₀T2+ξRD_tv0*σ))
            S2PrimeT = S2T*(1+βS2Td*RD+βS2Tlpr*(1-RD)^2*(1-T)+βσ1T*ψ*R)
        end
        fPrimeT = f₀T2+(fbT-f₀T2)*exp(αDiffT/S2PrimeT)
    end

    @pack! element.aero.BLflow = fN,fM,fT,fPrimeN,fPrimeM,fPrimeT

end


"""
BL_quasi_steady_separation_points(element::Element,Uₜ,Uₙ)

Computes the quasi-steady flow separation points of the modified Beddoes-Leishman model for an element

# Arguments
- element::Element
- Uₜ
- Uₙ
"""
function BL_quasi_steady_separation_points(element::Element,Uₜ,Uₙ)

    # Check tangential velocity
    if isapprox(abs(Uₜ),0)
        fN,fM,fT = 1.0,1.0,1.0
        return fN,fM,fT
    end

    @unpack α1₀N,α1₀M,α1₀T,f₀N,f₀M,f₀T,fbN,fbM,fbT,S1N,S1M,S1T,S2N,S2M,S2T = element.aero.airfoil.separatedFlowParameters

    # Angle of attack  
    α = atan(-Uₙ/Uₜ)
    absα = abs(α)

    # Quasi-steady separation points
    fN = absα <= α1₀N ? 1-(1-fbN)*exp((absα-α1₀N)/S1N) : f₀N+(fbN-f₀N)*exp((α1₀N-absα)/S2N)

    fM = absα <= α1₀M ? 1-(1-fbM)*exp((absα-α1₀M)/S1M) : f₀M+(fbM-f₀M)*exp((α1₀M-absα)/S2M)

    fT = absα <= α1₀T ? 1-(1-fbT)*exp((absα-α1₀T)/S1T) : f₀T+(fbT-f₀T)*exp((α1₀T-absα)/S2T)

    return fN,fM,fT

end


"""
BL_stall_time!(problem::Problem,element::Element)

Computes variables at the time of stall onset for the modified Beddoes-Leishman model

# Arguments
- problem::Problem
- element::Element
"""
function BL_stall_time!(problem::Problem,element::Element)

    @unpack timeNow,Δt = problem
    @unpack gᵥ,Tv = element.aero.airfoil.separatedFlowParameters
    @unpack f2primeN,RD = element.aero.BLstates
    @unpack qR,R = element.aero.BLkin
    @unpack stallOnsetRatio,upstroke,fN = element.aero.BLflow
    @unpack stallOnsetRatioPrev,fDiff_tv0P,qR_tv0P,R_tv0P,RD_tv0P,upstroke_tv0P,fDiff_tv0P2,RD_tv0P2,upstroke_tv0P2,fDiff_tv0N,qR_tv0N,R_tv0N,RD_tv0N,upstroke_tv0N,fDiff_tv0N2,RD_tv0N2,upstroke_tv0N2,tv0P,tv0N = element.aero.BLcompVars

    # Positive AoA
    #---------------------------------------------------------------------------
    # Primary vortex
    if stallOnsetRatio >= 1 && stallOnsetRatioPrev < 1 
        # Linear interpolation for time of vortex shedding 
        tv0P = interpolate([stallOnsetRatioPrev; stallOnsetRatio], [timeNow-Δt; timeNow], 1)
        # Variables at time of vortex shedding
        fDiff_tv0P = max(0, f2primeN-fN)
        qR_tv0P = qR
        R_tv0P = R
        RD_tv0P = RD
        upstroke_tv0P = upstroke
        # Reset flag to allow secondary vortices
        fDiff_tv0P2 = -2.0
    end
    # Time since primary vortex shedding 
    τvP = timeNow - tv0P
    # Secondary vortices
    if τvP >= (2+gᵥ)*Tv && fDiff_tv0P2 == -2.0
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
        tv0N = interpolate([-stallOnsetRatioPrev; -stallOnsetRatio], [timeNow-Δt; timeNow], 1)
        # Variables at time of vortex shedding
        fDiff_tv0N = max(0, f2primeN-fN)
        qR_tv0N = qR
        R_tv0N = R
        RD_tv0N = RD
        upstroke_tv0N = upstroke 
        # Reset flag to allow secondary vortices
        fDiff_tv0N2 = -2.0
    end
    # Time since primary vortex shedding 
    τvN = timeNow - tv0N
    # Secondary vortices
    if τvN >= (2+gᵥ)*Tv && fDiff_tv0N2 == -2.0
        # Set variables at the time of second vortex
        fDiff_tv0N2 = max(0, f2primeN-fN)
        RD_tv0N2 = RD
        upstroke_tv0N2 = upstroke
    end

    @pack! element.aero.BLcompVars = tv0P,fDiff_tv0P,qR_tv0P,R_tv0P,RD_tv0P,upstroke_tv0P,τvP,fDiff_tv0P2,RD_tv0P2,upstroke_tv0P2,tv0N,fDiff_tv0N,qR_tv0N,R_tv0N,RD_tv0N,upstroke_tv0N,τvN,fDiff_tv0N2,RD_tv0N2,upstroke_tv0N2

end


"""
BL_DSV_loads!(element::Element)

Computes DSV loads for the modified Beddoes-Leishman model

# Arguments
- element::Element
"""
function BL_DSV_loads!(element::Element)

    @unpack ϖMid = element.aero
    @unpack μv₂,ν₁,ν₂,ν₃,ν₄,ν₅,χu,χd,gᵥ,gᵥ₂,Tv,Tv₂,Vn₁,Vn₂,Vn₃,Vm,Vt = element.aero.airfoil.separatedFlowParameters
    @unpack qR,R = element.aero.BLkin
    @unpack stallOnsetRatio,upstroke,P = element.aero.BLflow
    @unpack fDiff_tv0P,qR_tv0P,R_tv0P,RD_tv0P,upstroke_tv0P,τvP,fDiff_tv0P2,RD_tv0P2,upstroke_tv0P2,fDiff_tv0N,qR_tv0N,R_tv0N,RD_tv0N,upstroke_tv0N,τvN,fDiff_tv0N2,RD_tv0N2,upstroke_tv0N2 = element.aero.BLcompVars

    # Positive AoA
    #---------------------------------------------------------------------------

    # First vortex: normal force and pitching moment
    #---------------------------------------------------------------------------
    if τvP <= 2*Tv 
        # Vortex strength 
        Vs₁ = Vn₁ * RD_tv0P^ν₁ * fDiff_tv0P^ν₂ * max(1, min(3, qR_tv0P^ν₃)) 
        # Vortex stretch factor
        χᵥ = upstroke_tv0P ? 1+χu*R_tv0P : 1+χd
        # Vortex shape function
        Vₓ = τvP/Tv <= 1/χᵥ ? sin(π/2*τvP/Tv*χᵥ)^(ν₄/χᵥ+RD_tv0P) : sin(π/2*(1+(χᵥ/(2*χᵥ-1))*(τvP/Tv-1/χᵥ)))^(ν₄*χᵥ+RD_tv0P)
        # Airloads coefficients
        cnV1P = Vs₁ * Vₓ
        cmV1P = -Vm * cnV1P
    else
        cnV1P = 0.0
        cmV1P = 0.0
    end

    # First vortex: tangential force 
    #---------------------------------------------------------------------------
    # Time at which vortex begins to affect tangential force
    tvT = Tv
    # Increment in vortex convection time at high pitch rates
    δTvc = min(1/2, qR_tv0P-R_tv0P)
    # Airload coefficient
    ccV1P = (τvP > Tv && τvP <= (2+δTvc)*Tv && upstroke_tv0P) ? -Vt * Vn₁ * RD_tv0P^ν₁ * fDiff_tv0P^ν₂ * min(5, qR_tv0P^ν₅) *((1-cos(2π*(τvP-tvT)/((1+δTvc)*Tv)))/2) : 0.0

    # Second vortex 
    #---------------------------------------------------------------------------
    # Time at which second vortex begins
    tv2 = (2+gᵥ)*Tv
    if τvP > tv2 && τvP <= tv2+2*Tv₂ && fDiff_tv0P2 > 0.05
        # Vortex strength
        Vs₂ = Vn₂ * RD_tv0P^ν₁ * fDiff_tv0P^ν₂ * max(1, min(3, qR_tv0P^ν₃)) * (1 + μv₂ * R_tv0P^2 * upstroke_tv0P2)
        # Vortex shape function
        Vx₂ = (1-cos(2π*(τvP-tv2)/(2*Tv₂)))/2 
        # Airloads coefficients
        cnV2P = Vs₂ * Vx₂
        cmV2P = -Vm * cnV2P * RD_tv0P2^3
    else
        cnV2P = 0.0
        cmV2P = 0.0
    end

    # Third vortex 
    #---------------------------------------------------------------------------
    # Time at which third vortex begins
    tv3 = tv2+(2+gᵥ₂)*Tv₂
    if τvP > tv3 && τvP <= tv3+2*Tv₂ && fDiff_tv0P2 > 0.1 && upstroke_tv0P2
        # Vortex strength
        Vs3 = Vn₃ * RD_tv0P2 * RD_tv0P^ν₁ * fDiff_tv0P^ν₂ * (1 + μv₂ * R_tv0P^2)
        # Vortex shape function 
        Vx₃ = (1-cos(2π*(τvP-tv3)/(2*Tv₂)))/2
        # Airloads coefficients
        cnV3P = Vs3 * Vx₃ 
        cmV3P = -Vm * cnV3P * RD_tv0P2^3
    else
        cnV3P = 0.0
        cmV3P = 0.0
    end 

    # Negative AoA
    #---------------------------------------------------------------------------

    # First vortex: normal force and pitching moment
    #---------------------------------------------------------------------------
    if τvN <= 2*Tv 
        # Vortex strength
        Vs₁ = Vn₁ * RD_tv0N^ν₁ * fDiff_tv0N^ν₂ * max(1, min(3, qR_tv0N^ν₃))
        # Vortex stretch factor
        χᵥ = upstroke_tv0N ? 1+χu*R_tv0N : 1+χd
        # Vortex shape function
        Vₓ = τvN/Tv <= 1/χᵥ ? sin(π/2*τvN/Tv*χᵥ)^(ν₄/χᵥ+RD_tv0N) : sin(π/2*(1+(χᵥ/(2*χᵥ-1))*(τvN/Tv-1/χᵥ)))^(ν₄*χᵥ+RD_tv0N)
        # Airloads coefficients
        cnV1N = -Vs₁ * Vₓ
        cmV1N = -Vm * cnV1N
    else
        cnV1N = 0.0
        cmV1N = 0.0
    end

    # First vortex: tangential force 
    #---------------------------------------------------------------------------
    # Time at which vortex begins to affect tangential force
    tvT = Tv
    # Increment in vortex convection time at high pitch rates
    δTvc = min(1/2, qR_tv0N-R_tv0N)
    # Airload coefficient
    ctV1N = (τvN > Tv && τvN <= (2+δTvc)*Tv && upstroke_tv0N) ? -Vt * Vn₁ * RD_tv0N^ν₁ * fDiff_tv0N^ν₂ * min(5, qR_tv0N^ν₅) *((1-cos(2π*(τvN-tvT)/((1+δTvc)*Tv)))/2) : 0.0

    # Second vortex 
    #---------------------------------------------------------------------------
    # Time at which second vortex begins 
    tv2 = (2+gᵥ)*Tv
    if τvN > tv2 && τvN <= tv2+2*Tv₂ && fDiff_tv0N2 > 0.05
        # Vortex strength
        Vs₂ = Vn₂ * RD_tv0N^ν₁ * fDiff_tv0N^ν₂ * max(1, min(3, qR_tv0N^ν₃)) * (1 + μv₂ * R_tv0N^2 * upstroke_tv0N2)
        # Vortex shape function 
        Vx₂ = (1-cos(2π*(τvN-tv2)/(2*Tv₂)))/2
        # Airloads coefficients
        cnV2N = -Vs₂ * Vx₂
        cmV2N = -Vm * cnV2N * RD_tv0N2^3
    else
        cnV2N = 0.0
        cmV2N = 0.0
    end

    # Third vortex 
    #---------------------------------------------------------------------------
    # Time at which third vortex begins
    tv3 = tv2+(2+gᵥ₂)*Tv₂
    if τvN > tv3 && τvN <= tv3+2*Tv₂ && fDiff_tv0N2 > 0.1 && upstroke_tv0N2
        # Vortex strength
        Vs3 = Vn₃ * RD_tv0N2 * RD_tv0N^ν₁ * fDiff_tv0N^ν₂ * (1 + μv₂*R_tv0N^2)
        # Vortex shape function 
        Vx₃ = (1-cos(2π*(τvN-tv3)/(2*Tv₂)))/2 
        # Airloads coefficients
        cnV3N = -Vs3 * Vx₃ 
        cmV3N = -Vm * cnV3N * RD_tv0N2^3
    else
        cnV3N = 0.0
        cmV3N = 0.0
    end

    # Total DSV loads
    #---------------------------------------------------------------------------
    cnV = cnV1P + cnV2P + cnV3P + cnV1N + cnV2N + cnV3N
    cmV = cmV1P + cmV2P + cmV3P + cmV1N + cmV2N + cmV3N
    ctV = ccV1P + ctV1N

    @pack! element.aero.aeroCoefficients = cnV,cmV,ctV

end


"""
BLi_cn!(element::Element,δNow)

Computes normal force coefficient for the modified incompressible Beddoes-Leishman model

# Arguments
- element::Element
- δNow
"""
function BLi_cn!(element::Element,δNow)

    @unpack flapLoadsSolver,flapped,b,δdotNow,δddotNow,ϖMid = element.aero
    @unpack αₑ = element.aero.flowAnglesAndRates
    @unpack Uᵢ,UₙdotMid = element.aero.flowVelocitiesAndRates
    @unpack ϵₙ,cnα,cnδ = element.aero.airfoil.separatedFlowParameters
    @unpack f2primeN = element.aero.BLstates
    @unpack cnV = element.aero.aeroCoefficients

    # Circulatory component - attached flow
    cnC = cnα * sin(αₑ)
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cnC += cnδ*δNow
    end

    # Circulatory component - separated flow
    cnF = cnC * ((1+sqrt(f2primeN))/2)^2

    # Inertial component
    cnI = π*b*UₙdotMid/Uᵢ^2
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cnI -= ϵₙ*b/Uᵢ^2*(Uᵢ*Th[4]*δdotNow+b*Th[1]*δddotNow)
    end

    # Total 
    cn = cnF + cnI + cnV

    # Scale by tip loss correction factor at element's midpoint
    cn,cnC,cnI,cnF,cnV = multiply_inplace!(ϖMid,cn,cnC,cnI,cnF,cnV)

    @pack! element.aero.aeroCoefficients = cn,cnC,cnI,cnF,cnV

end


"""
BLi_cm!(element::Element,δNow)

Computes the pitching moment aerodynamic coefficient at the attachment point (i.e., the beam reference line) for the modified incompressible Beddoes-Leishman model

# Arguments
- element::Element
- δNow
"""
function BLi_cm!(element::Element,δNow)

    @unpack flapLoadsSolver,flapped,b,normSparPos,normFlapPos,δdotNow,δddotNow,ϖMid = element.aero
    @unpack αₑ = element.aero.flowAnglesAndRates
    @unpack cnF = element.aero.aeroCoefficients
    @unpack Uᵢ,UₙdotMid,Ωₐ,Ωₐdot = element.aero.flowVelocitiesAndRates
    @unpack ϵₘ,κ₀,κ₁,κ₂,κ₃,cm₀,cmδ,cnα,K₀,K₁,K₂ = element.aero.airfoil.separatedFlowParameters
    @unpack R = element.aero.BLkin
    @unpack stallOnsetRatio,upstroke,S,P = element.aero.BLflow
    @unpack f2primeM,RD = element.aero.BLstates
    @unpack cmV = element.aero.aeroCoefficients

    # Center of pressure dynamic variables
    K1Prime = K₁*(1-κ₁*RD*(1-abs(stallOnsetRatio)))-κ₂*R*abs(stallOnsetRatio)*upstroke
    K2Prime = K₂*(1+κ₃*S*R^2*(1-P)*(!upstroke))                        

    # Center of pressure offset from quarter-chord
    δCP = K₀ + K1Prime*(1-f2primeM) + K2Prime*sin(π*f2primeM^κ₀)

    # Circulatory component - separated flow
    cmF = cm₀ + cnF/ϖMid*(normSparPos-(1/4-δCP))
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cmF += cmδ*δNow
    end

    # Inertial component
    cmI = -π*b/(2*Uᵢ^2)*((1-2*normSparPos)*UₙdotMid+b*Ωₐdot/8)
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cmI -= ϵₘ/(2*Uᵢ^2)*(Uᵢ^2*Th[14]*δNow+Uᵢ*b*Th[15]*δdotNow+b^2*Th[16]*δddotNow)
    end

    # Rotation-induced component (this is the increment due to the fact that UₙdotMid acts at midchord but the normal force is at the 3/4-chord)
    cmRot = -π/4*b*Ωₐ/Uᵢ

    # Total at attachment point
    cm = cmF + cmI + cmRot + cmV

    # Scale by tip loss correction factor at element's midpoint
    cm,cmF,cmI,cmRot,cmV = multiply_inplace!(ϖMid,cm,cmF,cmI,cmRot,cmV)

    @pack! element.aero.aeroCoefficients = cm,cmF,cmI,cmRot

end


"""
BLi_ct!(element::Element,δNow)

Computes the tangential force aerodynamic coefficient for the modified incompressible Beddoes-Leishman model

# Arguments
- element::Element
- δNow
"""
function BLi_ct!(element::Element,δNow)

    @unpack flapped,flapLoadsSolver,ϖMid = element.aero
    @unpack αₑ = element.aero.flowAnglesAndRates
    @unpack α1₀T,η,cd₀,cdδ,cnα,E₀,E₁ = element.aero.airfoil.separatedFlowParameters
    @unpack R = element.aero.BLkin
    @unpack stallOnsetRatio,upstroke,S = element.aero.BLflow
    @unpack αlag,f2primeT,RD = element.aero.BLstates
    @unpack ctV = element.aero.aeroCoefficients
    @unpack Ts = element.aero.BLcompVars

    # Ratio of lagged AoA to tangential steady breakpoint angle
    stallOnsetRatioT = abs(αlag)/α1₀T

    # Circulatory component - separated flow
    ctF = -cd₀/cos(αₑ) + (1-η*RD^2) * cnα * (αₑ)^2 * f2primeT^(1/2+ stallOnsetRatioT + E₀*RD*S*(1-Ts)*R*abs(stallOnsetRatio)^(1/2)*(!upstroke))
    if stallOnsetRatioT > 1
        ctF -= E₁ * (1-RD^3) * min(1, stallOnsetRatioT^3-1)
    end
    if flapped && typeof(flapLoadsSolver) == TableLookup
        ctF -= cdδ*abs(δNow)/cos(αₑ)
    end

    # Total
    ct = ctF + ctV

    # Scale by tip loss correction factor at element's midpoint
    ct,ctF,ctV = multiply_inplace!(ϖMid,ct,ctF,ctV)

    @pack! element.aero.aeroCoefficients = ct,ctF,ctV

end


"""
BLi_state_matrices!(element::Element,δNow)

Computes the aerodynamic state matrices for the modified incompressible Beddoes-Leishman model

# Arguments
- element::Element
- δNow
"""
function BLi_state_matrices!(element::Element,δNow)

    @unpack nTotalAeroStates,linearPitchPlungeStatesRange,nonlinearPitchPlungeStatesRange,flapStatesRange,linearGustStatesRange,nonlinearGustStatesRange,b = element.aero
    @unpack βₚ² = element.aero.flowParameters
    @unpack Ta,λbWMat = element.aero.airfoil.separatedFlowParameters
    @unpack Uₙ,Uᵢ,UₙdotTQC = element.aero.flowVelocitiesAndRates
    @unpack R = element.aero.BLkin
    @unpack Ta_SO,TfN,TfM,TfT,fPrimeN,fPrimeM,fPrimeT = element.aero.BLflow
    @unpack AW,bWMat = element.aero.solver
    @unpack cn = element.aero.aeroCoefficients

    # Initialize state matrices with appropriate types
    A = zeros(typeof(cn),nTotalAeroStates,nTotalAeroStates)
    B = zeros(typeof(cn),nTotalAeroStates)

    tmp = -Uᵢ/b*βₚ²
        
    # Get the rate of cnαUₙTQC
    cnαUₙTQCdot = cnαUₙTQC_rate(element)

    # Set entries of nonlinear states' matrices
    tmpA = Vector{typeof(cn)}(undef, 6)
    tmpB = Vector{typeof(cn)}(undef, 6)
    tmpA[1] = -1/Ta
    tmpB[1] = Uₙ/Ta
    tmpA[2] = -1/TfN
    tmpB[2] = fPrimeN/TfN
    tmpA[3] = -1/TfM
    tmpB[3] = fPrimeM/TfM
    tmpA[4] = -1/TfT
    tmpB[4] = fPrimeT/TfT
    tmpA[5] = -1/(3*Ta)
    tmpB[5] = R/(3*Ta)
    tmpA[6] = -1/Ta_SO
    tmpB[6] = R/Ta_SO
    
    # Pitch-plunge-induced flow state matrices
    A[linearPitchPlungeStatesRange,linearPitchPlungeStatesRange] .= tmp*bWMat.*λbWMat
    A[nonlinearPitchPlungeStatesRange,nonlinearPitchPlungeStatesRange] .= Diagonal(tmpA)
    B[linearPitchPlungeStatesRange] .= cnαUₙTQCdot*AW
    B[nonlinearPitchPlungeStatesRange] .= tmpB

    # Flap-induced flow state matrices
    if !isnothing(flapStatesRange)
        @unpack AWf,bWfMat = element.aero.flapLoadsSolver
        # Get flap normalwash rate
        wdotFlap = flap_normalwash_rate(element,δNow)
        # Set state matrices
        A[flapStatesRange,flapStatesRange] .= tmp*bWfMat
        B[flapStatesRange] .= wdotFlap*AWf
    end

    # Gust-induced flow state matrices
    if !isnothing(linearGustStatesRange)
        @unpack bGMat = element.aero.gustLoadsSolver
        @unpack UₙGust = element.aero.flowVelocitiesAndRates
        @unpack Tg = element.aero.airfoil.separatedFlowParameters
        # Linear part 
        A[linearGustStatesRange,linearGustStatesRange] .= tmp*bGMat
        B[linearGustStatesRange] .= UₙGust*ones(length(linearGustStatesRange))
        # Nonlinear part
        A[nonlinearGustStatesRange,nonlinearGustStatesRange] .= -1/Tg
        B[nonlinearGustStatesRange] .= UₙGust/Tg
    end

    @pack! element.aero = A,B

    return A,B
end


"""
update_initial_aero_states!(problem::Problem)

Updates the initial aerodynamic states assuming their rates are zero

# Arguments
- problem::Problem
"""
function update_initial_aero_states!(problem::Problem)

    @unpack x,model = problem

    # Loop over elements
    for element in model.elements

        @unpack DOF_χ = element

        # Skip elements without aerodynamic states
        if isempty(DOF_χ)
            continue 
        end

        @unpack a = model.atmosphere
        @unpack V,Ω,χ = element.states
        @unpack Vdot,Ωdot = element.statesRates
        @unpack solver,pitchPlungeStatesRange,flapStatesRange,pitchPlungeStatesRange,RwT,c,b,normSparPos = element.aero

        # Relative wind velocity and acceleration vectors at spar position, resolved in basis W
        U = -RwT*V
        Udot = -RwT*Vdot

        # Spanwise angular velocity and acceleration components, resolved in basis W
        Ωₐ = dot(RwT[1,:],Ω)
        Ωₐdot = dot(RwT[1,:],Ωdot)

        # Tangential and normal components of relative wind velocity and acceleration
        Uₜ,Uₙ = U[2],U[3]
        Uₜdot,Uₙdot = Udot[2],Udot[3]

        # In-plane relative wind velocity and acceleration 
        Uᵢ = sqrt(Uₜ^2+Uₙ^2)
        Uᵢdot = (Uₜ*Uₜdot+Uₙ*Uₙdot)/Uᵢ

        # Relative normal wind velocity and acceleration components at the airfoil's 1/2- and 3/4-chord
        UₙTQC = Uₙ - Ωₐ*c*(normSparPos-3/4)
        UₙdotTQC = Uₙdot - Ωₐdot*c*(normSparPos-3/4)

        # Mach, compressibility factor and its rate
        Ma = Uᵢ/a
        βₚ² = 1-Ma^2
        βₚ = sqrt(βₚ²)
        βₚdot = -Uᵢ*Uᵢdot/(a^2*βₚ)

        # Pitch-plunge states
        if typeof(solver) in [Indicial]
            @unpack AW,bW = solver
            @unpack cnα = element.aero.airfoil.attachedFlowParameters
            # Time derivative of cnα (assuming it scales with 1/βₚ)
            cnαdot = cnα*βₚdot/βₚ²
            # Time derivative of the product cnα * UₙTQC
            cnαUₙTQCdot = cnα*UₙdotTQC+cnαdot*UₙTQC
            # States
            χ[pitchPlungeStatesRange] = cnαUₙTQCdot*b*βₚ²/Uᵢ*AW./bW
        elseif typeof(solver) in [BLi]
            @unpack AW,bWMat = solver
            @unpack Ta,cnα,λbWMat = element.aero.airfoil.separatedFlowParameters
            @unpack R = element.aero.BLkin
            @unpack Ta_SO,TfN,TfM,TfT,fPrimeN,fPrimeM,fPrimeT = element.aero.BLflow
            # Time derivative of cnα (assuming it scales with 1/βₚ)
            cnαdot = cnα*βₚdot/βₚ²
            # Time derivative of the product cnα * UₙTQC
            cnαUₙTQCdot = cnα*UₙdotTQC+cnαdot*UₙTQC
            # State matrices
            A = -diagm([Uᵢ/b*βₚ²*bWMat[1,1]*λbWMat[1,1]; Uᵢ/b*βₚ²*bWMat[2,2]*λbWMat[2,2]; 1/Ta; 1/TfN; 1/TfM; 1/TfT; 1/(3*Ta); 1/Ta_SO])
            B = [cnαUₙTQCdot*AW[1]; cnαUₙTQCdot*AW[2]; Uₙ/Ta; fPrimeN/TfN; fPrimeM/TfM; fPrimeT/TfT; R/(3*Ta); R/Ta_SO]
            # States
            χ[pitchPlungeStatesRange] = -A\B
            # Overwrite separation points with quasi-steady values
            fN,fM,fT = BL_quasi_steady_separation_points(element,Uₜ,Uₙ)
            χ[pitchPlungeStatesRange[4:6]] = [fN; fM; fT]   
        elseif typeof(solver) == Inflow
            @unpack AₚInv,AₚInvcₚ = element.aero.solver
            # State matrices
            A = -Uᵢ/b*βₚ²*AₚInv
            B = UₙdotTQC*AₚInvcₚ
            # States
            χ[pitchPlungeStatesRange] = -A\B
        end

        # Flap-induced flow states
        if !isnothing(flapStatesRange)
            @unpack flapLoadsSolver,δNow,δdotNow,δddotNow = element.aero
            @unpack Th,AWf,bWf = flapLoadsSolver
            # Flap normalwash rate
            wdotFlap = Th[10]/π*(Uᵢ*δdotNow+Uᵢdot*δNow)+b*Th[11]/(2π)*δddotNow
            # Flap states
            χ[flapStatesRange] = wdotFlap*b/Uᵢ*AWf./bWf
        end

        # Update states
        x[DOF_χ] = χ
    end

    @pack! problem = x
end