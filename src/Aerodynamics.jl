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

    @pack! element.aero = Re,Ma,βₚ,βₚ²

end


"""
local_gust_velocity!(model::Model,element::Element,timeNow)

Computes the gust velocity in the local, deformed aerodynamic basis (basis W)

# Arguments
- model::Model
- element::Element
- timeNow
"""
function local_gust_velocity!(model::Model,element::Element,timeNow)

    @unpack R_AT = model
    @unpack R = element
    @unpack RwR0,UgustInertial = element.aero

    # Transform gust velocity vector from inertial basis to local, deformed aerodynamic basis 
    UgustLocal = (R*RwR0)'*R_AT*UgustInertial(timeNow)

    @pack! element.aero = UgustLocal

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
aero_coefficients!(element::Element,δNow)

Computes the aerodynamic coefficients

# Arguments
- element::Element
- δNow
"""
function aero_coefficients!(element::Element,δNow)

    @unpack solver = element.aero

    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        attached_flow_aero_coefficients!(element,δNow)
    elseif typeof(solver) == BLi
        BLi_aero_coefficients!(element)
    elseif typeof(solver) == BLc
        BLc_aero_coefficients!(element)    
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
        A,B = BLi_state_matrices!(element)
    elseif typeof(solver) == BLc
        A,B = BLc_state_matrices!(element)    
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
update_airfoil_parameters!(element::Element)

Updates the aerodynamic parameters of the airfoil according to current nondimensional flow parameters 

# Arguments
- element::Element
"""
function update_airfoil_parameters!(element::Element)

    @unpack airfoil,updateAirfoilParameters,Re,Ma,flapSiteID = element.aero
    @unpack name = airfoil

    # Skip if parameters are not to be updated
    if !updateAirfoilParameters
        return
    end
    
    # Get airfoil parameters
    attachedFlowParameters = AttachedFlowParameters(name,Re=Re,Ma=Ma,flapSiteID=flapSiteID)

    @pack! airfoil = attachedFlowParameters
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

    @unpack airfoil = element.aero
    @unpack Uₜ,UₜGust = element.aero.flowVelocitiesAndRates
    @unpack ϵₙ = airfoil.attachedFlowParameters

    # Effective normalwash: pitch-plunge-induced, flap-induced, gust-induced, and total
    wₑp = pitch_plunge_effective_normalwash(element,χ)
    wₑf = flap_effective_normalwash(element,δNow)
    wₑg = gust_effective_normalwash(element)
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
    cnC = cnα*αₑ
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
    cn = cnC+cnI

    # Scale by tip loss correction factor at element's midpoint
    cn,cnC,cnI = ϖMid*cn,ϖMid*cnC,ϖMid*cnI

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
    cmC = cm₀+cmα*αₑ+cnC*(normSparPos-1/4)
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cmC += cmδ*δNow
    end

    # Inertial component
    cmI = -π*b/(2*Uᵢ^2)*((1-2*normSparPos)*UₙdotMid+b*Ωₐdot/8)
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cmI -= ϵₘ/(2*Uᵢ^2)*(Uᵢ^2*(Th[4]+Th[10])*δNow+Uᵢ*b*(Th[1]-Th[8]-(normFlapPos-normSparPos)*Th[4]+Th[11]/2)*δdotNow+b^2*(Th[7]+Th[1]*(normFlapPos-normSparPos))*δddotNow)
    end

    # Rotation-induced component (this is the increment due to the fact that UₙdotMid acts at midchord but the normal force is at the 3/4-chord)
    cmRot = -π/4*b*Ωₐ/Uᵢ

    # Total at attachment point
    cm = cmC+cmI+cmRot

    # Scale by tip loss correction factor at element's midpoint
    cm,cmC,cmI,cmRot = ϖMid*cm,ϖMid*cmC,ϖMid*cmI,ϖMid*cmRot

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
    ct = -cd₀/cos(αₑ)+cnα*αₑ^2
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

    @unpack solver,pitchPlungeStatesRange = element.aero
    @unpack UₙTQC = element.aero.flowVelocitiesAndRates
    @unpack cnα = element.aero.airfoil.attachedFlowParameters

    if typeof(solver) == QuasiSteady
        wₑp = UₙTQC
    elseif typeof(solver) == Indicial
        wₑp = UₙTQC-sum(χ[pitchPlungeStatesRange])/cnα
    elseif typeof(solver) == Inflow
        @unpack bₚ = solver
        wₑp = UₙTQC-1/2*dot(bₚ,χ[pitchPlungeStatesRange])
    end

    return wₑp
end


"""
flap_effective_normalwash(element::Element,δNow)

Computes the effective (unsteady) flap-induced normalwash

# Arguments
- element::Element
- δNow
"""
function flap_effective_normalwash(element::Element,δNow)

    @unpack χ = element.states
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

    wFlap = Th[10]/π*Uᵢ*δNow+b*Th[11]/(2π)*δdotNow
  
    return wFlap
end


"""
gust_effective_normalwash(element::Element)

Computes the effective (unsteady) gust-induced normalwash

# Arguments
- element::Element
"""
function gust_effective_normalwash(element::Element)

    @unpack χ = element.states
    @unpack gustStatesRange = element.aero

    # Skip if there are not gust states
    if isnothing(gustStatesRange)
        return 0
    end

    @unpack b,βₚ² = element.aero
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates
    @unpack AGbG = element.aero.gustLoadsSolver

    # Effective gust-induced normalwash
    wₑg = Uᵢ/b*βₚ²*dot(AGbG,χ[gustStatesRange])

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

    wdotFlap = Th[10]/π*(Uᵢ*δdotNow+Uᵢdot*δNow)+b*Th[11]/(2π)*δddotNow

    return wdotFlap
end


"""
cnαUₙTQC_rate(element::Element)

Computes the time rate of the product of cnα by UₙTQC

# Arguments
- element::Element
"""
function cnαUₙTQC_rate(element::Element)

    @unpack Ma,βₚ,βₚ² = element.aero
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

    @unpack solver,nTotalAeroStates,pitchPlungeStatesRange,flapStatesRange,gustStatesRange,b,βₚ² = element.aero
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
        @unpack AW,bWDiag = element.aero.solver
        # Get the rate of cnαUₙTQC
        cnαUₙTQCdot = cnαUₙTQC_rate(element)
        # Set state matrices
        A[pitchPlungeStatesRange,pitchPlungeStatesRange] = tmp*bWDiag
        B[pitchPlungeStatesRange] = cnαUₙTQCdot*AW
    elseif typeof(solver) == Inflow
        @unpack AₚInv,AₚInvcₚ = element.aero.solver
        # Set state matrices
        A[pitchPlungeStatesRange,pitchPlungeStatesRange] = tmp*AₚInv
        B[pitchPlungeStatesRange] = UₙdotTQC*AₚInvcₚ
    end

    # Flap-induced flow states
    if !isnothing(flapStatesRange)
        @unpack AWf,bWfDiag = element.aero.flapLoadsSolver
        # Get flap normalwash rate
        wdotFlap = flap_normalwash_rate(element,δNow)
        # Set state matrices
        A[flapStatesRange,flapStatesRange] = tmp*bWfDiag
        B[flapStatesRange] = wdotFlap*AWf
    end

    # Gust-induced flow states
    if !isnothing(gustStatesRange)
        @unpack bGDiag = element.aero.gustLoadsSolver
        @unpack UₙGust = element.aero
        # Set state matrices
        A[gustStatesRange,gustStatesRange] = tmp*bGDiag
        B[gustStatesRange] = UₙGust*ones(length(gustStatesRange))
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
        @unpack solver,pitchPlungeStatesRange,flapStatesRange,RwT,c,b,normSparPos = element.aero

        # Relative wind velocity and acceleration vectors at spar position, resolved in basis W
        U = -RwT*V
        Udot = -RwT*Vdot

        # Spanwise angular velocity and acceleration components, resolved in basis W
        Ωₐ = dot(RwT[1,:],Ω)
        Ωₐdot = dot(RwT[1,:],Ωdot)

        # Tangential and normal components of relative wind velocity and acceleration
        Uₜ,Uₙ = U[1],U[2],U[3]
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
        if typeof(solver) == Indicial
            @unpack AW,bW = solver
            @unpack cnα = element.aero.airfoil.attachedFlowParameters
            # Time derivative of cnα (assuming it scales with 1/βₚ)
            cnαdot = cnα*βₚdot/βₚ²
            # Time derivative of the product cnα * UₙTQC
            cnαUₙTQCdot = cnα*UₙdotTQC+cnαdot*UₙTQC
            # States
            χ[pitchPlungeStatesRange] = cnαUₙTQCdot*b/Uᵢ*AW./bW
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