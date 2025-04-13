# Computes the steady aerodynamic kinematic variables
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


# Computes the unsteady aerodynamic kinematic variables
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


# Computes the nondimensional flow parameters
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

    # Flow inverse time scale (this is not nondimensional)
    Θ = Uᵢ/(c/2)*βₚ²

    # Chord-to-sound-speed time scale (this is not nondimensional)
    Tᵢ = c/a

    @pack! element.aero.flowParameters = Re,Ma,βₚ,βₚ²,Θ,Tᵢ

end


# Computes the gust velocity in the local, deformed aerodynamic basis (basis W)
function local_gust_velocity!(problem::Problem,model::Model,element::Element)

    @unpack timeNow = problem
    @unpack R_A,R_AT,gust = model
    @unpack R = element
    @unpack RwR0 = element.aero

    # Reset gust velocity vector
    UGust = zeros(3)

    # Gust defined over time
    if gust.isDefinedOverTime && (gust.initialTime <= timeNow <= gust.finalTime)
        @unpack UGustInertial = gust
        # Transform gust velocity vector from inertial basis to local deformed aerodynamic basis W
        UGust = (R*RwR0)'*R_AT * UGustInertial(timeNow)
    # Gust defined over space    
    elseif !gust.isDefinedOverTime
        @unpack UGustInertial = gust
        @unpack u_A = model
        @unpack r = element
        @unpack u = element.states
        # Get current position of element's midpoint, resolved in the inertial basis
        r⁰ = u_A(timeNow) + R_A * (r + u)
        # Transform gust velocity vector from inertial basis to local deformed aerodynamic basis W
        UGust = (R*RwR0)'*R_AT * UGustInertial(r⁰,timeNow)
    end

    # Set tangential and normal gust velocity components
    UₜGust,UₙGust = UGust[2],UGust[3]

    @pack! element.aero.flowVelocitiesAndRates = UₜGust,UₙGust

end


# Computes the current values of flap deflection rates
function flap_deflection_rates!(problem,element::Element)

    @unpack timeNow = problem
    @unpack δdot,δddot = element.aero
    
    δdotNow,δddotNow = δdot(timeNow),δddot(timeNow)

    @pack! element.aero = δdotNow,δddotNow

end


# Computes the aerodynamic coefficients
function aero_coefficients!(problem::Problem,element::Element,χ,δNow)

    @unpack solver = element.aero

    # Compute coefficients
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        attached_flow_aero_coefficients!(element,χ,δNow)
    elseif typeof(solver) == BLi
        BLi_aero_coefficients!(problem,element,χ,δNow)
    elseif typeof(solver) == BLo
        BLo_aero_coefficients!(problem,element,χ,δNow)    
    end

end


# Computes the aerodynamic state matrices
function aero_state_matrices!(element::Element,δNow)

    @unpack solver = element.aero

    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        A,B = attached_flow_state_matrices!(element,δNow)
    elseif typeof(solver) == BLi
        A,B = BLi_state_matrices!(element,δNow)
    elseif typeof(solver) == BLo
        A,B = BLo_state_matrices!(element,δNow)    
    end   
    
    return A,B

end


# Computes the aerodynamic nodal loads resultants
function aero_loads_resultants!(model::Model,element::Element)

    @unpack ρ = model.atmosphere
    @unpack Δℓ,R = element
    @unpack solver,c,RwR0,cosΛ,ϖ,ϖMid,hasTipCorrection,Λ = element.aero
    @unpack Uᵢ,Uᵢdot = element.aero.flowVelocitiesAndRates
    @unpack ct,cn,cm,ctNC,cnNC,cmNC = element.aero.aeroCoefficients

    # Aerodynamic loads per unit length, as a function of the local element coordinate (ζ)
    f(ζ) = 1/2*ρ*Uᵢ^2*(c/cosΛ)*Δℓ * (ϖ(ζ)/ϖMid * cosΛ^2 * [(ct-ctNC); (cn-cnNC); (c/cosΛ)*(cm-cmNC)] .+ [ctNC; cnNC; (c/cosΛ)*cmNC])

    # Shape function over element domain
    ϕ(n,ζ) = ifelse.(n==1, 1-ζ, ζ)

    # Loop over nodes - Set aerodynamic loads array, resolved in basis W
    T = promote_type(typeof(Uᵢ), typeof(Uᵢdot))
    F = zeros(T,12)
    for node=1:2  
        # Nodal loads array: perform integration only if tip correction is present, otherwise split equally among nodes
        if hasTipCorrection
            integrand(ζ) = f(ζ) .* ϕ(node,ζ) 
            F[6*node-4:6*node-2], = quadgk(integrand, 0, 1)
        else
            F[6*node-4:6*node-2] = f(1/2)/2
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


# Updates the aerodynamic parameters of the airfoil according to current nondimensional flow parameters 
function update_airfoil_parameters!(element::Element)

    @unpack solver,airfoil,updateAirfoilParameters,flapSiteID,b = element.aero
    @unpack name = airfoil
    @unpack Re,Ma = element.aero.flowParameters
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates

    # Skip if parameters are not to be updated
    if !updateAirfoilParameters
        return
    end
    
    # Get airfoil parameters
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        attachedFlowParameters = AttachedFlowParameters(name,Re=Re,Ma=Ma)
        @pack! airfoil = attachedFlowParameters
    elseif typeof(solver) == BLi
        parametersBLi = BLiParameters(name,Re=Re,Ma=Ma,U=Uᵢ,b=b)
        @pack! airfoil = parametersBLi
    elseif typeof(solver) == BLo
        parametersBLo = BLoParameters(name,Re=Re,Ma=Ma,U=Uᵢ,b=b)
        @pack! airfoil = parametersBLo  
    end

    @pack! element.aero = airfoil

end


# Computes the effective (unsteady) angle of attack
function effective_angle_of_attack!(element::Element,χ,δNow)

    @unpack solver,airfoil = element.aero
    @unpack Uₜ,UₜGust = element.aero.flowVelocitiesAndRates
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        @unpack α₀N,ϵₙ = airfoil.attachedFlowParameters
    elseif typeof(solver) in [BLi]
        @unpack α₀N,ϵₙ = airfoil.parametersBLi
    elseif typeof(solver) in [BLo]
        @unpack α₀N,ϵₙ = airfoil.parametersBLo    
    end

    # Effective normalwash: pitch-plunge-induced, flap-induced, gust-induced, and total
    wₑp = pitch_plunge_effective_normalwash(element,χ)
    wₑf = flap_effective_normalwash(element,χ,δNow)
    wₑg = gust_effective_normalwash(element,χ)
    wₑ = wₑp + ϵₙ*wₑf + wₑg

    # Effective circulatory AoA 
    αₑ = atan(-wₑ/(Uₜ+UₜGust)) - α₀N

    @pack! element.aero.flowAnglesAndRates = αₑ
    @pack! element.aero.flowVelocitiesAndRates = wₑp

end


# Computes the effective (unsteady) pitch-plunge-induced normalwash
function pitch_plunge_effective_normalwash(element::Element,χ)

    @unpack solver,circulatoryPitchPlungeStatesRange = element.aero
    @unpack UₙTQC = element.aero.flowVelocitiesAndRates
    @unpack Θ = element.aero.flowParameters
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        @unpack cnα = element.aero.airfoil.attachedFlowParameters
    elseif typeof(solver) in [BLi]
        @unpack cnα = element.aero.airfoil.parametersBLi
    elseif typeof(solver) == BLo
        @unpack cnα = element.aero.airfoil.parametersBLo    
    end

    # Effective pitch-plunge-induced normalwash
    if typeof(solver) == QuasiSteady
        wₑp = UₙTQC
    elseif typeof(solver) == Indicial
        wₑp = UₙTQC-sum(χ[circulatoryPitchPlungeStatesRange])/cnα
    elseif typeof(solver) == Inflow
        @unpack bₚ = solver
        wₑp = UₙTQC-1/2*dot(bₚ,χ[circulatoryPitchPlungeStatesRange])
    elseif typeof(solver) in [BLi,BLo]
        wₑp = UₙTQC-sum(χ[circulatoryPitchPlungeStatesRange])/cnα
    end

    return wₑp
end


# Computes the effective (unsteady) flap-induced normalwash
function flap_effective_normalwash(element::Element,χ,δNow)

    @unpack solver,flapStatesRange = element.aero
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        @unpack cnα = element.aero.airfoil.attachedFlowParameters
    elseif typeof(solver) in [BLi]
        @unpack cnα = element.aero.airfoil.parametersBLi
    elseif typeof(solver) == BLo
        @unpack cnα = element.aero.airfoil.parametersBLo    
    end

    # Skip if there are no flap states
    if isnothing(flapStatesRange)
        return 0
    end
    
    # Quasi-steady flap-induced normalwash
    wFlap = flap_normalwash(element,δNow)

    # Effective flap-induced normalwash
    wₑf = wFlap-sum(χ[flapStatesRange])/cnα

    return wₑf
end


# Computes the instantaneous (quasi-steady) flap-induced normalwash
function flap_normalwash(element::Element,δNow)

    @unpack b,δdotNow = element.aero
    @unpack Th = element.aero.flapLoadsSolver
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates

    wFlap = Th[17]*Uᵢ*δNow+b*Th[18]*δdotNow
  
    return wFlap
end


# Computes the time rate of the instantaneous flap-induced normalwash
function flap_normalwash_rate(element::Element,δNow)

    @unpack b,δdotNow,δddotNow = element.aero
    @unpack Th = element.aero.flapLoadsSolver
    @unpack Uᵢ,Uᵢdot = element.aero.flowVelocitiesAndRates

    wdotFlap = Th[17]*(Uᵢ*δdotNow+Uᵢdot*δNow)+b*Th[18]*δddotNow

    return wdotFlap
end


# Computes the effective (unsteady) gust-induced normalwash
function gust_effective_normalwash(element::Element,χ)

    @unpack linearGustStatesRange = element.aero

    # Skip if there are not gust states
    if isnothing(linearGustStatesRange)
        return 0
    end

    @unpack Θ = element.aero.flowParameters
    @unpack AGbG = element.aero.gustLoadsSolver

    # Effective gust-induced normalwash
    wₑg = Θ*dot(AGbG,χ[linearGustStatesRange])

    return wₑg
end


# Computes the time rate of the normal force coefficient slope, cnα
function cnα_rate(element::Element)

    @unpack solver = element.aero
    @unpack Ma,βₚ,βₚ² = element.aero.flowParameters
    @unpack Uᵢ,Uᵢdot = element.aero.flowVelocitiesAndRates
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        @unpack cnα = element.aero.airfoil.attachedFlowParameters
    elseif typeof(solver) in [BLi]
        @unpack cnα = element.aero.airfoil.parametersBLi
    elseif typeof(solver) == BLo
        @unpack cnα = element.aero.airfoil.parametersBLo    
    end

    # Time derivative of cnα (assuming it scales with 1/βₚ)
    βₚdot = -Uᵢ*Uᵢdot/((Uᵢ/Ma)^2*βₚ)
    cnαdot = cnα*βₚdot/βₚ^2

    return cnαdot
end


# Computes the time rate of the product of cnα by UₙTQC
function cnαUₙTQC_rate(element::Element)

    @unpack solver = element.aero
    @unpack UₙTQC,UₙdotTQC = element.aero.flowVelocitiesAndRates
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        @unpack cnα = element.aero.airfoil.attachedFlowParameters
    elseif typeof(solver) in [BLi]
        @unpack cnα = element.aero.airfoil.parametersBLi
    elseif typeof(solver) == BLo
        @unpack cnα = element.aero.airfoil.parametersBLo    
    end

    # Time derivative of cnα
    cnαdot = cnα_rate(element)

    # Time derivative of the product cnα * UₙTQC
    return cnα*UₙdotTQC+cnαdot*UₙTQC
end


# Computes the time rate of the product of cnα by wFlap
function cnαwFlap_rate(element::Element,δNow)

    @unpack solver = element.aero
    if typeof(solver) in [QuasiSteady,Indicial,Inflow]
        @unpack cnα = element.aero.airfoil.attachedFlowParameters
    elseif typeof(solver) in [BLi]
        @unpack cnα = element.aero.airfoil.parametersBLi
    elseif typeof(solver) == BLo
        @unpack cnα = element.aero.airfoil.parametersBLo    
    end

    # Quasi-steady flap-induced normalwash
    wFlap = flap_normalwash(element,δNow)

    # Time derivative of wFlap
    wdotFlap = flap_normalwash_rate(element,δNow)

    # Time derivative of cnα
    cnαdot = cnα_rate(element)

    # Time derivative of the product cnα * wFlap
    return cnα*wdotFlap+cnαdot*wFlap
end


# Computes the aerodynamic coefficients under attached flow conditions
function attached_flow_aero_coefficients!(element::Element,χ,δNow)

    # Update inertial parameters, if applicable
    if typeof(element.aero.solver) in [Indicial,Inflow] && !element.aero.solver.incompressibleInertialLoads
        update_inertial_parameters!(element)
    end

    # Normal force coefficient
    attached_flow_cn!(element,χ,δNow)

    # Pitching moment coefficient about the spar position
    attached_flow_cm!(element,χ,δNow)

    # Tangential force coefficient
    attached_flow_ct!(element,δNow)

end


# Computes the normal force aerodynamic coefficient for attached flow
function attached_flow_cn!(element::Element,χ,δNow)

    @unpack solver,flapLoadsSolver,flapped,b,aₕ,δdotNow,δddotNow,smallAngles,ϖMid,inertialPitchPlungeStatesRange = element.aero
    @unpack Ma,Tnα,TnM,Tnθ̇ = element.aero.flowParameters
    @unpack α,αₑ = element.aero.flowAnglesAndRates
    @unpack Uᵢ,UₙdotMid,Ωₐ = element.aero.flowVelocitiesAndRates
    @unpack ϵₙ,cnα = element.aero.airfoil.attachedFlowParameters
    @unpack cnδ = element.aero.airfoil.flapParameters
    if typeof(solver) in [Indicial,Inflow]
        @unpack incompressibleInertialLoads = element.aero.solver
    else
        incompressibleInertialLoads = true
    end

    # Circulatory component
    cnC = smallAngles ? ϖMid * cnα * αₑ : ϖMid * cnα * sin(αₑ) * cos(αₑ)
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cnC += ϖMid * cnδ * δNow
    end

    # Inertial component
    cnI = incompressibleInertialLoads ? π*b*UₙdotMid/Uᵢ^2 : 4/Ma * (α - χ[inertialPitchPlungeStatesRange[1]]/Tnα) + -4*aₕ*b/(Uᵢ*Ma) * (Ωₐ - χ[inertialPitchPlungeStatesRange[2]]/Tnθ̇) + 4*α/Ma^2 * (Ma - χ[inertialPitchPlungeStatesRange[6]]/TnM)
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cnI -= ϵₙ*b/Uᵢ^2*(Uᵢ*Th[4]*δdotNow+b*Th[1]*δddotNow)
    end

    # Total
    cn = cnC + cnI

    # Non-circulatory
    cnNC = cnI

    @pack! element.aero.aeroCoefficients = cn,cnC,cnI,cnNC
end


# Computes the pitching moment aerodynamic coefficient at the attachment point (i.e., the beam reference line) for attached flow
function attached_flow_cm!(element::Element,χ,δNow)

    @unpack solver,flapLoadsSolver,flapped,b,aₕ,normFlapPos,δdotNow,δddotNow,ϖMid,inertialPitchPlungeStatesRange = element.aero
    @unpack Ma,βₚ,Tmα,TmM,Tmθ̇ = element.aero.flowParameters
    @unpack α,αₑ = element.aero.flowAnglesAndRates
    @unpack cnC = element.aero.aeroCoefficients
    @unpack Uᵢ,UₙdotMid,Ωₐ,Ωₐdot = element.aero.flowVelocitiesAndRates
    @unpack α₀N,ϵₘ,cm₀,cmα,cnα = element.aero.airfoil.attachedFlowParameters
    @unpack cmδ = element.aero.airfoil.flapParameters
    if typeof(solver) in [Indicial,Inflow]
        @unpack incompressibleInertialLoads = element.aero.solver
        @unpack AI,bI = element.aero.solver
    else
        incompressibleInertialLoads = true
    end

    # Circulatory component
    cmC = cnC * (1/4+aₕ/2) + ϖMid * cmα * αₑ
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cmC += ϖMid * cmδ * δNow
    end

    # Inertial component
    cmI = incompressibleInertialLoads ? -π*b/(2*Uᵢ^2)*(-aₕ*UₙdotMid+b*Ωₐdot/8+Uᵢ*Ωₐ/2) : 2*aₕ/Ma * (α - χ[inertialPitchPlungeStatesRange[3]]*AI[1]/(bI[1]*Tmα) - χ[inertialPitchPlungeStatesRange[4]]*AI[2]/(bI[2]*Tmα)) + -2*b*(aₕ^2+1/3)/(Uᵢ*Ma) * (Ωₐ - χ[inertialPitchPlungeStatesRange[5]]/Tmθ̇) + 2*aₕ*α/Ma^2 * (Ma - χ[inertialPitchPlungeStatesRange[7]]*AI[1]/(bI[1]*TmM) - χ[inertialPitchPlungeStatesRange[8]]*AI[2]/(bI[2]*TmM)) - π*b*Ωₐ/(4*Uᵢ*βₚ)
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cmI -= ϵₘ/(2*Uᵢ^2)*(Uᵢ^2*Th[14]*δNow+Uᵢ*b*Th[15]*δdotNow+b^2*Th[16]*δddotNow)
    end

    # Non-circulatory
    cmNC = cmI + cm₀

    # Total at attachment point
    cm = cmC + cmNC

    @pack! element.aero.aeroCoefficients = cm,cmC,cmI,cmNC

end


# Computes the tangential force aerodynamic coefficient for attached flow
function attached_flow_ct!(element::Element,δNow)

    @unpack flapped,flapLoadsSolver,smallAngles,ϖMid,b = element.aero
    @unpack hasInducedDrag = element.parent.aeroSurface
    @unpack α,αₑ = element.aero.flowAnglesAndRates
    @unpack Uᵢ,Ωₐ = element.aero.flowVelocitiesAndRates
    @unpack cnC = element.aero.aeroCoefficients
    @unpack α₀N,η,cd₀ = element.aero.airfoil.attachedFlowParameters
    @unpack cdδ = element.aero.airfoil.flapParameters

    # Parasite component
    ct₀ = smallAngles ? -cd₀ : -cd₀/cos(α)

    # Circulatory component
    ctC = smallAngles ? η * cnC * (αₑ+b*Ωₐ/(2*Uᵢ)) : η * cnC * tan(αₑ+b*Ωₐ/(2*Uᵢ))
    if flapped && typeof(flapLoadsSolver) == TableLookup
        ctC -= ϖMid * cdδ * abs(δNow) / cos(α)
    end

    # Induced (drag) component
    if hasInducedDrag
        @unpack AR = element.parent.aeroSurface
        ctC -= cnC^2 / (π*AR) / cos(α)
    end

    # Total
    ct = ct₀ + ctC

    # Non-circulatory
    ctNC = ct₀

    @pack! element.aero.aeroCoefficients = ct,ctC,ctNC

end


# Computes the aerodynamic state matrices for the indicial and inflow solvers
function attached_flow_state_matrices!(element::Element,δNow)

    @unpack solver,nTotalAeroStates,pitchPlungeStatesRange,flapStatesRange,gustStatesRange,b = element.aero
    @unpack Ma,βₚ²,Θ,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = element.aero.flowParameters
    @unpack α = element.aero.flowAnglesAndRates
    @unpack Uᵢ,UₙdotTQC,Ωₐ = element.aero.flowVelocitiesAndRates
    if typeof(solver) in [Indicial,Inflow]
        @unpack incompressibleInertialLoads = element.aero.solver
    else
        incompressibleInertialLoads = true
    end

    # Initialize state matrices with appropriate types
    T = promote_type(typeof(α), typeof(UₙdotTQC))
    A = zeros(T, nTotalAeroStates,nTotalAeroStates)
    B = zeros(T, nTotalAeroStates)

    # Skip for quasi-steady solver
    if typeof(solver) == QuasiSteady
        return A,B
    end

    # Pitch-plunge-induced flow states
    if typeof(solver) == Indicial
        @unpack AC,bC,bCMat,bI = element.aero.solver
        # Get the rate of cnαUₙTQC
        cnαUₙTQCdot = cnαUₙTQC_rate(element)
        # Set state matrices
        A[pitchPlungeStatesRange,pitchPlungeStatesRange] .= incompressibleInertialLoads ? -Θ*bCMat : -Diagonal(vcat(Θ*bC, 1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM)))
        B[pitchPlungeStatesRange] .= incompressibleInertialLoads ? cnαUₙTQCdot*AC : vcat(cnαUₙTQCdot*AC, α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)
    elseif typeof(solver) == Inflow
        @unpack nInflowStates,nInertialStates,AₚInv,AₚInvcₚ,bI = element.aero.solver
        # Set state matrices
        A[pitchPlungeStatesRange,pitchPlungeStatesRange] .= incompressibleInertialLoads ? -Θ*AₚInv :
        [-Θ*AₚInv zeros(nInflowStates,nInertialStates); 
        zeros(nInertialStates,nInflowStates) -Diagonal([1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM)])]
        B[pitchPlungeStatesRange] .= incompressibleInertialLoads ? UₙdotTQC*AₚInvcₚ : vcat(UₙdotTQC*AₚInvcₚ, α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)
    end

    # Flap-induced flow states
    if !isnothing(flapStatesRange)
        @unpack ACf,bCfMat = element.aero.flapLoadsSolver
        # Get the rate of cnα*wFlap
        cnαwdotFlap = cnαwFlap_rate(element,δNow)
        # Set state matrices
        A[flapStatesRange,flapStatesRange] .= -Θ*bCfMat
        B[flapStatesRange] .= cnαwdotFlap*ACf
    end

    # Gust-induced flow states
    if !isnothing(gustStatesRange)
        @unpack bGMat = element.aero.gustLoadsSolver
        @unpack UₙGust = element.aero.flowVelocitiesAndRates
        # Set state matrices
        A[gustStatesRange,gustStatesRange] .= -Θ*bGMat
        B[gustStatesRange] .= UₙGust*ones(length(gustStatesRange))
    end

    @pack! element.aero = A,B

    return A,B
end


# Computes the aerodynamic coefficients according to the modified incompressible flow Beddoes-Leishman model
function BLi_aero_coefficients!(problem::Problem,element::Element,χ,δNow)

    # Name nonlinear states
    BLi_nonlinear_states!(element,χ)

    # Additional kinematics
    BLi_kinematics!(element)

    # Motion qualifiers
    BLi_motion_qualifiers!(problem,element)

    # Breakpoint angles
    BLi_breakpoint_angles!(element)

    # Time delay constants
    BLi_time_delays!(element)

    # Flow separation points
    BLi_separation_points!(element)

    # For dynamic problems
    if problem isa DynamicProblem

        # Variables at stall onset
        BLi_stall_time!(problem,element)

        # DSV loads
        BLi_DSV_loads!(element)
    end

    # Update inertial parameters, if applicable
    if !element.aero.solver.incompressibleInertialLoads
        update_inertial_parameters!(element)
    end

    # Normal force coefficient
    BLi_cn!(element,χ,δNow)

    # Pitching moment coefficient about the spar position
    BLi_cm!(element,χ,δNow)

    # Tangential flow coefficient
    BLi_ct!(element,δNow)

end


# Computes additional kinematics associated with the Beddoes-Leishman model
function BLi_kinematics!(element::Element)

    @unpack c = element.aero
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates
    @unpack αdot = element.aero.flowAnglesAndRates
    @unpack r₀ = element.aero.airfoil.parametersBLi

    # Non-dimensional pitch rate (round off to suppress noise)
    q = round_off!(αdot*c/Uᵢ,1e-8)

    # Reduced non-dimensional pitch rate
    r = q/2

    # Unsigned ratio of reduced pitch rate to critical pitch rate
    qR = abs(r)/r₀

    # Unsigned capped reduced pitch rate ratio 
    R = min(1, qR)

    @pack! element.aero.BLiKin = q,r,qR,R

end


# Sets the nonlinear states of the Beddoes-Leishman model
function BLi_nonlinear_states!(element::Element,χ)

    @unpack nonlinearPitchPlungeStatesRange,nonlinearGustStatesRange = element.aero
    @unpack Uₜ,UₜGust = element.aero.flowVelocitiesAndRates
    @unpack f₀N,f₀M,f₀T = element.aero.airfoil.parametersBLi

    if !isnothing(nonlinearGustStatesRange)
        αlag = abs(Uₜ+UₜGust) > 0 ? atan(-(χ[nonlinearPitchPlungeStatesRange[1]]+χ[nonlinearGustStatesRange[1]])/(Uₜ+UₜGust)) : 0.0
    else
        αlag = abs(Uₜ) > 0 ? atan(-χ[nonlinearPitchPlungeStatesRange[1]]/Uₜ) : 0.0
    end
    f2primeN = max(f₀N, χ[nonlinearPitchPlungeStatesRange[2]])
    f2primeM = max(f₀M, χ[nonlinearPitchPlungeStatesRange[3]])
    f2primeT = max(f₀T, χ[nonlinearPitchPlungeStatesRange[4]])
    RD = max(0, χ[nonlinearPitchPlungeStatesRange[5]])
    RD_stallOnsetRatio = max(0, χ[nonlinearPitchPlungeStatesRange[6]])

    @pack! element.aero.BLiStates = αlag,f2primeN,f2primeM,f2primeT,RD,RD_stallOnsetRatio

end


# Computes motion qualifiers of the modified Beddoes-Leishman model
function BLi_motion_qualifiers!(problem::Problem,element::Element)

    @unpack timeNow = problem
    @unpack αds₀,αₛₛ,γLS,Tf = element.aero.airfoil.parametersBLi
    @unpack r,qR = element.aero.BLiKin
    @unpack αlag,RD,RD_stallOnsetRatio = element.aero.BLiStates
    @unpack αlagPrev,stallOnsetRatioPrev,PPrev,qRPrev,upstrokePrev,maxStallOnsetRatio,minStallOnsetRatio,qRmax,tv0P,tv0N,RD_tv0P,RD_tv0N = element.aero.BLiCompVars

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
    P = S ? exp(-γLS*(abs(maxStallOnsetRatio)^4-1)) : PPrev 

    # Steady condition factor
    T = exp(-50*RD)

    # Update maximum stall onset ratio and corresponding delayed capped reduced pitch rate (update only for increments in the lagged AoA)
    if abs(stallOnsetRatio) > abs(stallOnsetRatioPrev) && abs(αlag) > abs(αlagPrev) 
        maxStallOnsetRatio = stallOnsetRatio
        @pack! element.aero.BLiCompVars = maxStallOnsetRatio
    end

    # Update minimum stall onset ratio (approximate as that at begin of upstroke)
    if upstrokeBeginning
        minStallOnsetRatio = stallOnsetRatio
        @pack! element.aero.BLiCompVars = minStallOnsetRatio
    end

    # Update maximum reduced pitch rate ratio
    if qR > qRPrev
        qRmax = qR
        @pack! element.aero.BLiCompVars = qRmax
    end

    # Time decay function for tangential force
    tv0 = max(tv0P, tv0N)
    lastRD_tv0 = tv0P > tv0N ? RD_tv0P : RD_tv0N
    Ts = exp(-5*(timeNow-tv0)*max(0.1,lastRD_tv0)^2/Tf)

    @pack! element.aero.BLiFlow = stallOnsetRatio,upstroke,S,P,T
    @pack! element.aero.BLiCompVars = Ts,lastRD_tv0

end


# Computes the breakpoint angles of the modified Beddoes-Leishman model
function BLi_breakpoint_angles!(element::Element)

    @unpack αds₀,αₛₛ,α1₀N,α1₀M,α1₀T,δα₀,δα₁,dt,dm,zm,ztd = element.aero.airfoil.parametersBLi
    @unpack RD = element.aero.BLiStates
    @unpack qR,R = element.aero.BLiKin
    @unpack upstroke,S,P = element.aero.BLiFlow
    @unpack qRmax = element.aero.BLiCompVars

    # Breakpoint angle offsets
    sqrtOfR = sqrt(R)
    if upstroke 
        fR = R/min(1, qRmax)*sqrtOfR
        δα1N = (αds₀-αₛₛ)*RD        
        δα1M = δα1N+dm*fR                                  
        δα1T = δα1N+dt*fR                                   
    else 
        δα1N = S ? -min(α1₀N/2, δα₀+δα₁*qR*(1+P)) : 0
        δα1M = δα1N*(1-zm*sqrtOfR*(1-P/2))                               
        δα1T = δα1N*(1-ztd*sqrtOfR*(1-P/4))    
    end

    # Unsteady breakpoint of separation angles
    α1N = α1₀N + δα1N 
    α1M = α1₀M + δα1M
    α1T = α1₀T + δα1T

    @pack! element.aero.BLiFlow = α1N,α1M,α1T

end


# Computes time delay variables for the modified Beddoes-Leishman model
function BLi_time_delays!(element::Element)

    @unpack Ta,Tf,λ₁,λ₂ = element.aero.airfoil.parametersBLi
    @unpack RD = element.aero.BLiStates
    @unpack qR,R = element.aero.BLiKin
    @unpack stallOnsetRatio,upstroke,P = element.aero.BLiFlow
    @unpack qRmax = element.aero.BLiCompVars

    # Stall supression time delay
    Ta_SO = Ta*(1/4+λ₁*RD^λ₂*(1-R^(1/4))*exp(-(abs(stallOnsetRatio)-1)^2/0.01))

    # Separation points time delays         
    TfN = Tf*(1+1/2*P^6*(!upstroke)+RD^3*(4*(1-qR/qRmax)*(!upstroke)+1/2*(qR/qRmax)^4*upstroke))
    TfM = Tf*(1-1/2*R*(1-P)*(!upstroke))
    TfT = TfN

    @pack! element.aero.BLiFlow = Ta_SO,TfN,TfM,TfT

end


# Computes the flow separation points of the modified Beddoes-Leishman model
function BLi_separation_points!(element::Element)

    @unpack α1₀N,α1₀M,α1₀T,βσ1N,βσ1T,βσ2N,βS2Nlpr,βS2Tlpr,βS1Nu,βS1Mu,βS1Tu,βS1Nd,βS1Md,βS1Td,βS2Nu,βS2Mu,βS2Tu,βS2Nd,βS2Md,βS2Td,ξ,f₀N,f₀M,f₀T,fbN,fbM,fbT,S1N,S1M,S1T,S2N,S2M,S2T = element.aero.airfoil.parametersBLi
    @unpack α = element.aero.flowAnglesAndRates
    @unpack αlag,RD = element.aero.BLiStates
    @unpack R = element.aero.BLiKin
    @unpack stallOnsetRatio,upstroke,S,P,T,α1N,α1M,α1T = element.aero.BLiFlow
    @unpack minStallOnsetRatio,lastRD_tv0 = element.aero.BLiCompVars

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
        αDiffN = α1N - absαlag
        if upstroke  
            S2PrimeN = S2N*(1+βS2Nu*RD+βS2Nlpr*(1-RD)^2*(1-T)+βσ2N*σ*RD*(1-RD))
            f₀Nprime = f₀N
        else                
            f₀Nprime = min(0.25, (1-sqrtOfR)*((S ? f₀N : fbN-0.2)+ξRD_tv0*σ))
            S2PrimeN = S2N*(1+βS2Nd*RD+βS2Nlpr*(1-RD)^2*(1-T)+βσ2N*σ*RD*(1-RD)+βσ1N*ψ*R)      
        end
        fPrimeN = f₀Nprime+(fbN-f₀Nprime)*exp(αDiffN/S2PrimeN)
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
            f₀Mprime = f₀M
        else                
            f₀Mprime = min(0.25, (1-sqrtOfR)*((S ? f₀M : fbM-0.2)+ξRD_tv0*σ))
            S2PrimeM = S2M*(1+βS2Md*RD)
        end
        fPrimeM = f₀Mprime+(fbM-f₀Mprime)*exp(αDiffM/S2PrimeM)
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
            f₀Tprime = f₀T
        else
            f₀Tprime = min(0.25, (1-sqrtOfR)*((S ? f₀T : fbT-0.2)+ξRD_tv0*σ))
            S2PrimeT = S2T*(1+βS2Td*RD+βS2Tlpr*(1-RD)^2*(1-T)+βσ1T*ψ*R)
        end
        fPrimeT = f₀Tprime+(fbT-f₀Tprime)*exp(αDiffT/S2PrimeT)
    end

    @pack! element.aero.BLiFlow = fN,fM,fT,fPrimeN,fPrimeM,fPrimeT

end


# Computes the quasi-steady flow separation points of the modified Beddoes-Leishman model for an element
function BLi_quasi_steady_separation_points(element::Element,Uₜ,Uₙ)

    # Check tangential velocity
    if isapprox(abs(Uₜ),0)
        fN,fM,fT = 1.0,1.0,1.0
        return fN,fM,fT
    end

    @unpack α1₀N,α1₀M,α1₀T,f₀N,f₀M,f₀T,fbN,fbM,fbT,S1N,S1M,S1T,S2N,S2M,S2T = element.aero.airfoil.parametersBLi

    # Angle of attack  
    α = atan(-Uₙ/Uₜ)
    absα = abs(α)

    # Quasi-steady separation points
    fN = absα <= α1₀N ? 1-(1-fbN)*exp((absα-α1₀N)/S1N) : f₀N+(fbN-f₀N)*exp((α1₀N-absα)/S2N)

    fM = absα <= α1₀M ? 1-(1-fbM)*exp((absα-α1₀M)/S1M) : f₀M+(fbM-f₀M)*exp((α1₀M-absα)/S2M)

    fT = absα <= α1₀T ? 1-(1-fbT)*exp((absα-α1₀T)/S1T) : f₀T+(fbT-f₀T)*exp((α1₀T-absα)/S2T)

    return fN,fM,fT

end


# Computes variables at the time of stall onset for the modified Beddoes-Leishman model
function BLi_stall_time!(problem::Problem,element::Element)

    @unpack timeNow,Δt = problem
    @unpack gᵥ,Tv = element.aero.airfoil.parametersBLi
    @unpack f2primeN,RD = element.aero.BLiStates
    @unpack qR,R = element.aero.BLiKin
    @unpack stallOnsetRatio,upstroke,fN = element.aero.BLiFlow
    @unpack stallOnsetRatioPrev,fDiff_tv0P,qR_tv0P,R_tv0P,RD_tv0P,upstroke_tv0P,fDiff_tv0P2,RD_tv0P2,upstroke_tv0P2,fDiff_tv0N,qR_tv0N,R_tv0N,RD_tv0N,upstroke_tv0N,fDiff_tv0N2,RD_tv0N2,upstroke_tv0N2,tv0P,tv0N = element.aero.BLiCompVars

    # Positive AoA
    #---------------------------------------------------------------------------
    # Primary vortex
    if stallOnsetRatio >= 1 && stallOnsetRatioPrev < 1 
        # Linear interpolation for time of vortex shedding 
        tv0P = LinearInterpolations.interpolate([stallOnsetRatioPrev; stallOnsetRatio], [timeNow-Δt; timeNow], 1)
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
        tv0N = LinearInterpolations.interpolate([-stallOnsetRatioPrev; -stallOnsetRatio], [timeNow-Δt; timeNow], 1)
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

    @pack! element.aero.BLiCompVars = tv0P,fDiff_tv0P,qR_tv0P,R_tv0P,RD_tv0P,upstroke_tv0P,τvP,fDiff_tv0P2,RD_tv0P2,upstroke_tv0P2,tv0N,fDiff_tv0N,qR_tv0N,R_tv0N,RD_tv0N,upstroke_tv0N,τvN,fDiff_tv0N2,RD_tv0N2,upstroke_tv0N2

end


# Computes DSV loads for the modified Beddoes-Leishman model
function BLi_DSV_loads!(element::Element)

    @unpack ϖMid = element.aero
    @unpack μv₂,ν₁,ν₂,ν₃,ν₄,ν₅,χu,χd,gᵥ,gᵥ₂,Tv,Tv₂,Vn₁,Vn₂,Vn₃,Vm,Vt = element.aero.airfoil.parametersBLi
    @unpack qR,R = element.aero.BLiKin
    @unpack stallOnsetRatio,upstroke,P = element.aero.BLiFlow
    @unpack fDiff_tv0P,qR_tv0P,R_tv0P,RD_tv0P,upstroke_tv0P,τvP,fDiff_tv0P2,RD_tv0P2,upstroke_tv0P2,fDiff_tv0N,qR_tv0N,R_tv0N,RD_tv0N,upstroke_tv0N,τvN,fDiff_tv0N2,RD_tv0N2,upstroke_tv0N2 = element.aero.BLiCompVars

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

    # Scale with tip correction factor
    cnV,cmV,ctV = multiply_inplace!(ϖMid, cnV,cmV,ctV)

    @pack! element.aero.aeroCoefficients = cnV,cmV,ctV

end


# Computes normal force coefficient for the modified incompressible Beddoes-Leishman model
function BLi_cn!(element::Element,χ,δNow)

    @unpack incompressibleInertialLoads = element.aero.solver
    @unpack flapLoadsSolver,flapped,b,δdotNow,δddotNow,c,aₕ,ϖMid,inertialPitchPlungeStatesRange = element.aero
    @unpack α,αₑ = element.aero.flowAnglesAndRates
    @unpack Uᵢ,UₙdotMid,Ωₐ = element.aero.flowVelocitiesAndRates
    @unpack α₀N,ϵₙ,cnα = element.aero.airfoil.parametersBLi
    @unpack cnδ = element.aero.airfoil.flapParameters
    @unpack f2primeN = element.aero.BLiStates
    @unpack cnV = element.aero.aeroCoefficients
    @unpack Ma = element.aero.flowParameters
    @unpack Tnα,TnM,Tnθ̇ = element.aero.flowParameters

    # Circulatory component - attached flow
    cnC = ϖMid * cnα * sin(αₑ)
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cnC += ϖMid * cnδ * δNow
    end

    # Circulatory component - separated flow
    cnF = cnC * ((1+sqrt(f2primeN))/2)^2

    # Inertial component
    cnI = incompressibleInertialLoads ? π*b*UₙdotMid/Uᵢ^2 : 4/Ma * (α - χ[inertialPitchPlungeStatesRange[1]]/Tnα) + -4*aₕ*b/(Uᵢ*Ma) * (Ωₐ - χ[inertialPitchPlungeStatesRange[2]]/Tnθ̇) + 4*α/Ma^2 * (Ma - χ[inertialPitchPlungeStatesRange[6]]/TnM)
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cnI -= ϵₙ*b/Uᵢ^2*(Uᵢ*Th[4]*δdotNow+b*Th[1]*δddotNow)
    end

    # Total 
    cn = cnF + cnI + cnV

    # Non-circulatory
    cnNC = cnI

    @pack! element.aero.aeroCoefficients = cn,cnC,cnI,cnF,cnV,cnNC

end


# Computes the pitching moment aerodynamic coefficient at the attachment point (i.e., the beam reference line) for the modified incompressible Beddoes-Leishman model
function BLi_cm!(element::Element,χ,δNow)

    @unpack incompressibleInertialLoads,AI,bI = element.aero.solver
    @unpack flapLoadsSolver,flapped,b,normSparPos,normFlapPos,δdotNow,δddotNow,aₕ,ϖMid,inertialPitchPlungeStatesRange = element.aero
    @unpack α = element.aero.flowAnglesAndRates
    @unpack cnF = element.aero.aeroCoefficients
    @unpack Uᵢ,UₙdotMid,Ωₐ,Ωₐdot = element.aero.flowVelocitiesAndRates
    @unpack ϵₘ,κ₀,κ₁,κ₂,κ₃,cm₀,cnα,K₀,K₁,K₂ = element.aero.airfoil.parametersBLi
    @unpack cmδ = element.aero.airfoil.flapParameters
    @unpack R = element.aero.BLiKin
    @unpack stallOnsetRatio,upstroke,S,P = element.aero.BLiFlow
    @unpack Tmα,TmM,Tmθ̇ = element.aero.flowParameters
    @unpack f2primeM,RD = element.aero.BLiStates
    @unpack cmV = element.aero.aeroCoefficients
    @unpack Ma,βₚ = element.aero.flowParameters

    # Center of pressure dynamic variables
    K1Prime = K₁*(1-κ₁*RD*(1-abs(stallOnsetRatio)))-κ₂*R*abs(stallOnsetRatio)*upstroke
    K2Prime = K₂*(1+κ₃*S*R^2*(1-P)*(!upstroke))                        

    # Center of pressure offset from quarter-chord
    δCP = K₀ + K1Prime*(1-f2primeM) + K2Prime*sin(π*f2primeM^κ₀)

    # Circulatory component - separated flow
    cmF = cnF * (normSparPos-(1/4-δCP))
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cmF += ϖMid * cmδ * δNow
    end

    # Inertial component
    cmI = incompressibleInertialLoads ? -π*b/(2*Uᵢ^2)*(-aₕ*UₙdotMid+b*Ωₐdot/8+Uᵢ*Ωₐ/2) : 2*aₕ/Ma * (α - χ[inertialPitchPlungeStatesRange[3]]*AI[1]/(bI[1]*Tmα) - χ[inertialPitchPlungeStatesRange[4]]*AI[2]/(bI[2]*Tmα)) + -2*b*(aₕ^2+1/3)/(Uᵢ*Ma) * (Ωₐ - χ[inertialPitchPlungeStatesRange[5]]/Tmθ̇) + 2*aₕ*α/Ma^2 * (Ma - χ[inertialPitchPlungeStatesRange[7]]*AI[1]/(bI[1]*TmM) - χ[inertialPitchPlungeStatesRange[8]]*AI[2]/(bI[2]*TmM)) - π*b*Ωₐ/(4*Uᵢ*βₚ)
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cmI -= ϵₘ/(2*Uᵢ^2)*(Uᵢ^2*Th[14]*δNow+Uᵢ*b*Th[15]*δdotNow+b^2*Th[16]*δddotNow)
    end

    # Non-circulatory
    cmNC = cmI + cm₀

    # Total at attachment point
    cm = cmF + cmV + cmNC

    @pack! element.aero.aeroCoefficients = cm,cmF,cmI,cmNC

end


# Computes the tangential force aerodynamic coefficient for the modified incompressible Beddoes-Leishman model
function BLi_ct!(element::Element,δNow)

    @unpack flapped,flapLoadsSolver,ϖMid,b = element.aero
    @unpack hasInducedDrag = element.parent.aeroSurface
    @unpack α,αₑ = element.aero.flowAnglesAndRates
    @unpack Uᵢ,Ωₐ = element.aero.flowVelocitiesAndRates
    @unpack α₀N,α1₀T,η,cd₀,cnα,E₀,E₁ = element.aero.airfoil.parametersBLi
    @unpack cdδ = element.aero.airfoil.flapParameters
    @unpack R = element.aero.BLiKin
    @unpack stallOnsetRatio,upstroke,S = element.aero.BLiFlow
    @unpack αlag,f2primeT,RD = element.aero.BLiStates
    @unpack cnF,ctV = element.aero.aeroCoefficients
    @unpack Ts = element.aero.BLiCompVars

    # Ratio of lagged AoA to tangential steady breakpoint angle
    stallOnsetRatioT = abs(αlag)/α1₀T

    # Parasite component
    ct₀ = -cd₀/cos(α)

    # Circulatory component - separated flow
    ctF = (1-η*RD^2) * (cnF/ϖMid) * tan(αₑ+b*Ωₐ/(2*Uᵢ)) * f2primeT^(1/2 + stallOnsetRatioT + E₀*RD*S*(1-Ts)*R*abs(stallOnsetRatio)^(1/2)*(!upstroke))
    if stallOnsetRatioT > 1
        ctF -= E₁ * (1-RD^3) * min(1, stallOnsetRatioT^3-1)
    end
    if flapped && typeof(flapLoadsSolver) == TableLookup
        ctF -= cdδ * abs(δNow) / cos(α)
    end
    ctF *= ϖMid

    # Induced (drag) component
    if hasInducedDrag
        @unpack AR = element.parent.aeroSurface
        ctF -= cnF^2 / (π*AR) / cos(α)
    end

    # Total
    ct = ct₀ + ctF + ctV

    # Non-circulatory
    ctNC = ct₀

    @pack! element.aero.aeroCoefficients = ct,ctF,ctV,ctNC

end


# Computes the aerodynamic state matrices for the modified incompressible Beddoes-Leishman model
function BLi_state_matrices!(element::Element,δNow)

    @unpack nTotalAeroStates,linearPitchPlungeStatesRange,nonlinearPitchPlungeStatesRange,flapStatesRange,linearGustStatesRange,nonlinearGustStatesRange,flapped,flapLoadsSolver = element.aero
    @unpack βₚ²,Θ,Ma,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = element.aero.flowParameters
    @unpack Ta,γbC,γbCMat = element.aero.airfoil.parametersBLi
    @unpack α = element.aero.flowAnglesAndRates
    @unpack Uₙ,Uᵢ,UₙdotTQC,Ωₐ = element.aero.flowVelocitiesAndRates
    @unpack R = element.aero.BLiKin
    @unpack Ta_SO,TfN,TfM,TfT,fPrimeN,fPrimeM,fPrimeT = element.aero.BLiFlow
    @unpack incompressibleInertialLoads,AC,bC,bI,bCMat = element.aero.solver

    # Initialize state matrices with appropriate types
    T = promote_type(typeof(α), typeof(UₙdotTQC))
    A = zeros(T, nTotalAeroStates,nTotalAeroStates)
    B = zeros(T, nTotalAeroStates)
        
    # Rate of cnαUₙTQC
    cnαUₙTQCdot = cnαUₙTQC_rate(element)

    # Quasi-steady flap-induced normalwash (considered only from thin-airfoil theory)
    wFlap = flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory ? flap_normalwash(element,δNow) : 0

    # Linear pitch-plunge-induced flow state matrices
    A[linearPitchPlungeStatesRange,linearPitchPlungeStatesRange] .= incompressibleInertialLoads ? -Θ * bCMat .* γbCMat : -Diagonal(vcat(Θ*bC.*γbC, 1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM)))
    B[linearPitchPlungeStatesRange] .= incompressibleInertialLoads ? cnαUₙTQCdot * AC : vcat(cnαUₙTQCdot*AC, α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)
    
    # Nonlinear pitch-plunge-induced flow state matrices
    A[nonlinearPitchPlungeStatesRange,nonlinearPitchPlungeStatesRange] .= -Diagonal([1/Ta, 1/TfN, 1/TfM, 1/TfT, 1/(3*Ta), 1/Ta_SO])
    B[nonlinearPitchPlungeStatesRange] .= [(Uₙ+wFlap)/Ta, fPrimeN/TfN, fPrimeM/TfM, fPrimeT/TfT, R/(3*Ta), R/Ta_SO]

    # Flap-induced flow states
    if !isnothing(flapStatesRange)
        @unpack ACf,bCfMat = element.aero.flapLoadsSolver
        # Get the rate of cnα*wFlap
        cnαwdotFlap = cnαwFlap_rate(element,δNow)
        # Set state matrices
        A[flapStatesRange,flapStatesRange] .= -Θ*bCfMat
        B[flapStatesRange] .= cnαwdotFlap*ACf
    end

    # Gust-induced flow state matrices
    if !isnothing(linearGustStatesRange)
        @unpack bGMat = element.aero.gustLoadsSolver
        @unpack UₙGust = element.aero.flowVelocitiesAndRates
        @unpack Tg = element.aero.airfoil.parametersBLi
        # Linear part 
        A[linearGustStatesRange,linearGustStatesRange] .= -Θ*bGMat
        B[linearGustStatesRange] .= UₙGust*ones(length(linearGustStatesRange))
        # Nonlinear part
        A[nonlinearGustStatesRange,nonlinearGustStatesRange] .= -1/Tg
        B[nonlinearGustStatesRange] .= UₙGust/Tg
    end

    @pack! element.aero = A,B

    return A,B
end


# Computes the aerodynamic coefficients according to the original Beddoes-Leishman model
function BLo_aero_coefficients!(problem::Problem,element::Element,χ,δNow)

    # Name nonlinear states
    BLo_nonlinear_states!(element,χ)

    # Motion qualifiers
    BLo_motion_qualifiers!(element)

    # Breakpoint angle
    BLo_breakpoint_angle!(element)

    # Flow separation points
    BLo_separation_points!(element)

    # For dynamic problems
    if problem isa DynamicProblem
        # Variables at stall onset
        BLo_stall_time!(problem,element)
    end

    # Time delay constants
    BLo_time_delays!(element)

    # Update inertial parameters, if applicable
    if !element.aero.solver.incompressibleInertialLoads
        update_inertial_parameters!(element)
    end

    # Normal force coefficient
    BLo_cn!(element,χ,δNow)

    # DSV accumulation rate
    BLo_vortex_accumulation_rate!(element,χ)

    # Pitching moment coefficient about the spar position
    BLo_cm!(element,χ,δNow)

    # Tangential flow coefficient
    BLo_ct!(element,δNow)

end


# Sets the nonlinear states of the original Beddoes-Leishman model
function BLo_nonlinear_states!(element::Element,χ)

    @unpack nonlinearPitchPlungeStatesRange,ϖMid = element.aero
    @unpack f₀ = element.aero.airfoil.parametersBLo

    cnPprime = χ[nonlinearPitchPlungeStatesRange[1]]
    f2Prime = min(max(χ[nonlinearPitchPlungeStatesRange[2]],f₀),1)
    fPrimeM = min(max(χ[nonlinearPitchPlungeStatesRange[3]],f₀),1)
    cnVP = χ[nonlinearPitchPlungeStatesRange[4]]
    cnVN = χ[nonlinearPitchPlungeStatesRange[5]]
    cnV = ϖMid * (cnVP + cnVN)

    @pack! element.aero.BLoStates = cnPprime,f2Prime,fPrimeM,cnVP,cnVN,cnV

end


# Computes motion qualifiers of the original Beddoes-Leishman model
function BLo_motion_qualifiers!(element::Element)

    @unpack c = element.aero
    @unpack Uᵢ = element.aero.flowVelocitiesAndRates
    @unpack αdot = element.aero.flowAnglesAndRates
    @unpack cn₁,cnα = element.aero.airfoil.parametersBLo
    @unpack cnPprime = element.aero.BLoStates

    # Non-dimensional pitch rate (round off to supress noise)
    q = round_off!(αdot*c/Uᵢ,1e-8)

    # Lagged angle of attack
    αlag = cnPprime/cnα

    # Dynamic stall onset ratio
    stallOnsetRatio = cnPprime/cn₁

    # TF for upstroke phase
    upstroke = stallOnsetRatio * q >= 0

    @pack! element.aero.BLoFlow = αlag,q,stallOnsetRatio,upstroke

end


# Computes the breakpoint angle of the original Beddoes-Leishman model
function BLo_breakpoint_angle!(element::Element)

    @unpack α1₀,δα = element.aero.airfoil.parametersBLo
    @unpack upstroke = element.aero.BLoFlow
    @unpack f2Prime = element.aero.BLoStates

    # Breakpoint angle offset
    δα1 = upstroke ? 0 : -(1-f2Prime)^(1/4)*δα

    # Unsteady breakpoint of separation angle
    α1 = α1₀ + δα1

    @pack! element.aero.BLoFlow = α1

end


# Computes the flow separation point of the original Beddoes-Leishman model
function BLo_separation_points!(element::Element)

    @unpack α1₀,cnα,f₀,fb,S1,S2 = element.aero.airfoil.parametersBLo
    @unpack α = element.aero.flowAnglesAndRates
    @unpack α1,αlag = element.aero.BLoFlow
    @unpack cnPprime = element.aero.BLoStates

    # Absolute values of angles
    absα = abs(α)
    absαlag = abs(αlag)

    # Quasi-steady separation point
    f = absα <= α1₀ ? 1-(1-fb)*exp((absα-α1₀)/S1) : f₀+(fb-f₀)*exp((α1₀-absα)/S2)

    # Unsteady lagged separation point based on lagged AoA
    fPrime = absαlag <= α1 ? 1-(1-fb)*exp((absαlag-α1)/S1) : f₀+(fb-f₀)*exp((α1-absαlag)/S2)

    @pack! element.aero.BLoFlow = f,fPrime

end


# Computes variables at the time of stall onset for the original Beddoes-Leishman model
function BLo_stall_time!(problem::Problem,element::Element)

    @unpack timeNow,Δt = problem
    @unpack stallOnsetRatio = element.aero.BLoFlow
    @unpack stallOnsetRatioPrev,tv0P,tv0N = element.aero.BLoCompVars

    # Positive AoA
    #---------------------------------------------------------------------------
    # Linear interpolation for time of vortex shedding
    if stallOnsetRatio >= 1 && stallOnsetRatioPrev < 1  
        tv0P = LinearInterpolations.interpolate([stallOnsetRatioPrev; stallOnsetRatio], [timeNow-Δt; timeNow], 1)
    end
    # Time since primary vortex shedding 
    τvP = timeNow - tv0P

    # Negative AoA
    #---------------------------------------------------------------------------
    # Linear interpolation for time of vortex shedding
    if stallOnsetRatio <= -1 && stallOnsetRatioPrev > -1 
        # Linear interpolation for time of vortex shedding 
        tv0N = LinearInterpolations.interpolate([-stallOnsetRatioPrev; -stallOnsetRatio], [timeNow-Δt; timeNow], 1)
    end
    # Time since primary vortex shedding 
    τvN = timeNow - tv0N

    @pack! element.aero.BLoCompVars = tv0P,τvP,tv0N,τvN

end


# Computes time delay variables for the original Beddoes-Leishman model
function BLo_time_delays!(element::Element)

    @unpack Tf₀,Tv₀,TvL,fb = element.aero.airfoil.parametersBLo
    @unpack stallOnsetRatio,upstroke = element.aero.BLoFlow
    @unpack f2Prime = element.aero.BLoStates
    @unpack τvP,τvN = element.aero.BLoCompVars

    τv = min(τvP,τvN)
    
    # Vortex shedding phase
    if abs(stallOnsetRatio) >= 1
        # Primary vortex
        if τv<=TvL && upstroke
            # Accelerate the rate of boundary layer dettachment during the vortex convection
            Tf = Tf₀*3/4
            # Nominal rate of vortex accumulation
            Tv = Tv₀
        elseif TvL<τv<=2*TvL && upstroke
            # Accelerate more the rate of dettachment after the vortex reaches the trailing edge 
            Tf = Tf₀/2
            # Increase the rate of decay of the vortex after it reaches the trailing edge 
            Tv = Tv₀/2
        elseif τv>2*TvL && upstroke
            # Maintain the rate of boundary layer dettachment after the vortex is totally shed
            Tf = Tf₀/2
            # Maintain a high rate of vortex lift decay after the vortex is totally shed
            Tv = Tv₀/2
        elseif τv<=2*TvL && !upstroke
            # Maintain the rate of dettachment if the rate of change of AoA changes during the vortex shedding 
            Tf = Tf₀/2
            # Increase the rate of decay of the vortex lift if the rate of change of AoA changes during the vortex shedding 
            Tv = Tv₀/2
        elseif τv>2*TvL && !upstroke
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

    @pack! element.aero.BLoFlow = Tf,Tv

end


# Compute the dynamic stall vortex accumulation rate for the original Beddoes-Leishman model
function BLo_vortex_accumulation_rate!(element::Element,χ)

    @unpack c,circulatoryPitchPlungeStatesRange = element.aero
    @unpack Θ = element.aero.flowParameters
    @unpack AC,bC = element.aero.solver
    @unpack Uₜ,UₙdotTQC = element.aero.flowVelocitiesAndRates
    @unpack cnα,TvL,γbC = element.aero.airfoil.parametersBLo
    @unpack stallOnsetRatio,Kf,fPrime,Tf = element.aero.BLoFlow
    @unpack f2Prime = element.aero.BLoStates
    @unpack τvP,τvN = element.aero.BLoCompVars
    @unpack cnC = element.aero.aeroCoefficients

    # Precomupte quantities, if vortex is over the chord
    if abs(stallOnsetRatio) >= 1 && (τvP <= TvL || τvN <= TvL)
        # Kirchhoff-Helmholtz factor time rate
        Kfdot = 1/4*(1+1/sqrt(f2Prime))*(fPrime-f2Prime)/Tf
        # Circulatory unsteady normal force coefficient time rate (assume cnC≈cnα*αₑ, αₑ≈wₑ/Uₜ, and neglect time rate of U and cnα)
        cnCdot = (-Θ*sum(bC.*γbC.*χ[circulatoryPitchPlungeStatesRange])+cnα*UₙdotTQC*(1-sum(AC)))/Uₜ
        # Vorticity coefficient rate
        cvdot = cnCdot*(1-Kf) - cnC*Kfdot
    end

    # Positive vorticity coefficient rate (if positive vortex is over the chord)
    cvdotP = (stallOnsetRatio >= 1 && τvP <= TvL) ? max(0,cvdot) : 0

    # Negative vorticity coefficient rate (if positive vortex is over the chord)
    cvdotN = (stallOnsetRatio <= -1 && τvN <= TvL) ? min(0,cvdot) : 0

    @pack! element.aero.BLoFlow = cvdotP,cvdotN

end


# Computes normal force coefficient for the original Beddoes-Leishman model
function BLo_cn!(element::Element,χ,δNow)

    @unpack incompressibleInertialLoads = element.aero.solver
    @unpack flapLoadsSolver,flapped,b,δdotNow,δddotNow,aₕ,ϖMid,inertialPitchPlungeStatesRange = element.aero
    @unpack α,αₑ = element.aero.flowAnglesAndRates
    @unpack Uᵢ,UₙdotMid,Ωₐ = element.aero.flowVelocitiesAndRates
    @unpack Ma,Tnα,TnM,Tnθ̇ = element.aero.flowParameters
    @unpack α₀N,ϵₙ,cnα = element.aero.airfoil.parametersBLo
    @unpack cnδ = element.aero.airfoil.flapParameters
    @unpack f2Prime,cnV = element.aero.BLoStates
    @unpack fPrime,Tf = element.aero.BLoFlow

    # Circulatory component - attached flow
    cnC = ϖMid * cnα * sin(αₑ)
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cnC += ϖMid * cnδ * δNow
    end

    # Kirchhoff-Helmholtz factor
    Kf = ((1+sqrt(f2Prime))/2)^2

    # Circulatory component - separated flow
    cnF = cnC * Kf

    # Inertial component
    cnI = incompressibleInertialLoads ? π*b*UₙdotMid/Uᵢ^2 : 4/Ma * (α - χ[inertialPitchPlungeStatesRange[1]]/Tnα) + -4*aₕ*b/(Uᵢ*Ma) * (Ωₐ - χ[inertialPitchPlungeStatesRange[2]]/Tnθ̇) + 4*α/Ma^2 * (Ma - χ[inertialPitchPlungeStatesRange[6]]/TnM)
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cnI -= ϵₙ*b/Uᵢ^2*(Uᵢ*Th[4]*δdotNow+b*Th[1]*δddotNow)
    end

    # Potential flow component
    cnP = cnC + cnI

    # Total 
    cn = cnF + cnI + cnV

    # Non-circulatory
    cnNC = cnI

    @pack! element.aero.aeroCoefficients = cn,cnC,cnI,cnF,cnP,cnV,cnNC
    @pack! element.aero.BLoFlow = Kf

end


# Computes the pitching moment aerodynamic coefficient at the attachment point (i.e., the beam reference line) for the original Beddoes-Leishman model
function BLo_cm!(element::Element,χ,δNow)

    @unpack flapLoadsSolver,flapped,normSparPos,normFlapPos,δdotNow,δddotNow,b,aₕ,ϖMid,inertialPitchPlungeStatesRange = element.aero
    @unpack incompressibleInertialLoads,AI,bI = element.aero.solver
    @unpack α = element.aero.flowAnglesAndRates
    @unpack Uᵢ,UₙdotMid,Ωₐ,Ωₐdot = element.aero.flowVelocitiesAndRates
    @unpack cnC,cnF = element.aero.aeroCoefficients
    @unpack ϵₘ,cm₀,K₀,K₁,K₂,TvL = element.aero.airfoil.parametersBLo
    @unpack cmδ = element.aero.airfoil.flapParameters
    @unpack Ma,βₚ,Tmα,TmM,Tmθ̇ = element.aero.flowParameters
    @unpack f2Prime,fPrimeM,cnVP,cnVN = element.aero.BLoStates
    @unpack τvP,τvN = element.aero.BLoCompVars

    # Circulatory unsteady - attached flow
    cmC = cnC * (normSparPos-(1/4-K₀))
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cmC += ϖMid * cmδ * δNow
    end

    # Lagged separation point and center of pressure offset from quarter-chord
    f2PrimeM = max(f2Prime, fPrimeM)
    δCP = K₀ + K₁*(1-f2PrimeM) + K₂*sin(π*f2PrimeM^2)

    # Circulatory unsteady - separated flow
    cmF = cnF * (normSparPos-(1/4-δCP))
    if flapped && typeof(flapLoadsSolver) == TableLookup
        cmF += ϖMid * cmδ * δNow
    end

    # Inertial
    cmI = incompressibleInertialLoads ? -π*b/(2*Uᵢ^2)*(-aₕ*UₙdotMid+b*Ωₐdot/8+Uᵢ*Ωₐ/2) : 2*aₕ/Ma * (α - χ[inertialPitchPlungeStatesRange[3]]*AI[1]/(bI[1]*Tmα) - χ[inertialPitchPlungeStatesRange[4]]*AI[2]/(bI[2]*Tmα)) + -2*b*(aₕ^2+1/3)/(Uᵢ*Ma) * (Ωₐ - χ[inertialPitchPlungeStatesRange[5]]/Tmθ̇) + 2*aₕ*α/Ma^2 * (Ma - χ[inertialPitchPlungeStatesRange[7]]*AI[1]/(bI[1]*TmM) - χ[inertialPitchPlungeStatesRange[8]]*AI[2]/(bI[2]*TmM)) - π*b*Ωₐ/(4*Uᵢ*βₚ)
    if flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory
        @unpack Th = flapLoadsSolver
        cmI -= ϵₘ/(2*Uᵢ^2)*(Uᵢ^2*Th[14]*δNow+Uᵢ*b*Th[15]*δdotNow+b^2*Th[16]*δddotNow)
    end

    # Vortex
    CPvP = τvP <= 2*TvL ? 0.25*(1-cos(π*τvP/TvL)) : 0
    CPvN = τvN <= 2*TvL ? 0.25*(1-cos(π*τvN/TvL)) : 0
    cmV = - ϖMid * (CPvP * cnVP + CPvN * cnVN)
    
    # Non-circulatory
    cmNC = cmI + cm₀

    # Total at attachment point
    cm = cmF + cmV + cmNC

    @pack! element.aero.aeroCoefficients = cm,cmC,cmF,cmI,cmV,cmNC

end


# Computes the tangential force aerodynamic coefficient for the original Beddoes-Leishman model
function BLo_ct!(element::Element,δNow)

    @unpack flapped,flapLoadsSolver,ϖMid,b = element.aero
    @unpack hasInducedDrag = element.parent.aeroSurface
    @unpack α,αₑ = element.aero.flowAnglesAndRates
    @unpack Uᵢ,Ωₐ = element.aero.flowVelocitiesAndRates
    @unpack α₀N,η,cd₀,cnα,cn₁,Df,E₀ = element.aero.airfoil.parametersBLo
    @unpack cdδ = element.aero.airfoil.flapParameters
    @unpack stallOnsetRatio = element.aero.BLoFlow
    @unpack cnPprime,f2Prime = element.aero.BLoStates
    @unpack cnF = element.aero.aeroCoefficients

    # Stall factor
    ϕ = abs(stallOnsetRatio) < 1 ? 0 : Df*(abs(cnPprime)-cn₁)

    # Parasite component
    ct₀ = -cd₀/cos(α)

    # Circulatory component - separated flow
    ctF = η * cnF * tan(αₑ+b*Ωₐ/(2*Uᵢ)) * (sqrt(f2Prime)*f2Prime^ϕ-E₀)
    if flapped && typeof(flapLoadsSolver) == TableLookup
        ctF -= ϖMid * cdδ * abs(δNow) / cos(α)
    end

    # Induced (drag) component
    if hasInducedDrag
        @unpack AR = element.parent.aeroSurface
        ctF -= cnF^2 / (π*AR) / cos(α)
    end
    
    # Total
    ct = ct₀ + ctF

    # Non-circulatory
    ctNC = ct₀

    @pack! element.aero.aeroCoefficients = ct,ctF,ctNC

end


# Computes the aerodynamic state matrices for the original Beddoes-Leishman model
function BLo_state_matrices!(element::Element,δNow)

    @unpack nTotalAeroStates,linearPitchPlungeStatesRange,nonlinearPitchPlungeStatesRange,flapStatesRange,linearGustStatesRange,nonlinearGustStatesRange = element.aero
    @unpack incompressibleInertialLoads,AC,bC,bI,bCMat = element.aero.solver
    @unpack βₚ²,Θ,Ma,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = element.aero.flowParameters
    @unpack Tf₀,Tp,γbC,γbCMat = element.aero.airfoil.parametersBLo
    @unpack α = element.aero.flowAnglesAndRates
    @unpack Uᵢ,Uᵢdot,UₙTQC,Ωₐ = element.aero.flowVelocitiesAndRates
    @unpack Tf,Tv,f,fPrime,cvdotP,cvdotN = element.aero.BLoFlow
    @unpack cnP = element.aero.aeroCoefficients

    # Initialize state matrices with appropriate types
    T = promote_type(typeof(α), typeof(Uᵢdot))
    A = zeros(T, nTotalAeroStates,nTotalAeroStates)
    B = zeros(T, nTotalAeroStates)

    # Get the rate of cnαUₙTQC
    cnαUₙTQCdot = cnαUₙTQC_rate(element)

    # Linear pitch-plunge-induced flow state matrices
    A[linearPitchPlungeStatesRange,linearPitchPlungeStatesRange] .= incompressibleInertialLoads ? -Θ * bCMat .* γbCMat : -Diagonal(vcat(Θ*bC.*γbC, 1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM)))
    B[linearPitchPlungeStatesRange] .= incompressibleInertialLoads ? cnαUₙTQCdot * AC : vcat(cnαUₙTQCdot*AC, α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)
    
    # Nonlinear pitch-plunge-induced flow state matrices
    A[nonlinearPitchPlungeStatesRange, nonlinearPitchPlungeStatesRange] .= -Diagonal([1/Tp, 1/Tf, 2/Tf₀, 1/Tv, 1/Tv])
    B[nonlinearPitchPlungeStatesRange] .= [cnP/Tp, fPrime/Tf, f/(Tf₀/2), cvdotP, cvdotN]

    # Flap-induced flow states
    if !isnothing(flapStatesRange)
        @unpack ACf,bCfMat = element.aero.flapLoadsSolver
        # Get the rate of cnα*wFlap
        cnαwdotFlap = cnαwFlap_rate(element,δNow)
        # Set state matrices
        A[flapStatesRange,flapStatesRange] .= -Θ*bCfMat
        B[flapStatesRange] .= cnαwdotFlap*ACf
    end

    # Gust-induced flow state matrices
    if !isnothing(linearGustStatesRange)
        @unpack bGMat = element.aero.gustLoadsSolver
        @unpack UₙGust = element.aero.flowVelocitiesAndRates
        A[linearGustStatesRange,linearGustStatesRange] .= -Θ*bGMat
        B[linearGustStatesRange] .= UₙGust*ones(length(linearGustStatesRange))
    end

    @pack! element.aero = A,B

    return A,B
end


# Updates the inertial indicial time scales
function update_inertial_parameters!(element::Element)

    @unpack solver = element.aero
    @unpack AC,bC,AI,bI = solver
    @unpack Ma,βₚ,Tᵢ = element.aero.flowParameters

    # Set/get airfoil parameters according to solver
    γbC = ones(length(bC))
    if typeof(solver) == BLi
        @unpack γbC = element.aero.airfoil.parametersBLi
    elseif typeof(solver) == BLo
        @unpack γbC = element.aero.airfoil.parametersBLo    
    end

    # Useful terms
    ΣγACbC = sum(AC.*bC.*γbC)

    # Inertial indicial time scale parameters
    Tnα = Tᵢ/(1-Ma+π*βₚ*Ma^2*ΣγACbC)
    TnM = Tᵢ/(1-Ma+π/βₚ*Ma^2*ΣγACbC)
    Tnθ̇ = Tᵢ/(1-Ma+2π*βₚ*Ma^2*ΣγACbC)
    Tmα = Tᵢ*(AI[1]*bI[2]+AI[2]*bI[1])/(bI[1]*bI[2]*(1-Ma))
    TmM = Tmα
    Tmθ̇ = 7*Tᵢ/(15*(1-Ma)+3π*βₚ*Ma^2*AI[3]*bI[3])

    # Apply correction factor (see Jose et al. - Unsteady Aerodynamic Modeling with Time-Varying Free-Stream Mach Numbers - 2006)
    λ = 0.75
    Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = multiply_inplace!(λ,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇)
    
    @pack! element.aero.flowParameters = Tᵢ,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇

end


# Updates the initial aerodynamic states assuming their rates are zero
function update_initial_aero_states!(problem::Problem;preInitialization::Bool=false)

    @unpack x = problem

    # Loop over elements
    for element in problem.model.elements

        @unpack DOF_χ = element

        # Skip elements without aerodynamic states
        if isnothing(element.aero) || typeof(element.aero.solver) == QuasiSteady
            continue 
        end

        @unpack solver = element.aero

        # Set aerodynamic quantities if in pre-initialization
        if preInitialization
            aero_steady_kinematics!(element,element.states.V,element.states.Ω)
            aero_unsteady_kinematics!(element,element.statesRates.Vdot,element.statesRates.Ωdot)
            nondimensional_flow_parameters!(problem.model,element)
        end

        @unpack incompressibleInertialLoads = element.aero.solver
        @unpack pitchPlungeStatesRange,flapStatesRange,b = element.aero
        @unpack χ = element.states
        @unpack Θ,Ma,Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = element.aero.flowParameters

        # Update inertial parameters, if applicable
        if !incompressibleInertialLoads
            update_inertial_parameters!(element)
        end

        # Pitch-plunge states
        if typeof(solver) == Indicial
            @unpack AC,bC,bCMat,bI = solver
            @unpack α = element.aero.flowAnglesAndRates
            @unpack Ωₐ = element.aero.flowVelocitiesAndRates
            @unpack Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = element.aero.flowParameters
            # Time derivative of the product cnα * UₙTQC
            cnαUₙTQCdot = cnαUₙTQC_rate(element)
            # State matrices
            A = incompressibleInertialLoads ? -Θ*bCMat : -Diagonal(vcat(Θ*bC, 1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM)))
            B = incompressibleInertialLoads ? cnαUₙTQCdot*AC : vcat(cnαUₙTQCdot*AC, α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)
            # States
            χ[pitchPlungeStatesRange] = -A\B
        elseif typeof(solver) == Inflow
            @unpack nInflowStates,nInertialStates,AₚInv,AₚInvcₚ,bI = element.aero.solver
            @unpack UₙdotTQC = element.aero.flowVelocitiesAndRates
            @unpack α = element.aero.flowAnglesAndRates
            @unpack Ωₐ = element.aero.flowVelocitiesAndRates
            @unpack Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = element.aero.flowParameters
            # State matrices
            A = incompressibleInertialLoads ? -Θ*AₚInv :
            [-Θ*AₚInv zeros(nInflowStates,nInertialStates); 
            zeros(nInertialStates,nInflowStates) -Diagonal([1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM)])]
            B = incompressibleInertialLoads ? UₙdotTQC*AₚInvcₚ : vcat(UₙdotTQC*AₚInvcₚ, α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma)
            # States
            χ[pitchPlungeStatesRange] = -A\B
        elseif typeof(solver) == BLi
            # Update kinematics
            BLi_kinematics!(element)
            # Unpack data
            @unpack bC,AC,bI = solver
            @unpack Ta,cnα,γbC = element.aero.airfoil.parametersBLi
            @unpack α = element.aero.flowAnglesAndRates
            @unpack R = element.aero.BLiKin
            @unpack Ta_SO,TfN,TfM,TfT = element.aero.BLiFlow
            @unpack Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = element.aero.flowParameters
            @unpack Uₜ,Uₙ,Ωₐ = element.aero.flowVelocitiesAndRates
            # Steady separation points
            fN,fM,fT = BLi_quasi_steady_separation_points(element,Uₜ,Uₙ)
            # Time derivative of the product cnα * UₙTQC
            cnαUₙTQCdot = cnαUₙTQC_rate(element)
            # State matrices
            A = incompressibleInertialLoads ? -Diagonal(vcat(Θ*bC.*γbC, 1/Ta, 1/TfN, 1/TfM, 1/TfT, 1/(3*Ta), 1/Ta_SO)) : -Diagonal(vcat(Θ*bC.*γbC, 1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM), 1/Ta, 1/TfN, 1/TfM, 1/TfT, 1/(3*Ta), 1/Ta_SO))
            B = incompressibleInertialLoads ? vcat(cnαUₙTQCdot*AC, Uₙ/Ta, fN/TfN, fM/TfM, fT/TfT, R/(3*Ta), R/Ta_SO) : vcat(cnαUₙTQCdot*AC, α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma, Uₙ/Ta, fN/TfN, fM/TfM, fT/TfT, R/(3*Ta), R/Ta_SO)
            # States
            χ[pitchPlungeStatesRange] = -A\B
        elseif typeof(solver) == BLo
            # Update kinematics
            BLo_motion_qualifiers!(element)
            # Unpack data
            @unpack bC,AC,bI = solver
            @unpack α = element.aero.flowAnglesAndRates
            @unpack Ωₐ = element.aero.flowVelocitiesAndRates
            @unpack α₀N,α1₀,cnα,f₀,fb,S1,S2,Tf₀,Tv₀,Tp,γbC = element.aero.airfoil.parametersBLo
            @unpack Tnα,TnM,Tnθ̇,Tmα,TmM,Tmθ̇ = element.aero.flowParameters
            # Steady separation point
            absα = abs(α)
            f = absα <= α1₀ ? 1-(1-fb)*exp((absα-α1₀)/S1) : f₀+(fb-f₀)*exp((α1₀-absα)/S2)
            # Time derivative of the product cnα * UₙTQC
            cnαUₙTQCdot = cnαUₙTQC_rate(element)
            # State matrices
            A = incompressibleInertialLoads ? -Diagonal(vcat(Θ*bC.*γbC, 1/Tp, 1/Tf₀, 1/(Tf₀/2), 1/Tv₀, 1/Tv₀)) : -Diagonal(vcat(Θ*bC.*γbC, 1/Tnα, 1/Tnθ̇, 1/(bI[1]*Tmα), 1/(bI[2]*Tmα), 1/Tmθ̇, 1/TnM, 1/(bI[1]*TmM), 1/(bI[2]*TmM), 1/Tp, 1/Tf₀, 1/(Tf₀/2), 1/Tv₀, 1/Tv₀))
            B = incompressibleInertialLoads ? vcat(cnαUₙTQCdot*AC, cnα*(α-α₀N)/Tp, f/Tf₀, f/(Tf₀/2), 0, 0) : vcat(cnαUₙTQCdot*AC, α, Ωₐ, α, α, Ωₐ, Ma, Ma, Ma, cnα*(α-α₀N)/Tp, f/Tf₀, f/(Tf₀/2), 0, 0)
            # States
            χ[pitchPlungeStatesRange] = -A\B
        end

        # Flap-induced flow states
        if !isnothing(flapStatesRange)
            # Unpack data
            @unpack δNow = element.aero
            @unpack ACf,bCf = element.aero.flapLoadsSolver
            # Get the rate of cnα*wFlap
            cnαwdotFlap = cnαwFlap_rate(element,δNow)
            # Flap states
            χ[flapStatesRange] = cnαwdotFlap*ACf./(Θ*bCf)
        end

        # Update states
        x[DOF_χ] = χ
    end

    @pack! problem = x
end