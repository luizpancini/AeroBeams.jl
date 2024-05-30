"""
mutable struct FlowAnglesAndRates

    FlowAnglesAndRates composite type

# Fields
- α = quasi-steady angle of attack
- β = angle of sideslip
- αdot = quasi-steady angle of attack rate
- αₑ = effective angle of attack at 3/4-chord
"""
mutable struct FlowAnglesAndRates
    
    # Fields
    α
    β
    αdot
    αₑ

    # Constructor
    function FlowAnglesAndRates() 

        α = NaN64
        β = NaN64
        αdot = 0.0
        αₑ = NaN64

        return new(α,β,αdot,αₑ)
    end
end


"""
mutable struct FlowVelocitiesAndRates

    FlowVelocitiesAndRates composite type

# Fields
- U = wind velocity vector at spar position, resolved in basis W
- U∞ = norm of U
- Uₛ = spanwise component of U
- Uₜ = tangential component of U
- Uₙ = normal component of U
- Uᵢ = in-plane component of U
- Ωₐ = spanwise angular velocity component, resolved in basis W
- UₙMid = relative normal wind velocity component at the airfoil's 1/2-chord
- UₙTQC = relative normal wind velocity component at the airfoil's 3/4-chord
- Udot = wind acceleration vector at spar position, resolved in basis W
- Uₜdot = tangential component of Udot
- Uₙdot = normal component of U
- Uᵢdot = in-plane component of U
- Ωₐdot = spanwise angular velocity component, resolved in basis W
- UₙdotMid = relative normal wind acceleration component at the airfoil's 1/2-chord
- UₙdotTQC = = relative normal wind acceleration component at the airfoil's 3/4-chord
- UₜGust = gust tangential velocity component
- UₙGust = gust normal velocity component
"""
mutable struct FlowVelocitiesAndRates
    
    # Fields
    U
    U∞
    Uₛ
    Uₜ
    Uₙ
    Uᵢ
    Ωₐ
    UₙMid
    UₙTQC
    Udot
    Uₜdot
    Uₙdot
    Uᵢdot
    Ωₐdot
    UₙdotMid
    UₙdotTQC
    UₜGust
    UₙGust

    # Constructor
    function FlowVelocitiesAndRates()

        U = zeros(3)
        U∞ = 0.0
        Uₛ = 0.0
        Uₜ = 0.0
        Uₙ = 0.0
        Uᵢ = 0.0
        Ωₐ = 0.0
        UₙMid = 0.0
        UₙTQC = 0.0
        Udot = zeros(3)
        Uₜdot = 0.0
        Uₙdot = 0.0
        Uᵢdot = 0.0
        Ωₐdot = 0.0
        UₙdotMid = 0.0
        UₙdotTQC = 0.0
        UₜGust = 0.0
        UₙGust = 0.0 

        return new(U,U∞,Uₛ,Uₜ,Uₙ,Uᵢ,Ωₐ,UₙMid,UₙTQC,Udot,Uₜdot,Uₙdot,Uᵢdot,Ωₐdot,UₙdotMid,UₙdotTQC,UₜGust,UₙGust)
    end
end


"""
mutable struct AeroCoefficients

    AeroCoefficients composite type

# Fields
- cn = normal force coefficient
- cm = pitching moment coefficient about spar position
- ct = tangential force coefficient
- cnC = circulatory component of cn
- cnI = inertial component of cn
- cmC = circulatory component of cm
- cmI = inertial component of cm
- cmRot = rotation-induced component of cm
"""
mutable struct AeroCoefficients
    
    # Fields
    cn 
    cm 
    ct 
    cnC 
    cnI 
    cmC 
    cmI 
    cmRot 

    # Constructor
    function AeroCoefficients() 

        cn = 0.0
        cm = 0.0
        ct = 0.0
        cnC = 0.0
        cnI = 0.0
        cmC = 0.0
        cmI = 0.0
        cmRot = 0.0

        return new(cn,cm,ct,cnC,cnI,cmC,cmI,cmRot)
    end
end


"""
mutable struct FlowVariables

    FlowVariables composite type

# Fields
- flowAnglesAndRates
- flowVelocitiesAndRates
- aeroCoefficients
"""
mutable struct FlowVariables
    
    # Fields
    α
    β
    αdot
    αₑ
    U
    U∞
    Uₛ
    Uₜ
    Uₙ
    Uᵢ
    Ωₐ
    UₙMid
    UₙTQC
    Udot
    Uₜdot
    Uₙdot
    Uᵢdot
    Ωₐdot
    UₙdotMid
    UₙdotTQC
    UₜGust
    UₙGust
    cn 
    cm 
    ct 
    cnC 
    cnI 
    cmC 
    cmI 
    cmRot

    # Constructor
    function FlowVariables(flowAnglesAndRates,flowVelocitiesAndRates,aeroCoefficients)
        
        @unpack α,β,αdot,αₑ = flowAnglesAndRates
        @unpack U,U∞,Uₛ,Uₜ,Uₙ,Uᵢ,Ωₐ,UₙMid,UₙTQC,Udot,Uₜdot,Uₙdot,Uᵢdot,Ωₐdot,UₙdotMid,UₙdotTQC,UₜGust,UₙGust = flowVelocitiesAndRates
        @unpack cn,cm,ct,cnC,cnI,cmC,cmI,cmRot = aeroCoefficients

        return new(α,β,αdot,αₑ,U,U∞,Uₛ,Uₜ,Uₙ,Uᵢ,Ωₐ,UₙMid,UₙTQC,Udot,Uₜdot,Uₙdot,Uᵢdot,Ωₐdot,UₙdotMid,UₙdotTQC,UₜGust,UₙGust,cn,cm,ct,cnC,cnI,cmC,cmI,cmRot)
    end
end


"""
@with_kw mutable struct AeroProperties

    AeroProperties composite type

# Fields

"""
@with_kw mutable struct AeroProperties

    # Aerodynamic solvers 
    solver::AeroSolver
    flapLoadsSolver::FlapAeroSolver
    gustLoadsSolver::GustAeroSolver
    # Total number of aerodynamic, flap and gust states, and respective ranges
    nTotalAeroStates::Int64
    nFlapStates::Int64
    nGustStates::Int64
    pitchPlungeStatesRange::UnitRange{Int64}
    flapStatesRange::Union{Nothing,UnitRange{Int64}}
    gustStatesRange::Union{Nothing,UnitRange{Int64}}
    # Aerodynamic derivatives calculation method
    derivationMethod::DerivationMethod
    # Geometry
    airfoil::Airfoil
    b::Number
    c::Number
    normSparPos::Float64
    Λ::Number
    Rw::Matrix{Float64}
    RwT::Matrix{Float64}
    RwR0::Matrix{Float64}
    RwR0T::Matrix{Float64}
    flapSiteID::Int64
    normFlapPos::Float64
    flapped::Bool
    # Flap deflection trim TF, values and its rates as functions of time
    δIsTrimVariable::Bool
    δ::Function
    δdot::Function
    δddot::Function
    # Current values of flap deflection and its rates
    δNow::Number
    δdotNow::Number
    δddotNow::Number
    # Flap deflection multiplier for slave surfaces
    δMultiplier::Number
    # TF to update airfoil parameters according to flow parameters
    updateAirfoilParameters::Bool
    # Tip loss correction factor as a function of the element's local coordinate (ζ), and value at midpoint 
    ϖ::Function
    ϖMid::Number = ϖ(1/2)
    # State matrices
    A = zeros(nTotalAeroStates,nTotalAeroStates)
    B = zeros(nTotalAeroStates)
    # Nondimensional flow parameters
    Re = 0.0
    Ma = 0.0
    βₚ = 1.0 
    βₚ² = 1.0
    # Flow angles and rates
    flowAnglesAndRates = FlowAnglesAndRates()
    # Flow velocities and rates
    flowVelocitiesAndRates = FlowVelocitiesAndRates()
    # Aerodynamic coefficients
    aeroCoefficients = AeroCoefficients()
    # Nodal aerodynamic loads resultants array
    F = zeros(12)
    # Aerodynamic derivatives w.r.t. elemental states
    f1χ_V::Matrix{Float64} = zeros(3,3)
    f2χ_V::Matrix{Float64} = zeros(3,3)
    f1χ_Ω::Matrix{Float64} = zeros(3,3)
    f2χ_Ω::Matrix{Float64} = zeros(3,3)
    m1χ_V::Matrix{Float64} = zeros(3,3)
    m2χ_V::Matrix{Float64} = zeros(3,3)
    m1χ_Ω::Matrix{Float64} = zeros(3,3)
    m2χ_Ω::Matrix{Float64} = zeros(3,3)
    f1χ_χ::Matrix{Float64} = zeros(3,nTotalAeroStates)
    f2χ_χ::Matrix{Float64} = zeros(3,nTotalAeroStates)
    m1χ_χ::Matrix{Float64} = zeros(3,nTotalAeroStates)
    m2χ_χ::Matrix{Float64} = zeros(3,nTotalAeroStates)
    f1χ_δ::Vector{Float64} = zeros(3)
    f2χ_δ::Vector{Float64} = zeros(3)
    m1χ_δ::Vector{Float64} = zeros(3)
    m2χ_δ::Vector{Float64} = zeros(3)
    F_χ_V::Matrix{Float64} = zeros(nTotalAeroStates,3)
    F_χ_Ω::Matrix{Float64} = zeros(nTotalAeroStates,3)
    F_χ_χ::Matrix{Float64} = initial_F_χ_χ(solver,nTotalAeroStates)
    # Aerodynamic derivatives w.r.t. elemental states' rates
    f1χ_Vdot::Matrix{Float64} = zeros(3,3)
    f2χ_Vdot::Matrix{Float64} = zeros(3,3)
    f1χ_Ωdot::Matrix{Float64} = zeros(3,3)
    f2χ_Ωdot::Matrix{Float64} = zeros(3,3)
    m1χ_Vdot::Matrix{Float64} = zeros(3,3)
    m2χ_Vdot::Matrix{Float64} = zeros(3,3)
    m1χ_Ωdot::Matrix{Float64} = zeros(3,3)
    m2χ_Ωdot::Matrix{Float64} = zeros(3,3)
    F_χ_Vdot::Matrix{Float64} = initial_F_χ_Vdot(solver,nTotalAeroStates,pitchPlungeStatesRange,airfoil.attachedFlowParameters.cnα)
    F_χ_Ωdot::Matrix{Float64} = initial_F_χ_Ωdot(solver,nTotalAeroStates,pitchPlungeStatesRange,c,normSparPos,airfoil.attachedFlowParameters.cnα)
    F_χ_χdot::Matrix{Float64} = Matrix(1.0*LinearAlgebra.I,nTotalAeroStates,nTotalAeroStates)

end

# Constructor
function AeroProperties(aeroSurface::AeroSurface,R0::Matrix{Float64},x1::Number,x1_norm::Float64,x1_n1_norm::Float64,x1_n2_norm::Float64)

    # Set geometric properties 
    airfoil = aeroSurface.airfoil
    c = aeroSurface.c isa Number ? aeroSurface.c : aeroSurface.c(x1)
    b = c/2
    normSparPos = aeroSurface.normSparPos isa Number ? aeroSurface.normSparPos : aeroSurface.normSparPos(x1)
    Λ = aeroSurface.Λ isa Number ? aeroSurface.Λ : aeroSurface.Λ(x1)
    Rw = rotation_tensor_E321([-Λ; 0; -airfoil.attachedFlowParameters.α₀N])
    RwT = Matrix(Rw')
    RwR0 = Rw*R0
    RwR0T = Matrix(RwR0')
    flapSiteID = !isnothing(aeroSurface.flapSiteID) ? aeroSurface.flapSiteID : 100
    if !isnothing(aeroSurface.normFlapSpan)
        normFlapPos = minimum(aeroSurface.normFlapSpan) < x1_norm < maximum(aeroSurface.normFlapSpan) ? aeroSurface.normFlapPos : 1.0
    else
        normFlapPos = 1.0
    end
    flapped = normFlapPos < 1 ? true : false

    # Set aerodynamic solvers and number of aerodynamic states
    solver = aeroSurface.solver
    flapLoadsSolver = aeroSurface.flapLoadsSolver
    gustLoadsSolver = aeroSurface.gustLoadsSolver

    # TF for flap states (the element is flapped, solver is thin-airfoil theory, and flap deflection is either an user input or a trim variable)
    hasFlapStates = flapped && typeof(flapLoadsSolver) == ThinAirfoilTheory && (aeroSurface.δIsInput || aeroSurface.δIsTrimVariable) ? true : false

    # Update Theodorsen's flap constants, if applicable
    if hasFlapStates
        flapLoadsSolver.Th = TheodorsenFlapConstants(normSparPos,normFlapPos)
    end

    # Total number of aerodynamic states (assume no gust states are active, update later upon model creation)
    nGustStates = 0
    nFlapStates = hasFlapStates ? flapLoadsSolver.nStates : 0
    nTotalAeroStates = solver.nStates + nFlapStates + nGustStates

    # Set aerodynamic states' ranges
    pitchPlungeStatesRange = 1:solver.nStates
    flapStatesRange = hasFlapStates ? (solver.nStates+1:solver.nStates+flapLoadsSolver.nStates) : nothing
    gustStatesRange = nothing

    # Set aerodynamic derivatives calculation method
    derivationMethod = aeroSurface.derivationMethod

    # Set flap deflection trim TF, value and rates
    δIsTrimVariable = aeroSurface.δIsTrimVariable
    δ = aeroSurface.δ
    δdot = aeroSurface.δdot
    δddot = aeroSurface.δddot
    δNow = δ(0)
    δdotNow = δdot(0)
    δddotNow = δddot(0)

    # Initialize flap deflection multiplier for slave surfaces (updated later on model assembly)
    δMultiplier = 1.0

    # Set TF to update airfoil parameters according to flow parameters
    updateAirfoilParameters = aeroSurface.updateAirfoilParameters

    # Set tip loss correction variables
    hasTipCorrection = aeroSurface.hasTipCorrection
    tipLossFunction = aeroSurface.tipLossFunction
    tipLossDecayFactor = aeroSurface.tipLossDecayFactor
    if !hasTipCorrection
        ϖ = ζ -> 1
    else
        if isnothing(tipLossFunction)
            ϖ = tipLossDecayFactor >= 0 ? ζ -> 1-exp(-tipLossDecayFactor*(1-(x1_n1_norm+ζ*(x1_n2_norm-x1_n1_norm)))) : ϖ = ζ -> 1-exp(-tipLossDecayFactor*(x1_n1_norm+ζ*(x1_n2_norm-x1_n1_norm)))
        else
            ϖ = ζ -> tipLossFunction(ζ)
        end
    end

    return AeroProperties(solver=solver,flapLoadsSolver=flapLoadsSolver,gustLoadsSolver=gustLoadsSolver,nTotalAeroStates=nTotalAeroStates,nFlapStates=nFlapStates,nGustStates=nGustStates,pitchPlungeStatesRange=pitchPlungeStatesRange,flapStatesRange=flapStatesRange,gustStatesRange=gustStatesRange,derivationMethod=derivationMethod,airfoil=airfoil,b=b,c=c,normSparPos=normSparPos,Λ=Λ,Rw=Rw,RwT=RwT,RwR0=RwR0,RwR0T=RwR0T,flapSiteID=flapSiteID,normFlapPos=normFlapPos,flapped=flapped,δIsTrimVariable=δIsTrimVariable,δ=δ,δdot=δdot,δddot=δddot,δNow=δNow,δdotNow=δdotNow,δddotNow=δddotNow,δMultiplier=δMultiplier,updateAirfoilParameters=updateAirfoilParameters,ϖ=ϖ)
end


"""
initial_F_χ_χ(solver::AeroSolver,nStates::Int64)

Computes the initial value of F_χ_χ (for zero relative airspeed) - theoretically, the ϵ value in the function should go to zero as the airspeed goes to zero, but numerical problems are found for the inflow solver

# Arguments
- solver::AeroSolver
- nStates::Int64
"""
function initial_F_χ_χ(solver::AeroSolver,nStates::Int64)

    # Calculate according to solver
    if typeof(solver) == QuasiSteady
        ϵ = 0.0 
    elseif typeof(solver) == Indicial  
        ϵ = 1e-4
    elseif typeof(solver) == Inflow
        ϵ = 1.0
    end

    return Matrix(-ϵ*LinearAlgebra.I,nStates,nStates)
end


"""
initial_F_χ_Vdot(solver::AeroSolver,nStates::Int64,pitchPlungeStatesRange::UnitRange{Int64},cnα::Number)

Computes the initial value of F_χ_Vdot (for zero relative airspeed)

# Arguments
- solver::AeroSolver
- nStates::Int64
- pitchPlungeStatesRange::UnitRange{Int64}
- cnα::Number
"""
function initial_F_χ_Vdot(solver::AeroSolver,nStates::Int64,pitchPlungeStatesRange::UnitRange{Int64},cnα::Number)

    F_χ_Vdot = zeros(nStates,3)

    # Calculate according to solver
    if typeof(solver) == Indicial  
        F_χ_Vdot[pitchPlungeStatesRange,3] = cnα*solver.AW
    elseif typeof(solver) == Inflow
        F_χ_Vdot[pitchPlungeStatesRange,3] = solver.AₚInvcₚ
    end

    return F_χ_Vdot
end


"""
initial_F_χ_Ωdot(solver::AeroSolver,nStates::Int64,pitchPlungeStatesRange::UnitRange{Int64},c::Number,normSparPos::Float64,cnα::Number)

Computes the initial value of F_χ_Ωdot (for zero relative airspeed)

# Arguments
- solver::AeroSolver
- nStates::Int64
- c::Number
- normSparPos::Float64
- cnα::Number
"""
function initial_F_χ_Ωdot(solver::AeroSolver,nStates::Int64,pitchPlungeStatesRange::UnitRange{Int64},c::Number,normSparPos::Float64,cnα::Number)

    F_χ_Ωdot = zeros(nStates,3)

    # Calculate according to solver
    if typeof(solver) == Indicial  
        F_χ_Ωdot[pitchPlungeStatesRange,1] = c*(normSparPos-3/4)*cnα*solver.AW
    elseif typeof(solver) == Inflow
        F_χ_Ωdot[pitchPlungeStatesRange,1] = c*(normSparPos-3/4)*solver.AₚInvcₚ
    end

    return F_χ_Ωdot
end