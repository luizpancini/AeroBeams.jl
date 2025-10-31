#
# mutable struct FlowParameters

#     FlowParameters composite type

# Fields
# - `Re`: Reynolds number
# - `Ma`: Mach number
# - `βₚ`: Prandtl-Glauert compressibility factor
# - `βₚ²`: Prandtl-Glauert compressibility factor squared
# - `Θ`: Uᵢ/b*βₚ², characteristic inverse time scale
# - `Tnα`: time scale for pitch-plunge-rate-induced inertial force
# - `TnM`: time scale for Mach-rate-induced normal inertial force
# - `Tnθ̇`: time scale for pitch-acceleration-induced normal inertial force
# - `Tmα`: time scale for pitch-plunge-rate-induced inertial moment
# - `TmM`: time scale for Mach-rate-induced normal inertial moment
# - `Tmθ̇`: time scale for pitch-acceleration-induced normal inertial moment
#
mutable struct FlowParameters
    
    # Fields
    Re
    Ma
    βₚ
    βₚ²
    Θ
    Tᵢ
    Tnα
    Tnθ̇
    TnM
    Tmα
    Tmθ̇
    TmM

    # Constructor
    function FlowParameters() 

        Re = 0.0
        Ma = 0.0
        βₚ = 1.0 
        βₚ² = 1.0
        Θ = 0.0
        Tᵢ = 1.0
        Tnα = 1.0
        Tnθ̇ = 1.0
        TnM = 1.0
        Tmα = 1.0
        Tmθ̇ = 1.0
        TmM = 1.0

        return new(Re,Ma,βₚ,βₚ²,Θ,Tᵢ,Tnα,Tnθ̇,TnM,Tmα,Tmθ̇,TmM)
    end
end


#
# mutable struct FlowAnglesAndRates

#     FlowAnglesAndRates composite type

# Fields
# - `α`: quasi-steady angle of attack
# - `β`: quasi-steady angle of sideslip
# - `αdot`: quasi-steady angle of attack rate
# - `αₑ`: effective (unsteady) angle of attack at 3/4-chord
#
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


#
# mutable struct FlowVelocitiesAndRates

#     FlowVelocitiesAndRates composite type

# Fields
# - `U`: wind velocity vector at spar position, resolved in basis W
# - `U∞`: norm of U, relative airspeed
# - `Uₛ`: spanwise component of U
# - `Uₜ`: tangential component of U
# - `Uₙ`: normal component of U
# - `Uᵢ`: in-plane component of U
# - `Ωₐ`: spanwise angular velocity component, resolved in basis W
# - `UₙMid`: relative normal wind velocity component at the airfoil's 1/2-chord
# - `UₙTQC`: relative normal wind velocity component at the airfoil's 3/4-chord
# - `Udot`: wind acceleration vector at spar position, resolved in basis W
# - `Uₜdot`: tangential component of Udot
# - `Uₙdot`: normal component of Udot
# - `Uᵢdot`: in-plane component of Udot
# - `Ωₐdot`: spanwise angular acce component, resolved in basis W
# - `UₙdotMid`: relative normal wind acceleration component at the airfoil's 1/2-chord
# - `UₙdotTQC`: = relative normal wind acceleration component at the airfoil's 3/4-chord
# - `UₜGust`: gust tangential velocity component
# - `UₙGust`: gust normal velocity component
# - `wₑp`: effective pitch-plunge-induced normal velocity component at the airfoil's 3/4-chord
#
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
    wₑp

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
        wₑp = 0.0

        return new(U,U∞,Uₛ,Uₜ,Uₙ,Uᵢ,Ωₐ,UₙMid,UₙTQC,Udot,Uₜdot,Uₙdot,Uᵢdot,Ωₐdot,UₙdotMid,UₙdotTQC,UₜGust,UₙGust,wₑp)
    end
end


#
# mutable struct AeroCoefficients

#     AeroCoefficients composite type

# Fields
# - `cn`: normal force coefficient
# - `cm`: pitching moment coefficient about spar position
# - `ct`: tangential force coefficient
# - `cnC`: circulatory component of cn
# - `cnI`: inertial component of cn
# - `cnP`: potential flow component of cn
# - `cmC`: circulatory component of cm
# - `cmI`: inertial component of cm
# - `cnF`: separated-flow circulatory component of cn
# - `cmF`: separated-flow circulatory component of cm
# - `ctC`: circulatory component of ct
# - `ctF`: separated-flow circulatory component of ct
# - `cnV`: DSV-induced component of cn
# - `cmV`: DSV-induced component of cm
# - `ctV`: DSV-induced component of ct
# - `cnNC`: Non-circulatory component of cn
# - `cmNC`: Non-circulatory component of cm
# - `ctNC`: Non-circulatory component of ct
#
mutable struct AeroCoefficients
    
    # Fields
    cn 
    cm 
    ct 
    cnC 
    cnI
    cnP 
    cmC 
    cmI 
    cnF 
    cmF
    ctC
    ctF
    cnV
    cmV
    ctV
    cnNC
    cmNC
    ctNC

    # Constructor
    function AeroCoefficients() 

        cn = 0.0
        cm = 0.0
        ct = 0.0
        cnC = 0.0
        cnI = 0.0
        cnP = 0.0
        cmC = 0.0
        cmI = 0.0
        cnF = 0.0
        cmF = 0.0
        ctC = 0.0
        ctF = 0.0
        cnV = 0.0
        cmV = 0.0
        ctV = 0.0
        cnNC = 0.0
        cmNC = 0.0
        ctNC = 0.0

        return new(cn,cm,ct,cnC,cnI,cnP,cmC,cmI,cnF,cmF,ctC,ctF,cnV,cmV,ctV,cnNC,cmNC,ctNC)
    end
end


#
# mutable struct BLiNamedStates

#     BLiNamedStates composite type

# Fields
# - `αlag`: lagged angle of attack
# - `f2primeN`: lagged separation point for cn
# - `f2primeM`: lagged separation point for cm
# - `f2primeT`: lagged separation point for ct
# - `RD`: lagged capped pitch rate ratio
# - `RD_stallOnsetRatio`: slightly lagged capped pitch rate ratio for stall onset
#
mutable struct BLiNamedStates
    
    # Fields
    αlag 
    f2primeN
    f2primeM
    f2primeT
    RD
    RD_stallOnsetRatio

    # Constructor
    function BLiNamedStates() 

        αlag = 0.0
        f2primeN = 1.0
        f2primeM = 1.0
        f2primeT = 1.0
        RD = 0.0
        RD_stallOnsetRatio = 0.0        

        return new(αlag,f2primeN,f2primeM,f2primeT,RD,RD_stallOnsetRatio)
    end
end


#
# mutable struct BLoNamedStates

#     BLoNamedStates composite type

# Fields
# - `cnPprime`: lagged cn
# - `f2Prime`: lagged separation point
# - `fPrimeM`: lagged separation point for cm on downstroke
# - `cnVP`: DSV-induced cn at positive angle of attack
# - `cnVN`: DSV-induced cn at negative angle of attack
# - `cnV`: sum of cnVP and cnVN
#
mutable struct BLoNamedStates
    
    # Fields
    cnPprime
    f2Prime
    fPrimeM
    cnVP
    cnVN
    cnV

    # Constructor
    function BLoNamedStates() 

        cnPprime = 0.0
        f2Prime = 1.0
        fPrimeM = 1.0
        cnVP = 0.0
        cnVN = 0.0
        cnV = 0.0     

        return new(cnPprime,f2Prime,fPrimeM,cnVP,cnVN,cnV)
    end
end


#
# mutable struct BLiKinematics

#     BLiKinematics composite type

# Fields
# - `q`: non-dimensional pitch rate
# - `r`: reduced non-dimensional pitch rate
# - `qR`: unsigned ratio of reduced pitch rate to critical pitch rate
# - `R`: unsigned capped reduced pitch rate ratio 
#
mutable struct BLiKinematics
    
    # Fields
    q
    r
    qR
    R

    # Constructor
    function BLiKinematics() 

        q = 0.0
        r = 0.0
        qR = 0.0
        R = 0.0      

        return new(q,r,qR,R)
    end
end


#
# mutable struct BLiFlowVariables

#     BLiFlowVariables composite type

# Fields
# - `stallOnsetRatio`
# - `upstroke`
# - `S`
# - `P`
# - `T`
# - `α1N`
# - `α1M`
# - `α1T`
# - `fN`
# - `fM`
# - `fT`
# - `fPrimeN`
# - `fPrimeM`
# - `fPrimeT`
# - `Ta_SO`
# - `TfN`
# - `TfM`
# - `TfT`
#
mutable struct BLiFlowVariables
    
    # Fields
    stallOnsetRatio
    upstroke
    S
    P
    T
    α1N
    α1M
    α1T
    fN
    fM
    fT
    fPrimeN
    fPrimeM
    fPrimeT
    Ta_SO
    TfN
    TfM
    TfT

    # Constructor
    function BLiFlowVariables() 

        stallOnsetRatio = 0.0
        upstroke = false
        S = false
        P = 0.0
        T = 0.0
        α1N = 0.0
        α1M = 0.0
        α1T = 0.0  
        fN = 1.0
        fM = 1.0
        fT = 1.0 
        fPrimeN = 1.0 
        fPrimeM = 1.0 
        fPrimeT = 1.0 
        Ta_SO = 1.0
        TfN = 1.0
        TfM = 1.0
        TfT = 1.0

        return new(stallOnsetRatio,upstroke,S,P,T,α1N,α1M,α1T,fN,fM,fT,fPrimeN,fPrimeM,fPrimeT,Ta_SO,TfN,TfM,TfT)
    end
end


#
# mutable struct BLoFlowVariables

#     BLoFlowVariables composite type

# Fields
# - `αlag`
# - `q`
# - `stallOnsetRatio`
# - `upstroke`
# - `α1`
# - `f`
# - `fPrime`
# - `Tf`
# - `Tv`
# - `Kf`
# - `cvdotP`
# - `cvdotN`
#
mutable struct BLoFlowVariables
    
    # Fields
    αlag
    q
    stallOnsetRatio
    upstroke
    α1
    f
    fPrime
    Tf
    Tv
    Kf
    cvdotP
    cvdotN

    # Constructor
    function BLoFlowVariables() 

        αlag = 0.0
        q = 0.0
        stallOnsetRatio = 0.0
        upstroke = false
        α1 = 0.0
        f = 1.0
        fPrime = 1.0
        Tf = 1.0
        Tv = 1.0
        Kf = 1.0
        cvdotP = 0.0
        cvdotN = 0.0

        return new(αlag,q,stallOnsetRatio,upstroke,α1,f,fPrime,Tf,Tv,Kf,cvdotP,cvdotN)
    end
end


#
# mutable struct BLiComplementaryVariables

#     BLiComplementaryVariables composite type

# Fields
# - `stallOnsetRatioPrev`
# - `αlagPrev`
# - `qRPrev`
# - `PPrev`
# - `upstrokePrev`
# - `maxStallOnsetRatio`
# - `minStallOnsetRatio` 
# - `qRmax`
# - `Ts`
# - `tv0P`
# - `fDiff_tv0P`
# - `qR_tv0P`
# - `R_tv0P`
# - `RD_tv0P`
# - `upstroke_tv0P`
# - `fDiff_tv0P2`
# - `RD_tv0P2`
# - `upstroke_tv0P2`
# - `tv0N`
# - `fDiff_tv0N`
# - `qR_tv0N`
# - `R_tv0N`
# - `RD_tv0N`
# - `upstroke_tv0N`
# - `fDiff_tv0N2`
# - `RD_tv0N2`
# - `upstroke_tv0N2`
# - `lastRD_tv0`
#
mutable struct BLiComplementaryVariables
    
    # Fields
    stallOnsetRatioPrev
    αlagPrev 
    qRPrev
    PPrev
    upstrokePrev 
    maxStallOnsetRatio
    minStallOnsetRatio 
    qRmax
    Ts
    tv0P
    fDiff_tv0P
    qR_tv0P
    R_tv0P
    RD_tv0P
    upstroke_tv0P
    τvP
    fDiff_tv0P2
    RD_tv0P2 
    upstroke_tv0P2
    tv0N
    fDiff_tv0N
    qR_tv0N
    R_tv0N
    RD_tv0N 
    upstroke_tv0N
    τvN
    fDiff_tv0N2
    RD_tv0N2  
    upstroke_tv0N2
    lastRD_tv0

    # Constructor
    function BLiComplementaryVariables() 

        stallOnsetRatioPrev = 0.0
        αlagPrev = 0.0
        qRPrev = 1.0
        PPrev = 0.0
        upstrokePrev = false  
        maxStallOnsetRatio = 1.0
        minStallOnsetRatio = 0.0 
        qRmax = 1.0
        Ts = 0.0
        tv0P = -Inf64
        fDiff_tv0P = 0.0
        qR_tv0P = 1.0
        R_tv0P = 1.0
        RD_tv0P = 1.0
        upstroke_tv0P = false
        τvP = Inf64
        fDiff_tv0P2 = -1.0
        RD_tv0P2 = 1.0
        upstroke_tv0P2 = false
        tv0N = -Inf64
        fDiff_tv0N = 0.0
        qR_tv0N = 1.0
        R_tv0N = 1.0
        RD_tv0N = 1.0
        upstroke_tv0N = false
        τvN = Inf64
        fDiff_tv0N2 = -1.0
        RD_tv0N2 = 1.0
        upstroke_tv0N2 = false
        lastRD_tv0 = 1.0

        return new(stallOnsetRatioPrev,αlagPrev,qRPrev,PPrev,upstrokePrev,maxStallOnsetRatio,minStallOnsetRatio,qRmax,Ts,tv0P,fDiff_tv0P,qR_tv0P,R_tv0P,RD_tv0P,upstroke_tv0P,τvP,fDiff_tv0P2,RD_tv0P2,upstroke_tv0P2,tv0N,fDiff_tv0N,qR_tv0N,R_tv0N,RD_tv0N,upstroke_tv0N,τvN,fDiff_tv0N2,RD_tv0N2,upstroke_tv0N2,lastRD_tv0)
    end
end


#
# mutable struct BLoComplementaryVariables

#     BLoComplementaryVariables composite type

# Fields
# - `stallOnsetRatioPrev`
# - `tv0P`
# - `tv0N`
# - `τvP`
# - `τvN`
#
mutable struct BLoComplementaryVariables
    
    # Fields
    stallOnsetRatioPrev
    tv0P
    tv0N
    τvP
    τvN

    # Constructor
    function BLoComplementaryVariables() 

        stallOnsetRatioPrev = 0.0
        tv0P = -Inf64
        tv0N = -Inf64
        τvP = Inf64
        τvN = Inf64

        return new(stallOnsetRatioPrev,tv0P,tv0N,τvP,τvN)
    end
end


#
# mutable struct AeroVariables

#     AeroVariables composite type

# Fields
# - `flowParameters::FlowParameters`
# - `flowAnglesAndRates::FlowAnglesAndRates`
# - `flowVelocitiesAndRates::FlowVelocitiesAndRates`
# - `aeroCoefficients::AeroCoefficients`
# - `BLiKin::BLiKinematics`
# - `BLiFlow::BLiFlowVariables`
# - `BLoFlow::BLoFlowVariables`
#
mutable struct AeroVariables
    
    # Fields
    flowParameters::FlowParameters
    flowAnglesAndRates::FlowAnglesAndRates
    flowVelocitiesAndRates::FlowVelocitiesAndRates
    aeroCoefficients::AeroCoefficients
    BLiKin::BLiKinematics
    BLiFlow::BLiFlowVariables
    BLoFlow::BLoFlowVariables

    # Constructor
    function AeroVariables(flowParameters,flowAnglesAndRates,flowVelocitiesAndRates,aeroCoefficients,BLiKin,BLiFlow,BLoFlow)
        return new(flowParameters,flowAnglesAndRates,flowVelocitiesAndRates,aeroCoefficients,BLiKin,BLiFlow,BLoFlow)
    end
end


# AeroProperties composite type: defines the aerodynamic properties of an element
@with_kw mutable struct AeroProperties

    # --- Immutable properties ---
    # Parent aerodynamic surface
    aeroSurface::AeroSurface
    # Geometry
    airfoil::Airfoil
    b::Real
    c::Real
    normSparPos::Float64
    aₕ::Float64
    Λ::Real
    φ::Real
    cosΛ::Real
    # Rotation tensors
    Rw::Matrix{Float64}
    RwT::Matrix{Float64}
    RwR0::Matrix{Float64}
    RwR0T::Matrix{Float64}
    # Flag for being flapped
    isFlapped::Bool
    # Theodorsen's flap constants
    Th::Union{Nothing,Vector{Float64}}
    # Tip loss correction factor as a function of the element's local coordinate (ζ)
    ϖ::Function
    # Total number of aerodynamic, flap and gust states, and respective ranges
    nTotalAeroStates::Int
    nFlapStates::Int
    nGustStates::Int
    pitchPlungeStatesRange::Union{Nothing,UnitRange{Int}}
    linearPitchPlungeStatesRange::Union{Nothing,UnitRange{Int}}
    nonlinearPitchPlungeStatesRange::Union{Nothing,UnitRange{Int}}
    circulatoryPitchPlungeStatesRange::Union{Nothing,UnitRange{Int}}
    inertialPitchPlungeStatesRange::Union{Nothing,UnitRange{Int}}
    flapStatesRange::Union{Nothing,UnitRange{Int}}
    gustStatesRange::Union{Nothing,UnitRange{Int}}
    linearGustStatesRange::Union{Nothing,UnitRange{Int}}
    nonlinearGustStatesRange::Union{Nothing,UnitRange{Int}}

    # --- Mutable properties ----
    # Tip loss correction factor at the element's midpoint
    ϖMid = ϖ(1/2)
    # State matrices
    A = zeros(nTotalAeroStates,nTotalAeroStates)
    B = zeros(nTotalAeroStates)
    # Nondimensional flow parameters
    flowParameters = FlowParameters()
    # Flow angles and rates
    flowAnglesAndRates = FlowAnglesAndRates()
    # Flow velocities and rates
    flowVelocitiesAndRates = FlowVelocitiesAndRates()
    # Aerodynamic coefficients
    aeroCoefficients = AeroCoefficients()
    # States of BL models
    BLiStates = BLiNamedStates()
    BLoStates = BLoNamedStates()
    # Flow kinematics of BL models
    BLiKin = BLiKinematics()
    # Flow variables of BL models
    BLiFlow = BLiFlowVariables()
    BLoFlow = BLoFlowVariables()
    # Complementary variables of BL models
    BLiCompVars = BLiComplementaryVariables()
    BLoCompVars = BLoComplementaryVariables()
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
    F_χ_χ::Matrix{Float64} = initial_F_χ_χ(aeroSurface.solver,nTotalAeroStates)
    # Aerodynamic derivatives w.r.t. elemental states' rates
    f1χ_Vdot::Matrix{Float64} = zeros(3,3)
    f2χ_Vdot::Matrix{Float64} = zeros(3,3)
    f1χ_Ωdot::Matrix{Float64} = zeros(3,3)
    f2χ_Ωdot::Matrix{Float64} = zeros(3,3)
    m1χ_Vdot::Matrix{Float64} = zeros(3,3)
    m2χ_Vdot::Matrix{Float64} = zeros(3,3)
    m1χ_Ωdot::Matrix{Float64} = zeros(3,3)
    m2χ_Ωdot::Matrix{Float64} = zeros(3,3)
    F_χ_Vdot::Matrix{Float64} = initial_F_χ_Vdot(aeroSurface.solver,nTotalAeroStates,circulatoryPitchPlungeStatesRange,airfoil.attachedFlowParameters.cnα)
    F_χ_Ωdot::Matrix{Float64} = initial_F_χ_Ωdot(aeroSurface.solver,nTotalAeroStates,circulatoryPitchPlungeStatesRange,c,normSparPos,airfoil.attachedFlowParameters.cnα)
    F_χ_χdot::Matrix{Float64} = Matrix(1.0*LinearAlgebra.I,nTotalAeroStates,nTotalAeroStates)

end

# AeroProperties constructor
function AeroProperties(aeroSurface::AeroSurface,rotationParametrization::String,R0::Matrix{Float64},x1::Real,x1_norm::Float64,x1_n1_norm::Float64,x1_n2_norm::Float64)

    # Possibly spanwise-variable geometric properties 
    airfoil = aeroSurface.airfoil
    c = aeroSurface.c isa Real ? aeroSurface.c : aeroSurface.c(x1)
    b = c/2
    normSparPos = aeroSurface.normSparPos isa Real ? aeroSurface.normSparPos : aeroSurface.normSparPos(x1)
    aₕ = 2*normSparPos .- 1
    Λ = aeroSurface.Λ isa Real ? aeroSurface.Λ : aeroSurface.Λ(x1)
    φ = aeroSurface.φ isa Real ? aeroSurface.φ : aeroSurface.φ(x1)
    cosΛ = cos(Λ)

    # Rotation tensors
    if rotationParametrization == "E321"
        Rw = rotation_tensor_E321([Λ; 0; φ])
    elseif rotationParametrization == "E213"
        Rw = rotation_tensor_E213([0; φ; Λ])
    elseif rotationParametrization == "E231"
        Rw = rotation_tensor_E231([0; Λ; φ])
    end
    RwT = Matrix(Rw')
    RwR0 = Rw*R0
    RwR0T = Matrix(RwR0')

    # Flag for being flapped (element's midpoint is inside the flapped portion of span)
    isFlapped = !isnothing(aeroSurface.normFlapSpan) && minimum(aeroSurface.normFlapSpan) < x1_norm < maximum(aeroSurface.normFlapSpan)

    # Set Theodorsen's flap constants, if applicable
    Th = isFlapped ? theodorsen_flap_constants(normSparPos,aeroSurface.normFlapPos) : nothing

    # Set tip loss correction variables
    hasTipCorrection = aeroSurface.hasTipCorrection
    tipLossFunction = aeroSurface.tipLossFunction
    tipLossDecayFactor = aeroSurface.tipLossDecayFactor
    if !hasTipCorrection
        # No tip correction -> ϖ is constant and equal to one
        ϖ = ζ -> 1
    else
        # Arclength coordinates
        s = ζ -> x1_n1_norm + ζ * (x1_n2_norm-x1_n1_norm)
        sFlip = ζ -> (1-x1_n1_norm) + ζ * ((1-x1_n2_norm)-(1-x1_n1_norm))
        if isnothing(tipLossFunction)
            # Set standard tip loss function with given tipLossDecayFactor
            ϖ = tipLossDecayFactor >= 0 ? ζ -> 1-exp(-tipLossDecayFactor*(1-s(ζ))) : ϖ = ζ -> 1-exp(tipLossDecayFactor*(1-sFlip(ζ)))
        else
            if aeroSurface.tipLossFunctionIsAirspeedDependent
                # Initialize tip loss function with airspeed assumed as zero (updated later on model creation)
                ϖ = ζ -> tipLossFunction(0)(s(ζ))
            else
                # Set given tip loss function
                ϖ = ζ -> tipLossFunction(s(ζ))
            end
        end
    end

    # Flag for flap states (solver is thin-airfoil theory, and flap deflection is either an user input or a trim variable)
    hasFlapStates = typeof(aeroSurface.flapLoadsSolver) == ThinAirfoilTheory && (aeroSurface.δIsInput || aeroSurface.δIsTrimVariable)

    # Total number of aerodynamic states (assume no gust states are active, update later upon model creation)
    solver = aeroSurface.solver
    nPitchPungeStates = solver.nStates
    nFlapStates = hasFlapStates ? aeroSurface.flapLoadsSolver.nStates : 0
    nGustStates = 0
    nTotalAeroStates = nPitchPungeStates + nFlapStates + nGustStates

    # Set aerodynamic states' ranges (assume no gust states are active, update later upon model creation)
    if typeof(solver) == QuasiSteady
        circulatoryPitchPlungeStatesRange = inertialPitchPlungeStatesRange = pitchPlungeStatesRange = linearPitchPlungeStatesRange = nonlinearPitchPlungeStatesRange = nothing
    elseif typeof(solver) == Indicial
        circulatoryPitchPlungeStatesRange = 1:solver.nCirculatoryStates
        inertialPitchPlungeStatesRange = solver.nCirculatoryStates+1:solver.nStates
        pitchPlungeStatesRange = linearPitchPlungeStatesRange = 1:solver.nStates
        nonlinearPitchPlungeStatesRange = nothing
    elseif typeof(solver) == Inflow
        circulatoryPitchPlungeStatesRange = 1:solver.nInflowStates
        inertialPitchPlungeStatesRange = solver.nInflowStates+1:solver.nStates
        pitchPlungeStatesRange = linearPitchPlungeStatesRange = 1:solver.nStates
        nonlinearPitchPlungeStatesRange = nothing
    elseif typeof(solver) in [BLi,BLo]
        circulatoryPitchPlungeStatesRange = 1:solver.nLinearCirculatoryStates
        inertialPitchPlungeStatesRange = solver.nLinearCirculatoryStates+1:solver.nLinearCirculatoryStates+solver.nInertialStates
        pitchPlungeStatesRange = 1:solver.nStates
        linearPitchPlungeStatesRange = 1:solver.nLinearCirculatoryStates+solver.nInertialStates
        nonlinearPitchPlungeStatesRange = solver.nLinearCirculatoryStates+solver.nInertialStates+1:solver.nStates
    end
    flapStatesRange = hasFlapStates ? (nPitchPungeStates+1:nPitchPungeStates+nFlapStates) : nothing
    nonlinearGustStatesRange = linearGustStatesRange = gustStatesRange = nothing

    return AeroProperties(aeroSurface=aeroSurface,airfoil=airfoil,b=b,c=c,normSparPos=normSparPos,aₕ=aₕ,Λ=Λ,φ=φ,cosΛ=cosΛ,Rw=Rw,RwT=RwT,RwR0=RwR0,RwR0T=RwR0T,isFlapped=isFlapped,Th=Th,ϖ=ϖ,nTotalAeroStates=nTotalAeroStates,nFlapStates=nFlapStates,nGustStates=nGustStates,pitchPlungeStatesRange=pitchPlungeStatesRange,linearPitchPlungeStatesRange=linearPitchPlungeStatesRange,nonlinearPitchPlungeStatesRange=nonlinearPitchPlungeStatesRange,circulatoryPitchPlungeStatesRange=circulatoryPitchPlungeStatesRange,inertialPitchPlungeStatesRange=inertialPitchPlungeStatesRange,flapStatesRange=flapStatesRange,gustStatesRange=gustStatesRange,linearGustStatesRange=linearGustStatesRange,nonlinearGustStatesRange=nonlinearGustStatesRange)
end


# Computes Theodorsen's flap constants
function theodorsen_flap_constants(normSparPos::Float64,normFlapPos::Float64)

    # Semichord-normalized spar position after midchord
    a = 2*normSparPos-1

    # Semichord-normalized flap hinge position after midchord
    d = 2*normFlapPos-1

    # Theodorsen's constants
    Th = zeros(18)
    Th[1] = -1/3*sqrt(1-d^2)*(2+d^2)+d*acos(d)
    Th[2] = d*(1-d^2)-sqrt(1-d^2)*(1+d^2)*acos(d)+d*acos(d)^2;
    Th[3] = -(1/8+d^2)*acos(d)^2+1/4*d*sqrt(1-d^2)*acos(d)*(7+2*d^2)-1/8*(1-d^2)*(5*d^2+4)
    Th[4] = -acos(d)+d*sqrt(1-d^2)
    Th[5] = -(1-d^2)-acos(d)^2+2*d*sqrt(1-d^2)*acos(d)
    Th[6] = Th[2]
    Th[7] = -(1/8+d^2)*acos(d)+1/8*d*sqrt(1-d^2)*(7+2*d^2)
    Th[8] = -1/3*sqrt(1-d^2)*(2*d^2+1)+d*acos(d)
    Th[9] = 1/2*(1/3*sqrt(1-d^2)^3+a*Th[4])
    Th[10] = sqrt(1-d^2)+acos(d)
    Th[11] = acos(d)*(1-2*d)+sqrt(1-d^2)*(2-d)
    Th[12] = sqrt(1-d^2)*(2+d)-acos(d)*(2*d+1)
    Th[13] = 1/2*(-Th[7]-(d-a)*Th[1])

    # Useful functions of the above
    Th[14] = Th[4]+Th[10]
    Th[15] = Th[1]-Th[8]-Th[4]*(d-a)+Th[11]/2
    Th[16] = Th[7]+Th[1]*(d-a)
    Th[17] = Th[10]/π
    Th[18] = Th[11]/(2π)

    return Th

end


# Computes the initial value of F_χ_χ (for zero relative airspeed) - theoretically, the ϵ value in the function should go to zero as the airspeed goes to zero, but numerical problems are found for the inflow solver
function initial_F_χ_χ(solver::AeroSolver,nStates::Int)

    # Calculate according to solver
    if typeof(solver) == QuasiSteady
        ϵ = 0.0 
    elseif typeof(solver) in [Indicial,BLi,BLo]  
        ϵ = 1e-4
    elseif typeof(solver) == Inflow
        ϵ = 1.0
    end

    return Matrix(-ϵ*LinearAlgebra.I,nStates,nStates)
end


# Computes the initial value of F_χ_Vdot (for zero relative airspeed)
function initial_F_χ_Vdot(solver::AeroSolver,nStates::Int,circulatoryPitchPlungeStatesRange::Union{Nothing,UnitRange{Int}},cnα::Real)

    F_χ_Vdot = zeros(nStates,3)

    # Calculate according to solver
    if typeof(solver) in [Indicial,BLi,BLo]
        F_χ_Vdot[circulatoryPitchPlungeStatesRange[1:length(solver.AC)],3] = cnα*solver.AC
    elseif typeof(solver) == Inflow
        F_χ_Vdot[circulatoryPitchPlungeStatesRange,3] = solver.AₚInvcₚ
    end

    return F_χ_Vdot
end


# Computes the initial value of F_χ_Ωdot (for zero relative airspeed)
function initial_F_χ_Ωdot(solver::AeroSolver,nStates::Int,circulatoryPitchPlungeStatesRange::Union{Nothing,UnitRange{Int}},c::Real,normSparPos::Float64,cnα::Real)

    F_χ_Ωdot = zeros(nStates,3)

    # Calculate according to solver
    if typeof(solver) in [Indicial,BLi,BLo]
        F_χ_Ωdot[circulatoryPitchPlungeStatesRange[1:length(solver.AC)],1] = c*(normSparPos-3/4)*cnα*solver.AC
    elseif typeof(solver) == Inflow
        F_χ_Ωdot[circulatoryPitchPlungeStatesRange,1] = c*(normSparPos-3/4)*solver.AₚInvcₚ  
    end

    return F_χ_Ωdot
end