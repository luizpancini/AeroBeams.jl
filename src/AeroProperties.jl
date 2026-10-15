## AEROPROPERTIES COMPOSITE TYPES
# ------------------------------------------------------------------------------
# AeroGeometry composite type
# Fields
# - `b`: semichord length of the element
# - `c`: chord length of the element
# - `normSparPos`: normalized spar position of the element
# - `aₕ`: normalized spar position in the range [-1,1]
# - `Λ`: sweep angle of the element
# - `φ`: twist angle of the element
# - `cosΛ`: cosine of the sweep angle
# - `cosΛ²`: cosine squared of the sweep angle
# - `csecΛ`: chord times secant of the sweep angle
# - `Rw`: rotation tensor from the undeformed beam basis b to the undeformed aerodynamic basis w
# - `RwT`: transpose of Rw
# - `RwR0`: rotation tensor from the undeformed beam basis b to the undeformed aerodynamic basis w, resolved in basis A
# - `isFlapped`: flag for being flapped (element's midpoint is inside the flapped portion of span)
# - `Th`: Theodorsen's flap constants
#
struct AeroGeometry

    b::Float64
    c::Float64
    normSparPos::Float64
    aₕ::Float64
    Λ::Float64
    φ::Float64

    cosΛ::Float64
    cosΛ²::Float64
    csecΛ::Float64

    Rw::SMatrix{3,3,Float64,9}
    RwT::SMatrix{3,3,Float64,9}
    RwR0::SMatrix{3,3,Float64,9}

    isFlapped::Bool
    Th::Union{Nothing,SVector{18,Float64}}
    
end


# StateRanges composite type
# Fields
# - `flapStatesRange`: range of flap states
# - `gustStatesRange`: range of gust states
#
struct StateRanges

    flapStatesRange::UnitRange{Int}
    gustStatesRange::UnitRange{Int}

end


# AeroLoadsDistribution composite type
# Fields
# - `ϖ`: tip loss function in the element's local coordinate (ζ)
# - `ϖMid`: tip loss correction factor at the element's midpoint
# - `Uref`: reference wind speed for the element's tip loss correction factor
#
mutable struct AeroLoadsDistribution{F<:Function}

    ϖ::F
    ϖMid::Float64
    Uref::Float64

end


# FlowParameters composite type
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
mutable struct FlowParameters{T<:Real}

    Re::T
    Ma::T
    βₚ::T
    βₚ²::T
    Θ::T
    Tᵢ::T

    Tnα::T
    Tnθ̇::T
    TnM::T
    Tmα::T
    Tmθ̇::T
    TmM::T

end
FlowParameters{T}() where {T<:Real} = FlowParameters{T}(
    zero(T), zero(T), one(T), one(T), zero(T), one(T),
    one(T), one(T), one(T), one(T), one(T), one(T), 
)
FlowParametersNaN() = FlowParameters{Float64}(
    ntuple(_ -> NaN, 12)...
)
Base.copy(x::FlowParameters{T}) where {T} = copy_fields(x)


# SteadyKinematics composite type
# Fields
# - `U`: wind velocity vector at spar position, resolved in basis W
# - `Uₛ`: spanwise component of U
# - `Uₜ`: tangential component of U
# - `Uₙ`: normal component of U
# - `Uᵢ`: in-plane component of U
# - `UₙMid`: normal wind velocity component at the airfoil's 1/2-chord
# - `UₙTQC`: normal wind velocity component at the airfoil's 3/4-chord
# - `UGust`: gust velocity vector at spar position, resolved in basis W
# - `UₜGust`: tangential component of UGust
# - `UₙGust`: normal component of UGust
# - `Ωₐ`: spanwise angular velocity component
# - `α`: quasi-steady angle of attack
# - `β`: quasi-steady angle of sideslip
# - `wf`: quasi-steady flap-induced normalwash
#
mutable struct SteadyKinematics{T<:Real}
        
    U::MVector{3,T}
    Uₛ::T
    Uₜ::T
    Uₙ::T
    Uᵢ::T
    UₙMid::T
    UₙTQC::T

    UGust::MVector{3,T}
    UₜGust::T
    UₙGust::T

    Ωₐ::T

    α::T
    β::T

    wf::T

end
SteadyKinematics{T}() where {T<:Real} = SteadyKinematics{T}(
    zeros(MVector{3,T}), ntuple(_ -> zero(T), 6)...,
    zeros(MVector{3,T}), ntuple(_ -> zero(T), 2)...,
    zero(T),
    ntuple(_ -> T(NaN), 2)...,
    zero(T),
)
SteadyKinematicsNaN() = SteadyKinematics{Float64}(
    fill(NaN, MVector{3,Float64}), ntuple(_ -> NaN, 6)...,
    fill(NaN, MVector{3,Float64}), ntuple(_ -> NaN, 2)...,
    NaN,
    ntuple(_ -> NaN, 2)...,
    NaN,
)
Base.copy(x::SteadyKinematics{T}) where {T} = copy_fields(x)


# UnsteadyKinematics composite type
# Fields
# - `Udot`: wind acceleration vector at spar position, resolved in basis W
# - `Uᵢdot`: in-plane component of Udot
# - `UₙdotMid`: normal wind acceleration component at the airfoil's 1/2-chord
# - `UₙdotTQC`: = normal wind acceleration component at the airfoil's 3/4-chord
# - `Ωₐdot`: spanwise angular acceleration component
# - `αₑ`: effective (unsteady) angle of attack at 3/4-chord
# - `αdot`: quasi-steady angle of attack rate
#
mutable struct UnsteadyKinematics{T<:Real}
    
    Udot::MVector{3,T}
    Uᵢdot::T
    UₙdotMid::T
    UₙdotTQC::T

    Ωₐdot::T

    αₑ::T
    αdot::T

end
UnsteadyKinematics{T}() where {T<:Real} = UnsteadyKinematics{T}(
    zeros(MVector{3,T}), ntuple(_ -> zero(T), 3)...,
    zero(T),
    ntuple(_ -> T(NaN), 2)...,
)
UnsteadyKinematicsNaN() = UnsteadyKinematics{Float64}(
    fill(NaN, MVector{3,Float64}), ntuple(_ -> NaN, 3)...,
    NaN,
    ntuple(_ -> NaN, 2)...,
)
Base.copy(x::UnsteadyKinematics{T}) where {T} = copy_fields(x)


# AeroCoefficients composite type
# Fields
# - `cn`: normal force coefficient
# - `cm`: pitching moment coefficient about spar position
# - `ct`: tangential force coefficient
# - `cnC`: circulatory component of cn
# - `cmC`: circulatory component of cm
# - `ctC`: circulatory component of ct
# - `cnNC`: non-circulatory component of cn
# - `cmNC`: non-circulatory component of cm
# - `ctNC`: non-circulatory component of ct
# - `cnI`: inertial component of cn
# - `cmI`: inertial component of cm
# - `cnP`: potential flow component of cn
# - `cnF`: separated-flow circulatory component of cn
# - `cmF`: separated-flow circulatory component of cm
# - `ctF`: separated-flow circulatory component of ct
# - `cnV`: DSV-induced component of cn
# - `cmV`: DSV-induced component of cm
# - `ctV`: DSV-induced component of ct
#
mutable struct AeroCoefficients{T<:Real}
    
    cn::T
    cm::T
    ct::T

    cnC::T
    cmC::T
    ctC::T

    cnNC::T
    cmNC::T
    ctNC::T

    cnI::T
    cmI::T

    cnP::T

    cnF::T
    cmF::T
    ctF::T

    cnV::T
    cmV::T
    ctV::T

end
AeroCoefficients{T}() where {T<:Real} = AeroCoefficients{T}(
    ntuple(_ -> zero(T), 3)...,
    ntuple(_ -> zero(T), 3)...,
    ntuple(_ -> zero(T), 3)...,
    ntuple(_ -> zero(T), 2)...,
    ntuple(_ -> zero(T), 1)...,
    ntuple(_ -> zero(T), 3)...,
    ntuple(_ -> zero(T), 3)...,
)
AeroCoefficientsNaN() = AeroCoefficients{Float64}(
    ntuple(_ -> NaN, 3)...,
    ntuple(_ -> NaN, 3)...,
    ntuple(_ -> NaN, 3)...,
    ntuple(_ -> NaN, 2)...,
    ntuple(_ -> NaN, 1)...,
    ntuple(_ -> NaN, 3)...,
    ntuple(_ -> NaN, 3)...,
)
Base.copy(x::AeroCoefficients{T}) where {T} = copy_fields(x)


# AttachedAeroVariables composite type
# Fields
# - `flowParameters`: flow parameters
# - `steadyKin`: steady kinematics
# - `unsteadyKin`: unsteady kinematics
# - `aeroCoefficients`: aerodynamic coefficients
#
struct AttachedAeroVariables{T<:Real}

    flowParameters::FlowParameters{T}
    steadyKin::SteadyKinematics{T}
    unsteadyKin::UnsteadyKinematics{T}
    aeroCoefficients::AeroCoefficients{T}

end
AttachedAeroVariables{T}() where {T<:Real} = AttachedAeroVariables{T}(
    FlowParameters{T}(), 
    SteadyKinematics{T}(), 
    UnsteadyKinematics{T}(),
    AeroCoefficients{T}(),
)
AttachedNaNAeroVariables() = AttachedAeroVariables{Float64}(
    FlowParametersNaN(),
    SteadyKinematicsNaN(),
    UnsteadyKinematicsNaN(),
    AeroCoefficientsNaN(),
)
Base.copy(x::AttachedAeroVariables{T}) where {T} = copy_fields(x)


# BLiNamedStates composite type
# Fields
# - `αlag`: lagged angle of attack
# - `f2primeN`: lagged separation point for cn
# - `f2primeM`: lagged separation point for cm
# - `f2primeT`: lagged separation point for ct
# - `RD`: lagged capped pitch rate ratio
# - `RD_stallOnsetRatio`: slightly lagged capped pitch rate ratio for stall onset
#
mutable struct BLiNamedStates{T<:Real}
    
    αlag::T
    f2primeN::T
    f2primeM::T
    f2primeT::T
    RD::T
    RD_stallOnsetRatio::T

end
BLiNamedStates{T}() where {T<:Real} = BLiNamedStates{T}(
    zero(T), one(T), one(T), one(T), zero(T), zero(T)
)
Base.copy(x::BLiNamedStates{T}) where {T} = copy_fields(x)


# BLiFlowVariables composite type
# Fields
# - `q`: non-dimensional pitch rate
# - `r`: reduced non-dimensional pitch rate
# - `qR`: unsigned ratio of reduced pitch rate to critical pitch rate
# - `R`: unsigned capped reduced pitch rate ratio 
# - `stallOnsetRatio`: ratio of current lagged angle to stall onset angle
# - `P`: lightness of stall factor (0-1, closer to 1 means lighter stall)
# - `SS`: steady condition indicator (0-1, closer to 1 means more steady)
# - `upstroke`: flag for motion on upstroke
# - `stalled`: flag for stalled condition on current cycle
# - `α1N`: unsteady breakpoint angle of attack for cn
# - `α1M`: unsteady breakpoint angle of attack for cm
# - `α1T`:unsteady breakpoint angle of attack for ct
# - `fN`: steady-state separation point for cn
# - `fM`: steady-state separation point for cm
# - `fT`: steady-state separation point for ct
# - `fPrimeN`: lagged separation point for cn
# - `fPrimeM`: lagged separation point for cm
# - `fPrimeT`: lagged separation point for ct
# - `Ta_SO`: time constant for stall onset
# - `TfN`: flow separation time constant for cn
# - `TfM`: flow separation time constant for cm
# - `TfT`: flow separation time constant for ct
#
mutable struct BLiFlowVariables{T<:Real}
    
    q::T
    r::T
    qR::T
    R::T
    stallOnsetRatio::T
    P::T
    SS::T

    upstroke::Bool
    stalled::Bool

    α1N::T
    α1M::T
    α1T::T

    fN::T
    fM::T
    fT::T

    fPrimeN::T
    fPrimeM::T
    fPrimeT::T

    Ta_SO::T
    TfN::T
    TfM::T
    TfT::T

end
BLiFlowVariables{T}() where {T<:Real} = BLiFlowVariables{T}(
    ntuple(_ -> zero(T), 7)...,
    ntuple(_ -> false, 2)...,
    ntuple(_ -> zero(T), 3)..., 
    ntuple(_ -> one(T), 3)...,
    ntuple(_ -> one(T), 3)...,
    ntuple(_ -> one(T), 4)...,
)
Base.copy(x::BLiFlowVariables{T}) where {T} = copy_fields(x)


# BLiComplementaryVariables composite type
# Fields
# - `stallOnsetRatioMax`: maximum value of stallOnsetRatio on current cycle
# - `stallOnsetRatioMin`: minimum value of stallOnsetRatio on current cycle
# - `stallOnsetRatioPrev`: previous value of stallOnsetRatio
# - `αlagPrev`: previous value of αlag
# - `qRPrev`: previous value of qR
# - `PPrev`: previous value of P
# - `upstrokePrev`: previous value of upstroke flag
# - `qRmax`: maximum value of qR on current cycle
# - `Ts`: Time decay function for tangential force
# - `tv0P`: time of last positive vortex shedding
# - `fDiff_tv0P`: f2primeN-fN @tv0P
# - `qR_tv0P`: qR @tv0P
# - `R_tv0P`: R @tv0P
# - `RD_tv0P`: RD @tv0P
# - `upstroke_tv0P`: upstroke flag @tv0P
# - `τvP`: nondimensional time since primary positive vortex shedding
# - `fDiff_tv0P2`: f2primeN-fN @tv0P2
# - `RD_tv0P2`: RD @tv0P2
# - `upstroke_tv0P2`: upstroke flag @tv0P2
# - `tv0N`: time of last negative vortex shedding
# - `fDiff_tv0N`: f2primeN-fN @tv0N
# - `qR_tv0N`: qR @tv0N
# - `R_tv0N`: R @tv0N
# - `RD_tv0N`: RD @tv0N
# - `upstroke_tv0N`: upstroke flag @tv0N
# - `τvN`: nondimensional time since primary negative vortex shedding
# - `fDiff_tv0N2`: f2primeN-fN @tv0N2
# - `RD_tv0N2`: RD @tv0N2
# - `upstroke_tv0N2`: upstroke flag @tv0N2
# - `lastRD_tv0`: last RD value @tv0P or @tv0N
#
mutable struct BLiComplementaryVariables{T<:Real}
    
    stallOnsetRatioMax::T
    stallOnsetRatioMin::T

    stallOnsetRatioPrev::T
    αlagPrev::T
    qRPrev::T
    PPrev::T
    upstrokePrev::Bool

    qRmax::T
    Ts::T

    tv0P::T
    fDiff_tv0P::T
    qR_tv0P::T
    R_tv0P::T
    RD_tv0P::T
    upstroke_tv0P::Bool
    secondPVortexYetToCome::Bool
    τvP::T
    fDiff_tv0P2::T
    RD_tv0P2::T 
    upstroke_tv0P2::Bool

    tv0N::T
    fDiff_tv0N::T
    qR_tv0N::T
    R_tv0N::T
    RD_tv0N::T 
    upstroke_tv0N::Bool
    secondNVortexYetToCome::Bool
    τvN::T
    fDiff_tv0N2::T
    RD_tv0N2::T  
    upstroke_tv0N2::Bool

    lastRD_tv0::T

end
BLiComplementaryVariables{T}() where {T<:Real} = BLiComplementaryVariables{T}(
    zero(T), one(T), 
    zero(T), zero(T), one(T), zero(T), false, 
    one(T), zero(T),
    T(-Inf), zero(T), one(T), one(T), one(T), false, false, T(Inf), -one(T), one(T), false,
    T(-Inf), zero(T), one(T), one(T), one(T), false, false, T(Inf), -one(T), one(T), false,
    one(T)
)
Base.copy(x::BLiComplementaryVariables{T}) where {T} = copy_fields(x)


# BLiAeroVariables composite type
# Fields
# - `states`: nonlinear states of BLi model
# - `flow`: flow variables of BLi model
# - `comp`: complementary variables of BLi model
#
struct BLiAeroVariables{T<:Real}

    states::BLiNamedStates{T}
    flow::BLiFlowVariables{T}
    comp::BLiComplementaryVariables{T}

end
BLiAeroVariables{T}() where {T<:Real} = BLiAeroVariables{T}(
    BLiNamedStates{T}(),
    BLiFlowVariables{T}(), 
    BLiComplementaryVariables{T}()
)
Base.copy(x::BLiAeroVariables{T}) where {T} = copy_fields(x)


# BLoNamedStates composite type
# Fields
# - `cnPprime`: lagged cnP
# - `f2Prime`: lagged separation point
# - `fPrimeM`: lagged separation point for cm on downstroke
# - `cnVP`: DSV-induced cn at positive angle of attack
# - `cnVN`: DSV-induced cn at negative angle of attack
#
mutable struct BLoNamedStates{T<:Real}

    cnPprime::T
    f2Prime::T
    fPrimeM::T
    cnVP::T
    cnVN::T

end
BLoNamedStates{T}() where {T<:Real} = BLoNamedStates{T}(
    zero(T), one(T), one(T), zero(T), zero(T)
)
Base.copy(x::BLoNamedStates{T}) where {T} = copy_fields(x)


# BLoFlowVariables composite type
# Fields
# - `αlag`: lagged angle of attack
# - `q`: non-dimensional pitch rate
# - `stallOnsetRatio`: ratio of current lagged angle to stall onset angle
# - `α1`: unsteady breakpoint angle of attack
# - `upstroke`: flag for motion on upstroke
# - `f`: steady-state separation point
# - `fPrime`: lagged separation point
# - `Tf`: flow separation time constant
# - `Tv`: vortex-shedding time constant
# - `Kf`: Kirchhoff-Helmholtz factor
# - `cvdotP`: positive vorticity coefficient rate
# - `cvdotN`: negative vorticity coefficient rate
#
mutable struct BLoFlowVariables{T<:Real}
    
    αlag::T
    q::T
    stallOnsetRatio::T
    α1::T
    
    upstroke::Bool

    f::T
    fPrime::T

    Tf::T
    Tv::T

    Kf::T
    cvdotP::T
    cvdotN::T

end
BLoFlowVariables{T}() where {T<:Real} = BLoFlowVariables{T}(
    ntuple(_ -> zero(T), 4)..., 
    false, 
    ntuple(_ -> one(T), 2)...,
    ntuple(_ -> one(T), 2)...,
    one(T), ntuple(_ -> zero(T), 2)...,
)
Base.copy(x::BLoFlowVariables{T}) where {T} = copy_fields(x)


# BLoComplementaryVariables composite type
# Fields
# - `stallOnsetRatioPrev`: previous value of stallOnsetRatio
# - `tv0P`: time of last positive vortex shedding
# - `tv0N`: time of last negative vortex shedding
# - `τvP`: nondimensional time since positive vortex shedding
# - `τvN`: nondimensional time since negative vortex shedding
#
mutable struct BLoComplementaryVariables{T<:Real}
    
    stallOnsetRatioPrev::T
    tv0P::T
    tv0N::T
    τvP::T
    τvN::T

end
BLoComplementaryVariables{T}() where {T<:Real} = BLoComplementaryVariables{T}(
    zero(T), T(-Inf), T(-Inf), T(Inf), T(Inf)
)
Base.copy(x::BLoComplementaryVariables{T}) where {T} = copy_fields(x)


# BLoAeroVariables composite type
# Fields
# - `states`: nonlinear states of BLo model
# - `flow`: flow variables of BLo model
# - `comp`: complementary variables of BLo model
#
struct BLoAeroVariables{T<:Real}

    states::BLoNamedStates{T}
    flow::BLoFlowVariables{T}
    comp::BLoComplementaryVariables{T}

end
BLoAeroVariables{T}() where {T<:Real} = BLoAeroVariables{T}(
    BLoNamedStates{T}(),
    BLoFlowVariables{T}(), 
    BLoComplementaryVariables{T}()
)
Base.copy(x::BLoAeroVariables{T}) where {T} = copy_fields(x)


# AeroLoadsResultants composite type
# Fields
# - `f1χ`: first node aerodynamic forces resultant, resolved in basis A
# - `f2χ`: second node aerodynamic forces resultant, resolved in basis A
# - `m1χ`: first node aerodynamic moments resultant, resolved in basis A
# - `m2χ`: second node aerodynamic moments resultant, resolved in basis A
# - `Fa`: aerodynamic forces and moments resultant vector, resolved in basis W
#
struct AeroLoadsResultants{T<:Real}
    
    f1χ::MVector{3,T}
    f2χ::MVector{3,T}
    m1χ::MVector{3,T}
    m2χ::MVector{3,T}

    Fa::MVector{12,T}

end
AeroLoadsResultants{T}() where {T<:Real} = AeroLoadsResultants{T}(
    ntuple(_ -> zeros(MVector{3,T}), 4)...,
    zeros(MVector{12,T}),
)


# AeroStateMatrices composite type
# Fields
# - `Ac`: aerodynamic state vector for linear circulatory states
# - `Ai`: aerodynamic state vector for inertial states
# - `An`: aerodynamic state vector for nonlinear states
# - `Af`: aerodynamic state vector for flap states
# - `Ag`: aerodynamic state vector for gust states
# - `A`: aerodynamic state matrix 
# - `B`: aerodynamic input matrix
#
struct AeroStateMatrices{T<:Real,N,Nc,Ni,Nn,Nf,Ng,NN}

    Ac::MVector{Nc,T}
    Ai::MVector{Ni,T}
    An::MVector{Nn,T}
    Af::MVector{Nf,T}
    Ag::MVector{Ng,T}

    A::MMatrix{N,N,T,NN}
    B::MVector{N,T}

    function AeroStateMatrices{T,N,Nc,Ni,Nn,Nf,Ng,NN}(Ac, Ai, An, Af, Ag, A, B) where {T,N,Nc,Ni,Nn,Nf,Ng,NN}
        new{T,N,Nc,Ni,Nn,Nf,Ng,NN}(Ac, Ai, An, Af, Ag, A, B)
    end

end
AeroStateMatrices{T,N,Nc,Ni,Nn,Nf,Ng}() where {T<:Real,N,Nc,Ni,Nn,Nf,Ng} =
    AeroStateMatrices{T,N,Nc,Ni,Nn,Nf,Ng,N*N}(
        zeros(MVector{Nc,T}), zeros(MVector{Ni,T}), zeros(MVector{Nn,T}), zeros(MVector{Nf,T}), zeros(MVector{Ng,T}),
        zeros(MMatrix{N,N,T,N*N}), zeros(MVector{N,T}),
    )


# AeroStateCache composite type
# Fields
# - `attachedVars`: attached-flow aerodynamic variables
# - `bliVars`: BLi model aerodynamic variables
# - `bloVars`: BLo model aerodynamic variables
# - `loadsResultants`: aerodynamic loads resultants
# - `stateMatrices`: aerodynamic state matrices
# - `statesBuffer`: buffer array for the input states
# - `statesRatesBuffer`: buffer array for the input states rates
# - `outBuffer`: buffer array for vcat(f1χ,f2χ,m1χ,m2χ,Ac,Ai,An,Af,Ag,B)
#
mutable struct AeroStateCache{T<:Real,N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,L,P}
    
    attachedVars::AttachedAeroVariables{T}
    bliVars::BLiAeroVariables{T}
    bloVars::BLoAeroVariables{T}
    loadsResultants::AeroLoadsResultants{T}
    stateMatrices::AeroStateMatrices{T,N,Nc,Ni,Nn,Nf,Ng,NN}

    statesBuffer::MVector{L,T}
    statesRatesBuffer::MVector{6,T}
    outBuffer::MVector{P,T}

end
AeroStateCache{T,N,Nc,Ni,Nn,Nf,Ng,Nδ}() where {T,N,Nc,Ni,Nn,Nf,Ng,Nδ} = 
    AeroStateCache{T,N,Nc,Ni,Nn,Nf,Ng,Nδ,N*N,9+N+Nδ,12+Nc+Ni+Nn+Nf+Ng+N}(
        AttachedAeroVariables{T}(),
        BLiAeroVariables{T}(),
        BLoAeroVariables{T}(),
        AeroLoadsResultants{T}(),
        AeroStateMatrices{T,N,Nc,Ni,Nn,Nf,Ng}(),
        zeros(MVector{9+N+Nδ,T}),
        zeros(MVector{6,T}),
        zeros(MVector{12+Nc+Ni+Nn+Nf+Ng+N,T}),
    )


# AeroMatricesDerivatives composite type
# Fields
# - `Ac_p`: ∂Ac/∂ṗ
# - `Ai_p`: ∂Ai/∂ṗ
# - `An_p`: ∂An/∂ṗ
# - `Af_p`: ∂Af/∂ṗ
# - `Ag_p`: ∂Ag/∂ṗ
# - `Ac_V`: ∂Ac/∂V
# - `Ai_V`: ∂Ai/∂V
# - `An_V`: ∂An/∂V
# - `Af_V`: ∂Af/∂V
# - `Ag_V`: ∂Ag/∂V
# - `Ac_Ω`: ∂Ac/∂Ω
# - `Ai_Ω`: ∂Ai/∂Ω
# - `An_Ω`: ∂An/∂Ω
# - `Af_Ω`: ∂Af/∂Ω
# - `Ag_Ω`: ∂Ag/∂Ω
# - `Ac_χ`: ∂Ac/∂χ
# - `Ai_χ`: ∂Ai/∂χ
# - `An_χ`: ∂An/∂χ
# - `Af_χ`: ∂Af/∂χ
# - `Ag_χ`: ∂Ag/∂χ
# - `A_p`: ∂A/∂ṗ
# - `A_V`: ∂A/∂V
# - `A_Ω`: ∂A/∂Ω
# - `A_χ`: ∂A/∂χ
# - `B_p`: ∂B/∂ṗ
# - `B_V`: ∂B/∂V
# - `B_Ω`: ∂B/∂Ω
# - `B_χ`: ∂B/∂χ
# - `Ac_Vdot`: ∂Ac/∂Vdot
# - `Ai_Vdot`: ∂Ai/∂Vdot
# - `An_Vdot`: ∂An/∂Vdot
# - `Af_Vdot`: ∂Af/∂Vdot
# - `Ag_Vdot`: ∂Ag/∂Vdot
# - `Ac_Ωdot`: ∂Ac/∂Ωdot
# - `Ai_Ωdot`: ∂Ai/∂Ωdot
# - `An_Ωdot`: ∂An/∂Ωdot
# - `Af_Ωdot`: ∂Af/∂Ωdot
# - `Ag_Ωdot`: ∂Ag/∂Ωdot
# - `A_Vdot`: ∂A/∂Vdot
# - `A_Ωdot`: ∂A/∂Ωdot
# - `B_Vdot`: ∂B/∂Vdot
# - `B_Ωdot`: ∂B/∂Ωdot
#
struct AeroMatricesDerivatives{N,Nc,Ni,Nn,Nf,Ng,N3,Nc3,Ni3,Nn3,Nf3,Ng3,NcN,NiN,NnN,NfN,NgN,NN,NN3,NNN}

    Ac_p::MMatrix{Nc,3,Float64,Nc3}
    Ai_p::MMatrix{Ni,3,Float64,Ni3}
    An_p::MMatrix{Nn,3,Float64,Nn3}
    Af_p::MMatrix{Nf,3,Float64,Nf3}
    Ag_p::MMatrix{Ng,3,Float64,Ng3}

    Ac_V::MMatrix{Nc,3,Float64,Nc3}
    Ai_V::MMatrix{Ni,3,Float64,Ni3}
    An_V::MMatrix{Nn,3,Float64,Nn3}
    Af_V::MMatrix{Nf,3,Float64,Nf3}
    Ag_V::MMatrix{Ng,3,Float64,Ng3}
    
    Ac_Ω::MMatrix{Nc,3,Float64,Nc3}
    Ai_Ω::MMatrix{Ni,3,Float64,Ni3}
    An_Ω::MMatrix{Nn,3,Float64,Nn3}
    Af_Ω::MMatrix{Nf,3,Float64,Nf3}
    Ag_Ω::MMatrix{Ng,3,Float64,Ng3}

    Ac_χ::MMatrix{Nc,N,Float64,NcN}
    Ai_χ::MMatrix{Ni,N,Float64,NiN}
    An_χ::MMatrix{Nn,N,Float64,NnN}
    Af_χ::MMatrix{Nf,N,Float64,NfN}
    Ag_χ::MMatrix{Ng,N,Float64,NgN}

    A_p::MArray{Tuple{N,N,3},Float64,3,NN3}
    A_V::MArray{Tuple{N,N,3},Float64,3,NN3}
    A_Ω::MArray{Tuple{N,N,3},Float64,3,NN3}
    A_χ::MArray{Tuple{N,N,N},Float64,3,NNN}
    B_p::MMatrix{N,3,Float64,N3}
    B_V::MMatrix{N,3,Float64,N3}
    B_Ω::MMatrix{N,3,Float64,N3}
    B_χ::MMatrix{N,N,Float64,NN}

    Ac_Vdot::MMatrix{Nc,3,Float64,Nc3}
    Ai_Vdot::MMatrix{Ni,3,Float64,Ni3}
    An_Vdot::MMatrix{Nn,3,Float64,Nn3}
    Af_Vdot::MMatrix{Nf,3,Float64,Nf3}
    Ag_Vdot::MMatrix{Ng,3,Float64,Ng3}

    Ac_Ωdot::MMatrix{Nc,3,Float64,Nc3}
    Ai_Ωdot::MMatrix{Ni,3,Float64,Ni3}
    An_Ωdot::MMatrix{Nn,3,Float64,Nn3}
    Af_Ωdot::MMatrix{Nf,3,Float64,Nf3}
    Ag_Ωdot::MMatrix{Ng,3,Float64,Ng3}
    
    A_Vdot::MArray{Tuple{N,N,3},Float64,3,NN3}
    A_Ωdot::MArray{Tuple{N,N,3},Float64,3,NN3}
    B_Vdot::MMatrix{N,3,Float64,N3}
    B_Ωdot::MMatrix{N,3,Float64,N3}

    function AeroMatricesDerivatives{N,Nc,Ni,Nn,Nf,Ng,N3,Nc3,Ni3,Nn3,Nf3,Ng3,NcN,NiN,NnN,NfN,NgN,NN,NN3,NNN}(Ac_p,Ai_p,An_p,Af_p,Ag_p,Ac_V,Ai_V,An_V,Af_V,Ag_V,Ac_Ω,Ai_Ω,An_Ω,Af_Ω,Ag_Ω,Ac_χ,Ai_χ,An_χ,Af_χ,Ag_χ,A_p,A_V,A_Ω,A_χ,B_p,B_V,B_Ω,B_χ,Ac_Vdot,Ai_Vdot,An_Vdot,Af_Vdot,Ag_Vdot,Ac_Ωdot,Ai_Ωdot,An_Ωdot,Af_Ωdot,Ag_Ωdot,A_Vdot,A_Ωdot,B_Vdot,B_Ωdot) where {N,Nc,Ni,Nn,Nf,Ng,N3,Nc3,Ni3,Nn3,Nf3,Ng3,NcN,NiN,NnN,NfN,NgN,NN,NN3,NNN}
        new{N,Nc,Ni,Nn,Nf,Ng,N3,Nc3,Ni3,Nn3,Nf3,Ng3,NcN,NiN,NnN,NfN,NgN,NN,NN3,NNN}(Ac_p,Ai_p,An_p,Af_p,Ag_p,Ac_V,Ai_V,An_V,Af_V,Ag_V,Ac_Ω,Ai_Ω,An_Ω,Af_Ω,Ag_Ω,Ac_χ,Ai_χ,An_χ,Af_χ,Ag_χ,A_p,A_V,A_Ω,A_χ,B_p,B_V,B_Ω,B_χ,Ac_Vdot,Ai_Vdot,An_Vdot,Af_Vdot,Ag_Vdot,Ac_Ωdot,Ai_Ωdot,An_Ωdot,Af_Ωdot,Ag_Ωdot,A_Vdot,A_Ωdot,B_Vdot,B_Ωdot)
    end

end
AeroMatricesDerivatives{N,Nc,Ni,Nn,Nf,Ng}() where {N,Nc,Ni,Nn,Nf,Ng} =
    AeroMatricesDerivatives{N,Nc,Ni,Nn,Nf,Ng,3*N,3*Nc,3*Ni,3*Nn,3*Nf,3*Ng,Nc*N,Ni*N,Nn*N,Nf*N,Ng*N,N*N,3*N*N,N*N*N}(
        zeros(MMatrix{Nc,3,Float64}), zeros(MMatrix{Ni,3,Float64}), zeros(MMatrix{Nn,3,Float64}), zeros(MMatrix{Nf,3,Float64}), zeros(MMatrix{Ng,3,Float64}),
        zeros(MMatrix{Nc,3,Float64}), zeros(MMatrix{Ni,3,Float64}), zeros(MMatrix{Nn,3,Float64}), zeros(MMatrix{Nf,3,Float64}), zeros(MMatrix{Ng,3,Float64}),
        zeros(MMatrix{Nc,3,Float64}), zeros(MMatrix{Ni,3,Float64}), zeros(MMatrix{Nn,3,Float64}), zeros(MMatrix{Nf,3,Float64}), zeros(MMatrix{Ng,3,Float64}),
        zeros(MMatrix{Nc,N,Float64}), zeros(MMatrix{Ni,N,Float64}), zeros(MMatrix{Nn,N,Float64}), zeros(MMatrix{Nf,N,Float64}), zeros(MMatrix{Ng,N,Float64}),
        zeros(MArray{Tuple{N,N,3},Float64,3,3*N*N}), zeros(MArray{Tuple{N,N,3},Float64,3,3*N*N}), zeros(MArray{Tuple{N,N,3},Float64,3,3*N*N}), zeros(MArray{Tuple{N,N,N},Float64,3,N*N*N}), zeros(MMatrix{N,3,Float64}), zeros(MMatrix{N,3,Float64}), zeros(MMatrix{N,3,Float64}), zeros(MMatrix{N,N,Float64}),
        zeros(MMatrix{Nc,3,Float64}), zeros(MMatrix{Ni,3,Float64}), zeros(MMatrix{Nn,3,Float64}), zeros(MMatrix{Nf,3,Float64}), zeros(MMatrix{Ng,3,Float64}),
        zeros(MMatrix{Nc,3,Float64}), zeros(MMatrix{Ni,3,Float64}), zeros(MMatrix{Nn,3,Float64}), zeros(MMatrix{Nf,3,Float64}), zeros(MMatrix{Ng,3,Float64}),
        zeros(MArray{Tuple{N,N,3},Float64,3,3*N*N}), zeros(MArray{Tuple{N,N,3},Float64,3,3*N*N}), zeros(MMatrix{N,3,Float64}), zeros(MMatrix{N,3,Float64}),
    )


# AeroProperties composite type
# Fields
# - `surface`: aerodynamic surface of the element's parent beam
# - `geometry`: aerodynamic geometric properties
# - `stateRanges`: ranges of element-dependent aerodynamic states (flap and gust)
# - `loadsDistribution`: aerodynamic loads distribution properties
# - `stateCache`: aerodynamic state cache
# - `stateCacheDiff`: pre-allocated aerodynamic state cache for Jacobian calculations
# - `matricesDerivatives`: aerodynamic matrices' derivatives
#
@with_kw mutable struct AeroProperties{N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,N3,Nc3,Ni3,Nn3,Nf3,Ng3,NcN,NiN,NnN,NfN,NgN,NN3,NNN,L,P,S1,S2,S3,F}

    # Surface
    surface::Union{Nothing,AeroSurface{S1,S2,S3}}

    # Geometry
    geometry::Union{Nothing,AeroGeometry}

    # State ranges
    stateRanges::Union{Nothing,StateRanges}

    # Loads distribution
    loadsDistribution::Union{Nothing,AeroLoadsDistribution{F}}

    # State cache
    stateCache::AeroStateCache{Float64,N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,L,P}
    
    # State cache for Jacobian calculations
    stateCacheDiff::Union{Nothing,AeroStateCache}

    # Aerodynamic matrices' derivatives
    matricesDerivatives::Union{Nothing,AeroMatricesDerivatives{N,Nc,Ni,Nn,Nf,Ng,N3,Nc3,Ni3,Nn3,Nf3,Ng3,NcN,NiN,NnN,NfN,NgN,NN,NN3,NNN}}

end


## PRIVATE API
# ------------------------------------------------------------------------------
# AeroProperties constructor
function AeroProperties(surface,R0,x1,x1_norm,x1_n1_norm,x1_n2_norm)

    # Aerodynamic geometric properties
    geometry = isnothing(surface) ?
               nothing :
               get_aero_geometric_properties(surface,R0,x1,x1_norm)

    # State ranges
    stateRanges = isnothing(surface) ?
                  nothing :
                  get_state_ranges(surface,geometry,false)              

    # Aerodynamic loads distribution properties
    loadsDistribution = isnothing(surface) ?
                        nothing :
                        get_aero_loads_distribution(surface,x1_n1_norm,x1_n2_norm)

    # Types and sizes
    if isnothing(surface)
        N = Nc = Ni = Nn = Nf = Ng = Nδ = 0
    else
        Nc = surface.solvers.ppSolver isa Inflow ? 
             surface.solvers.ppSolver.nInflowStates^2 : 
             surface.solvers.ppSolver.nCirculatoryStates
        Ni = surface.solvers.ppSolver.nInertialStates
        Nn = hasfield(typeof(surface.solvers.ppSolver), :nNonlinearStates) ? 
             surface.solvers.ppSolver.nNonlinearStates : 
             0
        Nf = geometry.isFlapped ? 
             surface.stateRanges.nFlapStates : 
             0
        Ng = surface.stateRanges.nGustStates
        Nδ = geometry.isFlapped && surface.assemblyData.δIsTrimVariable ? 
             1 :
             0
        N  =  surface.solvers.ppSolver.nStates + Nf + Ng
    end
    S1 = isnothing(surface) ? 
         NoAeroSolver : 
         typeof(surface.solvers.ppSolver)
    S2 = isnothing(surface) ? 
         NoFlapAeroSolver : 
         typeof(surface.solvers.flapSolver)
    S3 = isnothing(surface) ? 
         NoGustAeroSolver : 
         typeof(surface.solvers.gustSolver)
    F = isnothing(surface) ? 
        typeof(x -> 0) : 
        typeof(loadsDistribution.ϖ)      

    # Initalize state caches
    stateCache = AeroStateCache{Float64,N,Nc,Ni,Nn,Nf,Ng,Nδ}()
    stateCacheDiff = nothing

    # Initialize matrices derivatives
    matricesDerivatives = isnothing(surface) ?
                          nothing :
                          AeroMatricesDerivatives{N,Nc,Ni,Nn,Nf,Ng}()

    return AeroProperties{N,Nc,Ni,Nn,Nf,Ng,Nδ,N*N,3*N,3*Nc,3*Ni,3*Nn,3*Nf,3*Ng,Nc*N,Ni*N,Nn*N,Nf*N,Ng*N,3*N*N,N*N*N,9+N+Nδ,12+Nc+Ni+Nn+Nf+Ng+N,S1,S2,S3,F}(
        surface=surface,
        geometry=geometry,
        stateRanges=stateRanges,
        loadsDistribution=loadsDistribution,
        stateCache=stateCache,
        stateCacheDiff=stateCacheDiff,
        matricesDerivatives=matricesDerivatives,
        ), N, Nc, Ni, Nn, Nf, Ng, Nδ, S1, S2, S3, F

end


# AeroProperties re-constructor
function AeroProperties(element,N,Nδ,withGust)

    @unpack aero = element
    @unpack surface = aero

    # Get states sizes
    Nc = length(aero.stateCache.stateMatrices.Ac)
    Ni = length(aero.stateCache.stateMatrices.Ai)
    Nn = length(aero.stateCache.stateMatrices.An)
    Nf = length(aero.stateCache.stateMatrices.Af)
    Ng = withGust ? length(aero.stateCache.stateMatrices.Ag) : 0

    # Types
    S1 = typeof(aero.surface.solvers.ppSolver)
    S2 = typeof(aero.surface.solvers.flapSolver)
    S3 = typeof(aero.surface.solvers.gustSolver)
    F = typeof(aero.loadsDistribution.ϖ)

    # State ranges
    stateRanges = get_state_ranges(surface,aero.geometry,withGust)

    # Initalize state caches
    stateCache = AeroStateCache{Float64,N,Nc,Ni,Nn,Nf,Ng,Nδ}()
    stateCacheDiff = nothing

    # Initialize matrices derivatives
    matricesDerivatives = isnothing(surface) ?
                          nothing :
                          AeroMatricesDerivatives{N,Nc,Ni,Nn,Nf,Ng}()

    # Reconstruct AeroProperties with no gust states
    return AeroProperties{N,Nc,Ni,Nn,Nf,Ng,Nδ,N*N,3*N,3*Nc,3*Ni,3*Nn,3*Nf,3*Ng,Nc*N,Ni*N,Nn*N,Nf*N,Ng*N,3*N*N,N*N*N,9+N+Nδ,12+Nc+Ni+Nn+Nf+Ng+N,S1,S2,S3,F}(
        surface=aero.surface,
        geometry=aero.geometry,
        loadsDistribution=aero.loadsDistribution,
        stateRanges=stateRanges,
        stateCache=stateCache,
        stateCacheDiff=stateCacheDiff,
        matricesDerivatives=matricesDerivatives,
    ), Nc, Ni, Nn, Nf, Ng, Nδ, S1, S2, S3, F

end


# Computes the aerodynamic geometric properties of a surface at a given spanwise location
function get_aero_geometric_properties(surface,R0,x1,x1_norm)

    # Possibly spanwise-varying geometric properties 
    c = surface.geometry.c isa Float64 ? 
        surface.geometry.c : 
        surface.geometry.c(x1)
    b = c/2
    normSparPos = surface.geometry.normSparPos isa Float64 ? 
                  surface.geometry.normSparPos : 
                  surface.geometry.normSparPos(x1)
    aₕ = 2*normSparPos - 1
    Λ = surface.geometry.Λ isa Float64 ? 
        surface.geometry.Λ : 
        surface.geometry.Λ(x1)
    φ = surface.geometry.φ isa Float64 ? 
        surface.geometry.φ : 
        surface.geometry.φ(x1)
    cosΛ = cos(Λ)

    # Rotation tensor from basis b to basis w and its transpose
    Rw = SMatrix{3,3,Float64,9}(rotation_tensor_E321([Λ; 0; φ]))
    RwT = Rw'

    # Rotation tensor from basis A to basis w, resolved in basis A
    RwR0 = Rw * R0

    # Flag for being flapped
    isFlapped = !isnothing(surface.geometry.normFlapSpan) && minimum(surface.geometry.normFlapSpan) < x1_norm < maximum(surface.geometry.normFlapSpan)

    # Set Theodorsen's flap constants, if applicable
    Th = isFlapped && surface.solvers.flapSolver isa ThinAirfoilTheory ?
         theodorsen_flap_constants(normSparPos,surface.geometry.normFlapPos) : 
         nothing

    # Set geometry data 
    geometry = AeroGeometry(b,c,normSparPos,aₕ,Λ,φ,cosΛ,cosΛ^2,c/cosΛ,Rw,RwT,RwR0,isFlapped,Th)

    return geometry

end


# Computes Theodorsen's flap constants
function theodorsen_flap_constants(normSparPos,normFlapPos)

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
    Th[9] = 1/2*(1/3*(sqrt(1-d^2))^3+a*Th[4])
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

    return SVector{18,Float64}(Th)

end


# Computes the ranges of the element-dependent aerodynamic states (flap, gust)
function get_state_ranges(surface,geometry,withGust)

    # Unpack data
    @unpack nPitchPlungeStates,nFlapStates,nGustStates = surface.stateRanges
    @unpack isFlapped = geometry

    # Initialize offset for state indexing
    offset = nPitchPlungeStates

    # Flap states range
    if isFlapped
        flapStatesRange = (offset+1 : offset+nFlapStates)
        offset += nFlapStates
    else
        flapStatesRange = (1 : 0)
    end

    # Gust states ranges
    gustStatesRange = withGust ? 
                      (offset+1 : offset+nGustStates) : 
                      (1 : 0)
    
    # Set state ranges data
    stateRanges = StateRanges(flapStatesRange,gustStatesRange)

    return stateRanges

end


# Computes the aerodynamic loads distribution properites over the element's span
function get_aero_loads_distribution(surface,x1_n1_norm,x1_n2_norm)

    # Set local tip loss function
    if !surface.tipLoss.hasTipCorrection # No tip correction -> ϖ is constant and equal to one
        ϖ = (U,ζ) -> 1.0
    else # Tip correction is active -> ϖ is a function of the element's local coordinate (ζ)
        # Arclength coordinates
        s = ζ -> x1_n1_norm + ζ * (x1_n2_norm-x1_n1_norm)
        sFlip = ζ -> (1-x1_n1_norm) + ζ * ((1-x1_n2_norm)-(1-x1_n1_norm))
        if isnothing(surface.tipLoss.tipLossFunction) # No input
            # Set standard tip loss function with given tipLossDecayFactor
            ϖ = surface.tipLoss.tipLossDecayFactor >= 0 ? 
                (U,ζ) -> 1.0 - exp(-surface.tipLoss.tipLossDecayFactor*(1-s(ζ))) : 
                (U,ζ) -> 1.0 - exp(surface.tipLoss.tipLossDecayFactor*(1-sFlip(ζ)))
        else # Tip loss function is given
            ϖ = surface.tipLoss.tipLossFunctionIsAirspeedDependent ? 
                build_varpi_from_tlf_airspeed(surface.tipLoss.tipLossFunction, s) : 
                build_varpi_from_tlf(surface.tipLoss.tipLossFunction, s)
        end
    end

    # Initialize reference airspeed
    Uref = 0.0

    # Set loads distribution data
    loadsDistribution = AeroLoadsDistribution(ϖ,ϖ(Uref,1/2),Uref)

    return loadsDistribution

end


# Function barriers to concretize tipLossFunction type and avoid type instability in the loads distribution function ϖ
function build_varpi_from_tlf(tipLossFunction::T, s) where T
    return (U,ζ) -> tipLossFunction(s(ζ))
end
function build_varpi_from_tlf_airspeed(tipLossFunction::T, s) where T
    return (U,ζ) -> tipLossFunction(U)(s(ζ))
end


# Computes the initial value of χ (for zero relative airspeed)
function initial_χ(N,surface)

    # Initialize
    χ = zeros(MVector{N,Float64})

    # Return if the is no aerodynamic surface 
    if isnothing(surface)
        return χ
    end

    # Unpack surface data
    @unpack ppSolver = surface.solvers
    @unpack nonlinearStatesRange = surface.stateRanges

    # Set delayed separation states (f2Prime) to 1.0 for BL solvers
    if ppSolver isa BLi
        χ[nonlinearStatesRange[2:4]] .= 1.0
    elseif ppSolver isa BLo
        χ[nonlinearStatesRange[2:3]] .= 1.0
    end

    return χ

end


# Computes the initial value of F_χ_χ (for zero relative airspeed)
function initial_F_χ_χ(N,surface)

    # Return if the is no aerodynamic surface 
    if isnothing(surface)
        return zeros(MMatrix{N,N,Float64,N*N})
    end

    @unpack ppSolver = surface.solvers

    # Set according to ppSolver
    F_χ_χ = if ppSolver isa QuasiSteady
        zeros(MMatrix{N,N,Float64,N*N})
    elseif ppSolver isa IndicialSolver
        -sqrt(eps()) * MMatrix{N,N,Float64,N*N}(I)
    elseif ppSolver isa Inflow
        # -sqrt(eps()) * MMatrix{N,N,Float64,N*N}(ppSolver.Aₚ)
        -sqrt(eps()) * MMatrix{N,N,Float64,N*N}(I)
    end

    return F_χ_χ

end


# Computes the initial value of F_χ_Vdot (for zero relative airspeed)
function initial_F_χ_Vdot(N,surface)

    # Initialize
    F_χ_Vdot = zeros(MMatrix{N,3,Float64,3*N})

    # Return if the is no aerodynamic surface 
    if isnothing(surface)
        return F_χ_Vdot
    end

    @unpack airfoil = surface.geometry
    @unpack ppSolver = surface.solvers
    @unpack circulatoryPitchPlungeStatesRange = surface.stateRanges
    cnα = get_airfoil_state_field(airfoil, ppSolver, :cnα)

    # Calculate according to ppSolver
    if ppSolver isa IndicialSolver
        F_χ_Vdot[circulatoryPitchPlungeStatesRange[1:length(ppSolver.AC)],3] = cnα*ppSolver.AC
    elseif ppSolver isa Inflow
        F_χ_Vdot[circulatoryPitchPlungeStatesRange,3] = ppSolver.AₚInvcₚ
    end

    return F_χ_Vdot

end


# Computes the initial value of F_χ_Ωdot (for zero relative airspeed)
function initial_F_χ_Ωdot(N,surface,geometry)

    # Initialize
    F_χ_Ωdot = zeros(MMatrix{N,3,Float64,3*N})

    # Return if the is no aerodynamic surface 
    if isnothing(surface)
        return F_χ_Ωdot
    end

    @unpack c,normSparPos = geometry
    @unpack airfoil = surface.geometry
    @unpack ppSolver = surface.solvers
    @unpack circulatoryPitchPlungeStatesRange = surface.stateRanges
    cnα = get_airfoil_state_field(airfoil, ppSolver, :cnα)

    # Calculate according to ppSolver
    if ppSolver isa IndicialSolver
        F_χ_Ωdot[circulatoryPitchPlungeStatesRange[1:length(ppSolver.AC)],1] = c*(normSparPos-3/4)*cnα*ppSolver.AC
    elseif ppSolver isa Inflow
        F_χ_Ωdot[circulatoryPitchPlungeStatesRange,1] = c*(normSparPos-3/4)*ppSolver.AₚInvcₚ  
    end

    return F_χ_Ωdot

end


# Returns the state cache for the calculations of jacobians, creating a new AeroStateCache if necessary
function get_diff_state_cache!(aero::AeroProperties{N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,N3,Nc3,Ni3,Nn3,Nf3,Ng3,NcN,NiN,NnN,NfN,NgN,NN3,NNN,L,P,S1,S2,S3,F}, ::Type{T})::AeroStateCache{T,N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,L,P} where {T<:Real,N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,N3,Nc3,Ni3,Nn3,Nf3,Ng3,NcN,NiN,NnN,NfN,NgN,NN3,NNN,L,P,S1,S2,S3,F}
    cacheDiff = aero.stateCacheDiff
    ppSolver = aero.surface.solvers.ppSolver
    if cacheDiff isa AeroStateCache{T,N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,L,P}
        copy_aero_state!(cacheDiff, aero.stateCache, ppSolver)
        return cacheDiff
    else
        newCacheDiff = AeroStateCache{T,N,Nc,Ni,Nn,Nf,Ng,Nδ}()
        aero.stateCacheDiff = newCacheDiff
        return newCacheDiff
    end
end


# Copy function for AeroStateCache Float64 → Dual
@generated function copy_aero_state!(dst::AeroStateCache{T,N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,L,P}, src::AeroStateCache{S,N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,L,P}, ppSolver::AeroSolver) where {T<:Real,S<:Real,N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,L,P}
    exprs = Expr[]
    for f in fieldnames(dst)
        if f === :loadsResultants || f === :stateMatrices
            push!(exprs, :(copy_aero_state_struct_fields!(dst.$f, src.$f)))
        elseif f === :attachedVars
            push!(exprs, :(copy_aero_vars!(dst.$f, src.$f)))
        elseif f === :bliVars
            push!(exprs, :(ppSolver isa BLi && copy_aero_vars!(dst.$f, src.$f)))
        elseif f === :bloVars
            push!(exprs, :(ppSolver isa BLo && copy_aero_vars!(dst.$f, src.$f)))
        end
    end
    return quote
        @inbounds begin
            $(exprs...)
        end
        return dst
    end
end

# Copy function for aero variables structs
@generated function copy_aero_vars!(dst::Union{AttachedAeroVariables{T},BLiAeroVariables{T},BLoAeroVariables{T}}, src::Union{AttachedAeroVariables{S},BLiAeroVariables{S},BLoAeroVariables{S}}) where {T<:Real,S<:Real}
    fields = fieldnames(dst)
    exprs = map(fields) do f
        :(copy_aero_state_struct_fields!(getfield(dst, $(QuoteNode(f))), getfield(src, $(QuoteNode(f)))))
    end
    quote
        @inbounds begin
            $(exprs...)
        end
        return dst
    end
end

# Functional unit of copy_aero_state! for struct fields
@generated function copy_aero_state_struct_fields!(dst, src)
    fields = fieldnames(dst)
    exprs = map(fields) do f
        quote
            dst_f = getfield(dst, $(QuoteNode(f)))
            if dst_f isa AbstractArray
                copy_aero_state_array_fields!(dst_f, getfield(src, $(QuoteNode(f))))
            else
                setfield!(dst, $(QuoteNode(f)), convert(fieldtype(typeof(dst), $(QuoteNode(f))), getfield(src, $(QuoteNode(f)))))
            end
        end
    end
    quote
        @inbounds begin
            $(exprs...)
        end
        return dst
    end
end

# Functional unit of copy_aero_state! for array fields
copy_aero_state_array_fields!(dst::AbstractArray{T}, src::AbstractArray) where {T} = (dst .= T.(src); dst)