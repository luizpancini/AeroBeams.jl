#

    # AeroSurface composite type

#
@with_kw mutable struct AeroSurface

    # Primary (inputs to aero surface creation)
    # -----------------------------------------
    # Aerodynamic solver
    solver::AeroSolver
    # Flap loads solver
    flapLoadsSolver::FlapAeroSolver
    # Gust loads solver
    gustLoadsSolver::GustAeroSolver
    # Derivation method 
    derivationMethod::DerivationMethod
    # Airfoil
    airfoil::Airfoil
    # Airfoil chord (constant or a function of normalized span)
    c::Union{<:Function,Number}
    # Sweep angle (constant or a function of normalized span)
    Λ::Union{<:Function,Number} 
    # Spar position normalized by local chord (constant or a function of normalized span)
    normSparPos::Union{<:Function,Float64}
    # Trailing-edge flap span normalized by attached beam length
    normFlapSpan::Union{Nothing,Vector{<:Number}}
    # Trailing-edge flap chord normalized by local airfoil chord (only used if flapLoadsSolver is of type ThinAirfoilTheory)
    normFlapPos::Union{Nothing,Float64}
    # TF for flap deflection being a user input
    δIsInput::Bool
    # TF for flap deflection being a trim variable
    δIsTrimVariable::Bool
    # TF for flap deflection being zero over time
    δIsZero::Bool 
    # Flap deflection 
    δ::Union{Nothing,<:Function,Number}
    # Flap site ID (for flapLoadsSolver of type TableLookup)
    flapSiteID::Union{Nothing,Int64}
    # TF to update airfoil parameters according to flow parameters
    updateAirfoilParameters::Bool 
    # Tip loss correction variables
    hasTipCorrection::Bool
    tipLossFunction::Union{Nothing,<:Function}
    tipLossDecayFactor::Number
    # TF for small angle of attack approximations
    smallAngles::Bool

    # Secondary (outputs from aero surface creation)
    # ----------------------------------------------
    # Flap deflection rates
    δdot::Function
    δddot::Function
    # Flag for independent flap surface
    hasIndependentFlap::Bool

end
export AeroSurface


"""
    create_AeroSurface(; kwargs...)

Aerodynamic surface constructor

# Keyword arguments
- `solver::AeroSolver` = aerodynamic solver for pitch-plunge-induced loads
- `flapLoadsSolver::FlapAeroSolver` = aerodynamic solver for flap-induced loads
- `gustLoadsSolver::GustAeroSolver` = aerodynamic solver for gust-induced loads
- `derivationMethod::DerivationMethod` = method for calculation of aerodynamic derivatives
- `airfoil::Airfoil` = airfoil section
- `c::Union{<:Function,Number}` = chord
- `Λ::Union{<:Function,Number}` = sweep angle
- `normSparPos::Union{<:Function,Float64}` = normalized position of the spar (beam reference line) on the chord
- `normFlapSpan::Union{Nothing,Vector{<:Number}}` = normalized position of the trailing-edge flap along the span (beam arclength)
- `normFlapPos::Union{Nothing,Float64}` = normalized position of the trailing-edge flap hinge on the chord
- `δIsTrimVariable::Bool` = flag for trailing-edge deflection being a trim variable
- `δ::Union{Nothing,<:Function,Number}` = trailing-edge deflection [rad]
- `flapSiteID::Union{Nothing,Int64}` = trailing-edge flap site ID
- `updateAirfoilParameters::Bool` = flag to update airfoil parameters with local airspeed
- `hasTipCorrection::Bool` = flag to employ a tip correction on aerodynamic coefficients
- `tipLossFunction::String` = respective tip loss function
- `tipLossDecayFactor::Float64` = respective tip loss factor
- `smallAngles::Bool` = flag to employ small angles approximation on the calculation of the angle of attack
"""
function create_AeroSurface(; solver::AeroSolver=Indicial(),flapLoadsSolver::FlapAeroSolver=ThinAirfoilTheory(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),derivationMethod::DerivationMethod=AD(),airfoil::Airfoil,c::Union{<:Function,Number},Λ::Union{<:Function,Number}=0.0,normSparPos::Union{<:Function,Float64},normFlapSpan::Union{Nothing,Vector{<:Number}}=nothing,normFlapPos::Union{Nothing,Float64}=nothing,δIsTrimVariable::Bool=false,δ::Union{Nothing,<:Function,Number}=nothing,flapSiteID::Union{Nothing,Int64}=nothing,updateAirfoilParameters::Bool=true,hasTipCorrection::Bool=false,tipLossFunction::Union{Nothing,<:Function}=nothing,tipLossDecayFactor::Number=Inf64,smallAngles::Bool=false)

    # Validate
    if c isa Number
        @assert c > 0 "chord must be positive"
    end
    if Λ isa Number
        @assert -π/2 < Λ < π/2 "sweep angle too large (input must be in radians and smaller than π/2)"
    end
    if normSparPos isa Number
        @assert 0 < normSparPos < 1 "normSparPos must be between 0 and 1"
    end
    if !isnothing(normFlapSpan)
        @assert !isnothing(normFlapPos) "flap span was set, but position was not"
        @assert length(normFlapSpan) == 2 "normFlapSpan must be a vector with 2 entries"
        @assert first(normFlapSpan) >= 0 "first entry of normFlapSpan must be greater to or equal to 0"
        @assert last(normFlapSpan) <= 1 "second entry of normFlapSpan must be smaller to or equal to 1"
        @assert issorted(normFlapSpan) "set normFlapSpan in ascending order"
    end 
    if !isnothing(normFlapPos)
        @assert !isnothing(normFlapSpan) "flap position was set, but span was not"
        @assert 0.5 <= normFlapPos < 1 "normFlapPos must be between 0.5 and 1"
    end
    if !isnothing(δ)
        @assert !δIsTrimVariable "flap deflection cannot be a trim variable and an input"
    end
    @assert typeof(derivationMethod) in solver.availableDerivativesMethod

    # Set flap site ID
    if !isnothing(normFlapPos) && isnothing(flapSiteID)
        flapSiteID = round(Int,100*normFlapPos)
    end

    # Update airfoil flap parameters and Theodorsen parameters with known flap position
    if !isnothing(flapSiteID)
        airfoil.flapParameters = FlapParameters(airfoil.name; flapSiteID=flapSiteID)
        if typeof(flapLoadsSolver) == ThinAirfoilTheory
            flapLoadsSolver.Th = theodorsen_flap_constants(normSparPos,normFlapPos)
        end
    end

    # Set TF for δ being user input
    δIsInput = !isnothing(δ) 

    # Set TF for δ being zero over time
    δIsZero = isnothing(δ) && !δIsTrimVariable

    # Set flap deflection and rates as functions of time
    if isnothing(δ)
        δ = t -> 0
        δdot = t -> 0
        δddot = t -> 0
    elseif δ isa Number
        δconst = deepcopy(δ)
        δ = t -> δconst
        δdot = δIsTrimVariable ? t -> 0 : t -> ForwardDiff.derivative(δ, t)
        δddot = δIsTrimVariable ? t -> 0 : t -> ForwardDiff.derivative(δdot, t)
    elseif δ isa Function
        δdot = t -> ForwardDiff.derivative(δ, t)
        δddot = t -> ForwardDiff.derivative(δdot, t)
    end

    return AeroSurface(solver=solver,flapLoadsSolver=flapLoadsSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=c,Λ=Λ,normSparPos=normSparPos,normFlapSpan=normFlapSpan,normFlapPos=normFlapPos,δIsInput=δIsInput,δIsTrimVariable=δIsTrimVariable,δIsZero=δIsZero,δ=δ,flapSiteID=flapSiteID,updateAirfoilParameters=updateAirfoilParameters,hasTipCorrection=hasTipCorrection,tipLossFunction=tipLossFunction,tipLossDecayFactor=tipLossDecayFactor,δdot=δdot,δddot=δddot,smallAngles=smallAngles,hasIndependentFlap=true)

end
export create_AeroSurface

