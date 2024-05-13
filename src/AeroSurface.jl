"""
@with_kw mutable struct AeroSurface

    AeroSurface composite type

# Fields
- solver::AeroSolver
- flapLoadsSolver::FlapAeroSolver
- gustLoadsSolver::GustAeroSolver
- airfoil::Airfoil
- c::Union{<:Function,Number}
- Λ::Union{<:Function,Number}
- normSparPos::Union{<:Function,Float64}
- normFlapSpan::Union{Nothing,Vector{<:Number}}
- normFlapPos::Union{Nothing,Float64}
- flapSiteID::Union{Nothing,Int64}
- updateAirfoilParameters::Bool
- hasTipCorrection::Bool
- tipLossFunction::String
- tipLossDecayFactor::Float64
- αᵣ::Number
- δIsInput::Bool
- δdot::Function
- δddot::Function
"""
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

    # Secondary (outputs from aero surface creation)
    # ----------------------------------------------
    # TF for tip loss function being an user input
    tipLossFunctionIsInput::Bool
    # TF for flap deflection being an user input
    δIsInput::Bool
    # Flap deflection rates
    δdot::Function
    δddot::Function

end
export AeroSurface


# Constructor
function create_AeroSurface(;solver::AeroSolver=Indicial(),flapLoadsSolver::FlapAeroSolver=ThinAirfoilTheory(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),derivationMethod::DerivationMethod=AD(),airfoil::Airfoil,c::Union{<:Function,Number},Λ::Union{<:Function,Number}=0.0,normSparPos::Union{<:Function,Float64},normFlapSpan::Union{Nothing,Vector{<:Number}}=nothing,normFlapPos::Union{Nothing,Float64}=nothing,δ::Union{Nothing,<:Function,Number}=nothing,flapSiteID::Union{Nothing,Int64}=nothing,updateAirfoilParameters::Bool=true,hasTipCorrection::Bool=false,tipLossFunction::Union{Nothing,<:Function}=nothing,tipLossDecayFactor::Number=Inf64)

    # Validate
    x1 = LinRange(0,1,101)
    if c isa Number
        @assert c > 0 "chord must be positive"
    else
        @assert all(x-> x>0, c.(x1)) "chord must be positive over entire span"
    end
    if Λ isa Number
        @assert -π/2 < Λ < π/2 "sweep angle too large (input must be in radians and smaller than π/2)"
    else
        @assert all(x-> -π/2 < x < π/2, Λ.(x1)) "sweep angle too large (input must be in radians and smaller than π/2)"
    end
    if normSparPos isa Number
        @assert 0 < normSparPos < 1 "normSparPos must be between 0 and 1"
    else
        @assert all(x-> 0<x<1, normSparPos.(x1)) "normSparPos must be between 0 and 1 over entire span"
    end
    if !isnothing(normFlapSpan)
        @assert length(normFlapSpan) == 2 "normFlapSpan must be a vector with 2 entries"
        @assert first(normFlapSpan) >= 0 "first entry of normFlapSpan must be greater to or equal to 0"
        @assert last(normFlapSpan) <= 1 "second entry of normFlapSpan must be smaller to or equal to 1"
        @assert issorted(normFlapSpan) "set normFlapSpan in ascending order"
    end 
    if typeof(flapLoadsSolver) == TableLookup
        @assert !isnothing(flapSiteID) "define a flapSiteID for flap loads solver as TableLookup"
    end
    if !isnothing(normFlapPos)
        @assert 0.5 <= normFlapPos < 1 "normFlapPos must be between 0.5 and 1"
        flapSiteID = round(Int,100*normFlapPos)
    end
    @assert typeof(derivationMethod) in solver.availableDerivativesMethod

    # Update airfoil parameters with known flap position in the case of flap loads by table lookup
    if typeof(flapLoadsSolver) == TableLookup
        airfoil = create_flapped_Airfoil(name=airfoil,flapSiteID=flapSiteID)
    end

    # Set TF for tip loss function being an user input
    tipLossFunctionIsInput = !isnothing(tipLossFunction) ? true : false

    # Set TF for δ being user input
    δIsInput = !isnothing(δ) ? true : false

    # Set flap deflection and rates as functions of time
    if isnothing(δ)
        δ = t -> 0
        δdot = t -> 0
        δddot = t -> 0
    elseif δ isa Number
        δconst = deepcopy(δ)
        δ = t -> δconst
        δdot = t -> ForwardDiff.derivative(δ, t)
        δddot = t -> ForwardDiff.derivative(δdot, t)
    elseif δ isa Function
        δdot = t -> ForwardDiff.derivative(δ, t)
        δddot = t -> ForwardDiff.derivative(δdot, t)
    end

    return AeroSurface(solver=solver,flapLoadsSolver=flapLoadsSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=c,Λ=Λ,normSparPos=normSparPos,normFlapSpan=normFlapSpan,normFlapPos=normFlapPos,δ=δ,flapSiteID=flapSiteID,updateAirfoilParameters=updateAirfoilParameters,hasTipCorrection=hasTipCorrection,tipLossFunction=tipLossFunction,tipLossDecayFactor=tipLossDecayFactor,tipLossFunctionIsInput=tipLossFunctionIsInput,δIsInput=δIsInput,δdot=δdot,δddot=δddot)

end
export create_AeroSurface

