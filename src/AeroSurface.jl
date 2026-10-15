## AEROSURFACE COMPOSITE TYPES
# ------------------------------------------------------------------------------
# AeroSurfaceSolvers composite type
# Fields:
# - `ppSolver`: aerodynamic ppSolver for pitch-plunge-induced loads
# - `flapSolver`: aerodynamic ppSolver for flap-induced loads
# - `gustSolver`: aerodynamic ppSolver for gust-induced loads
# - `derivationMethod`: method for calculation of aerodynamic derivatives
#
struct AeroSurfaceSolvers{S1<:AeroSolver,S2<:FlapAeroSolver,S3<:GustAeroSolver}

    ppSolver::S1
    flapSolver::S2
    gustSolver::S3
    derivationMethod::DerivationMethod

end

# AeroSurface geometry composite type
# Fields:
# - `airfoil`: airfoil section
# - `γbC`: airfoil correction factor for circulatory indicial parameters in dynamic stall
# - `c`: chord
# - `Λ`: sweep angle [rad]
# - `φ`: twist angle [rad]
# - `normSparPos`: normalized position of the spar (beam reference line) on the chord
# - `normFlapPos`: normalized position of the trailing-edge flap hinge on the chord
# - `normFlapSpan`: normalized position of the trailing-edge flap along the span (beam arclength)
# - `flapParameters`: parameters of the trailing-edge flap
# - `updateAirfoilParameters`: flag to update airfoil parameters with local airspeed
# - `smallAngles`: flag indicating if small angle assumptions are used
# - `hasInducedDrag`: flag indicating if induced drag is considered
# - `hasSymmetricCounterpart`: flag indicating if the surface has a symmetric counterpart
#
struct AeroSurfaceGeometry

    airfoil::Airfoil
    γbC::Vector{Float64}

    c::Union{<:Function,Float64}
    Λ::Union{<:Function,Float64}
    φ::Union{<:Function,Float64}

    normSparPos::Union{<:Function,Float64}
    normFlapPos::Union{Nothing,Float64}
    normFlapSpan::Union{Nothing,SVector{2,Float64}}

    flapParameters::FlapParameters

    updateAirfoilParameters::Bool
    smallAngles::Bool
    hasInducedDrag::Bool
    hasSymmetricCounterpart::Bool

end

# AeroSurface tip loss composite type
# Fields:
# - `hasTipCorrection`: flag indicating if tip correction is applied
# - `tipLossFunctionIsAirspeedDependent`: flag indicating if the tip loss function depends on airspeed
# - `tipLossFunction`: tip loss function (if any)
# - `tipLossDecayFactor`: decay factor for the tip loss function
#
struct AeroSurfaceTipLoss

    hasTipCorrection::Bool
    tipLossFunctionIsAirspeedDependent::Bool
    tipLossFunction::Union{Nothing,<:Function}
    tipLossDecayFactor::Float64

end


# AeroSurface state ranges composite type
# Fields:
# - `nPitchPlungeStates`: number of pitch-plunge states
# - `nFlapStates`: number of flap states
# - `nGustStates`: number of gust states
# - `nTotalAeroStates`: total number of aerodynamic states
# - `pitchPlungeStatesRange`: range of pitch-plunge states
# - `linearPitchPlungeStatesRange`: range of linear pitch-plunge states
# - `nonlinearStatesRange`: range of nonlinear states
# - `circulatoryPitchPlungeStatesRange`: range of circulatory pitch-plunge states
# - `inertialPitchPlungeStatesRange`: range of inertial pitch-plunge states
#
struct AeroSurfaceStateRanges

    nPitchPlungeStates::Int
    nFlapStates::Int
    nGustStates::Int
    nTotalAeroStates::Int

    pitchPlungeStatesRange::Union{Nothing,UnitRange{Int}}
    circulatoryPitchPlungeStatesRange::Union{Nothing,UnitRange{Int}}
    inertialPitchPlungeStatesRange::Union{Nothing,UnitRange{Int}}
    linearPitchPlungeStatesRange::Union{Nothing,UnitRange{Int}}
    nonlinearStatesRange::Union{Nothing,UnitRange{Int}}

end


# AeroSurface assembly data composite type
# Fields
# - `DOF_δ`: global DOF of the flap trim state
# - `δFun`: trailing-edge deflection as a function of time [rad] (wrapped)
# - `δdotFun`: flap deflection rate as a function of time [rad/s] (wrapped)
# - `δddotFun`: flap deflection acceleration as a function of time [rad/s²] (wrapped)
# - `δFun_raw`: raw trailing-edge deflection function (for precompute without FunctionWrapper)
# - `δVec`: vector of flap deflection values over time [rad]
# - `δdotVec`: vector of flap deflection rate values over time [rad/s]
# - `δddotVec`: vector of flap deflection acceleration values over time [rad/s²]
# - `δIsTrimVariable`: flag for trailing-edge deflection being a trim variable
# - `δMultiplier`: flap deflection multiplier for slave surfaces
# - `area`: area of the aerodynamic surface
# - `AR`: aspect ratio of the aerodynamic surface
#
struct AeroSurfaceAssemblyData

    DOF_δ::Vector{Int}

    δFun::Function
    δdotFun::Function
    δddotFun::Function

    δVec::Vector{Float64}
    δdotVec::Vector{Float64}
    δddotVec::Vector{Float64}

    δIsTrimVariable::Bool

    δMultiplier::Float64

    area::Float64
    AR::Float64

end

# Generic constructor for any callable types
AeroSurfaceAssemblyData(δFun, δdotFun, δddotFun, δIsTrimVariable) = AeroSurfaceAssemblyData(
    Int[],
    δFun, δdotFun, δddotFun,
    [δFun(0)], [δdotFun(0)], [δddotFun(0)],
    δIsTrimVariable,
    1.0,
    0.0, 0.0
)

AeroSurfaceAssemblyData(surface, δFun, δdotFun, δddotFun, δMultiplier::Float64, δIsTrimVariable::Bool) = AeroSurfaceAssemblyData(
    surface.assemblyData.DOF_δ,
    δFun,
    δdotFun,
    δddotFun,
    surface.assemblyData.δVec, 
    surface.assemblyData.δdotVec,
    surface.assemblyData.δddotVec,
    δIsTrimVariable,
    δMultiplier,
    surface.assemblyData.area,
    surface.assemblyData.AR
)

AeroSurfaceAssemblyData(surface, area::Float64, AR::Float64) = AeroSurfaceAssemblyData(
    surface.assemblyData.DOF_δ,
    surface.assemblyData.δFun,
    surface.assemblyData.δdotFun,
    surface.assemblyData.δddotFun,
    surface.assemblyData.δVec, 
    surface.assemblyData.δdotVec,
    surface.assemblyData.δddotVec,
    surface.assemblyData.δIsTrimVariable,
    surface.assemblyData.δMultiplier,
    area,
    AR
)

AeroSurfaceAssemblyData(surface, DOF_δ::Vector{Int}) = AeroSurfaceAssemblyData(
    DOF_δ,
    surface.assemblyData.δFun,
    surface.assemblyData.δdotFun,
    surface.assemblyData.δddotFun,
    surface.assemblyData.δVec, 
    surface.assemblyData.δdotVec,
    surface.assemblyData.δddotVec,
    surface.assemblyData.δIsTrimVariable,
    surface.assemblyData.δMultiplier,
    surface.assemblyData.area,
    surface.assemblyData.AR
)

AeroSurfaceAssemblyData(surface, δVec::Vector{Float64},δdotVec::Vector{Float64},δddotVec::Vector{Float64}) = AeroSurfaceAssemblyData(
    surface.assemblyData.DOF_δ,
    surface.assemblyData.δFun,
    surface.assemblyData.δdotFun,
    surface.assemblyData.δddotFun,
    δVec, 
    δdotVec,
    δddotVec,
    surface.assemblyData.δIsTrimVariable,
    surface.assemblyData.δMultiplier,
    surface.assemblyData.area,
    surface.assemblyData.AR
)


# AeroSurface state composite type
# Fields
# - `δ`: flap deflection [rad]
# - `δdot`: flap deflection rate [rad/s]
# - `δddot`: flap deflection acceleration [rad/s²]
#
mutable struct AeroSurfaceState

    δ::Float64
    δdot::Float64
    δddot::Float64

end
AeroSurfaceState(assemblyData::AeroSurfaceAssemblyData) =
    AeroSurfaceState(assemblyData.δFun(0.0), assemblyData.δdotFun(0.0), assemblyData.δddotFun(0.0))


# AeroSurface composite type
# Fields


# - `assemblyData`: assembly data for the aerodynamic surface
# - `state`: state of the aerodynamic surface
#
@with_kw mutable struct AeroSurface{S1<:AeroSolver,S2<:FlapAeroSolver,S3<:GustAeroSolver}

    # Solvers and methods field
    solvers::AeroSurfaceSolvers{S1,S2,S3}

    # Geometry and configuration fields
    geometry::AeroSurfaceGeometry

    # Tip loss correction field
    tipLoss::AeroSurfaceTipLoss

    # State number and ranges field
    stateRanges::AeroSurfaceStateRanges   

    # Assembly-dependent fields (to be updated by Beam, Model and Problem)
    assemblyData::AeroSurfaceAssemblyData

    # State field
    state::AeroSurfaceState = AeroSurfaceState(assemblyData)

end
export AeroSurface


## PUBLIC API
# ------------------------------------------------------------------------------
"""
    create_AeroSurface(; kwargs...)

Aerodynamic surface constructor

# Keyword arguments
- `ppSolver::AeroSolver`: aerodynamic ppSolver for pitch-plunge-induced loads
- `flapSolver::FlapAeroSolver`: aerodynamic ppSolver for flap-induced loads
- `gustSolver::GustAeroSolver`: aerodynamic ppSolver for gust-induced loads
- `derivationMethod::DerivationMethod`: method for calculation of aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `c::Union{<:Function,Real}`: chord
- `Λ::Union{<:Function,Real}`: sweep angle [rad]
- `φ::Union{<:Function,Real}`: twist angle [rad]
- `normSparPos::Union{<:Function,Real}`: normalized position of the spar (beam reference line) on the chord
- `normFlapSpan::Union{Nothing,Vector{<:Real}}`: normalized position of the trailing-edge flap along the span (beam arclength)
- `normFlapPos::Union{Nothing,Real}`: normalized position of the trailing-edge flap hinge on the chord
- `δIsTrimVariable::Bool`: flag for trailing-edge deflection being a trim variable
- `δ::Union{Nothing,<:Function,Real}`: trailing-edge deflection [rad]
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with local airspeed
- `hasTipCorrection::Bool`: flag to employ a tip correction on 2D aerodynamic coefficients
- `tipLossFunctionIsAirspeedDependent::Bool`: flag for tip loss function being airspeed dependent
- `tipLossFunction::Union{Nothing,<:Function}`: corresponding tip loss function
- `tipLossDecayFactor::Real`: corresponding tip loss factor
- `smallAngles::Bool`: flag to employ small angles approximation on the calculation of the angle of attack
- `hasInducedDrag::Bool`: flag to include induced drag generated by the lifting surface
- `hasSymmetricCounterpart::Bool`: flag to indicate the surface has a symmetric counterpart (for the computation of aspect ratio and induced drag)
"""
function create_AeroSurface(; 
    ppSolver::AeroSolver=Indicial(),
    flapSolver::FlapAeroSolver=ThinAirfoilTheory(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(ppSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil,
    c::Union{<:Function,Real},
    Λ::Union{<:Function,Real}=0.0,
    φ::Union{<:Function,Real}=0.0,
    normSparPos::Union{<:Function,Real},
    normFlapSpan::Union{Nothing,Vector{<:Real}}=nothing,
    normFlapPos::Union{Nothing,Real}=nothing,
    δIsTrimVariable::Bool=false,
    δ::Union{Nothing,<:Function,Real}=nothing,updateAirfoilParameters::Bool=true,
    hasTipCorrection::Bool=false,
    tipLossFunctionIsAirspeedDependent::Bool=false,
    tipLossFunction::Union{Nothing,<:Function}=nothing,tipLossDecayFactor::Real=Inf,
    smallAngles::Bool=false,
    hasInducedDrag::Bool=false,
    hasSymmetricCounterpart::Bool=true
    )

    # Validate and set correct types of inputs
    c,Λ,φ,normSparPos,normFlapSpan,normFlapPos,tipLossDecayFactor,gustSolver = validate_aero_surface(c,Λ,φ,normSparPos,normFlapSpan,normFlapPos,tipLossFunctionIsAirspeedDependent,hasTipCorrection,tipLossFunction,δ,δIsTrimVariable,tipLossDecayFactor,ppSolver,gustSolver)

    # Get types
    S1 = typeof(ppSolver)
    S2 = typeof(flapSolver)
    S3 = typeof(gustSolver)

    # Set solvers and methods
    solvers = AeroSurfaceSolvers{S1,S2,S3}(ppSolver,flapSolver,gustSolver,derivationMethod)

    # Set geometry
    geometry = aero_surface_geometry(airfoil,c,Λ,φ,normSparPos,normFlapPos,normFlapSpan,updateAirfoilParameters,smallAngles,hasInducedDrag,hasSymmetricCounterpart,ppSolver)

    # Set tip loss correction data
    tipLoss = AeroSurfaceTipLoss(hasTipCorrection,tipLossFunctionIsAirspeedDependent,tipLossFunction,tipLossDecayFactor)

    # Set state ranges
    stateRanges = aero_surface_state_ranges(ppSolver,flapSolver,gustSolver)

    # Initialize assembly data
    assemblyData = aero_surface_assembly_data(δ,δIsTrimVariable)

    return AeroSurface{S1,S2,S3}(
        solvers=solvers,
        geometry=geometry,
        tipLoss=tipLoss,
        stateRanges=stateRanges,
        assemblyData=assemblyData,
        )

end
export create_AeroSurface


## PRIVATE API
# ------------------------------------------------------------------------------
# Validate the inputs for creating an aerodynamic surface
function validate_aero_surface(c,Λ,φ,normSparPos,normFlapSpan,normFlapPos,tipLossFunctionIsAirspeedDependent,hasTipCorrection,tipLossFunction,δ,δIsTrimVariable,tipLossDecayFactor,ppSolver,gustSolver)

    # Validate and set correct types of inputs
    if c isa Real
        @assert c > 0 "chord must be positive"
        c = Float64(c)
    end
    if Λ isa Real
        @assert -π/2 < Λ < π/2 "sweep angle too large (input modulus must be in radians and smaller than π/2)"
        Λ = Float64(Λ)
    end
    if φ isa Real
        @assert -π/4 <= φ <= π/4 "twist angle too large (input modulus must be in radians and smaller than or equal than π/4)"
        φ = Float64(φ)
    end
    if normSparPos isa Real
        @assert 0 < normSparPos < 1 "normalized spar position must be between 0 and 1"
        normSparPos = Float64(normSparPos)
    end
    if !isnothing(normFlapSpan)
        @assert !isnothing(normFlapPos) "flap span was set, but chord-normalized position was not"
        @assert length(normFlapSpan) == 2 "normFlapSpan must be a vector with 2 entries"
        @assert first(normFlapSpan) >= 0 "first entry of normFlapSpan must be greater to or equal to 0"
        @assert last(normFlapSpan) <= 1 "second entry of normFlapSpan must be smaller to or equal to 1"
        @assert issorted(normFlapSpan) "set normFlapSpan in ascending order"
        normFlapSpan = SVector(Float64.(normFlapSpan)...)
    end 
    if !isnothing(normFlapPos)
        @assert !isnothing(normFlapSpan) "chord-normalized flap position (normFlapPos) was set, but span (normFlapSpan) was not"
        @assert 0.5 <= normFlapPos <= 1 "normFlapPos must be between 0.5 and 1"
        normFlapPos = Float64(normFlapPos)
    end
    if tipLossFunctionIsAirspeedDependent
        @assert hasTipCorrection "tip loss function cannot be airspeed dependent if tip correction is disabled"
        @assert !isnothing(tipLossFunction) "tip loss function must be set if tip loss function is airspeed dependent"
    end
    if !isnothing(δ)
        @assert !δIsTrimVariable "flap deflection cannot be a trim variable and an input"
        if δ isa Function
            @assert applicable(δ, 0) "When input as a function, δ must be callable as a function of one argument (t)"
        end
    end
    tipLossDecayFactor = Float64(tipLossDecayFactor)

    # Enforce quasi-steady gust indicial solver if the pitch-plunge is quasi-steady
    if ppSolver isa QuasiSteady && gustSolver.indicialFunction != :QuasiSteady
        @warn "Enforcing quasi-steady gust indicial ppSolver because the pitch-plunge ppSolver is quasi-steady"
        gustSolver = IndicialGust(:QuasiSteady)
    end

    return c,Λ,φ,normSparPos,normFlapSpan,normFlapPos,tipLossDecayFactor,gustSolver

end

# Generate the geometry for an aerodynamic surface
function aero_surface_geometry(airfoil,c,Λ,φ,normSparPos,normFlapPos,normFlapSpan,updateAirfoilParameters,smallAngles,hasInducedDrag,hasSymmetricCounterpart,ppSolver)

    # Set flap site ID and flap parameters
    flapSiteID = !isnothing(normFlapPos) ? 
                 round(Int,100*normFlapPos) : 
                 100
    flapParameters = airfoil.database.flap.parameters[flapSiteID]

    # Airfoil-dependent correction factor for circulatory indicial parameters
    γbC = if (ppSolver isa IndicialSolver && ppSolver.nCirculatoryStates == 2) || (ppSolver isa Inflow && length(ppSolver.AC) == 2)
        airfoil.database.circulation.γbC2
    elseif (ppSolver isa IndicialSolver && ppSolver.nCirculatoryStates == 3) || (ppSolver isa Inflow && length(ppSolver.AC) == 3)
        airfoil.database.circulation.γbC3
    else
        zeros(0)
    end

    # Surface geometry data
    geometry = AeroSurfaceGeometry(airfoil,γbC,c,Λ,φ,normSparPos,normFlapPos,normFlapSpan,flapParameters,updateAirfoilParameters,smallAngles,hasInducedDrag,hasSymmetricCounterpart)

    return geometry

end


# Generate the state ranges for an aerodynamic surface
function aero_surface_state_ranges(ppSolver,flapSolver,gustSolver)

    # Number of aerodynamic states
    nPitchPlungeStates = ppSolver.nStates
    nFlapStates = flapSolver.nStates
    nGustStates = gustSolver.nStates
    nTotalAeroStates = nPitchPlungeStates + nFlapStates + nGustStates

    # Set aerodynamic states' ranges
    if ppSolver isa QuasiSteady
        circulatoryPitchPlungeStatesRange = inertialPitchPlungeStatesRange = pitchPlungeStatesRange = linearPitchPlungeStatesRange = nonlinearStatesRange = nothing
    elseif ppSolver isa Indicial
        circulatoryPitchPlungeStatesRange = 1:ppSolver.nCirculatoryStates
        inertialPitchPlungeStatesRange = ppSolver.nCirculatoryStates+1:ppSolver.nStates
        pitchPlungeStatesRange = linearPitchPlungeStatesRange = 1:ppSolver.nStates
        nonlinearStatesRange = nothing
    elseif ppSolver isa Inflow
        circulatoryPitchPlungeStatesRange = 1:ppSolver.nInflowStates
        inertialPitchPlungeStatesRange = ppSolver.nInflowStates+1:ppSolver.nStates
        pitchPlungeStatesRange = linearPitchPlungeStatesRange = 1:ppSolver.nStates
        nonlinearStatesRange = nothing
    elseif ppSolver isa DynamicStallSolver
        circulatoryPitchPlungeStatesRange = 1:ppSolver.nCirculatoryStates
        inertialPitchPlungeStatesRange = ppSolver.nCirculatoryStates+1:ppSolver.nCirculatoryStates+ppSolver.nInertialStates
        pitchPlungeStatesRange = 1:ppSolver.nStates
        linearPitchPlungeStatesRange = 1:ppSolver.nCirculatoryStates+ppSolver.nInertialStates
        nonlinearStatesRange = ppSolver.nCirculatoryStates+ppSolver.nInertialStates+1:ppSolver.nStates
    end

    # State ranges object
    stateRanges = AeroSurfaceStateRanges(
        nPitchPlungeStates,
        nFlapStates,
        nGustStates,
        nTotalAeroStates,
        pitchPlungeStatesRange,
        circulatoryPitchPlungeStatesRange,
        inertialPitchPlungeStatesRange,
        linearPitchPlungeStatesRange,
        nonlinearStatesRange,
    )

    return stateRanges

end


# Initialize the assembly data for an aerodynamic surface
function aero_surface_assembly_data(δ,δIsTrimVariable)

    # Initialize flap deflection and rates as functions of time (these may be updated later by the Model if current surface becomes a slave surface through a FlapLink)
    if isnothing(δ)
        δFun     = t -> 0.0
        δdotFun  = t -> 0.0
        δddotFun = t -> 0.0
    elseif δ isa Real
        δval = Float64(δ)
        δFun     = t -> δval
        δdotFun  = t -> 0.0
        δddotFun = t -> 0.0
    elseif δ isa Function
        δFun     = δ 
        δdotFun  = t -> ForwardDiff.derivative(δFun, t)
        δddotFun = t -> ForwardDiff.derivative(δdotFun, t)
    end

    # Initialize assembly data with flap deflection functions
    assemblyData = AeroSurfaceAssemblyData(δFun,δdotFun,δddotFun,δIsTrimVariable)

    return assemblyData

end