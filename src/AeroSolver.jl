abstract type AeroSolver end
abstract type FlapAeroSolver end
abstract type GustAeroSolver end
abstract type DerivationMethod end

"""

    AD (Automatic Differentiation) DerivationMethod composite type

"""
struct AD <: DerivationMethod

end
export AD


"""

    FD (Finite Differences) DerivationMethod composite type

# Fields
- method::FiniteDifferenceMethod
"""
@with_kw mutable struct FD <: DerivationMethod

    method::FiniteDifferenceMethod = central_fdm(3,1)

    # Constructor with 'nothing' argument
    function FD(null::Nothing)
        return new(central_fdm(3,1))
    end

    # Constructor with method and options
    function FD(FDFun::Function; nPoints::Int,adapt::Int=1)

        @assert typeof(FDFun) in [typeof(central_fdm),typeof(forward_fdm),typeof(backward_fdm)]
        @assert nPoints > 1
        @assert adapt >= 1

        method = FDFun(nPoints,1,adapt=adapt)

        return new(method)
    end
end
export FD


"""

    QuasiSteady AeroSolver composite type

"""
struct QuasiSteady <: AeroSolver

    name::String
    nStates::Int

    function QuasiSteady()

        # Name
        name = "QS"

        # Number of states per element
        nStates = 0

        return new(name,nStates)
    end

end
export QuasiSteady


"""

    Indicial AeroSolver composite type

"""
struct Indicial <: AeroSolver

    name::String
    circulatoryIndicialFunction::String
    incompressibleInertialLoads::Bool
    nCirculatoryStates::Int
    nInertialStates::Int
    nStates::Int
    AC::Vector{Float64}
    bC::Vector{Float64}
    AI::Vector{Float64}
    bI::Vector{Float64}
    bCMat::Matrix{Float64}

    function Indicial(; circulatoryIndicialFunction::String="Wagner",incompressibleInertialLoads::Bool=true)

        # Validate
        @assert circulatoryIndicialFunction in ["Beddoes","Wagner","Jose"]

        # Name
        name = "Indicial"

        # Circulatory indicial parameters
        AC,bC = circulatory_indicial_parameters(circulatoryIndicialFunction)
        bCMat = diagm(bC)

        # Number of states per element
        nCirculatoryStates = length(AC)
        nInertialStates = incompressibleInertialLoads ? 0 : 8
        nStates = nCirculatoryStates + nInertialStates

        # Inertial indicial parameters
        AI = [1.5; -0.5; 1.0]
        bI = [0.25; 0.1; 5.0]

        return new(name,circulatoryIndicialFunction,incompressibleInertialLoads,nCirculatoryStates,nInertialStates,nStates,AC,bC,AI,bI,bCMat)
    end

end
export Indicial


"""

    Incompressible modified Beddoes-Leishman AeroSolver composite type

"""
struct BLi <: AeroSolver

    name::String
    circulatoryIndicialFunction::String
    incompressibleInertialLoads::Bool
    nLinearCirculatoryStates::Int
    nInertialStates::Int
    nNonlinearCirculatoryStates::Int
    nStates::Int
    AC::Vector{Float64}
    bC::Vector{Float64}
    AI::Vector{Float64}
    bI::Vector{Float64}
    bCMat::Matrix{Float64}

    function BLi(; circulatoryIndicialFunction::String="Beddoes",incompressibleInertialLoads::Bool=true)

        # Validate
        @assert circulatoryIndicialFunction in ["Beddoes","Wagner","Jose"]

        # Name
        name = "BLi"

        # Circulatory indicial parameters
        AC,bC = circulatory_indicial_parameters(circulatoryIndicialFunction)
        bCMat = diagm(bC)

        # Number of states per element
        nLinearCirculatoryStates = length(AC)
        nInertialStates = incompressibleInertialLoads ? 0 : 8
        nNonlinearCirculatoryStates = 6
        nStates = nLinearCirculatoryStates + nInertialStates + nNonlinearCirculatoryStates

        # Inertial indicial parameters
        AI = [1.5; -0.5; 1.0]
        bI = [0.25; 0.1; 5.0]

        return new(name,circulatoryIndicialFunction,incompressibleInertialLoads,nLinearCirculatoryStates,nInertialStates,nNonlinearCirculatoryStates,nStates,AC,bC,AI,bI,bCMat)
    end

end
export BLi


"""

    Original Beddoes-Leishman AeroSolver composite type

"""
struct BLo <: AeroSolver

    name::String
    circulatoryIndicialFunction::String
    incompressibleInertialLoads::Bool
    nLinearCirculatoryStates::Int
    nInertialStates::Int
    nNonlinearCirculatoryStates::Int
    nStates::Int
    AC::Vector{Float64}
    bC::Vector{Float64}
    AI::Vector{Float64}
    bI::Vector{Float64}
    bCMat::Matrix{Float64}

    function BLo(; circulatoryIndicialFunction::String="Beddoes",incompressibleInertialLoads::Bool=false)

        # Validate
        @assert circulatoryIndicialFunction in ["Beddoes","Wagner","Jose"]

        # Name
        name = "BLo"

        # Circulatory indicial parameters
        AC,bC = circulatory_indicial_parameters(circulatoryIndicialFunction)
        bCMat = diagm(bC)

        # Number of states per element
        nLinearCirculatoryStates = length(AC)
        nInertialStates = incompressibleInertialLoads ? 0 : 8
        nNonlinearCirculatoryStates = 5
        nStates = nLinearCirculatoryStates + nInertialStates + nNonlinearCirculatoryStates

        # Inertial indicial parameters
        AI = [1.5; -0.5; 1.0]
        bI = [0.25; 0.1; 5.0]

        return new(name,circulatoryIndicialFunction,incompressibleInertialLoads,nLinearCirculatoryStates,nInertialStates,nNonlinearCirculatoryStates,nStates,AC,bC,AI,bI,bCMat)
    end

end
export BLo


"""

    Inflow AeroSolver composite type

"""
struct Inflow <: AeroSolver

    name::String
    circulatoryIndicialFunction::String
    incompressibleInertialLoads::Bool
    nInflowStates::Int
    nInertialStates::Int
    nStates::Int
    AₚInv::Matrix{Float64}
    AₚInvcₚ::Vector{Float64}
    bₚ::Vector{Float64}
    AC::Vector{Float64}
    bC::Vector{Float64}
    AI::Vector{Float64}
    bI::Vector{Float64}

    function Inflow(; circulatoryIndicialFunction::String="Jose", incompressibleInertialLoads::Bool=true, nInflowStates::Int=6)

        # Validate number of states per element
        @assert nInflowStates in [4,6,8] "select nInflowStates as 4, 6 or 8"

        # Name
        name = "Inflow"

        # Inflow arrays
        AₚInv,bₚ,cₚ = inflow_arrays(nInflowStates)
        AₚInvcₚ = AₚInv*cₚ

        # Number of states per element
        nInertialStates = incompressibleInertialLoads ? 0 : 8
        nStates = nInflowStates + nInertialStates

        # Circulatory indicial parameters (needed only to compute inertial response time scales)
        AC,bC = circulatory_indicial_parameters(circulatoryIndicialFunction)

        # Inertial indicial parameters
        AI = [1.5; -0.5; 1.0]
        bI = [0.25; 0.1; 5.0]

        return new(name,circulatoryIndicialFunction,incompressibleInertialLoads,nInflowStates,nInertialStates,nStates,AₚInv,AₚInvcₚ,bₚ,AC,bC,AI,bI)
    end

end
export Inflow


# Gets the circulatory indicial parameters given the indicial function
function circulatory_indicial_parameters(circulatoryIndicialFunction)
    ACbC = Dict(
        "Beddoes" => ([0.3; 0.7], [0.14; 0.53]),
        "Wagner"  => ([0.165; 0.335], [0.0455; 0.3]),
        "Jose"    => ([0.3493; 0.6507], [0.0984; 0.7759])
    )
    return get(ACbC, circulatoryIndicialFunction, nothing)
end


# Computes the fixed state arrays from Peters' inflow theory
function inflow_arrays(N::Int)

    # D matrix
    D = zeros(N,N)
    D[1,2] = -1/2
    D[end,end-1] = 1/(2*N)
    for n=2:N-1
        D[n,n+1] = -1/(2*n)
        D[n,n-1] = 1/(2*n)
    end

    # b vector
    b = zeros(N)
    for n=1:N-1
        F = factorial(N+n-1)/factorial(N-n-1)/(factorial(n))^2
        b[n] = (-1)^(n-1)*F
    end
    b[N] = (-1)^(N+1) 

    # c vector
    c = zeros(N)
    for n=1:N
        c[n] = 2/n
    end

    # d vector
    d = zeros(N)
    d[1] = 1/2

    # A matrix and its inverse
    A = D + d*b' + c*d' + 1/2*c*b'
    AInv = inv(A)

    return AInv,b,c
end


"""

    TableLookup FlapAeroSolver composite type

"""
struct TableLookup <: FlapAeroSolver

    name::String
    nStates::Int

    function TableLookup()

        # Name
        name = "TableLookup"

        # Number of states per element
        nStates = 0

        return new(name,nStates)
    end

end
export TableLookup


"""

    ThinAirfoilTheory FlapAeroSolver composite type

"""
mutable struct ThinAirfoilTheory <: FlapAeroSolver

    circulatoryIndicialFunction::String
    name::String
    nStates::Int
    ACf::Vector{Float64}
    bCf::Vector{Float64}
    bCfMat::Matrix{Float64}

    function ThinAirfoilTheory(circulatoryIndicialFunction::String="Wagner")

        # Name
        name = "ThinAirfoilTheory"

        # Circulatory indicial parameters
        ACf,bCf = circulatory_indicial_parameters(circulatoryIndicialFunction)
        bCfMat = diagm(bCf)

        # Number of states per element
        nStates = length(bCf)

        return new(circulatoryIndicialFunction,name,nStates,ACf,bCf,bCfMat)
    end

end
export ThinAirfoilTheory


"""

    IndicialGust GustAeroSolver composite type

"""
struct IndicialGust <: GustAeroSolver

    name::String
    indicialFunctionName::String
    nStates::Int
    AG::Vector{Float64}
    bG::Vector{Float64}
    AGbG::Vector{Float64}
    bGMat::Matrix{Float64}

    function IndicialGust(indicialFunctionName::String)

        # Gust solver name
        name = string("IndicialGust-",indicialFunctionName)

        # Set number of states and indicial parameters arrays according to indicial function
        if indicialFunctionName == "QuasiSteady"
            AG = [1.0]
            bG = [1e3]
        elseif indicialFunctionName == "Kussner"
            AG = [0.5; 0.5]
            bG = [0.13; 1.0]
        elseif indicialFunctionName == "Berci&Righi"
            AG = [0.3694; 0.0550; 0.2654; 0.1829; 0.0861; 0.0412]
            bG = [0.3733; 0.0179; 0.1096; 1.6003; 10.428; 170.93]
        else
            error("Indicial function unavailable")
        end

        # Number of gust states
        nStates = length(bG)

        # Useful functions of the indicial parameters arrays
        AGbG = AG .* bG
        bGMat = diagm(bG)

        return new(name,indicialFunctionName,nStates,AG,bG,AGbG,bGMat)
    end

end
export IndicialGust