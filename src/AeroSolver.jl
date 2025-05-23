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
    function FD(FDFun::Function; nPoints::Int64,adapt::Int64=1)

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

    nStates::Int64

    function QuasiSteady()

        # Number of states per element
        nStates = 0

        return new(nStates)
    end

end
export QuasiSteady


"""

    Indicial AeroSolver composite type

"""
struct Indicial <: AeroSolver

    circulatoryIndicialFunction::String
    incompressibleInertialLoads::Bool
    nCirculatoryStates::Int64
    nInertialStates::Int64
    nStates::Int64
    AC::Vector{Float64}
    bC::Vector{Float64}
    AI::Vector{Float64}
    bI::Vector{Float64}
    bCMat::Matrix{Float64}

    function Indicial(; circulatoryIndicialFunction::String="Wagner",incompressibleInertialLoads::Bool=true)

        # Validate
        @assert circulatoryIndicialFunction in ["Beddoes","Wagner","Jose"]

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

        return new(circulatoryIndicialFunction,incompressibleInertialLoads,nCirculatoryStates,nInertialStates,nStates,AC,bC,AI,bI,bCMat)
    end

end
export Indicial


"""

    Incompressible modified Beddoes-Leishman AeroSolver composite type

"""
struct BLi <: AeroSolver

    circulatoryIndicialFunction::String
    incompressibleInertialLoads::Bool
    nLinearCirculatoryStates::Int64
    nInertialStates::Int64
    nNonlinearCirculatoryStates::Int64
    nStates::Int64
    AC::Vector{Float64}
    bC::Vector{Float64}
    AI::Vector{Float64}
    bI::Vector{Float64}
    bCMat::Matrix{Float64}

    function BLi(; circulatoryIndicialFunction::String="Beddoes",incompressibleInertialLoads::Bool=true)

        # Validate
        @assert circulatoryIndicialFunction in ["Beddoes","Wagner","Jose"]

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

        return new(circulatoryIndicialFunction,incompressibleInertialLoads,nLinearCirculatoryStates,nInertialStates,nNonlinearCirculatoryStates,nStates,AC,bC,AI,bI,bCMat)
    end

end
export BLi


"""

    Original Beddoes-Leishman AeroSolver composite type

"""
struct BLo <: AeroSolver

    circulatoryIndicialFunction::String
    incompressibleInertialLoads::Bool
    nLinearCirculatoryStates::Int64
    nInertialStates::Int64
    nNonlinearCirculatoryStates::Int64
    nStates::Int64
    AC::Vector{Float64}
    bC::Vector{Float64}
    AI::Vector{Float64}
    bI::Vector{Float64}
    bCMat::Matrix{Float64}

    function BLo(; circulatoryIndicialFunction::String="Beddoes",incompressibleInertialLoads::Bool=false)

        # Validate
        @assert circulatoryIndicialFunction in ["Beddoes","Wagner","Jose"]

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

        return new(circulatoryIndicialFunction,incompressibleInertialLoads,nLinearCirculatoryStates,nInertialStates,nNonlinearCirculatoryStates,nStates,AC,bC,AI,bI,bCMat)
    end

end
export BLo


"""

    Inflow AeroSolver composite type

"""
struct Inflow <: AeroSolver

    circulatoryIndicialFunction::String
    incompressibleInertialLoads::Bool
    nInflowStates::Int64
    nInertialStates::Int64
    nStates::Int64
    AₚInv::Matrix{Float64}
    AₚInvcₚ::Vector{Float64}
    bₚ::Vector{Float64}
    AC::Vector{Float64}
    bC::Vector{Float64}
    AI::Vector{Float64}
    bI::Vector{Float64}

    function Inflow(; circulatoryIndicialFunction::String="Jose", incompressibleInertialLoads::Bool=true, nInflowStates::Int64=6)

        # Validate number of states per element
        @assert nInflowStates in [4,6,8] "select nInflowStates as 4, 6 or 8"

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

        return new(circulatoryIndicialFunction,incompressibleInertialLoads,nInflowStates,nInertialStates,nStates,AₚInv,AₚInvcₚ,bₚ,AC,bC,AI,bI)
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
function inflow_arrays(N::Int64)

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

    nStates::Int64

    function TableLookup()

        # Number of states per element
        nStates = 0

        return new(nStates)
    end

end
export TableLookup


"""

    ThinAirfoilTheory FlapAeroSolver composite type

"""
mutable struct ThinAirfoilTheory <: FlapAeroSolver

    nStates::Int64
    ACf::Vector{Float64}
    bCf::Vector{Float64}
    bCfMat::Matrix{Float64}
    Th::Vector{Float64} 

    function ThinAirfoilTheory()

        # Number of states per element
        nStates = 2

        # Circulatory indicial parameters (Wagner)
        ACf,bCf = circulatory_indicial_parameters("Wagner")
        bCfMat = diagm(bCf)

        # Initialize Theodorsen's flap constants assuming a pitch axis at 1/4-chord and a flap axis at 3/4-chord (updated later on model creation)
        Th = theodorsen_flap_constants(1/4,3/4)

        return new(nStates,ACf,bCf,bCfMat,Th)
    end

end
export ThinAirfoilTheory


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


"""

    IndicialGust GustAeroSolver composite type

"""
struct IndicialGust <: GustAeroSolver

    indicialFunctionName::String
    nStates::Int64
    AG::Vector{Float64}
    bG::Vector{Float64}
    AGbG::Vector{Float64}
    bGMat::Matrix{Float64}

    function IndicialGust(indicialFunctionName::String)

        # Set number of states and indicial parameters arrays according to indicial function
        if indicialFunctionName == "Kussner"
            nStates = 2
            AG = [0.5; 0.5]
            bG = [0.13; 1.0]
        elseif indicialFunctionName == "Berci&Righi"
            nStates = 6
            AG = [0.3694; 0.0550; 0.2654; 0.1829; 0.0861; 0.0412]
            bG = [0.3733; 0.0179; 0.1096; 1.6003; 10.428; 170.93]
        else
            error("Indicial function unavailable")
        end

        # Useful functions of the indicial parameters arrays
        AGbG = AG .* bG
        bGMat = diagm(bG)

        return new(indicialFunctionName,nStates,AG,bG,AGbG,bGMat)
    end

end
export IndicialGust