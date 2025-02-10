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
    availableDerivativesMethod::Vector{Type{<:DerivationMethod}}

    function QuasiSteady()
        nStates = 0
        availableDerivativesMethod = [AD,FD]
        return new(nStates,availableDerivativesMethod)
    end

end
export QuasiSteady


"""

    Indicial AeroSolver composite type

"""
struct Indicial <: AeroSolver

    nStates::Int64
    availableDerivativesMethod::Vector{Type{<:DerivationMethod}}
    AW::Vector{Float64}
    bW::Vector{Float64}
    bWMat::Matrix{Float64}

    function Indicial()
        nStates = 2
        availableDerivativesMethod = [AD,FD]
        AW = [0.165; 0.335]
        bW = [0.0455; 0.3]
        bWMat = diagm(bW)
        return new(nStates,availableDerivativesMethod,AW,bW,bWMat)
    end

end
export Indicial


"""

    Incompressible modified Beddoes-Leishman AeroSolver composite type

"""
struct BLi <: AeroSolver

    nStates::Int64
    availableDerivativesMethod::Vector{Type{<:DerivationMethod}}
    AW::Vector{Float64}
    bW::Vector{Float64}
    bWMat::Matrix{Float64}

    function BLi()
        nStates = 8
        availableDerivativesMethod = [AD,FD]
        AW = [0.165; 0.335]
        bW = [0.0455; 0.3]
        bWMat = diagm(bW)
        return new(nStates,availableDerivativesMethod,AW,bW,bWMat)
    end

end
export BLi


"""

    Original Beddoes-Leishman AeroSolver composite type

"""
struct BLo <: AeroSolver

    incompressibleInertialLoads::Bool
    nStates::Int64
    availableDerivativesMethod::Vector{Type{<:DerivationMethod}}
    a::Vector{Float64}
    b::Vector{Float64}
    b12Mat::Matrix{Float64}
    AW::Vector{Float64}
    bW::Vector{Float64}
    bWMat::Matrix{Float64}
    a1b1a2b2::Vector{Float64}

    function BLo(; incompressibleInertialLoads::Bool=false)
        nStates = incompressibleInertialLoads ? 7 : 13
        availableDerivativesMethod = [AD,FD]
        a = [0.3; 0.7; 1.5; -0.5]
        b = [0.14; 0.53; 0.25; 0.1; 0.5]
        bMat = diagm(b)
        AW = [0.165; 0.335]
        bW = [0.0455; 0.3]
        bWMat = diagm(bW)
        a1b1a2b2 = [a[1]*b[1]; a[2]*b[2]]
        return new(incompressibleInertialLoads,nStates,availableDerivativesMethod,a,b,bMat,AW,bW,bWMat,a1b1a2b2)
    end

end
export BLo


"""

    Inflow AeroSolver composite type

"""
struct Inflow <: AeroSolver

    nStates::Int64
    availableDerivativesMethod::Vector{Type{<:DerivationMethod}}
    AₚInv::Matrix{Float64}
    AₚInvcₚ::Vector{Float64}
    bₚ::Vector{Float64}

    function Inflow(nStates::Int64=6)
        @assert 2 <= nStates <= 8 "select between 2 and 8 inflow states"
        availableDerivativesMethod = [AD,FD]
        AₚInv,bₚ,cₚ = inflow_arrays(nStates)
        AₚInvcₚ = AₚInv*cₚ
        return new(nStates,availableDerivativesMethod,AₚInv,AₚInvcₚ,bₚ)
    end

end
export Inflow


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
    AWf::Vector{Float64}
    bWf::Vector{Float64}
    bWfMat::Matrix{Float64}
    Th::Vector{Float64} 

    function ThinAirfoilTheory()
        nStates = 2
        AWf = [0.165; 0.335]
        bWf = [0.0455; 0.3]
        bWfMat = diagm(bWf)
        Th = theodorsen_flap_constants(0.25,1.0)
        return new(nStates,AWf,bWf,bWfMat,Th)
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
    Th[15] = Th[1]-Th[8]-(normFlapPos-normSparPos)*Th[4]+Th[11]/2
    Th[16] = Th[7]+Th[1]*(normFlapPos-normSparPos)
    Th[17] = Th[11]/(2π)
    Th[18] = Th[10]/π

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
        AGbG = AG .* bG
        bGMat = diagm(bG)

        return new(indicialFunctionName,nStates,AG,bG,AGbG,bGMat)
    end

end
export IndicialGust