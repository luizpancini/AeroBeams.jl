"""
@with_kw mutable struct PointInertia

PointInertia composite type

# Fields
- elementLocalID::Int64
- mass::Float64 
- η::Vector{Float64} = position relative to elementLocalID's midpoint
- ηtilde::Matrix{Float64} = skew-symmetric matrix of η
- Ixx::Float64 
- Iyy::Float64 
- Izz::Float64 
- Ixy::Float64 
- Ixz::Float64 
- Iyz::Float64 
- inertiaMatrix::Union{Matrix{Float64},Matrix{Nothing}} 
"""
@with_kw mutable struct PointInertia
    # Fields
    elementLocalID::Int64
    mass::Float64
    η::Vector{Float64} = zeros(3)
    ηtilde::Matrix{Float64} = zeros(3,3)
    Ixx::Float64 = 0.0
    Iyy::Float64 = 0.0
    Izz::Float64 = 0.0
    Ixy::Float64 = 0.0
    Ixz::Float64 = 0.0
    Iyz::Float64 = 0.0
    inertiaMatrix::Union{Matrix{Float64},Matrix{Nothing}} = fill(nothing,3,3)

    # Constructor
    function PointInertia(elementLocalID::Int64,mass::Float64,η::Vector{Float64},ηtilde::Matrix{Float64},Ixx::Float64,Iyy::Float64,Izz::Float64,Ixy::Float64,Ixz::Float64,Iyz::Float64,inertiaMatrix::Union{Matrix{Float64},Matrix{Nothing}})

        # Validate inputs
        @assert mass > 0
        @assert length(η) == 3
        @assert size(inertiaMatrix) == (3, 3)

        # Get skew-symmetric matrix of position
        ηtilde = tilde(η)

        # Initialize inertia matrix if it was not input
        if all(x -> x === nothing, inertiaMatrix)
            inertiaMatrix = [Ixx Ixy Ixz; Ixy Iyy Iyz; Iyz Ixz Izz]
        end

        return new(elementLocalID,mass,η,ηtilde,Ixx,Iyy,Izz,Ixy,Ixz,Iyz,inertiaMatrix)

    end
end
export PointInertia