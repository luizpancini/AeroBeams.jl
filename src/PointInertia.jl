"""
@with_kw mutable struct PointInertia

PointInertia composite type

# Fields
- elementID::Int64
- mass::Number 
- η::Vector{Number} = position relative to element's midpoint
- Ixx::Number 
- Iyy::Number 
- Izz::Number 
- Ixy::Number 
- Ixz::Number 
- Iyz::Number 
- inertiaMatrix::Union{Matrix{Number},Matrix{Nothing}} 
"""
@with_kw mutable struct PointInertia

    elementID::Int64
    η::Vector{<:Number} = zeros(3)
    mass::Number = 0
    Iyy::Number = 0
    Izz::Number = 0
    Ixx::Number = Iyy+Izz
    Ixy::Number = 0
    Ixz::Number = 0
    Iyz::Number = 0
    inertiaMatrix::Union{Matrix{<:Number},Matrix{Nothing}} = fill(nothing,3,3)

    # Constructor
    function PointInertia(elementID::Int64,η::Vector{<:Number},mass::Number,Iyy::Number,Izz::Number,Ixx::Number,Ixy::Number,Ixz::Number,Iyz::Number,inertiaMatrix::Union{Matrix{<:Number},Matrix{Nothing}})

        # Validate inputs
        @assert mass >= 0
        @assert Ixx >= 0
        @assert Iyy >= 0
        @assert Izz >= 0
        @assert length(η) == 3
        @assert size(inertiaMatrix) == (3, 3)

        # Initialize inertia matrix if it was not input
        if all(x -> x === nothing, inertiaMatrix)
            inertiaMatrix = [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]
        end

        return new(elementID,η,mass,Iyy,Izz,Ixx,Ixy,Ixz,Iyz,inertiaMatrix)

    end
end
export PointInertia