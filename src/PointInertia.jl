"""
@with_kw mutable struct PointInertia

PointInertia composite type

# Fields
- `elementID::Int64` = local element ID to which the point inertia is attached
- `mass::Number` = mass
- `η::Vector{Number}` = position relative to element's midpoint's reference line
- `Ixx::Number` = mass moment of inertia about the x1-axis
- `Iyy::Number` = mass moment of inertia about the x2-axis
- `Izz::Number` = mass moment of inertia about the x3-axis
- `Ixy::Number` = mass product of inertia with respect to the x1-axis and x2-axis
- `Ixz::Number` = mass product of inertia with respect to the x1-axis and x3-axis
- `Iyz::Number` = mass product of inertia with respect to the x2-axis and x3-axis
- `inertiaMatrix::Union{Matrix{Number},Matrix{Nothing}}` = mass moment of inertia matrix (tensor)
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