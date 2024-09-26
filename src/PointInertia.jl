"""
    PointInertia composite type

# Fields
- `elementID::Int64` = local element ID to which the point inertia is attached
- `mass::Real` = mass
- `η::Vector{Real}` = position relative to element's midpoint's reference line
- `Ixx::Real` = mass moment of inertia about the x1-axis
- `Iyy::Real` = mass moment of inertia about the x2-axis
- `Izz::Real` = mass moment of inertia about the x3-axis
- `Ixy::Real` = mass product of inertia with respect to the x1-axis and x2-axis
- `Ixz::Real` = mass product of inertia with respect to the x1-axis and x3-axis
- `Iyz::Real` = mass product of inertia with respect to the x2-axis and x3-axis
- `inertiaMatrix::Union{Matrix{Real},Matrix{Nothing}}` = mass moment of inertia matrix (tensor)
"""
@with_kw mutable struct PointInertia

    elementID::Int64
    η::Vector{<:Real} = zeros(3)
    mass::Real = 0
    Iyy::Real = 0
    Izz::Real = 0
    Ixx::Real = Iyy+Izz
    Ixy::Real = 0
    Ixz::Real = 0
    Iyz::Real = 0
    inertiaMatrix::Union{Matrix{<:Real},Matrix{Nothing}} = fill(nothing,3,3)

    # Constructor
    function PointInertia(elementID::Int64,η::Vector{<:Real},mass::Real,Iyy::Real,Izz::Real,Ixx::Real,Ixy::Real,Ixz::Real,Iyz::Real,inertiaMatrix::Union{Matrix{<:Real},Matrix{Nothing}})

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