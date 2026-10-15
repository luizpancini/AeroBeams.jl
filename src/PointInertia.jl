# PointInertia composite type
# Fields
# - `mass`: mass of the point inertia
# - `η`: position vector of the point inertia relative to the element's midpoint, resolved in the local basis of the beam element
# - `inertiaMatrix`: mass moment of inertia matrix
#
@with_kw struct PointInertia

    mass::Float64
    η::SVector{3,Float64}
    inertiaMatrix::SMatrix{3,3,Float64,9}

end


"""
    create_PointInertia(; kwargs...)

Point inertia constructor (a point mass with mass moment of inertia to be attached to a beam element)

# Keyword arguments
- `mass::Real`: mass
- `η::AbstractVector{<:Real}`: position vector relative to element's midpoint's reference line
- `Ixx::Real`: mass moment of inertia about the x1-axis
- `Iyy::Real`: mass moment of inertia about the x2-axis
- `Izz::Real`: mass moment of inertia about the x3-axis
- `Ixy::Real`: mass product of inertia with respect to the x1-axis and x2-axis
- `Ixz::Real`: mass product of inertia with respect to the x1-axis and x3-axis
- `Iyz::Real`: mass product of inertia with respect to the x2-axis and x3-axis
- `inertiaMatrix::Union{Missing,AbstractMatrix{<:Real}}`: mass moment of inertia matrix
"""
function create_PointInertia(; 
    mass::Real=0,
    η::AbstractVector{<:Real}=zeros(3),
    Iyy::Real=0,
    Izz::Real=0,
    Ixx::Real=Iyy+Izz,
    Ixy::Real=0,
    Ixz::Real=0,
    Iyz::Real=0,
    inertiaMatrix::Union{Missing,AbstractMatrix{<:Real}}=missing,
    )

    # Validate inputs
    @assert mass >= 0 "mass must be non-negative"
    @assert all(x -> x >= 0, [Ixx, Iyy, Izz]) "mass moments of inertia must be non-negative"
    @assert length(η) == 3 "η must be a 3-element vector"

    # Set inertia matrix
    if ismissing(inertiaMatrix)
        inertiaMatrix = [
            Ixx Ixy Ixz
            Ixy Iyy Iyz
            Ixz Iyz Izz
        ]
    else
        @assert size(inertiaMatrix) == (3,3) "inertiaMatrix must be a 3x3 matrix"
    end

    # Convert to correct types
    mass = Float64(mass)
    η = SVector{3,Float64}(η)
    inertiaMatrix = SMatrix{3,3,Float64,9}(inertiaMatrix)

    return PointInertia(
        mass,
        η,
        inertiaMatrix
        )
end
export create_PointInertia