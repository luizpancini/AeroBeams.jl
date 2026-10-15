# CONSTANTS
# -----------------------------------------------------------------------------
# First vector of basis A, resolved in that basis
const a1 = @SVector [1, 0, 0]

# Second vector of basis A, resolved in that basis
const a2 = @SVector [0, 1, 0]

# Third vector of basis A, resolved in that basis
const a3 = @SVector [0, 0, 1]

# Empty 0-element StaticVector
const z0 = @SVector zeros(0)

# 3x1 zero vector
const z3 = @SVector zeros(3)

# 3x3 zero matrix
const Z3 = @SMatrix zeros(3,3)

# 6x6 zero matrix
const Z6 = @SMatrix zeros(6,6)

# 3x3 identity matrix
const I3 = SMatrix{3,3,Float64,9}(LinearAlgebra.I)

# 6x6 identity matrix
const I6 = SMatrix{6,6,Float64,36}(LinearAlgebra.I)

# Empty 3-element StaticVector and 3x3 StaticMatrix
const empty_SVector = SVector{0,Float64}()

const empty_SMatrix = SMatrix{0,0,Float64,0}()

# 3-point Gauss-Legendre quadrature node, points and (halved) weights
const ξ = 0.7745966692414834
const ζ1, ζ2, ζ3 = 0.5 - 0.5ξ, 0.5, 0.5 + 0.5ξ
const w1, w2, w3 = 5/18, 8/18, 5/18

# Gauss–Legendre quadrature nodes and weights for n=2 to n=7
const GL_POINTS = (
    SVector(-0.5773502691896257, 0.5773502691896257), # n=2
    SVector(-0.7745966692414834, 0.0, 0.7745966692414834), # n=3
    SVector(-0.8611363115940526, -0.3399810435848563,
            0.3399810435848563,  0.8611363115940526), # n=4
    SVector(-0.9061798459386640, -0.5384693101056831, 
            0.0,
            0.5384693101056831,  0.9061798459386640), # n=5
    SVector(-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
            0.2386191860831969,  0.6612093864662645,  0.9324695142031521), # n=6
    SVector(-0.9491079123427585, -0.7415311855993945, -0.4058451513773972,
            0.0,
            0.4058451513773972,  0.7415311855993945,  0.9491079123427585) # n=7
)
const GL_WEIGHTS = (
    SVector(1.0, 1.0), # n=2
    SVector(0.5555555555555556, 0.8888888888888888, 0.5555555555555556), # n=3
    SVector(0.3478548451374539, 0.6521451548625461,
            0.6521451548625461, 0.3478548451374539), # n=4
    SVector(0.2369268850561891, 0.4786286704993665, 
            0.5688888888888889,
            0.4786286704993665, 0.2369268850561891), # n=5
    SVector(0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
            0.4679139345726910, 0.3607615730481386, 0.1713244923791704), # n=6
    SVector(0.1294849661688697, 0.2797053914892766, 0.3818300505051189,
            0.4179591836734694,
            0.3818300505051189, 0.2797053914892766, 0.1294849661688697) # n=7
)


## PRIVATE API
# -----------------------------------------------------------------------------
# Gets the ith component of the 3x3 identity matrix
a_i(i::Int) = a[i]

# Divides the input variables in-place
divide_inplace!(divisor::Real, vars...) = 
    foreach(var -> var ./= divisor, vars)

# Multiplies the input variables in-place
multiply_inplace!(multiplier::Real, vars...) =
    foreach(var -> var .*= multiplier, vars)

# Multiplies the input variables in-place and returns them    
multiply_scalars!(multiplier::T, vars::Vararg{Any,N}) where {T<:Real,N} =
    ntuple(i -> vars[i] * multiplier, Val(N))

# Sets all elements of the input arrays to zero
zero_arrays!(arrays...) = foreach(x -> fill!(x, zero(eltype(x))), arrays)

# Computes the scalar product of a third-order tensor represented by matrices A1, A2, and A3 with the vector b
mul3(A1::AbstractMatrix{T},A2::AbstractMatrix{T},A3::AbstractMatrix{T},b::AbstractVector{T}) where {T<:Number} = SMatrix{3,3,T,9}(hcat(A1*b, A2*b, A3*b))

# Copies all fields of a struct
@generated function copy_fields(x::T) where T
    fields = fieldnames(T)
    exprs = [:(copy(getfield(x, $(QuoteNode(f))))) for f in fields]
    :(T($(exprs...)))
end
@generated function copy_fields!(dest::T, src::T) where T
    fields = fieldnames(T)
    exprs = map(fields) do f
        FT = fieldtype(T, f)
        if FT <: Union{MVector,MMatrix}
            :(getfield(dest, $(QuoteNode(f))) .= getfield(src, $(QuoteNode(f))))
        else
            :(copy_fields!(getfield(dest, $(QuoteNode(f))), getfield(src, $(QuoteNode(f)))))
        end
    end
    Expr(:block, exprs..., :dest)
end

# Safe norm function to avoid issues with very small numbers
function safe_norm(p)

    s = sum(abs2, p)
    δ = eps(Float64) * (1 + abs(ForwardDiff.value(real(s)))) * one(real(s))

    return sqrt(s + δ)

end

# Creates a mapping from the sparsity pattern of a sparse matrix J to the indices of a given set of rows and columns
function build_index_map(J, rows, cols)

    N = length(rows) * length(cols)

    # Initialize the index map with the correct size
    map = Vector{Int}(undef, N)

    # Loop over columns and rows to build the index map
    idx = 1
    @inbounds for j in cols
        colrange = J.colptr[j]:(J.colptr[j+1]-1)
        rv = view(J.rowval, colrange)
        for i in rows
            # Check if the index is valid and matches the expected row index 
            k = searchsortedfirst(rv, i)
            (k > length(rv) || rv[k] != i) && error("Sparsity pattern does not contain (i=$i, j=$j)")
            # Assign the index in the map
            map[idx] = colrange[k]
            idx += 1
        end
    end

    return SVector{N, Int}(map)
    
end

# Sets a N = n x m block in the nzval indices of a sparse matrix using the provided index map
@inline function set_block!(nzval, map::SVector{N, Int}, block) where {N}
    @inbounds @simd for i in eachindex(map)
        nzval[map[i]] = block[i]
    end
end

# Adds to a N = n x m block in the nzval indices of a sparse matrix using the provided index map
@inline function add_block!(nzval, map::SVector{N, Int}, block) where {N}
    @inbounds for i in eachindex(map)
        nzval[map[i]] += block[i]
    end
end

# Function barrier for allocation-free evaluation of motion functions
@inline function eval_motion(f::Function, t::Float64)::SVector{3,Float64}
    return f(t)::SVector{3,Float64}
end

# Computes the beam length from the curvature vector and the endpoint position vector
function beam_length(r,k)

    # Norm of the endpoint position vector and its square
    rn = norm(r)
    rn2 = rn^2

    # Norm of the curvature vector
    knorm = norm(k)

    # Set length according to curvature vector
    if iszero(knorm) || iszero(k[2]) && iszero(k[3]) # Uncurved, but potentially twisted beam
        ℓ = rn
    elseif iszero(k[1]) # Untwisted, curved beam
        ℓ = 2/knorm * asin(knorm/2 * rn)
    else # Twisted and curved beam
        f(Δℓ) = rn2 - ( (k[1]*knorm*Δℓ)^2 + 
                        4*(knorm^2-k[1]^2) * sin(knorm*Δℓ/2)^2 
                      ) / knorm^4
        ℓ = Roots.find_zero(f, rn)
    end

    return ℓ

end

# Computes some useful functions of the curvature vector
function curvature_quantities(k)

    # External self-product 
    kkT = k * k'

    # Skew-symmetric matrix 
    ktilde = tilde(k)

    # Norm and its square
    knorm2 = dot(k,k)
    knorm = sqrt(knorm2)

    return kkT, ktilde, knorm, knorm2

end

# Computes the position vector at an arclength value
function position_vector(R0,k,x1)

    # Curvature quantities
    kkT, ktilde, knorm, knorm2 = curvature_quantities(k)

    # Set according to knorm
    r = knorm > 0 ?
        R0 * (1/knorm * (I3 - kkT/knorm2) * sin(knorm*x1) +
              ktilde/knorm2 * (1 - cos(knorm*x1)) + 
              kkT/knorm2 * x1
             ) * a1 :
        x1 * R0[:,1]

    return MVector{3,Float64}(r)

end

# Computes the rotation tensor at an arclength position
function rotation_tensor_from_curvature(R0,k,x1) 

    # Curvature quantities
    kkT, ktilde, knorm, knorm2 = curvature_quantities(k)

    # Set according to knorm
    R = knorm > 0 ?
        R0 * ( (I3 - kkT/knorm2) * cos(knorm*x1) + 
               ktilde/knorm * sin(knorm*x1) + 
               kkT/knorm2
             ) :
        R0

    # Round off
    round_off!(MMatrix{3,3}(R))

    return SMatrix{3,3}(R)

end

# Private version of function to compute the skew-symmetric matrix of a 3-element StaticVector
function tilde(v::StaticVector{3,T}) where {T<:Number}
    z = zero(T)
    return SMatrix{3,3,T,9}(
        z, v[3], -v[2],
        -v[3], z, v[1],
        v[2], -v[1], z
    )
end

# No-overhead, in-place allocation, private version of the Gauss–Legendre quadrature function with 3 Gauss points for a function weighted by the element shape functions (1-ζ) and ζ, where ζ is the local coordinate in the interval [0, 1]
function gauss_legendre3_shaped_fun!(I, f)

    # Gauss points and weights
    ξ = (-0.7745966692414834, 0.0, 0.7745966692414834)
    w = (5/9, 8/9, 5/9)

    # Element local coordinates at Gauss points
    ζ = (ξ .+ 1) ./ 2

    # Compute integrals at nodes (nodal resultants)
    fill!(I, 0.0)
    @inbounds for i in 1:3
        val = f(ζ[i])
        I[1,:] .+= w[i] .* val .* (1 - ζ[i])
        I[2,:] .+= w[i] .* val .* ζ[i]
    end
    I .*= 0.5

    return I

end

# In-place version of the Wiener-Milenkovic rotation tensor function
function R_WM!(R,p)

    Θ,υ²,ps,ps0,pspsT,ps_p = R_core(p)

    R .= υ² * Θ

    return ps,ps0,υ²,pspsT,ps_p

end

# Computes the Wiener-Milenkovic rotation tensor's "core" and associated quantities
function R_core(p)

    # Scaling factor and scaled rotation parameters
    λ,ps,_,pNorm,psNorm = rotation_parameter_scaling(p)

    # Scaled parameters derivatives w.r.t. extended rotation parameters
    ps_p = ∂ps_∂p(λ,p,pNorm)

    # Shorthand for components of scaled rotation parameters 
    ps1 = ps[1]
    ps2 = ps[2]
    ps3 = ps[3]

    # Useful functions of the scaled rotation parameters
    pspsT = SMatrix{3,3}(ps*ps')
    ps0 = 2 - psNorm^2/8
    υ = 1/(4-ps0)
    υ² = υ^2
    ps0s = ps0^2
    ps1s = ps1^2
    ps2s = ps2^2
    ps3s = ps3^2
    ps0ps1 = ps0*ps1
    ps0ps2 = ps0*ps2
    ps0ps3 = ps0*ps3
    ps1ps2 = ps1*ps2
    ps2ps3 = ps2*ps3
    ps1ps3 = ps1*ps3

    # Rotation tensor's "core"
    Θ = @SMatrix [ps0s+ps1s-ps2s-ps3s    2*(ps1ps2-ps0ps3)    2*(ps1ps3+ps0ps2);
                    2*(ps1ps2+ps0ps3)  ps0s-ps1s+ps2s-ps3s    2*(ps2ps3-ps0ps1);
                    2*(ps1ps3-ps0ps2)    2*(ps2ps3+ps0ps1)  ps0s-ps1s-ps2s+ps3s]

    return Θ,υ²,ps,ps0,pspsT,ps_p

end

# Scales the Wiener-Milenkovic rotation parameters
function rotation_parameter_scaling(p::AbstractVector{T}) where {T<:Number}

    # Initialize scaling factor and number of odd half rotations
    λ = one(T)
    halfRotations = 0

    # Norm of the extended parameters vector
    pNorm = safe_norm(p)
    
    # Scale according to norm
    if pNorm > 4
        halfRotations = round(Int,pNorm/8)
        λ -= 8 * halfRotations / pNorm
    end
    
    # Scaled parameters
    ps = λ * SVector{3}(p)
    
    # Norm of scaled parameters vector
    psNorm = abs(λ) * pNorm

    return λ,ps,halfRotations,pNorm,psNorm
    
end

# Computes the derivatives of the scaled rotation parameters with respect to the extended rotation parameters, ∂ps/∂p = λ*I + (1-λ)*p*p'/|p|²
function ∂ps_∂p(λ::T,p,pNorm) where T<:Number

    p = SVector{3}(p)

    return λ == one(T) ? 
           SMatrix{3,3,T,9}(LinearAlgebra.I) : 
           SMatrix{3,3,T,9}(λ*I3 + (1-λ)/pNorm^2 * p*p')

end

# Computes the transpose of the tangent operator tensor w.r.t. to the extended Wiener-Milenkovic rotation parameters
function HT_WM(p::AbstractVector{<:Real})

    # Scaling factor and scaled rotation parameters
    λ,ps,_,pNorm,psNorm = rotation_parameter_scaling(p)

    # Scaled parameters derivatives w.r.t. extended rotation parameters
    ps_p = ∂ps_∂p(λ,p,pNorm)

    # Useful functions of the scaled rotation parameters 
    ps0 = 2 - psNorm^2/8
    υ² = (4-ps0)^-2
    pspsT = SMatrix{3,3}(ps*ps')

    # Tangent operator transpose
    return HT_WM(ps,ps0,υ²,pspsT,ps_p)

end
function HT_WM(ps,ps0,υ²,pspsT,ps_p)

    # ω = HT(p)*pdot = HT(ps)*psdot -> HT(p) = HT(ps)*∂ps/∂p
    return 2*υ²*(ps0*I3 + 1/4*pspsT - tilde(ps)) * ps_p

end
function HT_WM!(HT,ps,ps0,υ²,pspsT,ps_p)

    HT .= HT_WM(ps,ps0,υ²,pspsT,ps_p)

end

# Computes the inverse of the transpose of the tangent operator tensor w.r.t the extended Wiener-Milenkovic parameters
function HTinv_WM(p)

    # Scaling factor and scaled rotation parameters
    λ,ps,_,pNorm,psNorm = rotation_parameter_scaling(p)

    # Scaled parameters derivatives w.r.t. extended rotation parameters
    ps_p = ∂ps_∂p(λ,p,pNorm)

    # Useful functions of the scaled rotation parameters 
    ps0 = 2 - psNorm^2/8
    pspsT = SMatrix{3,3}(ps*ps')

    # Tangent operator inverse transpose
    return HTinv_WM(ps,ps0,pspsT,ps_p)

end
function HTinv_WM(ps,ps0,pspsT,ps_p)

    # pdot = (∂ps/∂p)^-1 * HTinv(ps) * ω -> HTinv(p) = (∂ps/∂p)^-1 * HTinv(ps)
    return inv(ps_p) * 1/2*(ps0*I3 + 1/4*pspsT + tilde(ps))

end
function HTinv_WM!(HTinv,ps,ps0,pspsT,ps_p)

    HTinv .= HTinv_WM(ps,ps0,pspsT,ps_p)

end

# Computes the time derivative of the rotation tensor from basis A to basis B, resolved in basis A ( Rdot = R * tilde(HT * pdot) )
function RdotR0!(RdotR0,R,R0,HT,pdot)

    # Static versions for performance
    R = SMatrix{3,3}(R)
    HT = SMatrix{3,3}(HT)
    pdot = SVector{3}(pdot)
    
    RdotR0 .= R * tilde(HT * pdot) * R0

end

# Computes the derivatives of the rotation tensor with respect to the extended rotation parameters
function ∂R_∂p(p)

    # Jacobian of the vectorized rotation tensor, as an array indexed by (row, column, parameter)
    J = SArray{Tuple{3,3,3}}(Tuple(ForwardDiff.jacobian(R_vectorized, SVector{3}(p))))

    return J[:,:,1], J[:,:,2], J[:,:,3]

end
R_vectorized(p) = vec(rotation_tensor_WM(p))
HTinv_vectorized(p) = vec(HTinv_WM(p))

# Computes the derivatives of the rotation tensor, the tangent tensor's transpose and its inverse with respect to the extended rotation parameters
function ∂R_∂p_and_∂HT_∂p_and_∂HTinv_∂p!(R_p1,R_p2,R_p3,HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3,p)

    # Jacobian of the stacked tensors, as an array indexed by (row, column, tensor, parameter)
    J = SArray{Tuple{3,3,3,3}}(Tuple(ForwardDiff.jacobian(R_HT_HTinv_stacked, SVector{3}(p))))

    # Derivatives of the rotation tensor
    R_p1 .= J[:,:,1,1]
    R_p2 .= J[:,:,1,2]
    R_p3 .= J[:,:,1,3]

    # Derivatives of the tangent tensor's transpose
    HT_p1 .= J[:,:,2,1]
    HT_p2 .= J[:,:,2,2]
    HT_p3 .= J[:,:,2,3]

    # Derivatives of the tangent tensor's transpose inverse
    HTinv_p1 .= J[:,:,3,1]
    HTinv_p2 .= J[:,:,3,2]
    HTinv_p3 .= J[:,:,3,3]

end
R_HT_HTinv_stacked(p) = vcat(vec(rotation_tensor_WM(p)), vec(HT_WM(p)), vec(HTinv_WM(p)))

# Computes the derivative of the time derivative of the rotation tensor with respect to the extended rotation parameters, from Rdot = R*tilde(HT*pdot) with pdot = 2/Δt*p - const
function ∂Rdot_∂p!(Rdot_p1,Rdot_p2,Rdot_p3,R,HT,R_p1,R_p2,R_p3,HT_p1,HT_p2,HT_p3,pdot,Δt)

    # Static versions for performance
    R = SMatrix{3,3}(R)
    HT = SMatrix{3,3}(HT)
    pdot = SVector{3}(pdot)
    R_p1 = SMatrix{3,3}(R_p1)
    R_p2 = SMatrix{3,3}(R_p2)
    R_p3 = SMatrix{3,3}(R_p3)
    HT_p1 = SMatrix{3,3}(HT_p1)
    HT_p2 = SMatrix{3,3}(HT_p2)
    HT_p3 = SMatrix{3,3}(HT_p3)

    # Angular velocity component
    ω = HT * pdot

    # Derivatives w.r.t. each extended rotation parameter
    Rdot_p1 .= R_p1*tilde(ω) + R*tilde(HT_p1*pdot + 2/Δt*HT[:,1])
    Rdot_p2 .= R_p2*tilde(ω) + R*tilde(HT_p2*pdot + 2/Δt*HT[:,2])
    Rdot_p3 .= R_p3*tilde(ω) + R*tilde(HT_p3*pdot + 2/Δt*HT[:,3])

end

## PUBLIC API
# ------------------------------------------------------------------------------
"""
    round_off(x,tol)

Rounds the array to input tolerance

# Arguments
- `x::AbstractArray{<:Number}`: array
- `tol::Real`: tolerance (defaults to machine epsilon)
"""
function round_off(x, tol::Real=eps(eltype(x)))
    return @. round(x/tol)*tol
end
export round_off


"""
    round_off!(x,tol)

Rounds the array to input tolerance in-place

# Arguments
- `x::AbstractArray{<:Number}`: array
- `tol::Real`: tolerance (defaults to machine epsilon)
"""
function round_off!(x::AbstractArray{<:Number}, tol::Real=eps(eltype(x)))
    x .= round_off(x, tol)
    return x
end
export round_off!


"""
    average(x)

Computes the arithmetic mean of an array `x`.

# Arguments
- `x::AbstractArray{<:Real}`: array
"""
function average(x::AbstractArray{T}) where {T<:Real}
    return float(T)(sum(x) / length(x))
end


"""
    average(args...)

Computes the arithmetic mean of a list of numbers.

# Arguments
- `args...`: list of numbers
"""
function average(args...)
    return sum(args)/length(args)
end


"""
    rms(x)

Computes the root mean square of an array `x`

# Arguments
- `x::AbstractArray{<:Real}`: array
"""
function rms(x::AbstractArray{T}) where {T<:Real}
    return float(T)(sqrt(average(x .^ 2)))
end
export rms


"""
    tilde(v)

Computes the skew-symmetric matrix associated with a three-element vector `v`

# Arguments
- `v::AbstractVector{<:Number}`: 3-element vector
"""
function tilde(v::AbstractVector{T}) where {T<:Number}
    @assert length(v) == 3 "Input vector must have length 3"
    z = zero(T)
    return T[  z  -v[3]   v[2];
             v[3]    z   -v[1];
            -v[2]  v[1]      z]
end
export tilde


"""
    axial(R)

Computes the axial part of a 3x3 matrix

# Arguments
- `R::AbstractMatrix{<:Number}`: 3x3 matrix
"""
function axial(R::AbstractMatrix{T}) where {T<:Number}
    @assert size(R) == (3,3) "Input matrix must be 3x3"
    return 1/2 * T[R[3,2]-R[2,3]; 
                   R[1,3]-R[3,1];
                   R[2,1]-R[1,2]
                   ]
end
export axial


"""
    gauss_legendre(f, a, b; n)

Computes the Gauss–Legendre quadrature (integral) for the vector-valued function `f` in the interval `[a, b]`, using `n` Gauss points (`2 ≤ n ≤ 7`)

# Arguments
- `f::Function`: function to integrate
- `a::Real`: lower limit
- `b::Real`: upper limit
- `n::Int`: number of Gauss points
"""
function gauss_legendre(f::Function, a::Real, b::Real; n::Int=3)

    # Gauss points and weights
    ξ = GL_POINTS[n-1]
    w = GL_WEIGHTS[n-1]

    # Change of interval from [-1, 1] to [a, b]
    mid = (a + b) / 2
    half = (b - a) / 2

    # Weighted sum
    s = 0.0
    @inbounds @simd for i in eachindex(w)
        s += w[i] * f(mid + half*ξ[i])
    end

    return half * s
    
end
export gauss_legendre


"""
    isotropic_stiffness_matrix(; kwargs...)

Creates a 6x6 sectional stiffness matrix for a cross-section made of isotropic material

# Keyword Arguments
- `∞::Real`: value for rigid properties
- `EA::Union{Missing,Real}`: axial stiffness
- `GAy::Union{Missing,Real}`: shear stiffness in the x2 direction
- `GAz::Union{Missing,Real}`: shear stiffness in the x3 direction
- `GJ::Union{Missing,Real}`: torsional stiffness
- `EIy::Union{Missing,Real}`: bending stiffness in the x2 direction
- `EIz::Union{Missing,Real}`: bending stiffness in the x3 direction
- `E::Union{Missing,Real}`: Young's modulus, used to compute `EA`, `EIy` and `EIz` when those are not provided explicitly
- `G::Union{Missing,Real}`: shear modulus, used to compute `GAy`, `GAz` and `GJ` when those are not provided explicitly
- `ν::Union{Missing,Real}`: Poisson's ratio, used to compute `G = E / (2(1 + ν))` when `G` is not provided explicitly
- `A::Union{Missing,Real}`: area, used to compute `EA` and shear stiffnesses when those are not provided explicitly
- `Iy::Union{Missing,Real}`: second moment of area about x2, used to compute `EIy` when it is not provided explicitly
- `Iz::Union{Missing,Real}`: second moment of area about x3, used to compute `EIz` when it is not provided explicitly
- `J::Union{Missing,Real}`: torsion constant, used to compute `GJ` when it is not provided explicitly
- `Ksy::Real`: shear correction factor for x2, used to compute `GAy = G*A*Ksy`
- `Ksz::Real`: shear correction factor for x3, used to compute `GAz = G*A*Ksz`
- `Kt::Real`: torsion correction factor, used to compute `GJ = G*J*Kt`
- `t2::Real`: offset from reference line to tension center in the x2 direction
- `t3::Real`: offset from reference line to tension center in the x3 direction
- `s2::Real`: offset from reference line to shear center in the x2 direction
- `s3::Real`: offset from reference line to shear center in the x3 direction
"""
function isotropic_stiffness_matrix(; 
    ∞::Real=1e16,
    EA::Union{Missing,Real}=missing,
    GAy::Union{Missing,Real}=missing,
    GAz::Union{Missing,Real}=missing,
    GJ::Union{Missing,Real}=missing,
    EIy::Union{Missing,Real}=missing,
    EIz::Union{Missing,Real}=missing,
    E::Union{Missing,Real}=missing,
    G::Union{Missing,Real}=missing,
    ν::Union{Missing,Real}=missing,
    A::Union{Missing,Real}=missing,
    Iy::Union{Missing,Real}=missing,
    Iz::Union{Missing,Real}=missing,
    J::Union{Missing,Real}=missing,
    Ksy::Real=1,
    Ksz::Real=1,
    Kt::Real=1,
    t2::Real=0,
    t3::Real=0,
    s2::Real=0,
    s3::Real=0,
    )

    # Resolve elastic properties from the inputs when needed
    E = ismissing(E) ? 
        missing : 
        E
    G = ismissing(G) ? 
            (ismissing(E) || ismissing(ν) ? 
            missing : 
            E / (2 * (1 + ν))) : 
        G

    EA = ismissing(EA) ? 
            (ismissing(E) || ismissing(A) ? 
            missing : 
            E * A) : 
         EA
    GAy = ismissing(GAy) ? 
            (ismissing(G) || ismissing(A) ? 
            missing : 
            G * A * Ksy) : 
          GAy
    GAz = ismissing(GAz) ? 
            (ismissing(G) || ismissing(A) ? 
            missing : 
            G * A * Ksz) : 
          GAz
    GJ = ismissing(GJ) ? 
            (ismissing(G) || ismissing(J) ? 
            missing : 
            G * J * Kt) : 
         GJ         
    EIy = ismissing(EIy) ? 
            (ismissing(E) || ismissing(Iy) ? 
            missing : 
            E * Iy) : 
          EIy
    EIz = ismissing(EIz) ? 
            (ismissing(E) || ismissing(Iz) ? 
            missing : 
            E * Iz) :
          EIz

    # Preserve rigid defaults when no derived values are available
    EA = ismissing(EA) ? 
         ∞ : 
         EA
    GAy = ismissing(GAy) ?
          ∞ : 
          GAy
    GAz = ismissing(GAz) ? 
          ∞ : 
          GAz
    GJ = ismissing(GJ) ? 
         ∞ : 
         GJ
    EIy = ismissing(EIy) ? 
          ∞ : 
          EIy
    EIz = ismissing(EIz) ? 
          ∞ : 
          EIz

    # Validate
    @assert all(x->x>0,[EA,GAy,GAz,GJ,EIy,EIz]) "stiffness properties must be positive"

    # Set sectional stiffness matrix (See Hodges' book Eqs. 4.114 - 4.118)
    z = [0 t3 -t2; -s3 0 0; s2 0 0]

    A = diagm([EA,GAy,GAz])

    Tinv = diagm([GJ,EIy,EIz])

    B = A * z

    D = Tinv + z' * A * z

    S = [A B; B' D]

    return S

end
export isotropic_stiffness_matrix


"""
    anisotropic_stiffness_matrix(; kwargs...)

Creates a 6x6 sectional stiffness matrix for a cross-section made of anisotropic material

# Keyword Arguments
- `∞::Real`: value for rigid properties
- `EA::Real`: axial stiffness
- `GAy::Real`: shear stiffness in the x2 direction
- `GAz::Real`: shear stiffness in the x3 direction
- `GJ::Real`: torsional stiffness
- `EIy::Real`: bending stiffness in the x2 direction
- `EIz::Real`: bending stiffness in the x3 direction
- `c12::Real`: coupling stiffness between axial and shear in the x2 direction
- `c13::Real`: coupling stiffness between axial and shear in the x3 direction
- `c14::Real`: coupling stiffness between axial and torsion
- `c15::Real`: coupling stiffness between axial and bending in the x2 direction
- `c16::Real`: coupling stiffness between axial and bending in the x3 direction
- `c23::Real`: coupling stiffness between shear in the x2 and x3 directions
- `c24::Real`: coupling stiffness between shear in the x2 direction and torsion
- `c25::Real`: coupling stiffness between shear in the x2 direction and bending in the x2 direction
- `c26::Real`: coupling stiffness between shear in the x2 direction and bending in the x3 direction
- `c34::Real`: coupling stiffness between shear in the x3 direction and torsion
- `c35::Real`: coupling stiffness between shear in the x3 direction and bending in the x2 direction
- `c36::Real`: coupling stiffness between shear in the x3 direction and bending in the x3 direction
- `c45::Real`: coupling stiffness between torsion and bending in the x2 direction
- `c46::Real`: coupling stiffness between torsion and bending in the x3 direction
- `c56::Real`: coupling stiffness between bending in the x2 and x3 directions
"""
function anisotropic_stiffness_matrix(; 
    ∞::Real=1e16,
    EA::Real=∞,
    GAy::Real=∞,
    GAz::Real=∞,
    GJ::Real=∞,
    EIy::Real=∞,
    EIz::Real=∞,
    c12::Real=0,
    c13::Real=0,
    c14::Real=0,
    c15::Real=0,
    c16::Real=0,
    c23::Real=0,
    c24::Real=0,
    c25::Real=0,
    c26::Real=0,
    c34::Real=0,
    c35::Real=0,
    c36::Real=0,
    c45::Real=0,
    c46::Real=0,
    c56::Real=0,
    )

    # Validate
    @assert all(x->x>0,[EA,GAy,GAz,GJ,EIy,EIz]) "diagonal stiffness properties must be positive"

    # Set matrix
    S = [ EA   c12  c13  c14  c15  c16;
         c12   GAy  c23  c24  c25  c26;
         c13   c23  GAz  c34  c35  c36;
         c14   c24  c34   GJ  c45  c46;
         c15   c25  c35  c45  EIy  c56;
         c16   c26  c36  c46  c56  EIz]

    return S

end
export anisotropic_stiffness_matrix


"""
    inertia_matrix(; kwargs...)

Creates a 6x6 sectional inertia matrix

# Keyword Arguments
- `ρA::Union{Missing,Real}`: mass per unit length
- `ρIy::Union{Missing,Real}`: mass moment of inertia per unit length about the x2-axis
- `ρIz::Union{Missing,Real}`: mass moment of inertia per unit length about the x3-axis
- `ρIyz::Union{Missing,Real}`: mass product of inertia per unit length
- `ρ::Union{Missing,Real}`: material density, used to compute `ρA`, `ρIy`, `ρIz` and `ρIyz` when those are not provided explicitly
- `A::Union{Missing,Real}`: area, used to compute `ρA` when it is not provided explicitly
- `Iy::Union{Missing,Real}`: second moment of area about x2, used to compute `ρIy` when it is not provided explicitly
- `Iz::Union{Missing,Real}`: second moment of area about x3, used to compute `ρIz` when it is not provided explicitly
- `Iyz::Union{Missing,Real}`: product of area inertia, used to compute `ρIyz` when it is not provided explicitly
- `e2::Real`: offset from reference line to center of gravity in the x2 direction
- `e3::Real`: offset from reference line to center of gravity in the x3 direction
"""
function inertia_matrix(; 
    ρA::Union{Missing,Real}=missing,
    ρIy::Union{Missing,Real}=missing,
    ρIz::Union{Missing,Real}=missing,
    ρIyz::Union{Missing,Real}=missing,
    ρ::Union{Missing,Real}=missing,
    A::Union{Missing,Real}=missing,
    Iy::Union{Missing,Real}=missing,
    Iz::Union{Missing,Real}=missing,
    Iyz::Union{Missing,Real}=missing,
    e2::Real=0,
    e3::Real=0
    )

    # Resolve sectional inertia properties from material/geometry inputs when needed
    ρA = ismissing(ρA) ?
            (ismissing(ρ) || ismissing(A) ?
            missing :
            ρ * A) :
         ρA
    ρIy = ismissing(ρIy) ?
            (ismissing(ρ) || ismissing(Iy) ?
            missing :
            ρ * Iy) :
          ρIy
    ρIz = ismissing(ρIz) ?
            (ismissing(ρ) || ismissing(Iz) ?
            missing :
            ρ * Iz) :
          ρIz
    ρIyz = ismissing(ρIyz) ?
            (ismissing(ρ) || ismissing(Iyz) ?
            missing :
            ρ * Iyz) :
           ρIyz

    # Preserve previous defaults when no values are available
    ρA = ismissing(ρA) ? 
         0 : 
         ρA
    ρIy = ismissing(ρIy) ? 
          0 : 
          ρIy
    ρIz = ismissing(ρIz) ? 
          0 : 
          ρIz
    ρIyz = ismissing(ρIyz) ? 
           0 : 
           ρIyz

    # Validate
    @assert all(x->x>=0,[ρA,ρIy,ρIz]) "mass per unit length and mass moments of inertia must be non-negative"

    # Skew-symmetric operator of the offset vector from reference line to center of mass
    ηtilde = tilde([0; ρA*e2; ρA*e3])

    # Sectional inertia matrix
    I = diagm([ρA, ρA, ρA, ρIy+ρIz, ρIy, ρIz])

    I[5,6] = I[6,5] = ρIyz
    I[1:3,4:6] = -ηtilde
    I[4:6,1:3] = ηtilde

    return I

end
export inertia_matrix


"""
    rotation_tensor_E321(p)

Computes the rotation tensor according to Euler parameters sequence 3-2-1

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters
"""
function rotation_tensor_E321(p::AbstractVector{<:Number})

    # Euler angles 3-2-1 sequence (yaw, roll and pitch angles), respective sines and cosines
    yaw,roll,pitch = p
    sy,cy = sincos(yaw)
    sr,cr = sincos(roll)
    sp,cp = sincos(pitch)

    # Rotation tensor that brings the reference basis to the final basis
    R = [cr*cy cy*sp*sr-cp*sy sp*sy+cp*cy*sr;
         cr*sy cy*cp+sp*sr*sy cp*sr*sy-cy*sp;
           -sr          cr*sp          cr*cp]

    round_off!(R)       

    return R

end
export rotation_tensor_E321


"""
    rotation_tensor_E213(p)

Computes the rotation tensor according to Euler parameters sequence 2-1-3.

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters
"""
function rotation_tensor_E213(p::AbstractVector{<:Number})

    # Euler angles 2-1-3 sequence (roll, pitch, yaw), respective sines and cosines
    roll,pitch,yaw = p
    sr,cr = sincos(roll)
    sp,cp = sincos(pitch)
    sy,cy = sincos(yaw)

    # Rotation tensor
    R = [cy*cr+sy*sp*sr cy*sp*sr-sy*cr cp*sr;
                  sy*cp          cy*cp   -sp;
         sy*sp*cr-cy*sr sy*sr+cy*sp*cr cp*cr]

    round_off!(R)

    return R

end
export rotation_tensor_E213


"""
    rotation_tensor_E231(p)

Computes the rotation tensor according to Euler parameters sequence 2-3-1.

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters
"""
function rotation_tensor_E231(p::AbstractVector{<:Number})

    # Euler angles 2-3-1 sequence (roll, yaw, pitch), respective sines and cosines
    roll,yaw,pitch = p
    sr,cr = sincos(roll)
    sy,cy = sincos(yaw)
    sp,cp = sincos(pitch)

    # Rotation tensor
    R = [ cr*cy sp*sr-sy*cp*cr sr*cp+cr*sp*sy;
             sy          cp*cy         -sp*cy;
         -sr*cy cr*sp+cp*sy*sr cp*cr-sp*sr*sy]

    round_off!(R)

    return R

end
export rotation_tensor_E231


"""
    rotation_tensor_E313(p)

Computes the rotation tensor according to Euler parameters sequence 3-1-3

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters
"""
function rotation_tensor_E313(p::AbstractVector{<:Number})

    # Euler angles 3-1-3 sequence (precession=ϕ, nutation=θ and spin=ψ angles), respective sines and cosines
    ϕ,θ,ψ = p
    sϕ,cϕ = sincos(ϕ)
    sθ,cθ = sincos(θ)
    sψ,cψ = sincos(ψ)
       
    # Rotation tensor that brings the reference basis to the final basis
    R = [-sϕ*cθ*sψ+cϕ*cψ -sϕ*cθ*cψ-cϕ*sψ  sϕ*sθ;
          cϕ*cθ*sψ+sϕ*cψ  cϕ*cθ*cψ-sϕ*sψ -cϕ*sθ;
                   sθ*sψ           sθ*cψ     cθ]

    round_off!(R)

    return R

end
export rotation_tensor_E313


"""
    rotation_tensor_quaternion(q)

Computes the rotation tensor according to quaternion parameters

# Arguments
- `q::AbstractVector{<:Number}`: quaternion parameters
"""
function rotation_tensor_quaternion(q::AbstractVector{<:Number})

    # Validate inputs
    @assert length(q) == 4 "quaternion must be a 4-element vector"

    # Quaternion components
    q0,q1,q2,q3 = q

    # Rotation tensor that brings the reference basis to the final basis
    R = [q0^2+q1^2-q2^2-q3^2    2*(q1*q2-q0*q3)      2*(q1*q3+q0*q2);
           2*(q1*q2+q0*q3)    q0^2-q1^2+q2^2-q3^2    2*(q2*q3-q0*q1);
           2*(q1*q3-q0*q2)      2*(q2*q3+q0*q1)    q0^2-q1^2-q2^2+q3^2]

    round_off!(R)

    return R

end
export rotation_tensor_quaternion


""" 
    transform_to_WM(p; from)

Transforms rotation parameters from a given parametrization to Wiener-Milenkovic parameters

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters

# Keyword arguments
- `from::Symbol`: rotation parametrization, either :E321 (Euler angles 3-2-1), :E213 (Euler angles 2-1-3), :E231 (Euler angles 2-3-1), :E313 (Euler angles 3-1-3) or :quaternion
"""
function transform_to_WM(p::AbstractVector{<:Number}; from::Symbol)

    # Validate inputs
    @assert length(p) == 3 "Input rotation parameters must be a 3-element vector"
    @assert from in [:E321, :E213, :E231, :E313, :quaternion] "from must be one of :E321, :E213, :E231, :E313, :quaternion"

    # Transform to Wiener-Milenkovic parameters
    p_WM = if from == :E321
        rotation_parameters_WM(rotation_tensor_E321(p))
    elseif from == :E213
        rotation_parameters_WM(rotation_tensor_E213(p))
    elseif from == :E231
        rotation_parameters_WM(rotation_tensor_E231(p))
    elseif from == :E313
        rotation_parameters_WM(rotation_tensor_E313(p))
    elseif from == :quaternion
        rotation_parameters_WM(rotation_tensor_quaternion(p))
    end

    return p_WM

end
export transform_to_WM


"""
    rotation_tensor_WM(p)

Computes the rotation tensor according to Wiener-Milenkovic parameters

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters
"""
function rotation_tensor_WM(p::AbstractVector{<:Number})

    # Core
    Θ,υ²,_ = R_core(p)
            
    # Rotation tensor
    R = υ² * Θ

    return R

end
export rotation_tensor_WM


"""
    scaled_rotation_parameters(p)

Returns the scaled (extended) Wiener-Milenkovic rotation parameters

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters
"""
function scaled_rotation_parameters(p::AbstractVector{<:Number})

    # Perform scaling
    _,ps,_ = rotation_parameter_scaling(p)

    return ps
    
end
export scaled_rotation_parameters


"""
    twist_angle(p; rotationParametrization, outputUnit)

Computes the twist angle given the rotation parameters

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters

# Keyword arguments
- `rotationParametrization::Symbol`: rotation parametrization, either :WM (Wiener-Milenkovic), :E321 (Euler angles 3-2-1), :E213 (Euler angles 2-1-3) or :E313 (Euler angles 3-1-3)
- `outputUnit::Symbol`: output unit for the twist angle, either :deg (degrees) or :rad (radians)
"""
function twist_angle(p::AbstractVector{<:Number}; rotationParametrization::Symbol=:WM, outputUnit::Symbol=:deg)

    # Validate inputs
    @assert length(p) == 3 "Input rotation parameters must be a 3-element vector"
    @assert rotationParametrization in [:WM, :E321, :E213, :E313] "rotationParametrization must be one of :WM, :E321, :E213, :E313"
    @assert outputUnit in [:deg, :rad] "outputUnit must be either :deg or :rad"

    # Rotation tensor
    R = if rotationParametrization == :WM
        rotation_tensor_WM(p)
    elseif rotationParametrization == :E321
        rotation_tensor_E321(p)
    elseif rotationParametrization == :E213
        rotation_tensor_E213(p)
    elseif rotationParametrization == :E313
        rotation_tensor_E313(p)
    end

    # Rotation axis
    Δ = R*a2

    # Twist angle
    twist = outputUnit == :deg ? 
            asind(Δ[3]) : 
            asin(Δ[3])

    return twist

end
export twist_angle


"""
    rotation_angle(p)

Computes the rotation angle given the Wiener-Milenkovic rotation parameters

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters
"""
function rotation_angle(p::AbstractVector{<:Number})

    # Scale rotation parameters and get number of half rotations
    _,_,halfRotations,pNorm,_ = rotation_parameter_scaling(p)

    # Highest absolute value component
    greatestComp = argmax(SVector{3}(abs.(p)))

    # Sign of rotation angle
    signal = sign(p[greatestComp])

    return signal * (4*atan((pNorm-8*halfRotations)/4) + 2π*halfRotations)

end
export rotation_angle


"""
    rotation_angle_limited(p)

Computes the rotation angle (in the range -360 to 360 degrees) given the Wiener-Milenkovic rotation parameters

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters
"""
function rotation_angle_limited(p::AbstractVector{<:Number})

    # Scale rotation parameters and get number of half rotations
    _,_,halfRotations,pNorm,_ = rotation_parameter_scaling(p)

    # Highest absolute value component
    greatestComp = argmax(abs.(p))

    # Sign of rotation angle
    signal = sign(p[greatestComp])

    return signal * (4*atan((pNorm-8*halfRotations)/4))

end
export rotation_angle_limited


"""
    rotation_parameters_WM(R)

Computes the Wiener-Milenkovic rotation parameters given a rotation tensor

# Arguments
- `R::AbstractMatrix{<:Number}`: rotation tensor
"""
function rotation_parameters_WM(R::AbstractMatrix{<:Number})

    # Get quaternion
    q = quaternion_from_rotation_tensor(R)

    # Vector part of quaternion
    e = view(q, 2:4)

    # Angle of rotation
    ϕ = 2*asin(safe_norm(e))
    
    # Bauchau's ν parameter for Wiener-Milenkovic parametrization
    ν = cos(ϕ/4)^2

    # Wiener-Milenkovic rotation parameters
    p = 2/ν * SVector{3}(e)

    return p

end
export rotation_parameters_WM


"""
    rotation_parameters_Rodrigues(R)

Computes the Rodrigues rotation parameters given a rotation tensor

# Arguments
- `R::AbstractMatrix{<:Number}`: rotation tensor
"""
function rotation_parameters_Rodrigues(R::AbstractMatrix{<:Number})

    # Get quaternion
    q = quaternion_from_rotation_tensor(R)

    # Vector part of quaternion
    e = view(q, 2:4)

    # Angle of rotation
    ϕ = 2*asin(safe_norm(e))
    
    # Bauchau's ν parameter for Rodrigues parametrization
    ν = cos(ϕ/2)

    # Rodrigues rotation parameters
    p = 2/ν * SVector{3}(e)
    
    return p

end
export rotation_parameters_Rodrigues


"""
    yrp_from_rotation_tensor(R; ϵround=1e-10,ϵsingularity=1e-3,assumeNullYawInSingularity=true)

Computes the Euler angles from the sequence 3-2-1 (yaw, roll, pitch) given the rotation tensor

# Arguments
- `R::AbstractMatrix{<:Number}`: rotation tensor

# Keyword arguments
- `ϵround::Real`: tolerance for rounding off elements of R to zero
- `ϵsingularity::Real`: tolerance to consider singular case
- `assumeNullYawInSingularity::Bool`: flag to assume zero yaw angle in singularity
"""
function yrp_from_rotation_tensor(R::AbstractMatrix{<:Number}; 
    ϵround::Real=1e-10,
    ϵsingularity::Real=1e-3,
    assumeNullYawInSingularity::Bool=true
    )

    # Round-off
    round_off!(R, ϵround)
    
    # Roll angle
    roll = asin(-R[3,1])

    # Check if near singularity (roll ≈ ±π/2)
    if abs(cos(roll)) < ϵsingularity
        # Singular case - only sum of yaw and roll can be determined: set angles according to assumption of either one being zero
        yaw_plus_pitch = atan(R[1,2]/R[1,3])
        if assumeNullYawInSingularity
            yaw = 0
            pitch = yaw_plus_pitch
        else
            pitch = 0
            yaw = yaw_plus_pitch
        end
    else
        # Non-singular case
        yaw = atan(R[2,1]/R[1,1])
        pitch = atan(R[3,2]/R[3,3])
    end

    return yaw, roll, pitch

end
export yrp_from_rotation_tensor


"""
    quaternion_from_rotation_tensor(R)

Computes the quaternion (Euler parameters) given a rotation tensor. Derived from Bauchau's book section 13.3.4. 

# Arguments
- `R::AbstractMatrix{<:Number}`: rotation tensor
"""
function quaternion_from_rotation_tensor(R::AbstractMatrix{<:Number})
    
    # Quaternion component product matrix (T = [1+tr(R) 2*axial(R)'; 2*axial(R) (1-tr(R))*I3+2*Symmetric(R)])
    T = @SMatrix [1+R[1,1]+R[2,2]+R[3,3] R[3,2]-R[2,3] R[1,3]-R[3,1] R[2,1]-R[1,2];
    R[3,2]-R[2,3] 1+R[1,1]-R[2,2]-R[3,3] R[1,2]+R[2,1] R[1,3]+R[3,1];
    R[1,3]-R[3,1] R[2,1]+R[1,2] 1-R[1,1]+R[2,2]-R[3,3] R[2,3]+R[3,2];
    R[2,1]-R[1,2] R[1,3]+R[3,1] R[2,3]+R[3,2] 1-R[1,1]-R[2,2]+R[3,3];
    ]

    # Maximum value along diagonal and position
    m = argmax(diag(T))
    Tmm = T[m,m]

    # Largest quaternion component
    qm = 1/2*sqrt(Tmm)

    # Other components
    q = if m == 1
        SVector(qm, T[1,2]/(4qm), T[1,3]/(4qm), T[1,4]/(4qm))
    elseif m == 2
        SVector(T[2,1]/(4qm), qm, T[2,3]/(4qm), T[2,4]/(4qm))
    elseif m == 3
        SVector(T[3,1]/(4qm), T[3,2]/(4qm), qm, T[3,4]/(4qm))
    else
        SVector(T[4,1]/(4qm), T[4,2]/(4qm), T[4,3]/(4qm), qm)
    end

    # Ensure the scalar part of the quaternion is non-negative
    return q[1] < zero(q[1]) ? -q : q
    
end
export quaternion_from_rotation_tensor


"""
    WM_to_yrp(p)

Transforms Wiener-Milenkovic parameters to Euler parameters of sequence 3-2-1 (yaw, roll, pitch)

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters
"""
function WM_to_yrp(p::AbstractVector{<:Number})

    R,_ = rotation_tensor_WM(p)

    yaw,roll,pitch = yrp_from_rotation_tensor(R)
    
    return [yaw; roll; pitch]

end
export WM_to_yrp


"""
    yrp_to_WM(p)

Transforms Euler parameters of sequence 3-2-1 (yaw, roll, pitch) to Wiener-Milenkovic parameters

# Arguments
- `p::AbstractVector{<:Number}`: rotation parameters
"""
function yrp_to_WM(p::AbstractVector{<:Number})
    return rotation_parameters_WM(rotation_tensor_E321(p))
end
export yrp_to_WM


"""
    rotation_between_WM(p1,p2)

Computes the Wiener-Milenkovic parameters describing the rotation from `p1` to `p2`, i.e., `p12` such that `R2(p2) = R12(p12)*R1(p1)`

# Arguments
- `p1::AbstractVector{<:Number}`: initial rotation parameters
- `p2::AbstractVector{<:Number}`: final rotation parameters
"""
function rotation_between_WM(p1::AbstractVector{<:Number},p2::AbstractVector{<:Number})

    # Rotation tensors
    R1 = rotation_tensor_WM(p1)
    R2 = rotation_tensor_WM(p2)
    R12 = R2*R1'

    # Rotation parameters vector and its negative
    p12 = rotation_parameters_WM(R12)
    mp12 = - SVector{3}(p12)

    # Check consistency: true rotation parameters might actually be -p12 (since rotations of ϕ about axis n or -ϕ about axis -n are equal)
    R12Check1 = rotation_tensor_WM(p12)
    R12Check2 = rotation_tensor_WM(mp12)
    check1 = norm(R12Check1-R12)
    check2 = norm(R12Check2-R12)

    # Select correct rotation parameters vector (the one for which norm(R12-R12_check) = 0)
    if check2 < check1
        p12 .= mp12
    end

    return p12

end
export rotation_between_WM


"""
    mode_tracking(sols; method,α,β,γ,maxCost)

Applies mode tracking based on eigenvectors, frequencies and dampings match.

# Arguments
- `sols::AbstractVector{<:EigenSolution}`: vector of eigensolutions to be tracked

# Keyword Arguments
- `method::Symbol`: mode assignment method
- `α::Real`: weight for eigenvector similarity in hungarian method
- `β::Real`: weight for frequency continuity in hungarian method
- `γ::Real`: weight for damping continuity in hungarian method
- `maxCost::Real`: maximum allowable cost for mode matching in the hungarian method
"""
function mode_tracking(sols::AbstractVector{<:EigenSolution}; method::Symbol=:hungarian,α::Real=1.0,β::Real=0.01,γ::Real=0.01,maxCost::Real=Inf)

    # Validate inputs
    @assert method in (:hungarian, :greedy) "Select mode tracking `method` as either :hungarian or :greedy"
    @assert α > 0 && β >= 0 && γ >= 0 "Weights must satisfy α > 0, β ≥ 0, γ ≥ 0"
    N = length(sols)
    @assert N > 0 "Empty vector of eigensolutions"

    # trackIdx[k][c] = index (in sols[c]) of the mode belonging to track k, or 0 if absent
    trackIdx = Vector{Vector{Int}}()

    # Reference (last seen) eigenvalue and eigenvector of each track
    refλ = ComplexF64[]
    refV = nothing

    # Loop over each solution in the vector of eigensolutions
    for c in 1:N

        # Current solution's eigenvalues and eigenvectors
        λc = sols[c].λ
        Vc = sols[c].V
        nc = length(λc)

        # Skip undefined solutions
        if nc == 0 || any(isnan, λc) || any(isnan, Vc)
            continue
        end

        # Assignment of tracks to modes of the current solution (0 = unassigned)
        nTracks = length(trackIdx)
        assignment = zeros(Int, nTracks)
        if nTracks > 0
            # Compute the cost matrix for matching current modes to existing tracks
            C = eigenmatch_cost_matrix(refV, Vc, refλ, λc; α=α, β=β, γ=γ)
            # Perform mode assignment using the selected method
            assignment = method == :greedy ?
                         greedy_assignment(C) :
                         first(Hungarian.hungarian(C))
            # Reject matches that are too expensive
            for k in 1:nTracks
                if assignment[k] != 0 && C[k,assignment[k]] > maxCost
                    assignment[k] = 0
                end
            end
        end

        # Update existing tracks
        for k in 1:nTracks
            j = assignment[k]
            if j != 0
                trackIdx[k][c] = j
                refλ[k] = λc[j]
                refV[:,k] = Vc[:,j]
            end
        end

        # Start new tracks for unmatched modes
        used = falses(nc)
        for j in assignment
            j != 0 && (used[j] = true)
        end
        for j in findall(!, used)
            idx = zeros(Int, N)
            idx[c] = j
            push!(trackIdx, idx)
            push!(refλ, λc[j])
            refV = isnothing(refV) ? Vc[:,j:j] : hcat(refV, Vc[:,j])
        end

    end

    # Tracked eigenvalues, corresponding frequencies and dampings, and mode shapes
    λT, modeShapesT = tracked_eigenvalues_and_modeshapes(sols, trackIdx)
    freqsT = imag.(λT)
    dampingsT = real.(λT)

    return freqsT, dampingsT, modeShapesT, trackIdx

end
export mode_tracking

# Greedy assignment for rectangular cost matrix (lowest cost first). Returns assignment[row] = col or 0
function greedy_assignment(C)
    nr, nc = size(C)
    assignment = zeros(Int, nr)
    rowFree = trues(nr)
    colFree = trues(nc)
    for lin in sortperm(vec(C))
        i, j = Tuple(CartesianIndices(C)[lin])
        if rowFree[i] && colFree[j]
            assignment[i] = j
            rowFree[i] = false
            colFree[j] = false
        end
    end
    return assignment
end

# Computes the cost matrix (nTracks × nModes)
function eigenmatch_cost_matrix(Vref, Vcurr, λref, λcurr; α, β, γ)
    mac = abs.(Vref' * Vcurr)
    C = zeros(length(λref), length(λcurr))
    for j in eachindex(λcurr), i in eachindex(λref)
        Δλ = λref[i] - λcurr[j]
        C[i,j] = -α * mac[i,j] + β * abs(imag(Δλ)) + γ * abs(real(Δλ))
    end
    return C
end

# Extracts the tracked eigenvalues and corresponding mode shapes based on the tracking indices
function tracked_eigenvalues_and_modeshapes(sols, trackIdx)
    nT, N = length(trackIdx), length(sols)
    λT = fill(complex(NaN, NaN), nT, N)
    modeShapesT = Matrix{Union{ModeShape,Nothing}}(nothing, nT, N)
    for k in 1:nT, c in 1:N
        j = trackIdx[k][c]
        j != 0 && (λT[k,c] = sols[c].λ[j])
        j != 0 && (modeShapesT[k,c] = sols[c].modeShapes[j])
    end
    return λT, modeShapesT
end


"""
    filter_tracked_modes(freqs, damps, modeShapes; fmin, fmax, maxTracks, whichTracks, trackwise, dropEmpty, dropOutsideBounds, tolFreq, tolDamp)

Filters tracked modes based on frequency range and removes empty tracks.

# Arguments
- `freqs::AbstractVecOrMat{<:Real}`: matrix of frequencies (nTracks x N)
- `damps::AbstractVecOrMat{<:Real}`: matrix of damping ratios (nTracks x N)
- `modeShapes::AbstractVecOrMat{Union{ModeShape,Nothing}}`: matrix of mode shapes (nTracks x N)

# Keyword Arguments
- `fmin::Real`: minimum frequency to keep
- `fmax::Real`: maximum frequency to keep
- `maxTracks::Int`: maximum number of tracks to keep
- `whichTracks::Symbol`: whether to keep the first or last `maxTracks` tracks (`:first` or `:last`)
- `trackwise::Bool`: whether to keep a track entirely if it lies within the frequency limits at any point
- `dropEmpty::Bool`: whether to drop tracks that are entirely NaN
- `dropOutsideBounds::Bool`: whether to drop tracks that have frequencies outside the [fmin, fmax] range
- `tolFreq::Real`: tolerance for considering frequencies as zero
- `tolDamp::Real`: tolerance for considering dampings as zero
"""
function filter_tracked_modes(freqs::AbstractVecOrMat{<:Real}, damps::AbstractVecOrMat{<:Real}=zeros(size(freqs)), modeShapes::AbstractVecOrMat{Union{ModeShape,Nothing}}=Array{Union{Nothing,AeroBeams.ModeShape}}(nothing,size(freqs)); fmin::Real=0, fmax::Real=Inf, maxTracks::Int=typemax(Int), whichTracks::Symbol=:first, trackwise::Bool=false, dropEmpty::Bool=true, dropOutsideBounds::Bool=false, tolFreq::Real=eps(), tolDamp::Real=eps())

    # Validate
    @assert size(freqs) == size(damps) == size(modeShapes) "freqs, damps, and modeShapes must have the same dimensions"
    @assert 0 <= fmin <= fmax "0 <= fmin <= fmax must hold"
    @assert maxTracks > 0 "maxTracks must be positive"
    @assert whichTracks == :first || whichTracks == :last "whichTracks must be either :first or :last"

    # Apply tolerance filtering to frequencies and dampings
    f = ifelse.(abs.(freqs) .<= tolFreq, 0.0, float.(freqs))
    d = ifelse.(abs.(damps) .<= tolDamp, 0.0, float.(damps))

    # Keep non-negative frequencies within [fmin, fmax] (comparisons with NaN are false, so absent points are dropped)
    keep = @. fmin <= f <= fmax

    # Optionally keep a track entirely if it lies within the limits at any point
    if trackwise
        anyIn = vec(any(keep, dims=2))
        keep = @. anyIn && !isnan(f) && f >= 0
    end

    f = ifelse.(keep, f, NaN)
    d = ifelse.(keep, d, NaN)

    # Filter mode shapes based on the keep matrix, setting entries to `nothing` if they are not kept
    m = [ifelse(keep[k,c], modeShapes[k,c], nothing) for k in axes(modeShapes,1), c in axes(modeShapes,2)]

    # Indices of the original tracks that remain (dropEmpty drops tracks that are entirely NaN and dropOutsideBounds drops tracks that have any NaN values, i.e., are not entirely within the frequency bounds)
    keptTracks = dropEmpty ?
                 findall(k -> any(!isnan, @view f[k,:]), 1:size(f,1)) :
                 collect(1:size(f,1))              
    keptTracks = dropOutsideBounds ?
                 findall(k -> all(!isnan, @view f[k,:]), 1:size(f,1)) : 
                 keptTracks

    # Limit the number of tracks to maxTracks
    if length(keptTracks) > maxTracks
        keptTracks = whichTracks == :first ? 
                     keptTracks[1:maxTracks] : 
                     keptTracks[end-maxTracks+1:end]
    end
                 
    return f[keptTracks,:], d[keptTracks,:], m[keptTracks,:], keptTracks

end
export filter_tracked_modes


"""
    find_flutter(controlParam,freqs,damps,extras...)

Finds the flutter onset and offset points for each mode given the frequency and damping matrices and the control parameter vector.

# Arguments
- `controlParam::AbstractVector{<:Real}`: vector of control parameter values (length N)
- `freqs::AbstractMatrix{<:Real}`: matrix of frequencies (nTracks x N)
- `damps::AbstractMatrix{<:Real}`: matrix of damping ratios (nTracks x N)
- `extras::AbstractVector{<:Real}...`: additional variables of the same length as `controlParam`

# Keyword Arguments
- `sortby::Symbol=:CP`: determines the sorting criterion for the flutter points
"""
function find_flutter(controlParam::AbstractVector{<:Real},freqs::AbstractMatrix{<:Real},damps::AbstractMatrix{<:Real},extras::AbstractVector{<:Real}...; sortby::Symbol=:CP)

    # Validate
    @assert size(freqs,2) == size(damps,2) == length(controlParam) "Number of columns in freqs and damps must match length of controlParam"
    @assert size(freqs,1) == size(damps,1) "Number of rows in freqs must match number of rows in damps"
    @assert all(length(e) == length(controlParam) for e in extras) "All extra variables must have the same length as controlParam"
    @assert sortby in (:CP, :freq, :mode) "sortby must be either :CP, :freq or :mode"

    # Number of modes and control parameter points
    nModes, N = size(freqs)

    # Initialize outputs
    flutterOnsetCP = Float64[]
    flutterOnsetFreqs = Float64[]
    flutterOnsetModes = Int[]
    flutterOnsetExtras = [Float64[] for _ in extras]
    flutterOffsetCP = Float64[]
    flutterOffsetFreqs = Float64[]
    flutterOffsetModes = Int[]
    flutterOffsetExtras = [Float64[] for _ in extras]

    # Loop over modes
    for mode in 1:nModes
        # Flutter onset
        # ----------------------------------------------------------------------
        # Instability onset indices
        iOnset = 1 .+ findall(i -> damps[mode,i] < 0 && damps[mode,i+1] > 0 && !iszero(freqs[mode,i]) && !iszero(freqs[mode,i+1]), 1:N-1)
        # Skip if onset is at the first index or there are no valid onsets
        (isempty(iOnset) || isempty(filter(x->x!=1,iOnset))) && continue
        # Loop over onset indices
        for i in iOnset
            push!(flutterOnsetCP,LinearInterpolations.interpolate(damps[mode,i-1:i],controlParam[i-1:i],0))
            push!(flutterOnsetFreqs,LinearInterpolations.interpolate(damps[mode,i-1:i],freqs[mode,i-1:i],0))
            push!(flutterOnsetModes,mode)
            for (out,e) in zip(flutterOnsetExtras,extras)
                push!(out,LinearInterpolations.interpolate(damps[mode,i-1:i],e[i-1:i],0))
            end
        end
        # Flutter offset
        # ----------------------------------------------------------------------
        # Instability offset indices
        iOffset = 1 .+ findall(i -> damps[mode,i] > 0 && damps[mode,i+1] < 0 && !iszero(freqs[mode,i]) && !iszero(freqs[mode,i+1]), 1:N-1)
        # Skip if there are no valid offsets
        isempty(iOffset) && continue
        # Loop over offset indices
        for i in iOffset
            push!(flutterOffsetCP,LinearInterpolations.interpolate(-damps[mode,i-1:i],controlParam[i-1:i],0))
            push!(flutterOffsetFreqs,LinearInterpolations.interpolate(-damps[mode,i-1:i],freqs[mode,i-1:i],0))
            push!(flutterOffsetModes,mode)
            for (out,e) in zip(flutterOffsetExtras,extras)
                push!(out,LinearInterpolations.interpolate(-damps[mode,i-1:i],e[i-1:i],0))
            end
        end
    end

    # Sort in ascending order of control parameter values
    if !isempty(flutterOnsetCP)
        sortIdx = if sortby == :CP
            sortperm(flutterOnsetCP)
        elseif sortby == :freq
            sortperm(flutterOnsetFreqs)
        elseif sortby == :mode
            sortperm(flutterOnsetModes)
        end
        flutterOnsetCP = flutterOnsetCP[sortIdx]
        flutterOnsetFreqs = flutterOnsetFreqs[sortIdx]
        flutterOnsetModes = flutterOnsetModes[sortIdx]
        flutterOnsetExtras = map(v -> v[sortIdx], flutterOnsetExtras)
    end
    if !isempty(flutterOffsetCP)
        sortIdx = if sortby == :CP
            sortperm(flutterOffsetCP)
        elseif sortby == :freq
            sortperm(flutterOffsetFreqs)
        elseif sortby == :mode
            sortperm(flutterOffsetModes)
        end
        flutterOffsetCP = flutterOffsetCP[sortIdx]
        flutterOffsetFreqs = flutterOffsetFreqs[sortIdx]
        flutterOffsetModes = flutterOffsetModes[sortIdx]
        flutterOffsetExtras = map(v -> v[sortIdx], flutterOffsetExtras)
    end

    return flutterOnsetCP,flutterOnsetFreqs,flutterOnsetModes,flutterOnsetExtras...,flutterOffsetCP,flutterOffsetFreqs,flutterOffsetModes,flutterOffsetExtras...
    
end
export find_flutter


"""
    find_divergence(controlParam,freqs,damps; method,atol)

Finds the control parameter values and corresponding modes in damps for which an instability occurs (positive damping), indicating divergence.

# Arguments
- `controlParam::AbstractVector{<:Real}`: vector of control parameter values (length N)
- `freqs::AbstractMatrix{<:Real}`: matrix of frequencies (nTracks x N)
- `damps::AbstractMatrix{<:Real}`: matrix of damping ratios (nTracks x N)

# Keyword arguments
- `method::Symbol`: method to use for finding divergence, either `:regular` or `:backExtrap`
- `atol::Real`: absolute tolerance for zero frequency comparison
"""
function find_divergence(controlParam::AbstractVector{<:Real},freqs::AbstractMatrix{<:Real},damps::AbstractMatrix{<:Real}; method::Symbol=:regular,atol::Real=eps())

    # Validate
    @assert length(controlParam) == size(freqs,2) == size(damps,2) "controlParam, freqs, and damps must have the same length"
    @assert method in (:regular, :backExtrap) "method must be one of :regular or :backExtrap"
    @assert atol >= 0 "atol must be non-negative"

    # Get number of modes and number of control parameter points
    nModes, N = size(freqs)
    
    # Initialize instability values of control parameter and corresponding modes
    divergenceCP = Float64[]
    divergenceModes = Int[]
    divergenceInds = Int[]
    
    # Loop over modes
    for mode in 1:nModes
        # Find the interval [i, i+1] where the zero crossing happens
        if method == :regular
            # Instability: damping goes from negative to positive (with zero frequency)
            k = findfirst(i -> damps[mode,i] < 0 && damps[mode,i+1] > 0 && isapprox(freqs[mode,i],0, atol=atol) && isapprox(freqs[mode,i+1],0), 1:N-1)
            isnothing(k) && continue
            idx = k:k+1
            y = view(damps, mode, idx)
        else #if method == :backExtrap
            # Instability: estimated damping has opposite sign to the actual one (with zero frequency)
            est = backward_extrapolation(damps[mode,:])
            k = findfirst(i -> damps[mode,i] * est[i] < 0 && isapprox(freqs[mode,i],0, atol=atol), 2:N)
            isnothing(k) && continue
            idx = k:k+1
            y = view(est, idx)
        end
        push!(divergenceCP, LinearInterpolations.interpolate(y, view(controlParam, idx), 0))
        push!(divergenceModes, mode)
        push!(divergenceInds, k)
    end

    return divergenceCP, divergenceModes, divergenceInds

end
export find_divergence


"""
    get_FFT_and_PSD(t,y; tol)

Computes the FFT and PSD of signal y(t)

# Arguments
- `t::AbstractVector{<:Real}`: time signal
- `y::Vector{<:Real}`: quantity signal

# Keyword arguments
- `tol::Real`: tolerance for time signal being equally spaced
"""
function get_FFT_and_PSD(t::AbstractVector{<:Real},y::Vector{<:Real}; tol::Real=1e-9)
    
    @assert length(t) > 1
    @assert maximum(abs.(diff(diff(t)))) < tol "t must be an evenly spaced vector"
    @assert length(y) == length(t) "t and y must have the same length"

    # Frequency vector
    # --------------------------------------------------------------------------
    # Duration of time vector
    τ = t[end] - t[1]
    # Time step
    Δt = t[2] - t[1]
    # Minimum detectable frequency [Hz]
    fmin = 1/(τ+Δt)
    # Maximum detectable frequency (Nyquist frequency) [Hz]
    fmax = 1/(2*Δt)
    # Frequency vector [Hz]
    f = collect(0:fmin:fmax)

    # FFT and PSD
    # --------------------------------------------------------------------------
    # Length of time vector
    N = length(t)
    # FFT
    yFFT = FFTW.fft(y)
    # Two-sided spectrum
    ySP2 = abs.(yFFT/N)
    # Single-sided spectrum
    ySP1 = 2*ySP2[1:floor(Int,N/2)+1]
    # PSD
    yPSD = ySP1.^2 * τ/2
    # Adjust single-sided spectrum after calculating the PSD (constant component does not get multiplied by 2)
    ySP1[1] /= 2

    return f,yFFT,yPSD

end
export get_FFT_and_PSD


"""
    blockdiag_dense(mats::AbstractMatrix...) -> Matrix

Creates a block diagonal matrix from the given matrices
"""
function blockdiag_dense(mats::AbstractMatrix...)

    # Set type of output matrix
    T = promote_type(eltype.(mats)...)

    # Number of rows and columns of each matrix
    rows = size.(mats, 1)
    cols = size.(mats, 2)

    # Number of rows and columns of output matrix
    R, C = sum(rows), sum(cols)

    # Initialize output matrix
    out = zeros(T, R, C)

    # Initialize row and column counter
    ro, co = 0, 0

    # Loop over matrices and fill the output matrix
    for (m, r, c) in zip(mats, rows, cols)
        out[ro+1:ro+r, co+1:co+c] .= m
        ro += r
        co += c
    end

    return out

end
export blockdiag_dense


"""
    moving_average(v, n)

Computes the moving average of a series

# Arguments
- `v::AbstractVector{<:Real}`: series
- `n::Int`: moving average window size
"""
function moving_average(v::AbstractVector{<:Real}, n::Int)

    if n == 0
        return v
    end

    @assert iseven(n) "n must be even"
    @assert 2 ≤ n ≤ length(v) "n must be between 2 and length(v)"

    k = div(n, 2)
    result = similar(v, length(v))

    for i in 1:length(v)
        start_idx = max(1, i - k + 1)
        end_idx = min(length(v), i + k)
        result[i] = average(v[start_idx:end_idx])
    end

    return result

end
export moving_average


"""
    local_extrema(x)

Computes the local extrema of a series, and corresponding indices

# Arguments
- `x::AbstractVector{<:Real}`: series
"""
function local_extrema(x::AbstractVector{<:Real})

    n = length(x)
    maximaIdx = Int[]
    minimaIdx = Int[]
    
    for i in 2:(n-1)
        if x[i] > x[i-1] && x[i] > x[i+1]
            push!(maximaIdx, i)
        elseif x[i] < x[i-1] && x[i] < x[i+1]
            push!(minimaIdx, i)
        end
    end
    
    return (maxIdx = maximaIdx, 
            minIdx = minimaIdx, 
            maxVal = x[maximaIdx], 
            minVal = x[minimaIdx])

end
export local_extrema


"""
    backward_extrapolation(v; order)

Estimates the values of an array from backward finite difference extrapolation

# Arguments
- `v::AbstractVector{<:Real}`: vector

# Keyword arguments
- `order::Int`: order of the approximation
"""
function backward_extrapolation(v::AbstractVector{<:Real}; order::Int=5)

    # Size of the vector
    n = length(v)

    # Fix Inf values
    v[v .== Inf .|| v .== -Inf] .= 0

    # Initialize estimated vector
    vEst = copy(v)
    
    # Compute estimated vector
    if order == 3
        for i in 4:n
            vEst[i] = 3v[i-1] - 3v[i-2] + v[i-3]
        end
    elseif order == 4
        for i in 5:n
            vEst[i] = 4v[i-1] - 6v[i-2] + 4v[i-3] - v[i-4]
        end
    elseif order == 5
        for i in 6:n
            vEst[i] = 5v[i-1] - 10v[i-2] + 10v[i-3] - 5v[i-4] + v[i-5]
        end
    else
        error("Set order between 3 and 5")
    end
    
    return vEst

end
export backward_extrapolation


"""
    count_of_exceedance(time, series, datum, marker)

Computes the count of exceedance of the quantities in the marker vector above the datum for the given time series

# Keyword arguments
- `time::AbstractVector{<:Real}`: time vector
- `series::Vector{<:Real}`: time series
- `datum::Real`: datum value for the time series
- `marker::Vector{<:Real}`: markers for the count of exceedance computation
"""
function count_of_exceedance(; time::AbstractVector{<:Real}, series::Vector{<:Real}, datum::Real, marker::Vector{<:Real})

    # Validate
    @assert length(time) == length(series) "time and series vectors must be of the same length"
    @assert issorted(marker) "Marker vector must be in ascending order"

    # Initialize
    n = length(series)
    counts = zeros(Int, length(marker))
    aboveDatum = series[1] > datum
    regionStart = 1

    # Sweep time 
    for i in 2:n
        isAbove = series[i] > datum
        if isAbove != aboveDatum
            region = series[regionStart:i-1]
            if aboveDatum
                peak = maximum(region)
                for (j, m) in enumerate(marker)
                    if m < peak
                        counts[j] += 1
                    end
                end
            end
            regionStart = i
            aboveDatum = isAbove
        end
    end

    # Handle last region
    region = series[regionStart:end]
    if aboveDatum
        peak = maximum(region)
        for (j, m) in enumerate(marker)
            if m < peak
                counts[j] += 1
            end
        end
    end

    return counts

end
export count_of_exceedance


"""
    frequency_of_exceedance(timeLength, count)

Computes the frequency of exceedance

# Keyword arguments
- `count::Vector{<:Real}`: count of exceedance
- `timeLength::Real`: length of time
"""
function frequency_of_exceedance(; count::Vector{<:Real}, timeLength::Real)
    return count / timeLength
end
export frequency_of_exceedance


"""
    gcd_floats_approx(v; tol)

Computes the approximate greatest common divisor of a vector of floats

# Arguments
- `v::Vector{Float64}`: vector

# Keyword arguments
- `tol::Real`: tolerance
"""
function gcd_floats_approx(v::Vector{Float64}; tol::Real=1e-8)
    scale = 10.0 ^ ceil(Int, -log10(tol))
    ints = round.(Int, v .* scale)
    g = reduce(gcd, ints)
    return g / scale
end
export gcd_floats_approx


"""
    spectrogram(t, x; window, overlap, db)

Compute the spectrogram (Short-Time Fourier Transform) of signal `x` sampled at times `t`

# Arguments
- `t::AbstractVector{<:Real}`: time vector (must be equally spaced)
- `x::AbstractVector{<:Real}`: signal values

# Keyword arguments
- `window::Int`: window length in samples (default: 256)
- `overlap::Real`: fraction of window overlap (0–1, default: 0.5)
- `db::Bool`: if true, return power in decibels (default: true)

# Outputs
- `f::Vector`: frequency vector
- `tt::Vector`: time vector (centers of analysis windows)
- `S::Matrix`: spectrogram matrix (|STFT|²), size (length(f), length(tt))
"""
function spectrogram(t::AbstractVector{<:Real}, x::AbstractVector{<:Real}; window::Int=256, overlap::Real=0.5, db::Bool=true)

    # Validate
    @assert length(t) == length(x) "time vector and signal must have the same length"
    dt = diff(t)
    @assert maximum(abs.(dt .- dt[1])) < 1e-8*dt[1] "time vector must be equally spaced"
    @assert 0 < overlap < 1 "overlap must be between 0 and 1"

    # Sampling frequency
    fs = 1/(t[2]-t[1])

    # Length of overlap (hop size)
    hop = round(Int, window*(1-overlap))

    # Compute STFT
    X = DSP.stft(x, window, hop; fs=fs, window=DSP.hanning, onesided=true)

    # Frequency and time vectors
    nfft = window
    f = collect(0:fs/nfft:fs/2)
    tt = collect(0:hop:(size(X,2)-1)*hop) ./ fs

    # Power spectrogram
    S = abs.(X).^2
    if db
        S = 10 .* log10.(S .+ eps())  # dB scale
    end

    return f, tt, S

end
export spectrogram


"""
    NASA_frames_loader(frame)

Loads dynamic stall data from a frame of the NASA report (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19830003778.pdf)

# Arguments
- `frame::Union{Int,String}`: frame number
"""
function NASA_frames_loader(frame::Union{Int,String})

    # Atmosphere
    altitude = 0
    atmosphere = standard_atmosphere(altitude)

    # Frame data lookup (airfoil, mean angle, angle amplitude, airfoil semichord, reduced frequency, Mach number)
    frames = Dict(
        10_022 => (airfoil=deepcopy(NACA0012), a₀=12π/180, a₁=9.9π/180, b=0.305, k=0.098, Ma=0.301)
    )

    # Retrieve frame data
    frame_data = get(frames, frame, nothing)
    isnothing(frame_data) && error("Frame $frame is unavailable")

    # Extract variables
    airfoil,a₀,a₁,b,k,Ma = frame_data

    # Airspeed
    U = Ma*atmosphere.a

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,Ma)

    return airfoil,a₀,a₁,b,k,Ma,U

end
export NASA_frames_loader


"""
    GU_frames_loader(frame)

Loads dynamic stall data from a frame of the University of Glasgow (DOI:10.5525/gla.researchdata.464)

# Arguments
- `frame::Union{Int,String}`: frame number
"""
function GU_frames_loader(frame::Union{Int,String})

    if frame isa String
        frame = parse(Int, replace(frame, "_" => ""))
    end

    # Frame data lookup (airfoil, mean angle, angle amplitude, airfoil semichord, reduced frequency, Mach number, airspeed)
    frames = Dict(
        02_000_101 => (airfoil=deepcopy(NACA23012A), a₀=0.24817, a₁=0.26967, b=0.275, k=0.001, Ma=0.11931, U=41.570),
        02_010_151 => (airfoil=deepcopy(NACA23012A), a₀=0.18100, a₁=0.10401, b=0.275, k=0.12489, Ma=0.11635, U=40.605),
        02_010_351 => (airfoil=deepcopy(NACA23012A), a₀=0.17899, a₁=0.17741, b=0.275, k=0.174, Ma=0.11694, U=40.815)
    )

    # Retrieve frame data
    frame_data = get(frames, frame, nothing)
    isnothing(frame_data) && error("Frame $frame is unavailable")

    # Extract variables
    airfoil,a₀,a₁,b,k,Ma,U = frame_data

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,Ma)

    return airfoil,a₀,a₁,b,k,Ma,U

end
export GU_frames_loader