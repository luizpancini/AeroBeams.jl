using Revise, Parameters, BenchmarkTools, SparseArrays, LinearAlgebra, Statistics, Plots, QuadGK, LinearInterpolations, ForwardDiff

"""
const a1 = [1.0; 0.0; 0.0]

First vector of basis A, resolved in that basis
"""
const a1 = [1.0; 0.0; 0.0]


"""
I3 = Matrix(1.0*LinearAlgebra.I,3,3)

3x3 identity matrix
"""
const I3 = Matrix(1.0*LinearAlgebra.I,3,3)


"""
I6 = Matrix(1.0*LinearAlgebra.I,6,6)

6x6 identity matrix
"""
const I6 = Matrix(1.0*LinearAlgebra.I,6,6)


"""
round_off!(x)

Rounds the array to input tolerance (defaults to machine epsilon)

# Arguments
- x
"""
function round_off!(x,tol=eps())

    x[abs.(x).<tol] .= 0.0

    return x
end


"""
divide_inplace(divisor::Number, vars...)

Divides the input variables in-place

# Arguments
- divisor::Number = divisor
- vars... = variables to be divided
"""
function divide_inplace(divisor::Number, vars...)
    return (var ./ divisor for var in vars)
end


"""
tilde(v::Vector{<:Number})

Gets the skew-symmetric matrix associated with a vector

# Arguments
- v::Vector{<:Number} = three-element vector
"""
function tilde(v::Vector{<:Number})
    return [  0.0  -v[3]   v[2]
             v[3]    0.0  -v[1]
            -v[2]   v[1]   0.0]
end
export tilde


"""
mul3(A1::Matrix{<:Number},A2::Matrix{<:Number},A3::Matrix{<:Number},b::Vector{<:Number})

Gets the scalar product of a third-order tensor represented by matrices A1, A2, and A3 with the vector b

# Arguments
- A1::Matrix{<:Number}
- A2::Matrix{<:Number}
- A3::Matrix{<:Number}
- b::Vector{<:Number}
"""
function mul3(A1::Matrix{<:Number},A2::Matrix{<:Number},A3::Matrix{<:Number},b::Vector{<:Number})
    
    return hcat(A1*b, A2*b, A3*b)
    
end

"""
curvature_quantities(k::Vector{Float64})

Gets some useful functions of the curvature vector

# Arguments
- k::Vector{Float64} = curvature vector
"""
function curvature_quantities(k::Vector{Float64})

    # External product with itself
    kkT = k * k'

    # Skew-symmetric matrix 
    ktilde = tilde(k)

    # Norm 
    knorm = sqrt(dot(k,k))

    return kkT, ktilde, knorm
end

"""
position_vector_from_curvature(R0::Matrix{Float64}, k::Vector{Float64}, x1::Float64)

Gets the position vector at an arclength value

# Arguments
- R0::Matrix{Float64} = initial rotation tensor (at arclength position zero)
- k::Vector{Float64} = curvature vector
- x1::Float64 = arclength position
"""
function position_vector_from_curvature(R0::Matrix{Float64}, k::Vector{Float64}, x1::Float64) 

    kkT, ktilde, knorm = curvature_quantities(k)

    if knorm > 0
        return R0 * (1.0/knorm*(I3 .- kkT/knorm^2)*sin(knorm*x1) .+ ktilde/knorm^2*(1.0-cos(knorm*x1)) .+ kkT/knorm^2*x1) * a1
    else
        return x1*R0[:,1]
    end
end


"""
rotation_tensor_from_curvature(R0::Matrix{Float64}, k::Vector{Float64}, x1::Float64)

Gets the rotation tensor at an arclength position 

# Arguments
- R0::Matrix{Float64} = initial rotation tensor (at arclength position zero)
- k::Vector{Float64} = curvature vector
- x1::Float64 = arclength position
"""
function rotation_tensor_from_curvature(R0::Matrix{Float64},k::Vector{Float64},x1::Float64) 

    kkT, ktilde, knorm = curvature_quantities(k)

    # Set according to |k|
    if knorm > 0.0
        R = R0 * ((I3 .- kkT/knorm^2)*cos(knorm*x1) .+ ktilde/knorm*sin(knorm*x1) .+ kkT/knorm^2)
    else
        R = R0
    end

    # Round off
    round_off!(R)

    return R
end


"""
rotation_tensor_E321(p::Vector{<:Union{Float64,ComplexF64}})

Gets the rotation tensor according to Euler parameters sequence 3-2-1

# Arguments
- p::Vector{<:Union{Float64,ComplexF64}} = rotation parameters
"""
function rotation_tensor_E321(p::Vector{<:Union{Float64,ComplexF64}})

    # Euler angles 3-2-1 sequence (yaw, pitch and roll angles), respective sines and cosines
    yaw = p[1]
    cy = cos(yaw)
    sy = sin(yaw)
    pitch = p[2]
    cp = cos(pitch)
    sp = sin(pitch)
    roll = p[3]
    cr = cos(roll)
    sr = sin(roll)

    # Rotation tensor that brings the reference basis to the final basis
    R = [cp*cy cy*sr*sp - cr*sy sr*sy + cr*cy*sp;
         cp*sy cr*cy + sr*sp*sy cr*sp*sy - cy*sr;
           -sp            cp*sr            cr*cp]

    round_off!(R)       

    return R

end


"""
rotation_tensor_E321(p::::Vector{<:Union{Float64,ComplexF64}})

Gets the rotation tensor according to Euler parameters sequence 3-1-3

# Arguments
- p::Vector{<:Union{Float64,ComplexF64}} = rotation parameters
"""
function rotation_tensor_E313(p::Vector{<:Union{Float64,ComplexF64}})

    # Euler angles 3-1-3 sequence (precession=ϕ, nutation=θ and spin=ψ angles), respective sines and cosines
    ϕ = p[1]
    cϕ = cos(ϕ)
    sϕ = sin(ϕ)   
    θ = p[2]
    cθ = cos(θ)
    sθ = sin(θ)  
    ψ = p[3]
    cψ = cos(ψ)
    sψ = sin(ψ)
       
    # Rotation tensor that brings the reference basis to the final basis
    R = [-sϕ*cθ*sψ+cϕ*cψ -sϕ*cθ*cψ-cϕ*sψ  sψ*sθ;
          cϕ*cθ*sψ+sϕ*cψ  cϕ*cθ*cψ-sϕ*sψ -cϕ*sθ;
                   sθ*sψ           sθ*cψ     cθ]

    round_off!(R)

    return R

end


"""
rotation_tensor_WM(p::Vector{<:Union{Float64,ComplexF64}})

Gets the rotation tensor according to Wiener-Milenkovic parameters

# Arguments
- p::Vector{<:Union{Float64,ComplexF64}} = rotation parameters
"""
function rotation_tensor_WM(p::Vector{<:Union{Float64,ComplexF64}})

    # Scaling factor and scaled rotation parameters
    λ, ps, _ , pNorm, psNorm = rotation_parameter_scaling(p)

    # Shorthand for scaled rotation parameters 
    ps1 = ps[1]
    ps2 = ps[2]
    ps3 = ps[3]

    # Useful functions of the rotation parameters 
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
    Θ = [ps0s+ps1s-ps2s-ps3s    2*(ps1ps2-ps0ps3)    2*(ps1ps3+ps0ps2);
           2*(ps1ps2+ps0ps3)  ps0s-ps1s+ps2s-ps3s    2*(ps2ps3-ps0ps1);
           2*(ps1ps3-ps0ps2)    2*(ps2ps3+ps0ps1)  ps0s-ps1s-ps2s+ps3s]
        
    # Rotation tensor from Wiener-Milenkovic parameters 
    R = υ² * Θ

    return R,Θ,pNorm,λ,ps,ps1,ps2,ps3,ps0,υ,υ²,ps1s,ps2s,ps3s,ps1ps2,ps2ps3,ps1ps3

end


"""
rotation_parameter_scaling(p::Vector{<:Union{Float64,ComplexF64}})

Scales the Wiener-Milenkovic rotation parameters

# Arguments
- p::Vector{<:Union{Float64,ComplexF64}} = unscaled rotation parameters
"""
function rotation_parameter_scaling(p::Vector{<:Union{Float64,ComplexF64}})

    # Initialize scaling factor and number of odd half rotations
    λ = 1.0
    halfRotations = 0
    
    # Norm of the unscaled/extended parameters vector
    pNorm = norm(p)
    
    # Scale according to norm
    if pNorm > 4.0  
        halfRotations = round(pNorm/8.0)
        λ = 1.0 - 8.0 * halfRotations / pNorm
    end
    
    # Scaled parameters
    ps = λ * p
    
    # Norm of scaled parameters vector
    psNorm = abs(λ) * pNorm

    return λ,ps,halfRotations,pNorm,psNorm
    
end


"""
rotation_tensor_derivatives_scaled_parameters(ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²)

Gets the derivatives of the rotation tensor with respect to the scaled rotation parameters

# Arguments
- ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²
"""
function rotation_tensor_derivatives_scaled_parameters(ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²)

    # Derivatives of υ² w.r.t. scaled rotation parameters
    υ²_p = -1/2*υ^3*ps
    υ²_ps1 = υ²_p[1]
    υ²_ps2 = υ²_p[2]
    υ²_ps3 = υ²_p[3]

    # Local temporary constants
    c1 = 1-ps0/4
    c2 = 1+ps0/4
    p1so4 = ps1s/4
    p2so4 = ps2s/4
    p3so4 = ps3s/4
    p1p2o4 = ps1ps2/4
    p1p3o4 = ps1ps3/4
    p2p3o4 = ps2ps3/4

    # Rotation tensor's core' derivatives w.r.t. scaled rotation parameters
    Θ_ps1 = 2*[     ps1*c1   ps2+p1p3o4    ps3-p1p2o4;
                ps2-p1p3o4      -ps1*c2  -(ps0-p1so4);
                ps3+p1p2o4    ps0-p1so4       -ps1*c2]
    Θ_ps2 = 2*[    -ps2*c2   ps1+p2p3o4   ps0-p2so4;
                ps1-p2p3o4       ps2*c1  ps3+p1p2o4;
              -(ps0-p2so4)   ps3-p1p2o4     -ps2*c2]
    Θ_ps3 = 2*[    -ps3*c2  -(ps0-p3so4)  ps1-p2p3o4;
                 ps0-p3so4      -ps3*c2   ps2+p1p3o4;
                ps1+p2p3o4   ps2-p1p3o4       ps3*c1]

    # Rotation tensor derivatives w.r.t. scaled rotation parameters         
    R_ps1 = υ²_ps1*Θ .+ υ²*Θ_ps1
    R_ps2 = υ²_ps2*Θ .+ υ²*Θ_ps2
    R_ps3 = υ²_ps3*Θ .+ υ²*Θ_ps3

    return R_ps1,R_ps2,R_ps3,υ²_ps1,υ²_ps2,υ²_ps3,Θ_ps1,Θ_ps2,Θ_ps3
    
end


"""
scaling_derivatives_extended_parameters(λ::Float64,p::Vector{Float64},pNorm::Float64)

Gets the derivatives of the scaling factor with respect to the extended rotation parameters

# Arguments
- λ::Float64 = scaling factor
- p::Vector{Float64} = rotation parameters
- pNorm::Float64 = norm of rotation parameters
"""
function scaling_derivatives_extended_parameters(λ::Float64,p::Vector{Float64},pNorm::Float64)

    if λ == 1.0
        return zeros(3)
    else 
        return (1-λ)*p/pNorm^2
    end

end


"""
rotation_tensor_derivatives_extended_parameters(p,pNorm,λ,ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²)

Gets the derivatives of the rotation tensor with respect to the extended rotation parameters

# Arguments
- p,pNorm,λ,ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²
"""
function rotation_tensor_derivatives_extended_parameters(p,pNorm,λ,ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²)

    # Rotation tensor derivatives w.r.t. scaled rotation parameters
    R_ps1,R_ps2,R_ps3,υ²_ps1,υ²_ps2,υ²_ps3,Θ_ps1,Θ_ps2,Θ_ps3 = rotation_tensor_derivatives_scaled_parameters(ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²)

    # Scaling derivatives w.r.t. extended rotation parameters
    λ_p = scaling_derivatives_extended_parameters(λ,p,pNorm)

    # Scaled rotation parameters derivatives w.r.t. extended rotation parameters 
    ps_p = λ_p*p' .+ λ*I3
    ps1_p1 = ps_p[1,1]
    ps2_p1 = ps_p[2,1]
    ps3_p1 = ps_p[3,1]
    ps1_p2 = ps_p[1,2]
    ps2_p2 = ps_p[2,2]
    ps3_p2 = ps_p[3,2]
    ps1_p3 = ps_p[1,3]
    ps2_p3 = ps_p[2,3]
    ps3_p3 = ps_p[3,3]

    # Rotation tensor derivatives w.r.t. extended rotation parameters
    R_p1 = R_ps1*ps1_p1 .+ R_ps2*ps2_p1 .+ R_ps3*ps3_p1
    R_p2 = R_ps1*ps1_p2 .+ R_ps2*ps2_p2 .+ R_ps3*ps3_p2
    R_p3 = R_ps1*ps1_p3 .+ R_ps2*ps2_p3 .+ R_ps3*ps3_p3

    return R_p1,R_p2,R_p3,R_ps1,R_ps2,R_ps3,υ²_ps1,υ²_ps2,υ²_ps3,Θ_ps1,Θ_ps2,Θ_ps3,ps_p,ps1_p1,ps2_p1,ps3_p1,ps1_p2,ps2_p2,ps3_p2,ps1_p3,ps2_p3,ps3_p3
end


"""
rotation_tensor_derivatives_extended_parameters(p::Vector{Float64})

Gets the derivatives of the rotation tensor with respect to the extended rotation parameters

# Arguments
- p::Vector{Float64}
"""
function rotation_tensor_derivatives_extended_parameters(p::Vector{Float64})

    _,Θ,pNorm,λ,ps,ps1,ps2,ps3,ps0,υ,υ²,ps1s,ps2s,ps3s,ps1ps2,ps2ps3,ps1ps3 = rotation_tensor_WM(p)

    # Rotation tensor derivatives w.r.t. scaled rotation parameters
    R_ps1,R_ps2,R_ps3,_ = rotation_tensor_derivatives_scaled_parameters(ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²)

    # Scaling derivatives w.r.t. extended rotation parameters
    λ_p = scaling_derivatives_extended_parameters(λ,p,pNorm)

    # Scaled rotation parameters derivatives w.r.t. extended rotation parameters 
    ps_p = λ_p*p' .+ λ*I3
    ps1_p1 = ps_p[1,1]
    ps2_p1 = ps_p[2,1]
    ps3_p1 = ps_p[3,1]
    ps1_p2 = ps_p[1,2]
    ps2_p2 = ps_p[2,2]
    ps3_p2 = ps_p[3,2]
    ps1_p3 = ps_p[1,3]
    ps2_p3 = ps_p[2,3]
    ps3_p3 = ps_p[3,3]

    # Rotation tensor derivatives w.r.t. extended rotation parameters
    R_p1 = R_ps1*ps1_p1 .+ R_ps2*ps2_p1 .+ R_ps3*ps3_p1
    R_p2 = R_ps1*ps1_p2 .+ R_ps2*ps2_p2 .+ R_ps3*ps3_p2
    R_p3 = R_ps1*ps1_p3 .+ R_ps2*ps2_p3 .+ R_ps3*ps3_p3

    return R_p1,R_p2,R_p3
end


"""
rotation_tensor_time_derivative(R_ps1,R_ps2,R_ps3,ps_p,pdot)

Gets the time derivative of the rotation tensor 

# Arguments
- R_ps1,R_ps2,R_ps3,ps_p,pdot
"""
function rotation_tensor_time_derivative(R_ps1,R_ps2,R_ps3,ps_p,pdot)

    # Scaled rotation parameters time derivative
    psdot = ps_p*pdot
    ps1dot = psdot[1]
    ps2dot = psdot[2]
    ps3dot = psdot[3]

    # Rotation tensor time derivative 
    Rdot = R_ps1*ps1dot .+ R_ps2*ps2dot .+ R_ps3*ps3dot

    return Rdot,ps1dot,ps2dot,ps3dot

end


"""
rotation_tensor_time_derivative(p::Vector{Float64},pdot::Vector{Float64})

Gets the derivatives of the time derivative of the rotation tensor, given the rotation parameters and their rates 

# Arguments
- p::Vector{Float64}
- pdot::Vector{Float64}
"""
function rotation_tensor_time_derivative(p::Vector{Float64},pdot::Vector{Float64})

    # Rotation parameters' variables
    _,Θ,pNorm,λ,ps,ps1,ps2,ps3,ps0,υ,υ²,ps1s,ps2s,ps3s,ps1ps2,ps2ps3,ps1ps3 = rotation_tensor_WM(p)

    # Rotation tensor and rotation parameters' derivatives w.r.t. extended rotation parameters
    _,_,_,R_ps1,R_ps2,R_ps3,_,_,_,_,_,_,ps_p,_,_,_,_,_,_,_,_,_ = rotation_tensor_derivatives_extended_parameters(p,pNorm,λ,ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²)

    # Rotation tensor time derivative 
    Rdot,_ = rotation_tensor_time_derivative(R_ps1,R_ps2,R_ps3,ps_p,pdot)

    return Rdot
end


"""
rotation_tensor_derivatives_time_extended_parameters(ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,Θ_ps1,Θ_ps2,Θ_ps3,υ,υ²,υ²_ps1,υ²_ps2,υ²_ps3,ps1dot,ps2dot,ps3dot,ps1_p1,ps2_p1,ps3_p1,ps1_p2,ps2_p2,ps3_p2,ps1_p3,ps2_p3,ps3_p3)

Gets the derivatives of the time derivative of the rotation tensor with respect to the extended rotation parameters 

# Arguments
- ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,Θ_ps1,Θ_ps2,Θ_ps3,υ,υ²,υ²_ps1,υ²_ps2,υ²_ps3,ps1dot,ps2dot,ps3dot,ps1_p1,ps2_p1,ps3_p1,ps1_p2,ps2_p2,ps3_p2,ps1_p3,ps2_p3,ps3_p3
"""
function rotation_tensor_derivatives_time_extended_parameters(ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,Θ_ps1,Θ_ps2,Θ_ps3,υ,υ²,υ²_ps1,υ²_ps2,υ²_ps3,ps1dot,ps2dot,ps3dot,ps1_p1,ps2_p1,ps3_p1,ps1_p2,ps2_p2,ps3_p2,ps1_p3,ps2_p3,ps3_p3)

    # υ² second derivatives w.r.t. scaled rotation parameters
    υ²_psps = -1/2*υ^3*(I3-3/4*υ*ps*(ps'))

    # Local temporary constants
    c1 = 1-ps0/4
    c2 = 1+ps0/4
    ps1o4 = ps1/4
    ps2o4 = ps2/4
    ps3o4 = ps3/4
    ps1so16 = ps1s/16
    ps2so16 = ps2s/16
    ps3so16 = ps3s/16
    ps1ps2o16 = ps1ps2/16
    ps1ps3o16 = ps1ps3/16
    ps2ps3o16 = ps2ps3/16

    # Rotation tensor core's second derivatives w.r.t. scaled rotation parameters
    Θ_ps1ps1 = 2*[ps1so16+c1      ps3o4     -ps2o4;
                      -ps3o4 ps1so16-c2    3*ps1o4;
                       ps2o4   -3*ps1o4 ps1so16-c2]                 
    Θ_ps2ps2 = 2*[ps2so16-c2      ps3o4   -3*ps2o4;
                      -ps3o4 ps2so16+c1      ps1o4;
                     3*ps2o4     -ps1o4 ps2so16-c2]   
    Θ_ps3ps3 = 2*[ps3so16-c2    3*ps3o4     -ps2o4;
                    -3*ps3o4 ps3so16-c2      ps1o4;
                       ps2o4     -ps1o4 ps3so16+c1]
    Θ_ps1ps2 = 2*[ps1ps2o16           1     -ps1o4;
                          1   ps1ps2o16      ps2o4;
                      ps1o4      -ps2o4  ps1ps2o16]   
    Θ_ps1ps3 = 2*[ps1ps3o16       ps1o4          1;
                     -ps1o4   ps1ps3o16      ps3o4;
                          1      -ps3o4  ps1ps3o16]  
    Θ_ps2ps3 = 2*[ps2ps3o16       ps2o4     -ps3o4;
                     -ps2o4   ps2ps3o16          1;
                      ps3o4           1  ps2ps3o16]

    # Rotation tensor second derivatives w.r.t scaled rotation parameters           
    R_ps1ps1 = υ²_psps[1,1]*Θ .+ 2*υ²_ps1*Θ_ps1 .+ υ²*Θ_ps1ps1
    R_ps2ps2 = υ²_psps[2,2]*Θ .+ 2*υ²_ps2*Θ_ps2 .+ υ²*Θ_ps2ps2
    R_ps3ps3 = υ²_psps[3,3]*Θ .+ 2*υ²_ps3*Θ_ps3 .+ υ²*Θ_ps3ps3
    R_ps2ps1 = υ²_psps[1,2]*Θ .+ υ²_ps2*Θ_ps1 .+ υ²_ps1*Θ_ps2 .+ υ²*Θ_ps1ps2
    R_ps3ps1 = υ²_psps[1,3]*Θ .+ υ²_ps3*Θ_ps1 .+ υ²_ps1*Θ_ps3 .+ υ²*Θ_ps1ps3
    R_ps3ps2 = υ²_psps[2,3]*Θ .+ υ²_ps3*Θ_ps2 .+ υ²_ps2*Θ_ps3 .+ υ²*Θ_ps2ps3
    R_ps1ps2 = R_ps2ps1
    R_ps1ps3 = R_ps3ps1
    R_ps2ps3 = R_ps3ps2

    # Rotation tensor time derivatives' derivative w.r.t. scaled rotation parameters
    Rdot_ps1 = R_ps1ps1*ps1dot .+ R_ps2ps1*ps2dot .+ R_ps3ps1*ps3dot
    Rdot_ps2 = R_ps1ps2*ps1dot .+ R_ps2ps2*ps2dot .+ R_ps3ps2*ps3dot
    Rdot_ps3 = R_ps1ps3*ps1dot .+ R_ps2ps3*ps2dot .+ R_ps3ps3*ps3dot

    # Rotation tensor time derivatives' derivative w.r.t. extended rotation parameters
    Rdot_p1 = Rdot_ps1*ps1_p1 .+ Rdot_ps2*ps2_p1 .+ Rdot_ps3*ps3_p1
    Rdot_p2 = Rdot_ps1*ps1_p2 .+ Rdot_ps2*ps2_p2 .+ Rdot_ps3*ps3_p2
    Rdot_p3 = Rdot_ps1*ps1_p3 .+ Rdot_ps2*ps2_p3 .+ Rdot_ps3*ps3_p3

    return Rdot_p1,Rdot_p2,Rdot_p3

end


"""
tangent_operator_transpose_WM(p::Vector{Float64})

Gets the transpose of tangent operator tensor according to Wiener-Milenkovic parameters

# Arguments
- p::Vector{Float64} = rotation parameters
"""
function tangent_operator_transpose_WM(p::Vector{Float64})

    # Scaling factor and scaled rotation parameters
    _, ps, _, _, psNorm = rotation_parameter_scaling(p)

    # Useful functions of the scaled rotation parameters 
    ps0 = 2 - psNorm^2/8
    υ² = (4-ps0)^-2

    # Tangent operator transpose
    return 2*υ²*(ps0*I3 .+ 1/4*ps*(ps') .- tilde(ps))

end


"""
tangent_operator_transpose_WM(ps::Vector{Float64},ps0::Float64,υ²::Float64)

Gets the transpose of tangent operator tensor according to Wiener-Milenkovic parameters

# Arguments
- ps::Vector{Float64} = scaled rotation parameters
- ps0::Float64
- υ²::Float64
"""
function tangent_operator_transpose_WM(ps::Vector{Float64},ps0::Float64,υ²::Float64)

    return 2*υ²*(ps0*I3 .+ 1/4*ps*(ps') .- tilde(ps))

end


"""
tangent_operator_transpose_inverse_WM(p::Vector{<:Union{Float64,ComplexF64}})

Gets the inverse of the transpose of the tangent operator tensor according to Wiener-Milenkovic parameters

# Arguments
- p::Vector{<:Union{Float64,ComplexF64}} = rotation parameters
"""
function tangent_operator_transpose_inverse_WM(p::Vector{<:Union{Float64,ComplexF64}})

    # Scaling factor and scaled rotation parameters
    _, ps, _, _, psNorm = rotation_parameter_scaling(p)

    # Useful functions of the scaled rotation parameters 
    ps0 = 2 - psNorm^2/8

    # Tangent operator inverse transpose
    return 1/2*(ps0*I3 .+ 1/4*ps*(ps') .+ tilde(ps))

end


"""
tangent_operator_transpose_inverse_WM(ps::Vector{Float64},ps0::Float64)

Gets the inverse of the transpose of the tangent operator tensor according to Wiener-Milenkovic parameters

# Arguments
- ps::Vector{Float64} = scaled rotation parameters
- ps0::Float64
"""
function tangent_operator_transpose_inverse_WM(ps::Vector{Float64},ps0::Float64)

    return 1/2*(ps0*I3 .+ 1/4*ps*(ps') .+ tilde(ps))

end


"""
tangent_tensor_functions_derivatives_extended_parameters(HT,ps1,ps2,ps3,υ²,υ²_ps1,υ²_ps2,υ²_ps3,ps_p)

Gets the derivatives of the tangent tensor's transpose and its inverse with respect to the extended rotation parameters

# Arguments
- HT,ps1,ps2,ps3,υ²,υ²_ps1,υ²_ps2,υ²_ps3,ps_p
"""
function tangent_tensor_functions_derivatives_extended_parameters(HT,ps1,ps2,ps3,υ²,υ²_ps1,υ²_ps2,υ²_ps3,ps_p)

    # Tangent operator transpose derivatives w.r.t. scaled rotation parameters
    HT_ps1 = υ²_ps1*HT/υ² .+ υ²*1/2*[ps1 ps2 ps3; ps2 -ps1 4; ps3 -4 -ps1]
    HT_ps2 = υ²_ps2*HT/υ² .+ υ²*1/2*[-ps2 ps1 -4; ps1 ps2 ps3; 4 ps3 -ps2]
    HT_ps3 = υ²_ps3*HT/υ² .+ υ²*1/2*[-ps3 4 ps1; -4 -ps3 ps2; ps1 ps2 ps3]

    # Tangent operator transpose inverse derivatives w.r.t. scaled rotation parameters
    HTinv_ps1 = 1/8*[ps1 ps2 ps3; ps2 -ps1 -4; ps3 4 -ps1]
    HTinv_ps2 = 1/8*[-ps2 ps1 4; ps1 ps2 ps3; -4 ps3 -ps2]
    HTinv_ps3 = 1/8*[-ps3 -4 ps1; 4 -ps3 ps2; ps1 ps2 ps3]

    # Tangent operator transpose derivatives w.r.t. extended rotation parameters
    HT_p1 = HT_ps1*ps_p[1,1] .+ HT_ps2*ps_p[2,1] .+ HT_ps3*ps_p[3,1]
    HT_p2 = HT_ps1*ps_p[1,2] .+ HT_ps2*ps_p[2,2] .+ HT_ps3*ps_p[3,2]
    HT_p3 = HT_ps1*ps_p[1,3] .+ HT_ps2*ps_p[2,3] .+ HT_ps3*ps_p[3,3]

    # Tangent operator transpose inverse derivatives w.r.t. extended rotation parameters
    HTinv_p1 = HTinv_ps1*ps_p[1,1] .+ HTinv_ps2*ps_p[2,1] .+ HTinv_ps3*ps_p[3,1]
    HTinv_p2 = HTinv_ps1*ps_p[1,2] .+ HTinv_ps2*ps_p[2,2] .+ HTinv_ps3*ps_p[3,2]
    HTinv_p3 = HTinv_ps1*ps_p[1,3] .+ HTinv_ps2*ps_p[2,3] .+ HTinv_ps3*ps_p[3,3]

    return HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3

end


"""
tangent_tensor_transpose_derivatives_extended_parameters(p::Vector{Float64})

Gets the derivatives of the tangent tensor's transpose with respect to the extended rotation parameters

# Arguments
- p::Vector{Float64}
"""
function tangent_tensor_transpose_derivatives_extended_parameters(p::Vector{Float64})

    # Scaling factor and scaled rotation parameters
    λ,ps,_,pNorm,psNorm = rotation_parameter_scaling(p)

    # Scaling derivatives w.r.t. extended rotation parameters
    λ_p = scaling_derivatives_extended_parameters(λ,p,pNorm)

    # Scaled parameters derivatives w.r.t. extended rotation parameters
    ps_p = λ_p*p' + λ*I3
    
    # Useful functions of the scaled rotation parameters 
    ps0 = 2 - psNorm^2/8
    ps1 = ps[1]
    ps2 = ps[2]
    ps3 = ps[3]
    υ = 1/(4-ps0)
    υ² = υ^2
    υ²_p = -1/2*υ^3*ps
    υ²_ps1 = υ²_p[1]
    υ²_ps2 = υ²_p[2]
    υ²_ps3 = υ²_p[3]

    # Tangent operator transpose
    HT = tangent_operator_transpose_WM(ps,ps0,υ²)
    
    # Tangent operator transpose derivatives w.r.t. scaled rotation parameters
    HT_ps1 = υ²_ps1*HT/υ² .+ υ²*1/2*[ps1 ps2 ps3; ps2 -ps1 4; ps3 -4 -ps1]
    HT_ps2 = υ²_ps2*HT/υ² .+ υ²*1/2*[-ps2 ps1 -4; ps1 ps2 ps3; 4 ps3 -ps2]
    HT_ps3 = υ²_ps3*HT/υ² .+ υ²*1/2*[-ps3 4 ps1; -4 -ps3 ps2; ps1 ps2 ps3]

    # Tangent operator transpose derivatives w.r.t. extended rotation parameters
    HT_p1 = HT_ps1*ps_p[1,1] .+ HT_ps2*ps_p[2,1] .+ HT_ps3*ps_p[3,1]
    HT_p2 = HT_ps1*ps_p[1,2] .+ HT_ps2*ps_p[2,2] .+ HT_ps3*ps_p[3,2]
    HT_p3 = HT_ps1*ps_p[1,3] .+ HT_ps2*ps_p[2,3] .+ HT_ps3*ps_p[3,3]

    return HT_p1,HT_p2,HT_p3

end

"""
force_scaling(S::Vector{Matrix{Float64}})

Gets the appropriate force scaling for the linear system of equations

# Arguments
- S::Vector{Matrix{Float64}} = array of elemental sectional compliance matrices
"""
function force_scaling(S::Vector{Matrix{Float64}})

    # Create array with compliance entries
    complianceEntries = vcat(vcat(S...)...)
            
    # Get all nonzero entries
    complianceNonzeroIndices = findall( x -> x != 0.0, complianceEntries)
    
    # Set force scaling based on number of nonzero compliance matrix entries
    forceScaling = 1.0
    if !isempty(complianceNonzeroIndices)
        
        nonzeroComplianceEntries = complianceEntries[complianceNonzeroIndices]

        α = length(nonzeroComplianceEntries)/sum(nonzeroComplianceEntries)/100

        if α > 2
            forceScaling = nextpow(2,α)
        else
            forceScaling = 2/nextpow(2,1/α)
        end
    end
    
    return forceScaling
end


"""
rotation_angle(p::Vector{<:Union{Float64,ComplexF64}})

Gets the rotation angle

# Arguments
- p::Vector{<:Union{Float64,ComplexF64}} = rotation parameters
"""
function rotation_angle(p::Vector{<:Union{Float64,ComplexF64}})

    # Scale rotation parameters and get number of half rotations
    _,_,halfRotations,pNorm,_ = rotation_parameter_scaling(p)

    if pNorm > 0 && count(!iszero,p/pNorm) == 1 && any(x -> x == -1.0, real.(p/pNorm))
        θ_sign = -1
    else
        θ_sign = 1
    end

    return θ_sign * (4*atan((pNorm-8*halfRotations)/4) + 2*π*halfRotations)

end


"""
rotation_parameters_WM(R::Matrix{<:Number})

Computes the Wiener-Milenkovic rotation parameters given a rotation tensor

# Arguments
- R::Matrix{<:Number} = rotation tensor
"""
function rotation_parameters_WM(R::Matrix{<:Number})

    # Get quaternion
    q = quaternion_from_rotation_tensor(R)

    # Vector part of quaternion
    e = q[2:4]

    # Angle of rotation
    ϕ = 2*asin(norm(e))
    
    # Bauchau's ν parameter for Wiener-Milenkovic parametrization
    ν = cos(ϕ/4)^2

    # Wiener-Milenkovic rotation parameters
    p = 2/ν*e
    
    return p
end


"""
quaternion_from_rotation_tensor(R::Matrix{<:Number})

Computes the quaternion (Euler parameters) given a rotation tensor

# Arguments
- R::Matrix{<:Number} = rotation tensor
"""
function quaternion_from_rotation_tensor(R::Matrix{<:Number})

    # Quaternion component product matrix
    T = [1+R[1,1]+R[2,2]+R[3,3] R[3,2]-R[2,3] R[1,3]-R[3,1] R[2,1]-R[1,2];
         R[3,2]-R[2,3] 1+R[1,1]-R[2,2]-R[3,3] R[1,2]+R[2,1] R[1,3]+R[3,1];
         R[1,3]-R[3,1] R[2,1]+R[1,2] 1-R[1,1]+R[2,2]-R[3,3] R[2,3]+R[3,2];
         R[2,1]-R[1,2] R[1,3]+R[3,1] R[2,3]+R[3,2] 1-R[1,1]-R[2,2]+R[3,3]]
    
    # Maximum value and position
    i = argmax(tr(T))
    Tii = T[i]
    
    # Initialize quaternion
    q = zeros(4)

    # Largest quaternion component
    q[i] = 1/2*sqrt(Tii)

    # Other components
    ind = [1,2,3,4]
    popat!(ind,i)
    q[ind] = T[i,ind]/(4*q[i])

    return q
    
end