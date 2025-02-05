# First vector of basis A, resolved in that basis
const a1 = [1.0; 0.0; 0.0]


# Second vector of basis A, resolved in that basis
const a2 = [0.0; 1.0; 0.0]


# Third vector of basis A, resolved in that basis
const a3 = [0.0; 0.0; 1.0]


# 3x3 identity matrix
const I3 = Matrix(1.0*LinearAlgebra.I,3,3)


# 6x6 identity matrix
const I6 = Matrix(1.0*LinearAlgebra.I,6,6)


"""
    round_off!(x)

Rounds the array or number to input tolerance (defaults to machine epsilon)

# Arguments
- `x`: array or number
- `tol`: tolerance
"""
function round_off!(x,tol=eps())

    if x isa Number
        x = round(x/tol)*tol
    else
        for i in eachindex(x)
            x[i] = round(x[i]/tol)*tol
        end
    end

    return x
end
export round_off!


"""
    rms(x)

Computes the root mean square of an array

"""
function rms(x)
    return sqrt(mean(x .^ 2))
end
export rms


"""
    tilde(v)

Computes the skew-symmetric matrix associated with a three-element vector

# Arguments
- `v`: vector
"""
function tilde(v)
    return [  0.0  -v[3]   v[2]
             v[3]    0.0  -v[1]
            -v[2]   v[1]   0.0]
end
export tilde


"""
    axial(R)

Computes the axial part of a 3x3 matrix

# Arguments
- `R`: matrix
"""
function axial(R)
    return 1/2*[R[3,2]-R[2,3]; R[1,3]-R[3,1]; R[2,1]-R[1,2]]
end
export axial


# Divides the input variables in-place
function divide_inplace!(divisor, vars...)
    return (var ./ divisor for var in vars)
end


# Multiplies the input variables in-place
function multiply_inplace!(multiplier, vars...)
    return (var .* multiplier for var in vars)
end


# Computes the scalar product of a third-order tensor represented by matrices A1, A2, and A3 with the vector b
function mul3(A1,A2,A3,b)
    
    return hcat(A1*b, A2*b, A3*b)
    
end

# Gets the ith component of the 3x3 identity matrix
function a_i(i)

    @assert i in [1,2,3]

    if i==1
        return a1
    elseif i==2
        return a2
    elseif i==3
        return a3
    end

end


"""
    isotropic_stiffness_matrix(; kwargs...)

Creates a 6x6 sectional stiffness matrix for a cross-section made of isotropic material

# Arguments
- `∞`: value for rigid properties
- `EA`: axial stiffness
- `GAy`: shear stiffness in the x2 direction
- `GAz`: shear stiffness in the x3 direction
- `GJ`: torsional stiffness
- `EIy`: bending stiffness in the x2 direction
- `EIz`: bending stiffness in the x3 direction
- `t2`: offset from reference line to tension center in the x2 direction
- `t3`: offset from reference line to tension center in the x3 direction
- `s2`: offset from reference line to shear center in the x2 direction
- `s3`: offset from reference line to shear center in the x3 direction
"""
function isotropic_stiffness_matrix(; ∞=1e16,EA=∞,GAy=∞,GAz=∞,GJ=∞,EIy=∞,EIz=∞,t2=0,t3=0,s2=0,s3=0)

    # Validate
    @assert all(x->x>0,[EA,GAy,GAz,GJ,EIy,EIz])

    # See Hodges' book Eqs. 4.114 - 4.118
    z = [0 t3 -t2; -s3 0 0; s2 0 0]

    A = diagm([EA,GAy,GAz])

    Tinv = diagm([GJ,EIy,EIz])

    B = A*z

    D = Tinv + z'*A*z

    S = [A B; B' D]

    return S

end
export isotropic_stiffness_matrix


"""
    inertia_matrix(; kwargs...)

Creates a 6x6 sectional inertia matrix

# Arguments
- `ρA`: mass per unit length
- `ρIy`: mass moment of inertia per unit length about the x2-axis
- `ρIz`: mass moment of inertia per unit length about the x3-axis
- `ρIs`: mass moment of inertia per unit length about the x1-axis
- `ρIyz`: mass product of inertia per unit length
- `e2`: offset from reference line to center of gravity in the x2 direction
- `e3`: offset from reference line to center of gravity in the x3 direction
"""
function inertia_matrix(; ρA=0,ρIy=0,ρIz=0,ρIs=ρIy+ρIz,ρIyz=0,e2=0,e3=0)

    # Validate
    @assert all(x->x>=0,[ρA,ρIy,ρIz,ρIs])

    ηtilde = tilde([0;ρA*e2;ρA*e3])
    I = diagm([ρA,ρA,ρA,ρIs,ρIy,ρIz])
    I[5,6] = I[6,5] = ρIyz
    I[1:3,4:6] = -ηtilde
    I[4:6,1:3] = ηtilde

    return I

end
export inertia_matrix


# Computes some useful functions of the curvature vector
function curvature_quantities(k)

    # External self-product 
    kkT = k * k'

    # Skew-symmetric matrix 
    ktilde = tilde(k)

    # Norm 
    knorm = sqrt(dot(k,k))

    return kkT, ktilde, knorm
end


# Computes the position vector at an arclength value
function position_vector_from_curvature(R0,k,x1) 

    kkT, ktilde, knorm = curvature_quantities(k)

    # Set according to |k|
    if knorm > 0
        return R0 * (1.0/knorm*(I3 .- kkT/knorm^2)*sin(knorm*x1) .+ ktilde/knorm^2*(1.0-cos(knorm*x1)) .+ kkT/knorm^2*x1) * a1
    else
        return x1*R0[:,1]
    end
end


# Computes the rotation tensor at an arclength position
function rotation_tensor_from_curvature(R0,k,x1) 

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
    rotation_tensor_E321(p)

Computes the rotation tensor according to Euler parameters sequence 3-2-1

# Arguments
- `p`: rotation parameters
"""
function rotation_tensor_E321(p)

    # Validate
    @assert length(p) == 3

    # Euler angles 3-2-1 sequence (yaw, pitch and roll angles), respective sines and cosines
    yaw,pitch,roll = p[1],p[2],p[3]
    sy,cy = sincos(yaw)
    sp,cp = sincos(pitch)
    sr,cr = sincos(roll)

    # Rotation tensor that brings the reference basis to the final basis
    R = [cp*cy cy*sr*sp-cr*sy sr*sy+cr*cy*sp;
         cp*sy cy*cr+sr*sp*sy cr*sp*sy-cy*sr;
           -sp          cp*sr          cp*cr]

    round_off!(R)       

    return R

end
export rotation_tensor_E321


"""
    rotation_tensor_E313(p)

Computes the rotation tensor according to Euler parameters sequence 3-1-3

# Arguments
- `p`: rotation parameters
"""
function rotation_tensor_E313(p)

    # Validate
    @assert length(p) == 3

    # Euler angles 3-1-3 sequence (precession=ϕ, nutation=θ and spin=ψ angles), respective sines and cosines
    ϕ,θ,ψ = p[1],p[2],p[3]
    sϕ,cϕ = sincos(ϕ)
    sθ,cθ = sincos(θ)
    sψ,cψ = sincos(ψ)
       
    # Rotation tensor that brings the reference basis to the final basis
    R = [-sϕ*cθ*sψ+cϕ*cψ -sϕ*cθ*cψ-cϕ*sψ  sψ*sθ;
          cϕ*cθ*sψ+sϕ*cψ  cϕ*cθ*cψ-sϕ*sψ -cϕ*sθ;
                   sθ*sψ           sθ*cψ     cθ]

    round_off!(R)

    return R

end
export rotation_tensor_E313


"""
    rotation_tensor_WM(p)

Computes the rotation tensor according to Wiener-Milenkovic parameters

# Arguments
- `p`: rotation parameters
"""
function rotation_tensor_WM(p)

    # Scaling factor and scaled rotation parameters
    λ,ps,_,pNorm,psNorm = rotation_parameter_scaling(p)

    # Shorthand for components of scaled rotation parameters 
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
        
    # Rotation tensor
    R = υ² * Θ

    return R,Θ,pNorm,λ,ps,ps1,ps2,ps3,ps0,υ,υ²,ps1s,ps2s,ps3s,ps1ps2,ps2ps3,ps1ps3

end
export rotation_tensor_WM


# Scales the Wiener-Milenkovic rotation parameters
function rotation_parameter_scaling(p)

    # Initialize scaling factor and number of odd half rotations
    λ = 1.0
    halfRotations = 0
    
    # Norm of the extended parameters vector
    pNorm = norm(p)
    
    # Scale according to norm
    if pNorm > 4
        halfRotations = round(pNorm/8)
        λ = 1 - 8 * halfRotations / pNorm
    end
    
    # Scaled parameters
    ps = λ * p
    
    # Norm of scaled parameters vector
    psNorm = abs(λ) * pNorm

    return λ,ps,halfRotations,pNorm,psNorm
    
end


"""
    scaled_rotation_parameters(p)

Returns the scaled (extended) Wiener-Milenkovic rotation parameters

# Arguments
- `p`: rotation parameters
"""
function scaled_rotation_parameters(p)

    # Perform scaling
    _,ps,_ = rotation_parameter_scaling(p)

    return ps
    
end
export scaled_rotation_parameters


# Computes the derivatives of the rotation tensor with respect to the scaled rotation parameters
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
    R_ps1 = υ²_ps1*Θ + υ²*Θ_ps1
    R_ps2 = υ²_ps2*Θ + υ²*Θ_ps2
    R_ps3 = υ²_ps3*Θ + υ²*Θ_ps3

    return R_ps1,R_ps2,R_ps3,υ²_ps1,υ²_ps2,υ²_ps3,Θ_ps1,Θ_ps2,Θ_ps3
    
end


# Computes the derivatives of the scaling factor with respect to the extended rotation parameters
function scaling_derivatives_extended_parameters(λ,p,pNorm)

    λ_p = λ == 1.0 ? zeros(3) : (1-λ)*p/pNorm^2

    return λ_p

end


# Computes the derivatives of the rotation tensor with respect to the extended rotation parameters
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
    R_p1 = R_ps1*ps1_p1 + R_ps2*ps2_p1 + R_ps3*ps3_p1
    R_p2 = R_ps1*ps1_p2 + R_ps2*ps2_p2 + R_ps3*ps3_p2
    R_p3 = R_ps1*ps1_p3 + R_ps2*ps2_p3 + R_ps3*ps3_p3

    return R_p1,R_p2,R_p3,R_ps1,R_ps2,R_ps3,υ²_ps1,υ²_ps2,υ²_ps3,Θ_ps1,Θ_ps2,Θ_ps3,ps_p,ps1_p1,ps2_p1,ps3_p1,ps1_p2,ps2_p2,ps3_p2,ps1_p3,ps2_p3,ps3_p3
end

function rotation_tensor_derivatives_extended_parameters(p)

    _,Θ,pNorm,λ,ps,ps1,ps2,ps3,ps0,υ,υ²,ps1s,ps2s,ps3s,ps1ps2,ps2ps3,ps1ps3 = rotation_tensor_WM(p)

    # Rotation tensor derivatives w.r.t. scaled rotation parameters
    R_ps1,R_ps2,R_ps3,_ = rotation_tensor_derivatives_scaled_parameters(ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²)

    # Scaling derivatives w.r.t. extended rotation parameters
    λ_p = scaling_derivatives_extended_parameters(λ,p,pNorm)

    # Scaled rotation parameters derivatives w.r.t. extended rotation parameters 
    ps_p = λ_p*p' + λ*I3
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
    R_p1 = R_ps1*ps1_p1 + R_ps2*ps2_p1 + R_ps3*ps3_p1
    R_p2 = R_ps1*ps1_p2 + R_ps2*ps2_p2 + R_ps3*ps3_p2
    R_p3 = R_ps1*ps1_p3 + R_ps2*ps2_p3 + R_ps3*ps3_p3

    return R_p1,R_p2,R_p3
end


# Computes the time derivative of the rotation tensor 
function rotation_tensor_time_derivative(R_ps1,R_ps2,R_ps3,ps_p,pdot)

    # Scaled rotation parameters time derivative
    psdot = ps_p*pdot
    ps1dot = psdot[1]
    ps2dot = psdot[2]
    ps3dot = psdot[3]

    # Rotation tensor time derivative 
    Rdot = R_ps1*ps1dot + R_ps2*ps2dot + R_ps3*ps3dot

    return Rdot,ps1dot,ps2dot,ps3dot

end

function rotation_tensor_time_derivative(p,pdot)

    # Rotation parameters' variables
    _,Θ,pNorm,λ,ps,ps1,ps2,ps3,ps0,υ,υ²,ps1s,ps2s,ps3s,ps1ps2,ps2ps3,ps1ps3 = rotation_tensor_WM(p)

    # Rotation tensor and rotation parameters' derivatives w.r.t. extended rotation parameters
    _,_,_,R_ps1,R_ps2,R_ps3,_,_,_,_,_,_,ps_p,_,_,_,_,_,_,_,_,_ = rotation_tensor_derivatives_extended_parameters(p,pNorm,λ,ps,ps0,ps1,ps2,ps3,ps1s,ps2s,ps3s,ps1ps2,ps1ps3,ps2ps3,Θ,υ,υ²)

    # Rotation tensor time derivative 
    Rdot,_ = rotation_tensor_time_derivative(R_ps1,R_ps2,R_ps3,ps_p,pdot)

    return Rdot
end


# Computes the derivatives of the time derivative of the rotation tensor with respect to the extended rotation parameters 
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
    R_ps1ps1 = υ²_psps[1,1]*Θ + 2*υ²_ps1*Θ_ps1 + υ²*Θ_ps1ps1
    R_ps2ps2 = υ²_psps[2,2]*Θ + 2*υ²_ps2*Θ_ps2 + υ²*Θ_ps2ps2
    R_ps3ps3 = υ²_psps[3,3]*Θ + 2*υ²_ps3*Θ_ps3 + υ²*Θ_ps3ps3
    R_ps2ps1 = υ²_psps[1,2]*Θ + υ²_ps2*Θ_ps1 + υ²_ps1*Θ_ps2 + υ²*Θ_ps1ps2
    R_ps3ps1 = υ²_psps[1,3]*Θ + υ²_ps3*Θ_ps1 + υ²_ps1*Θ_ps3 + υ²*Θ_ps1ps3
    R_ps3ps2 = υ²_psps[2,3]*Θ + υ²_ps3*Θ_ps2 + υ²_ps2*Θ_ps3 + υ²*Θ_ps2ps3
    R_ps1ps2 = R_ps2ps1
    R_ps1ps3 = R_ps3ps1
    R_ps2ps3 = R_ps3ps2

    # Rotation tensor time derivatives' derivative w.r.t. scaled rotation parameters
    Rdot_ps1 = R_ps1ps1*ps1dot + R_ps2ps1*ps2dot + R_ps3ps1*ps3dot
    Rdot_ps2 = R_ps1ps2*ps1dot + R_ps2ps2*ps2dot + R_ps3ps2*ps3dot
    Rdot_ps3 = R_ps1ps3*ps1dot + R_ps2ps3*ps2dot + R_ps3ps3*ps3dot

    # Rotation tensor time derivatives' derivative w.r.t. extended rotation parameters
    Rdot_p1 = Rdot_ps1*ps1_p1 + Rdot_ps2*ps2_p1 + Rdot_ps3*ps3_p1
    Rdot_p2 = Rdot_ps1*ps1_p2 + Rdot_ps2*ps2_p2 + Rdot_ps3*ps3_p2
    Rdot_p3 = Rdot_ps1*ps1_p3 + Rdot_ps2*ps2_p3 + Rdot_ps3*ps3_p3

    return Rdot_p1,Rdot_p2,Rdot_p3

end


# Computes the transpose of tangent operator tensor according to Wiener-Milenkovic parameters
function tangent_operator_transpose_WM(p)

    # Scaling factor and scaled rotation parameters
    _,ps,_,_,psNorm = rotation_parameter_scaling(p)

    # Useful functions of the scaled rotation parameters 
    ps0 = 2 - psNorm^2/8
    υ² = (4-ps0)^-2

    # Tangent operator transpose
    return 2*υ²*(ps0*I3 + 1/4*ps*(ps') - tilde(ps))

end
export tangent_operator_transpose_WM

function tangent_operator_transpose_WM(ps,ps0,υ²)

    return 2*υ²*(ps0*I3 + 1/4*ps*(ps') - tilde(ps))

end


# Computes the inverse of the transpose of the tangent operator tensor according to Wiener-Milenkovic parameters
function tangent_operator_transpose_inverse_WM(p)

    # Scaling factor and scaled rotation parameters
    _,ps,_,_,psNorm = rotation_parameter_scaling(p)

    # Useful functions of the scaled rotation parameters 
    ps0 = 2 - psNorm^2/8

    # Tangent operator inverse transpose
    return 1/2*(ps0*I3 + 1/4*ps*(ps') + tilde(ps))

end

function tangent_operator_transpose_inverse_WM(ps,ps0)

    return 1/2*(ps0*I3 + 1/4*ps*(ps') + tilde(ps))

end


# Computes the derivatives of the tangent tensor's transpose and its inverse with respect to the extended rotation parameters
function tangent_tensor_functions_derivatives_extended_parameters(HT,ps1,ps2,ps3,υ²,υ²_ps1,υ²_ps2,υ²_ps3,ps_p)

    # Tangent operator transpose derivatives w.r.t. scaled rotation parameters
    HT_ps1 = υ²_ps1*HT/υ² + υ²*1/2*[ps1 ps2 ps3; ps2 -ps1 4; ps3 -4 -ps1]
    HT_ps2 = υ²_ps2*HT/υ² + υ²*1/2*[-ps2 ps1 -4; ps1 ps2 ps3; 4 ps3 -ps2]
    HT_ps3 = υ²_ps3*HT/υ² + υ²*1/2*[-ps3 4 ps1; -4 -ps3 ps2; ps1 ps2 ps3]

    # Tangent operator transpose inverse derivatives w.r.t. scaled rotation parameters
    HTinv_ps1 = 1/8*[ps1 ps2 ps3; ps2 -ps1 -4; ps3 4 -ps1]
    HTinv_ps2 = 1/8*[-ps2 ps1 4; ps1 ps2 ps3; -4 ps3 -ps2]
    HTinv_ps3 = 1/8*[-ps3 -4 ps1; 4 -ps3 ps2; ps1 ps2 ps3]

    # Tangent operator transpose derivatives w.r.t. extended rotation parameters
    HT_p1 = HT_ps1*ps_p[1,1] + HT_ps2*ps_p[2,1] + HT_ps3*ps_p[3,1]
    HT_p2 = HT_ps1*ps_p[1,2] + HT_ps2*ps_p[2,2] + HT_ps3*ps_p[3,2]
    HT_p3 = HT_ps1*ps_p[1,3] + HT_ps2*ps_p[2,3] + HT_ps3*ps_p[3,3]

    # Tangent operator transpose inverse derivatives w.r.t. extended rotation parameters
    HTinv_p1 = HTinv_ps1*ps_p[1,1] + HTinv_ps2*ps_p[2,1] + HTinv_ps3*ps_p[3,1]
    HTinv_p2 = HTinv_ps1*ps_p[1,2] + HTinv_ps2*ps_p[2,2] + HTinv_ps3*ps_p[3,2]
    HTinv_p3 = HTinv_ps1*ps_p[1,3] + HTinv_ps2*ps_p[2,3] + HTinv_ps3*ps_p[3,3]

    return HT_p1,HT_p2,HT_p3,HTinv_p1,HTinv_p2,HTinv_p3

end


# Computes the derivatives of the tangent tensor's transpose with respect to the extended rotation parameters
function tangent_tensor_transpose_derivatives_extended_parameters(p)

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
    HT_ps1 = υ²_ps1*HT/υ² + υ²*1/2*[ps1 ps2 ps3; ps2 -ps1 4; ps3 -4 -ps1]
    HT_ps2 = υ²_ps2*HT/υ² + υ²*1/2*[-ps2 ps1 -4; ps1 ps2 ps3; 4 ps3 -ps2]
    HT_ps3 = υ²_ps3*HT/υ² + υ²*1/2*[-ps3 4 ps1; -4 -ps3 ps2; ps1 ps2 ps3]

    # Tangent operator transpose derivatives w.r.t. extended rotation parameters
    HT_p1 = HT_ps1*ps_p[1,1] + HT_ps2*ps_p[2,1] + HT_ps3*ps_p[3,1]
    HT_p2 = HT_ps1*ps_p[1,2] + HT_ps2*ps_p[2,2] + HT_ps3*ps_p[3,2]
    HT_p3 = HT_ps1*ps_p[1,3] + HT_ps2*ps_p[2,3] + HT_ps3*ps_p[3,3]

    return HT_p1,HT_p2,HT_p3

end
export tangent_tensor_transpose_derivatives_extended_parameters


# Computes the appropriate force scaling for the linear system of equations
function force_scaling(S)

    # Create array with compliance entries
    complianceEntries = vcat(vcat(S...)...)
            
    # Get all nonzero entries
    complianceNonzeroIndices = findall( x -> x != 0.0, complianceEntries)
    
    # Set force scaling based on number of nonzero compliance matrix entries
    forceScaling = 1.0
    if !isempty(complianceNonzeroIndices)
        
        nonzeroComplianceEntries = complianceEntries[complianceNonzeroIndices]

        α = length(nonzeroComplianceEntries)/sum(nonzeroComplianceEntries)/100

        forceScaling = α > 2 ? nextpow(2,α) : 2/nextpow(2,1/α)
    end
    
    return forceScaling
end


"""
    rotation_angle(p)

Computes the rotation angle given the Wiener-Milenkovic rotation parameters

# Arguments
- `p`: rotation parameters
"""
function rotation_angle(p)

    # Scale rotation parameters and get number of half rotations
    _,_,halfRotations,pNorm,_ = rotation_parameter_scaling(p)

    # Highest absolute value component
    greatestComp = argmax(abs.(p))

    # Sign of rotation angle
    signal = sign(p[greatestComp])

    return signal * (4*atan((pNorm-8*halfRotations)/4) + 2π*halfRotations)

end
export rotation_angle


"""
    rotation_angle_limited(p)

Computes the rotation angle (in the range -360 to 360 degrees) given the Wiener-Milenkovic rotation parameters

# Arguments
- `p`: rotation parameters
"""
function rotation_angle_limited(p)

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
- `R`: rotation tensor
"""
function rotation_parameters_WM(R)

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
export rotation_parameters_WM


"""
    rotation_parameters_Rodrigues(R)

Computes the Rodrigues rotation parameters given a rotation tensor

# Arguments
- `R`: rotation tensor
"""
function rotation_parameters_Rodrigues(R)

    # Get quaternion
    q = quaternion_from_rotation_tensor(R)

    # Vector part of quaternion
    e = q[2:4]

    # Angle of rotation
    ϕ = 2*asin(norm(e))
    
    # Bauchau's ν parameter for Rodrigues parametrization
    ν = cos(ϕ/2)

    # Rodrigues rotation parameters
    p = 2/ν*e
    
    return p
end
export rotation_parameters_Rodrigues


"""
    ypr_from_rotation_tensor(R; ϵround=1e-10,ϵsingularity=1e-3,assumeNullYawInSingularity=true)

Computes the Euler angles from the sequence 3-2-1 (yaw, pitch, roll) given the rotation tensor

# Arguments
- `R`: rotation tensor

# Keyword arguments
- `ϵround`: tolerance for rounding off elements of R to zero
- `ϵsingularity`: tolerance to consider singular case 
- `assumeNullYawInSingularity`: flag to assume zero yaw angle in singularity
"""
function ypr_from_rotation_tensor(R; ϵround=1e-10,ϵsingularity=1e-3,assumeNullYawInSingularity=true)

    # Round-off
    round_off!(R,ϵround)
    
    # Pitch angle
    pitch = asin(-R[3,1])

    # Check if near singularity (pitch ≈ ±π/2)
    if abs(cos(pitch)) < ϵsingularity
        # Singular case - only sum of yaw and roll can be determined: set angles according to assumption of either one being zero
        yaw_plus_roll = atan(R[1,2]/R[1,3])
        if assumeNullYawInSingularity
            yaw = 0
            roll = yaw_plus_roll
        else
            roll = 0
            yaw = yaw_plus_roll
        end
    else
        # Non-singular case
        yaw = atan(R[2,1]/R[1,1])
        roll = atan(R[3,2]/R[3,3])
    end

    return yaw, pitch, roll

end
export ypr_from_rotation_tensor


"""
    quaternion_from_rotation_tensor(R)

Computes the quaternion (Euler parameters) given a rotation tensor. Derived from Bauchau's book section 13.3.4. 

# Arguments
- `R`: rotation tensor
"""
function quaternion_from_rotation_tensor(R)

    # Initialize quaternion
    q = zeros(eltype(R),4)
    
    # Quaternion component product matrix (T = [1+tr(R) 2*axial(R)'; 2*axial(R) (1-tr(R))*I3+2*Symmetric(R)], but using this form yields innacuracies when computing derivatives with ForwardDiff)
    T = [1+R[1,1]+R[2,2]+R[3,3] R[3,2]-R[2,3] R[1,3]-R[3,1] R[2,1]-R[1,2];
         R[3,2]-R[2,3] 1+R[1,1]-R[2,2]-R[3,3] R[1,2]+R[2,1] R[1,3]+R[3,1];
         R[1,3]-R[3,1] R[2,1]+R[1,2] 1-R[1,1]+R[2,2]-R[3,3] R[2,3]+R[3,2];
         R[2,1]-R[1,2] R[1,3]+R[3,1] R[2,3]+R[3,2] 1-R[1,1]-R[2,2]+R[3,3]]
    
    # Maximum value along diagonal and position
    m = argmax(diag(T))
    Tmm = T[m,m]

    # Largest quaternion component
    q[m] = 1/2*sqrt(Tmm)

    # Other components
    ind = [1,2,3,4]
    popat!(ind,m)
    q[ind] = T[m,ind]/(4*q[m])

    return q 
end
export quaternion_from_rotation_tensor


"""
    WM_to_ypr(p)

Transforms Wiener-Milenkovic parameters to Euler parameters of sequence 3-2-1 (yaw, pitch, roll)

# Arguments
- `p`: rotation parameters
"""
function WM_to_ypr(p)

    R,_ = rotation_tensor_WM(p)

    yaw,pitch,roll = ypr_from_rotation_tensor(R)
    
    return [yaw; pitch; roll]
end
export WM_to_ypr


"""
    ypr_to_WM(p)

Transforms Euler parameters of sequence 3-2-1 (yaw, pitch, roll) to Wiener-Milenkovic parameters

# Arguments
- `p`: rotation parameters
"""
function ypr_to_WM(p)

    return rotation_parameters_WM(rotation_tensor_E321(p))
end
export ypr_to_WM


"""
    rotation_between_WM(p1,p2)

Computes the Wiener-Milenkovic parameters describing the rotation from p1 to p2, i.e., p12 such that R2(p2) = R12(p12)*R1(p1)

# Arguments
- `p1`: initial rotation parameters
- `p2`: final rotation parameters
"""
function rotation_between_WM(p1,p2)

    # Rotation tensors
    R1,_ = rotation_tensor_WM(p1)
    R2,_ = rotation_tensor_WM(p2)
    R12 = R2*R1'

    # Rotation parameters vector
    p12 = rotation_parameters_WM(R12)

    # Check consistency: true rotation parameters might actually be -p12 (since rotations of ϕ about axis n or -ϕ about axis -n are equal)
    R12_check1,_ = rotation_tensor_WM(p12)
    R12_check2,_ = rotation_tensor_WM(-p12)
    check1 = norm(R12_check1.-R12)
    check2 = norm(R12_check2.-R12)

    # Select correct rotation parameters vector (the one for which norm(R12.-R12_check) = 0)
    if check2 < check1
        p12 = -p12
    end

    return p12
end
export rotation_between_WM


"""
    mode_tracking(controlParam::Vector{<:Real},freqs::Array{Vector{Float64}},damps::Array{Vector{Float64}},eigenvectors::Array{Matrix{ComplexF64}})

Applies mode tracking based on eigenvectors match

# Arguments
- `controlParam::Vector{<:Real}`: vector of control parameter
- `freqs::Array{Vector{Float64}}`: frequencies vector
- `damps::Array{Vector{Float64}}`: dampings vector
- `eigenvectors::Array{Matrix{ComplexF64}}`: complex-valued eigenvectors
"""
function mode_tracking(controlParam::Vector{<:Real},freqs::Array{Vector{Float64}},damps::Array{Vector{Float64}},eigenvectors::Array{Matrix{ComplexF64}})

    # Validate inputs
    @assert length(controlParam) == length(freqs) == length(damps) == length(eigenvectors)
    for i in eachindex(freqs)
        @assert length(freqs[i]) == length(damps[i]) == size(eigenvectors[i],2)
    end

    # Length of control parameter range
    N = length(controlParam)
 
    # Number of modes
    nModes = length(freqs[1])

    # Initialize array of matched modes
    matchedModes = Array{Vector{Int64}}(undef,N)
    matchedModes[1] = collect(1:nModes)

    # Loop over control parameter
    for c = 2:N
        # Skip if eigenvectors are undefined
        if all(isnan,eigenvectors[c]) || all(isnan,eigenvectors[c-1])
            continue
        end
        # Mode matching index matrix: eigenvector alignment
        I = [abs(dot(eigenvectors[c-1][:,m1],eigenvectors[c][:,m2])) for m1 in 1:nModes, m2 in 1:nModes]
        # Set matched modes by highest index
        matchedModes[c] = highest_in_rowcol(I)
        # Skip if matched modes are in order
        if matchedModes[c] == collect(1:nModes)
            continue
        end
        # Set frequencies, dampings and eigenvectors in new order of matched modes
        freqs[c:end] = [freqs[i][matchedModes[c]] for i=c:N]
        damps[c:end] = [damps[i][matchedModes[c]] for i=c:N]
        eigenvectors[c:end] = [eigenvectors[i][:,matchedModes[c]] for i=c:N]
    end

    return freqs,damps,eigenvectors,matchedModes
end
export mode_tracking


# Finds the unrepeated columns of the n highest values in matrix O, where n is the size of O
function highest_in_rowcol(O)

    # Size of the original matrix
    n = size(O, 1)

    # Initialize output vector
    v = fill(0,n)

    # Array of leftover rows and columns of original matrix
    rowO = collect(1:n)
    colO = collect(1:n)

    # Copy of original matrix
    C = deepcopy(O)

    # Loop
    for _ in 1:n
        # Cartesian index of highest value in copy matrix
        ind = argmax(C)
        # Corresponding row and column
        rowC = ind[1]
        colC = ind[2]
        # ith entry is the corresponding column of the original matrix 
        v[rowO[rowC]] = colO[colC]
        # Remove row and column from copy matrix
        C = C[setdiff(1:end, colC), setdiff(1:end, colC)]
        # Remove row and column from leftover arrays
        popat!(rowO,rowC)
        popat!(colO,colC)
    end

    return v
end


"""
    get_FFT_and_PSD(t::Vector{<:Real},y::Vector{<:Real}; tol::Float=1e3*eps())

Computes the FFT and PSD of signal y(t)

# Arguments
- `t::Vector{<:Real}`: time signal
- `y::Vector{<:Real}`: quantity signal

# Keyword arguments
- `tol::AbstractFloat`: tolerance for time signal being equally spaced
"""
function get_FFT_and_PSD(t::Vector{<:Real},y::Vector{<:Real}; tol::AbstractFloat=1e3*eps())
    
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
    yFFT = fft(y)
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
    Newton_solver(f::Function, x0::AbstractArray; absTol::Real=1e-9, relTol::Real=1e-9, maxIter::Int=50)

Solves a nonlinear algebraic system of equations \$( f(x) = 0 \$) using the Newton-Raphson method.

# Arguments
- `f::Function`: the nonlinear system of equations to solve. Must return a vector of residuals for a given `x`.
- `x0::AbstractArray`: initial guess for the solution.

# Keyword arguments
- `absTol::Real`: absolute tolerance for convergence (default: `1e-9`).
- `relTol::Real`: relative tolerance for convergence (default: `1e-9`).
- `maxIter::Int`: maximum number of iterations allowed (default: `50`).

# Returns
- `x::AbstractArray`: the solution vector, if convergence is achieved.
- `converged::Bool`: `true` if the solution converged within the given tolerances, `false` otherwise.

"""
function Newton_solver(f::Function, x0::AbstractArray; absTol::Real=1e-9, relTol::Real=1e-9, maxIter::Int=50)

    # Validate inputs
    @assert length(x0) == length(f(x0)) "The size of x0 must match the output dimension of f(x)"

    # Define the Jacobian function using ForwardDiff
    J(x) = ForwardDiff.jacobian(f, x)

    # Initialize variables
    x = deepcopy(x0)               # Current solution estimate
    converged = false              # Convergence flag
    iter = 0                       # Iteration counter
    ϵabs = ϵrel = Inf              # Initial error values

    # Main Newton-Raphson iteration loop
    while !converged && iter < maxIter
        # Update number of iterations
        iter += 1
        # Compute residual and Jacobian
        res = f(x)
        jac = J(x)
        # Check if Jacobian is singular
        if det(jac) ≈ 0
            error("Jacobian matrix is singular at iteration $iter.")
        end
        # Solve for the Newton step
        Δx = -jac \ res
        # Update solution
        x += Δx
        # Compute convergence metrics
        ϵabs = norm(res, Inf)  # Absolute residual norm
        ϵrel = norm(Δx ./ x, Inf)  # Relative solution change
        # Check convergence
        converged = ϵabs < absTol || ϵrel < relTol
    end

    # Return solution and convergence status
    return x, converged
end
export Newton_solver