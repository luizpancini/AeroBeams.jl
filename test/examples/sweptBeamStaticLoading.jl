using AeroBeams, LinearInterpolations

# Sweep angle range
ΛRange = π/180*collect(0:10:40)

# Discretization
nElem = 20

# Beam properties
L = 1
EIy = 1
GJ = 1
ρA = 1
∞ = 1e16

# Tip force magnitude
F = 1

# Pre-allocate memory and initialize output arrays
problem = Array{SteadyProblem}(undef,length(ΛRange))
tip_u3 = Array{Float64}(undef,length(ΛRange))
tip_twist = Array{Float64}(undef,length(ΛRange))

# Sweep angle of sweep
for (i,Λ) in enumerate(ΛRange)
    println("Solving for Λ = $(round(Λ*180/π)) deg")
    # Beam
    beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,GJ=GJ/cos(Λ),EIy=EIy*cos(Λ)^2)],I=[inertia_matrix(ρA=ρA)],rotationParametrization="E321",p0=[-Λ;0;0])
    # beam.S[1][4,5] = beam.S[1][5,4] = EIy*GJ*sin(Λ)*cos(Λ)^2
    # BCs
    clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
    tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,types=["F3A"],values=[F])
    # Model
    model = create_Model(name="beamModel",beams=[beam],BCs=[clamp,tipForce])
    # Create and solve eigenproblem
    problem[i] = create_SteadyProblem(model=model)
    solve!(problem[i])
    # Outputs
    tip_u3[i] = problem[i].nodalStatesOverσ[end][nElem].u_n2[3]
    tip_p = problem[i].nodalStatesOverσ[end][nElem].p_n2
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*[0; 1; 0]
    tip_twist[i] = asind(Δ[3])
end

println("Finished sweptBeamStaticLoading.jl")