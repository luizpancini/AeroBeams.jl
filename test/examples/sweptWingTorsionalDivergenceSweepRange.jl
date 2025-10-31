using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = Indicial()

# Sweep angle range
ΛRange = π/180*vcat(-60:5:60)

# Pitch angle (the divergence speed can be found even with θ=0)
θ = 0*π/180

# Discretization
nElem = 20

# Number of non-oscillatory modes of interest
nNOModes = 5

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)
ρ = atmosphere.ρ

# Aerodynamic surface properties
airfoil = deepcopy(flatPlate)
c = 0.2
normSparPos = 1/2
normCGPos = 1/2
a = 2*normSparPos - 1
d = 2*normCGPos - 1
e = c * (normSparPos - 1/4)
e2 = -(d-a)*c/2

# Wing spar properties
L = 1
GJ = 10
ρA = 1
∞ = 1e16

# Airspeed range (considering minimum and maximum analytical divergence speeds)
Umin = sqrt(2*(π/2)^2*GJ/(ρ*e*c*2π*L^2*cos(minimum(abs.(ΛRange)))^2))
Umax = sqrt(2*(π/2)^2*GJ/(ρ*e*c*2π*L^2*cos(maximum(abs.(ΛRange)))^2))
URange = collect(floor(Umin-5):0.5:ceil(Umax+1))

# Pre-allocate memory and initialize output arrays
dampingsNonOscillatory = Array{Vector{Float64}}(undef,length(ΛRange),length(URange))
modeDampings = Array{Vector{Float64}}(undef,length(ΛRange),nNOModes)
modeDampingsEst = Array{Vector{Float64}}(undef,length(ΛRange),nNOModes)
UD = Array{Float64}(undef,length(ΛRange))

# Sweep angle of sweep
for (i,Λ) in enumerate(ΛRange)
    # Aero surface
    surf = create_AeroSurface(solver=aeroSolver,airfoil=airfoil,c=c,normSparPos=normSparPos,Λ=Λ,updateAirfoilParameters=false,smallAngles=true)
    # Wing
    wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,GJ=GJ)],I=[inertia_matrix(ρA=ρA,e2=e2)],aeroSurface=surf,rotationParametrization="E321",p0=[-Λ;0;θ])
    # BCs
    clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
    # Model
    model = create_Model(name="wingModel",beams=[wing],BCs=[clamp],atmosphere=atmosphere)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for Λ = $(round(Λ*180/π)) deg, U = $U m/s")
        # Update airspeed on model
        set_motion_basis_A!(model=model,v_A=[0;U;0])
        # Create and solve eigenproblem
        problem = create_EigenProblem(model=model,nModes=1)
        solve!(problem)
        # Damping of non-oscillatory modes
        dampingsNonOscillatory[i,j] = problem.dampingsNonOscillatory[1:nNOModes]
    end
    # Get divergence speed
    UDvec,_ = find_non_oscillatory_instability(URange,dampingsNonOscillatory[i,:])
    UD[i] = minimum(UDvec)
    println("UD = $(UD[i]) m/s at Λ = $(round(Λ*180/π)) deg")
end

# Approximate divergence speed - Eq. 4.103 of Hodges and Pierce
UDRef = @. sqrt(2*(π/2)^2*GJ/(ρ*e*c*2π*L^2*cos(ΛRange)^2))

println("Finished sweptWingTorsionalDivergenceSweepRange.jl")