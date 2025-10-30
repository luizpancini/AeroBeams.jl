using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = Indicial()

# Number of elements for the wing
nElem = 2

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)
ρ = atmosphere.ρ

# Aerodynamic surface
airfoil = deepcopy(flatPlate)
c = 0.2
normSparPos = 1/2
normCGPos = 1/2
a = 2*normSparPos - 1
e = 2*normCGPos - 1
e2 = -(e-a)*c/2
surf = create_AeroSurface(solver=aeroSolver,airfoil=airfoil,c=c,normSparPos=normSparPos,updateAirfoilParameters=false,smallAngles=true)

# Rigid wing
L = 1
GJ = 30
ρA = 1
∞ = 1e15
wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=ρA,e2=e2)],aeroSurface=surf)

# Attachment spring in pitch
kα = GJ/L
spring = create_Spring(elementsIDs=[1],nodesSides=[2],kp=kα*[1;0;0])
add_springs_to_beam!(beam=wing,springs=[spring])

# BCs
journal1 = create_BC(name="journal-1",beam=wing,node=1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])
journal2 = create_BC(name="journal-2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])

# Model
typicalSectionDivergence = create_Model(name="typicalSectionDivergence",beams=[wing],BCs=[journal1,journal2],atmosphere=atmosphere)

# Airspeed range
URange = collect(.1:.1:40)

# Number of non-oscillatory modes of interest
nNOModes = 2

# Pre-allocate memory and initialize output arrays
problem = Array{EigenProblem}(undef,length(URange))
dampingsNonOscillatory = Array{Vector{Float64}}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update airspeed on model
    set_motion_basis_A!(model=typicalSectionDivergence,v_A=[0;U;0])
    # Create and solve eigenproblem
    problem[i] = create_EigenProblem(model=typicalSectionDivergence,nModes=1)
    solve!(problem[i])
    # Damping of non-oscillatory modes
    dampingsNonOscillatory[i] = problem[i].dampingsNonOscillatory[1:nNOModes]
end

# Get divergence speed
UDvec,_ = find_non_oscillatory_instability(URange,dampingsNonOscillatory)
UD = minimum(UDvec)

# Analytical divergence speed - Eq. 4.8 of Hodges and Pierce
UDRef = sqrt(kα/(1/2*ρ*c*L*2π*c*(normSparPos-1/4)))

# Relative error in divergence speed
ϵrel = UD/UDRef - 1
println("Relative error in divergence speed = $ϵrel")

println("Finished typicalSectionDivergence.jl")