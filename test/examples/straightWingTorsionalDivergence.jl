using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = Indicial()

# Airspeed range
URange = collect(0:0.25:30)

# Pitch angle (the divergence speed can be found even with θ=0)
θ = 0.1*π/180

# Discretization
nElem = 20

# Number of non-oscillatory modes of interest
nNOModes = 3

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
d = 2*normCGPos - 1
e = c * (normSparPos - 1/4)
e2 = -(d-a)*c/2
surf = create_AeroSurface(solver=aeroSolver,airfoil=airfoil,c=c,normSparPos=normSparPos,updateAirfoilParameters=false,smallAngles=true)

# Wing
L = 1
GJ = 10
ρA = 1
∞ = 1e16
wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,GJ=GJ)],I=[inertia_matrix(ρA=ρA,e2=e2)],aeroSurface=surf,rotationParametrization="E321",p0=[0;0;θ])

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
straightWingTorsionalDivergence = create_Model(name="straightWingTorsionalDivergence",beams=[wing],BCs=[clamp],atmosphere=atmosphere)

# Elemental and nodal arclength coordinate
x1_e = [straightWingTorsionalDivergence.elements[e].x1 for e in 1:nElem]
x1_n = vcat([vcat(straightWingTorsionalDivergence.elements[e].r_n1[1],straightWingTorsionalDivergence.elements[e].r_n2[1]) for e in 1:nElem]...)

# Pre-allocate memory and initialize output arrays
problem = Array{EigenProblem}(undef,length(URange))
dampingsNonOscillatory = Array{Vector{Float64}}(undef,length(URange))
AoA = Array{Vector{Float64}}(undef,length(URange))
tipAoA = Array{Float64}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update airspeed on model
    set_motion_basis_A!(model=straightWingTorsionalDivergence,v_A=[0;U;0])
    # Create and solve eigenproblem
    problem[i] = create_EigenProblem(model=straightWingTorsionalDivergence,nModes=1)
    solve!(problem[i])
    # Damping of non-oscillatory modes
    dampingsNonOscillatory[i] = problem[i].dampingsNonOscillatory[1:nNOModes]
    # Tip angle of attack
    AoA[i] = [problem[i].model.elements[e].aero.flowAnglesAndRates.αₑ*180/π for e in 1:nElem]
    tipAoA[i] = AoA[i][end]
end

# Get divergence speed, dynamic pressure and λ
UDvec,_ = find_non_oscillatory_instability(URange,dampingsNonOscillatory)
UD = minimum(UDvec)
qD = 1/2*ρ*UD^2
λD = sqrt(qD*c*2π*e/GJ)
iD = findlast(x -> x<UD, URange)

# Analytical angle of attack and twist - Eq. 4.52 of Hodges and Pierce
x = 0:0.1:L
q = 1/2*ρ*URange.^2
λ = sqrt.(q*c*2π*e/GJ)
αRef = @. tan(λD*L)*sin(λD*x)+cos(λD*x) - 1
if αRef[end] < 0
    αRef ./= minimum(αRef)
else
    αRef ./= maximum(αRef)
end
αTipRef = @. 180/π*θ*(sec(λ*L)-1) + 180/π*θ
θTipRef = @. 180/π*θ*(sec(λ*L)-1)

# Analytical divergence speed - Eq. 4.55 of Hodges and Pierce
qDRef = GJ/(e*c*2π)*(π/(2L))^2
UDRef = sqrt(2*qDRef/ρ)

# Relative error in divergence speed
ϵrel = UD/UDRef - 1
println("UDRef = $UDRef, UD = $UD, ϵrel = $ϵrel")

println("Finished straightWingTorsionalDivergence.jl")