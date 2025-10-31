using AeroBeams

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor (for the structure)
λ = 1

# Altitude [m]
h = 20e3

# Gravity 
g = 9.80665

# Airspeed [m/s]
U = 30

# Wing and stabilizers parasite drag
wingCd0 = stabsCd0 = 0

# Twist curvature
k1 = 0.002

# Discretization
nElemWing = 30

# System solver
σ0 = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=1.0)

# Model
conventionalHALEclampedSteady,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,g=g,altitude=h,airspeed=U,k1=k1,nElemWing=nElemWing,wingCd0=wingCd0,stabsCd0=stabsCd0,includeVS=false)

# Set clamped BC
clamp = create_BC(name="clamp",beam=conventionalHALEclampedSteady.beams[2],node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
push!(conventionalHALEclampedSteady.BCs,clamp)
update_model!(conventionalHALEclampedSteady)

# Create and solve steady problem
problem = create_SteadyProblem(model=conventionalHALEclampedSteady,systemSolver=NR)
solve!(problem)

# Undeformed nodal positions
x1_0_l = vcat([vcat(conventionalHALEclampedSteady.beams[1].elements[e].r_n1[1],conventionalHALEclampedSteady.beams[1].elements[e].r_n2[1]) for e in 1:div(nElemWing,2)]...)
x1_0_r = vcat([vcat(conventionalHALEclampedSteady.beams[2].elements[e].r_n1[1],conventionalHALEclampedSteady.beams[2].elements[e].r_n2[1]) for e in 1:div(nElemWing,2)]...)
x1_0 = vcat(x1_0_l,x1_0_r)
x3_0_l = vcat([vcat(conventionalHALEclampedSteady.beams[1].elements[e].r_n1[3],conventionalHALEclampedSteady.beams[1].elements[e].r_n2[3]) for e in 1:div(nElemWing,2)]...)
x3_0_r = vcat([vcat(conventionalHALEclampedSteady.beams[2].elements[e].r_n1[3],conventionalHALEclampedSteady.beams[2].elements[e].r_n2[3]) for e in 1:div(nElemWing,2)]...)
x3_0 = vcat(x3_0_l,x3_0_r)

# Undeformed elemental positions
x1_e_l = reverse(-[conventionalHALEclampedSteady.beams[1].elements[e].x1 for e in 1:div(nElemWing,2)])
x1_e_r = [conventionalHALEclampedSteady.beams[2].elements[e].x1 for e in 1:div(nElemWing,2)]
x1_e = vcat(x1_e_l,x1_e_r)

# Extract outputs
# Displacements over span
u1_of_x1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElemWing]...)
u3_of_x1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElemWing]...)
# Deformed nodal positions
x1_def = x1_0 .+ u1_of_x1
x3_def = x3_0 .+ u3_of_x1
# Angle of attack over span
α_of_x1 = [problem.aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:nElemWing]

println("Finished conventionalHALEclampedSteady.jl") #src