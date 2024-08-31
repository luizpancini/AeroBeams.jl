using AeroBeams

# Aerodynamic solver and derivatives method
aeroSolver = Indicial()
derivationMethod = AD()

# Altitude
h = 20e3

# Gravity
g = 9.80665

# Pitch angle
θ = 2*π/180

# Discretization
nElem = 32

# Model
SMWSteady,L = create_SMW(aeroSolver=aeroSolver,derivationMethod=derivationMethod,θ=θ,nElem=nElem,altitude=h,g=g)

# Set system solver options (limit initial load factor)
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Undeformed nodal and midpoint positions
x1_0 = vcat([vcat(SMWSteady.beams[1].elements[e].r_n1[1],SMWSteady.beams[1].elements[e].r_n2[1]) for e in 1:nElem]...)
x3_0 = vcat([vcat(SMWSteady.beams[1].elements[e].r_n1[3],SMWSteady.beams[1].elements[e].r_n2[3]) for e in 1:nElem]...)

# Set airspeed range, and initialize outputs
URange = collect(0:1:30)
x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))
tip_u3 = Array{Float64}(undef,length(URange))
tip_twist = Array{Float64}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Update velocity of basis A (and update model)
    set_motion_basis_A!(model=SMWSteady,v_A=[0;U;0])
    # Create and solve problem
    global problem = create_SteadyProblem(model=SMWSteady,systemSolver=NR)
    solve!(problem)
    # Displacements over span
    u1_of_x1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
    u3_of_x1 = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
    # Deformed nodal positions
    x1_def[i] = x1_0 .+ u1_of_x1
    x3_def[i] = x3_0 .+ u3_of_x1
    # Tip OOP displacement
    tip_u3[i] = problem.nodalStatesOverσ[end][nElem].u_n2[3]
    # Tip twist
    tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*[0; 1; 0]
    tip_twist[i] = asind(Δ[3])
end

println("Finished SMWSteady.jl")