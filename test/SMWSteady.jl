using AeroBeams, LinearAlgebra, LinearInterpolations, Plots, ColorSchemes, DelimitedFiles

# Wing surface
chord = 1.0
normSparPos = 0.5
aeroSolver = Indicial()
derivationMethod = AD()
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=flatPlate,c=chord,normSparPos=normSparPos)

# Wing beam
θ = π/180*2
L = 16
GJ,EIy,EIz = 1e4,2e4,4e6
ρA,ρIs = 0.75,0.1
ρIy = ρIs*EIy/EIz
ρIz = ρIs-ρIy
nElem = 32
∞ = 1e12
wing = create_Beam(name="wing",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs)],rotationParametrization="E321",p0=[0;0;θ],aeroSurface=surf)

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
g = 9.80665
h = 20e3
SMWSteady = create_Model(name="SMWSteady",beams=[wing],BCs=[clamp],gravityVector=[0;0;-g],altitude=h)

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
    problem = create_SteadyProblem(model=SMWSteady,systemSolver=NR)
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

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 3
# Normalized deformed wingspan
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$ [% semispan]", xlims=[0,1], ylims=[-20,60])
for (i,U) in enumerate(URange)
    plot!(x1_def[i]/L, x3_def[i]/L*100, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/SMWSteady_1.pdf"))
# Tip OOP disp vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp [% semispan]", xlims=[0,30], ylims=[-20,60])
plot!(URange, tip_u3/L*100, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/SMWSteady_2.pdf"))
# Tip twist vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,30], ylims=[-1,3])
plot!(URange, tip_twist, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/SMWSteady_3.pdf"))

println("Finished SMWSteady.jl")