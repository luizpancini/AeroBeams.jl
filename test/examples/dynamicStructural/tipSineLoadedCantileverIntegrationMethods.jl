using AeroBeams

# Beam
L,b,H = 1.0, 0.1, 0.1
E,ρ = 200e9, 7.8e3
A,Iy,Iz, = b*H, b*H^3/12, H*b^3/12
∞ = 1e12
nElements = 20
beam = create_Beam(name="beam",length=L,nElements=nElements,S=[isotropic_stiffness_matrix(EA=E*A,EIy=E*Iy,∞=∞)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*Iy,ρIz=ρ*Iz)])

# BCs
F = 1e3
ω = 200
force = create_BC(name="force",beam=beam,node=nElements+1,conditions=Dict(:F3A=>t->F*sin(ω*t)))
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
tipSineLoadedCantileverIntegrationMethods = create_Model(name="tipSineLoadedCantileverIntegrationMethods",beams=[beam],BCs=[clamp,force])

# Time variables
T = 2π/ω
cycles = 20
tf = cycles * T
Δt = T/100

# Time integration methods: trapezoidal rule (Newmark's method with γ = 1/2), BDF2 and generalized-α
integrationMethods = [:newmark, :bdf2, GeneralizedAlpha(ρ∞=0.8)]
methodLabels = ["Trapezoidal", "BDF2", "Generalized-α"]

# Initialize outputs
t = Vector{Vector{Float64}}(undef,length(integrationMethods))
u3_tip = Vector{Vector{Float64}}(undef,length(integrationMethods))
M2_root = Vector{Vector{Float64}}(undef,length(integrationMethods))

# Loop integration methods
for (i,integrationMethod) in enumerate(integrationMethods)
    # Create and solve the problem
    problem = create_DynamicProblem(model=tipSineLoadedCantileverIntegrationMethods,finalTime=tf,Δt=Δt,integrationMethod=integrationMethod)
    solve!(problem)
    # Unpack numerical solution
    sol = problem.solutionHistory
    t[i] = sol.t
    u3_tip[i] = [sol.states[j][end].nodal.u2[3] for j in 1:length(sol.t)]
    M2_root[i] = [sol.states[j][1].nodal.M1[2] for j in 1:length(sol.t)]
end

## Plots
# Set paths
relPath = "/test/outputs/figures/tipSineLoadedCantileverIntegrationMethods"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(integrationMethods), categorical=true)
gr()
lw = 1

# Tip u3
plt1 = plot(xlabel="\$t\$ [s]", ylabel="Tip \$u_3\$ [mm]")
for i in eachindex(integrationMethods)
    plot!(t[i], u3_tip[i]*1e3, c=colors[i], lw=lw, label=methodLabels[i])
end
display(plt1)
savefig(string(absPath,"/tipSineLoadedCantileverIntegrationMethods_u3.pdf"))

# Root M2
plt2 = plot(xlabel="\$t\$ [s]", ylabel="Root \$M_2^*\$ [Nm]")
for i in eachindex(integrationMethods)
    plot!(t[i], M2_root[i], c=colors[i], lw=lw, label=methodLabels[i])
end
display(plt2)
savefig(string(absPath,"/tipSineLoadedCantileverIntegrationMethods_M2.pdf"))

println("Finished tipSineLoadedCantileverIntegrationMethods.jl")