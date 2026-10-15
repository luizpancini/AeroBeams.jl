using AeroBeams

# Stiffness factor for the aircraft's wing
λ = 1e0

# Altitude [m]
h = 0e3

# Discretization
nElem = 20

# Landing gear height
LGheight = 1

# Set bending curvature range
k2Range = range(-0.015,0.045,5)

# System solver
systemSolver = create_NLsolver()

# Initialize outputs
x1_def = Array{Vector{Float64}}(undef,length(k2Range))
x3_def = Array{Vector{Float64}}(undef,length(k2Range))

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Tip roller
    ρ = 1/k2
    θ = 16/ρ
    z = abs(k2) > 0 ? ρ*(1-cos(θ)) : 0
    dummyBeam = create_Beam(length=1,nElements=nElem,S=[AeroBeams.I6])
    roller = create_BC(name="roller",beam=dummyBeam,node=nElem+1,conditions=Dict(:u3A => z-LGheight))
    # Initialize model
    SMW,L = create_SMW(nElem=nElem,stiffnessFactor=λ,altitude=h,k2=k2,additionalBCs=[roller])
    # Undeformed nodal positions
    x1_0 = reduce(vcat, 
                [SMW.assemblyData.elements[e].assemblyData.r_n1[1], 
                SMW.assemblyData.elements[e].assemblyData.r_n2[1]] 
            for e in 1:nElem)
    x3_0 = reduce(vcat, 
                [SMW.assemblyData.elements[e].assemblyData.r_n1[3], 
                SMW.assemblyData.elements[e].assemblyData.r_n2[3]] 
            for e in 1:nElem)
    # Create and solve steady problem
    problem = create_SteadyProblem(model=SMW,systemSolver=systemSolver)
    solve!(problem)
    # Final states
    states = problem.loadstepHistory.states[end]
    # Displacements over span
    u1_of_x1 = reduce(vcat,[states[e].nodal.u1[1],
                            states[e].nodal.u2[1]] for e in 1:nElem)
    u3_of_x1 = reduce(vcat,[states[e].nodal.u1[3],
                            states[e].nodal.u2[3]] for e in 1:nElem)
    # Deformed nodal positions
    x1_def[i] = x1_0 + u1_of_x1
    x3_def[i] = x3_0 + u3_of_x1
end

## Plots
# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_on_ground"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold])
lw = 2
labels = ["\$k_2 = $(k2) \$" for k2 in k2Range]
gr()

# Normalized deformed span
plt1 = plot(xlabel="Normalized spanwise direction", ylabel="OOP direction [m]", xlims=[0,1], ylims=[-LGheight,Inf64])
for (i,k2) in enumerate(k2Range)
    plot!(x1_def[i]/L, x3_def[i], c=colors[i], lw=lw, label=labels[i])
end
display(plt1)
savefig(string(absPath,"/cHALEwing_on_ground.pdf"))

println("Finished cHALEwing_on_ground.jl")