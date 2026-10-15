using AeroBeams

# Airspeed range [m/s]
URange = vcat(1e-2,1:1:40)

# Altitude
h = 20e3

# Pitch angle [deg]
θ = deg2rad(2)

# Discretization
nElem = 32

# Aerodynamic solver
aeroSolver = Indicial()

# Model
SMWSteady,L = create_SMW(aeroSolver=aeroSolver,pitchAngle=θ,nElem=nElem,altitude=h)

# Undeformed nodal positions
x1_0 = reduce(vcat, 
             [SMWSteady.assemblyData.elements[e].assemblyData.r_n1[1], 
              SMWSteady.assemblyData.elements[e].assemblyData.r_n2[1]] 
        for e in 1:nElem)
x3_0 = reduce(vcat, 
             [SMWSteady.assemblyData.elements[e].assemblyData.r_n1[3], 
              SMWSteady.assemblyData.elements[e].assemblyData.r_n2[3]] 
        for e in 1:nElem)

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(URange))
x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))
tip_u3 = Array{Float64}(undef,length(URange))
tip_twist = Array{Float64}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Update velocity of basis A
    update_v_A!(SMWSteady; v_A=[0,U,0])
    # Set initial guess solution for this airspeed
    x0 = i == 1 ? zeros(0) : problem[i-1].systemArrays.x
    # Create and solve problem
    problem[i] = create_SteadyProblem(model=SMWSteady,x0=x0)
    solve!(problem[i])
    # Solution history and final states
    hist = problem[i].loadstepHistory
    states = hist.states[end]
    tipNodeFinalState = states[nElem].nodal
    tipElemFinalAero = hist.aeroVariables[end][nElem]
    # Displacements over span
    u1_of_x1 = reduce(vcat,[states[e].nodal.u1[1],
                            states[e].nodal.u2[1]] for e in 1:nElem)
    u3_of_x1 = reduce(vcat,[states[e].nodal.u1[3],
                            states[e].nodal.u2[3]] for e in 1:nElem)
    # Deformed nodal positions
    x1_def[i] = x1_0 .+ u1_of_x1
    x3_def[i] = x3_0 .+ u3_of_x1
    # Tip OOP displacement
    tip_u3[i] = tipNodeFinalState.u2[3]
    # Tip twist
    tip_twist[i] = twist_angle(tipNodeFinalState.p2_b)
end

## Plots
# Set paths
relPath = "/test/outputs/figures/SMWSteady"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot deformed shape
deformationPlot = plot_steady_deformation(problem[end],save=true,savePath=string(relPath,"/SMWSteady_deformation.pdf"))
display(deformationPlot)

# Plot configurations
gr()
lw = 2
ms = 3

# Normalized deformed wingspan
plt1 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$ [% semispan]", xlims=[0,1], ylims=[-20,70])
for (i,U) in enumerate(URange)
    plot!(x1_def[i]/L, x3_def[i]/L*100, lz=U, c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [m/s]")
end
display(plt1)
savefig(string(absPath,"/SMWSteady_disp.pdf"))

# Tip OOP disp vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp [% semispan]", xlims=[0,30], ylims=[-20,60])
plot!(URange, tip_u3/L*100, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/SMWSteady_OOP.pdf"))

# Tip twist vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,30], ylims=[-1,3])
plot!(URange, tip_twist, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/SMWSteady_twist.pdf"))

println("Finished SMWSteady.jl")