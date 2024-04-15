using AeroBeams, LinearAlgebra, Plots

# Select initial conditions as either "displacement", "rotation" or "both"
initialConditions = "displacement"

# Select boundary conditions as either "ss-ss" (simple support on both sides) or "ss-roller" (simple support on the left and roller on the right)
BCType = "ss-roller"

# Initial conditions: sinusoidal velocity V₃ and/or angular velocity Ω₂
σ = 1e-1
udot₃ = x1 -> σ*sin.(2*π*x1/L)
θdot₂ = x1 -> -σ*2*π/L*cos.(2*π*x1/L)
pdot₂ = x1 -> sec.(0/4).^2 .* θdot₂(x1)

# Beam
L = 1.0
EIy = 1.0
ρA = 1.0
Φ = 1e4
nElem = 48
stiffnessMatrix = diagm([Φ,Φ,Φ,Φ,EIy,Φ])
inertiaMatrix = diagm([ρA,ρA,ρA,0,0,0])
beam = Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])
if initialConditions == "displacement"
    beam.udot0_of_x1=x1->[0; 0; udot₃(x1)]
elseif initialConditions == "rotation"
    beam.pdot0_of_x1=x1->[0; pdot₂(x1); 0]
elseif initialConditions == "both"
    beam.udot0_of_x1=x1->[0; 0; udot₃(x1)]
    beam.pdot0_of_x1=x1->[0; pdot₂(x1); 0]
else
    error("Wrong initialConditions")
end
update_beam!(beam)

# BCs
ss1 = create_BC(name="simple-support-1",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
ss2 = create_BC(name="simple-support-2",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p3A"],values=[0,0,0,0,0])
roller = create_BC(name="roller",beam=beam,node=nElem+1,types=["u3A"],values=[0])
if BCType == "ss-ss"
    bcs = [ss1,ss2]
elseif BCType == "ss-roller"
    bcs = [ss1,roller]
end

# Model
initialVelocityBeam = Model(name="initialVelocityBeam",beams=[beam],BCs=bcs)

# plot_undeformed_assembly(initialVelocityBeam)

# Time and frequency variables
ω₂ = (2*π/L)^2*sqrt(EIy/ρA)
T = 2*π/ω₂
cycles = 1
tf = cycles*T
Δt = 2.5e-4

# Initial velocities update options
if initialConditions == "displacement"
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-4, displayProgress=true, relaxFactor = 0.5, Δt = 2.5e-6)
elseif initialConditions == "rotation"
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=100,tol=1e-4, displayProgress=true)
elseif initialConditions == "both"
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-4, displayProgress=true, relaxFactor = 0.5, Δt = 2.5e-6)
end

# Create and solve the problem
problem = DynamicProblem(model=initialVelocityBeam,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tNorm = t/T
u₃_quarter = [problem.nodalStatesOverTime[i][div(nElem,4)].u_n2[3] for i in 1:length(tNorm)]
V₃_quarter = [(problem.elementalStatesOverTime[i][div(nElem,4)].V[3]+problem.elementalStatesOverTime[i][div(nElem,4)+1].V[3])/2.0 for i in 1:length(tNorm)]
Vdot₃_quarter = [(problem.elementalStatesRatesOverTime[i][div(nElem,4)].Vdot[3]+problem.elementalStatesRatesOverTime[i][div(nElem,4)+1].Vdot[3])/2.0 for i in 1:length(tNorm)]
θ₂_root = [problem.nodalStatesOverTime[i][1].θ_n1 for i in 1:length(tNorm)]
Ω₂_mid = [(problem.elementalStatesOverTime[i][div(nElem,2)].Ω[2]+problem.elementalStatesOverTime[i][div(nElem,2)+1].Ω[2])/2.0 for i in 1:length(tNorm)]
Ωdot₂_mid = [(problem.elementalStatesRatesOverTime[i][div(nElem,2)].Ωdot[2]+problem.elementalStatesRatesOverTime[i][div(nElem,2)+1].Ωdot[2])/2.0 for i in 1:length(tNorm)]

# Compute analytical values
u₃_quarter_analytic = σ/ω₂*sin.(ω₂*t)*sin(2*π*1/4)
V₃_quarter_analytic = σ*cos.(ω₂*t)*sin(2*π*1/4)
Vdot₃_quarter_analytic = -σ*ω₂*sin.(ω₂*t)*sin(2*π*1/4)
θ₂_root_analytic = -σ*2*π/L/ω₂*sin.(ω₂*t)*cos(2*π*0)
Ω₂_mid_analytic = -σ*2*π/L*cos.(ω₂*t)*cos(2*π*1/2)
Ωdot₂_mid_analytic = σ*2*π/L*ω₂*sin.(ω₂*t)*cos(2*π*1/2)

# Plots
# ------------------------------------------------------------------------------
# # Normalized displacement at quarter-length
# Normalized displacement at quarter-length
plt1 = Plots.plot()
Plots.plot!(tNorm,u₃_quarter/σ, c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$u_3/\\delta\$ at \$x_1=L/4\$", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],u₃_quarter_analytic[1:20:end]/σ, c=:blue, markersize=3, label="Analytical", show=true)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/initialVelocityBeam_1.pdf"))
# Normalized velocity at quarter-length
plt2 = Plots.plot()
Plots.plot!(tNorm,V₃_quarter/σ, c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$V_3/\\delta\$ at \$x_1=L/4\$ [\$1\$/s]", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],V₃_quarter_analytic[1:20:end]/σ, c=:blue, markersize=3, label="Analytical", show=true)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/initialVelocityBeam_2.pdf"))
# Normalized acceleration at quarter-length
plt3 = Plots.plot()
Plots.plot!(tNorm,Vdot₃_quarter/σ, c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$\\dot{V}_3/\\delta\$ at \$x_1=L/4\$ [\$1\$/\$s^2\$]", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],Vdot₃_quarter_analytic[1:20:end]/σ, c=:blue, markersize=3, label="Analytical", show=true)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/initialVelocityBeam_3.pdf"))
# Normalized rotation at root
plt4 = Plots.plot()
Plots.plot!(tNorm,θ₂_root/(2*π)/σ, c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$\\theta/(2\\pi\\delta)\$ at \$x_1=0\$", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],θ₂_root_analytic[1:20:end]/(2*π)/σ, c=:blue, markersize=3, label="Analytical", show=true)
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/initialVelocityBeam_4.pdf"))
# Angular velocity at mid-length
plt5 = Plots.plot()
Plots.plot!(tNorm,Ω₂_mid, c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$\\Omega_2\$ at \$x_1=L/2\$ [rad/s]", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],Ω₂_mid_analytic[1:20:end], c=:blue, markersize=3, label="Analytical", show=true)
display(plt5)
savefig(string(pwd(),"/test/outputs/figures/initialVelocityBeam_5.pdf"))
# Angular acceleration at mid-length
plt6 = Plots.plot()
Plots.plot!(tNorm,Ωdot₂_mid, c=:black, linewidth=2, xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_2\$ at \$x_1=L/2\$ [rad/\$s^2\$]", label="Numerical", show=true)
Plots.scatter!(tNorm[1:20:end],Ωdot₂_mid_analytic[1:20:end], c=:blue, markersize=3, label="Analytical", show=true)
display(plt6)
savefig(string(pwd(),"/test/outputs/figures/initialVelocityBeam_6.pdf"))

println("Finished initialVelocityBeam.jl")