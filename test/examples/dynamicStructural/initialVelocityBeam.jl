using AeroBeams

# Select initial conditions as either :displacement, :rotation or :both
initialConditions = :both

# Select boundary conditions as either "ss-ss" (simple support on both sides) or "ss-roller" (simple support on the left and roller on the right)
BCType = "ss-roller"

# Initial conditions: sinusoidal velocity V3 and/or angular velocity Ω2
σ = 1e-1
udot3 = x1 -> σ*sin.(2π*x1/L)
θdot2 = x1 -> -σ*2π/L*cos.(2π*x1/L)
pdot2 = x1 -> sec.(0/4).^2 .* θdot2(x1)

# Beam
L = 1
EI = 1
ρA = 1
∞ = 1e4
nElem = 40
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,EIy=EI)],I=[inertia_matrix(ρA=ρA)])

# Set initial conditions
if initialConditions == :displacement
    add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:udot0_of_x1 => x1->[0; 0; udot3(x1)]))
elseif initialConditions == :rotation
    add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:pdot0_of_x1 => x1->[0; pdot2(x1); 0]))
elseif initialConditions == :both
    add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:udot0_of_x1 => x1->[0; 0; udot3(x1)], :pdot0_of_x1 => x1->[0; pdot2(x1); 0]))
else
    error("initialConditions must be either :displacement, :rotation or :both")
end

# BCs
ss1 = create_BC(name="simple-support-1",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
ss2 = create_BC(name="simple-support-2",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
roller = create_BC(name="roller",beam=beam,node=nElem+1,conditions=Dict(:u3A=>0))
if BCType == "ss-ss"
    bcs = [ss1,ss2]
elseif BCType == "ss-roller"
    bcs = [ss1,roller]
end

# Model
initialVelocityBeam = create_Model(name="initialVelocityBeam",beams=[beam],BCs=bcs)

# Time and frequency variables
ω2 = (2π/L)^2 * sqrt(EI/ρA)
T = 2π/ω2
cycles = 1
tf = cycles * T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=initialVelocityBeam,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
Nt = length(t)
tNorm = t/T
midElem = div(nElem,2)
quarterElem = div(nElem,4)
u3_quarter = [states[i][quarterElem].nodal.u2[3] for i in 1:Nt]
V3_quarter = [(states[i][quarterElem].midpoint.V[3]+states[i][quarterElem+1].midpoint.V[3])/2 for i in 1:Nt]
Vdot3_quarter = [(states[i][quarterElem].rates.Vdot[3]+states[i][quarterElem+1].rates.Vdot[3])/2 for i in 1:Nt]
θ2_root = [states[i][1].nodal.θ[2] for i in 1:Nt]
Ω2_mid = [(states[i][midElem].midpoint.Ω[2]+states[i][midElem+1].midpoint.Ω[2])/2 for i in 1:Nt]
Ωdot2_mid = [(states[i][midElem].rates.Ωdot[2]+states[i][midElem+1].rates.Ωdot[2])/2 for i in 1:Nt]

# Compute analytical values
u3_quarter_analytic = σ/ω2*sin.(ω2*t)*sin(2π*1/4)
V3_quarter_analytic = σ*cos.(ω2*t)*sin(2π*1/4)
Vdot3_quarter_analytic = -σ*ω2*sin.(ω2*t)*sin(2π*1/4)
θ2_root_analytic = -σ*2π/L/ω2*sin.(ω2*t)*cos(2π*0)
Ω2_mid_analytic = -σ*2π/L*cos.(ω2*t)*cos(2π*1/2)
Ωdot2_mid_analytic = σ*2π/L*ω2*sin.(ω2*t)*cos(2π*1/2)

## Plots
# Set paths
relPath = "/test/outputs/figures/initialVelocityBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,scale=1/δ/10,scalePos=[0.15;-0.05;0],timeStampPos=[0.5;0.05;0],plotFrequency=1,plotLimits=([0,L],[-L/2,L/2],[-L/3,L/3]),save=true,savePath=string(relPath,"/initialVelocityBeam_deformation.gif"),displayProgress=true)
display(anim)

# Plot configurations
lw = 2
ms = 5
gr()

# Normalized displacement at quarter-length
plt1 = plot(xlabel="\$t/T\$", ylabel="\$u_3/\\delta\$ at \$x_1=L/4\$")
plot!(tNorm,u3_quarter/σ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],u3_quarter_analytic[1:5:end]/σ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt1)
savefig(string(absPath,"/initialVelocityBeam_disp.pdf"))

# Normalized velocity at quarter-length
plt2 = plot(xlabel="\$t/T\$", ylabel="\$V_3/\\delta\$ at \$x_1=L/4\$ [\$1\$/s]")
plot!(tNorm,V3_quarter/σ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],V3_quarter_analytic[1:5:end]/σ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt2)
savefig(string(absPath,"/initialVelocityBeam_vel.pdf"))

# Normalized acceleration at quarter-length
plt3 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3/\\delta\$ at \$x_1=L/4\$ [\$1\$/\$s^2\$]")
plot!(tNorm,Vdot3_quarter/σ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Vdot3_quarter_analytic[1:5:end]/σ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt3)
savefig(string(absPath,"/initialVelocityBeam_acc.pdf"))

# Normalized rotation at root
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\theta/(2\\pi\\delta)\$ at \$x_1=0\$")
plot!(tNorm,θ2_root/(2π)/σ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],θ2_root_analytic[1:5:end]/(2π)/σ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt4)
savefig(string(absPath,"/initialVelocityBeam_rot.pdf"))

# Angular velocity at mid-length
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_2\$ at \$x_1=L/2\$ [rad/s]")
plot!(tNorm,Ω2_mid, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ω2_mid_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt5)
savefig(string(absPath,"/initialVelocityBeam_angVel.pdf"))

# Angular acceleration at mid-length
plt6 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_2\$ at \$x_1=L/2\$ [rad/\$s^2\$]")
plot!(tNorm,Ωdot2_mid, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ωdot2_mid_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt6)
savefig(string(absPath,"/initialVelocityBeam_angAcc.pdf"))

println("Finished initialVelocityBeam.jl")