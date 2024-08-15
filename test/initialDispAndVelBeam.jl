using AeroBeams, LinearAlgebra, Plots
 
# Select initial conditions as either "displacement", "rotation" or "both"
initialConditions = "displacement"

# Select boundary conditions as either "ss-ss" (simple support on both sides) or "ss-roller" (simple support on the left and roller on the right)
BCType = "ss-roller"

# Initial conditions: sinusoidal displacement of u₃ and/or rotation of θ₂
δ = 1e-3
σ = 1e-1
u₃ = x1 -> δ*sin.(2*π*x1/L)
θ₂ = x1 -> -δ*2*π/L*cos.(2*π*x1/L)
p₂ = x1 -> 4*tan.(θ₂(x1)/4)
udot₃ = x1 -> σ*sin.(2*π*x1/L)
θdot₂ = x1 -> -σ*2*π/L*cos.(2*π*x1/L)
pdot₂ = x1 -> sec.(0/4).^2 .* θdot₂(x1)

# Beam
L = 1
EIy = 1
ρA = 1
Φ = 1e4
nElem = 48
stiffnessMatrix = diagm([Φ,Φ,Φ,Φ,EIy,Φ])
inertiaMatrix = diagm([ρA,ρA,ρA,0,0,0])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])
if initialConditions == "displacement"
    beam.u0_of_x1=x1->[0; 0; u₃(x1)]
    beam.udot0_of_x1=x1->[0; 0; udot₃(x1)]
elseif initialConditions == "rotation"
    beam.p0_of_x1=x1->[0; p₂(x1); 0]
    beam.pdot0_of_x1=x1->[0; pdot₂(x1); 0]
elseif initialConditions == "both"
    beam.u0_of_x1=x1->[0; 0; u₃(x1)]
    beam.p0_of_x1=x1->[0; p₂(x1); 0]
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
initialDispAndVelBeam = create_Model(name="initialDispAndVelBeam",beams=[beam],BCs=bcs)

# Time and frequency variables
ω₂ = (2*π/L)^2*sqrt(EIy/ρA)
T = 2*π/ω₂
cycles = 1
tf = cycles*T
Δt = T/100

# Initial velocities update options
if initialConditions == "displacement"
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,displayProgress=true, relaxFactor=0.5, Δt=Δt/1e3)
elseif initialConditions == "rotation"
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2, tol=1e-8, displayProgress=true, relaxFactor=0.5, Δt=Δt/1e3)
elseif initialConditions == "both"
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2, displayProgress=true, relaxFactor=0.5, Δt=Δt/1e3)
end

# Create and solve the problem
problem = create_DynamicProblem(model=initialDispAndVelBeam,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
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
u₃_quarter_analytic = δ*cos.(ω₂*t)*sin(2*π*1/4) + σ/ω₂*sin.(ω₂*t)*sin(2*π*1/4)
V₃_quarter_analytic = -δ*ω₂*sin.(ω₂*t)*sin(2*π*1/4) + σ*cos.(ω₂*t)*sin(2*π*1/4)
Vdot₃_quarter_analytic = -δ*ω₂^2*cos.(ω₂*t)*sin(2*π*1/4) + -σ*ω₂*sin.(ω₂*t)*sin(2*π*1/4)
θ₂_root_analytic = -δ*2*π/L*cos.(ω₂*t)*cos(2*π*0) + -σ*2*π/L/ω₂*sin.(ω₂*t)*cos(2*π*0)
Ω₂_mid_analytic = δ*2*π/L*ω₂*sin.(ω₂*t)*cos(2*π*1/2) + -σ*2*π/L*cos.(ω₂*t)*cos(2*π*1/2)
Ωdot₂_mid_analytic = δ*2*π/L*ω₂^2*cos.(ω₂*t)*cos(2*π*1/2) + σ*2*π/L*ω₂*sin.(ω₂*t)*cos(2*π*1/2)

# Plots
# --------------------------------------------------------------------------------
lw = 2
ms = 5
relPath = "/test/outputs/figures/initialDispAndVelBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,scale=1/δ/10,scalePos=[0.15;-0.05;0],timeStampPos=[0.5;-0.05;0],plotFrequency=1,plotLimits=[(0,L),(-L/2,L/2),(-L,L)],save=true,savePath=string(relPath,"/initialDispAndVelBeam_deformation.gif"),displayProgress=true)
# Displacement at quarter-length
gr()
plt1 = plot(xlabel="\$t/T\$", ylabel="\$u_3\$ at \$x_1=L/4\$")
plot!(tNorm,u₃_quarter, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],u₃_quarter_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt1)
savefig(string(absPath,"/initialDispAndVelBeam_disp.pdf"))
# Velocity at quarter-length
plt2 = plot(xlabel="\$t/T\$", ylabel="\$V_3\$ at \$x_1=L/4\$ [\$1\$/s]")
plot!(tNorm,V₃_quarter, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],V₃_quarter_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt2)
savefig(string(absPath,"/initialDispAndVelBeam_vel.pdf"))
# Acceleration at quarter-length
plt3 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3\$ at \$x_1=L/4\$ [\$1\$/\$s^2\$]")
plot!(tNorm,Vdot₃_quarter, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Vdot₃_quarter_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt3)
savefig(string(absPath,"/initialDispAndVelBeam_acc.pdf"))
# Rotation at root
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\theta/(2\\pi\\delta)\$ at \$x_1=0\$")
plot!(tNorm,θ₂_root/(2*π), c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],θ₂_root_analytic[1:5:end]/(2*π), c=:blue, ms=ms, msw=0, label="Analytical")
display(plt4)
savefig(string(absPath,"/initialDispAndVelBeam_rot.pdf"))
# Angular velocity at mid-length
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_2\$ at \$x_1=L/2\$ [rad/s]")
plot!(tNorm,Ω₂_mid, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ω₂_mid_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt5)
savefig(string(absPath,"/initialDispAndVelBeam_angVel.pdf"))
# Angular acceleration at mid-length
plt6 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_2\$ at \$x_1=L/2\$ [rad/\$s^2\$]")
plot!(tNorm,Ωdot₂_mid, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ωdot₂_mid_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt6)
savefig(string(absPath,"/initialDispAndVelBeam_angAcc.pdf"))

println("Finished initialDispAndVelBeam.jl")