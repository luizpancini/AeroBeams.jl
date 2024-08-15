using AeroBeams, LinearAlgebra, Plots
 
# Select initial conditions as either "displacement", "rotation" or "both"
initialConditions = "displacement"

# Initial conditions: displacement of u₃ and/or rotation of θ₂
δ = 1e-3
u₃ = x1 -> δ * (sin.(π*x1/L) - π*x1/L.*(1.0 .- x1/L))
θ₂ = x1 -> -δ * (π/L*cos.(π*x1/L) - π/L*(1.0 .- 2*x1/L))
p₂ = x1 -> 4*tan.(θ₂(x1)/4)

# Beam
L = 1.0
EA,GA,GJ,EIy,EIz = 1e6,1e4,1e9,1.0,1e9
ρA,ρI = 1.0,0
nElem = 40
stiffnessMatrix = diagm([EA,GA,GA,GJ,EIy,EIz])
inertiaMatrix = diagm([ρA,ρA,ρA,ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])
if initialConditions == "displacement"
    beam.u0_of_x1=x1->[0; 0; u₃(x1)]
elseif initialConditions == "rotation"
    beam.p0_of_x1=x1->[0; p₂(x1); 0]
elseif initialConditions == "both"
    beam.u0_of_x1=x1->[0; 0; u₃(x1)]
    beam.p0_of_x1=x1->[0; p₂(x1); 0]
else
    error("Wrong initialConditions")
end
update_beam!(beam)

# BCs
clamp1 = create_BC(name="clamp1",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=beam,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
biclampedBeam = create_Model(name="biclampedBeam",beams=[beam],BCs=[clamp1,clamp2])

# Time and frequency variables
ω = 2*π*3.5653
T = 2*π/ω
cycles = 2
tf = cycles*T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=biclampedBeam,finalTime=tf,Δt=Δt)
solve!(problem)
# @time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tNorm = t/T
midElem = div(nElem,2)
quarterElem = div(nElem,4)
u₃_mid = [problem.nodalStatesOverTime[i][midElem].u_n2[3] for i in 1:length(tNorm)]
V₃_mid = [(problem.elementalStatesOverTime[i][midElem].V[3]+problem.elementalStatesOverTime[i][midElem+1].V[3])/2.0 for i in 1:length(tNorm)]
Vdot₃_mid = [(problem.elementalStatesRatesOverTime[i][midElem].Vdot[3]+problem.elementalStatesRatesOverTime[i][midElem+1].Vdot[3])/2.0 for i in 1:length(tNorm)]
θ₂_quarter = [problem.nodalStatesOverTime[i][quarterElem].θ_n2 for i in 1:length(tNorm)]
Ω₂_quarter = [(problem.elementalStatesOverTime[i][quarterElem].Ω[2]+problem.elementalStatesOverTime[i][quarterElem+1].Ω[2])/2.0 for i in 1:length(tNorm)]
Ωdot₂_quarter = [(problem.elementalStatesRatesOverTime[i][quarterElem].Ωdot[2]+problem.elementalStatesRatesOverTime[i][quarterElem+1].Ωdot[2])/2.0 for i in 1:length(tNorm)]

# Compute analytical values
u₃_mid_analytic = δ*cos.(ω*t)*(sin(π/2)-π/4)
V₃_mid_analytic = -δ*ω*sin.(ω*t)*(sin(π/2)-π/4)
Vdot₃_mid_analytic = -δ*ω^2*cos.(ω*t)*(sin(π/2)-π/4)
θ₂_quarter_analytic = -δ*π/L*cos.(ω*t)*(cos(π/4)-1/2)
Ω₂_quarter_analytic = δ*ω*π/L*sin.(ω*t)*(cos(π/4)-1/2)
Ωdot₂_quarter_analytic = δ*ω^2*π/L*cos.(ω*t)*(cos(π/4)-1/2)

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 6
relPath = "/test/outputs/figures/biclampedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,scale=1/δ,plotLimits=[(0,L),(-0.25,0.25),(0,1)],save=true,savePath=string(relPath,"/biclampedBeam_deformation.gif"))
# Normalized displacement at mid-length
gr()
plt1 = plot(xlabel="\$t/T\$", ylabel="\$u_3/\\delta\$ at \$x_1=L/2\$")
plot!(tNorm,u₃_mid/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],u₃_mid_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt1)
savefig(string(absPath,"/biclampedBeam_midDisp.pdf"))
# Normalized velocity at mid-length 
plt2 = plot(xlabel="\$t/T\$", ylabel="\$V_3/\\delta\$ at \$x_1=L/2\$ [\$1\$/s]")
plot!(tNorm,V₃_mid/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],V₃_mid_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt2)
savefig(string(absPath,"/biclampedBeam_midVel.pdf"))
# Normalized acceleration at mid-length
plt3 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3/\\delta\$ at \$x_1=L/2\$ [\$1\$/\$s^2\$]")
plot!(tNorm,Vdot₃_mid/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Vdot₃_mid_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt3)
savefig(string(absPath,"/biclampedBeam_midAcc.pdf"))
# Normalized rotation at quarter-length
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\theta/(2\\pi\\delta)\$ at \$x_1=L/4\$")
plot!(tNorm,θ₂_quarter/(2*π)/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],θ₂_quarter_analytic[1:5:end]/(2*π)/δ, c=:blue, ms=ms, msw=0, label="Analytical")
display(plt4)
savefig(string(absPath,"/biclampedBeam_quarterRot.pdf"))
# Angular velocity at quarter-length
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_2\$ at \$x_1=L/4\$ [rad/s]")
plot!(tNorm,Ω₂_quarter, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ω₂_quarter_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt5)
savefig(string(absPath,"/biclampedBeam_quarterAngVel.pdf"))
# Angular acceleration at quarter-length
plt6 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_2\$ at \$x_1=L/4\$ [rad/\$s^2\$]")
plot!(tNorm,Ωdot₂_quarter, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ωdot₂_quarter_analytic[1:5:end], c=:blue, ms=ms, msw=0, label="Analytical")
display(plt6)
savefig(string(absPath,"/biclampedBeam_quarterAngAcc.pdf"))

println("Finished biclampedBeam.jl")