using AeroBeams
 
# Select initial conditions input as either :displacement, :rotation or :both
initialConditions = :both

# Initial conditions: displacement of u3 and/or rotation of θ2
δ = 1e-3
u3 = x1 -> δ * (sin.(π*x1/L) - π*x1/L.*(1.0 .- x1/L))
θ2 = x1 -> -δ * (π/L*cos.(π*x1/L) - π/L*(1.0 .- 2*x1/L))
p2 = x1 -> 4*tan.(θ2(x1)/4)

# Beam
L = 1
EA = 1e6
GA = 1e4
EIy = 1
∞ = 1e9
ρA = 1
nElem = 40
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=EA,GAy=GA,GAz=GA,EIy=EIy,∞=∞)],I=[inertia_matrix(ρA=ρA)])
if initialConditions == :displacement
    add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:u0_of_x1 => x1->[0; 0; u3(x1)]))
elseif initialConditions == :rotation
    add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:p0_of_x1 => x1->[0; p2(x1); 0]))
elseif initialConditions == :both
    add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:u0_of_x1 => x1->[0; 0; u3(x1)], :p0_of_x1 => x1->[0; p2(x1); 0]))
else
    error("initialConditions must be either :displacement, :rotation or :both")
end

# BCs
clamp1 = create_BC(name="clamp1",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
clamp2 = create_BC(name="clamp2",beam=beam,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
biclampedBeam = create_Model(name="biclampedBeam",beams=[beam],BCs=[clamp1,clamp2])

# Time and frequency variables
ω = Hz2rads(3.5653)
T = 2π/ω
cycles = 2
tf = cycles * T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=biclampedBeam,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
states = sol.states
Nt = length(t)
tNorm = t/T
midElem = div(nElem,2)
quarterElem = div(nElem,4)
u3_mid = [states[i][midElem].nodal.u2[3] for i in 1:Nt]
V3_mid = [(states[i][midElem].midpoint.V[3]+states[i][midElem+1].midpoint.V[3])/2 for i in 1:Nt]
Vdot3_mid = [(states[i][midElem].rates.Vdot[3]+states[i][midElem+1].rates.Vdot[3])/2 for i in 1:Nt]
θ2_quarter = [states[i][quarterElem].nodal.θ[2] for i in 1:Nt]
Ω2_quarter = [(states[i][quarterElem].midpoint.Ω[2]+states[i][quarterElem+1].midpoint.Ω[2])/2 for i in 1:Nt]
Ωdot2_quarter = [(states[i][quarterElem].rates.Ωdot[2]+states[i][quarterElem+1].rates.Ωdot[2])/2 for i in 1:Nt]

# Compute analytical values
u3_mid_analytic = δ*cos.(ω*t)*(sin(π/2)-π/4)
V3_mid_analytic = -δ*ω*sin.(ω*t)*(sin(π/2)-π/4)
Vdot3_mid_analytic = -δ*ω^2*cos.(ω*t)*(sin(π/2)-π/4)
θ2_quarter_analytic = -δ*π/L*cos.(ω*t)*(cos(π/4)-1/2)
Ω2_quarter_analytic = δ*ω*π/L*sin.(ω*t)*(cos(π/4)-1/2)
Ωdot2_quarter_analytic = δ*ω^2*π/L*cos.(ω*t)*(cos(π/4)-1/2)

## Plots
# Set paths
relPath = "/test/outputs/figures/biclampedBeam"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,scale=1/δ,plotLimits=([0,L],[-0.25,0.25],[-0.25,0.25]),save=true,savePath=string(relPath,"/biclampedBeam_deformation.gif"))
display(anim)

# Plot configurations
lw = 2
ms = 6
msw = 0
gr()

# Normalized displacement at mid-length
plt1 = plot(xlabel="\$t/T\$", ylabel="\$u_3/\\delta\$ at \$x_1=L/2\$")
plot!(tNorm,u3_mid/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],u3_mid_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt1)
savefig(string(absPath,"/biclampedBeam_midDisp.pdf"))

# Normalized velocity at mid-length 
plt2 = plot(xlabel="\$t/T\$", ylabel="\$V_3/\\delta\$ at \$x_1=L/2\$ [\$1\$/s]")
plot!(tNorm,V3_mid/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],V3_mid_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt2)
savefig(string(absPath,"/biclampedBeam_midVel.pdf"))

# Normalized acceleration at mid-length
plt3 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3/\\delta\$ at \$x_1=L/2\$ [\$1\$/\$s^2\$]")
plot!(tNorm,Vdot3_mid/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Vdot3_mid_analytic[1:5:end]/δ, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt3)
savefig(string(absPath,"/biclampedBeam_midAcc.pdf"))

# Normalized rotation at quarter-length
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\theta/(2\\pi\\delta)\$ at \$x_1=L/4\$")
plot!(tNorm,θ2_quarter/(2*π)/δ, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],θ2_quarter_analytic[1:5:end]/(2*π)/δ, c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt4)
savefig(string(absPath,"/biclampedBeam_quarterRot.pdf"))

# Angular velocity at quarter-length
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_2\$ at \$x_1=L/4\$ [rad/s]")
plot!(tNorm,Ω2_quarter, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ω2_quarter_analytic[1:5:end], c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt5)
savefig(string(absPath,"/biclampedBeam_quarterAngVel.pdf"))

# Angular acceleration at quarter-length
plt6 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_2\$ at \$x_1=L/4\$ [rad/\$s^2\$]")
plot!(tNorm,Ωdot2_quarter, c=:black, lw=lw, label="Numerical")
scatter!(tNorm[1:5:end],Ωdot2_quarter_analytic[1:5:end], c=:blue, ms=ms, msw=msw, label="Analytical")
display(plt6)
savefig(string(absPath,"/biclampedBeam_quarterAngAcc.pdf"))

println("Finished biclampedBeam.jl")