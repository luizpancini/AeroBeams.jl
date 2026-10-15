```@meta
EditURL = "../../../test/examples/dynamicStructural/initialDispAndVelBeam.jl"
```

# Motion of a simply supported beam under initial conditions
This example illustrates how to set up a dynamic analysis with initial conditions. To do that, we take a simply supported beam under a specified combination of sinusoidal initial displacement and velocity.

!!! tip
    The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/dynamicStructural/initialDispAndVelBeam.jl).

### Initial conditions
The setup of our problem begins by specifying the initial conditions, defined by the variable `initialConditions`. These can be defined by the displacement (`u3`), the angle of rotation (`θ2`), or both (in this example, each one with its time rate). The package makes the initial states of the dynamic problem consistent with whichever data are specified, so the three options lead to practically the same solution. The displacement and rotation (and their time rates) correspond to those of a simply supported beam's second bending mode shape. Notice that `θ2` is the angle of rotation over the beam's arclength coordinate (`x1`) and `p2` is the corresponding Wiener-Milenkovic rotation parameter.

````@example initialDispAndVelBeam
using AeroBeams

# Select initial conditions as either :displacement, :rotation or :both
initialConditions = :displacement

# Initial conditions: sinusoidal displacement and velocity of u3 and/or rotation of θ2
δ = 1e-3
σ = 1e-1
u3 = x1 -> δ*sin.(2π*x1/L)
θ2 = x1 -> -δ*2π/L*cos.(2π*x1/L)
p2 = x1 -> 4*tan.(θ2(x1)/4)
udot3 = x1 -> σ*sin.(2π*x1/L)
θdot2 = x1 -> -σ*2π/L*cos.(2π*x1/L)
pdot2 = x1 -> sec.(θ2(x1)/4).^2 .* θdot2(x1)
nothing #hide
````

### Beam
We create a simple isotropic beam for this example. The only relevant properties are the beam's length (`L`), bending stiffness (`EI`) and mass per unit length (`ρA`). The initial conditions are appended to the beam accordingly, with the function [`add_initial_conditions_to_beam!`](@ref add_initial_conditions_to_beam!).

````@example initialDispAndVelBeam
# Beam
L = 1
EI = 1
ρA = 1
∞ = 1e4
nElem = 40
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,EIy=EI)],I=[inertia_matrix(ρA=ρA)])

# Set initial conditions
if initialConditions == :displacement
    add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:u0_of_x1 => x1->[0; 0; u3(x1)], :udot0_of_x1 => x1->[0; 0; udot3(x1)]))
elseif initialConditions == :rotation
    add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:p0_of_x1 => x1->[0; p2(x1); 0], :pdot0_of_x1 => x1->[0; pdot2(x1); 0]))
elseif initialConditions == :both
    add_initial_conditions_to_beam!(beam=beam, conditions=Dict(:u0_of_x1 => x1->[0; 0; u3(x1)], :p0_of_x1 => x1->[0; p2(x1); 0], :udot0_of_x1 => x1->[0; 0; udot3(x1)], :pdot0_of_x1 => x1->[0; pdot2(x1); 0]))
else
    error("initialConditions must be either :displacement, :rotation or :both")
end
nothing #hide
````

### Boundary conditions
The beam is simply-supported: pinned at one end and with a roller at the other.

````@example initialDispAndVelBeam
# BCs
pin = create_BC(name="pin",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p3A=>0))
roller = create_BC(name="roller",beam=beam,node=nElem+1,conditions=Dict(:u3A=>0))
nothing #hide
````

### Model
Our model consists of the beam and the associated boundary conditions. Notice that gravity is neglected.

````@example initialDispAndVelBeam
# Model
initialDispAndVelBeam = create_Model(name="initialDispAndVelBeam",beams=[beam],BCs=[pin,roller])
nothing #hide
````

### Creating the dynamic problem
The time and frequency variables of our problem are defined by the beam's second bending mode frequency, `ω2`, whose analytical value is known. We'll run the simulation for one cycle, and with a time step one hundredth of that.

````@example initialDispAndVelBeam
# Time and frequency variables
ω2 = (2π/L)^2 * sqrt(EI/ρA)
T = 2π/ω2
cycles = 1
tf = cycles * T
Δt = T/100
nothing #hide
````

In initial conditions problems, in general we do not know the initial values of all the states of the formulation (displacements, rotations, internal loads and sectional velocities) and of their time rates. Thus, they have to be determined together with the initial step, which is automatically done by the dynamic problem: the unspecified displacements and rotations are those of the static equilibrium, the unspecified velocities follow the specified ones as in a quasi-static motion, and the accelerations satisfy the equations of motion at the initial time. This procedure does not depend on any tuning parameters. Now we can create the problem and solve it.

````@example initialDispAndVelBeam
# Create and solve the problem
using Suppressor #hide
problem = create_DynamicProblem(model=initialDispAndVelBeam,finalTime=tf,Δt=Δt)
@suppress begin #hide
solve!(problem)
end #hide
nothing #hide
````

### Post-processing
First, we unpack the numerical solution of our variables of interest and compute analytical values.

````@example initialDispAndVelBeam
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
u3_quarter_analytic = δ*cos.(ω2*t)*sin(2π*1/4) + σ/ω2*sin.(ω2*t)*sin(2π*1/4)
V3_quarter_analytic = -δ*ω2*sin.(ω2*t)*sin(2π*1/4) + σ*cos.(ω2*t)*sin(2π*1/4)
Vdot3_quarter_analytic = -δ*ω2^2*cos.(ω2*t)*sin(2π*1/4) + -σ*ω2*sin.(ω2*t)*sin(2π*1/4)
θ2_root_analytic = -δ*2π/L*cos.(ω2*t)*cos(2π*0) + -σ*2π/L/ω2*sin.(ω2*t)*cos(2π*0)
Ω2_mid_analytic = δ*2π/L*ω2*sin.(ω2*t)*cos(2π*1/2) + -σ*2π/L*cos.(ω2*t)*cos(2π*1/2)
Ωdot2_mid_analytic = δ*2π/L*ω2^2*cos.(ω2*t)*cos(2π*1/2) + σ*2π/L*ω2*sin.(ω2*t)*cos(2π*1/2)
nothing #hide
````

Now let's plot the displacement, rotation, velocities and accelerations at the appropriate points of the beam. The correspondance with the analytical reference is very good.

````@example initialDispAndVelBeam
gr()
ENV["GKSwstype"] = "100" #hide

# Displacement at quarter-length
plt1 = plot(xlabel="\$t/T\$", ylabel="\$u_3\$ at \$x_1=L/4\$ [m]")
plot!(tNorm,u3_quarter, c=:black, lw=2, label="AeroBeams")
scatter!(tNorm[1:5:end],u3_quarter_analytic[1:5:end], c=:blue, ms=5, msw=0, label="Analytical")
savefig("initialDispAndVelBeam_disp.svg") #hide

# Velocity at quarter-length
plt2 = plot(xlabel="\$t/T\$", ylabel="\$V_3\$ at \$x_1=L/4\$ [m/s]")
plot!(tNorm,V3_quarter, c=:black, lw=2, label="AeroBeams")
scatter!(tNorm[1:5:end],V3_quarter_analytic[1:5:end], c=:blue, ms=5, msw=0, label="Analytical")
savefig("initialDispAndVelBeam_vel.svg") #hide

# Acceleration at quarter-length
plt3 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3\$ at \$x_1=L/4\$ [m/\$s^2\$]")
plot!(tNorm,Vdot3_quarter, c=:black, lw=2, label="AeroBeams")
scatter!(tNorm[1:5:end],Vdot3_quarter_analytic[1:5:end], c=:blue, ms=5, msw=0, label="Analytical")
savefig("initialDispAndVelBeam_acc.svg") #hide

# Rotation at root
plt4 = plot(xlabel="\$t/T\$", ylabel="\$\\theta/(2\\pi)\$ at \$x_1=0\$")
plot!(tNorm,θ2_root/(2π), c=:black, lw=2, label="AeroBeams")
scatter!(tNorm[1:5:end],θ2_root_analytic[1:5:end]/(2π), c=:blue, ms=5, msw=0, label="Analytical")
savefig("initialDispAndVelBeam_rot.svg") #hide

# Angular velocity at mid-length
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_2\$ at \$x_1=L/2\$ [rad/s]")
plot!(tNorm,Ω2_mid, c=:black, lw=2, label="AeroBeams")
scatter!(tNorm[1:5:end],Ω2_mid_analytic[1:5:end], c=:blue, ms=5, msw=0, label="Analytical")
savefig("initialDispAndVelBeam_angVel.svg") #hide

# Angular acceleration at mid-length
plt6 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_2\$ at \$x_1=L/2\$ [rad/\$s^2\$]")
plot!(tNorm,Ωdot2_mid, c=:black, lw=2, label="AeroBeams")
scatter!(tNorm[1:5:end],Ωdot2_mid_analytic[1:5:end], c=:blue, ms=5, msw=0, label="Analytical")
savefig("initialDispAndVelBeam_angAcc.svg") #hide
nothing #hide
````

![](initialDispAndVelBeam_disp.svg)
![](initialDispAndVelBeam_vel.svg)
![](initialDispAndVelBeam_acc.svg)
![](initialDispAndVelBeam_rot.svg)
![](initialDispAndVelBeam_angVel.svg)
![](initialDispAndVelBeam_angAcc.svg)

Finally, we can visualize the motion of the beam using the function [`plot_dynamic_deformation`](@ref plot_dynamic_deformation) with the appropriate inputs.

````@example initialDispAndVelBeam
# Animation
@suppress begin #hide
plot_dynamic_deformation(problem,scale=1/δ/10,scalePos=[0.15;-0.05;0],timeStampPos=[0.5;-0.05;0],plotFrequency=1,plotLimits=([0,L],[-L/2,L/2],[-L/3,L/3]),save=true,savePath="/docs/build/literate/initialDispAndVelBeam_motion.gif")
nothing #hide
end #hide
````

![](initialDispAndVelBeam_motion.gif)

````@example initialDispAndVelBeam
nothing #hide
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

