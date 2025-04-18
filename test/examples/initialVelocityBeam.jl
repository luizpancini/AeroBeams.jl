using AeroBeams, LinearAlgebra

# Select initial conditions as either "displacement", "rotation" or "both"
initialConditions = "displacement"

# Select boundary conditions as either "ss-ss" (simple support on both sides) or "ss-roller" (simple support on the left and roller on the right)
BCType = "ss-roller"

# Initial conditions: sinusoidal velocity V3 and/or angular velocity Ω2
σ = 1e-1
udot3 = x1 -> σ*sin.(2*π*x1/L)
θdot2 = x1 -> -σ*2*π/L*cos.(2*π*x1/L)
pdot2 = x1 -> sec.(0/4).^2 .* θdot2(x1)

# Beam
L = 1
EIy = 1
ρA = 1
ρI = 0
∞ = 1e4
nElem = 48
stiffnessMatrix = diagm([∞,∞,∞,∞,EIy,∞])
inertiaMatrix = diagm([ρA,ρA,ρA,ρI,ρI,0])
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[stiffnessMatrix],I=[inertiaMatrix])
if initialConditions == "displacement"
    beam.udot0_of_x1=x1->[0; 0; udot3(x1)]
elseif initialConditions == "rotation"
    beam.pdot0_of_x1=x1->[0; pdot2(x1); 0]
elseif initialConditions == "both"
    beam.udot0_of_x1=x1->[0; 0; udot3(x1)]
    beam.pdot0_of_x1=x1->[0; pdot2(x1); 0]
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
initialVelocityBeam = create_Model(name="initialVelocityBeam",beams=[beam],BCs=bcs)

# Time and frequency variables
ω2 = (2*π/L)^2*sqrt(EIy/ρA)
T = 2*π/ω2
cycles = 1
tf = cycles*T
Δt = T/100

# Initial velocities update options
if initialConditions == "displacement"
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=true, relaxFactor = 0.5, Δt=Δt/1e2)
elseif initialConditions == "rotation"
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=100,tol=1e-8, displayProgress=true)
elseif initialConditions == "both"
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-10, displayProgress=true, relaxFactor = 0.5, Δt=Δt/1e2)
end

# Create and solve the problem
problem = create_DynamicProblem(model=initialVelocityBeam,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tNorm = t/T
u3_quarter = [problem.nodalStatesOverTime[i][div(nElem,4)].u_n2[3] for i in 1:length(tNorm)]
V3_quarter = [(problem.elementalStatesOverTime[i][div(nElem,4)].V[3]+problem.elementalStatesOverTime[i][div(nElem,4)+1].V[3])/2.0 for i in 1:length(tNorm)]
Vdot3_quarter = [(problem.elementalStatesRatesOverTime[i][div(nElem,4)].Vdot[3]+problem.elementalStatesRatesOverTime[i][div(nElem,4)+1].Vdot[3])/2.0 for i in 1:length(tNorm)]
θ2_root = [problem.nodalStatesOverTime[i][1].θ_n1 for i in 1:length(tNorm)]
Ω2_mid = [(problem.elementalStatesOverTime[i][div(nElem,2)].Ω[2]+problem.elementalStatesOverTime[i][div(nElem,2)+1].Ω[2])/2.0 for i in 1:length(tNorm)]
Ωdot2_mid = [(problem.elementalStatesRatesOverTime[i][div(nElem,2)].Ωdot[2]+problem.elementalStatesRatesOverTime[i][div(nElem,2)+1].Ωdot[2])/2.0 for i in 1:length(tNorm)]

# Compute analytical values
u3_quarter_analytic = σ/ω2*sin.(ω2*t)*sin(2*π*1/4)
V3_quarter_analytic = σ*cos.(ω2*t)*sin(2*π*1/4)
Vdot3_quarter_analytic = -σ*ω2*sin.(ω2*t)*sin(2*π*1/4)
θ2_root_analytic = -σ*2*π/L/ω2*sin.(ω2*t)*cos(2*π*0)
Ω2_mid_analytic = -σ*2*π/L*cos.(ω2*t)*cos(2*π*1/2)
Ωdot2_mid_analytic = σ*2*π/L*ω2*sin.(ω2*t)*cos(2*π*1/2)

println("Finished initialVelocityBeam.jl")