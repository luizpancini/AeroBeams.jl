using AeroBeams, LinearAlgebra
 
# Select initial conditions as either "displacement", "rotation" or "both"
initialConditions = "displacement"

# Initial conditions: displacement of u3 and/or rotation of θ2
δ = 1e-3
u3 = x1 -> δ * (sin.(π*x1/L) - π*x1/L.*(1.0 .- x1/L))
θ2 = x1 -> -δ * (π/L*cos.(π*x1/L) - π/L*(1.0 .- 2*x1/L))
p2 = x1 -> 4*tan.(θ2(x1)/4)

# Beam
L = 1.0
EA,GA,GJ,EIy,EIz = 1e6,1e4,1e9,1.0,1e9
ρA,ρI = 1.0,0
nElem = 40
stiffnessMatrix = diagm([EA,GA,GA,GJ,EIy,EIz])
inertiaMatrix = diagm([ρA,ρA,ρA,ρI,ρI,ρI])
beam = create_Beam(name="beam",length=L,nElements=nElem,C=[stiffnessMatrix],I=[inertiaMatrix])
if initialConditions == "displacement"
    beam.u0_of_x1=x1->[0; 0; u3(x1)]
elseif initialConditions == "rotation"
    beam.p0_of_x1=x1->[0; p2(x1); 0]
elseif initialConditions == "both"
    beam.u0_of_x1=x1->[0; 0; u3(x1)]
    beam.p0_of_x1=x1->[0; p2(x1); 0]
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
ω = 2π*3.5653
T = 2π/ω
cycles = 2
tf = cycles*T
Δt = T/100

# Create and solve the problem
problem = create_DynamicProblem(model=biclampedBeam,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tNorm = t/T
midElem = div(nElem,2)
quarterElem = div(nElem,4)
u3_mid = [problem.nodalStatesOverTime[i][midElem].u_n2[3] for i in 1:length(tNorm)]
V3_mid = [(problem.elementalStatesOverTime[i][midElem].V[3]+problem.elementalStatesOverTime[i][midElem+1].V[3])/2.0 for i in 1:length(tNorm)]
Vdot3_mid = [(problem.elementalStatesRatesOverTime[i][midElem].Vdot[3]+problem.elementalStatesRatesOverTime[i][midElem+1].Vdot[3])/2.0 for i in 1:length(tNorm)]
θ2_quarter = [problem.nodalStatesOverTime[i][quarterElem].θ_n2 for i in 1:length(tNorm)]
Ω2_quarter = [(problem.elementalStatesOverTime[i][quarterElem].Ω[2]+problem.elementalStatesOverTime[i][quarterElem+1].Ω[2])/2.0 for i in 1:length(tNorm)]
Ωdot2_quarter = [(problem.elementalStatesRatesOverTime[i][quarterElem].Ωdot[2]+problem.elementalStatesRatesOverTime[i][quarterElem+1].Ωdot[2])/2.0 for i in 1:length(tNorm)]

# Compute analytical values
u3_mid_analytic = δ*cos.(ω*t)*(sin(π/2)-π/4)
V3_mid_analytic = -δ*ω*sin.(ω*t)*(sin(π/2)-π/4)
Vdot3_mid_analytic = -δ*ω^2*cos.(ω*t)*(sin(π/2)-π/4)
θ2_quarter_analytic = -δ*π/L*cos.(ω*t)*(cos(π/4)-1/2)
Ω2_quarter_analytic = δ*ω*π/L*sin.(ω*t)*(cos(π/4)-1/2)
Ωdot2_quarter_analytic = δ*ω^2*π/L*cos.(ω*t)*(cos(π/4)-1/2)

println("Finished biclampedBeam.jl")