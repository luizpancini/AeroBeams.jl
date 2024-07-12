using AeroBeams, LinearAlgebra, Plots, DelimitedFiles

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Airspeed
Ma = 0.5
U = Ma*atmosphere.a

# Wing surface data
airfoil = flatPlate
chord = 0.18
normSparPos = 0.25
normFlapPos = 0.75
normFlapSpan = [0; 1]

# Flap deflection profile
k = 0.098
A = 2.5*π/180
ω = k*U/(chord/2)
δ = t -> A*sin(ω*t)

# Create wing surface
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos,normFlapPos=normFlapPos,normFlapSpan=normFlapSpan,δ=δ)

# Wing beam
L = 1
nElem = 1
∞ = 1e10
wing = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],aeroSurface=surf)

# BCs
clamp1 = create_BC(name="clamp1",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
flapOscillation = create_Model(name="flapOscillation",beams=[wing],BCs=[clamp1,clamp2],atmosphere=atmosphere,v_A=[0;U;0])

# Time variables
T = 2π/ω
cycles = 10
tf = cycles*T
Δt = T/100

# Create and solve problem
problem = create_DynamicProblem(model=flapOscillation,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
t = problem.timeVector
cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:length(t)]
cm = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cm for i in 1:length(t)]

# Load reference data by TIJDEMAN & SCHIPPERS (1973) and LEISHMAN (2006)
cnExp = readdlm(string(pwd(),"/test/referenceData/flapOscillation/cnVsDeltaExp.txt"))
cmExp = readdlm(string(pwd(),"/test/referenceData/flapOscillation/cmVsDeltaExp.txt"))
cnRefMod = readdlm(string(pwd(),"/test/referenceData/flapOscillation/cnVsDeltaRefMod.txt"))
cmRefMod = readdlm(string(pwd(),"/test/referenceData/flapOscillation/cmVsDeltaRefMod.txt"))

# Plots
# ------------------------------------------------------------------------------
rangeLastCycle = findfirst(x->x==tf-T,t):length(t)
lw = 2
ms = 3
# cn vs δ
plt1 = plot(xlabel="\$\\delta\$ [deg]", ylabel="\$c_n/\\pi\$", xlims=[-3,3], ylims=[-0.075,0.075])
plot!(δ.(t[rangeLastCycle])*180/π, cn[rangeLastCycle]/π, c=:black, lw=lw, label="AeroBeams")
plot!(cnRefMod[1,:], cnRefMod[2,:]/π, c=:black, ls=:dash, lw=lw, label="Incompressible thoery by Leishman (2006)")
scatter!(cnExp[1,:], cnExp[2,:]/π, c=:black, ms=ms, msw=0, label="Experiment by Tijdeman & Schippers (1973)")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/flapOscillation_cn.pdf"))
# cm vs δ
plt2 = plot(xlabel="\$\\delta\$ [deg]", ylabel="\$-2c_m/\\pi\$", xlims=[-3,3], ylims=[-0.03,0.03])
plot!(δ.(t[rangeLastCycle])*180/π, -2*cm[rangeLastCycle]/π, c=:black, lw=lw, label="AeroBeams")
plot!(cmRefMod[1,:], -2*cmRefMod[2,:]/π, c=:black, ls=:dash, lw=lw, label="Incompressible thoery by Leishman (2006)")
scatter!(cmExp[1,:], -2*cmExp[2,:]/π, c=:black, ms=ms, msw=0, label="Experiment by Tijdeman & Schippers (1973)")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/flapOscillation_cm.pdf"))

println("Finished flapOscillation.jl")