using AeroBeams

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Airspeed
Ma = 0.5
U = Ma*atmosphere.a

# Wing surface data
aeroSolver = Indicial()
flapSolver = ThinAirfoilTheory()
airfoil = deepcopy(FlatPlate)
chord = 0.1
normSparPos = 0.25
normFlapPos = 0.75
normFlapSpan = [0; 1]

# Flap deflection profile
k = 0.098
A = deg2rad(2.5)
ω = k*U/(chord/2)
δ = t -> A*sin(ω*t)

# Update airfoil parameters
update_all_airfoil_parameters!(airfoil,Ma)

# Create wing surfaces
surf1 = create_AeroSurface(ppSolver=aeroSolver,flapSolver=flapSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,normFlapPos=normFlapPos,normFlapSpan=normFlapSpan,δ=δ)
surf2 = create_AeroSurface(ppSolver=aeroSolver,flapSolver=flapSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,normFlapPos=normFlapPos,normFlapSpan=normFlapSpan)

# Wing beams
L = 1
nElem = 1
wing1 = create_Beam(name="wing1",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e16)],I=[inertia_matrix(ρA=1)],surface=surf1)
wing2 = create_Beam(name="wing2",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e16)],I=[inertia_matrix(ρA=1)],surface=surf2)

# Flap links
m = -1
flapLink = create_FlapLink(masterBeam=wing1,slaveBeams=[wing2],δMultipliers=[m])

# BCs
clamp1 = create_BC(name="clamp1",beam=wing1,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
clamp2 = create_BC(name="clamp2",beam=wing2,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
flapOscillationLinked = create_Model(name="flapOscillationLinked",beams=[wing1,wing2],BCs=[clamp1,clamp2],flapLinks=[flapLink],atmosphere=atmosphere,v_A=[0;U;0])

# Time variables
T = 2π/ω
cycles = 10
tf = cycles*T
Δt = T/100

# Create and solve problem
problem = create_DynamicProblem(model=flapOscillationLinked,finalTime=tf,Δt=Δt)
solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
aeroVars = sol.aeroVariables
tNorm = t/T
cnMaster = [aeroVars[i][1].aeroCoefficients.cn for i in eachindex(t)]
cmMaster = [aeroVars[i][1].aeroCoefficients.cm for i in eachindex(t)]
cnSlave = [aeroVars[i][2].aeroCoefficients.cn for i in eachindex(t)]
cmSlave = [aeroVars[i][2].aeroCoefficients.cm for i in eachindex(t)]

## Plots
# Set paths
relPath = "/test/outputs/figures/flapOscillationLinked"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
anim = plot_dynamic_deformation(problem,follow=:basisA,plotFrequency=10,plotLimits=([0,2*L],[-L,L],[-L,L]),save=true,savePath=string(relPath,"/flapOscillationLinked_deformation.gif"),displayProgress=true,plotAeroSurf=false)
display(anim)

# Plot configurations
lw = 2
gr()

# cn and cm vs time
plt11 = plot(ylabel="\$c_n\$", xlims=[0,cycles])
plot!(tNorm, [cnMaster, cnSlave], lw=lw, label=["Master" "Slave"])
plt12 = plot(xlabel="\$t/T\$", ylabel="\$c_m\$", xlims=[0,cycles])
plot!(tNorm, [cmMaster, cmSlave], lw=lw, label=false)
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/flapOscillationLinked.pdf"))

println("Finished flapOscillationLinked.jl")