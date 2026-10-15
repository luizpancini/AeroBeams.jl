using AeroBeams

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Airspeed
Ma = 0.5
U = Ma*atmosphere.a

# Wing surface data
aeroSolver = BLi()
flapSolver = ThinAirfoilTheory()
airfoil = deepcopy(FlatPlate)
chord = 0.18
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

# Create wing surface
surf = create_AeroSurface(ppSolver=aeroSolver,flapSolver=flapSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,normFlapPos=normFlapPos,normFlapSpan=normFlapSpan,δ=δ)

# Wing beam
L = 1
nElem = 1
wing = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e16)],I=[inertia_matrix(ρA=1)],surface=surf)

# BCs
clamp1 = create_BC(name="clamp1",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
flapOscillation = create_Model(name="flapOscillation",beams=[wing],BCs=[clamp1,clamp2],v_A=[0;U;0],atmosphere=atmosphere)

# Time variables
T = 2π/ω
cycles = 10
tf = cycles*T
Δt = T/100

# Create and solve problem
problem = create_DynamicProblem(model=flapOscillation,finalTime=tf,Δt=Δt)
@time solve!(problem)
# @profview_allocs solve!(problem) sample_rate=1e-1
# using Profile
# Profile.clear()
# Profile.init(n=10^8,delay=1e-6)
# @profview solve!(problem)

# Unpack numerical solution
sol = problem.solutionHistory
t = sol.t
aeroVars = sol.aeroVariables
cn = [aeroVars[i][1].aeroCoefficients.cn for i in eachindex(t)]
cm = [aeroVars[i][1].aeroCoefficients.cm for i in eachindex(t)]

# Load reference data by TIJDEMAN & SCHIPPERS (1973) and LEISHMAN (2006)
cnExp = readdlm(pkgdir(AeroBeams,"test","referenceData","flapOscillation", "cnVsDeltaExp.txt"))
cmExp = readdlm(pkgdir(AeroBeams,"test","referenceData","flapOscillation", "cmVsDeltaExp.txt"))
cnRefMod = readdlm(pkgdir(AeroBeams,"test","referenceData","flapOscillation", "cnVsDeltaRefMod.txt"))
cmRefMod = readdlm(pkgdir(AeroBeams,"test","referenceData","flapOscillation", "cmVsDeltaRefMod.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/flapOscillation"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
rangeLastCycle = findfirst(x->x>=tf-T,t):length(t)
ts = 10
fs = 16
lw = 2
ms = 3
gr()

# cn vs δ
plt1 = plot(xlabel="\$\\delta\$ [deg]", ylabel="\$c_n/\\pi\$", xlims=[-3,3], ylims=[-0.075,0.075], tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
plot!(δ.(t[rangeLastCycle])*180/π, cn[rangeLastCycle]/π, c=:black, lw=lw, label="AeroBeams")
plot!(cnRefMod[1,:], cnRefMod[2,:]/π, c=:black, ls=:dashdot, lw=lw, label="Model - Leishman (2006)")
scatter!(cnExp[1,:], cnExp[2,:]/π, c=:black, ms=ms, msw=0, label="Exp. - Tijdeman & Schippers (1973)")
display(plt1)
savefig(string(absPath,"/flapOscillation_cn.pdf"))

# cm vs δ
plt2 = plot(xlabel="\$\\delta\$ [deg]", ylabel="\$-2c_m/\\pi\$", xlims=[-3,3], ylims=[-0.03,0.03], tickfont=font(ts), guidefont=font(fs))
plot!(δ.(t[rangeLastCycle])*180/π, -2*cm[rangeLastCycle]/π, c=:black, lw=lw, label=false)
plot!(cmRefMod[1,:], -2*cmRefMod[2,:]/π, c=:black, ls=:dashdot, lw=lw, label=false)
scatter!(cmExp[1,:], -2*cmExp[2,:]/π, c=:black, ms=ms, msw=0, label=false)
display(plt2)
savefig(string(absPath,"/flapOscillation_cm.pdf"))

println("Finished flapOscillation.jl")