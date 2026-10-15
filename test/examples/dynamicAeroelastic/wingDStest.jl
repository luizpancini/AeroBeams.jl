using AeroBeams, ForwardDiff

# Flag to apply a transition function to airspeed and pitch profiles (suppresses noise from accelerations)
applyTransitionFunction = true

# Aerodynamic solvers
aeroSolvers = [Indicial(), BLo(), BLi()]

# Frame data
frame = 10_022
airfoil,a₀,a₁,b,k,Ma,U = NASA_frames_loader(frame)

# Pitch profile
ω = k*U/b
τ = 2π/ω
t₀ = -τ/4
θ = applyTransitionFunction ? 
    t -> a₁ * (1+sin(ω*(t+t₀))) * (t/τ)^10 / (1 + (t/τ)^10) : 
    t -> a₁ * (1+sin(ω*(t+t₀)))
p = t -> 4 * tan(θ(t)/4)
pdot = t -> ForwardDiff.derivative(p,t)

# Wing beam properties
L = 10*2b
EIy,GJ = 1e6,1e8
ρA,ρIy,ρIz = 10,1e-2,1e-2
nElem = 10
∞ = 1e12

# Time variables
nCycles = 4
Δt = τ/1000
tf = nCycles*τ

# Elements where output data is to be gathered (root and tip)
elemRangePlot = [1,nElem]

# System solver
systemSolver = create_NewtonRaphson(allowAdvanceThroughUnconvergedAeroStates=true,displayStatus=false)

# Initialize outputs
dynProblem = Array{DynamicProblem}(undef,length(aeroSolvers))
t = Array{Vector{Float64}}(undef,length(aeroSolvers))
α = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
αₑ = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
cn = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
cm = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
ct = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
cl = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
cdrag = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))

# Loop aero solvers
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Display progress
    println("Solving for $(aeroSolvers[i].name)")
    # Aerodynamic surface
    surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=2*b,normSparPos=0.25,updateAirfoilParameters=true)
    # Wing
    wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=10*EIy)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz)],p0=transform_to_WM([0,0,a₀-a₁],from=:E321),surface=surf)
    # Set initial conditions
    add_initial_conditions_to_beam!(beam=wing, conditions=Dict(:pdot0_of_x1 => x1->[pdot(0),0,0]))
    # BCs
    driver = create_BC(name="driver",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>t->p(t),:p2A=>0,:p3A=>0))
    journal = create_BC(name="journal",beam=wing,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p2A=>0,:p3A=>0))
    # Model
    wingDStest = create_Model(name="wingDStest",beams=[wing],BCs=[driver],v_A=[0,U,0])
    # Create and solve steady problem for initial solution
    steadyProblem = create_SteadyProblem(model=wingDStest)
    solve!(steadyProblem)
    # Create and solve dynamic problem
    dynProblem[i] = create_DynamicProblem(model=wingDStest,systemSolver=systemSolver,finalTime=tf,Δt=Δt,x0=steadyProblem.systemArrays.x)
    solve!(dynProblem[i])
    # Unpack numerical solution
    sol = dynProblem[i].solutionHistory
    t[i] = sol.t
    aeroVars = sol.aeroVariables
    α[i] = [[rad2deg.(aeroVars[j][elem].steadyKin.α) for j in eachindex(t[i])] for elem in elemRangePlot]
    αₑ[i] = [[rad2deg.(aeroVars[j][elem].unsteadyKin.αₑ) for j in eachindex(t[i])] for elem in elemRangePlot]
    cn[i] = [[aeroVars[j][elem].aeroCoefficients.cn for j in eachindex(t[i])] for elem in elemRangePlot]
    cm[i] = [[aeroVars[j][elem].aeroCoefficients.cm for j in eachindex(t[i])] for elem in elemRangePlot]
    ct[i] = [[aeroVars[j][elem].aeroCoefficients.ct for j in eachindex(t[i])] for elem in elemRangePlot]
    cl[i] = [cn[i][e].*cosd.(α[i][e]) .+ ct[i][e].*sind.(α[i][e]) for e in eachindex(elemRangePlot)]
    cdrag[i] = [cn[i][e].*sind.(α[i][e]) .- ct[i][e].*cosd.(α[i][e]) for e in eachindex(elemRangePlot)]
end

# Load reference data from McAlister et al.
clRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cl.txt")
cmRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cm.txt")
cdRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cd.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/wingDStest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ind2plot = 1
lw = 2
gr()

# Animation
anim = plot_dynamic_deformation(dynProblem[end],follow=:basisA,view=(30,15),plotDistLoads=false,plotFrequency=round(Int,1e-2/Δt),plotLimits=([-L,L],[-L,L],[-L,L]),save=true,savePath=string(relPath,"/wingDStest_deformation_BLi.gif"),displayProgress=true)
display(anim)

# Range of indices for the last cycle
rangeLastCycle = t -> (t .>= (nCycles-1)*τ) .& (t .<= nCycles*τ)

# Aerodynamic coefficients over time
plt_cl_t = plot(xlim=(0,nCycles), ylim=(0,2.5), ylabel="\$c_l\$")
for i in eachindex(aeroSolvers)
    plot!(plt_cl_t, t[i]/τ, cl[i][ind2plot], label=aeroSolvers[i].name)
end
plt_cm_t = plot(xlim=(0,nCycles), ylim=(-0.4,0.1), ylabel="\$c_m\$")
for i in eachindex(aeroSolvers)
    plot!(plt_cm_t, t[i]/τ, cm[i][ind2plot], label=false)
end
plt_cd_t = plot(xlim=(0,nCycles), ylim=(-0.05,0.65), xlabel="Cycles", ylabel="\$c_d\$")
for i in eachindex(aeroSolvers)
    plot!(plt_cd_t, t[i]/τ, cdrag[i][ind2plot], label=false)
end
plt1 = plot(plt_cl_t, plt_cm_t, plt_cd_t, layout=(3,1))
display(plt1)
savefig(plt1, string(absPath,"/wingDStest_coefficients_time.png"))

# Aerodynamic coefficients over angle of attack
plt_cl_AoA = plot(xlim=(0,25), ylim=(0,2.5), xlabel="\$\\alpha\$ [deg]", ylabel="\$c_l\$")
scatter!(plt_cl_AoA, clRef[1,:], clRef[2,:], c=:black, ms=3, label="Exp", color=:black)
for i in eachindex(aeroSolvers)
    plot!(plt_cl_AoA, α[i][ind2plot][rangeLastCycle(t[i])], cl[i][ind2plot][rangeLastCycle(t[i])], label=aeroSolvers[i].name)
end
plt_cm_AoA = plot(xlim=(0,25), ylim=(-0.4,0.1), xlabel="\$\\alpha\$ [deg]", ylabel="\$c_m\$")
scatter!(plt_cm_AoA, cmRef[1,:], cmRef[2,:], c=:black, ms=3, label=false, color=:black)
for i in eachindex(aeroSolvers)
    plot!(plt_cm_AoA, α[i][ind2plot][rangeLastCycle(t[i])], cm[i][ind2plot][rangeLastCycle(t[i])], label=false)
end
plt_cd_AoA = plot(xlim=(0,25), ylim=(-0.05,0.65), xlabel="\$\\alpha\$ [deg]", ylabel="\$c_d\$")
scatter!(plt_cd_AoA, cdRef[1,:], cdRef[2,:], c=:black, ms=3, label=false, color=:black)
for i in eachindex(aeroSolvers)
    plot!(plt_cd_AoA, α[i][ind2plot][rangeLastCycle(t[i])], cdrag[i][ind2plot][rangeLastCycle(t[i])], label=false)
end
plt2 = plot(plt_cl_AoA, plt_cm_AoA, plt_cd_AoA, layout=@layout [ a [b; c] ])
display(plt2)
savefig(plt2, string(absPath,"/wingDStest_coefficients_AoA.png"))

println("Finished wingDStest.jl")