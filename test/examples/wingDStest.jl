using AeroBeams, ForwardDiff, DelimitedFiles

# Aerodynamic solvers
aeroSolvers = [Indicial(); BLo(); BLi()]

# Frame data
frame = 10_022
airfoil,a₀,a₁,b,k,Ma,U = NASA_frames_loader(frame)

# Pitch profile
ω = k*U/b
τ = 2π/ω
t₀ = -τ/4
θ = t -> a₁*(1+sin(ω*(t+t₀)))
p = t -> 4*tan(θ(t)/4)
pdot = t -> ForwardDiff.derivative(p,t)

# Flag to update airfoil parameters
updateAirfoilParameters = true

# Wing beam properties
L = 10*2b
EIy,GJ = 1e6,1e8
ρA,ρIy,ρIz = 10,1e-2,1e-2
nElem = 20
∞ = 1e12

# Set system solver options
maxIter = 50
NR = create_NewtonRaphson(maximumIterations=maxIter,displayStatus=false,allowAdvanceThroughUnconvergedAeroStates=true)

# Time variables
nCycles = 4
Δt = τ/1000
tf = nCycles*τ

# Elements where output data is to be gathered (root and tip)
elemRangePlot = [1,nElem]

# Initialize outputs
dynProblem = Array{DynamicProblem}(undef,length(aeroSolvers))
t = Array{Vector{Float64}}(undef,length(aeroSolvers))
α = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
cn = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
cm = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
ct = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
cl = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))
cdrag = Array{Vector{Vector{Float64}}}(undef,length(aeroSolvers))

# Loop aero solvers
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Display progress
    println("Solving for $(aeroSolver.name) solver")
    # Aerodynamic surface
    surf = create_AeroSurface(solver=aeroSolver,airfoil=airfoil,c=2*b,normSparPos=0.25,updateAirfoilParameters=updateAirfoilParameters)
    # Wing
    wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=10*EIy)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz)],rotationParametrization="E321",p0=[0;0;a₀-a₁],aeroSurface=surf,pdot0_of_x1=x1->[pdot(0);0;0])
    # BCs
    driver = create_BC(name="driver",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t->p(t),0,0])
    # Model
    wingDStest = create_Model(name="wingDStest",beams=[wing],BCs=[driver],v_A=[0;U;0])
    # Create and solve steady problem for initial solution
    steadyProblem = create_SteadyProblem(model=wingDStest,systemSolver=NR)
    solve!(steadyProblem)
    # Create and solve dynamic problem
    dynProblem[i] = create_DynamicProblem(model=wingDStest,finalTime=tf,Δt=Δt,systemSolver=NR,x0=steadyProblem.x)
    solve!(dynProblem[i])
    # Unpack numerical solution
    t[i] = dynProblem[i].savedTimeVector
    α[i] = [[dynProblem[i].aeroVariablesOverTime[j][elem].flowAnglesAndRates.α for j in 1:length(t[i])] for elem in elemRangePlot]
    cn[i] = [[dynProblem[i].aeroVariablesOverTime[j][elem].aeroCoefficients.cn for j in 1:length(t[i])] for elem in elemRangePlot]
    cm[i] = [[dynProblem[i].aeroVariablesOverTime[j][elem].aeroCoefficients.cm for j in 1:length(t[i])] for elem in elemRangePlot]
    ct[i] = [[dynProblem[i].aeroVariablesOverTime[j][elem].aeroCoefficients.ct for j in 1:length(t[i])] for elem in elemRangePlot]
    cl[i] = [cn[i][e].*cos.(α[i][e]) .+ ct[i][e].*sin.(α[i][e]) for e in eachindex(elemRangePlot)]
    cdrag[i] = [cn[i][e].*sin.(α[i][e]) .- ct[i][e].*cos.(α[i][e]) for e in eachindex(elemRangePlot)]
end

# Load reference data from McAlister et al.
clRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cl.txt")
cmRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cm.txt")
cdRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/NASAframes/"*string(frame)*"_cd.txt")

println("Finished wingDStest.jl")
