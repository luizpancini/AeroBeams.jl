using AeroBeams, LinearAlgebra, ForwardDiff, Plots, ColorSchemes, DelimitedFiles, BenchmarkTools

# Aerodynamic solver
aeroSolver = BLi()

# Derivation method
derivationMethod = AD()

# Frame data
a₀ = 0.20944
a₁ = 0.17279
b = 0.305
k = 0.098
U = 102.34

# Pitch profile
ω = k*U/b
τ = 2π/ω
t₀ = -τ/4
θ = t -> a₁*(1+sin(ω*(t+t₀)))
p = t -> 4*tan(θ(t)/4)
pdot = t -> ForwardDiff.derivative(p,t)

# Aerodynamic surface
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=deepcopy(NACA0012),c=2*b,normSparPos=0.25,updateAirfoilParameters=false)

# Wing
L = 1.0
EI,GJ = 2e10,2e10
nElem = 1
∞ = 1e12
wing = create_Beam(name="wing",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EI,EIz=EI)],I=[inertia_matrix(ρA=10,ρIy=1e-2,ρIz=1e-2)],rotationParametrization="E321",p0=[0;0;a₀-a₁],aeroSurface=surf,pdot0_of_x1=x1->[pdot(0); 0.0; 0.0])

# BCs
driver = create_BC(name="driver",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t -> p(t),0,0])
journal = create_BC(name="journal",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])

# Model
DSModelTest = create_Model(name="DSModelTest",beams=[wing],BCs=[driver,journal],v_A=[0;U;0])

# Set system solver options
σ0 = 1.0
maxIter = 20
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=true)

# Time variables
nCycles = 4
Δt = τ/500
tf = nCycles*τ

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=true, relaxFactor=0.5, Δt=Δt/1e3)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=DSModelTest,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions,adaptableΔt=false,minΔt=Δt/2^2)
# solve!(problem)
@time solve!(problem)
# @profview solve!(problem)

# Unpack numerical solution
t = problem.timeVector
α = [problem.aeroVariablesOverTime[i][1].flowAnglesAndRates.α for i in 1:length(t)]
cn = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cn for i in 1:length(t)]
cm = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.cm for i in 1:length(t)]
ct = [problem.aeroVariablesOverTime[i][1].aeroCoefficients.ct for i in 1:length(t)]
cl = @. cn*cos(α) + ct*sin(α)
cd = @. cn*sin(α) - ct*cos(α)
Vdot₃ = [problem.elementalStatesRatesOverTime[i][1].Vdot[3] for i in 1:length(t)]
χ = [problem.elementalStatesOverTime[i][1].χ for i in 1:length(t)]
χdot = [problem.elementalStatesRatesOverTime[i][1].χdot for i in 1:length(t)]

# Load reference data from McAlister et al (frame 10022)
clRef = readdlm(string(pwd(),"/test/referenceData/DSModelTest/cl.txt"))
cmRef = readdlm(string(pwd(),"/test/referenceData/DSModelTest/cm.txt"))
cdRef = readdlm(string(pwd(),"/test/referenceData/DSModelTest/cd.txt"))

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 3
# Pitch angle
plt1 = plot(xlabel="Time [s]", ylabel="Pitch angle [deg]")
plot!(t, α*180/π, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/DSModelTest_a.pdf"))
# Normal relative wind acceleration
plt0 = plot(xlabel="Time [s]", ylabel="\$\\dot{V}_3\$ [m/s^2]")
plot!(t, Vdot₃, color=:black, lw=lw, label=false)
display(plt0)
# cn vs time
plt2 = plot(xlabel="Time [s]", ylabel="\$c_n\$")
plot!(t, cn, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/DSModelTest_cnt.pdf"))
# cm vs time
plt3 = plot(xlabel="Time [s]", ylabel="\$c_m\$")
plot!(t, cm, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/DSModelTest_cmt.pdf"))
# ct vs time
plt4 = plot(xlabel="Time [s]", ylabel="\$c_t\$")
plot!(t, ct, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/DSModelTest_ctt.pdf"))
# cl vs α
plt7 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_l\$")
plot!(α*180/π, cl, color=:black, lw=lw, label="AeroBeams")
scatter!(clRef[1,:], clRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
display(plt7)
savefig(string(pwd(),"/test/outputs/figures/DSModelTest_cla.pdf"))
# cm vs α
plt8 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_m\$")
plot!(α*180/π, cm, color=:black, lw=lw, label="AeroBeams")
scatter!(cmRef[1,:], cmRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
display(plt8)
savefig(string(pwd(),"/test/outputs/figures/DSModelTest_cma.pdf"))
# cd vs α
plt9 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_d\$")
plot!(α*180/π, cd, color=:black, lw=lw, label="AeroBeams")
scatter!(cdRef[1,:], cdRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
display(plt9)
savefig(string(pwd(),"/test/outputs/figures/DSModelTest_cda.pdf"))
# Aero states at 3/4-span
nTotalAeroStates = problem.model.elements[1].aero.nTotalAeroStates
colors = get(colorschemes[:rainbow], LinRange(0, 1, nTotalAeroStates))
χ_ = Array{Vector{Float64}}(undef,nTotalAeroStates)
for i in 1:nTotalAeroStates
    χ_[i] = [χ[tt][i] for tt in 1:length(t)]
end
plt6 = plot(xlabel="Time [s]", ylabel="")
for i in 1:nTotalAeroStates
    plot!(t, χ_[i], c=colors[i], lw=lw, label="\$\\chi $(i)\$")
end
display(plt6)
savefig(string(pwd(),"/test/outputs/figures/DSModelTest_6.pdf"))
# Aero states' rates at 3/4-span
χdot_ = Array{Vector{Float64}}(undef,nTotalAeroStates)
for i in 1:nTotalAeroStates
    χdot_[i] = [χdot[tt][i] for tt in 1:length(t)]
end
plt7 = plot(xlabel="Time [s]", ylabel="")
for i in 1:nTotalAeroStates
    plot!(t, χdot_[i], c=colors[i], lw=lw, label="\$\\dot{\\chi} $(i)\$")
end
display(plt7)
savefig(string(pwd(),"/test/outputs/figures/DSModelTest_7.pdf"))

println("Finished DSModelTest.jl")
