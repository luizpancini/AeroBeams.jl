using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil
airfoil = deepcopy(flatPlate)

# Derivation method
derivationMethod = AD()

# Root pitch angle [rad]
θ = 5*π/180

# Airspeed
U = 50

# Gust solver
gustLoadsSolver = IndicialGust("Kussner")

# Pazy wing
wing,L,nElem,chord,normSparPos,airfoil,surf = create_Pazy(aeroSolver=aeroSolver,gustLoadsSolver=gustLoadsSolver,airfoil=airfoil,derivationMethod=derivationMethod,p0=[0;-π/2;θ])

# Gust (defined such that it begins at time t0 and last for τ seconds)
t0 = 0.5
τ = 0.5
Ug = -U*1/10
gustLength = U*τ
gustWidth = 2*L
pg = [0;-π/2;0]
c0 = [0; U*t0; 0]
gust = create_DARPAGust(length=gustLength,width=gustWidth,verticalVelocity=Ug,c0=c0,p=pg)

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Set tip loss function at specified airspeed and root angle
surf.tipLossDecayFactor = Pazy_tip_loss_factor(θ*180/π,U)
update_beam!(wing)

# Model
PazyWingDARPAGust = create_Model(name="PazyWingDARPAGust",beams=[wing],BCs=[clamp],gravityVector=[0;0;-9.80665],v_A=[0;U;0],gust=gust)

# Set system solver options
σ0 = 1.0
maxIter = 100
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

# Time variables
Δt = τ/500
tf = 10*τ

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=true, relaxFactor=0.5, Δt=Δt/10)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=PazyWingDARPAGust,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions,adaptableΔt=false)
solve!(problem)
# @profview solve!(problem)
# @time solve!(problem)

# Unpack numerical solution
t = problem.timeVector
tipAoA = [problem.aeroVariablesOverTime[i][nElem].flowAnglesAndRates.αₑ for i in 1:length(t)]
tipOOP = -[problem.nodalStatesOverTime[i][nElem].u_n2[1] for i in 1:length(t)]
tqSpan_cn = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.cn for i in 1:length(t)]
tqSpan_cm = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.cm for i in 1:length(t)]
tqSpan_ct = [problem.aeroVariablesOverTime[i][12].aeroCoefficients.ct for i in 1:length(t)]
tqsχ = [problem.elementalStatesOverTime[i][12].χ for i in 1:length(t)]

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 3
# Tip displacement
plt1 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, tipOOP/L*100, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/PazyWingDARPAGust_1.pdf"))
# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/PazyWingDARPAGust_2.pdf"))
# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/PazyWingDARPAGust_3.pdf"))
# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/PazyWingDARPAGust_4.pdf"))
# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(pwd(),"/test/outputs/figures/PazyWingDARPAGust_5.pdf"))
# Aero states at 3/4-span
nAeroStates = problem.model.elements[1].aero.nTotalAeroStates
colors = get(colorschemes[:rainbow], LinRange(0, 1, nAeroStates))
tqsχ_ = Array{Vector{Float64}}(undef,nAeroStates)
for i in 1:nAeroStates
    tqsχ_[i] = [tqsχ[tt][i] for tt in 1:length(t)]
end
plt6 = plot(xlabel="Time [s]", ylabel="")
for i in 1:nAeroStates
    plot!(t, tqsχ_[i], c=colors[i], lw=lw, label="\$\\chi $(i)\$")
end
display(plt6)
savefig(string(pwd(),"/test/outputs/figures/PazyWingDARPAGust_6.pdf"))

println("Finished PazyWingDARPAGust.jl")
