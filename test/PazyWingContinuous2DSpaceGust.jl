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

# Gust (defined such that it begins at time t0 and lasts for τ seconds)
spectrum = "vK"
t0 = 0.5
τ = 1.0
σ = U/15
c0 = [0;t0*U;0]
pg = [0;-π/2;0]
gust = create_Continuous2DSpaceGust(spectrum=spectrum,length=τ*U,width=2*L,Nx=101,Ny=101,σ=σ,c0=c0,p=pg)

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Set tip loss function at specified airspeed and root angle
surf.tipLossDecayFactor = Pazy_tip_loss_factor(θ*180/π,U)
update_beam!(wing)

# Model
PazyWingContinuous2DSpaceGust = create_Model(name="PazyWingContinuous2DSpaceGust",beams=[wing],BCs=[clamp],gravityVector=[0;0;-9.80665],v_A=[0;U;0],gust=gust)

# Set system solver options
σ0 = 1.0
maxIter = 100
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter,displayStatus=false,alwaysUpdateJacobian=false,minConvRateAeroJacUpdate=1.2,minConvRateJacUpdate=1.2)

# Time variables
Δt = τ/500
tf = 5*t0 + τ

# Initial velocities update options
initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=2,tol=1e-8, displayProgress=true, relaxFactor=0.5, Δt=Δt/10)

# Create and solve dynamic problem
problem = create_DynamicProblem(model=PazyWingContinuous2DSpaceGust,finalTime=tf,Δt=Δt,systemSolver=NR,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions,adaptableΔt=false)
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
savefig(string(pwd(),"/test/outputs/figures/PazyWingContinuousSpaceGust_1.pdf"))
# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/PazyWingContinuousSpaceGust_2.pdf"))
# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/PazyWingContinuousSpaceGust_3.pdf"))
# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/PazyWingContinuousSpaceGust_4.pdf"))
# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(pwd(),"/test/outputs/figures/PazyWingContinuousSpaceGust_5.pdf"))
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
savefig(string(pwd(),"/test/outputs/figures/PazyWingContinuousSpaceGust_6.pdf"))
# "Vertical" gust velocity
W = t -> t0<t<t0+τ ? gust.W([0; t*U; 0]) : 0
plt7 = plot(xlabel="Time [s]", ylabel="Gust velocity [m/s]")
plot!(t, W.(t), color=:black, lw=lw, label=false)
display(plt7)
savefig(string(pwd(),"/test/outputs/figures/PazyWingContinuousSpaceGust_7.pdf"))

println("Finished PazyWingContinuous2DSpaceGust.jl")
