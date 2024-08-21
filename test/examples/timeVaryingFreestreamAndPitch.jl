using AeroBeams, LinearAlgebra, ForwardDiff, Plots, ColorSchemes, DelimitedFiles

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Mean airspeed and reduced frequency of airspeed oscillation
kᵤ = 0.2
Ma = 0.3
U₀ = Ma*atmosphere.a

# Wing surface and solver
aeroSolver = Indicial()
airfoil = create_Airfoil(name="flatPlate",Ma=Ma)
chord = 0.1
normSparPos = 1/4
surf = create_AeroSurface(solver=aeroSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,updateAirfoilParameters=false)

# Wing pitch and rates as functions of time
kₚ = 0.2
ωₚ = kₚ*U₀/(chord/2)
τₚ = 2π/ωₚ
θ₀ = 1*π/180
Δθ = 1*π/180
θ = t -> @. θ₀ + Δθ*sin(ωₚ*t)*((t/τₚ)^10/(1+(t/τₚ)^10))
p = t -> 4*tan(θ(t)/4)
Δp = t -> p(t) - p(0)
θdot = iszero(Δθ) ? t -> 0 : t -> ForwardDiff.derivative(θ,t)
θddot = iszero(Δθ) ? t -> 0 : t -> ForwardDiff.derivative(θdot,t)
pdot = iszero(Δθ) ? t -> 0 : t -> ForwardDiff.derivative(p,t)
pddot = iszero(Δθ) ? t -> 0 : t -> ForwardDiff.derivative(pdot,t)

# Wing beam
L = 1
nElem = 1
∞ = 1e12
wing = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],rotationParametrization="E321",p0=[0;0;0],pdot0_of_x1=[pdot(0);0;0],aeroSurface=surf)

# BCs
journal1 = create_BC(name="journal1",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t->p(t),0,0])
journal2 = create_BC(name="journal2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])

# Model
timeVaryingFreestreamAndPitch = create_Model(name="timeVaryingFreestreamAndPitch",beams=[wing],BCs=[journal1,journal2],atmosphere=atmosphere)

# Set normalized airspeed amplitude range and initialize outputs
λᵤRange = collect(0.2:0.2:0.8)
t = Array{Vector{Float64}}(undef,length(λᵤRange))
tNorm = Array{Vector{Float64}}(undef,length(λᵤRange))
αₑ = Array{Vector{Float64}}(undef,length(λᵤRange))
cn = Array{Vector{Float64}}(undef,length(λᵤRange))
cm = Array{Vector{Float64}}(undef,length(λᵤRange))
V2 = Array{Vector{Float64}}(undef,length(λᵤRange))
V3 = Array{Vector{Float64}}(undef,length(λᵤRange))
Ω1 = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot2 = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot3 = Array{Vector{Float64}}(undef,length(λᵤRange))
Ωdot1 = Array{Vector{Float64}}(undef,length(λᵤRange))
V2Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
V3Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot2Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot3Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
rangeLastCycle = Array{StepRange{Int64,Int64}}(undef,length(λᵤRange))

# Loop normalized airspeed amplitude
for (i,λᵤ) in enumerate(λᵤRange)
    # Update airspeed and its rate as functions of time
    ΔU = λᵤ*U₀
    ωᵤ = kᵤ*U₀/(chord/2)
    U = t -> U₀ + ΔU*sin.(ωᵤ*t)
    Udot = t -> ForwardDiff.derivative(U,t)
    # Update velocity of basis A (and update model)
    set_motion_basis_A!(model=timeVaryingFreestreamAndPitch,v_A=t->[0;U(t);0])
    # Time variables
    T = 2π/ωᵤ
    cycles = 20
    tf = cycles*T
    Δt = T/200
    # Initial velocities update options
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=10,displayProgress=true, relaxFactor=0.5, Δt=Δt/1e3)
    # Create and solve problem
    global problem = create_DynamicProblem(model=timeVaryingFreestreamAndPitch,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
    solve!(problem)
    # Unpack numerical solution
    t[i] = problem.timeVector
    tNorm[i] = t[i]/T
    αₑ[i] = [problem.aeroVariablesOverTime[j][1].flowAnglesAndRates.αₑ for j in 1:length(t[i])]
    cn[i] = [problem.aeroVariablesOverTime[j][1].aeroCoefficients.cn for j in 1:length(t[i])]
    cm[i] = [problem.aeroVariablesOverTime[j][1].aeroCoefficients.cm for j in 1:length(t[i])]
    V2[i] = [problem.elementalStatesOverTime[j][1].V[2] for j in 1:length(t[i])]
    V3[i] = [problem.elementalStatesOverTime[j][1].V[3] for j in 1:length(t[i])]
    Ω1[i] = [problem.elementalStatesOverTime[j][1].Ω[1] for j in 1:length(t[i])]
    Vdot2[i] = [problem.elementalStatesRatesOverTime[j][1].Vdot[2] for j in 1:length(t[i])]
    Vdot3[i] = [problem.elementalStatesRatesOverTime[j][1].Vdot[3] for j in 1:length(t[i])]
    Ωdot1[i] = [problem.elementalStatesRatesOverTime[j][1].Ωdot[1] for j in 1:length(t[i])]
    rangeLastCycle[i] = ceil(Int,(tf-T)/Δt):length(t[i])
    # Analytical solution for relative wind speed and acceleration 
    V2Analytical[i] = @. U(t[i])*cos.(θ(t[i]))
    V3Analytical[i] = @. - U(t[i])*sin(θ(t[i]))
    Vdot2Analytical[i] = @. Udot(t[i])*cos.(θ(t[i])) - U(t[i])*sin(θ(t[i]))*θdot(t[i])
    Vdot3Analytical[i] = @. - (Udot(t[i])*sin(θ(t[i])) + U(t[i])*cos(θ(t[i]))*θdot(t[i]))
end

# Analytical solution for angular velocity and acceleration
HT = t -> AeroBeams.tangent_operator_transpose_WM([p(t);0;0])
HT_p = t -> AeroBeams.tangent_tensor_transpose_derivatives_extended_parameters([p(t);0;0])
HT_p1,HT_p2,HT_p3 = t -> HT_p(t)[1], t -> HT_p(t)[2], t -> HT_p(t)[3]
ΩAnalytical = t -> HT(t)*[pdot(t);0;0]
ΩdotAnalytical = t -> AeroBeams.mul3(HT_p1(t),HT_p2(t),HT_p3(t),[pdot(t);0;0])*[pdot(t);0;0] + HT(t)*[pddot(t);0;0]
Ω1Analytical = t -> ΩAnalytical(t)[1]
Ωdot1Analytical = t -> ΩdotAnalytical(t)[1]

# Load reference data JOSE (2006)
cnCFDLambda0_2 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestreamAndPitch/cnCFDLambda0_2.txt"))
cmCFDLambda0_2 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestreamAndPitch/cmCFDLambda0_2.txt"))
cnCFDLambda0_4 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestreamAndPitch/cnCFDLambda0_4.txt"))
cmCFDLambda0_4 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestreamAndPitch/cmCFDLambda0_4.txt"))
cnCFDLambda0_6 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestreamAndPitch/cnCFDLambda0_6.txt"))
cmCFDLambda0_6 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestreamAndPitch/cmCFDLambda0_6.txt"))
cnCFDLambda0_8 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestreamAndPitch/cnCFDLambda0_8.txt"))
cmCFDLambda0_8 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestreamAndPitch/cmCFDLambda0_8.txt"))

cnCFD = Array{Matrix{Float64}}(undef,4)
cmCFD = Array{Matrix{Float64}}(undef,4)
cnCFD[1] = cnCFDLambda0_2
cmCFD[1] = cmCFDLambda0_2
cnCFD[2] = cnCFDLambda0_4
cmCFD[2] = cmCFDLambda0_4
cnCFD[3] = cnCFDLambda0_6
cmCFD[3] = cmCFDLambda0_6
cnCFD[4] = cnCFDLambda0_8
cmCFD[4] = cmCFDLambda0_8

# Plots
# ------------------------------------------------------------------------------
cn_qs = t -> 2π*θ₀
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λᵤRange)))
lw = 2
ms = 3
relPath = "/test/outputs/figures/timeVaryingFreestreamAndPitch"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,showScale=false,plotAeroSurf=false,plotLimits=[(0,L),(-L/2,L/2),(-L/2,L/2)],save=true,savePath=string(relPath,"/timeVaryingFreestreamAndPitch_deformation.gif"),displayProgress=true)
# Ratio of unsteady to quasi-steady cn over cycle
gr()
plt1 = plot(xlabel="\$t/T\$", ylabel="\$c_n/c_{n_{QS}}\$", xlims=[0,1], legend=:topleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cn[i][rangeLastCycle[i]]./cn_qs(t[i][rangeLastCycle[i]]), c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cnCFD[i][1,:], cnCFD[i][2,:], c=colors[i], ms=ms, msw=0, label=false)
end
display(plt1)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_cn.pdf"))
# cm over cycle
plt2 = plot(xlabel="\$t/T\$", ylabel="\$c_m\$", xlims=[0,1], legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cm[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cmCFD[i][1,:], cmCFD[i][2,:], c=colors[i], ms=ms, msw=0, label=false)
end
display(plt2)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_cm.pdf"))
# Relative wind velocity over cycle
plt31 = plot(xlabel="\$t/T\$", ylabel="\$V_2\$", xlims=[0,1])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="Analytical")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], V2[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], V2Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=0, label=false)
end
plt32 = plot(xlabel="\$t/T\$", ylabel="\$V_3\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], V3[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], V3Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=0, label=false)
end
plt33 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_1\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Ω1[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Ω1Analytical.(t[i][rangeLastCycle[i][1:10:end]]), c=colors[i], ms=ms, msw=0, label=false)
end
plt3 = plot(plt31,plt32,plt33, layout=(3,1))
display(plt3)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_velocities.pdf"))
# Relative wind acceleration over cycle
plt41 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_2\$", xlims=[0,1])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="Analytical")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot2[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot2Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=0, label=false)
end
plt42 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot3[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot3Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=0, label=false)
end
plt43 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_1\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Ωdot1[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Ωdot1Analytical.(t[i][rangeLastCycle[i][1:10:end]]), c=colors[i], ms=ms, msw=0, label=false)
end
plt4 = plot(plt41,plt42,plt43, layout=(3,1))
display(plt4)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_accelerations.pdf"))
# Effective angle of attack over cycle
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\alpha_E\$ [deg]", xlims=[0,1], legend=:bottomleft)
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], αₑ[i][rangeLastCycle[i]]*180/π, c=colors[i], lw=lw, label="λ = $λᵤ")
end
display(plt5)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_alpha.pdf"))

println("Finished timeVaryingFreestreamAndPitch.jl")