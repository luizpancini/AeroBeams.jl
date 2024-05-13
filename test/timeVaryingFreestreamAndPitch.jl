using AeroBeams, LinearAlgebra, ForwardDiff, Plots, ColorSchemes, DelimitedFiles

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Mean airspeed and reduced frequency of airspeed oscillation
kᵤ = 0.2
Ma = 0.3
U₀ = Ma*atmosphere.a

# Wing surface
airfoil = create_Airfoil(name="flatPlate",Ma=Ma)
chord = 0.1
normSparPos = 1/4
surf = create_AeroSurface(airfoil=airfoil,c=chord,normSparPos=normSparPos,updateAirfoilParameters=false)

# Wing pitch as a function of time
kₚ = 0.2
ωₚ = kₚ*U₀/(chord/2)
θ₀ = 1*π/180
Δθ = 1*π/180
θ = t -> θ₀ + Δθ*sin.(ωₚ*t)
p = t -> 4*tan(θ(t)/4)
Δp = t -> p(t) - p(0)
θdot = iszero(Δθ) ? t-> 0 : t -> ForwardDiff.derivative(θ,t)
pdot = iszero(Δθ) ? t-> 0 : t -> ForwardDiff.derivative(p,t)

# Wing beam
L = 1
nElem = 1
∞ = 1e12
wing = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],rotationParametrization="E321",p0=[0;0;0],aeroSurface=surf)

# BCs
journal1 = create_BC(name="journal1",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,t->p(t),0,0])
journal2 = create_BC(name="journal2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p2A","p3A"],values=[0,0,0,0,0])

# Model
timeVaryingFreestreamAndPitch = create_Model(name="timeVaryingFreestreamAndPitch",beams=[wing],BCs=[journal1,journal2],atmosphere=atmosphere)

# Set normalized airspeed amplitude range and initialize outputs
λᵤRange = collect(0.2:0.2:0.8)
t = Array{Vector{Float64}}(undef,length(λᵤRange))
tNorm = Array{Vector{Float64}}(undef,length(λᵤRange))
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
Ω1Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot2Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot3Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
Ωdot1Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
rangeLastCycle = Array{StepRange{Int64,Int64}}(undef,length(λᵤRange))

# Loop normalized airspeed amplitude
for (i,λᵤ) in enumerate(λᵤRange)
    # Update airspeed as a function of time
    ΔU = λᵤ*U₀
    ωᵤ = kᵤ*U₀/(chord/2)
    U = t -> U₀ + ΔU*sin.(ωᵤ*t)
    # Update velocity of basis A (and update model)
    set_motion_basis_A!(model=timeVaryingFreestreamAndPitch,v_A=t->[0;U(t);0],ω_A=t->[0;0;0])
    # Time variables
    T = 2π/ωᵤ
    cycles = 10
    tf = cycles*T
    Δt = T/200
    # Initial velocities update options
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=20,displayProgress=true, relaxFactor=0.5, Δt=Δt/1e3)
    # Create and solve problem
    problem = create_DynamicProblem(model=timeVaryingFreestreamAndPitch,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
    solve!(problem)
    # Unpack numerical solution
    t[i] = problem.timeVector
    tNorm[i] = t[i]/T
    cn[i] = [problem.flowVariablesOverTime[j][1].cn for j in 1:length(t[i])]
    cm[i] = [problem.flowVariablesOverTime[j][1].cm for j in 1:length(t[i])]
    V2[i] = [problem.elementalStatesOverTime[j][1].V[2] for j in 1:length(t[i])]
    V3[i] = [problem.elementalStatesOverTime[j][1].V[3] for j in 1:length(t[i])]
    Ω1[i] = [problem.elementalStatesOverTime[j][1].Ω[1] for j in 1:length(t[i])]
    Vdot2[i] = [problem.elementalStatesRatesOverTime[j][1].Vdot[2] for j in 1:length(t[i])]
    Vdot3[i] = [problem.elementalStatesRatesOverTime[j][1].Vdot[3] for j in 1:length(t[i])]
    Ωdot1[i] = [problem.elementalStatesRatesOverTime[j][1].Ωdot[1] for j in 1:length(t[i])]
    rangeLastCycle[i] = ceil(Int,(tf-T)/Δt):length(t[i])
    # Analytical solution for relative wind speed and acceleration 
    Udot = t -> ForwardDiff.derivative(U,t)
    θddot = t -> ForwardDiff.derivative(θdot,t)
    V2Analytical[i] = @. U(t[i])*cos.(θ(t[i]))
    V3Analytical[i] = @. - U(t[i])*sin(θ(t[i]))
    Ω1Analytical[i] = @. θdot(t[i])
    Vdot2Analytical[i] = @. Udot(t[i])*cos.(θ(t[i])) - U(t[i])*sin(θ(t[i]))*θdot(t[i])
    Vdot3Analytical[i] = @. - (Udot(t[i])*sin(θ(t[i])) + U(t[i])*cos(θ(t[i]))*θdot(t[i]))
    Ωdot1Analytical[i] = @. θddot(t[i])
end

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
cn_qs = t -> 2π*θ.(t)
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λᵤRange)))
lw = 2
ms = 3
# Ratio of unsteady to quasi-steady cn over cycle
plt1 = plot(xlabel="\$t/T\$", ylabel="\$c_n/c_{n_{QS}}\$", xlims=[0,1], legend=:topleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cn[i][rangeLastCycle[i]]./cn_qs(t[i][rangeLastCycle[i]]), c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cnCFD[i][1,:], cnCFD[i][2,:], c=colors[i], ms=ms, msw=0, label=false)
end
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/timeVaryingFreestreamAndPitch_cn.pdf"))
# cm over cycle
plt2 = plot(xlabel="\$t/T\$", ylabel="\$c_m\$", xlims=[0,1], legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cm[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cmCFD[i][1,:], cmCFD[i][2,:], c=colors[i], ms=ms, msw=0, label=false)
end
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/timeVaryingFreestreamAndPitch_cm.pdf"))
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
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Ω1Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=0, label=false)
end
plt3 = plot(plt31,plt32,plt33, layout=(3,1))
display(plt3)
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
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Ωdot1Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=0, label=false)
end
plt4 = plot(plt41,plt42,plt43, layout=(3,1))
display(plt4)

println("Finished timeVaryingFreestreamAndPitch.jl")