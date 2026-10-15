using AeroBeams, ForwardDiff

# Flag to apply a transition function to airspeed and pitch profiles
applyTransitionFunction = true

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)

# Normalized airspeed amplitude range
λᵤRange = collect(0.2:0.2:0.8)

# Mean airspeed and reduced frequency of airspeed oscillation
kᵤ = 0.2
Ma = 0.3
U₀ = Ma*atmosphere.a

# Aerodynamic solver
aeroSolver = Indicial(circulatoryIndicialFunction=:Jose,incompressibleInertialLoads=false)

# Wing surface and solver
chord = 0.1
normSparPos = 1/4
airfoil = deepcopy(FlatPlate)
update_all_airfoil_parameters!(airfoil,Ma)
surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,updateAirfoilParameters=false)

# Number of cycles to run for each case
cycles = 5

# Wing pitch and rates as functions of time
kₚ = 0.2
ωₚ = kₚ*U₀/(chord/2)
τₚ = 2π/ωₚ
θ₀ = deg2rad(1)
Δθ = deg2rad(1)
θ = applyTransitionFunction ?
    t -> θ₀ + Δθ*sin(ωₚ*t) *  (t/τₚ)^10 / (1 + (t/τₚ)^10) :
    t -> θ₀ + Δθ*sin(ωₚ*t)
p = t -> [4*tan(θ(t)/4), 0, 0]
θdot = iszero(Δθ) ? t -> 0 : t -> ForwardDiff.derivative(θ,t)
θddot = iszero(Δθ) ? t -> 0 : t -> ForwardDiff.derivative(θdot,t)
pdot = iszero(Δθ) ? t -> 0 : t -> ForwardDiff.derivative(p,t)
pddot = iszero(Δθ) ? t -> 0 : t -> ForwardDiff.derivative(pdot,t)

# Wing beam
L = 1
nElem = 1
wing = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e16)],I=[inertia_matrix(ρA=1)],surface=surf)

# Set initial conditions
add_initial_conditions_to_beam!(beam=wing,conditions=Dict(:pdot0_of_x1=>x1->pdot(0)))

# BCs
journal1 = create_BC(name="journal1",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>t->p(t)[1],:p2A=>t->p(t)[2],:p3A=>t->p(t)[3]))
journal2 = create_BC(name="journal2",beam=wing,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p2A=>0,:p3A=>0))

# Model
timeVaryingFreestreamAndPitch = create_Model(name="timeVaryingFreestreamAndPitch",beams=[wing],BCs=[journal1,journal2],atmosphere=atmosphere)

# Initialize outputs
problem = Array{DynamicProblem}(undef,length(λᵤRange))
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
    T = 2π/ωᵤ
    U = applyTransitionFunction ? 
        t -> U₀ + ΔU*sin.(ωᵤ*t) * (t/T)^10 / (1 + (t/T)^10) :
        t -> U₀ + ΔU*sin.(ωᵤ*t)
    Udot = t -> ForwardDiff.derivative(U,t)
    # Time variables
    tf = cycles*T
    Δt = T/200
    # Update velocity of basis A
    update_v_A!(timeVaryingFreestreamAndPitch,v_A=t->[0,U(t),0])
    # Create and solve problem
    problem[i] = create_DynamicProblem(model=timeVaryingFreestreamAndPitch,finalTime=tf,Δt=Δt)
    solve!(problem[i])
    # Unpack numerical solution
    sol = problem[i].solutionHistory
    states = sol.states
    aeroVars = sol.aeroVariables    
    t[i] = sol.t
    tNorm[i] = t[i]/T
    αₑ[i] = [aeroVars[j][1].unsteadyKin.αₑ for j in eachindex(t[i])]
    cn[i] = [aeroVars[j][1].aeroCoefficients.cn for j in eachindex(t[i])]
    cm[i] = [aeroVars[j][1].aeroCoefficients.cm for j in eachindex(t[i])]
    V2[i] = [states[j][1].midpoint.V[2] for j in eachindex(t[i])]
    V3[i] = [states[j][1].midpoint.V[3] for j in eachindex(t[i])]
    Ω1[i] = [states[j][1].midpoint.Ω[1] for j in eachindex(t[i])]
    Vdot2[i] = [states[j][1].rates.Vdot[2] for j in eachindex(t[i])]
    Vdot3[i] = [states[j][1].rates.Vdot[3] for j in eachindex(t[i])]
    Ωdot1[i] = [states[j][1].rates.Ωdot[1] for j in eachindex(t[i])]
    rangeLastCycle[i] = ceil(Int,(tf-T)/Δt):length(t[i])
    # Analytical solution for relative wind speed and acceleration 
    V2Analytical[i] = @. U(t[i])*cos.(θ(t[i]))
    V3Analytical[i] = @. - U(t[i])*sin(θ(t[i]))
    Vdot2Analytical[i] = @. Udot(t[i])*cos.(θ(t[i])) - U(t[i])*sin(θ(t[i]))*θdot(t[i])
    Vdot3Analytical[i] = @. - (Udot(t[i])*sin(θ(t[i])) + U(t[i])*cos(θ(t[i]))*θdot(t[i]))
end

# Quasi-steady normal force coefficient
cn_qs = t -> 2π * θ₀

# Analytical solution for angular velocity and acceleration
HT = t -> AeroBeams.HT_WM(p(t))
HT_p = t -> ForwardDiff.jacobian(p -> AeroBeams.HT_WM(p), p(t))
HT_p1,HT_p2,HT_p3 = t -> HT_p(t)[1:3,:], t -> HT_p(t)[4:6,:], t -> HT_p(t)[7:9,:]
ΩAnalytical = t -> HT(t)*pdot(t)
ΩdotAnalytical = t -> AeroBeams.mul3(HT_p1(t),HT_p2(t),HT_p3(t),pdot(t))*pdot(t) + HT(t)*pddot(t)
Ω1Analytical = t -> ΩAnalytical(t)[1]
Ωdot1Analytical = t -> ΩdotAnalytical(t)[1]

# Load reference data
cnCFDLambda0_2 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestreamAndPitch", "cnCFDLambda0_2.txt"))
cmCFDLambda0_2 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestreamAndPitch", "cmCFDLambda0_2.txt"))
cnCFDLambda0_4 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestreamAndPitch", "cnCFDLambda0_4.txt"))
cmCFDLambda0_4 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestreamAndPitch", "cmCFDLambda0_4.txt"))
cnCFDLambda0_6 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestreamAndPitch", "cnCFDLambda0_6.txt"))
cmCFDLambda0_6 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestreamAndPitch", "cmCFDLambda0_6.txt"))
cnCFDLambda0_8 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestreamAndPitch", "cnCFDLambda0_8.txt"))
cmCFDLambda0_8 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestreamAndPitch", "cmCFDLambda0_8.txt"))

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

## Plots
# Set paths
relPath = "/test/outputs/figures/timeVaryingFreestreamAndPitch"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(λᵤRange), categorical=true)
ts = 10
fs = 16
lw = 2
ms = 3
msw = 0
gr()

# Ratio of unsteady to quasi-steady cn over cycle
plt1 = plot(xlabel="\$t/T\$", ylabel="\$c_n/c_{n_{QS}}\$", xlims=[0,1], legend=:topleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="CFD - Jose et al. (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cn[i][rangeLastCycle[i]]./cn_qs(t[i][rangeLastCycle[i]]), c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cnCFD[i][1,:], cnCFD[i][2,:], c=colors[i], ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_cn.pdf"))

# cm over cycle
plt2 = plot(xlabel="\$t/T\$", ylabel="\$c_m\$", xlims=[0,1], legend=:bottomleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
plot!([NaN], [NaN], c=:black, lw=lw, label=false)
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label=false)
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cm[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(cmCFD[i][1,:], cmCFD[i][2,:], c=colors[i], ms=ms, msw=msw, label=false)
end
display(plt2)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_cm.pdf"))

# Relative wind velocity over cycle
plt31 = plot(xlabel="\$t/T\$", ylabel="\$V_2^\\star\$", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Analytical")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], V2[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], V2Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt32 = plot(xlabel="\$t/T\$", ylabel="\$V_3^\\star\$", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], V3[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], V3Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt33 = plot(xlabel="\$t/T\$", ylabel="\$\\Omega_1^\\star\$", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Ω1[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Ω1Analytical.(t[i][rangeLastCycle[i][1:10:end]]), c=colors[i], ms=ms, msw=msw, label=false)
end
plt3 = plot(plt31,plt32,plt33, layout=(3,1))
display(plt3)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_velocities.pdf"))

# Relative wind acceleration over cycle
plt41 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_2^\\star\$", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
plot!([NaN], [NaN], c=:black, lw=lw, label=false)
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label=false)
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot2[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot2Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt42 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3^\\star\$", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot3[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot3Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt43 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{\\Omega}_1^\\star\$", xlims=[0,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Ωdot1[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Ωdot1Analytical.(t[i][rangeLastCycle[i][1:10:end]]), c=colors[i], ms=ms, msw=msw, label=false)
end
plt4 = plot(plt41,plt42,plt43, layout=(3,1))
display(plt4)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_accelerations.pdf"))

# Effective angle of attack over cycle
plt5 = plot(xlabel="\$t/T\$", ylabel="\$\\alpha_E\$ [deg]", xlims=[0,1], legend=:bottomleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], αₑ[i][rangeLastCycle[i]]*180/π, c=colors[i], lw=lw, label="λ = $λᵤ")
end
display(plt5)
savefig(string(absPath,"/timeVaryingFreestreamAndPitch_alpha.pdf"))

println("Finished timeVaryingFreestreamAndPitch.jl")
