using AeroBeams, ForwardDiff

# Flag to apply a transition function to airspeed profile
applyTransitionFunction = false

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
aeroSolver = Indicial(circulatoryIndicialFunction=:Jose,incompressibleInertialLoads=true)

# Wing surface
chord = 0.1
normSparPos = 1/4
airfoil = deepcopy(FlatPlate)
update_all_airfoil_parameters!(airfoil,Ma)
surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,updateAirfoilParameters=false)

# Number of cycles to run for each case
cycles = 10

# Wing pitch
θ = deg2rad(1)

# Wing beam
L = 1
nElem = 1
wing = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=1e16)],I=[inertia_matrix(ρA=1)],p0=transform_to_WM([0,0,θ],from=:E321),surface=surf)

# BCs
clamp1 = create_BC(name="clamp1",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
timeVaryingFreestream = create_Model(name="timeVaryingFreestream",beams=[wing],BCs=[clamp1,clamp2],atmosphere=atmosphere)

# Initialize outputs
problem = Array{DynamicProblem}(undef,length(λᵤRange))
t = Array{Vector{Float64}}(undef,length(λᵤRange))
tNorm = Array{Vector{Float64}}(undef,length(λᵤRange))
cn = Array{Vector{Float64}}(undef,length(λᵤRange))
cm = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot2 = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot3 = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot2Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
Vdot3Analytical = Array{Vector{Float64}}(undef,length(λᵤRange))
rangeLastCycle = Array{StepRange{Int64,Int64}}(undef,length(λᵤRange))

# Loop normalized airspeed amplitude
for (i,λᵤ) in enumerate(λᵤRange)
    # Update airspeed as a function of time
    ΔU = λᵤ*U₀
    ωᵤ = kᵤ*U₀/(chord/2)
    T = 2π/ωᵤ
    U = applyTransitionFunction ?
        t -> U₀ + ΔU*sin.(ωᵤ*t) * (t/T)^10 / (1 + (t/T)^10) :
        t -> U₀ + ΔU*sin.(ωᵤ*t)
    # Time variables
    tf = cycles*T
    Δt = T/200
    # Update velocity of basis A
    update_v_A!(timeVaryingFreestream,v_A=t->[0,U(t),0])
    # Create and solve problem
    problem[i] = create_DynamicProblem(model=timeVaryingFreestream,finalTime=tf,Δt=Δt)
    solve!(problem[i])
    # Unpack numerical solution
    sol = problem[i].solutionHistory
    states = sol.states
    aeroVars = sol.aeroVariables    
    t[i] = sol.t
    tNorm[i] = t[i]/T
    cn[i] = [aeroVars[j][1].aeroCoefficients.cn for j in eachindex(t[i])]
    cm[i] = [aeroVars[j][1].aeroCoefficients.cm for j in eachindex(t[i])]
    Vdot2[i] = [states[j][1].rates.Vdot[2] for j in eachindex(t[i])]
    Vdot3[i] = [states[j][1].rates.Vdot[3] for j in eachindex(t[i])]
    rangeLastCycle[i] = ceil(Int,(tf-T)/Δt):length(t[i])
    # Analytical solution for relative wind acceleration 
    Udot = t -> ForwardDiff.derivative(U,t)
    Vdot2Analytical[i] = Udot.(t[i]).*cos(θ)
    Vdot3Analytical[i] = -Udot.(t[i]).*sin(θ)
end

# Quasi-steady normal force coefficient
cn_qs = t -> 2π * θ₀

# Load reference data
cnCFDLambda0_2 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestream", "cnCFDLambda0_2.txt"))
cmCFDLambda0_2 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestream", "cmCFDLambda0_2.txt"))
cnCFDLambda0_4 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestream", "cnCFDLambda0_4.txt"))
cmCFDLambda0_4 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestream", "cmCFDLambda0_4.txt"))
cnCFDLambda0_6 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestream", "cnCFDLambda0_6.txt"))
cmCFDLambda0_6 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestream", "cmCFDLambda0_6.txt"))
cnCFDLambda0_8 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestream", "cnCFDLambda0_8.txt"))
cmCFDLambda0_8 = readdlm(pkgdir(AeroBeams,"test","referenceData","timeVaryingFreestream", "cmCFDLambda0_8.txt"))

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
relPath = "/test/outputs/figures/timeVaryingFreestream"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(λᵤRange), categorical=true)
lw = 2
ms = 3
msw = 0
gr()

# Ratio of unsteady to quasi-steady cn over cycle
plt1 = plot(xlabel="\$t/T\$", ylabel="\$c_n/c_{n_{QS}}\$", xlims=[0,1], legend=:topleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cn[i][rangeLastCycle[i]]/(2π*θ), c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cnCFD[i][1,:], cnCFD[i][2,:], c=colors[i], ms=ms, msw=msw, label=false)
end
display(plt1)
savefig(string(absPath,"/timeVaryingFreestream_cn.pdf"))

# cm over cycle
plt2 = plot(xlabel="\$t/T\$", ylabel="\$c_m\$", xlims=[0,1], legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cm[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cmCFD[i][1,:], cmCFD[i][2,:], c=colors[i], ms=ms, msw=msw, label=false)
end
display(plt2)
savefig(string(absPath,"/timeVaryingFreestream_cm.pdf"))

# Relative wind acceleration over cycle
plt31 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_2\$", xlims=[0,1])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Analytical")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot2[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot2Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt32 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot3[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot3Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=msw, label=false)
end
plt3 = plot(plt31,plt32, layout=(2,1))
display(plt3)
savefig(string(absPath,"/timeVaryingFreestream_Vdot.pdf"))

println("Finished timeVaryingFreestream.jl")
