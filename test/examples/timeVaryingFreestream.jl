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

# Wing pitch
θ = 1*π/180

# Wing beam
L = 1
nElem = 1
∞ = 1e10
wing = create_Beam(name="beam",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=1)],rotationParametrization="E321",p0=[0;0;θ],aeroSurface=surf)

# BCs
clamp1 = create_BC(name="clamp1",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
clamp2 = create_BC(name="clamp2",beam=wing,node=nElem+1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
timeVaryingFreestream = create_Model(name="timeVaryingFreestream",beams=[wing],BCs=[clamp1,clamp2],atmosphere=atmosphere)

# Set normalized airspeed amplitude range and initialize outputs
λᵤRange = collect(0.2:0.2:0.8)
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
    U = t -> U₀ + ΔU*sin.(ωᵤ*t)
    # Update velocity of basis A (and update model)
    set_motion_basis_A!(model=timeVaryingFreestream,v_A=t->[0;U(t);0])
    # Time variables
    T = 2π/ωᵤ
    cycles = 10
    tf = cycles*T
    Δt = T/200
    # Initial velocities update options
    initialVelocitiesUpdateOptions = InitialVelocitiesUpdateOptions(maxIter=10,displayProgress=true, relaxFactor=0.5, Δt=Δt/1e3)
    # Create and solve problem
    global problem = create_DynamicProblem(model=timeVaryingFreestream,finalTime=tf,Δt=Δt,initialVelocitiesUpdateOptions=initialVelocitiesUpdateOptions)
    solve!(problem)
    # Unpack numerical solution
    t[i] = problem.timeVector
    tNorm[i] = t[i]/T
    cn[i] = [problem.aeroVariablesOverTime[j][1].aeroCoefficients.cn for j in 1:length(t[i])]
    cm[i] = [problem.aeroVariablesOverTime[j][1].aeroCoefficients.cm for j in 1:length(t[i])]
    Vdot2[i] = [problem.elementalStatesRatesOverTime[j][1].Vdot[2] for j in 1:length(t[i])]
    Vdot3[i] = [problem.elementalStatesRatesOverTime[j][1].Vdot[3] for j in 1:length(t[i])]
    rangeLastCycle[i] = ceil(Int,(tf-T)/Δt):length(t[i])
    # Analytical solution for relative wind acceleration 
    Udot = t -> ForwardDiff.derivative(U,t)
    Vdot2Analytical[i] = Udot.(t[i]).*cos(θ)
    Vdot3Analytical[i] = -Udot.(t[i]).*sin(θ)
end

# Load reference data JOSE (2006)
cnCFDLambda0_2 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestream/cnCFDLambda0_2.txt"))
cmCFDLambda0_2 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestream/cmCFDLambda0_2.txt"))
cnCFDLambda0_4 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestream/cnCFDLambda0_4.txt"))
cmCFDLambda0_4 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestream/cmCFDLambda0_4.txt"))
cnCFDLambda0_6 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestream/cnCFDLambda0_6.txt"))
cmCFDLambda0_6 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestream/cmCFDLambda0_6.txt"))
cnCFDLambda0_8 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestream/cnCFDLambda0_8.txt"))
cmCFDLambda0_8 = readdlm(string(pwd(),"/test/referenceData/timeVaryingFreestream/cmCFDLambda0_8.txt"))

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
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λᵤRange)))
lw = 2
ms = 3
relPath = "/test/outputs/figures/timeVaryingFreestream"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,showScale=false,plotAeroSurf=false,plotLimits=[(0,L),(-L/2,L/2),(-L/2,L/2)],save=true,savePath=string(relPath,"/timeVaryingFreestream_deformation.gif"),displayProgress=true)
# Ratio of unsteady to quasi-steady cn over cycle
gr()
plt1 = plot(xlabel="\$t/T\$", ylabel="\$c_n/c_{n_{QS}}\$", xlims=[0,1], legend=:topleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cn[i][rangeLastCycle[i]]/(2π*θ), c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cnCFD[i][1,:], cnCFD[i][2,:], c=colors[i], ms=ms, msw=0, label=false)
end
display(plt1)
savefig(string(absPath,"/timeVaryingFreestream_cn.pdf"))
# cm over cycle
plt2 = plot(xlabel="\$t/T\$", ylabel="\$c_m\$", xlims=[0,1], legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="CFD - Jose (2006)")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], cm[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(cmCFD[i][1,:], cmCFD[i][2,:], c=colors[i], ms=ms, msw=0, label=false)
end
display(plt2)
savefig(string(absPath,"/timeVaryingFreestream_cm.pdf"))
# Relative wind acceleration over cycle
plt31 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_2\$", xlims=[0,1])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="Analytical")
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot2[i][rangeLastCycle[i]], c=colors[i], lw=lw, label="λ = $λᵤ")
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot2Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=0, label=false)
end
plt32 = plot(xlabel="\$t/T\$", ylabel="\$\\dot{V}_3\$", xlims=[0,1])
for (i,λᵤ) in enumerate(λᵤRange)
    plot!(tNorm[i][rangeLastCycle[i]].-tNorm[i][rangeLastCycle[i][1]], Vdot3[i][rangeLastCycle[i]], c=colors[i], lw=lw, label=false)
    scatter!(tNorm[i][rangeLastCycle[i][1:10:end]].-tNorm[i][rangeLastCycle[i][1]], Vdot3Analytical[i][rangeLastCycle[i][1:10:end]], c=colors[i], ms=ms, msw=0, label=false)
end
plt3 = plot(plt31,plt32, layout=(2,1))
display(plt3)
savefig(string(absPath,"/timeVaryingFreestream_Vdot.pdf"))

println("Finished timeVaryingFreestream.jl")