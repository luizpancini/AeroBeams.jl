using AeroBeams, DelimitedFiles, Plots, ColorSchemes

# Sweep angle [rad]
Λ = 30*π/180

# Root pitch angle range
θRange = π/180*vcat(0,3,5,7)

# Airspeed range
URange = collect(1:1:100)

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-undef"

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(flatPlate)

# Flag for upright position
upright = false

# Gravity
g = 0

# System solver
σ0 = 0.5
σstep = 0.5
maxIter = 50
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep,maximumIterations=maxIter)

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(θRange),length(URange))
tipTwist = Array{Float64}(undef,length(θRange),length(URange))
tipAoA = Array{Float64}(undef,length(θRange),length(URange))
AoA = Array{Vector{Float64}}(undef,length(θRange),length(URange))
cn = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u1_of_x1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u2_of_x1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(θ*180/π)) deg, U = $U m/s")
        # Model
        model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,Λ=Λ,θ=θ,airspeed=U,g=g,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepStructuralCorrections=sweepStructuralCorrections)
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=model,systemSolver=NR)
        solve!(problem[i,j])
        # Outputs
        tipOOP[i,j] = problem[i,j].nodalStatesOverσ[end][nElem].u_n2_b[3]
        tip_p = problem[i,j].nodalStatesOverσ[end][nElem].p_n2_b
        R = first(rotation_tensor_WM(tip_p))
        Δ = R*AeroBeams.a2
        tipTwist[i,j] = asind(Δ[3])
        tipAoA[i,j] = problem[i,j].model.elements[end].aero.flowAnglesAndRates.αₑ*180/π
        cn[i,j] = [problem[i,j].model.elements[e].aero.aeroCoefficients.cn for e in 1:nElem]
        AoA[i,j] = [problem[i,j].model.elements[e].aero.flowAnglesAndRates.αₑ for e in 1:nElem]
        u1_of_x1[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[1],problem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
        u2_of_x1[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[2],problem[i,j].nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
        u3_of_x1[i,j] = vcat([vcat(problem[i,j].nodalStatesOverσ[end][e].u_n1[3],problem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
    end
end

# Undeformed nodal and elemental positions
x1_0 = vcat([vcat(problem[1].model.elements[e].r_n1[1],problem[1].model.elements[e].r_n2[1]) for e in 1:nElem]...)
x2_0 = vcat([vcat(problem[1].model.elements[e].r_n1[2],problem[1].model.elements[e].r_n2[2]) for e in 1:nElem]...)
x3_0 = vcat([vcat(problem[1].model.elements[e].r_n1[3],problem[1].model.elements[e].r_n2[3]) for e in 1:nElem]...)
x1_e = getindex.(getfield.(problem[1].model.elements, :x1), 1)

# Deformed nodal positions
x1_def = [x1_0 .+ u1_of_x1[i,j] for i in eachindex(θRange), j in eachindex(URange)]
x2_def = [x2_0 .+ u2_of_x1[i,j] for i in eachindex(θRange), j in eachindex(URange)]
x3_def = [x3_0 .+ u3_of_x1[i,j] for i in eachindex(θRange), j in eachindex(URange)]

# Load reference data (from AePW4 meetings)
dispΛ30θ5U60_Nastran = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda30AoA5U60_Nastran.txt")

# Set paths
relPath = "/dev/sweptPazy/General/outputs/sweptPazySteadyPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed state of selected condition (AoA = 5 deg and U = 60 m/s)
deformationPlot = plot_steady_deformation(problem[3,60],view=(30,15),plotDistLoads=false,save=true,savePath="/dev/sweptPazy/General/outputs/sweptPazySteadyPitchRange/sweptPazySteadyPitchRange_deformation.pdf")
display(deformationPlot)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0
gr()

# Deformed position at θ=5 deg, U=60 m/s
x = upright ? x3_def[3,60] : x1_def[3,60]
y = upright ? -x1_def[3,60] : x3_def[3,60]
plt_defPos = plot(xlabel="\$x_1\$ [m]", ylabel="\$x_3\$ [m]", xlims=[0,0.5], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!(x, y, c=:black, lw=lw, ls=:solid, label="AeroBeams")
scatter!(dispΛ30θ5U60_Nastran[1,:],dispΛ30θ5U60_Nastran[2,:], c=:black, ms=ms, msw=msw, label="Technion (Nastran)")
display(plt_defPos)
savefig(string(absPath,"/sweptPazySteadyPitchRange_defPos.pdf"))

# Tip OOP displacement vs. airspeed
plt_tipOOP = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,100], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
for (i,θ) in enumerate(θRange)
    plot!(URange, tipOOP[i,:]/L*100, c=colors[i], lw=lw, ls=:solid, label="\$\\theta = $(round(Int,θ*180/pi)) ^\\circ\$")
end
display(plt_tipOOP)
savefig(string(absPath,"/sweptPazySteadyPitchRange_tipOOP.pdf"))

# Tip twist vs. airspeed
plt_tipTwist = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,100], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tipTwist[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_tipTwist)
savefig(string(absPath,"/sweptPazySteadyPitchRange_tipTwist.pdf"))

# Tip AoA vs. airspeed
plt_tipAOA = plot(xlabel="Airspeed [m/s]", ylabel="Tip angle of attack [deg]", xlims=[0,100], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tipAoA[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_tipAOA)
savefig(string(absPath,"/sweptPazySteadyPitchRange_tipAoA.pdf"))

# Normal force coefficient slope over span at selected condition (AoA = 5 deg and U = 60 m/s)
plt_cna = plot(xlabel="Normalized span", ylabel="\$c_{n_\\alpha}\$ [1/rad]", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
plot!(x1_e/L, cn[3,60]./AoA[3,60], c=:black, lw=lw, label=false)
display(plt_cna)
savefig(string(absPath,"/sweptPazySteadyPitchRange_cna.pdf"))

println("Finished sweptPazySteadyPitchRange.jl")