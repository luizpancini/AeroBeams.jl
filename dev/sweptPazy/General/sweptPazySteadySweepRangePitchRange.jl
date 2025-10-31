using AeroBeams, DelimitedFiles, Plots, ColorSchemes

# --- Configurations ---
# Flag for tip correction
hasTipCorrectionConfig = [false,true,true,true]
# Tip correction function type
tipLossTypeConfig = ["None","Exponential","VLM-undef","VLM-def"]

# Sweep angle range
ΛRange = π/180*vcat(0,10,20,30)

# Root pitch angle range
θRange = π/180*vcat(0,3,5,7)

# Airspeed range
URange = collect(30:1:70)

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = false

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
maxIter = 50
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumIterations=maxIter)

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
tipTwist = Array{Float64}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
tipAoA = Array{Float64}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
AoA = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
cn = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
u1_of_x1 = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
u2_of_x1 = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))

# Sweep configurations
for (c,hasTipCorrection,tipLossType) in zip(1:length(tipLossTypeConfig),hasTipCorrectionConfig,tipLossTypeConfig)
    # Sweep angle of sweep
    for (i,Λ) in enumerate(ΛRange)
        # Sweep pitch angle
        for (j,θ) in enumerate(θRange)
            # Sweep airspeed
            for (k,U) in enumerate(URange)
                # Display progress
                println("Solving for configuration $c, Λ = $(round(Λ*180/π)) deg, θ = $(round(θ*180/π)) deg, U = $U m/s")
                # Model
                model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,Λ=Λ,θ=θ,airspeed=U,g=g,hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepStructuralCorrections=sweepStructuralCorrections)
                # Create and solve problem
                problem[c,i,j,k] = create_SteadyProblem(model=model,systemSolver=NR)
                solve!(problem[c,i,j,k])
                # Outputs
                tipOOP[c,i,j,k] = problem[c,i,j,k].nodalStatesOverσ[end][nElem].u_n2_b[3]
                tip_p = problem[c,i,j,k].nodalStatesOverσ[end][nElem].p_n2_b
                R = first(rotation_tensor_WM(tip_p))
                Δ = R*AeroBeams.a2
                tipTwist[c,i,j,k] = asind(Δ[3])
                tipAoA[c,i,j,k] = problem[c,i,j,k].model.elements[end].aero.flowAnglesAndRates.αₑ*180/π
                cn[c,i,j,k] = [problem[c,i,j,k].model.elements[e].aero.aeroCoefficients.cn for e in 1:nElem]
                AoA[c,i,j,k] = [problem[c,i,j,k].model.elements[e].aero.flowAnglesAndRates.αₑ for e in 1:nElem]
                u1_of_x1[c,i,j,k] = vcat([vcat(problem[c,i,j,k].nodalStatesOverσ[end][e].u_n1[1],problem[c,i,j,k].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
                u2_of_x1[c,i,j,k] = vcat([vcat(problem[c,i,j,k].nodalStatesOverσ[end][e].u_n1[2],problem[c,i,j,k].nodalStatesOverσ[end][e].u_n2[2]) for e in 1:nElem]...)
                u3_of_x1[c,i,j,k] = vcat([vcat(problem[c,i,j,k].nodalStatesOverσ[end][e].u_n1[3],problem[c,i,j,k].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
            end
        end
    end
end

# Undeformed and deformed nodal positions
x1_0 = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
x2_0 = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
x3_0 = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
x1_def = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
x2_def = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
x3_def = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
for c in 1:length(tipLossTypeConfig)
    for i in eachindex(ΛRange)
        for j in eachindex(θRange)
            for k in eachindex(URange)
                x1_0[c,i,j,k] = vcat([vcat(problem[c,i,j,k].model.elements[e].r_n1[1],problem[c,i,j,k].model.elements[e].r_n2[1]) for e in 1:nElem]...)
                x2_0[c,i,j,k] = vcat([vcat(problem[c,i,j,k].model.elements[e].r_n1[2],problem[c,i,j,k].model.elements[e].r_n2[2]) for e in 1:nElem]...)
                x3_0[c,i,j,k] = vcat([vcat(problem[c,i,j,k].model.elements[e].r_n1[3],problem[c,i,j,k].model.elements[e].r_n2[3]) for e in 1:nElem]...)
                x1_def[c,i,j,k] = x1_0[c,i,j,k] .+ u1_of_x1[c,i,j,k]
                x2_def[c,i,j,k] = x2_0[c,i,j,k] .+ u2_of_x1[c,i,j,k]
                x3_def[c,i,j,k] = x3_0[c,i,j,k] .+ u3_of_x1[c,i,j,k]
            end
        end
    end
end

# Load reference data (from AePW4 meetings)
dispΛ0θ7U70_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda0AoA7U70_Sharpy.txt")
dispΛ10θ7U70_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda10AoA7U70_Sharpy.txt")
dispΛ20θ7U70_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda20AoA7U70_Sharpy.txt")
dispΛ30θ7U70_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda30AoA7U70_Sharpy.txt")

# Set paths
relPath = "/dev/sweptPazy/General/outputs/sweptPazySteadySweepRangePitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
xΛ = [ 0.85, 0.95, 0.95, 0.95]
yΛ = [0.65,  0.43, 0.27, 0.15]
colors = cgrad(:rainbow, length(ΛRange), categorical=true)
ls = [:solid, :dash, :dot, :dashdot]
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 3
msw = 0
gr()

# Deformed position at θ=7 deg, U=70 m/s
plt_defPos = plot(xlabel="Normalized horizontal position", ylabel="Normalized vertical position", title=string("\$\\alpha_r=7^\\circ\$",", \$U=70\$ m/s"), xlims=[0,1], ylims=[0,0.8], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topleft)
plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[1], label="No loss")
plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[2], label="Exponential loss")
plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[3], label="VLM - undeformed")
plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[4], label="VLM - deformed")
scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Sharpy (VLM)")
for c in 1:length(tipLossTypeConfig)
    for i in eachindex(ΛRange)
        plot!(x1_def[c,i,end,end]/L, x3_def[c,i,end,end]/L, c=colors[i], lw=lw, ls=ls[c], label=false)
        annotate!([xΛ[i]],[yΛ[i]], text("\$\\Lambda=$(round(Int,ΛRange[i]*180/π)) ^\\circ\$", 12, colors[i]))
        if i == 1
            scatter!(dispΛ0θ7U70_Sharpy[1,:], dispΛ0θ7U70_Sharpy[2,:], c=colors[i], ms=ms, msw=msw, label=false)
        elseif i == 2
            scatter!(dispΛ10θ7U70_Sharpy[1,:], dispΛ10θ7U70_Sharpy[2,:], c=colors[i], ms=ms, msw=msw, label=false)
        elseif i == 3
            scatter!(dispΛ20θ7U70_Sharpy[1,:], dispΛ20θ7U70_Sharpy[2,:], c=colors[i], ms=ms, msw=msw, label=false)
        elseif i == 4
            scatter!(dispΛ30θ7U70_Sharpy[1,:], dispΛ30θ7U70_Sharpy[2,:], c=colors[i], ms=ms, msw=msw, label=false)    
        end
    end
end
display(plt_defPos)
savefig(string(absPath,"/sweptPazySteadySweepRangePitchRange_defPos.pdf"))

println("Finished sweptPazySteadySweepRangePitchRange.jl")