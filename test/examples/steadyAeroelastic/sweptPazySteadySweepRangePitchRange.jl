using AeroBeams

# --- Configurations ---
# Flag for tip correction
hasTipCorrectionConfig = [false,true,true,true]
# Tip correction function type
tipLossTypeConfig = ["None","Exponential","VLM-undef","VLM-def"]

# Sweep angle range [rad]
ΛRange = deg2rad.(vcat(0,10,20,30))

# Root pitch angle range [rad]
θRange = deg2rad.(vcat(0,3,5,7))

# Airspeed range [m/s]
URange = collect(30:1:70)

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = false

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(FlatPlate)

# Flag for upright position
upright = false

# Gravity
g = 0

# Geometric properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Undeformed nodal positions (vary only with sweep angle and upright position)
x1_0 = Array{Vector{Float64}}(undef,length(ΛRange))
x2_0 = Array{Vector{Float64}}(undef,length(ΛRange))
x3_0 = Array{Vector{Float64}}(undef,length(ΛRange))
for (i,Λ) in enumerate(ΛRange)
    model,_ = create_Pazy(upright=upright,sweepAngle=Λ)
    wing = model.beams[1]
    x1_0[i] = vcat([vcat(wing.elements[e].assemblyData.r_n1[1],wing.elements[e].assemblyData.r_n2[1]) for e in 1:nElem]...)
    x2_0[i] = vcat([vcat(wing.elements[e].assemblyData.r_n1[2],wing.elements[e].assemblyData.r_n2[2]) for e in 1:nElem]...)
    x3_0[i] = vcat([vcat(wing.elements[e].assemblyData.r_n1[3],wing.elements[e].assemblyData.r_n2[3]) for e in 1:nElem]...)
end

# Initialize outputs
tipOOP = Array{Float64}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
tipTwist = Array{Float64}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
tipAoA = Array{Float64}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
AoA = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
cn = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
x1_def = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
x2_def = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))
x3_def = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange),length(θRange),length(URange))

# Sweep configurations
for (c,hasTipCorrection,tipLossType) in zip(1:length(tipLossTypeConfig),hasTipCorrectionConfig,tipLossTypeConfig)
    # Sweep angle of sweep
    for (i,Λ) in enumerate(ΛRange)
        # Sweep pitch angle
        for (j,θ) in enumerate(θRange)
            # Initialize guess solution
            x0 = zeros(0)
            # Sweep airspeed
            for (k,U) in enumerate(URange)
                # Display progress
                println("Solving for configuration $c, Λ = $(round(Int,rad2deg(Λ))) deg, θ = $(round(Int,rad2deg(θ))) deg, U = $U m/s")
                # Model
                model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,sweepAngle=Λ,pitchAngle=θ,airspeed=U,gravityVector=[0,0,-g],hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepStructuralCorrections=sweepStructuralCorrections)
                # Create and solve problem
                problem = create_SteadyProblem(model=model,x0=x0)
                solve!(problem)
                # Final states
                state = problem.loadstepHistory.states[end]
                aeroVars = problem.loadstepHistory.aeroVariables[end]
                tipElemFinalAero = aeroVars[nElem]
                tipNodeFinalState = state[nElem].nodal
                # Outputs
                tipTwist[c,i,j,k] = twist_angle(tipNodeFinalState.p2_b)
                tipOOP[c,i,j,k] = tipNodeFinalState.u2_b[3]
                tipAoA[c,i,j,k] = rad2deg(tipElemFinalAero.unsteadyKin.αₑ)
                cn[c,i,j,k] = [aeroVars[e].aeroCoefficients.cn for e in 1:nElem]
                AoA[c,i,j,k] = rad2deg.([aeroVars[e].unsteadyKin.αₑ for e in 1:nElem])
                u1 = vcat([vcat(state[e].nodal.u1[1],state[e].nodal.u2[1]) for e in 1:nElem]...)
                u2 = vcat([vcat(state[e].nodal.u1[2],state[e].nodal.u2[2]) for e in 1:nElem]...)
                u3 = vcat([vcat(state[e].nodal.u1[3],state[e].nodal.u2[3]) for e in 1:nElem]...)
                x1_def[c,i,j,k] = x1_0[i] .+ u1
                x2_def[c,i,j,k] = x2_0[i] .+ u2
                x3_def[c,i,j,k] = x3_0[i] .+ u3
                # Set initial guess solution for next iteration
                x0 = problem.systemArrays.x
            end
        end
    end
end

# Load reference data (from AePW4 meetings)
dispΛ0θ7U70_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda0AoA7U70_Sharpy.txt")
dispΛ10θ7U70_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda10AoA7U70_Sharpy.txt")
dispΛ20θ7U70_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda20AoA7U70_Sharpy.txt")
dispΛ30θ7U70_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda30AoA7U70_Sharpy.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptPazySteadySweepRangePitchRange"
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