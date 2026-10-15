using AeroBeams

# Root pitch angle range [rad]
θRange = deg2rad.(vcat(0,3,5,7))

# Airspeed range [m/s]
URange = vcat(1e-3,1:1:100)

# Sweep angle [rad]
Λ = deg2rad(30)

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = true

# Flag for tip correction
hasTipCorrection = true

# Tip correction function type
tipLossType = "VLM-undef"

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(FlatPlate)

# Flag for upright position
upright = false

# Gravity
g = 0

# Initialize model
model,nElem,L,chord,normSparPos = create_Pazy(upright=upright)

# Undeformed nodal and midpoint positions
wing = model.beams[1]
x1_0 = vcat([vcat(wing.elements[e].assemblyData.r_n1[1],wing.elements[e].assemblyData.r_n2[1]) for e in 1:nElem]...)
x2_0 = vcat([vcat(wing.elements[e].assemblyData.r_n1[2],wing.elements[e].assemblyData.r_n2[2]) for e in 1:nElem]...)
x3_0 = -vcat([vcat(wing.elements[e].assemblyData.r_n1[3],wing.elements[e].assemblyData.r_n2[3]) for e in 1:nElem]...)
x1_e = [wing.elements[e].geometry.x1 for e in 1:nElem]

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))
tipOOP = Array{Float64}(undef,length(θRange),length(URange))
tipTwist = Array{Float64}(undef,length(θRange),length(URange))
tipAoA = Array{Float64}(undef,length(θRange),length(URange))
AoA = Array{Vector{Float64}}(undef,length(θRange),length(URange))
cn = Array{Vector{Float64}}(undef,length(θRange),length(URange))
x1_def = Array{Vector{Float64}}(undef,length(θRange),length(URange))
x2_def = Array{Vector{Float64}}(undef,length(θRange),length(URange))
x3_def = Array{Vector{Float64}}(undef,length(θRange),length(URange))

# Sweep pitch angle
for (i,θ) in enumerate(θRange)
    # Model
    model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,sweepAngle=Λ,pitchAngle=θ,gravityVector=[0,0,-g],hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepStructuralCorrections=sweepStructuralCorrections)
    # Initialize guess solution
    x0 = zeros(0)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $(round(rad2deg(θ))) deg, U = $U m/s")
        # Update airspeed
        update_v_A!(model, v_A=[0,U,0])
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=model,x0=x0)
        solve!(problem[i,j])
        # Final states
        state = problem[i,j].loadstepHistory.states[end]
        aeroVars = problem[i,j].loadstepHistory.aeroVariables[end]
        tipElemFinalAero = aeroVars[nElem]
        tipNodeFinalState = state[nElem].nodal
        # Outputs
        tipTwist[i,j] = twist_angle(tipNodeFinalState.p2_b)
        tipOOP[i,j] = tipNodeFinalState.u2_b[3]
        tipAoA[i,j] = rad2deg(tipElemFinalAero.unsteadyKin.αₑ)
        cn[i,j] = [aeroVars[e].aeroCoefficients.cn for e in 1:nElem]
        AoA[i,j] = rad2deg.([aeroVars[e].unsteadyKin.αₑ for e in 1:nElem])
        u1 = vcat([vcat(state[e].nodal.u1[1],state[e].nodal.u2[1]) for e in 1:nElem]...)
        u2 = vcat([vcat(state[e].nodal.u1[2],state[e].nodal.u2[2]) for e in 1:nElem]...)
        u3 = vcat([vcat(state[e].nodal.u1[3],state[e].nodal.u2[3]) for e in 1:nElem]...)
        x1_def[i,j] = x1_0 .+ u1
        x2_def[i,j] = x2_0 .+ u2
        x3_def[i,j] = x3_0 .+ u3
        # Update initial guess for next airspeed
        x0 = problem[i,j].systemArrays.x
    end
end

# Load reference data (from AePW4 meetings)
dispΛ30θ5U60_Nastran = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/dispLambda30AoA5U60_Nastran.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptPazySteadyPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed state of selected condition (AoA = 5 deg and U = 60 m/s)
deformationPlot = plot_steady_deformation(problem[3,61],view=(30,15),plotDistLoads=false,save=true,savePath="/test/outputs/figures/sweptPazySteadyPitchRange/sweptPazySteadyPitchRange_deformation.pdf")
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