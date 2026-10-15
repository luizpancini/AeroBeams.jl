using AeroBeams

# Root pitch angle range [rad]
θRange = deg2rad.([2.5, 5, 7.5, 10])

# Airspeed range [m/s]
URange = vcat(1e-3,0.5:0.5:30)

# Wingtip configuration
wingtipTab = false

# Discretization
nElementsInner = 16
nElementsOuter = 6
nElementsFFWT = 8

# Tip loss options (tipLossDecayFactor=missing uses custom tip loss function dependent on airspeed and root pitch angle, made to match Healy's data)
hasTipCorrection = true
tipLossDecayFactor = missing

# System solver
systemSolver = create_NLsolver(method=:newton,linesearch=true)

# Initialize model
model = create_HealyLCOFFWT(wingtipTab=wingtipTab,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsOuter=nElementsOuter,nElementsFFWT=nElementsFFWT)

# Fixed properties of the model
nElem = nElementsInner + nElementsOuter + nElementsFFWT
beam = model.beams[1]
L = beam.geometry.ℓ
r_n1 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[1], 
                    beam.elements[e].assemblyData.r_n2[1]] for e in 1:nElem)
r_n2 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[2], 
                    beam.elements[e].assemblyData.r_n2[2]] for e in 1:nElem)
r_n3 = reduce(vcat,[beam.elements[e].assemblyData.r_n1[3], 
                    beam.elements[e].assemblyData.r_n2[3]] for e in 1:nElem)
x1_n = reduce(vcat,[beam.elements[e].geometry.x1_n1[1],
                    beam.elements[e].geometry.x1_n2[1]] for e in 1:nElem)
x1_e = [beam.elements[e].geometry.x1 for e in 1:nElem]

# Initialize outputs
u1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u2 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
u3 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
p1 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
p2 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
p3 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
M2 = Array{Vector{Float64}}(undef,length(θRange),length(URange))
α = Array{Vector{Float64}}(undef,length(θRange),length(URange))
cn = Array{Vector{Float64}}(undef,length(θRange),length(URange))
pHinge = Array{Vector{Float64}}(undef,length(θRange),length(URange))
ϕHinge = Array{Float64}(undef,length(θRange),length(URange))

# Loop root pitch angle
for (i,θ) in enumerate(θRange)
    # Initialize guess solution
    x0 = zeros(0)
    # Loop airspeed
    for (j,U) in enumerate(URange)
        println("Solving for θ = $(round(rad2deg(θ),digits=1)) deg, U = $U m/s")
        # Update model
        model = create_HealyLCOFFWT(wingtipTab=wingtipTab,airspeed=U,pitchAngle=θ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor,nElementsInner=nElementsInner,nElementsOuter=nElementsOuter,nElementsFFWT=nElementsFFWT)
        # Create and solve problem
        problem = create_SteadyProblem(model=model,systemSolver=systemSolver,x0=x0)
        solve!(problem)
        # Get outputs
        state = problem.loadstepHistory.states[end]
        hacState = problem.loadstepHistory.hacStates[end]
        aeroVars = problem.loadstepHistory.aeroVariables[end]
        u1[i,j] = reduce(vcat,[state[e].nodal.u1[1],
                            state[e].nodal.u2[1]] for e in 1:nElem)
        u2[i,j] = reduce(vcat,[state[e].nodal.u1[2],
                            state[e].nodal.u2[2]] for e in 1:nElem)
        u3[i,j] = reduce(vcat,[state[e].nodal.u1[3],
                            state[e].nodal.u2[3]] for e in 1:nElem)
        p1[i,j] = reduce(vcat,[state[e].nodal.p1[1],
                            state[e].nodal.p2[1]] for e in 1:nElem)
        p2[i,j] = reduce(vcat,[state[e].nodal.p1[2],
                            state[e].nodal.p2[2]] for e in 1:nElem)
        p3[i,j] = reduce(vcat,[state[e].nodal.p1[3],
                            state[e].nodal.p2[3]] for e in 1:nElem)
        M2[i,j] = reduce(vcat,[state[e].nodal.M1[2],
                            state[e].nodal.M2[2]] for e in 1:nElem)
        α[i,j] = [rad2deg(aeroVars[e].unsteadyKin.αₑ) for e in 1:nElem]
        cn[i,j] = [aeroVars[e].aeroCoefficients.cn for e in 1:nElem]
        pHinge[i,j] = hacState[1].pH
        ϕHinge[i,j] = rad2deg(hacState[1].ϕ)
        # Set initial guess for next iteration
        x0 = problem.systemArrays.x
    end
end

# Load reference data
aoa_25_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_2.5_exp.txt")
aoa_50_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_5.0_exp.txt")
aoa_75_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_7.5_exp.txt")
aoa_10_exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_10.0_exp.txt")
aoa_25_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_2.5_num.txt")
aoa_50_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_5.0_num.txt")
aoa_75_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_7.5_num.txt")
aoa_10_num = readdlm(pkgdir(AeroBeams)*"/test/referenceData/HealyLCOFFWTsteadyURangeAoARangeCoast/aoa_10.0_num.txt")

## Plots
# Set paths
relPath = "/test/outputs/figures/HealyLCOFFWTsteadyURangeAoARangeCoast"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 12
fs = 16
lfs = 9
lw = 2
ms = 6
msw = 0
gr()
lstyle = [:dash, :solid]

# Select airspeed indices for plots
kPlot = [1, length(URange)]
plt_OOP = plot(xlims=[0,1], xlabel="Normalized spanwise position", ylabel="Normalized OOP position", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:bottomleft)
plot!([NaN], [NaN], lw=lw, ls=lstyle[1], c=:black, label="\$U = $(URange[kPlot[1]])\$ m/s")
plot!([NaN], [NaN], lw=lw, ls=lstyle[2], c=:black, label="\$U = $(URange[kPlot[2]])\$ m/s")
for (j,θ) in enumerate(θRange)
    plot!([NaN], [NaN], lw=lw, c=colors[j], label="\$\\theta = $(round(rad2deg(θ),digits=1)) ^\\circ\$")
end
plt_p1 = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$p_1\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_p2 = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$p_2\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_p3 = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$p_3\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_M2 = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$M_2\$ [N.m]", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_AOA = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$\\alpha\$ [deg]", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plt_cn = plot(xlims=[0,1], xlabel="\$x_1/L\$", ylabel="\$c_n\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (j,θ) in enumerate(θRange)
    for (kk, k) in enumerate(kPlot)
        U = URange[k]
        # OOP displacement
        plot!(plt_OOP, (r_n1.+u1[j,k])/x1_n[end], (r_n3.+u3[j,k])/x1_n[end], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # p1
        plot!(plt_p1, x1_n/x1_n[end], p1[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # p2
        plot!(plt_p2, x1_n/x1_n[end], p2[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # p3
        plot!(plt_p3, x1_n/x1_n[end], p3[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # Bending moment
        plot!(plt_M2, x1_n/x1_n[end], M2[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # Angle of attack
        plot!(plt_AOA, x1_e/x1_n[end], α[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
        # cn
        plot!(plt_cn, x1_e/x1_n[end], cn[j,k], lw=lw, ls=lstyle[kk], c=colors[j], label=false)
    end
end

# Display and save figures
display(plt_OOP)
display(plt_p1)
display(plt_p2)
display(plt_p3)
display(plt_M2)
display(plt_AOA)
display(plt_cn)
savefig(string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_OOP.pdf"))
savefig(string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_p1.pdf"))
savefig(string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_p2.pdf"))
savefig(string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_p3.pdf"))
savefig(string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_M2.pdf"))
savefig(string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_AOA.pdf"))
savefig(string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_cn.pdf"))

# Coast angle vs airspeed for each AoA
plt = plot(xlabel="Airspeed [m/s]", ylabel="Coast angle [deg]", xlims=[0,30], ylims=[-90,60], yticks=-90:15:60, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topleft)
plot!([NaN],[NaN], c=:black, ls=:solid, lw=lw, label="AeroBeams")
scatter!([NaN],[NaN], c=:black, ms=ms, msw=msw, label="Healy (2023) - Exp.")
plot!([NaN],[NaN], c=:black, ls=:dash, lw=lw, label="Healy (2023) - Num.")
for (j,θ) in enumerate(θRange)
    plot!(URange, -ϕHinge[j,:], lw=lw, ls=:solid, c=colors[j], label="\$\\theta = $(round(rad2deg(θ),digits=1)) \\degree\$")
    if θ == deg2rad(2.5)
        scatter!(aoa_25_exp[1,:], aoa_25_exp[2,:], c=colors[j], ms=ms, msw=msw, label=false)
        plot!(aoa_25_num[1,:], aoa_25_num[2,:], c=colors[j], ls=:dash, lw=lw, label=false)
    elseif θ == deg2rad(5.0)
        scatter!(aoa_50_exp[1,:], aoa_50_exp[2,:], c=colors[j], ms=ms, msw=msw, label=false)
        plot!(aoa_50_num[1,:], aoa_50_num[2,:], c=colors[j], ls=:dash, lw=lw, label=false)
    elseif θ == deg2rad(7.5)
        scatter!(aoa_75_exp[1,:], aoa_75_exp[2,:], c=colors[j], ms=ms, msw=msw, label=false)
        plot!(aoa_75_num[1,:], aoa_75_num[2,:], c=colors[j], ls=:dash, lw=lw, label=false)
    elseif θ == deg2rad(10)
        scatter!(aoa_10_exp[1,:], aoa_10_exp[2,:], c=colors[j], ms=ms, msw=msw, label=false)
        plot!(aoa_10_num[1,:], aoa_10_num[2,:], c=colors[j], ls=:dash, lw=lw, label=false)
    end
end
display(plt)
savefig(string(absPath,"/HealyLCOFFWTsteadyURangeAoARangeCoast_coastAngle.pdf"))

println("Finished HealyLCOFFWTsteadyURangeAoARangeCoast.jl")