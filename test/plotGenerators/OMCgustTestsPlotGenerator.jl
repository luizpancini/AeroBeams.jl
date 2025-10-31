using Plots, ColorSchemes

# Run the script
include("../examples/OMCgustTests.jl")

# Plot configurations
colors = cgrad(:rainbow, length(aeroSolvers), categorical=true)
linestyles = [:solid :dash :dot :dashdot :dashdotdot]
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 5
gr()

# Loop aerodynamic solver
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Loop gust solver
    for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
        # Loop test cases
        for k in tests
            # Set paths
            relPath = string("/test/outputs/figures/OMCgustTests/",aeroSolver.name,"_",gustLoadsSolver.indicialFunctionName,"_test",k)
            absPath = string(pwd(),relPath)
            mkpath(absPath)
            # Lift coefficient increment over time
            plt = plot(xlabel="\$\\tau\$ [semichords]", ylabel="\$\\Delta c_l\$", title="Test $k, $(aeroSolver.name), $(gustLoadsSolver.indicialFunctionName)")
            plot!(τ[i,j,k], Δcl[i,j,k], color=:black, lw=lw, label="AeroBeams")
            scatter!(ΔclRef[i,j,k][1,:], ΔclRef[i,j,k][2,:], color=:black, ms=ms, label="Mallik & Raveh (2019)")
            display(plt)
            savefig(string(absPath,"/OMCgustTests_dcl.pdf"))
        end
    end
end

# Set path
absPath = string(pwd(),"/test/outputs/figures/OMCgustTests")
mkpath(absPath)

# Test 2 - with Kussner indicial function
plt12 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="\$\\Delta c_l\$", xlims=[0,30], ylims=[0.0,0.25], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
scatter!(ΔclRef[1,1,2][1,:], ΔclRef[1,1,2][2,:], color=:black, ms=ms, label="CFD - Mallik & Raveh (2019)")
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(τ[i,1,2], Δcl[i,1,2], color=colors[i], lw=lw, ls=linestyles[i], label=aeroSolver.name)
end
display(plt12)
savefig(string(absPath,"/OMCgustTests_test2K.pdf"))

# Test 6 - with Berci and Righi indicial function
plt26 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="\$\\Delta c_l\$", xlims=[0,80], ylims=[-0.06,0.25], tickfont=font(ts), guidefont=font(fs))
scatter!(ΔclRef[1,1,6][1,:], ΔclRef[1,1,6][2,:], color=:black, ms=ms, label=false)
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(τ[i,2,6], Δcl[i,2,6], color=colors[i], lw=lw, ls=linestyles[i], label=false)
end
display(plt26)
savefig(string(absPath,"/OMCgustTests_test6BR.pdf"))

println("Finished OMCgustTestsPlotGenerator.jl")