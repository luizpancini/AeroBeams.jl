using Plots, ColorSchemes

# Run the script
include("../examples/OMCgustTests.jl")

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(aeroSolvers)))
labels = ["QS" "Indicial" "Inflow" "BLi" "BLo"]
linestyles = [:solid :dash :dot :dashdot :dashdotdot]
ts = 10
fs = 16
lw = 2
ms = 5
gr()

# Loop aerodynamic solver
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Loop gust solver
    for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
        # Loop test cases
        for k in tests
            # Aerodynamic solver name
            if typeof(aeroSolver) == QuasiSteady
                aeroSolverName = "QS"
            elseif typeof(aeroSolver) == Indicial
                aeroSolverName = "Indicial"
            elseif typeof(aeroSolver) == Inflow
                aeroSolverName = "Inflow"
            elseif typeof(aeroSolver) == BLi
                aeroSolverName = "BLi"
            elseif typeof(aeroSolver) == BLo
                aeroSolverName = "BLo"
            end
            # Gust indicial solver name
            gustSolverName = gustLoadsSolver.indicialFunctionName
            # Set paths
            relPath = string("/test/outputs/figures/OMCgustTests/",aeroSolverName,"_",gustSolverName,"_test",k)
            absPath = string(pwd(),relPath)
            mkpath(absPath)
            # Lift coefficient increment over time
            plt1 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="\$\\Delta c_l\$")
            plot!(τ[i,j,k], Δcl[i,j,k], color=:black, lw=lw, label="AeroBeams")
            scatter!(ΔclRef[i,j,k][1,:], ΔclRef[i,j,k][2,:], color=:black, ms=ms, label="Mallik & Raveh (2019)")
            display(plt1)
            savefig(string(absPath,"/OMCgustTests_dcl.pdf"))
        end
    end
end

# Set path
absPath = string(pwd(),"/test/outputs/figures/OMCgustTests")
mkpath(absPath)

## Test 2 - with Kussner indicial function
plt21 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="\$\\Delta c_l\$", xlims=[0,30], ylims=[0.0,0.14], tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
scatter!(ΔclRef[1,1,2][1,:], ΔclRef[1,1,2][2,:], color=:black, ms=ms, label="CFD - Mallik & Raveh (2019)")
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(τ[i,1,2], Δcl[i,1,2], color=colors[i], lw=lw, ls=linestyles[i], label=labels[i])
end
display(plt21)
savefig(string(absPath,"/OMCgustTests_test2K.pdf"))

## Test 6 - with Berci and Righi indicial function
plt62 = plot(xlabel="\$\\tau\$ [semichords]", ylabel="\$\\Delta c_l\$", xlims=[0,80], ylims=[-0.05,0.2], tickfont=font(ts), guidefont=font(fs))
scatter!(ΔclRef[1,1,6][1,:], ΔclRef[1,1,6][2,:], color=:black, ms=ms, label="CFD - Mallik & Raveh (2019)")
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(τ[i,2,6], Δcl[i,2,6], color=colors[i], lw=lw, ls=linestyles[i], label=labels[i])
end
display(plt62)
savefig(string(absPath,"/OMCgustTests_test6BR.pdf"))

println("Finished OMCgustTestsPlotGenerator.jl")