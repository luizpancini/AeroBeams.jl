using Plots

# Run the script
include("../examples/OMCgustTests.jl")

# Plot configurations
lw = 2
ms = 5
gr()

# Loop aerodynamic solver
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Loop gust solver
    for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
        # Loop test cases
        for (k,testCase) in enumerate(1:3)
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
            relPath = string("/test/outputs/figures/SEgustTests/OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase)
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

println("Finished OMCgustTestsPlotGenerator.jl")