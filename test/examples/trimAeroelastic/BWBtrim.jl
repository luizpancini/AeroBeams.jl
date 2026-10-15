using AeroBeams

# Airspeed range [m/s]
URange = collect(30:5:160)

# Altitude [m]
h = feet2meters(0e3)

# Attachment spring stiffness
μ = 1e-1

# Aerodynamic solver
aeroSolver = Indicial()

# Model and its beams
BWB = first(create_BWB(aeroSolver=aeroSolver,δElevIsTrimVariable=true,thrustIsTrimVariable=true,altitude=h,withAttachmentSprings=true,μ=μ))

# Root element's global ID
rootElemGlobalID = BWB.beams[3].elements[1].assemblyData.globalID

# System solver
solver = create_NLsolver(trimSolverMethod=:LM,displayStatus=false)

# Initialize outputs
problem = Array{TrimProblem}(undef,length(URange))
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Trimming for airspeed U = $U m/s")
    # Update airspeed on model
    update_v_A!(BWB,v_A=[0,U,0])
    # Set initial guess for the trim problem
    x0 = i==1 ? zeros(0) : problem[i-1].systemArrays.x
    # Create and solve trim problem
    problem[i] = create_TrimProblem(model=BWB,systemSolver=solver,x0=x0)
    solve!(problem[i])
    # Trim results
    trimSol = problem[i].trimSolution
    aeroVars = problem[i].loadstepHistory.aeroVariables
    trimAoA[i] = rad2deg(aeroVars[end][rootElemGlobalID].unsteadyKin.αₑ)
    trimThrust[i] = trimSol.trimLoads[1]
    trimδ[i] = rad2deg(trimSol.trimδ[1])
    println("AoA = $(round(trimAoA[i], digits=2)) deg, T = $(round(trimThrust[i], digits=2)) N, δ = $(round(trimδ[i], digits=2)) deg")
end

# Load reference solution
trimAoARef = readdlm(pkgdir(AeroBeams,"test","referenceData","BWB", "trimAoA.txt"))
trimThrustRef = readdlm(pkgdir(AeroBeams,"test","referenceData","BWB", "trimThrust.txt"))
trimδRef = readdlm(pkgdir(AeroBeams,"test","referenceData","BWB", "trimDelta.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/BWBtrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem[end],view=(45,30),element2centralize=11,save=true,savePath=string(relPath,"/BWBtrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
ms = 4
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[URange[1],URange[end]])
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="UM/NAST")
plot!(URange, trimAoA, c=:black, lw=lw, label=false)
scatter!(trimAoARef[1,:],trimAoARef[2,:], mc=:black, ms=ms, label=false)
display(plt1)
savefig(string(absPath,"/BWBtrim_AoA.pdf"))

# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[URange[1],URange[end]])
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="UM/NAST")
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
scatter!(trimThrustRef[1,:],trimThrustRef[2,:], mc=:black, ms=ms, label=false)
display(plt2)
savefig(string(absPath,"/BWBtrim_thrust.pdf"))

# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[URange[1],URange[end]])
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="UM/NAST")
plot!(URange, trimδ, c=:black, lw=lw, label=false)
scatter!(trimδRef[1,:],trimδRef[2,:], mc=:black, ms=ms, label=false)
display(plt3)
savefig(string(absPath,"/BWBtrim_delta.pdf"))

println("Finished BWBtrim.jl")