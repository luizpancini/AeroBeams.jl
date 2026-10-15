using AeroBeams

# Altitude [m]
h = 20e3

# Airspeed range [m/s]
URange = collect(20:1:35)

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1e0

# Option for stabilizers
stabilizersAero = true
includeVS = true

# Parasite drag coefficients
wingCd0 = 1e-2
stabsCd0 = 1e-2

# Discretization
nElemWing = 20
nElemTailBoom = 5
nElemHorzStabilizer = 4

# Attachment spring stiffness
μ = 1e-1

# Model and its beams
cHALE,leftWing,rightWing,tailboom,_ = create_cHALE(aeroSolver=aeroSolver,altitude=h,stiffnessFactor=λ,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,withAttachmentSprings=true,μ=μ)

# Root element's global ID
rootElemGlobalID = rightWing.elements[1].assemblyData.globalID

# System solver 
NR = create_NewtonRaphson(trimSolverMethod=:LM,displayStatus=false)

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
    update_v_A!(cHALE,v_A=[0,U,0])
    # Set initial guess solution for the trim problem
    x0 = i==1 ? zeros(0) : problem[i-1].systemArrays.x
    # Create and solve trim problem
    problem[i] = create_TrimProblem(model=cHALE,systemSolver=NR,x0=x0)
    solve!(problem[i])
    # Trim results
    trimSol = problem[i].trimSolution
    aeroVars = problem[i].loadstepHistory.aeroVariables
    trimAoA[i] = rad2deg(aeroVars[end][rootElemGlobalID].unsteadyKin.αₑ)
    trimThrust[i] = trimSol.trimLoads[1]
    trimδ[i] = stabilizersAero ? 
               rad2deg(trimSol.trimδ[1]) : 
               NaN
    println("AoA = $(round(trimAoA[i], digits=2)) deg, T = $(round(trimThrust[i], digits=2)) N, δ = $(round(trimδ[i], digits=2)) deg")
end

## Plots
# Set paths
relPath = "/test/outputs/figures/cHALEfullTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformation plot
deformationPlot = plot_steady_deformation(problem[end],save=true,savePath=string(relPath,"/cHALEtrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
lw = 2
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[20,35], ylims=[0,15])
plot!(URange, trimAoA, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/cHALEfullTrim_AoA.pdf"))

# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[20,35])
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/cHALEfullTrim_thrust.pdf"))

# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[20,35])
plot!(URange, trimδ, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/cHALEfullTrim_delta.pdf"))

println("Finished cHALEfullTrim.jl")