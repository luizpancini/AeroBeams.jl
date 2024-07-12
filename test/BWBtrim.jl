using AeroBeams, LinearAlgebra, Plots

# Stiffness factor
λ = 1e0

# Aerodynamic solver
aeroSolver = Indicial()

# Tip loss option
hasTipCorrection = true

# Model 
BWB = create_BWB(aeroSolver=aeroSolver,stiffnessFactor=λ,δElevIsTrimVariable=true,thrustIsTrimVariable=true,hasTipCorrection=hasTipCorrection)
# plt = plot_undeformed_assembly(BWB)
# display(plt)

# Set NR system solver 
relaxFactor = 0.5
displayStatus = true
maxiter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxiter,displayStatus=displayStatus)

# Attachment springs
μ = 1e-2
ku = μ*[1; 1; 1]
kp = ku
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)
spring2 = create_Spring(elementsIDs=[3],nodesSides=[2],ku=ku,kp=kp)
add_springs_to_beam!(beam=BWB.beams[2],springs=[spring1])
add_springs_to_beam!(beam=BWB.beams[3],springs=[spring2])

# Set airspeed range and initialize outputs
URange = collect(40:5:120)
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Trimming for airspeed U = $U m/s")
    # Update airspeed on model
    set_motion_basis_A!(model=BWB,v_A=[0;U;0])
    # Set initial guess solution as previous known solution
    x0 = i == 1 ? zeros(0) : problem.x
    # Create and solve trim problem
    global problem = create_TrimProblem(model=BWB,systemSolver=NR,x0=x0)
    solve!(problem)
    # Trim results
    trimAoA[i] = problem.flowVariablesOverσ[end][BWB.beams[3].elementRange[1]].αₑ*180/π
    trimThrust[i] = problem.x[end-1]*problem.model.forceScaling 
    trimδ[i] = problem.x[end]*180/π 
    println("AoA = $(trimAoA[i]), T = $(trimThrust[i]), δ = $(trimδ[i])")
end

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 3
# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[40,120])
plot!(URange, trimAoA, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/BWBtrim_AoA.pdf"))
# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[40,120])
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/BWBtrim_thrust.pdf"))
# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[40,120])
plot!(URange, trimδ, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/BWBtrim_delta.pdf"))

println("Finished BWBtrim.jl")