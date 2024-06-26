using AeroBeams, LinearAlgebra, Plots

# Stiffness factor
λ = 1e0

# Aerodynamic solver
aeroSolver = Indicial()

# Option for stabilizers
stabilizersAero = true
includeVS = true

# Parasite drag coefficients
wingCd0 = 1e-2
stabsCd0 = 1e-2

# Model and its beams
conventionalHALE,leftWing,rightWing,tailboom,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true)

# Set NR system solver 
relaxFactor = 0.5
displayStatus = false
NR = create_NewtonRaphson(ρ=relaxFactor,displayStatus=displayStatus)

# Set airspeed range and initialize outputs
URange = collect(20:1:35)
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))

# Add attachment springs
μ = 1e-1
ku = μ*[1; 1; 1]
kp = ku
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)
spring2 = create_Spring(elementsIDs=[nElemTailBoom],nodesSides=[2],ku=ku,kp=kp)
add_springs_to_beam!(beam=tailboom,springs=[spring1,spring2])
update_model!(conventionalHALE)

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Trimming for airspeed U = $U m/s")
    # Update airspeed on model
    set_motion_basis_A!(model=conventionalHALE,v_A=[0;U;0])
    # Set initial guess solution as previous known solution
    x0 = i == 1 ? zeros(0) : problem.x
    # Create and solve trim problem
    global problem = create_TrimProblem(model=conventionalHALE,systemSolver=NR,x0=x0)
    solve!(problem)
    # Trim results
    trimAoA[i] = problem.flowVariablesOverσ[end][rightWing.elementRange[1]].αₑ*180/π
    trimThrust[i] = stabilizersAero ? problem.x[end-1]*problem.model.forceScaling : problem.x[end]*problem.model.forceScaling
    trimδ[i] = stabilizersAero ? problem.x[end]*180/π : 0
    println("AoA = $(trimAoA[i]), T = $(trimThrust[i]), δ = $(trimδ[i])")
end

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 3
# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[20,35], ylims=[0,15])
plot!(URange, trimAoA, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEfullTrim_AoA.pdf"))
# Trim propeller force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[20,35])
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEfullTrim_thrust.pdf"))
# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[20,35])
plot!(URange, trimδ, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEfullTrim_delta.pdf"))

println("Finished conventionalHALEfullTrim.jl")