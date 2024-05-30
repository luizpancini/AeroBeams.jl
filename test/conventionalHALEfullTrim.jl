using AeroBeams, LinearAlgebra, Plots

# Stiffness factor
λ = 1

# Model and its beams
conventionalHALE,leftWing,rightWing,_ = create_conventional_HALE(aeroSolver=Inflow(),nElemWing=20,stiffnessFactor=λ,stabilizersAero=true,includeVS=true,wingCd0=1e-2,stabsCd0=1e-2,δElevIsTrimVariable=true,thrustIsTrimVariable=true)
# plt = plot_undeformed_assembly(conventionalHALE)
# display(plt)

# Set NR system solver 
relaxFactor = 0.5
displayStatus = true
NR = create_NewtonRaphson(ρ=relaxFactor,displayStatus=displayStatus)

# Set airspeed range and initialize outputs
URange = collect(20:0.5:35)
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))

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
    trimThrust[i] = problem.x[end-1]*problem.model.forceScaling
    trimδ[i] = problem.x[end]*180/π
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