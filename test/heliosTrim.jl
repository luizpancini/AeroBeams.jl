using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Wing airfoil
# wingAirfoil = NACA23012A
wingAirfoil = HeliosWingAirfoil

# Option for reduced chord
reducedChord = false

# TF to include beam pods 
beamPods = true

# Option to set payload on wing
payloadOnWing = false

# Aerodynamic solver
aeroSolver = Indicial()

# Set NR system solver 
relaxFactor = 0.5
displayStatus = false
NR = create_NewtonRaphson(ρ=relaxFactor,displayStatus=displayStatus)

# Airspeed
U = 40*0.3048

# Set stiffness factor and payload ranges, and initialize outputs
λRange = [1,50]
PRange = collect(0:20:500)
trimAoA = Array{Float64}(undef,length(λRange),length(PRange))
trimThrust = Array{Float64}(undef,length(λRange),length(PRange))
trimδ = Array{Float64}(undef,length(λRange),length(PRange))

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Sweep payload
    for (j,P) in enumerate(PRange)
        # Display progress
        println("Trimming for λ = $λ payload = $P lb")
        # Model and its beams
        helios,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,reducedChord=reducedChord,wingAirfoil=wingAirfoil,payloadOnWing=payloadOnWing,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
        # plt = plot_undeformed_assembly(helios)
        # display(plt)
        # Set initial guess solution as previous known solution
        x0 = (j==1) ? zeros(0) : problem.x
        # Create and solve trim problem
        global problem = create_TrimProblem(model=helios,systemSolver=NR,x0=x0)
        solve!(problem)
        # Trim results
        trimAoA[i,j] = problem.flowVariablesOverσ[end][midSpanElem].αₑ*180/π
        trimThrust[i,j] = problem.x[end-1]*problem.model.forceScaling
        trimδ[i,j] = problem.x[end]*180/π
        println("AoA = $(trimAoA[i,j]), T = $(trimThrust[i,j]), δ = $(trimδ[i,j])")
    end
end

# Load reference data
αFlexibleRef = readdlm(string(pwd(),"/test/referenceData/Helios/trim_AoA_flexible.txt"))
αRigidRef = readdlm(string(pwd(),"/test/referenceData/Helios/trim_AoA_rigid.txt"))
δFlexibleRef = readdlm(string(pwd(),"/test/referenceData/Helios/trim_delta_flexible.txt"))
δRigidRef = readdlm(string(pwd(),"/test/referenceData/Helios/trim_delta_rigid.txt"))

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(λRange)))
lw = 2
ms = 3
labels = ["Flexible" "Rigid"]
# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Payload [lb]", ylabel="Trim root AoA [deg]", xlims=[0,500], ylims=[0,5])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, label="Patil & Hodges (2006)")
for (i,λ) in enumerate(λRange)
    plot!(PRange, trimAoA[i,:], c=colors[i], lw=lw, label=labels[i])
    if i==1
        scatter!(αFlexibleRef[1,:], αFlexibleRef[2,:], c=colors[i], ms=ms, msw=0, label=false)
    else
        scatter!(αRigidRef[1,:], αRigidRef[2,:], c=colors[i], ms=ms, msw=0, label=false)
    end
end
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/heliosTrim_AoA.pdf"))
# Trim propeller force vs airspeed
plt2 = plot(xlabel="Payload [lb]", ylabel="Trim thrust per motor [N]", xlims=[0,500])
for (i,λ) in enumerate(λRange)
    plot!(PRange, trimThrust[i,:], c=colors[i], lw=lw, label=labels[i])
end
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/heliosTrim_thrust.pdf"))
# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Payload [lb]", ylabel="Trim elevator deflection [deg]", xlims=[0,500], ylims=[0,10])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, label="Patil & Hodges (2006)")
for (i,λ) in enumerate(λRange)
    plot!(PRange, trimδ[i,:], c=colors[i], lw=lw, label=labels[i])
    if i==1
        scatter!(δFlexibleRef[1,:], δFlexibleRef[2,:], c=colors[i], ms=ms, msw=0, label=false)
    else
        scatter!(δRigidRef[1,:], δRigidRef[2,:], c=colors[i], ms=ms, msw=0, label=false)
    end
end
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/heliosTrim_delta.pdf"))

println("Finished heliosTrim.jl")