using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Wing precurvature
k2 = 0.0

# NR system solver 
maxit = 100
displayStatus = true
NR = create_NewtonRaphson(maximumIterations=maxit,displayStatus=displayStatus)

# Set stiffness factor and airspeed ranges, and initialize outputs
λRange = [1; 50]
URange = collect(20:1:35)
trimAoA = Array{Float64}(undef,length(λRange),length(URange))
trim_u1 = Array{Vector{Float64}}(undef,length(λRange),length(URange))
trim_u3 = Array{Vector{Float64}}(undef,length(λRange),length(URange))

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Model
    conventionalHALE,leftWing,rightWing,_ = create_conventional_HALE(aeroSolver=Inflow(),nElemWing=20,stiffnessFactor=λ,stabilizersAero=false,includeVS=false,k2=k2)
    # plt = plot_undeformed_assembly(conventionalHALE,(0,0))
    # display(plt)
    # Get element ranges and nodal arclength positions of right wing
    global x1_0 = vcat([vcat(rightWing.elements[e].r_n1[1],rightWing.elements[e].r_n2[1]) for e in 1:rightWing.nElements]...)
    global x3_0 = vcat([vcat(rightWing.elements[e].r_n1[3],rightWing.elements[e].r_n2[3]) for e in 1:rightWing.nElements]...)
    rWGlobalElemRange = rightWing.elementRange[1]:rightWing.elementRange[end]
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Trimming for λ = $λ, U = $U m/s")
        # Update airspeed on model
        set_motion_basis_A!(model=conventionalHALE,v_A=[0;U;0])
        # Set initial guess solution as previous known solution
        x0 = j == 1 ? zeros(0) : problem.x
        # Create and solve trim problem
        global problem = create_TrimProblem(model=conventionalHALE,systemSolver=NR,x0=x0)
        solve!(problem)
        # Trim results
        trimAoA[i,j] = problem.flowVariablesOverσ[end][rWGlobalElemRange[1]].αₑ*180/π
        trim_u1[i,j] = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[1],problem.nodalStatesOverσ[end][e].u_n2[1]) for e in rWGlobalElemRange]...)
        trim_u3[i,j] = vcat([vcat(problem.nodalStatesOverσ[end][e].u_n1[3],problem.nodalStatesOverσ[end][e].u_n2[3]) for e in rWGlobalElemRange]...)
        println("AoA = $(trimAoA[i,j]) deg")
    end
end

# Load reference solutions
trimAoAERef = readdlm(string(pwd(),"/test/referenceData/conventionalHALE/trimAoAVsAirspeedElastic.txt"))
trimAoARRef = readdlm(string(pwd(),"/test/referenceData/conventionalHALE/trimAoAVsAirspeedRigid.txt"))
trimDispRef = readdlm(string(pwd(),"/test/referenceData/conventionalHALE/trimDispAtU25.txt"))

# Plots
# ------------------------------------------------------------------------------
colorsU = get(colorschemes[:rainbow], LinRange(0, 1, length(URange)))
colors = get(colorschemes[:rainbow], LinRange(0, 1, 2))
labels = ["Elastic" "Rigid"]
lw = 2
ms = 3
# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root angle of attack [deg]", xlims=[URange[1],URange[end]], ylims=[0,20])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, msw=0, label="Patil et al. (2001)")
for (i,λ) in enumerate(λRange)
    plot!(URange, trimAoA[i,:], c=colors[i], lw=lw, label=labels[i])
    if i==1
        scatter!(trimAoAERef[1,:], trimAoAERef[2,:], c=colors[i], ms=ms, msw=0, label=false)
    else
        scatter!(trimAoARRef[1,:], trimAoARRef[2,:], c=colors[i], ms=ms, msw=0, label=false)
    end
end
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEtrim_AoA.pdf"))
display(plt1)
# Trim deflected wingspan at U = 25 m/s
U2plot = 25.0
indU = findfirst(x->x==U2plot,URange)
if !isnothing(indU)
    plt2 = plot(xlabel="Spanwise length [m]", ylabel="Vertical displacement [m]", xlims=[0,16], ylims=[0,16], xticks=collect(0:4:16), yticks=collect(0:4:16))
    plot!(x1_0.+(trim_u1[1,indU].-trim_u1[1,indU][1]), x3_0.+trim_u3[1,indU].-trim_u3[1,indU][1], c=:black, lw=lw, label="AeroBeams")
    scatter!(trimDispRef[1,:], trimDispRef[2,:], c=:black, ms=ms, msw=0, label="Patil et al. (2001)")
    display(plt2)
    savefig(string(pwd(),"/test/outputs/figures/conventionalHALEtrim_disp.pdf"))
end
# Trim deflected wingspan over airspeed
plt3 = plot(xlabel="Normalized spanwise length", ylabel="Vertical displacement [% semispan]", xlims=[0,1], ylims=[0,100])
for (i,U) in enumerate(URange)
    plot!((x1_0.+(trim_u1[1,i].-trim_u1[1,i][1]))/16, (x3_0.+(trim_u3[1,i].-trim_u3[1,i][1]))/16*100, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/conventionalHALEtrim_u3OverU.pdf"))

println("Finished conventionalHALEtrim.jl")