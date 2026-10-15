using AeroBeams

# Stiffness factor range
λRange = [1; 50]

# Airspeed range [m/s] 
URange = collect(20:1:35)

# Altitude [m]
h = 20e3

# Aerodynamic solver
aeroSolver = Indicial()

# Option for stabilizers
stabilizersAero = false
includeVS = false

# Parasite drag coefficients
wingCd0 = stabsCd0 = 0

# Discretization
nElemWing = 20
nElemTailBoom = 5
nElemHorzStabilizer = 4

# Attachment spring stiffness
μ = 1e-1

# System solver 
NR = create_NewtonRaphson(trimSolverMethod=:LQ,displayStatus=false)

# Initialize model to get nodal arclength positions and element ranges of right wing
_,_,rightWing,_ = create_cHALE(nElemWing=nElemWing)
x1_0 = vcat([vcat(rightWing.elements[e].assemblyData.r_n1[1],rightWing.elements[e].assemblyData.r_n2[1]) for e in 1:rightWing.discretization.nElements]...)
x3_0 = vcat([vcat(rightWing.elements[e].assemblyData.r_n1[3],rightWing.elements[e].assemblyData.r_n2[3]) for e in 1:rightWing.discretization.nElements]...)
rWGlobalElemRange = rightWing.assemblyData.elementsGlobalIDs

# Initialize outputs
problem = Array{TrimProblem}(undef,length(λRange),length(URange))
trimAoA = Array{Float64}(undef,length(λRange),length(URange))
trim_u1 = Array{Vector{Float64}}(undef,length(λRange),length(URange))
trim_u3 = Array{Vector{Float64}}(undef,length(λRange),length(URange))

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Model
    cHALE,_ = create_cHALE(aeroSolver=aeroSolver,stiffnessFactor=λ,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,withAttachmentSprings=true,μ=μ)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Trimming for λ = $λ, U = $U m/s")
        # Update airspeed on model
        update_v_A!(cHALE,v_A=[0,U,0])
        # Set initial guess solution for the trim problem
        x0 = j==1 ? zeros(0) : problem[i,j-1].systemArrays.x
        # Create and solve trim problem
        problem[i,j] = create_TrimProblem(model=cHALE,systemSolver=NR,x0=x0)
        solve!(problem[i,j])
        # Trim results
        states = problem[i,j].loadstepHistory.states
        aeroVars = problem[i,j].loadstepHistory.aeroVariables
        trimAoA[i,j] = rad2deg(aeroVars[end][rWGlobalElemRange[1]].unsteadyKin.αₑ)
        trim_u1[i,j] = vcat([vcat(states[end][e].nodal.u1[1],states[end][e].nodal.u2[1]) for e in rWGlobalElemRange]...)
        trim_u3[i,j] = vcat([vcat(states[end][e].nodal.u1[3],states[end][e].nodal.u2[3]) for e in rWGlobalElemRange]...)
        println("AoA = $(round(trimAoA[i,j], digits=2)) deg")
    end
end

# Load reference solutions
trimAoAERef = readdlm(pkgdir(AeroBeams,"test","referenceData","cHALE", "trimAoAVsAirspeedElastic.txt"))
trimAoARRef = readdlm(pkgdir(AeroBeams,"test","referenceData","cHALE", "trimAoAVsAirspeedRigid.txt"))
trimDispRef = readdlm(pkgdir(AeroBeams,"test","referenceData","cHALE", "trimDispAtU25.txt"))

## Plots
# Set paths
relPath = "/test/outputs/figures/cHALEtrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformation plot
deformationPlot = plot_steady_deformation(problem[end],element2centralize=div(nElemWing,2)+1,save=true,savePath=string(relPath,"/cHALEtrim_deformation.pdf"))
display(deformationPlot)

# Plot configurations
colorsU = cgrad(:rainbow, length(URange), categorical=true)
colors = cgrad(:rainbow, 2, categorical=true)
labels = ["Elastic" "Rigid"]
lw = 2
ms = 3
msw = 0
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root angle of attack [deg]", xlims=[URange[1],URange[end]], ylims=[0,20])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Patil et al. (2001)")
for (i,λ) in enumerate(λRange)
    plot!(URange, trimAoA[i,:], c=colors[i], lw=lw, label=labels[i])
    if i==1
        scatter!(trimAoAERef[1,:], trimAoAERef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
    else
        scatter!(trimAoARRef[1,:], trimAoARRef[2,:], c=colors[i], ms=ms, msw=msw, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/cHALEtrim_AoA.pdf"))

# Trim deflected wingspan at U = 25 m/s
U2plot = 25.0
indU = findfirst(x->x==U2plot,URange)
if !isnothing(indU)
    plt2 = plot(xlabel="Spanwise length [m]", ylabel="Vertical displacement [m]", xlims=[0,16], ylims=[0,16], xticks=collect(0:4:16), yticks=collect(0:4:16))
    plot!(x1_0.+(trim_u1[1,indU].-trim_u1[1,indU][1]), x3_0.+trim_u3[1,indU].-trim_u3[1,indU][1], c=:black, lw=lw, label="AeroBeams")
    scatter!(trimDispRef[1,:], trimDispRef[2,:], c=:black, ms=ms, msw=msw, label="Patil et al. (2001)")
    display(plt2)
    savefig(string(absPath,"/cHALEtrim_disp.pdf"))
end

# Trim deflected wingspan over airspeed
plt3 = plot(xlabel="Normalized spanwise length", ylabel="Vertical displacement [% semispan]", xlims=[0,1], ylims=[0,70])
for (i,U) in enumerate(URange)
    plot!((x1_0.+(trim_u1[1,i].-trim_u1[1,i][1]))/16, (x3_0.+(trim_u3[1,i].-trim_u3[1,i][1]))/16*100, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt3)
savefig(string(absPath,"/cHALEtrim_u3OverU.pdf"))

println("Finished cHALEtrim.jl")