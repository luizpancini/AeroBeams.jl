using AeroBeams, DelimitedFiles

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = false
includeVS = false

# Parasite drag coefficients
wingCd0 = stabsCd0 = 0

# Option to include induced drag
hasInducedDrag = false

# Discretization
nElemWing = 40
nElemTailBoom = 10
nElemHorzStabilizer = 10

# Wing precurvature
k2 = 0.0

# NR system solver
relaxFactor = 0.5
maxIter = 100
displayStatus = false
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,displayStatus=displayStatus)

# Set stiffness factor and airspeed ranges
λRange = [1; 50]
URange = collect(20:1:35)

# Initialize outputs
trimProblem = Array{TrimProblem}(undef,length(λRange),length(URange))
trimAoA = Array{Float64}(undef,length(λRange),length(URange))
trim_u1 = Array{Vector{Float64}}(undef,length(λRange),length(URange))
trim_u3 = Array{Vector{Float64}}(undef,length(λRange),length(URange))

# Sweep stiffness factor
for (i,λ) in enumerate(λRange)
    # Model
    conventionalHALE,leftWing,rightWing,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,hasInducedDrag=hasInducedDrag,k2=k2)
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
        x0 = j == 1 ? zeros(0) : trimProblem[i,j-1].x
        # Create and solve trim problem
        trimProblem[i,j] = create_TrimProblem(model=conventionalHALE,systemSolver=NR,x0=x0)
        solve!(trimProblem[i,j])
        # Trim results
        trimAoA[i,j] = trimProblem[i,j].aeroVariablesOverσ[end][rWGlobalElemRange[1]].flowAnglesAndRates.αₑ*180/π
        trim_u1[i,j] = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[1],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in rWGlobalElemRange]...)
        trim_u3[i,j] = vcat([vcat(trimProblem[i,j].nodalStatesOverσ[end][e].u_n1[3],trimProblem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in rWGlobalElemRange]...)
        println("AoA = $(trimAoA[i,j]) deg")
    end
end

# Load reference solutions
trimAoAERef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/conventionalHALE/trimAoAVsAirspeedElastic.txt")
trimAoARRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/conventionalHALE/trimAoAVsAirspeedRigid.txt")
trimDispRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/conventionalHALE/trimDispAtU25.txt")

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/outputs/figures/cHALE_trim_compare_Patil"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, 2, categorical=true)
labels = ["Elastic" "Rigid"]
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0
L = 16
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root angle of attack [deg]", xlims=[URange[1],URange[end]], ylims=[0,15], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
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
savefig(string(absPath,"/cHALE_trim_compare_Patil_AoA.pdf"))

# Trim deflected wingspan at U = 25 m/s
U2plot = 25.0
indU = findfirst(x->x==U2plot,URange)
if !isnothing(indU)
    plt2 = plot(xlabel="Spanwise length [m]", ylabel="Vertical displacement [m]", xlims=[0,16], ylims=[0,16], xticks=collect(0:4:16), yticks=collect(0:4:16), tickfont=font(ts), guidefont=font(fs))
    plot!(x1_0.+(trim_u1[1,indU].-trim_u1[1,indU][1]), x3_0.+trim_u3[1,indU].-trim_u3[1,indU][1], c=:black, lw=lw, label=false)
    scatter!(trimDispRef[1,:], trimDispRef[2,:], c=:black, ms=ms, msw=msw, label=false)
    display(plt2)
    savefig(string(absPath,"/cHALE_trim_compare_Patil_disp.pdf"))
end

# Trim deflected wingspan over airspeed
plt3 = plot(xlabel="Normalized spanwise length", ylabel="Vertical displacement [% semispan]", xlims=[0,1], ylims=[0,100], tickfont=font(ts), guidefont=font(fs))
for (i,U) in enumerate(URange)
    plot!((x1_0.+(trim_u1[1,i].-trim_u1[1,i][1]))/L, (x3_0.+(trim_u3[1,i].-trim_u3[1,i][1]))/L*100, lz=U, c=:rainbow, lw=lw, label=false, colorbar_title="Airspeed [m/s]")
end
display(plt3)
savefig(string(absPath,"/cHALE_trim_compare_Patil_u3OverU.pdf"))

println("Finished cHALE_trim_compare_Patil.jl")