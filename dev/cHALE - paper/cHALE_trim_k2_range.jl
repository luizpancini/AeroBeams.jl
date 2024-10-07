using AeroBeams

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Stiffness factor for the aircraft's wing
λ = 1e0

# Discretization
nElemWing = 20
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# System solver for trim problem
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
σstep = 0.5
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,displayStatus=false)

# Set bending curvature and airspeed ranges
k2Range = range(-0.015, 0.045, 5)
URange = collect(20:1:35)

# Initialize outputs
trimAoA = Array{Float64}(undef,length(k2Range),length(URange))
trimThrust = Array{Float64}(undef,length(k2Range),length(URange))
trimδ = Array{Float64}(undef,length(k2Range),length(URange))

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model for trim problem
        cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stiffnessFactor=λ,∞=1e12,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag)
        # Update model
        cHALEtrim.skipValidationMotionBasisA = true
        update_model!(cHALEtrim)
        # Set initial guess solution as previous known solution
        x0Trim = j == 1 ? zeros(0) : trimProblem.x
        # Create and trim problem
        global trimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NR,x0=x0Trim)
        solve!(trimProblem)
        # Extract trim variables
        trimAoA[i,j] = trimProblem.aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
        trimThrust[i,j] = stabilizersAero ? trimProblem.x[end-1]*trimProblem.model.forceScaling : trimProblem.x[end]*trimProblem.model.forceScaling
        trimδ[i,j] = stabilizersAero ? trimProblem.x[end] : 0
        println("Trim AoA = $(trimAoA[i,j]*180/π), trim thrust = $(trimThrust[i,j]), trim δ = $(trimδ[i,j]*180/π)")
    end
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/outputs/figures/cHALE_trim_k2_range"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], range(0, 1, length(k2Range)))
lw = 2
ms = 3
msw = 0
mshape = [:circle, :star, :utriangle, :pentagon]
labels = ["\$k_2 = $(k2) \$" for k2 in k2Range]
L = 16
gr()

# Trim root angle of attack
pltT1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim \$\\alpha_r\$ [deg]", xlims=[URange[1],URange[end]], ylims=[0,15])
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimAoA[i,:]*180/π, c=colors[i], lw=lw, label=labels[i])
end
display(pltT1)
savefig(string(absPath,"/cHALE_flutter_k2_range_trimAoA.pdf"))

# Trim thrust
pltT2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[URange[1],URange[end]], ylims=[0,40])
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimThrust[i,:], c=colors[i], lw=lw, label=labels[i])
end
display(pltT2)
savefig(string(absPath,"/cHALE_flutter_k2_range_trimThrust.pdf"))

# Trim elevator deflection
pltT3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim \$\\delta_f\$ [deg]", xlims=[URange[1],URange[end]], ylims=[-40,5])
for (i,k2) in enumerate(k2Range)
    plot!(URange, trimδ[i,:]*180/π, c=colors[i], lw=lw, label=labels[i])
end
display(pltT3)
savefig(string(absPath,"/cHALE_flutter_k2_range_trimDelta.pdf"))

println("Finished cHALE_flutter_k2_range.jl")