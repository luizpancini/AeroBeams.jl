using AeroBeams, DelimitedFiles

# Wing airfoil
wingAirfoil = deepcopy(NACA23012A)

# Option for reduced chord
reducedChord = false

# Flag to include beam pods 
beamPods = true

# Flag to set payload on wing
payloadOnWing = false

# Circulatory indicial function
circulatoryIndicialFunction = "Wagner"

# Airspeed [m/s]
U = 40*0.3048

# Aerodynamic solvers
aeroSolvers = [Indicial(circulatoryIndicialFunction=circulatoryIndicialFunction); BLi(circulatoryIndicialFunction=circulatoryIndicialFunction)]

# Stiffness range
λRange = [1, 1e4]

# Payload range
PRange = collect(0:10:500)

# Set NR system solver
relTol = 1e-11
absTol = 1e-11
relaxFactor = 0.5
maxIter = 100
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,displayStatus=false,relativeTolerance=relTol,absoluteTolerance=absTol)

# Initialize outputs
problem = Array{TrimProblem}(undef,length(aeroSolvers),length(λRange),length(PRange))
trimAoA = Array{Float64}(undef,length(aeroSolvers),length(λRange),length(PRange))
trimThrust = Array{Float64}(undef,length(aeroSolvers),length(λRange),length(PRange))
trimδ = Array{Float64}(undef,length(aeroSolvers),length(λRange),length(PRange))
tip_OOP = Array{Float64}(undef,length(aeroSolvers),length(λRange),length(PRange))

# Sweep aero solvers
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Sweep stiffness factor
    for (j,λ) in enumerate(λRange)
        # Sweep payload
        for (k,P) in enumerate(PRange)
            # Display progress
            println("Trimming for aeroSolver $i, λ = $λ, payload = $P lb")
            # Model
            model,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,wingAirfoil=wingAirfoil,payloadOnWing=payloadOnWing,reducedChord=reducedChord,stiffnessFactor=λ,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
            # Set initial guess solution as previous known solution
            x0 = (k==1) ? zeros(0) : problem[i,j,k-1].x
            # Create and solve trim problem
            problem[i,j,k] = create_TrimProblem(model=model,systemSolver=NR,x0=x0)
            solve!(problem[i,j,k])
            # Trim results
            trimAoA[i,j,k] = (problem[i,j,k].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.α-wingAirfoil.attachedFlowParameters.α₀N)*180/π
            trimThrust[i,j,k] = problem[i,j,k].x[end-1]*problem[i,j,k].model.forceScaling
            trimδ[i,j,k] = problem[i,j,k].x[end]*180/π
            println("AoA = $(trimAoA[i,j,k]), T = $(trimThrust[i,j,k]), δ = $(trimδ[i,j,k])")
            # Other outputs
            tip_OOP[i,j,k] = problem[i,j,k].nodalStatesOverσ[end][1].u_n1[3] - problem[i,j,k].nodalStatesOverσ[end][midSpanElem].u_n2[3]
        end
    end
end

# Set paths
relPath = "/dev/helios/figures/heliosTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 4
msw = 0
labels = ["Attached flow - Elastic" "Attached flow - Rigid"; "Dynamic stall - Elastic" "Dynamic stall - Rigid"]
colors = cgrad(:rainbow, length(aeroSolvers), categorical=true)
ls = [:solid :dash]

# Trim root angle of attack
plt_AoA = plot(xlabel="Payload [lb]", ylabel="Trim root AoA [deg]", xlims=[0,502], ylims=[0,6], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:bottomright)
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(PRange, trimAoA[i,j,:], c=colors[i], lw=lw, ls=ls[j], label=labels[i,j])
    end
end
display(plt_AoA)
savefig(string(absPath,"/heliosTrim_AoA.pdf"))

# Trim elevator deflection
plt_delta = plot(xlabel="Payload [lb]", ylabel="Trim elevator deflection [deg]", xlims=[0,502], ylims=[0,12], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(PRange, trimδ[i,j,:], c=colors[i], lw=lw, ls=ls[j], label=false)
    end
end
display(plt_delta)
savefig(string(absPath,"/heliosTrim_delta.pdf")) 

# Trim thrust per motor
plt_thrust = plot(xlabel="Payload [lb]", ylabel="Trim thrust per motor [N]", xlims=[0,502], ylims=[0,20], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(PRange, trimThrust[i,j,:], c=colors[i], lw=lw, ls=ls[j], label=false)
    end
end
display(plt_thrust)
savefig(string(absPath,"/heliosTrim_thrust.pdf")) 

# Trim tip OOP deflection
L = 0.3048*40*(2+cosd(10))
plt_OOP = plot(xlabel="Payload [lb]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,502], ylims=[0,30], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,λ) in enumerate(λRange)
        plot!(PRange, tip_OOP[i,j,:]/L*100, c=colors[i], lw=lw, ls=ls[j], label=false)
    end
end
display(plt_OOP)
savefig(string(absPath,"/heliosTrim_OOP.pdf")) 

println("Finished heliosTrim.kl")