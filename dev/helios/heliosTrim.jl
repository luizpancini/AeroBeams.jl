using AeroBeams, DelimitedFiles

# Wing airfoil
wingAirfoil = deepcopy(NACA23012A)

# Option for reduced chord
reducedChord = true

# Flag to include beam pods 
beamPods = true

# Flag to set payload on wing
payloadOnWing = false

# Circulatory indicial function
circulatoryIndicialFunction = "Wagner"

# Aerodynamic solvers
aeroSolvers = [Indicial(circulatoryIndicialFunction=circulatoryIndicialFunction); BLi(circulatoryIndicialFunction=circulatoryIndicialFunction)]

# Set NR system solver 
relaxFactor = 0.5
maxIter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,displayStatus=false)

# Airspeed [m/s]
U = 40*0.3048

# Set payload range
PRange = collect(0:20:500)

# Initialize outputs
problem = Array{TrimProblem}(undef,length(aeroSolvers),length(PRange))
trimAoA = Array{Float64}(undef,length(aeroSolvers),length(PRange))
trimThrust = Array{Float64}(undef,length(aeroSolvers),length(PRange))
trimδ = Array{Float64}(undef,length(aeroSolvers),length(PRange))

# Sweep aero solvers
for (i,aeroSolver) in enumerate(aeroSolvers)
    # Sweep payload
    for (j,P) in enumerate(PRange)
        # Display progress
        println("Trimming for aeroSolver $i, payload = $P lb")
        # Model
        model,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,wingAirfoil=wingAirfoil,payloadOnWing=payloadOnWing,reducedChord=reducedChord,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
        # Set initial guess solution as previous known solution
        x0 = (j==1) ? zeros(0) : problem[i,j-1].x
        # Create and solve trim problem
        problem[i,j] = create_TrimProblem(model=model,systemSolver=NR,x0=x0)
        solve!(problem[i,j])
        # Trim results
        trimAoA[i,j] = problem[i,j].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.αₑ*180/π
        trimThrust[i,j] = problem[i,j].x[end-1]*problem[i,j].model.forceScaling
        trimδ[i,j] = problem[i,j].x[end]*180/π
        println("AoA = $(trimAoA[i,j]), T = $(trimThrust[i,j]), δ = $(trimδ[i,j])")
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
labels = ["Linear aero" "Dynamic stall"]
colors = cgrad(:rainbow, length(aeroSolvers), categorical=true)

# Trim root angle of attack
plt_AoA = plot(xlabel="Payload [lb]", ylabel="Trim root AoA [deg]", xlims=[0,500], ylims=[0,7], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(PRange, trimAoA[i,:], c=colors[i], lw=lw, label=labels[i])
end
display(plt_AoA)
savefig(string(absPath,"/heliosTrim_AoA.pdf"))

# Trim elevator deflection
plt_delta = plot(xlabel="Payload [lb]", ylabel="Trim elevator deflection [deg]", xlims=[0,500], ylims=[0,12], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(PRange, trimδ[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_delta)
savefig(string(absPath,"/heliosTrim_delta.pdf")) 

# Trim thrust per motor
plt_thrust = plot(xlabel="Payload [lb]", ylabel="Trim thrust per motor [N]", xlims=[0,500], ylims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,aeroSolver) in enumerate(aeroSolvers)
    plot!(PRange, trimThrust[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_thrust)
savefig(string(absPath,"/heliosTrim_thrust.pdf")) 

println("Finished heliosTrim.jl")