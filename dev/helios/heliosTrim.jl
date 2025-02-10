using AeroBeams, DelimitedFiles

# Wing airfoil
wingAirfoil = deepcopy(NACA23012A)

# Option for reduced chord
reducedChord = false

# Flag to include beam pods 
beamPods = true

# Flag to set payload on wing
payloadOnWing = false

# Aerodynamic solver
aeroSolver = Indicial()

# Set NR system solver 
relaxFactor = 0.5
NR = create_NewtonRaphson(ρ=relaxFactor,displayStatus=false)

# Airspeed [m/s]
U = 40*0.3048

# Set payload range
PRange = collect(0:20:500)

# Initialize outputs
problem = Array{TrimProblem}(undef,length(PRange))
trimAoA = Array{Float64}(undef,length(PRange))
trimThrust = Array{Float64}(undef,length(PRange))
trimδ = Array{Float64}(undef,length(PRange))
 
# Sweep payload
for (i,P) in enumerate(PRange)
    # Display progress
    println("Trimming for payload = $P lb")
    # Model
    helios,midSpanElem,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,wingAirfoil=wingAirfoil,payloadOnWing=payloadOnWing,reducedChord=reducedChord,payloadPounds=P,airspeed=U,δIsTrimVariable=true,thrustIsTrimVariable=true)
    # Set initial guess solution as previous known solution
    x0 = (i==1) ? zeros(0) : problem[i-1].x
    # Create and solve trim problem
    problem[i] = create_TrimProblem(model=helios,systemSolver=NR,x0=x0)
    solve!(problem[i])
    # Trim results
    trimAoA[i] = problem[i].aeroVariablesOverσ[end][midSpanElem].flowAnglesAndRates.αₑ*180/π
    trimThrust[i] = problem[i].x[end-1]*problem[i].model.forceScaling
    trimδ[i] = problem[i].x[end]*180/π
    println("AoA = $(trimAoA[i]), T = $(trimThrust[i]), δ = $(trimδ[i])")
end
 
# Load reference data
αFlexibleRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Helios/trim_AoA_flexible.txt")
αRigidRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Helios/trim_AoA_rigid.txt")
δFlexibleRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Helios/trim_delta_flexible.txt")
δRigidRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Helios/trim_delta_rigid.txt")
TRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Helios/trim_thrust.txt") 

# Set paths
relPath = "/dev/helios/figures/heliosTrim"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 16
lw = 2
ms = 4
msw = 0

# Trim root angle of attack
plt1 = plot(xlabel="Payload [lb]", ylabel="Trim root AoA [deg]", xlims=[0,500], ylims=[0,6], tickfont=font(ts), guidefont=font(fs), legendfontsize=12)
plot!(PRange, trimAoA, c=:black, lw=lw, label="AeroBeams")
scatter!(αFlexibleRef[1,:], αFlexibleRef[2,:], c=:black, ms=ms, msw=msw, label="Patil & Hodges (2006)")
display(plt1)
savefig(string(absPath,"/heliosTrim_AoA.pdf"))

# Trim elevator deflection
plt2 = plot(xlabel="Payload [lb]", ylabel="Trim elevator deflection [deg]", xlims=[0,500], ylims=[0,6], tickfont=font(ts), guidefont=font(fs))
plot!(PRange, trimδ, c=:black, lw=lw, label=false)
scatter!(δFlexibleRef[1,:], δFlexibleRef[2,:], c=:black, ms=ms, msw=msw, label=false)
display(plt2)
savefig(string(absPath,"/heliosTrim_delta.pdf")) 

# Trim thrust per motor
plt3 = plot(xlabel="Payload [lb]", ylabel="Trim thrust per motor [N]", xlims=[0,500], ylims=[0,40], tickfont=font(ts), guidefont=font(fs))
plot!(PRange, trimThrust, c=:black, lw=lw, label=false)
scatter!(TRef[1,:], TRef[2,:], c=:black, ms=ms, msw=msw, label=false)
display(plt3)
savefig(string(absPath,"/heliosTrim_thrust.pdf")) 

println("Finished heliosTrim.jl")