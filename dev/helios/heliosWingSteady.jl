using AeroBeams, DelimitedFiles

# Airspeed range [m/s]
URange = 0.3048*vcat(0:2:30)

# Wing airfoil
wingAirfoil = deepcopy(NACA23012A)

# Option for reduced chord
reducedChord = false

# Flag to include beam pods 
beamPods = false

# Flag to set payload on wing
payloadOnWing = false

# Aerodynamic solver
aeroSolver = BLi()

# Root pitch angle
θ = 0*π/180

# Discretization
nElemStraightSemispan = 10
nElemDihedralSemispan = 5
nElem = nElemStraightSemispan + nElemDihedralSemispan

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(URange))
x1_def = Array{Vector{Float64}}(undef,length(URange))
x3_def = Array{Vector{Float64}}(undef,length(URange))
tipOOP = Array{Float64}(undef,length(URange))
tipTwist = Array{Float64}(undef,length(URange))
tipAoA = Array{Float64}(undef,length(URange))

# Initialize model
_,_,model,_ = create_Helios(aeroSolver=aeroSolver,beamPods=beamPods,wingAirfoil=wingAirfoil,payloadOnWing=payloadOnWing,reducedChord=reducedChord,airspeed=0,wingRootAoA=θ,nElemStraightSemispan=nElemStraightSemispan,nElemDihedralSemispan=nElemDihedralSemispan)

# Undeformed nodal positions
x1_0 = [model.r_n[n][1] for n in 1:nElem+1]
x3_0 = [model.r_n[n][3] for n in 1:nElem+1]

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $(U/.3048) ft/s")
    # Update airspeed
    set_motion_basis_A!(model=model,v_A=[0;U;0])
    # Set initial guess solution as previous known solution
    x0 = (i==1) ? zeros(0) : problem[i-1].x
    # Create and solve steady problem
    problem[i] = create_SteadyProblem(model=model,x0=x0)
    solve!(problem[i])
    # Outputs
    u1 = vcat([problem[i].nodalStatesOverσ[end][e].u_n1[1] for e in 1:nElem],problem[i].nodalStatesOverσ[end][end].u_n2[1])
    u3 = vcat([problem[i].nodalStatesOverσ[end][e].u_n1[3] for e in 1:nElem],problem[i].nodalStatesOverσ[end][end].u_n2[3])
    x1_def[i] = x1_0 .+ u1
    x3_def[i] = x3_0 .+ u3
    tip_p = model.elements[end].nodalStates.p_n2
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*AeroBeams.a2
    tipTwist[i] = asind(Δ[3])
    tipOOP[i] = model.elements[end].nodalStates.u_n2[3]
    tipAoA[i] = model.elements[end].aero.flowAnglesAndRates.αₑ*180/π
end
 
# Set paths
relPath = "/dev/helios/figures/heliosWingSteady"
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

# Normalized deformed wingspan
plt1 = plot(xlabel="\$x_1\$ [m]", ylabel="\$x_3\$ [m]", xlims=[0,40], ylims=[-20,20])
for (i,U) in enumerate(URange)
    plot!(x1_def[i], x3_def[i], lz=(U/.3048), c=:rainbow, lw=lw, label=false,  colorbar_title="Airspeed [ft/s]")
end
display(plt1)
savefig(string(absPath,"/heliosWingSteady_disp.pdf"))

# Tip OOP disp vs airspeed
plt2 = plot(xlabel="Airspeed [ft/s]", ylabel="Tip OOP disp. [m]", xlims=[0,30])
plot!(URange/.3048, tipOOP, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/heliosWingSteady_OOP.pdf"))

# Tip twist vs airspeed
plt3 = plot(xlabel="Airspeed [ft/s]", ylabel="Tip twist [deg]", xlims=[0,30])
plot!(URange/.3048, tipTwist, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/heliosWingSteady_twist.pdf"))

# Tip AoA vs airspeed
plt4 = plot(xlabel="Airspeed [ft/s]", ylabel="Tip AoA [deg]", xlims=[0,30])
plot!(URange/.3048, tipAoA, c=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/heliosWingSteady_aoa.pdf"))

println("Finished heliosWingSteady.jl")