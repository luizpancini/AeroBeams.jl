using AeroBeams

# Sweep angle range [rad]
ΛRange = deg2rad.(collect(0:10:40))

# Discretization
nElem = 40

# Beam properties (correction factors from FREUND & KARAKOC - Shear and torsion correction factors of Timoshenko beam model for generic cross sections (2016))
L = 0.5
b, h = 2e-2, 5e-3
A, Iy, Iz = b*h, b*h^3/12, h*b^3/12
J = Iy+Iz
E, G, ρ = 2e11, 7.6923e10, 7.85e3
Ksy = Ksz = 5/6
Kt = 0.198
∞ = 1e16

# Tip force magnitude
F = 1e2

# Initialize output arrays
problem = Array{SteadyProblem}(undef,length(ΛRange))
tip_u3 = Array{Float64}(undef,length(ΛRange))
tip_twist = Array{Float64}(undef,length(ΛRange))

# Sweep angle of sweep
for (i,Λ) in enumerate(ΛRange)
    println("Solving for Λ = $(round(rad2deg(Λ))) deg")
    # Beam
    beam = create_Beam(name="beam",length=L,nElements=nElem,S=[anisotropic_stiffness_matrix(∞=∞,EA=E*A,GAy=G*A*Ksy,GAz=G*A*Ksz,GJ=G*J*Kt/cos(Λ),EIy=E*Iy*cos(Λ)^3,EIz=E*Iz*cos(Λ)^3)],I=[inertia_matrix(ρA=ρ*A,ρIy=ρ*Iy,ρIz=ρ*Iz)],p0=transform_to_WM([-Λ;0;0],from=:E321))
    # beam.S[1][4,5] = beam.S[1][5,4] = EIy*GJ*sin(Λ)*cos(Λ)^2
    # BCs
    clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
    tipForce = create_BC(name="tipForce",beam=beam,node=nElem+1,conditions=Dict(:F3A=>F))
    # Model
    model = create_Model(name="beamModel",beams=[beam],BCs=[clamp,tipForce])
    # Create and solve eigenproblem
    problem[i] = create_SteadyProblem(model=model)
    solve!(problem[i])
    # Outputs
    tipNodeFinalState = problem[i].loadstepHistory.states[end][nElem].nodal
    tip_u3[i] = tipNodeFinalState.u2[3]
    tip_twist[i] = twist_angle(tipNodeFinalState.p2)
end

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptBeamStaticLoading"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot displaced beams
plt_beams = plot_steady_deformations(problem,backendSymbol=:plotlyjs,plotBCs=true,plotUndeformed=false,plotDistLoads=false,plotAxes=true,plotGrid=true,legendEntries=["\$Λ = $(round(Int,rad2deg(Λ))) ^\\circ \$" for Λ in ΛRange],legendPos=:right,view=(80,30),plotLimits=([0,L],[-L*3/4,L/4],[-L/2,L/2]),save=true,savePath=string(relPath,"/beam_shapes.pdf"))
display(plt_beams)

# Plot configurations
lw = 2
ms = 4
gr()

# Tip displacement vs. sweep angle
plt_u3 = plot(xlabel="Sweep angle [deg]", ylabel="Tip \$u_3/L\$")
plot!(rad2deg.(ΛRange), tip_u3/L, c=:black, lw=lw, label=false)
display(plt_u3)
savefig(string(absPath,"/sweptBeamStaticLoadingTipDisp.pdf"))

# Tip twist vs. sweep angle
plt_twist = plot(xlabel="Sweep angle [deg]", ylabel="Tip twist [deg]")
plot!(rad2deg.(ΛRange), tip_twist, c=:black, lw=lw, label=false)
display(plt_twist)
savefig(string(absPath,"/sweptBeamStaticLoadingTipTwist.pdf"))

println("Finished sweptBeamStaticLoading.jl")