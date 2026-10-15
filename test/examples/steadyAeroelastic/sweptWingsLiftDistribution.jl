using AeroBeams

# Sweep angle range [rad]
ΛRange = deg2rad.([-30; 0; 30])

# Corresponding pitch angle range (in order to obtain approximately same lift)
θRange = deg2rad.([0.95; 1; 1.35])

# Airspeed [m/s]
U = 30

# Wing span and chord
L = 1
chord = 0.15

# Wing spar properties
EIy = 100
GJ = 1*EIy
∞ = 1e12
S = isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy)

# Tip loss options
hasTipCorrection = true
tipLossDecayFactor = 8

# Discretization
nElem = 100

# Initialize outputs
x1_e = Array{Vector{Float64}}(undef,length(ΛRange))
cn = Array{Vector{Float64}}(undef,length(ΛRange))

# Loop over sweep angles
for (i,Λ,θ) in zip(1:length(ΛRange),ΛRange,θRange)
    # Wing surface
    wingSurf = create_AeroSurface(ppSolver=Indicial(),airfoil=deepcopy(FlatPlate),c=chord,normSparPos=1/4,Λ=Λ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor)
    # Wing beam
    beam = create_Beam(name="beam",length=L,nElements=nElem,S=[S],surface=wingSurf,p0=transform_to_WM([-Λ,0,θ],from=:E321))
    # BCs
    clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
    # Wing model
    wing = create_Model(name="wing",beams=[beam],BCs=[clamp],v_A=[0,U,0])
    # Create and solve steady problem
    problem = create_SteadyProblem(model=wing)
    solve!(problem)
    # Undeformed elemental positions
    x1_e[i] = [wing.beams[1].elements[e].geometry.x1 for e in 1:nElem]
    # Normal force over the span
    cn[i] = [problem.loadstepHistory.aeroVariables[end][e].aeroCoefficients.cn for e in 1:nElem]
end

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptWingsLiftDistribution"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
colors = cgrad(:rainbow, length(ΛRange), categorical=true)
labels = ["Forward swept" "Unswept" "Backward swept"]
gr()

# Normalized normal force coefficient over span (lift per unit span) - compare to Fig. 4.22 of Hodges & Pierce
plt_cn = plot(xlabel="\$x_1/L\$", ylabel="Normalized lift per unit span", xlims=[0,1], ylims=[0,Inf])
for (i,Λ) in enumerate(ΛRange)
    plot!(x1_e[i]/L, cn[i]/cn[2][1], c=colors[i], lw=lw, label=labels[i])
end
display(plt_cn)
savefig(string(absPath,"/sweptWingsLiftDistribution_cn.pdf"))

println("Finished sweptWingsLiftDistribution.jl")