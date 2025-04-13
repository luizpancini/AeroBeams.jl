using AeroBeams

# Sweep angle range
ΛRange = π/180*[-30; 0; 30]

# Corresponding pitch angle range (in order to obtain approximately same lift)
θRange = π/180*[0.95; 1; 1.35]

# Airspeed
U = 30

# Wing span and chord
L = 1
chord = 0.15

# Wing spar properties
EIy = 100
GJ = 1*EIy
∞ = 1e12
S = isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy)

# Tip loss option
hasTipCorrection = true
tipLossDecayFactor = 8

# Discretization
nElem = 100

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(ΛRange))
x1_e = Array{Vector{Float64}}(undef,length(ΛRange))
cn = Array{Vector{Float64}}(undef,length(ΛRange))

# Loop over sweep angles
for (i,Λ,θ) in zip(1:length(ΛRange),ΛRange,θRange)
    # Wing surface
    wingSurf = create_AeroSurface(solver=Indicial(),airfoil=deepcopy(flatPlate),c=chord,normSparPos=1/4,Λ=Λ,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor)
    # Wing beam
    beam = create_Beam(name="beam",length=L,nElements=nElem,S=[S],aeroSurface=wingSurf,rotationParametrization="E321",p0=[-Λ;0;θ])
    # BCs
    clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
    # Wing model
    wing = create_Model(name="wing",beams=[beam],BCs=[clamp],v_A=[0;U;0])
    # Create and solve steady problem
    problem[i] = create_SteadyProblem(model=wing)
    solve!(problem[i])
    # Undeformed elemental positions
    x1_e[i] = [wing.beams[1].elements[e].x1 for e in 1:nElem]
    # Normal force over the span
    cn[i] = [problem[i].model.elements[e].aero.aeroCoefficients.cn for e in 1:nElem]
end

using Plots

# Set paths
relPath = "/test/outputs/figures/sweptWingsLiftDistribution"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
gr()
lw = 2
colors = cgrad(:rainbow, length(ΛRange), categorical=true)
labels = ["Forward swept" "Unswept" "Backward swept"]

# Normalized normal force coefficient over span (lift per unit span) - compare to Fig. 4.22 of Hodges & Pierce
plt_cn = plot(xlabel="\$x_1/L\$", ylabel="Normalized lift per unit span", xlims=[0,1], ylims=[0,Inf])
for (i,Λ) in enumerate(ΛRange)
    plot!(x1_e[i]/L, cn[i]/cn[2][1], c=colors[i], lw=lw, label=labels[i])
end
display(plt_cn)
savefig(string(absPath,"/sweptWingsLiftDistribution_cn.pdf"))

println("Finished sweptWingsLiftDistribution.jl")