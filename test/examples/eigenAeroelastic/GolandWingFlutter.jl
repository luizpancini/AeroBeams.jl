# Goland wing parameters taken from PALACIOS & EPUREANU - An Intrinsic Description of the Nonlinear Aeroelasticity of Very Flexible Wings - (2011) (https://doi.org/10.2514/6.2011-1917)

using AeroBeams

# Airspeed range[m/s]
URange = vcat(1e-2,5:5:160)

# Altitude [m] and atmosphere
altitude = feet2meters(6.125e3)
atmosphere = standard_atmosphere(altitude)

# Pitch angle
θ = deg2rad(0)

# Gravity [m/s²]
g = 0

# Tip loss factor
τ = Inf

# Aerodynamic solver
aeroSolver = Inflow()

# Aerodynamic surface
airfoil = deepcopy(FlatPlate)
chord = feet2meters(6)
normSparPos = 0.33
normCGPos = 0.43
surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,hasTipCorrection=true,tipLossDecayFactor=τ,updateAirfoilParameters=true)

# Wing beam
L = feet2meters(20)
EIy = 9.77e6
GJ = 0.99e6
ρA = 35.71
ρIx = 8.64
σ = 0.1
ρIy,ρIz = σ*ρIx, (1-σ)*ρIx
e2 = -(normCGPos-normSparPos)*chord
nElem = 20
beam = create_Beam(name="wingBeam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(GJ=GJ,EIy=EIy)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,e2=e2)],p0=transform_to_WM([0,0,θ],from=:E321),surface=surf)

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
GolandWing = create_Model(name="GolandWing",beams=[beam],BCs=[clamp],gravityVector=[0;0;g],atmosphere=atmosphere)

# Number of modes to track
nModes = 2

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
nev = 80

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(URange))

# Initial guess solution
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Update airspeed on model
    update_v_A!(GolandWing,v_A=[0,U,0])
    # Create and solve problem
    problem = create_EigenProblem(model=GolandWing,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(problem)
    # Eigensolution
    eigensols[i] = problem.eigensolution
    # Set initial guess solution for next iteration
    global x0 = problem.systemArrays.x
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(eigensols)

# Get frequencies and dampings of oscillatory modes
freqsOsc,dampsOsc,_ = filter_tracked_modes(freqs,damps; fmin=eps(), fmax=100, maxTracks=nModes, dropOutsideBounds=true)

# Flutter speed and frequency
flutterSpeeds,flutterFreqs,_ = find_flutter(URange,freqsOsc,dampsOsc)
Uf = flutterSpeeds[1]
ωf = flutterFreqs[1]

# Reference results
UfRef = 140
ωfRef = 69.9

# Relative errors (%)
ϵU = 100*(Uf/UfRef - 1)
ϵF = 100*(ωf/ωfRef - 1)

# Pretty table of results
data = [
    "Flutter speed"      Uf  UfRef  ϵU
    "Flutter frequency"  ωf  ωfRef  ϵF
]
pretty_table(
    data;
    column_labels = ["", "AeroBeams", "Palacios & Epureanu", "Δ [%]"],
    title         = "Goland Wing Flutter",
    formatters    = [
        (v, i, j) -> (v isa Number) ? round(v, digits=1) : v,
        (v, i, j) -> (j == 4 && v isa Number) ? string(round(v, digits=2), " %") : v,
    ],
    source_notes  = "Reference: PALACIOS & EPUREANU - An Intrinsic Description of the Nonlinear Aeroelasticity of Very Flexible Wings - (2011)",
)

## Plots
# Set paths
relPath = "/test/outputs/figures/GolandWingFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
lw = 2
ms = 2
gr()

# V-g-f
plt11 = plot(ylabel="Frequency [rad/s]")
for mode in 1:nModes
    plot!(URange, freqsOsc[mode,:], c=modeColors[mode], lw=lw,  label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.25,0.1],legend=:topleft)
for mode in 1:nModes
    plot!(URange, dampsOsc[mode,:]./freqsOsc[mode,:], c=modeColors[mode], lw=lw, label=false)
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(absPath,"/GolandWingFlutter_Vgf.pdf"))

println("Finished GolandWingFlutter.jl")