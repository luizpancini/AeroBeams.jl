using AeroBeams

# Airspeed range [m/s]
URange = collect(0.1:0.1:40)

# Aerodynamic solver
aeroSolver = Indicial()

# Number of elements for the wing
nElem = 2

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)
ρ = atmosphere.ρ

# Aerodynamic surface
airfoil = deepcopy(FlatPlate)
c = 0.2
normSparPos = 1/2
normCGPos = 1/2
a = 2*normSparPos - 1
e = 2*normCGPos - 1
e2 = -(e-a)*c/2
surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=c,normSparPos=normSparPos,updateAirfoilParameters=false,smallAngles=true)

# Rigid wing
L = 1
GJ = 30
ρA = 1
∞ = 1e9
wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=ρA,e2=e2)],surface=surf)

# Attachment spring in pitch
kα = GJ/L
spring = create_Spring(kp=kα)
add_spring_to_beam!(beam=wing,spring=spring,nodeID=div(nElem,2)+1)

# BCs
journal1 = create_BC(name="journal-1",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p2A=>0,:p3A=>0))
journal2 = create_BC(name="journal-2",beam=wing,node=nElem+1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p2A=>0,:p3A=>0))

# Model
typicalSectionDivergence = create_Model(name="typicalSectionDivergence",beams=[wing],BCs=[journal1,journal2],atmosphere=atmosphere)

# Number of modes to track
nModes = 10

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 16

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(URange))

# Initialize guess solution
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update airspeed on model
    update_v_A!(typicalSectionDivergence,v_A=[0,U,0])
    # Create and solve eigenproblem
    problem = create_EigenProblem(model=typicalSectionDivergence,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(problem)
    # Eigensolution
    eigensols[i] = problem.eigensolution
    # Set initial guess solution for next iteration
    global x0 = problem.systemArrays.x
end

# Apply mode tracking
freqs,damps,modeShapes,_ = mode_tracking(eigensols)

# Get frequencies and dampings of non-oscillatory modes
freqsNonOsc,dampsNonOsc,_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=0, fmax=0, maxTracks=nModes)

# Get divergence speed
divergenceSpeeds,_ = find_divergence(URange,freqsNonOsc,dampsNonOsc)
UD = minimum(divergenceSpeeds)

# Analytical divergence speed - Eq. 4.8 of Hodges and Pierce
UDRef = sqrt(kα/(1/2*ρ*c*L*2π*c*(normSparPos-1/4)))

# Relative error in divergence speed
ϵrel = UD/UDRef - 1
println("Relative error in divergence speed = $(round(ϵrel, sigdigits=3))")

## Plots
# Set paths
relPath = "/test/outputs/figures/typicalSectionDivergence"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
gr()

# Damping plot
plt_damp = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", xlims=extrema(URange), ylims=[-50,5],legend=:topleft)
for mode in 1:3
    plot!(URange, dampsNonOsc[mode,:], c=:black, lw=lw, label=false)
end
display(plt_damp)
savefig(string(absPath,"/typicalSectionDivergence.pdf"))

println("Finished typicalSectionPitchDivergence.jl")