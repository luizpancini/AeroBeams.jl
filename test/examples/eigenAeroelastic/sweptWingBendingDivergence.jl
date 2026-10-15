using AeroBeams

# Airspeed range [m/s]
URange = collect(0.2:0.2:20)

# Pitch angle (the divergence speed can be found even with θ=0)
θ = deg2rad(0)

# Sweep angle [rad]
Λ = deg2rad(-30)

# Discretization
nElem = 10

# Aerodynamic solver
aeroSolver = Indicial()

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
d = 2*normCGPos - 1
e = c * (normSparPos - 1/4)
e2 = -(d-a)*c/2
surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=c,normSparPos=normSparPos,Λ=Λ,updateAirfoilParameters=false,smallAngles=true)

# Wing
L = 1
EIy = 10
ρA = 1
∞ = 1e9
wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,EIy=EIy)],I=[inertia_matrix(ρA=ρA,e2=e2)],surface=surf,p0=transform_to_WM([-Λ,0,θ],from=:E321))

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
sweptWingBendingDivergence = create_Model(name="sweptWingBendingDivergence",beams=[wing],BCs=[clamp],atmosphere=atmosphere)

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 16

# Initialize outputs
problem = Array{EigenProblem}(undef,length(URange))
eigensols = Array{EigenSolution}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update airspeed on model
    update_v_A!(sweptWingBendingDivergence,v_A=[0,U,0])
    # Create and solve eigenproblem
    problem[i] = create_EigenProblem(model=sweptWingBendingDivergence,eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(eigensols)

# Get frequencies and dampings of non-oscillatory modes
freqsNonOsc,dampsNonOsc,_ = filter_tracked_modes(freqs,damps; fmin=0, fmax=0)

# Get divergence speed
divergenceSpeeds,_ = find_divergence(URange,freqsNonOsc,dampsNonOsc)
UD = minimum(divergenceSpeeds)

# Approximate divergence speed - Eq. 4.104 of Hodges and Pierce
βD = -19/3
qDRef = βD*EIy/(c*2π*L^3*sin(Λ)*cos(Λ))
UDRef = sqrt(2*qDRef/ρ)

# Relative error in divergence speed
ϵrel = UD/UDRef - 1
println("UDRef = $(round(UDRef,digits=2)), UD = $(round(UD,digits=2)), ϵrel = $(round(100ϵrel,sigdigits=2))%")

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptWingBendingDivergence"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
ms = 4
gr()

# Damping plot
plt_damp = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", xlims=extrema(URange), ylims=[-50,5],legend=:topleft)
for mode in axes(dampsNonOsc,1)
    plot!(URange, dampsNonOsc[mode,:], c=:black, lw=lw, label=false)
end
display(plt_damp)
savefig(string(absPath,"/sweptWingBendingDivergence.pdf"))

println("Finished sweptWingBendingDivergence.jl")