using AeroBeams

# Wing spar properties
L = 1
GJ = 5
EIy = 10
ρA = 1
∞ = 1e9

# Sweep angle needed to avoid infinite divergence speed [deg]
Λ∞ = atand(76*EIy*e/(3*π^2*GJ*L))

# Sweep angle range [rad]
ΛRange = deg2rad.(-60:5:floor(Λ∞-5))

# Airspeed range [m/s]
URange = collect(0.5:0.5:30)

# Pitch angle (the divergence speed can be found even with θ=0)
θ = deg2rad(0)

# Discretization
nElem = 20

# Aerodynamic solver
aeroSolver = Indicial()

# Atmosphere 
altitude = 0
atmosphere = standard_atmosphere(altitude)
ρ = atmosphere.ρ

# Aerodynamic surface properties
airfoil = deepcopy(FlatPlate)
c = 0.2
normSparPos = 1/2
normCGPos = 1/2
a = 2*normSparPos - 1
d = 2*normCGPos - 1
e = c * (normSparPos - 1/4)
e2 = -(d-a)*c/2

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 16

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(ΛRange),length(URange))
UD = Array{Float64}(undef,length(ΛRange))

# Sweep angle of sweep
for (i,Λ) in enumerate(ΛRange)
    # Aero surface
    surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=c,normSparPos=normSparPos,Λ=Λ,updateAirfoilParameters=false,smallAngles=true)
    # Wing
    wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy)],I=[inertia_matrix(ρA=ρA,e2=e2)],surface=surf,p0=transform_to_WM([-Λ,0,θ],from=:E321))
    # BCs
    clamp = create_BC(name="clamp",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
    # Model
    model = create_Model(name="wingModel",beams=[wing],BCs=[clamp],atmosphere=atmosphere)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for Λ = $(round(Int,rad2deg(Λ))) deg, U = $U m/s")
        # Update airspeed on model
        update_v_A!(model,v_A=[0,U,0])
        # Create and solve eigenproblem
        problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev)
        solve!(problem)
        # Eigensolution
        eigensols[i,j] = problem.eigensolution
    end
    # Apply mode tracking
    freqs,damps,_ = mode_tracking(eigensols[i,:])
    # Get frequencies and dampings of non-oscillatory modes
    freqsNonOsc,dampsNonOsc,_ = filter_tracked_modes(freqs,damps; fmin=0, fmax=0)
    # Get divergence speed
    divergenceSpeeds,_ = find_divergence(URange,freqsNonOsc,dampsNonOsc)
    UD[i] = minimum(divergenceSpeeds)
    println("UD = $(round(UD[i],digits=2)) m/s at Λ = $(round(Int,rad2deg(Λ))) deg")
end

# Approximate divergence speed - Eq. 4.113 of Hodges and Pierce
UDRef = @. sqrt(2/ρ*19/3*EIy/(2π*c*L^3)*(1+tan(ΛRange)^2)/(tand(Λ∞)-tan(ΛRange)))

# Maximum relative error
ϵmax = maximum(abs.(UD - UDRef)./UDRef)
println("Maximum relative error in divergence speed: $(round(100ϵmax, sigdigits=2))%")

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptWingBendingTorsionalDivergenceSweepRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 4
gr()

# Divergence speed vs. sweep angle
plt_UD = plot(xlabel="Sweep angle [deg]", ylabel="Divergence speed [m/s]", xlims=[-60,Λ∞], ylims=[0,Inf], xticks=[-60,-45,-30,-15,0], legend=:top, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!(rad2deg.(ΛRange), UD, c=:black, lw=lw, ls=:solid, label="AeroBeams")
scatter!(rad2deg.(ΛRange), UDRef, c=:black, ms=ms, label="Hodges and Pierce (2011)")
display(plt_UD)
savefig(string(absPath,"/sweptWingBendingTorsionalDivergenceSweepRange_UD.pdf"))

# Divergence speed error vs. sweep angle
plt_eps = plot(xlabel="Sweep angle [deg]", ylabel="Error in divergence speed [%]", xlims=[-60,Λ∞], xticks=[-60,-45,-30,-15,0], tickfont=font(ts), guidefont=font(fs))
plot!(rad2deg.(ΛRange), (UD./UDRef .- 1)*100, c=:black, lw=lw, ls=:solid, label=false)
display(plt_eps)
savefig(string(absPath,"/sweptWingBendingTorsionalDivergenceSweepRange_error.pdf"))

println("Finished sweptWingBendingTorsionalDivergenceSweepRange.jl")