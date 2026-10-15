using AeroBeams

# Airspeed range [m/s]
URange = vcat(1e-2,0.25:0.25:30)

# Pitch angle (the divergence speed can be found even with θ=0)
θ = deg2rad(0.1)

# Discretization
nElem = 20

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
surf = create_AeroSurface(ppSolver=aeroSolver,airfoil=airfoil,c=c,normSparPos=normSparPos,updateAirfoilParameters=false,smallAngles=true)

# Wing
L = 1
GJ = 10
ρA = 1
∞ = 1e16
wing = create_Beam(name="wing",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(∞=∞,GJ=GJ)],I=[inertia_matrix(ρA=ρA,e2=e2)],surface=surf,p0=transform_to_WM([0;0;θ],from=:E321))

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Model
straightWingTorsionalDivergence = create_Model(name="straightWingTorsionalDivergence",beams=[wing],BCs=[clamp],atmosphere=atmosphere)

# Elemental and nodal arclength coordinate
elems = [straightWingTorsionalDivergence.assemblyData.elements[e] for e in 1:nElem]
x1_e = [elem.geometry.x1 for elem in elems]
x1_n = vcat([vcat(elem.assemblyData.r_n1[1],elem.assemblyData.r_n2[1]) for elem in elems]...)

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Initialize outputs
eigensols = Array{EigenSolution}(undef,length(URange))
AoA = Array{Vector{Float64}}(undef,length(URange))
tipAoA = Array{Float64}(undef,length(URange))

# Initialize guess solution
x0 = zeros(0)

# Sweep airspeed
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s")
    # Update airspeed on model
    update_v_A!(straightWingTorsionalDivergence,v_A=[0,U,0])
    # Create and solve eigenproblem
    problem = create_EigenProblem(model=straightWingTorsionalDivergence,eigenmethod=eigenmethod,nev=nev,x0=x0)
    solve!(problem)
    # Eigensolution
    eigensols[i] = problem.eigensolution
    # Angle of attack
    AoA[i] = [rad2deg(problem.loadstepHistory.aeroVariables[end][e].unsteadyKin.αₑ) for e in 1:nElem]
    tipAoA[i] = AoA[i][end]
    # Set initial guess solution for next iteration
    global x0 = problem.systemArrays.x
end

# Apply mode tracking
freqs,damps,_ = mode_tracking(eigensols)

# Get frequencies and dampings of non-oscillatory modes
freqsNonOsc,dampsNonOsc,_ = filter_tracked_modes(freqs,damps; fmin=0, fmax=0)

# Get divergence speed, dynamic pressure and λ
divergenceSpeeds,_,divergenceInds = find_divergence(URange,freqsNonOsc,dampsNonOsc)
UD = divergenceSpeeds[1]
qD = 1/2*ρ*UD^2
λD = sqrt(qD*c*2π*e/GJ)

# Analytical angle of attack and twist - Eq. 4.52 of Hodges and Pierce
x = 0:0.1:L
q = 1/2*ρ*URange.^2
λ = sqrt.(q*c*2π*e/GJ)
αRef = @. tan(λD*L)*sin(λD*x)+cos(λD*x) - 1
if αRef[end] < 0
    αRef ./= minimum(αRef)
else
    αRef ./= maximum(αRef)
end
αTipRef = @. 180/π*θ*(sec(λ*L)-1) + 180/π*θ
θTipRef = @. 180/π*θ*(sec(λ*L)-1)

# Analytical divergence speed - Eq. 4.55 of Hodges and Pierce
qDRef = GJ/(e*c*2π)*(π/(2L))^2
UDRef = sqrt(2*qDRef/ρ)

# Relative error in divergence speed
ϵrel = UD/UDRef - 1
println("UDRef = $(round(UDRef, digits=2)), UD = $(round(UD, digits=2)), ϵrel = $(round(100ϵrel, sigdigits=2))%")

## Plots
# Plot configurations
lw = 2
ms = 4
gr()

# AoA distribution at divergence iminence
iD = divergenceInds[1]
plt_AoA = plot(xlabel="Normalized arclength", ylabel="Normalized angle of attack", title="\$U=$(round(URange[iD]/UD,digits=4))U_D\$", xlims=[0,1])
plot!(x1_e/L, AoA[iD]/maximum(AoA[iD]), c=:black, ls=:solid, lw=lw, label="AeroBeams")
scatter!(x, αRef, c=:black, marker=:circle, ms=ms, label="Hodges and Pierce (2011)")
display(plt_AoA)
savefig(string(absPath,"/straightWingTorsionalDivergence_AoAatUD.pdf"))

# Tip AoA vs U/U_D
plt_tipAoA = plot(xlabel="\$U/U_D\$", ylabel="Tip angle of attack [deg]", xlims=[0,1], ylims=[0,5])
plot!(URange[1:iD]/UD,tipAoA[1:iD], c=:black, ls=:solid, lw=lw, label="AeroBeams")
scatter!(URange[1:iD]/UDRef,αTipRef[1:iD], c=:black, marker=:circle, ms=ms, label="Hodges and Pierce (2011)")
display(plt_tipAoA)
savefig(string(absPath,"/straightWingTorsionalDivergence_AoA.pdf"))

# Tip twist angle vs q/q_D
plt_tipTwist = plot(xlabel="\$q/q_D\$", ylabel="Tip twist angle [deg]", xlims=[0,1], ylims=[0,5], legend=:topleft)
plot!((URange[1:iD]/UD).^2,tipAoA[1:iD].-θ*180/π, c=:black, ls=:solid, lw=lw, label="AeroBeams")
scatter!((URange[1:iD]/UDRef).^2,θTipRef[1:iD], c=:black, marker=:circle, ms=ms, label="Hodges and Pierce (2011)")
display(plt_tipTwist)
savefig(string(absPath,"/straightWingTorsionalDivergence_twist.pdf"))

# Damping plot
plt_damp = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", ylims=[-20,5],legend=:topleft)
for mode in axes(dampsNonOsc,1)
    plot!(URange, dampsNonOsc[mode,:], c=:black, lw=lw, label=false)
end
display(plt_damp)
savefig(string(absPath,"/straightWingTorsionalDivergence_damp.pdf"))

println("Finished straightWingTorsionalDivergence.jl")