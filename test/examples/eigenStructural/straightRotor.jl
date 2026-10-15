using AeroBeams

# Geometric properties
L = 37.5
b,H = 1, 0.117
A,Iy,Iz = b*H, b*H^3/12, H*b^3/12
J = Iy + Iz
Ksy,Ksz,Kt = 5/6, 1/14.625, 1/65.852
r0 = [2.5; 0.0; 0.0]

# Material properties
E11,E22,E33 = 2.059e7, 1.42e6, 1.42e6
G12,G13,G23 = 8.9e5, 8.9e5, 8e5
ν12,ν13,ν23 = 0.42, 0.42, 0.54
ρ = 1.44e-4

# Discretization variables
nElem = 30

# Range of angular velocity [rad/s]
ωRPM = collect(0:50:750)
ωRange = rpm2rads.(ωRPM)

# Beam
beam = create_Beam(name="beam",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=E11*A,GAy=G12*A*Ksy,GAz=G13*A*Ksz,GJ=G23*J*Kt,EIy=E11*Iy,EIz=E11*Iz)],I=[inertia_matrix(ρ=ρ,A=A,Iy=Iy,Iz=Iz)])

# BCs
clamp = create_BC(name="clamp",beam=beam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

# Create model
straightRotor = create_Model(name="straightRotor",beams=[beam],BCs=[clamp],firstNodePosition=r0,units=create_UnitsSystem(length="in",force="lbf",frequency="Hz"))

# Number of modes to track
nModes = 6

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 20

# Initialize outputs
problem = Array{EigenProblem}(undef,length(ωRange))
eigensols = Array{EigenSolution}(undef,length(ωRange))

# Loop over angular velocities
for (i,ω) in enumerate(ωRange)
    # Update angular velocity of basis A
    update_ω_A!(straightRotor,ω_A=[0, 0, ω])
    # Create and solve eigenproblem
    problem[i] = create_EigenProblem(model=straightRotor,eigenmethod=eigenmethod,nev=nev)
    solve!(problem[i])
    # Eigensolution
    eigensols[i] = problem[i].eigensolution
end

# Apply mode tracking
freqs,_ = mode_tracking(eigensols)

# Get frequencies of oscillatory modes
freqs,_ = filter_tracked_modes(rads2Hz.(freqs); fmin=eps(), maxTracks=nModes)

# Load experimental values from Epps & Chandra (1996)
ωRPMExp = [0, 250, 500, 750]
expFreqs = [5.2374e+00   7.4239e+00   1.1018e+01   1.5820e+01;
            3.2999e+01   3.4582e+01   4.0188e+01   4.7001e+01;
            9.0439e+01   9.1112e+01   9.7733e+01   1.0558e+02;
            1.7237e+02   1.7532e+02   1.7967e+02   1.8874e+02]

## Plots
# Set paths
relPath = "/test/outputs/figures/straightRotor"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes 
modesPlot = plot_mode_shapes(problem[end],scale=5,view=(30,30),nModes=nModes,save=true,savePath=string(relPath,"/straightRotor_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
gr()
lw = 2
ms = 5
msw = 0
colors = cgrad(:darkrainbow, 4, categorical=true)
labels = ["1B", "2B", "3B", "4B"]

# Plot frequency over angular velocity for several modes
plt1 = plot(xlabel="Angular velocity [rpm]", ylabel="Frequency [Hz]", title="Bending modes", xticks=ωRPMExp, yticks=collect(0:20:200))
modes = [1,2,5,6]
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Epps & Chandra (1996)")
for (m,mode) in enumerate(modes)  
    plot!(ωRPM, freqs[mode,:], lc=colors[m], lw=lw, label=false)
    scatter!(ωRPMExp,expFreqs[m,:], mc=colors[m], ms=ms, msw=msw, label=false)
    plot!([NaN], [NaN], lc=colors[m], m=colors[m], lw=lw, ms=ms, msw=msw, label=labels[m])
end
plot!(legend=(0.15,0.8))
display(plt1)
savefig(string(absPath,"/straightRotor_freqVsOmega.pdf"))

println("Finished straightRotor.jl")