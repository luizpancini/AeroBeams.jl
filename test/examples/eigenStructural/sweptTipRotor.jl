# # Modal analysis of a rotating beam
# This example illustrates how to set up a modal (eigen) analysis, using a swept-tip rotor. We will investigate the vibration modes of this rotor at several rotating speeds and several angles of the tip. This problem was proposed by [Hodges, Shang and Cesnik](https://doi.org/10.2514/6.1995-1174), based on the experiments of [Epps and Chandra](https://doi.org/10.4050/JAHS.41.29).

#md # ![](../assets/sweptTipRotor.png)

#md # *Geometry of the swept-tip rotor* by [Hodges, Shang and Cesnik](https://doi.org/10.2514/6.1995-1174)

#md # !!! tip
#md #     The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/eigenStructural/sweptTipRotor.jl).

# ### Problem setup
# To set up the problem, we need to define the geometric and material properties of the beams, to later create them. We also define the number of elements for each beam, the number of vibration modes to extract and the ranges for the angular velocity and the tip angle. Notice that the beam is clamped at a position `r0` relative to the axis of rotation, which is input to the function [`create_Model`](@ref create_Model) through the argument `firstNodePosition`. Moreover, the model assembly assumes by default that all the beams given by the `beams` argument are connected in sequence (*e.g.*, the first node of `tipBeam` is the same as the last of `baseBeam`).
using AeroBeams

## Geometric properties
L1,L2 = 31.5, 6.0
b,H = 1, 0.063
A,Iy,Iz = b*H, b*H^3/12, H*b^3/12
J = Iy + Iz
Ksy,Ksz,Kt = 5/6, 1/14.625, 1/65.852
r0 = [2.5; 0.0; 0.0]

## Material properties (aluminum)
E = 1.06e7
G = E/(2*(1+0.325))
ρ = 2.51e-4

## Sectional stiffness and inertia matrices
S = isotropic_stiffness_matrix(EA=E*A,GAy=G*A*Ksy,GAz=G*A*Ksz,GJ=G*J*Kt,EIy=E*Iy,EIz=E*Iz)
I = inertia_matrix(ρ=ρ,A=A,Iy=Iy,Iz=Iz)

## Discretization variables
nElem1 = 20
nElem2 = 10

## Range of angular velocity [rad/s]
ωRange = rpm2rads.([0,500,750])

## Range of tip beam angles [rad]
tipAngleRange = deg2rad.(0:2.5:45)

## Number of modes to track
nModes = 8

## Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :fulldecomp
nev = 30

## Beam frame
baseBeam = create_Beam(name="baseBeam",length=L1,nElements=nElem1,S=[S],I=[I])
tipBeam = create_Beam(name="tipBeam",length=L2,nElements=nElem2,S=[S],I=[I])

## BCs
clamp = create_BC(name="clamp",beam=baseBeam,node=1,conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))

## Create model
sweptTipRotor = create_Model(name="sweptTipRotor",beams=[baseBeam,tipBeam],BCs=[clamp],firstNodePosition=r0,units=create_UnitsSystem(length="in",force="lbf",frequency="Hz"))

## Initialize outputs
problem = Array{EigenProblem}(undef,length(ωRange),length(tipAngleRange))
eigensols = Array{EigenSolution}(undef,length(ωRange),length(tipAngleRange))
freqs = Array{Matrix{Float64}}(undef,length(ωRange))
#md nothing #hide

# ### Solving the problem
# We now solve an eigenproblem at each combination of tip angle and angular velocity, and extract the modal frequencies. Notice that we set the angular velocity of the entire model through the angular velocity of basis ``\mathcal{A}``, using the variable `ω_A`. AeroBeams will automatically solve the steady problem to find the deformed state of the model, and then perform an eigen-analysis around that state.
## Loop over sweep variables
for (i,ω) in enumerate(ωRange)
    ## Set angular velocity of basis A 
    update_ω_A!(sweptTipRotor,ω_A=[0; 0; ω])
    ## Loop tip angles
    for (j,tipAngle) in enumerate(tipAngleRange)
        ## Display progress #src
        display("Solving for ω = $(round(Int,rads2rpm(ω))) rpm and tip angle = $(round(rad2deg(tipAngle),digits=1)) deg") #src
        ## Update tipBeam orientation angle
        update_beam_orientation!(sweptTipRotor; beam=tipBeam,p0=transform_to_WM([-tipAngle,0,0],from=:E321))
        ## Create and solve eigenproblem
        problem[i,j] = create_EigenProblem(model=sweptTipRotor,eigenmethod=eigenmethod,nev=nev)
        solve!(problem[i,j])
        ## Eigensolution
        eigensols[i,j] = problem[i,j].eigensolution
    end
    ## Apply mode tracking
    freqs[i],_ = mode_tracking(eigensols[i,:])
    ## Get frequencies of oscillatory modes in Hz
    freqs[i],_ = filter_tracked_modes(rads2Hz.(freqs[i]); fmin=eps(), maxTracks=nModes)
end
#md nothing #hide

# ### Post-processing
# We begin by loading the experimental values by [Epps and Chandra](https://doi.org/10.4050/JAHS.41.29).
## Experimental values
expTipAngles = [0, 15, 30, 45]
expFreqs1 = [1.4 1.8 1.7 1.6; 10.2 10.1 10.2 10.2; 14.8 14.4 14.9 14.7]
expFreqs2 = [10.3 10.2 10.4 10.4; 25.2 25.2 23.7 21.6; 36.1 34.8 30.7 26.1]
expFreqs3 = [27.7 27.2 26.6 24.8; 47.0 44.4 39.3 35.1; 62.9 55.9 48.6 44.8]
expFreqs4 = [95.4 87.5 83.7 78.8; 106.6 120.1 122.6 117.7; 132.7 147.3 166.2 162.0]
#md nothing #hide

# Now we plot the evolution of the modal frequencies with the tip angle for the range of angular velocities. The numerical results agree well with the experimental ones, and match those described by [Hodges, Shang and Cesnik](https://doi.org/10.2514/6.1995-1174) and [Yu and Blair](https://doi.org/10.1016/j.compstruct.2012.04.007).
#md gr()
#md ENV["GKSwstype"] = "100" #hide
#md colors = cgrad(:darkrainbow, length(ωRange), categorical=true)
#md lw = 2
#md ms = 5
#md msw = 0

## 1st bending mode
#md plt1 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]", title="1st bending", legend=(0.15,0.4))
#md plot!([NaN],[NaN], lc=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Epps & Chandra (1996)")
#md for (i,ω) in enumerate(ωRange) 
#md     mode = 1
#md     plot!(rad2deg.(tipAngleRange),freqs[i][mode,:], lc=colors[i], lw=lw, label=false)
#md     scatter!(expTipAngles,expFreqs1[i,:], mc=colors[i], ms=ms, msw=msw, label=false)
#md     plot!([NaN], [NaN], lc=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label="\$\\omega\$ = $(round(Int,rads2rpm(ω))) rpm")
#md end
#md savefig("sweptTipRotor_1B.svg") #hide
#md nothing #hide

## 2nd bending
#md plt2 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]", title="2nd bending", legend=(0.15,0.4))
#md plot!([NaN],[NaN], lc=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Epps & Chandra (1996)")
#md for (i,ω) in enumerate(ωRange) 
#md     mode = i < 3 ? 2 : 3
#md     plot!(rad2deg.(tipAngleRange),freqs[i][mode,:], lc=colors[i], lw=lw, label=false)
#md     scatter!(expTipAngles,expFreqs2[i,:], mc=colors[i], ms=ms, msw=msw, label=false)
#md     plot!([NaN], [NaN], lc=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label="\$\\omega\$ = $(round(Int,rads2rpm(ω))) rpm")
#md end
#md savefig("sweptTipRotor_2B.svg") #hide
#md nothing #hide

## 3rd bending mode
#md plt3 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]", title="3rd bending", legend=:topright)
#md plot!([NaN],[NaN], lc=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Epps & Chandra (1996)")
#md for (i,ω) in enumerate(ωRange) 
#md     mode = 4
#md     plot!(rad2deg.(tipAngleRange),freqs[i][mode,:], lc=colors[i], lw=lw, label=false)
#md     scatter!(expTipAngles,expFreqs3[i,:], mc=colors[i], ms=ms, msw=msw, label=false)
#md     plot!([NaN], [NaN], lc=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label="\$\\omega\$ = $(round(Int,rads2rpm(ω))) rpm")
#md end
#md savefig("sweptTipRotor_3B.svg") #hide
#md nothing #hide

## Coupled bending-torsion modes at ω = 750 rpm
#md plt4 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]", title="Coupled bending-torsion at \$\\omega\$ = 750 rpm", legend=(0.1,0.99))
#md plot!([NaN],[NaN], lc=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Epps & Chandra (1996)")
#md modes = [5,6,7]
#md modeLabels = ["1T/5B","5B/1T","4B/1T"]
#md for (i,m) in enumerate(modes)
#md     plot!(rad2deg.(tipAngleRange),freqs[3][m,:], c=colors[i], lw=lw, label=false)
#md     scatter!(expTipAngles,expFreqs4[i,:], c=colors[i], ms=ms, msw=msw, label=false)
#md     plot!([NaN], [NaN], lc=colors[i], mc=colors[i], lw=lw, ms=ms, msw=msw, label=modeLabels[i])
#md end
#md savefig("sweptTipRotor_TB.svg") #hide
#md nothing #hide

#md # ![](sweptTipRotor_1B.svg)
#md # ![](sweptTipRotor_2B.svg)
#md # ![](sweptTipRotor_3B.svg)
#md # ![](sweptTipRotor_TB.svg)

## Plots
# Set paths
relPath = "/test/outputs/figures/sweptTipRotor"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot mode shapes for the final configuration
modesPlot = plot_mode_shapes(problem[end],nModes=nModes,scale=5,view=(30,30),legendPos=:best,save=true,savePath=string(relPath,"/sweptTipRotor_modeShapes.pdf"))
display(modesPlot)

# Plot configurations
colors = cgrad(:darkrainbow, length(ωRange), categorical=true)
lw = 2
ms = 5
msw = 0
gr()

# Plot 1st bending mode frequency over tip angle for several angular velocities
mode = 1
plt1 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]", title="1st bending", legend=(0.15,0.4))
plot!([NaN],[NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Epps & Chandra (1996)")
for (i,ω) in enumerate(ωRange)
    plot!(rad2deg.(tipAngleRange),freqs[i][mode,:], lc=colors[i], lw=lw, label=false)
    scatter!(expTipAngles,expFreqs1[i,:], mc=colors[i], ms=ms, msw=msw, label=false)
    plot!([NaN], [NaN], lc=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label="\$\\omega\$ = $(round(Int,rads2rpm(ω))) rpm")
end
display(plt1)
savefig(string(absPath,"/sweptTipRotor_1B.pdf"))

# Plot 2nd bending mode frequency over tip angle for several angular velocities
plt2 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]", title="2nd bending", legend=(0.15,0.4))
plot!([NaN],[NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Epps & Chandra (1996)")
for (i,ω) in enumerate(ωRange)
    mode = i < 3 ? 2 : 3
    plot!(rad2deg.(tipAngleRange),freqs[i][mode,:], lc=colors[i], lw=lw, label=false)
    scatter!(expTipAngles,expFreqs2[i,:], mc=colors[i], ms=ms, msw=msw, label=false)
    plot!([NaN], [NaN], lc=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label="\$\\omega\$ = $(round(Int,rads2rpm(ω))) rpm")
end
display(plt2)
savefig(string(absPath,"/sweptTipRotor_2B.pdf"))

# Plot 3rd bending mode frequency over tip angle for several angular velocities
mode = 4
plt3 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]", title="3rd bending", legend=:topright)
plot!([NaN],[NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Epps & Chandra (1996)")
for (i,ω) in enumerate(ωRange) 
    plot!(rad2deg.(tipAngleRange),freqs[i][mode,:], lc=colors[i], lw=lw, label=false)
    scatter!(expTipAngles,expFreqs3[i,:], mc=colors[i], ms=ms, msw=msw, label=false)
    plot!([NaN], [NaN], lc=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label="\$\\omega\$ = $(round(Int,rads2rpm(ω)) ) rpm")
end
display(plt3)
savefig(string(absPath,"/sweptTipRotor_3B.pdf"))

# Plot frequency over tip angle for coupled bending-torsion modes at ω = 750 rpm
plt4 = plot(xlabel="Tip sweep angle [deg]", ylabel="Frequency [Hz]", title="Coupled bending-torsion at \$\\omega\$ = 750 rpm", legend=(0.1,0.99))
plot!([NaN],[NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Epps & Chandra (1996)")
modes = [5,6,7]
modeLabels = ["1T/5B","5B/1T","4B/1T"]
for (i,m) in enumerate(modes)
    plot!(rad2deg.(tipAngleRange),freqs[3][m,:], c=colors[i], lw=lw, label=false)
    scatter!(expTipAngles,expFreqs4[i,:], c=colors[i], ms=ms, msw=msw, label=false)
    plot!([NaN], [NaN], lc=colors[i], mc=colors[i], lw=lw, ms=ms, msw=msw, label=modeLabels[i])
end
display(plt4)
savefig(string(absPath,"/sweptTipRotor_TB.pdf"))

println("Finished sweptTipRotor.jl") #src