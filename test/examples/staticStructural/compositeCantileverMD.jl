# # Static structural analysis of composite beams
# This example illustrates how to set up a static structural analysis of composite beams. We will investigate the deflection of several composite beams under tip loads caused by attached weights. This problem was proposed by [Minguet and Dugundji](https://doi.org/10.2514/3.25255):

#md # ![](../assets/compositeCantileverMD.png)

#md # *Beam geometry* by [Minguet and Dugundji](https://doi.org/10.2514/3.25255)

#md # ![](../assets/compositeCantileverMDTable2.png)

#md # *Laminate properties* by [Minguet and Dugundji](https://doi.org/10.2514/3.25255)

#md # !!! tip
#md #     The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/staticStructural/compositeCantileverMD.jl).

# ### Problem setup
# Let's begin by setting up the variables of our problem. We create three beams, each of a different composite laminate, clamp them ar the root, and apply the tip load. Notice that only beam 1 has a diagonal sectional stiffness matrix, whereas the others have coupling terms.
using AeroBeams

## Properties common to all beams (length, in meters, and number of elements)
L = 0.55
nElem = 22

## Root rotation angles [rad]
θRange = π/4*[0, 1]

## Gravity [m/s²]
g = standard_atmosphere(0).g

## Tip load mass [kg] and weight
m = 0.5
W = m*g

## Beam 1
beam1 = create_Beam(name="beam1",length=L,nElements=nElem,S=[isotropic_stiffness_matrix(EA=3.7e6,GAy=2.6e5,GAz=2.9e5,GJ=0.183,EIy=0.707,EIz=276)])

## Beam 2
beam2 = create_Beam(name="beam2",length=L,nElements=nElem,S=[anisotropic_stiffness_matrix(EA=4e6,GAy=2.6e5,GAz=5.5e5,GJ=0.368,EIy=0.522,EIz=298,c12=-2.7e5,c45=-0.102)])

## Beam 3
beam3 = create_Beam(name="beam3",length=L,nElements=nElem,S=[anisotropic_stiffness_matrix(EA=3.9e6,GAy=1.1e6,GAz=1.2e5,GJ=1.18,EIy=0.983,EIz=290,c14=-522)])

## BCs
clamp1 = create_BC(name="clamp1",beam=beam1,node=1,conditions=Dict(:u1A=>0, :u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipLoad1 = create_BC(name="tipLoad1",beam=beam1,node=nElem+1,conditions=Dict(:F3A=>-W))
clamp2 = create_BC(name="clamp2",beam=beam2,node=1,conditions=Dict(:u1A=>0, :u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipLoad2 = create_BC(name="tipLoad2",beam=beam2,node=nElem+1,conditions=Dict(:F3A=>-W))
clamp3 = create_BC(name="clamp3",beam=beam3,node=1,conditions=Dict(:u1A=>0, :u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0))
tipLoad3 = create_BC(name="tipLoad3",beam=beam3,node=nElem+1,conditions=Dict(:F3A=>-W))

## Models
compositeCantilever1 = create_Model(name="compositeCantilever1",beams=[beam1],BCs=[clamp1,tipLoad1])
compositeCantilever2 = create_Model(name="compositeCantilever2",beams=[beam2],BCs=[clamp2,tipLoad2])
compositeCantilever3 = create_Model(name="compositeCantilever3",beams=[beam3],BCs=[clamp3,tipLoad3])

## Set system solver options
σ0 = 0
σstep = 0.1
NR = create_NLsolver(initLoadFac=σ0,maxLoadFacStep=σstep)

## Initialize outputs: displacements at x1 = 500 mm
σVec = Array{Vector{Float64}}(undef,3,2)
u1_500mm = Array{Vector{Float64}}(undef,3,2)
u2_500mm = Array{Vector{Float64}}(undef,3,2)
u3_500mm = Array{Vector{Float64}}(undef,3,2)
#md nothing #hide

# ### Problem solution
# We can sweep the range of root rotation angle, `θRange`, to find the solution at each value for each beam.
## Loop root angles
for (j,θ) in enumerate(θRange)
    ## Update beams' rotation
    update_beam_orientation!(compositeCantilever1,beam=beam1,p0=transform_to_WM([0;0;θ];from=:E321))
    update_beam_orientation!(compositeCantilever2,beam=beam2,p0=transform_to_WM([0;0;θ];from=:E321))
    update_beam_orientation!(compositeCantilever3,beam=beam3,p0=transform_to_WM([0;0;θ];from=:E321))
    ## Create and solve the problems
    problem1 = create_SteadyProblem(model=compositeCantilever1,systemSolver=NR)
    solve!(problem1)
    problem2 = create_SteadyProblem(model=compositeCantilever2,systemSolver=NR)
    solve!(problem2)
    problem3 = create_SteadyProblem(model=compositeCantilever3,systemSolver=NR)
    solve!(problem3)
    ## Save outputs
    hist1 = problem1.loadstepHistory
    hist2 = problem2.loadstepHistory
    hist3 = problem3.loadstepHistory
    stt1 = hist1.states
    stt2 = hist2.states
    stt3 = hist3.states
    σVec[1,j] = hist1.σ
    σVec[2,j] = hist2.σ
    σVec[3,j] = hist3.σ
    u1_500mm[1,j] = [stt1[i][nElem-1].nodal.u1[1] for i in 1:length(σVec[1,j])]
    u1_500mm[2,j] = [stt2[i][nElem-1].nodal.u1[1] for i in 1:length(σVec[2,j])]
    u1_500mm[3,j] = [stt3[i][nElem-1].nodal.u1[1] for i in 1:length(σVec[3,j])]
    u2_500mm[1,j] = [stt1[i][nElem-1].nodal.u1[2] for i in 1:length(σVec[1,j])]
    u2_500mm[2,j] = [stt2[i][nElem-1].nodal.u1[2] for i in 1:length(σVec[2,j])]
    u2_500mm[3,j] = [stt3[i][nElem-1].nodal.u1[2] for i in 1:length(σVec[3,j])]
    u3_500mm[1,j] = [stt1[i][nElem-1].nodal.u1[3] for i in 1:length(σVec[1,j])]
    u3_500mm[2,j] = [stt2[i][nElem-1].nodal.u1[3] for i in 1:length(σVec[2,j])]
    u3_500mm[3,j] = [stt3[i][nElem-1].nodal.u1[3] for i in 1:length(σVec[3,j])]
end
#md nothing #hide

# ### Post-processing
# Let's load the reference experimental solution by [Minguet and Dugundji](https://doi.org/10.2514/3.25255).
## Load reference solution
u1_b1_th0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b1_th0_u1.txt")
u2_b1_th0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b1_th0_u2.txt")
u3_b1_th0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b1_th0_u3.txt")
u1_b1_th45_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b1_th45_u1.txt")
u2_b1_th45_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b1_th45_u2.txt")
u3_b1_th45_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b1_th45_u3.txt")
u1_b2_th0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b2_th0_u1.txt")
u2_b2_th0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b2_th0_u2.txt")
u3_b2_th0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b2_th0_u3.txt")
u1_b2_th45_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b2_th45_u1.txt")
u2_b2_th45_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b2_th45_u2.txt")
u3_b2_th45_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b2_th45_u3.txt")
u1_b3_th0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b3_th0_u1.txt")
u2_b3_th0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b3_th0_u2.txt")
u3_b3_th0_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b3_th0_u3.txt")
u1_b3_th45_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b3_th45_u1.txt")
u2_b3_th45_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b3_th45_u2.txt")
u3_b3_th45_ref = readdlm(pkgdir(AeroBeams)*"/test/referenceData/compositeCantileverMD/b3_th45_u3.txt")
#md nothing #hide

#md # Finally, we may plot the displacements at ``x_1 = 500`` mm for each beam in each configuration. The correlation with the experimental data is good, and matches the numerical results of [Minguet and Dugundji](https://doi.org/10.2514/3.25255) as well (not shown).
#md gr()
#md ENV["GKSwstype"] = "100" #hide
#md nothing #hide

## Plot configurations
#md labels = ["\$-u_1\$" "\$u_2\$" "\$-u_3\$"]
#md xLabel = "\$-u_1, u_2, -u_3\$ [m]"
#md yLabel = "Load [grams]"
#md colors = [:blue,:orange,:green]
#md lw = 2
#md ms = 5
#md msw = 0

## Beam 1, θ=0⁰
#md plt11 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 1, \$\\theta=0^{\\degree}\$")
#md plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN], [NaN], mc=:black, ms=ms, label=" Experimental")
#md for i=1:3
#md     plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
#md end
#md for (x, c) in zip([-u1_500mm[1,1], u2_500mm[1,1], -u3_500mm[1,1]], colors)
#md     plot!(x, σVec[1,1]*W/g*1e3, c=c, lw=lw, label=false)
#md end
#md for (x, y, c) in zip([u1_b1_th0_ref[1,:], u2_b1_th0_ref[1,:], u3_b1_th0_ref[1,:]], [u1_b1_th0_ref[2,:], u2_b1_th0_ref[2,:], u3_b1_th0_ref[2,:]], colors)
#md     scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
#md end
#md savefig("compositeCantileverMD_b1th0.svg") #hide
#md nothing #hide

## Beam 1, θ=45⁰
#md plt12 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 1, \$\\theta=45^{\\degree}\$")
#md plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN], [NaN], mc=:black, ms=ms, label="  Experimental")
#md for i=1:3
#md     plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
#md end
#md for (x, c) in zip([-u1_500mm[1,2], u2_500mm[1,2], -u3_500mm[1,2]], colors)
#md     plot!(x, σVec[1,2]*W/g*1e3, c=c, lw=lw, label=false)
#md end
#md for (x, y, c) in zip([u1_b1_th45_ref[1,:], u2_b1_th45_ref[1,:], u3_b1_th45_ref[1,:]], [u1_b1_th45_ref[2,:], u2_b1_th45_ref[2,:], u3_b1_th45_ref[2,:]], colors)
#md     scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
#md end
#md savefig("compositeCantileverMD_b1th45.svg") #hide
#md nothing #hide

## Beam 2, θ=0⁰
#md plt21 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 2, \$\\theta=0^{\\degree}\$")
#md plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN], [NaN], mc=:black, ms=ms, label="  Experimental")
#md for i=1:3
#md     plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
#md end
#md for (x, c) in zip([-u1_500mm[2,1], u2_500mm[2,1], -u3_500mm[2,1]], colors)
#md     plot!(x, σVec[2,1]*W/g*1e3, c=c, lw=lw, label=false)
#md end
#md for (x, y, c) in zip([u1_b2_th0_ref[1,:], u2_b2_th0_ref[1,:], u3_b2_th0_ref[1,:]], [u1_b2_th0_ref[2,:], u2_b2_th0_ref[2,:], u3_b2_th0_ref[2,:]], colors)
#md     scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
#md end
#md savefig("compositeCantileverMD_b2th0.svg") #hide
#md nothing #hide

## Beam 2, θ=45⁰
#md plt22 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 2, \$\\theta=45^{\\degree}\$", legend=:bottomright)
#md plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN], [NaN], mc=:black, ms=ms, label="  Experimental")
#md for i=1:3
#md     plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
#md end
#md for (x, c) in zip([-u1_500mm[2,2], u2_500mm[2,2], -u3_500mm[2,2]], colors)
#md     plot!(x, σVec[2,2]*W/g*1e3, c=c, lw=lw, label=false)
#md end
#md for (x, y, c) in zip([u1_b2_th45_ref[1,:], u2_b2_th45_ref[1,:], u3_b2_th45_ref[1,:]], [u1_b2_th45_ref[2,:], u2_b2_th45_ref[2,:], u3_b2_th45_ref[2,:]], colors)
#md     scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
#md end
#md savefig("compositeCantileverMD_b2th45.svg") #hide
#md nothing #hide

## Beam 3, θ=0⁰
#md plt31 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 3, \$\\theta=0^{\\degree}\$")
#md plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN], [NaN], mc=:black, ms=ms, label="  Experimental")
#md for i=1:3
#md     plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
#md end
#md for (x, c) in zip([-u1_500mm[3,1], u2_500mm[3,1], -u3_500mm[3,1]], colors)
#md     plot!(x, σVec[3,1]*W/g*1e3, c=c, lw=lw, label=false)
#md end
#md for (x, y, c) in zip([u1_b3_th0_ref[1,:], u2_b3_th0_ref[1,:], u3_b3_th0_ref[1,:]], [u1_b3_th0_ref[2,:], u2_b3_th0_ref[2,:], u3_b3_th0_ref[2,:]], colors)
#md     scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
#md end
#md savefig("compositeCantileverMD_b3th0.svg") #hide
#md nothing #hide

## Beam 3, θ=45⁰
#md plt32 = plot(xlabel=xLabel, ylabel=yLabel, title="Beam 3, \$\\theta=45^{\\degree}\$")
#md plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
#md scatter!([NaN], [NaN], mc=:black, ms=ms, label="  Experimental")
#md for i=1:3
#md     plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
#md end
#md for (x, c) in zip([-u1_500mm[3,2], u2_500mm[3,2], -u3_500mm[3,2]], colors)
#md     plot!(x, σVec[3,2]*W/g*1e3, c=c, lw=lw, label=false)
#md end
#md for (x, y, c) in zip([u1_b3_th45_ref[1,:], u2_b3_th45_ref[1,:], u3_b3_th45_ref[1,:]], [u1_b3_th45_ref[2,:], u2_b3_th45_ref[2,:], u3_b3_th45_ref[2,:]], colors)
#md     scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
#md end
#md savefig("compositeCantileverMD_b3th45.svg") #hide
#md nothing #hide

#md # ![](compositeCantileverMD_b1th0.svg)
#md # ![](compositeCantileverMD_b1th45.svg)
#md # ![](compositeCantileverMD_b2th0.svg)
#md # ![](compositeCantileverMD_b2th45.svg)
#md # ![](compositeCantileverMD_b3th0.svg)
#md # ![](compositeCantileverMD_b3th45.svg)

## Plots
# Set paths
relPath = "/test/outputs/figures/compositeCantileverMD"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
labels = ["\$-u_1\$" "\$u_2\$" "\$-u_3\$"]
xLabel = "Displacements [m]"
yLabel = "Load [grams]"
colors = [:blue,:orange,:green]
ts = 10
fs = 16
lw = 2
ms = 5
msw = 0
gr()

# Beam 1, θ=0⁰
plt11 = plot(xlabel=xLabel, ylabel=yLabel, tickfont=font(ts), guidefont=font(fs), legendfontsize=9)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Exp. - Minguet & Dugundji (1990)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-u1_500mm[1,1], u2_500mm[1,1], -u3_500mm[1,1]], colors)
    plot!(x, σVec[1,1]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b1_th0_ref[1,:], u2_b1_th0_ref[1,:], u3_b1_th0_ref[1,:]], [u1_b1_th0_ref[2,:], u2_b1_th0_ref[2,:], u3_b1_th0_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt11)
savefig(string(absPath,"/compositeCantileverMD_b1th0.pdf"))

# Beam 1, θ=45⁰
plt12 = plot(xlabel=xLabel, ylabel=yLabel, tickfont=font(ts), guidefont=font(fs), legendfontsize=9)
plot!([NaN], [NaN], lc=:black, lw=lw, label=false)
scatter!([NaN], [NaN], mc=:black, ms=ms, label=false)
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=false)
end
for (x, c) in zip([-u1_500mm[1,2], u2_500mm[1,2], -u3_500mm[1,2]], colors)
    plot!(x, σVec[1,2]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b1_th45_ref[1,:], u2_b1_th45_ref[1,:], u3_b1_th45_ref[1,:]], [u1_b1_th45_ref[2,:], u2_b1_th45_ref[2,:], u3_b1_th45_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt12)
savefig(string(absPath,"/compositeCantileverMD_b1th45.pdf"))

# Beam 2, θ=0⁰
plt21 = plot(xlabel=xLabel, ylabel=yLabel, tickfont=font(ts), guidefont=font(fs), legendfontsize=9)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Exp. - Minguet & Dugundji (1990)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-u1_500mm[2,1], u2_500mm[2,1], -u3_500mm[2,1]], colors)
    plot!(x, σVec[2,1]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b2_th0_ref[1,:], u2_b2_th0_ref[1,:], u3_b2_th0_ref[1,:]], [u1_b2_th0_ref[2,:], u2_b2_th0_ref[2,:], u3_b2_th0_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt21)
savefig(string(absPath,"/compositeCantileverMD_b2th0.pdf"))

# Beam 2, θ=45⁰
plt22 = plot(xlabel=xLabel, ylabel=yLabel, tickfont=font(ts), guidefont=font(fs), legendfontsize=9)
plot!([NaN], [NaN], lc=:black, lw=lw, label=false)
scatter!([NaN], [NaN], mc=:black, ms=ms, label=false)
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=false)
end
for (x, c) in zip([-u1_500mm[2,2], u2_500mm[2,2], -u3_500mm[2,2]], colors)
    plot!(x, σVec[2,2]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b2_th45_ref[1,:], u2_b2_th45_ref[1,:], u3_b2_th45_ref[1,:]], [u1_b2_th45_ref[2,:], u2_b2_th45_ref[2,:], u3_b2_th45_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
plot!(legend=:bottomright)
display(plt22)
savefig(string(absPath,"/compositeCantileverMD_b2th45.pdf"))

# Beam 3, θ=0⁰
plt31 = plot(xlabel=xLabel, ylabel=yLabel, tickfont=font(ts), guidefont=font(fs), legendfontsize=9)
plot!([NaN], [NaN], lc=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Exp. - Minguet & Dugundji (1990)")
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=labels[i])
end
for (x, c) in zip([-u1_500mm[3,1], u2_500mm[3,1], -u3_500mm[3,1]], colors)
    plot!(x, σVec[3,1]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b3_th0_ref[1,:], u2_b3_th0_ref[1,:], u3_b3_th0_ref[1,:]], [u1_b3_th0_ref[2,:], u2_b3_th0_ref[2,:], u3_b3_th0_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt31)
savefig(string(absPath,"/compositeCantileverMD_b3th0.pdf"))

# Beam 3, θ=45⁰
plt32 = plot(xlabel=xLabel, ylabel=yLabel, tickfont=font(ts), guidefont=font(fs), legendfontsize=9)
plot!([NaN], [NaN], lc=:black, lw=lw, label=false)
scatter!([NaN], [NaN], mc=:black, ms=ms, label=false)
for i=1:3
    plot!([NaN], [NaN], c=colors[i], m=colors[i], lw=lw, ms=ms, msw=msw, label=false)
end
for (x, c) in zip([-u1_500mm[3,2], u2_500mm[3,2], -u3_500mm[3,2]], colors)
    plot!(x, σVec[3,2]*W/g*1e3, c=c, lw=lw, label=false)
end
for (x, y, c) in zip([u1_b3_th45_ref[1,:], u2_b3_th45_ref[1,:], u3_b3_th45_ref[1,:]], [u1_b3_th45_ref[2,:], u2_b3_th45_ref[2,:], u3_b3_th45_ref[2,:]], colors)
    scatter!(x, y, c=c, ms=ms, msw=msw, label=false)
end
display(plt32)
savefig(string(absPath,"/compositeCantileverMD_b3th45.pdf"))

println("Finished compositeCantileverMD.jl") #src