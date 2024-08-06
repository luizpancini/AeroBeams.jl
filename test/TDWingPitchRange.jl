using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Wing beam
L = 0.4508
chord = 0.0508
GJ,EIy,EIz = 0.9539,0.4186,18.44
ρA,ρIs,ρIy = 0.2351,0.2056e-4,1e-6
ρIz = ρIy*EIz/EIy
e3 = 1e-2*chord
nElem = 20
∞ = 1e6
wing = create_Beam(name="wing",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs,e3=e3)],rotationParametrization="E321")

# Wing's tip store
tipMass = 0.0417
tipIyy = 0.3783e-5
tipIzz = 0.9753e-4
tipStore = PointInertia(elementID=nElem,η=[L/nElem/2;0;0],mass=tipMass,Iyy=tipIyy,Izz=tipIzz)
add_point_inertias_to_beam!(wing,inertias=[tipStore])

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
g = 9.807
TDWingPitchRange = create_Model(name="TDWingPitchRange",beams=[wing],BCs=[clamp],gravityVector=[0,0,-g])

# Loop configurations
θRange = LinRange(0,90,21) 
freqs = Array{Vector{Float64}}(undef,length(θRange))
tip_u2 = Array{Float64}(undef,length(θRange))
tip_u3 = Array{Float64}(undef,length(θRange))
tip_twist = Array{Float64}(undef,length(θRange))
for (i,θ) in enumerate(θRange)
    display("Solving for root angle = $θ deg")
    # Update beam root rotation angle
    wing.p0 = [0;0;θ*π/180]
    update_beam!(wing)
    update_model!(TDWingPitchRange)
    # Create and solve problem
    global problem = create_EigenProblem(model=TDWingPitchRange,nModes=4,frequencyFilterLimits=[0.1,Inf64],normalizeModeShapes=true)
    solve!(problem)
    # Get outputs
    freqs[i] = problem.frequenciesOscillatory/(2π)
    tip_u2[i] = problem.nodalStatesOverσ[end][nElem].u_n2_b[2]
    tip_u3[i] = problem.nodalStatesOverσ[end][nElem].u_n2_b[3]
    tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*[0; 1; 0]
    tip_twist[i] = asind(Δ[3])
end

# Load reference solutions
u2_exp = readdlm(string(pwd(),"/test/referenceData/TDWingPitchRange/u2_exp.txt"))
u3_exp = readdlm(string(pwd(),"/test/referenceData/TDWingPitchRange/u3_exp.txt"))
th_exp = readdlm(string(pwd(),"/test/referenceData/TDWingPitchRange/th_exp.txt"))
u2_num = readdlm(string(pwd(),"/test/referenceData/TDWingPitchRange/u2_num.txt"))
u3_num = readdlm(string(pwd(),"/test/referenceData/TDWingPitchRange/u3_num.txt"))
th_num = readdlm(string(pwd(),"/test/referenceData/TDWingPitchRange/th_num.txt"))
freqs_exp = readdlm(string(pwd(),"/test/referenceData/TDWingPitchRange/freqs_exp.txt"))
freq1_num = readdlm(string(pwd(),"/test/referenceData/TDWingPitchRange/freq1_num.txt"))
freq2_num = readdlm(string(pwd(),"/test/referenceData/TDWingPitchRange/freq2_num.txt"))
freq4_num = readdlm(string(pwd(),"/test/referenceData/TDWingPitchRange/freq4_num.txt"))
freqs_num = Array{Matrix{Float64}}(undef,3)
freqs_num[1] = freq1_num
freqs_num[2] = freq2_num
freqs_num[3] = freq4_num

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:darkrainbow], LinRange(0, 1, 3))
lw = 2
ms = 3
# Tip flapwise displacement
plt1 = plot(xlabel="Root angle [deg]", ylabel="Tip flapwise displacement [m]", xticks=collect(0:15:90))
plot!(θRange,-tip_u3, c=:black, lw=2, ls=:solid, label="AeroBeams")
plot!(u3_num[1,:],u3_num[2,:], c=:black, lw=2, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!(u3_exp[1,:],u3_exp[2,:], mc=:black, ms=ms, label="Tang & Dowell (2001) - Exp.")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/TDWingPitchRange/TDWingPitchRange_u3.pdf"))
# Tip chordwise displacement
plt2 = plot(xlabel="Root angle [deg]", ylabel="Tip chordwise displacement [m]", xticks=collect(0:15:90))
plot!(θRange,-tip_u2, c=:black, lw=2, ls=:solid, label="AeroBeams")
plot!(u2_num[1,:],u2_num[2,:], c=:black, lw=2, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!(u2_exp[1,:],u2_exp[2,:], mc=:black, ms=ms, label="Tang & Dowell (2001) - Exp.")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/TDWingPitchRange/TDWingPitchRange_u2.pdf"))
# Tip twist
plt3 = plot(xlabel="Root angle [deg]", ylabel="Tip twist [deg]", xticks=collect(0:15:90))
plot!(θRange,-tip_twist, c=:black, lw=2, ls=:solid, label="AeroBeams")
plot!(th_num[1,:],th_num[2,:], c=:black, lw=2, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!(th_exp[1,:],th_exp[2,:], mc=:black, ms=ms, label="Tang & Dowell (2001) - Exp.")
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/TDWingPitchRange/TDWingPitchRange_twist.pdf"))
# Structural frequencies 
plt4 = plot(xlabel="Root angle [deg]", ylabel="Frequency [Hz]", legend=:outertop, xticks=collect(0:15:90))
modeLabels = ["Flapwise bending" "Chordwise bending" "Torsion"]
plot!([NaN],[NaN], c=:black, lw=2, ls=:solid, label="AeroBeams")
plot!([NaN],[NaN], c=:black, lw=2, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!([NaN],[NaN], mc=:black, ms=ms, label="Tang & Dowell (2001) - Exp.")
for (m,mode) in enumerate([1,2,4])
    freqsMode = [freqs[j][mode] for j in 1:length(freqs)]
    plot!(θRange,freqsMode, c=colors[m], lw=2, ls=:solid, label=false)
    plot!(freqs_num[m][1,:],freqs_num[m][2,:], c=colors[m], lw=2, ls=:dash, label=false)
    scatter!(freqs_exp[1,:],freqs_exp[1+m,:], mc=colors[m], ms=ms, label=false)
    annotate!(30, freqsMode[1], text(modeLabels[m], :bottom, colors[m]))
end
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/TDWingPitchRange/TDWingPitchRange_freqs.pdf"))

println("Finished TDWingPitchRange.jl")