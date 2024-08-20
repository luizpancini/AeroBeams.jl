using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Wing surface
chord = 0.0508
normSparPos = 0.5
aeroSolver = Indicial()
derivationMethod = AD()
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=flatPlate,c=chord,normSparPos=normSparPos)

# Wing beam
L = 0.4508
GJ,EIy,EIz = 0.9539,0.4186,18.44
ρA,ρIs,ρIy = 0.2351,0.2056e-4,1e-6
ρIz = ρIy*EIz/EIy
e3 = 1e-2*chord
nElem = 20
∞ = 1e6
wing = create_Beam(name="wing",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs,e3=e3)],rotationParametrization="E321",aeroSurface=surf)

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
h = 0
TDWingAirspeedRange = create_Model(name="TDWingAirspeedRange",beams=[wing],BCs=[clamp],gravityVector=[0,0,-g],altitude=h)

# Loop configurations
θRange = π/180*[1.0; 2.2]
URange = LinRange(0,40,41)
freqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
tip_u3 = Array{Float64}(undef,length(θRange),length(URange))
tip_twist = Array{Float64}(undef,length(θRange),length(URange))
for (i,θ) in enumerate(θRange)
    # Update beam root rotation angle
    wing.p0 = [0;0;θ]
    update_beam!(wing)
    # Loop airspeeds
    for (j,U) in enumerate(URange)
        display("Solving for root angle = $(θ*180/π) deg, U = $(round(Int,U)) m/s")
        # Update velocity of basis A 
        set_motion_basis_A!(model=TDWingAirspeedRange,v_A=[0;U;0])
        # Create and solve problem
        global problem = create_EigenProblem(model=TDWingAirspeedRange,nModes=4,frequencyFilterLimits=[0.1,Inf64],normalizeModeShapes=true)
        solve!(problem)
        # Get outputs
        freqs[i,j] = problem.frequenciesOscillatory/(2π)
        tip_u3[i,j] = problem.nodalStatesOverσ[end][nElem].u_n2[3]
        tip_p = problem.nodalStatesOverσ[end][nElem].p_n2
        R,_ = rotation_tensor_WM(tip_p)
        Δ = R*[0; 1; 0]
        tip_twist[i,j] = asind(Δ[3])
    end
end

# Load reference solutions
u3_1deg_exp = readdlm(string(pwd(),"/test/referenceData/TDWingAirspeedRange/u3_1deg_exp.txt"))
u3_1deg_num = readdlm(string(pwd(),"/test/referenceData/TDWingAirspeedRange/u3_1deg_num.txt"))
u3_2_2deg_exp = readdlm(string(pwd(),"/test/referenceData/TDWingAirspeedRange/u3_2_2deg_exp.txt"))
u3_2_2deg_num = readdlm(string(pwd(),"/test/referenceData/TDWingAirspeedRange/u3_2_2deg_num.txt"))
th_1deg_exp = readdlm(string(pwd(),"/test/referenceData/TDWingAirspeedRange/th_1deg_exp.txt"))
th_1deg_num = readdlm(string(pwd(),"/test/referenceData/TDWingAirspeedRange/th_1deg_num.txt"))
th_2_2deg_exp = readdlm(string(pwd(),"/test/referenceData/TDWingAirspeedRange/th_2_2deg_exp.txt"))
th_2_2deg_num = readdlm(string(pwd(),"/test/referenceData/TDWingAirspeedRange/th_2_2deg_num.txt"))
freqs_ref = readdlm(string(pwd(),"/test/referenceData/TDWingAirspeedRange/freqs.txt"))

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:darkrainbow], LinRange(0, 1, length(θRange)))
colors2 = get(colorschemes[:darkrainbow], LinRange(0, 1, 3))
lw = 2
ms = 3
labels = ["\\theta_{r} = 1.0 deg" "\\theta_{r} = 2.2 deg"]
relPath = "/test/outputs/figures/TDWingAirspeedRange"
absPath = string(pwd(),relPath)
mkpath(absPath)
# Deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/TDWingAirspeedRange_deformation.pdf"))
display(deformationPlot)
# Tip flapwise displacement
gr()
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Tip flapwise displacement [m]")
plot!([NaN], [NaN], lc=:black,  lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN], [NaN], lc=:black,  lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Tang & Dowell (2001) - Exp.")
for i=eachindex(θRange)
    plot!(URange,tip_u3[i,:], c=colors[i], lw=2, ls=:solid, label=false)
    if i==1
        plot!(u3_1deg_num[1,:],u3_1deg_num[2,:], c=colors[i], lw=2, ls=:dash, label=false)
        scatter!(u3_1deg_exp[1,:],u3_1deg_exp[2,:], mc=colors[i], ms=ms, label=false)
        annotate!(30, tip_u3[i,2], text(labels[i], :bottom, :left, colors[i]))
    else
        plot!(u3_2_2deg_num[1,:],u3_2_2deg_num[2,:], c=colors[i], lw=2, ls=:dash, label=false)
        scatter!(u3_2_2deg_exp[1,:],u3_2_2deg_exp[2,:], mc=colors[i], ms=ms, label=false)
        annotate!(35, tip_u3[i,end-1], text(labels[i], :top, :right, colors[i]))
    end
end
display(plt1)
savefig(string(absPath,"/test/outputs/figures/TDWingAirspeedRange/TDWingAirspeedRange_disp.pdf"))
# Tip twist
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]")
plot!([NaN], [NaN], lc=:black,  lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN], [NaN], lc=:black,  lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!([NaN], [NaN ], mc=:black, ms=ms, label="Tang & Dowell (2001) - Exp.")
for i=eachindex(θRange)
    plot!(URange,tip_twist[i,:], c=colors[i], lw=2, ls=:solid, label=false)
    if i==1
        plot!(th_1deg_num[1,:],th_1deg_num[2,:], c=colors[i], lw=2, ls=:dash, label=false)
        scatter!(th_1deg_exp[1,:],th_1deg_exp[2,:], mc=colors[i], ms=ms, label=false)
        annotate!(30, tip_twist[i,2], text(labels[i], :bottom, :left, colors[i]))
    else
        plot!(th_2_2deg_num[1,:],th_2_2deg_num[2,:], c=colors[i], lw=2, ls=:dash, label=false)
        scatter!(th_2_2deg_exp[1,:],th_2_2deg_exp[2,:], mc=colors[i], ms=ms, label=false)
        annotate!(33, tip_twist[i,end-1], text(labels[i], :top, :right, colors[i]))
    end
end
display(plt2)
savefig(string(absPath,"/TDWingAirspeedRange_twist.pdf"))
# Aeroelastic frequencies at root angle of 1 deg
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]",legend=:outertop)
plot!([NaN], [NaN], lc=:black,  lw=lw, ls=:solid, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Arena & Lacarbonara (2013) - Num.")
modeLabels = ["Flapwise bending" "Chordwise bending" "Torsion"]
for (m,mode) in enumerate([1,2,4])
    freqsMode = [freqs[1,j][mode] for j in 1:size(freqs, 2)]
    plot!(URange,freqsMode, c=colors2[m], lw=2, ls=:solid, label=false)
    scatter!(freqs_ref[1,:],freqs_ref[1+m,:], mc=colors2[m], ms=ms, label=false)
    annotate!(10, freqsMode[1], text(modeLabels[m], :bottom, colors2[m]))
end
display(plt3)
savefig(string(absPath,"/TDWingAirspeedRange_freqs.pdf"))

println("Finished TDWingAirspeedRange.jl")