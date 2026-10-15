using AeroBeams, JLD2

# Bending curvature and airspeed ranges
k2Range = range(-0.015,0.045,5)
URange = vcat(1e-3,0.5:0.5:50)

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1

# Altitude [m]
h = 20e3

# No gravity
g = 0

# Root pitch angle [rad]
θ = deg2rad(0)

# Option to include induced drag
hasInducedDrag = true

# Parasite drag
cd0 = 1e-2

# Discretization
nElem = 40

# System solver
systemSolver = create_NLsolver()

# Number of modes to track
nModes = 5

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 80

# Initialize outputs
eigenProblem = Array{EigenProblem}(undef,length(k2Range),length(URange))
eigensols = Array{EigenSolution}(undef,length(k2Range),length(URange))
freqsOsc = Array{Matrix{Float64}}(undef, length(k2Range))
dampsOsc = Array{Matrix{Float64}}(undef, length(k2Range))
modeShapesOsc = Array{Matrix{Union{ModeShape, Nothing}}}(undef, length(k2Range))
flutterOnsetSpeeds = Array{Vector{Float64}}(undef, length(k2Range))
flutterOnsetFreqs = Array{Vector{Float64}}(undef, length(k2Range))
flutterOnsetModes = Array{Vector{Int}}(undef, length(k2Range))
flutterOnsetTipDisp = Array{Vector{Float64}}(undef, length(k2Range))
flutterOffsetSpeeds = Array{Vector{Float64}}(undef, length(k2Range))
flutterOffsetFreqs = Array{Vector{Float64}}(undef, length(k2Range))
flutterOffsetModes = Array{Vector{Int}}(undef, length(k2Range))
flutterOffsetTipDisp = Array{Vector{Float64}}(undef, length(k2Range))
x1_0 = Array{Vector{Float64}}(undef, length(k2Range))
x3_0 = Array{Vector{Float64}}(undef, length(k2Range))
x1_n = Array{Vector{Float64}}(undef, length(k2Range))
x1_e = Array{Vector{Float64}}(undef, length(k2Range))
x1_def = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
x3_def = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
twist = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
tip_u3 = Array{Float64}(undef, length(k2Range),length(URange))
α = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
cl = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
cdrag = Array{Vector{Float64}}(undef, length(k2Range),length(URange))

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model
        model,L = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElem=nElem,altitude=h,g=g,pitchAngle=θ,k2=k2,cd0=cd0,hasInducedDrag=hasInducedDrag)
        # Set initial guess solution as previous known solution
        x0 = j == 1 ? zeros(0) : eigenProblem[i,j-1].systemArrays.x
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=model,systemSolver=systemSolver,eigenmethod=eigenmethod,nev=nev,x0=x0)
        solve!(eigenProblem[i,j])
        # Eigensolution
        eigensols[i,j] = eigenProblem[i,j].eigensolution
        # Final states
        states = eigenProblem[i,j].loadstepHistory.states[end]
        aeroVars = eigenProblem[i,j].loadstepHistory.aeroVariables[end]
        # Undeformed jig-shape properties
        if j == 1
            # Undeformed nodal positions
            x1_0[i] = reduce(vcat, 
                        [model.assemblyData.elements[e].assemblyData.r_n1[1], 
                         model.assemblyData.elements[e].assemblyData.r_n2[1]] 
                    for e in 1:nElem)
            x3_0[i] = reduce(vcat, 
                        [model.assemblyData.elements[e].assemblyData.r_n1[3], 
                         model.assemblyData.elements[e].assemblyData.r_n2[3]] 
                    for e in 1:nElem)
            # Nodal and elemental arclength positions
            x1_n[i] = vcat([vcat(model.assemblyData.elements[e].geometry.x1_n1,model.assemblyData.elements[e].geometry.x1_n2) for e in 1:nElem]...)
            x1_e[i] = [model.assemblyData.elements[e].geometry.x1 for e in 1:nElem]
        end
        # Displacements over span
        u1_of_x1 = reduce(vcat,[states[e].nodal.u1[1],
                                states[e].nodal.u2[1]] for e in 1:nElem)
        u3_of_x1 = reduce(vcat,[states[e].nodal.u1[3],
                                states[e].nodal.u2[3]] for e in 1:nElem)
        u1_of_x1 .-= u1_of_x1[1]
        u3_of_x1 .-= u3_of_x1[1]
        # Deformed nodal positions
        x1_def[i,j] = x1_0[i] + u1_of_x1
        x3_def[i,j] = x3_0[i] + u3_of_x1
        # Angle of twist and tip OOP displacement
        twist[i,j] = vcat([vcat(twist_angle(states[e].nodal.p1_b), twist_angle(states[e].nodal.p2_b)) for e in 1:nElem]...)
        tip_u3[i,j] = states[nElem].nodal.u2[3]
        # Angle of attack and force coefficients over span
        α[i,j] = [rad2deg(aeroVars[e].unsteadyKin.αₑ) for e in 1:nElem]
        cn = [aeroVars[e].aeroCoefficients.cn for e in 1:nElem]
        ct = [aeroVars[e].aeroCoefficients.ct for e in 1:nElem]
        cl[i,j] = @. cn*cosd(α[i,j]) + ct*sind(α[i,j])
        cdrag[i,j] = @. cn*sind(α[i,j]) - ct*cosd(α[i,j])
    end
    # Apply mode tracking
    freqs,damps,modeShapes,_ = mode_tracking(eigensols[i,:])
    # Frequencies, dampings and shapes of oscillatory modes
    freqsOsc[i],dampsOsc[i],modeShapesOsc[i],_ = filter_tracked_modes(freqs,damps,modeShapes; fmin=1e-2, maxTracks=nModes)
    # Flutter variables
    flutterOnsetSpeeds[i],flutterOnsetFreqs[i],flutterOnsetModes[i],flutterOnsetTipDisp[i],flutterOffsetSpeeds[i],flutterOffsetFreqs[i],flutterOffsetModes[i],flutterOffsetTipDisp[i] = find_flutter(URange,freqsOsc[i],dampsOsc[i],tip_u3[i,:])
end

# Lowest flutter onset and offset indices
flutterOnsetInds = [searchsortedfirst.(Ref(URange), U) for U in flutterOnsetSpeeds]
flutterOffsetInds = [searchsortedfirst.(Ref(URange), U) for U in flutterOffsetSpeeds]

## Plots and saves
# Save flutter data
flutterOnsetSpeedWingg0 = flutterOnsetSpeeds
flutterOnsetFreqWingg0 = flutterOnsetFreqs
flutterOnsetTipOOPWingg0 = flutterOnsetTipDisp
@save absPathData*string("wing_g0_lambda",λ,"_flutterSpeed.jld2") flutterOnsetSpeedWingg0
@save absPathData*string("wing_g0_lambda",λ,"_flutterFreq.jld2") flutterOnsetFreqWingg0
@save absPathData*string("wing_g0_lambda",λ,"_flutterOOP.jld2") flutterOnsetTipOOPWingg0

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_g0_flutter_fixedPitch_k2_range"
absPath = string(pwd(),relPath)
mkpath(absPath)
relPathData = "/dev/cHALE/Flexible/outputs/data/cHALEwing_g0_flutter_fixedPitch_k2_range/"
absPathData = string(pwd(),relPathData)
mkpath(absPathData)

# Plot configurations
colors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold])
ts = 10
fs = 16
lfs = 10
lw = 2
ms = 3
msw = 0
tsz = 9
DPI = 300
fps = 10
mshape = [:circle, :star, :utriangle, :pentagon, :diamond]
labels = ["\$k_2 = $(k2) \$" for k2 in k2Range]
L = 16
gr()

# Undeformed jig-shapes
plt_shapes = plot(xlabel="Normalized spanwise direction", ylabel="Normalized vertical direction", xlims=[0,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:bottomleft, dpi=DPI)
for (i,k2) in enumerate(k2Range)
    plot!(x1_0[i]/L, x3_0[i]/L, c=colors[i], lw=lw, label=labels[i])
end
display(plt_shapes)
savefig(string(absPath,"/cHALEwing_undef_k2range.pdf"))

# Root locus
freqLim = [0,50]
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,2], ylims=freqLim, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([0; 0], freqLim, c=:black, ls=:dash, lw=2, label=false)
if λ == 1
    OOPtext = "OOP"
    TIP1text = "1st T-IP"
    TIP2text = "2nd T-IP"
    Ttext = "T"
    IPtext = "IP"
    OOPtextPos = [-8, 20]
    TIP1textPos = [1, 5]
    TIP2textPos = [0.75, 45]
    TtextPos = [0.25, 30]
    IPtextPos = [0.25, 32]
    annotate!(OOPtextPos[1], OOPtextPos[2], text(OOPtext, tsz))
    annotate!(TIP1textPos[1], TIP1textPos[2], text(TIP1text, tsz))
    annotate!(TIP2textPos[1], TIP2textPos[2], text(TIP2text, tsz))
    annotate!(TtextPos[1], TtextPos[2], text(Ttext, tsz))
    annotate!(IPtextPos[1], IPtextPos[2], text(IPtext, tsz))
    quiver!([OOPtextPos[1],OOPtextPos[1]+0.5,OOPtextPos[1]+0.5], [OOPtextPos[2]-2,OOPtextPos[2]-0.5,OOPtextPos[2]+1], quiver=([0,1.75,3.75], [-10,-3.5,17]), arrow=:closed, linecolor=:black)
end
plt_RL_base = deepcopy(plt_RL)
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label="\$U_{\\infty} = $(round(Int, URange[1]))\\;\\mathrm{m/s}\$")
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        plot!(dampsOsc[i][mode,:], freqsOsc[i][mode,:], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        plot!([dampsOsc[i][mode,1]], [freqsOsc[i][mode,1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
display(plt_RL)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_rootlocus_lambda",λ,".pdf")))

# Animated root locus plot
plt_RL_anim = plot(plt_RL_base)
plot!(dpi=DPI)
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
end
anim = @animate for (j,U) in enumerate(URange)
    title!("\$U_{\\infty} = $(round(Int, U))\\;\\mathrm{m/s}\$")
    for (i,k2) in enumerate(k2Range)
        for mode in 1:nModes
            plot!([dampsOsc[i][mode,j]], [freqsOsc[i][mode,j]], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        end
    end
end
gif_handle = gif(anim, string(absPath,"/cHALEwing_g0_flutter_fixedPitch_k2_range_RL.gif"), fps=fps)
display(gif_handle)

# V-g-f
Vgf_k2_ind = [2,5]
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,50], tickfont=font(ts), guidefont=font(12))
for (i,k2) in enumerate(k2Range)
    if !(i in Vgf_k2_ind)
        continue
    end
    for mode in 1:nModes
        plot!(URange, freqsOsc[i][mode,:], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.2,0.2], tickfont=font(ts), guidefont=font(12), legendfontsize=10, legend=:topleft)
plot!(plt_Vg,URange,zeros(length(URange)), c=:gray, lw=lw, ls=:dash, label=false)
for (i,k2) in enumerate(k2Range)
    if !(i in Vgf_k2_ind)
        continue
    end
    for mode in 1:nModes
        labelNow = mode == 1 ? labels[i] : false
        plot!(URange, dampsOsc[i][mode,:]./freqsOsc[i][mode,:], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=labelNow)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_Vgf_lambda",λ,".pdf")))

# Normalized deformed span at iminence of flutter onset
plt_disp = plot(xlabel="Normalized spanwise direction", ylabel="Normalized vertical direction", xlims=[0,1], ylims=[-1,1], yticks=-1:0.2:1, tickfont=font(ts), guidefont=font(fs), legendfontsize=8, legend_position=:topright)
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dot, label="Undeformed")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dash, label="At 1st flutter condition")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:solid, label="At 2nd flutter condition")
for (i,k2) in enumerate(k2Range)
    plot!(x1_0[i]/L, x3_0[i]/L, c=colors[i], lw=lw, ls=:dot, label=false)
    if !isempty(flutterOnsetInds[i])
        plot!(x1_def[i,flutterOnsetInds[i][1]]/L, x3_def[i,flutterOnsetInds[i][1]]/L, c=colors[i], lw=lw, ls=:dash, label=false)
        if length(flutterOnsetInds[i]) > 1
            plot!(x1_def[i,flutterOnsetInds[i][2]]/L, x3_def[i,flutterOnsetInds[i][2]]/L, c=colors[i], lw=lw, ls=:solid, label=labels[i])
        end
    end
end
display(plt_disp)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_disp_lambda",λ,".pdf")))

# AoA at iminence of flutter onset
plt_AoA = plot(xlabel="Normalized arclength", ylabel="Angle of attack [deg]", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    if !isempty(flutterOnsetInds[i])
        plot!(x1_e[i]/L, α[i,flutterOnsetInds[i][1]], c=colors[i], lw=lw, ls=:dash, label=false)
        if length(flutterOnsetInds[i]) > 1
            plot!(x1_e[i]/L, α[i,flutterOnsetInds[i][2]], c=colors[i], lw=lw, ls=:solid, label=false)
        end
    end
end
display(plt_AoA)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_AoA_lambda",λ,".pdf")))

# Angle of twist at iminence of flutter onset
plt_twist = plot(xlabel="Normalized arclength", ylabel="Angle of twist [deg]", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    if !isempty(flutterOnsetInds[i])
        plot!(x1_n[i]/L, twist[i,flutterOnsetInds[i][1]], c=colors[i], lw=lw, ls=:dash, label=false)
        if length(flutterOnsetInds[i]) > 1
            plot!(x1_n[i]/L, twist[i,flutterOnsetInds[i][2]], c=colors[i], lw=lw, ls=:solid, label=false)
        end
    end
end
display(plt_twist)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_twist_lambda",λ,".pdf")))

# Flutter onset speeds vs k2
plt_Uf = plot(xlabel="\$k_2\$ [1/m]", ylabel="Flutter speed [m/s]", xlims=[-0.02,0.05], ylims=[0,35], xticks=k2Range, tickfont=font(ts), guidefont=font(fs))
plot!(k2Range,[flutterOnsetSpeeds[i][1] for i in eachindex(k2Range)], marker=(:circle, 10), msw=msw, c=:black, lw=lw, label=false)
display(plt_Uf)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_speedOn.pdf")))

# Flutter onset speeds vs tip OOP position
plt_Uf = plot(xlabel="Tip OOP position [% semispan]", ylabel="Flutter speed [m/s]", xlims=[-40,20], ylims=[0,40], xticks=vcat(-40:10:20), tickfont=font(ts), guidefont=font(fs))
plot!([flutterOnsetTipDisp[i][1] for i in eachindex(k2Range)]*100/L,[flutterOnsetSpeeds[i][1] for i in eachindex(k2Range)], marker=(:circle, 10), msw=msw, c=:black, lw=lw, label=false)
display(plt_Uf)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_speedOn_OOP_lambda",λ,".pdf")))

##
# Mode shapes of straight wing at lowest airspeed
modesPlot = plot_mode_shapes(eigenProblem[2,1],modeShapes=modeShapesOsc[2][:,1],scale=5,view=(60,15),legendPos=(0,0.99),nullDampingTol=1e-4,modalColorScheme=colors,modeLabels=["1st OOP","2nd OOP","T","IP","3rd OOP"],save=true,savePath=string(relPath,"/cHALEwing_modeShapes_k0U1_lambda",λ,".pdf"))
display(modesPlot)

# Mode shapes of straight wing at highest airspeed
modesPlot = plot_mode_shapes(eigenProblem[2,end],modeShapes=modeShapesOsc[2][:,end],scale=5,view=(60,15),legendPos=(0,0.99),nullDampingTol=1e-4,modalColorScheme=colors,modeLabels=["1st OOP","2nd OOP","T","IP","3rd OOP"],save=true,savePath=string(relPath,"/cHALEwing_modeShapes_k0U35_lambda",λ,".pdf"))
display(modesPlot)

# Mode shapes of wing with k_2=0.045 at lowest airspeed
modesPlot = plot_mode_shapes(eigenProblem[5,1],modeShapes=modeShapesOsc[5][:,1],scale=5,view=(60,15),legendPos=(0,0.99),nullDampingTol=1e-4,modalColorScheme=colors,modeLabels=["1st OOP","1st T-IP","2nd OOP","3rd OOP","2nd T-IP"],save=true,savePath=string(relPath,"/cHALEwing_modeShapes_k0045U1_lambda",λ,".pdf"))
display(modesPlot)

# Mode shapes of wing with k_2=0.045 at highest airspeed
modesPlot = plot_mode_shapes(eigenProblem[5,end],modeShapes=modeShapesOsc[5][:,end],scale=5,view=(60,15),legendPos=(0,0.99),nullDampingTol=1e-4,modalColorScheme=colors,modeLabels=["1st OOP","1st T-IP","2nd OOP","2nd T-IP","3rd OOP"],save=true,savePath=string(relPath,"/cHALEwing_modeShapes_k0045U35_lambda",λ,".pdf"))
display(modesPlot)

# Mode shapes animation of wing with k_2=0.0 at lowest airspeed
modesAnim = plot_mode_shapes_animation(eigenProblem[2,1],modeShapes=modeShapesOsc[2][3:4,1],scale=10,view=(60,15),legendPos=(0.4,0.9),modeLabels=["T","IP"],nFramesPerCycle=21,plotSteady=true,plotBCs=true,plotAxes=false,modalColorScheme=colors,legendFontSize=10,nullDampingTol=1e-4,plotLimits=([0,L],[-L/2,L/2],[-L/2,L/2]),save=true,savePath=string(relPath,"/cHALEwing_modeShapes_k0U1_lambda",λ,".gif"),displayProgress=true)
display(modesAnim)

# Mode shapes animation of wing with k_2=0.045 at lowest airspeed
modesAnim = plot_mode_shapes_animation(eigenProblem[end,1],modeShapes=modeShapesOsc[5][[2,5],1],scale=10,view=(60,15),legendPos=(0.4,0.9),modeLabels=["T-IP-1","T-IP-2"],nFramesPerCycle=101,plotSteady=true,plotBCs=true,plotAxes=false,modalColorScheme=colors,legendFontSize=10,nullDampingTol=1e-4,save=true,savePath=string(relPath,"/cHALEwing_modeShapes_k0045U1_lambda",λ,".gif"),displayProgress=true)
display(modesAnim)

# Steady plot of unloaded straight wing
plt_steady = plot_steady_deformation(eigenProblem[2,1],plotUndeformed=false,plotBCs=true,plotDistLoads=false,plotLegend=false,plotAxes=false,plotGrid=false,view=(30,30),plotLimits=([-L*1/4,L*3/4],[-L*1/2,L*1/2],[-L*1/2,L*1/2]),save=true,savePath=string(relPath,"/cHALEwing_steady_k0U0_lambda",λ,".pdf"))
display(plt_steady)

println("Finished cHALEwing_g0_flutter_fixedPitch_k2_range.jl")