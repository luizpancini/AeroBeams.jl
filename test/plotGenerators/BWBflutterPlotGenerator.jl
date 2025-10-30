using Plots, ColorSchemes

# Run the script
include("../examples/BWBflutter.jl")

# Set paths
relPath = "/test/outputs/figures/BWBflutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at highest airspeed
modesPlot = plot_mode_shapes(eigenProblem[end],scale=2,view=(45,30),element2centralize=11,legendPos=:outertop,modalColorScheme=:rainbow,save=true,savePath=string(relPath,"/BWBflutter_modeShapes.pdf"))
display(modesPlot)

# Flutter mode shape animation at flutter speed
flutterModeAnim = plot_mode_shapes_animation(eigenProblem[flutterSpeedInd],modes2plot=[4],nFramesPerCycle=21,fps=20,scale=2,view=(45,30),element2centralize=11,legendPos=:top,plotBCs=false,plotAxes=false,displayProgress=true,save=true,savePath=string(relPath,"/BWBflutter_flutterMode.gif"))
display(flutterModeAnim)

# Each mode shape animation at lowest flutter speed
for mode in 1:nModes
    modeAnim = plot_mode_shapes_animation(eigenProblem[flutterSpeedInd],modes2plot=[mode],nFramesPerCycle=21,fps=20,scale=2,view=(45,30),element2centralize=11,legendPos=:top,plotBCs=false,plotAxes=false,displayProgress=true,save=true,savePath=string(relPath,"/BWBflutter_mode",mode,".gif"))
    display(modeAnim)
end

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0
fps = 5
DPI = 300
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[URange[1],URange[end]+1], tickfont=font(ts), guidefont=font(fs), legend=:topright, legendfontsize=lfs, xticks=30:10:160)
plot!(URange, trimAoA*180/π, c=:black, lw=lw, label="AeroBeams")
scatter!(trimAoARef[1,:],trimAoARef[2,:], c=:black, ms=ms, label="UM/NAST")
display(plt1)
savefig(string(absPath,"/BWBflutter_AoA.pdf"))

# Trim thrust force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[URange[1],URange[end]+1], tickfont=font(ts), guidefont=font(fs), xticks=30:10:160)
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
scatter!(trimThrustRef[1,:],trimThrustRef[2,:], c=:black, ms=ms, label=false)
display(plt2)
savefig(string(absPath,"/BWBflutter_thrust.pdf"))

# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[URange[1],URange[end]+1], tickfont=font(ts), guidefont=font(fs), xticks=30:10:160)
plot!(URange, trimδ*180/π, c=:black, lw=lw, label=false)
scatter!(trimδRef[1,:],trimδRef[2,:], c=:black, ms=ms, label=false)
display(plt3)
savefig(string(absPath,"/BWBflutter_delta.pdf"))

# Root locus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-35,5],ylims=[0,125], tickfont=font(ts), guidefont=font(fs), legend_position=(0.2,0.5), legendfontsize=10)
scatter!([NaN],[NaN], c=:black, shape=:circle, ms=ms, msw=msw, label="AeroBeams")
scatter!([NaN],[NaN], c=:black, shape=:utriangle, ms=ms, msw=msw, label="UM/NAST")
plt_RL_base = deepcopy(plt_RL)
for mode in 1:nModes
    scatter!(dampsRef[mode+1,:], 2π*freqsRef[mode+1,:], c=:black, shape=:utriangle, ms=ms, msw=msw, label=false)
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
end
display(plt_RL)
savefig(string(absPath,"/BWBflutter_rootlocus.pdf"))

# Animated root locus plot
plt_RL_anim = plot(plt_RL_base)
plot!(dpi=DPI)
anim = @animate for (j,U) in enumerate(URange)
    title!("\$U_{\\infty} = $U\$ m/s")
    for mode in 1:nModes
        scatter!([dampsRef[mode+1,j]], [2π*freqsRef[mode+1,j]], c=:black, shape=:utriangle, ms=ms, msw=msw, label=false)
        scatter!([modeDampings[mode][j]], [modeFrequencies[mode][j]], c=modeColors[mode], ms=ms, msw=msw, label=false)
    end
end
gif_handle = gif(anim, string(absPath,"/BWBflutter_rootlocus.gif"), fps=fps)
display(gif_handle)

# V-g-f
plt51 = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]+1], ylims=[0,125], tickfont=font(ts), guidefont=font(12), xticks=30:10:160)
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt52 = plot(xlabel="Airspeed [m/s]", ylabel="Damping ratio", xlims=[URange[1],URange[end]+1], ylims=[-0.3,0.2], tickfont=font(ts), guidefont=font(12), xticks=30:10:160, yticks=-0.3:0.1:0.2)
for mode in 1:nModes
    plot!(URange, modeDampings[mode]./modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt5 = plot(plt51,plt52, layout=(2,1))
display(plt5)
savefig(string(absPath,"/BWBflutter_Vgf.pdf"))

println("Finished BWBflutterPlotGenerator.jl")