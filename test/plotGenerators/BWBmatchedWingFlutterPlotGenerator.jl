using Plots, ColorSchemes

# Run the script
include("../examples/BWBmatchedWingFlutter.jl")

# Set paths
relPath = "/test/outputs/figures/BWBmatchedWingFlutter"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Each mode shape animation at lowest flutter speed
for mode in 1:nModes
    modeAnim = plot_mode_shapes_animation(eigenProblem[flutterSpeedInd],modes2plot=[mode],nFramesPerCycle=21,fps=20,scale=2,view=(45,30),element2centralize=11,legendPos=:top,plotBCs=false,plotAxes=false,displayProgress=true,save=true,savePath=string(relPath,"/BWBmatchedWingFlutter_mode",mode,".gif"))
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
gr()

# Trim root angle of attack vs airspeed
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=extrema(URange), tickfont=font(ts), guidefont=font(fs))
plot!(URange, trimAoA*180/π, c=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/BWBmatchedWingFlutter_AoA.pdf"))

# Trim thrust force vs airspeed
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=extrema(URange), tickfont=font(ts), guidefont=font(fs))
plot!(URange, trimThrust, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/BWBmatchedWingFlutter_thrust.pdf"))

# Trim elevator deflection vs airspeed
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=extrema(URange), tickfont=font(ts), guidefont=font(fs))
plot!(URange, trimδ*180/π, c=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/BWBmatchedWingFlutter_delta.pdf"))

# Root locus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-30,5],ylims=[0,125], tickfont=font(ts), guidefont=font(fs))
for mode in 1:nModes
    scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], ms=ms, msw=msw, label=false)
end
display(plt_RL)
savefig(string(absPath,"/BWBmatchedWingFlutter_rootlocus.pdf"))

# V-g-f
plt51 = plot(ylabel="Frequency [rad/s]", xlims=extrema(URange), ylims=[0,100], tickfont=font(ts), guidefont=font(12))
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt52 = plot(xlabel="Airspeed [m/s]", ylabel="Damping ratio", xlims=extrema(URange), ylims=[-0.3,0.2], tickfont=font(ts), guidefont=font(12), yticks=-0.3:0.1:0.2)
for mode in 1:nModes
    plot!(URange, modeDampings[mode]./modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=ms, msw=msw, label=false)
end
plt5 = plot(plt51,plt52, layout=(2,1))
display(plt5)
savefig(string(absPath,"/BWBmatchedWingFlutter_Vgf.pdf"))

println("Finished BWBmatchedWingFlutterPlotGenerator.jl")