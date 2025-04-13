using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazySteadyPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazySteadyPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed state of selected condition
deformationPlot = plot_steady_deformation(problem[3,60], interactive=true,view=(30,15),plotDistLoads=false,save=true,savePath="/test/outputs/figures/sweptPazySteadyPitchRange/sweptPazySteadyPitchRange_deformation.pdf")
display(deformationPlot)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 3
msw = 0
gr()

# Deformed position at θ=5 deg, U=60 m/s
x = upright ? x3_def[3,60] : x1_def[3,60]
y = upright ? -x1_def[3,60] : x3_def[3,60]
plt_defPos = plot(xlabel="\$x_1\$ [m]", ylabel="\$x_3\$ [m]", xlims=[0,0.5], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!(x, y, c=:black, lw=lw, ls=:solid, label="AeroBeams")
scatter!(dispRef[1,:],dispRef[2,:], c=:black, ms=ms, msw=msw, label="Technion (Nastran)")
display(plt_defPos)
savefig(string(absPath,"/sweptPazySteadyPitchRange_defPos.pdf"))

# Tip OOP displacement vs. airspeed
plt_tipOOP = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,100], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
for (i,θ) in enumerate(θRange)
    plot!(URange, tipOOP[i,:]/L*100, c=colors[i], lw=lw, ls=:solid, label="\$\\theta = $(round(Int,θ*180/pi)) ^\\circ\$")
end
display(plt_tipOOP)
savefig(string(absPath,"/sweptPazySteadyPitchRange_tipOOP.pdf"))

# Tip twist vs. airspeed
plt_tipTwist = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,100], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tipTwist[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_tipTwist)
savefig(string(absPath,"/sweptPazySteadyPitchRange_tipTwist.pdf"))

# Tip AoA vs. airspeed
plt_tipAOA = plot(xlabel="Airspeed [m/s]", ylabel="Tip angle of attack [deg]", xlims=[0,100], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tipAoA[i,:], c=colors[i], lw=lw, label=false)
end
display(plt_tipAOA)
savefig(string(absPath,"/sweptPazySteadyPitchRange_tipAoA.pdf"))

# Normal force coefficient slope over span at selected condition
plt_cna = plot(xlabel="Normalized span", ylabel="\$c_{n_\\alpha}\$ [1/rad]", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
plot!(x1_e/L, cn[3,60]./AoA[3,60], c=:black, lw=lw, label=false)
display(plt_cna)
savefig(string(absPath,"/sweptPazySteadyPitchRange_cna.pdf"))

println("Finished sweptPazySteadyPitchRangePlotGenerator.jl")