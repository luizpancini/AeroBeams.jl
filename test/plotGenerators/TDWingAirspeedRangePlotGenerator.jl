using Plots, ColorSchemes

# Run the script
include("../examples/TDWingAirspeedRange.jl")

# Set paths
relPath = "/test/outputs/figures/TDWingAirspeedRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed shape
deformationPlot = plot_steady_deformation(problem,save=true,savePath=string(relPath,"/TDWingAirspeedRange_deformation.pdf"))
display(deformationPlot)

# Plot configurations
colors = get(colorschemes[:darkrainbow], LinRange(0, 1, length(θRange)))
colors2 = get(colorschemes[:darkrainbow], LinRange(0, 1, 3))
lw = 2
ms = 4
msw = 0
labels = ["\\theta_{r} = 1.0 deg" "\\theta_{r} = 2.2 deg"]
gr()

# Tip flapwise displacement
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Tip flapwise displacement [m]")
plot!([NaN], [NaN], lc=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN], [NaN], lc=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
for i=eachindex(θRange)
    plot!(URange,tip_u3[i,:], c=colors[i], lw=lw, ls=:solid, label=false)
    if i==1
        plot!(u3_1deg_num[1,:],u3_1deg_num[2,:], c=colors[i], lw=lw, ls=:dash, label=false)
        scatter!(u3_1deg_exp[1,:],u3_1deg_exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        annotate!(30, tip_u3[i,2], text(labels[i], :bottom, :left, colors[i]))
    else
        plot!(u3_2_2deg_num[1,:],u3_2_2deg_num[2,:], c=colors[i], lw=lw, ls=:dash, label=false)
        scatter!(u3_2_2deg_exp[1,:],u3_2_2deg_exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        annotate!(35, tip_u3[i,end-1], text(labels[i], :top, :right, colors[i]))
    end
end
display(plt1)
savefig(string(absPath,"/TDWingAirspeedRange_disp.pdf"))

# Tip twist
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]")
plot!([NaN], [NaN], lc=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN], [NaN], lc=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!([NaN], [NaN ], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
for i=eachindex(θRange)
    plot!(URange,tip_twist[i,:], c=colors[i], lw=lw, ls=:solid, label=false)
    if i==1
        plot!(th_1deg_num[1,:],th_1deg_num[2,:], c=colors[i], lw=lw, ls=:dash, label=false)
        scatter!(th_1deg_exp[1,:],th_1deg_exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        annotate!(30, tip_twist[i,2], text(labels[i], :bottom, :left, colors[i]))
    else
        plot!(th_2_2deg_num[1,:],th_2_2deg_num[2,:], c=colors[i], lw=lw, ls=:dash, label=false)
        scatter!(th_2_2deg_exp[1,:],th_2_2deg_exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        annotate!(33, tip_twist[i,end-1], text(labels[i], :top, :right, colors[i]))
    end
end
display(plt2)
savefig(string(absPath,"/TDWingAirspeedRange_twist.pdf"))

# Aeroelastic frequencies at root angle of 1 deg
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]",legend=:outertop)
plot!([NaN], [NaN], lc=:black, lw=lw, ls=:solid, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=msw, label="Arena & Lacarbonara (2013) - Num.")
modeLabels = ["Flapwise bending" "Chordwise bending" "Torsion"]
for (m,mode) in enumerate([1,2,4])
    freqsMode = [freqs[1,j][mode] for j in 1:size(freqs, 2)]
    plot!(URange,freqsMode, c=colors2[m], lw=lw, ls=:solid, label=false)
    scatter!(freqs_ref[1,:],freqs_ref[1+m,:], mc=colors2[m], ms=ms, msw=msw, label=false)
    annotate!(10, freqsMode[1], text(modeLabels[m], :bottom, colors2[m]))
end
display(plt3)
savefig(string(absPath,"/TDWingAirspeedRange_freqs.pdf"))

println("Finished TDWingAirspeedRangePlotGenerator.jl")