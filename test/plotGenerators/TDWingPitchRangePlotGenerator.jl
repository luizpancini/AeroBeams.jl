using Plots, ColorSchemes

# Run the script
include("../examples/TDWingPitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/TDWingPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:darkrainbow], LinRange(0, 1, 3))
lw = 2
ms = 5
msw = 0

# Tip flapwise displacement
plt1 = plot(xlabel="Root angle [deg]", ylabel="Tip flapwise displacement [m]", xticks=collect(0:15:90))
plot!(θRange*180/π,-tip_u3, c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!(u3_num[1,:],u3_num[2,:], c=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!(u3_exp[1,:],u3_exp[2,:], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
display(plt1)
savefig(string(absPath,"/TDWingPitchRange_u3.pdf"))

# Tip chordwise displacement
plt2 = plot(xlabel="Root angle [deg]", ylabel="Tip chordwise displacement [m]", xticks=collect(0:15:90))
plot!(θRange*180/π,-tip_u2, c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!(u2_num[1,:],u2_num[2,:], c=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!(u2_exp[1,:],u2_exp[2,:], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
display(plt2)
savefig(string(absPath,"/TDWingPitchRange_u2.pdf"))

# Tip twist
plt3 = plot(xlabel="Root angle [deg]", ylabel="Tip twist [deg]", xticks=collect(0:15:90), legend=:bottom)
plot!(θRange*180/π,-tip_twist, c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!(th_num[1,:],th_num[2,:], c=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!(th_exp[1,:],th_exp[2,:], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
display(plt3)
savefig(string(absPath,"/TDWingPitchRange_twist.pdf"))

# Structural frequencies 
plt4 = plot(xlabel="Root angle [deg]", ylabel="Frequency [Hz]", legend=:outertop, xticks=collect(0:15:90))
modeLabels = ["Flapwise bending" "Chordwise bending" "Torsion"]
plot!([NaN],[NaN], c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN],[NaN], c=:black, lw=lw, ls=:dash, label="Tang & Dowell (2001) - Num.")
scatter!([NaN],[NaN], mc=:black, ms=ms, msw=msw, label="Tang & Dowell (2001) - Exp.")
for (m,mode) in enumerate([1,2,4])
    freqsMode = [freqs[j][mode] for j in 1:length(freqs)]
    plot!(θRange*180/π,freqsMode, c=colors[m], lw=lw, ls=:solid, label=false)
    plot!(freqs_num[m][1,:],freqs_num[m][2,:], c=colors[m], lw=lw, ls=:dash, label=false)
    scatter!(freqs_exp[1,:],freqs_exp[1+m,:], mc=colors[m], ms=ms, msw=msw, label=false)
    annotate!(30, freqsMode[1], text(modeLabels[m], :bottom, colors[m]))
end
display(plt4)
savefig(string(absPath,"/TDWingPitchRange_freqs.pdf"))

println("Finished TDWingPitchRangePlotGenerator.jl")