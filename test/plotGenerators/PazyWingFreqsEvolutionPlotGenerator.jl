using Plots, ColorSchemes

# Run the script
include("../examples/PazyWingFreqsEvolution.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingFreqsEvolution"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
modeColors = cgrad(:rainbow, nModes, categorical=true)
ts = 10
fs = 16
lfs = 10
lw = 2
ms = 3
msw = 0
gr()

# Frequencies vs. airspeed
pltU = plot(xlabel="Airspeed [m/s]", ylabel="Frequency [Hz]", title="\$\\alpha_r = 7^\\circ\$", xlims=[0,45], ylims=[0,120], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dash, label="Riso & Cesnik (2023)")
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode]/(2π), c=modeColors[mode], lw=lw, label=false)
    plot!(freqsVsU_alpha7[1,:], freqsVsU_alpha7[2,:], lw=lw, ls=:dash, c=:black, label=false)
end
display(pltU)
savefig(string(absPath,"/PazyWingFreqsEvolution_U.pdf"))

# Frequencies vs. tip OOP displacement 
pltDisp = plot(xlabel="Tip displacement [% semispan]", ylabel="Frequency [Hz]", title="\$\\alpha_r = 7^\\circ\$", xlims=[0,45], ylims=[0,120], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dash, label="Riso & Cesnik (2023)")
for mode in 1:nModes
    plot!(tipOOP*100/L, modeFrequencies[mode]/(2π), c=modeColors[mode], lw=lw, label=false)
    plot!(freqsVsDisp_alpha7[1,:], freqsVsDisp_alpha7[2,:], lw=lw, ls=:dash, c=:black, label=false)
end
display(pltDisp)
savefig(string(absPath,"/PazyWingFreqsEvolution_disp.pdf"))

println("Finished PazyWingFreqsEvolutionPlotGenerator.jl")