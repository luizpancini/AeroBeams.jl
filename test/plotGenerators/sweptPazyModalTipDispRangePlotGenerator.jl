using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazyModalTipDispRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazyModalTipDispRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, nModes, categorical=true)
lfs = 12
lw = 2
ls = [:solid :dash :dot]
gr()

# Frequencies vs sweep angle
plt1 = plot(xlabel="Normalized tip displacement", ylabel="Frequency [Hz]", xlims=[0,0.5], ylims=[0,50], legendfontsize=lfs, colorbar_title="Mode")
for (i,Λ) in enumerate(ΛRange)
    plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[i], label="\$\\Lambda=\$$(round.(Int,Λ*180/π))\$^\\circ\$")
end
for (i,Λ) in enumerate(ΛRange)
    for mode in 1:nModes
        plot!(tipOOP[i,:]/L, modeFrequencies[i,mode], c=colors, lw=lw, ls=ls[i], lz=mode, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/sweptPazyModalTipDispRange.pdf"))

# Tip displacement vs. sweep angle for largest load
tip_u3_maxload = tip_u3[:,end]
plt_u3 = plot(xlabel="Sweep angle [deg]", ylabel="Tip \$u_3/L\$")
for i in eachindex(ΛRange)
    plot!(ΛRange*180/π, tip_u3_maxload/L, c=:black, lw=lw, label=false)
end
display(plt_u3)
savefig(string(absPath,"/sweptPazyModalTipDispRangeDisp.pdf"))

# Tip twist vs. sweep angle for largest load
tip_twist_maxload = tip_twist[:,end]
plt_twist = plot(xlabel="Sweep angle [deg]", ylabel="Tip twist [deg]")
for i in eachindex(ΛRange)
    plot!(ΛRange*180/π, tip_twist_maxload, c=:black, lw=lw, label=false)
end
display(plt_twist)
savefig(string(absPath,"/sweptPazyModalTipDispRangeTwist.pdf"))

println("Finished sweptPazyModalTipDispRangePlotGenerator.jl")