using Plots, ColorSchemes

# Run the script
include("../examples/SMWFlutterStructuralCouplingRange.jl")

# Set paths
relPath = "/test/outputs/figures/SMWFlutterStructuralCouplingRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(ΨRange)))
lw = 2
ms = 3
gr()

# Flutter speed vs tip load for several structural couplings
plt1 = plot(xlabel="Tip Load [N]", ylabel="Flutter speed [m/s]", xlims=[0,35], ylims=[0,35])
plot!([NaN], [NaN], c=:black, lw=lw, legend=:bottomleft, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Patil et al. (2001)")
for (i,Ψ) in enumerate(ΨRange)
    plot!(F3Range, flutterSpeed[i,:], c=colors[i], lw = lw, label="\\Psi=$Ψ")
    if Ψ == -0.2
        scatter!(flutterSpeedVsTipLoadΨm02[1,:], flutterSpeedVsTipLoadΨm02[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif Ψ == 0.0
        scatter!(flutterSpeedVsTipLoadΨ0[1,:], flutterSpeedVsTipLoadΨ0[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif Ψ == +0.2
        scatter!(flutterSpeedVsTipLoadΨp02[1,:], flutterSpeedVsTipLoadΨp02[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/SMWFlutterStructuralCouplingRange.pdf"))

println("Finished SMWFlutterStructuralCouplingRangePlotGenerator.jl")