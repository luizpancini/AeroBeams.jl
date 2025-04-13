using Plots, ColorSchemes

# Run the script
include("../examples/straightWingTorsionalDivergence.jl")

# Plot configurations
lw = 2
ms = 4
gr()

# AoA distribution at divergence iminence
plt_AoA = plot(xlabel="Normalized arclength", ylabel="Angle of attack [deg]", title="\$U=$(round(URange[iD-1]/UD,digits=2))U_D\$", xlims=[0,1])
plot!(x1_e/L,AoA[iD-1], c=:black, ls=:solid, lw=lw, label="AeroBeams")
plot!(x,αRef, c=:black, ls=:dot, lw=lw, marker=:circle, ms=ms, label="Hodges and Pierce (2011)")
display(plt_AoA)
savefig(string(absPath,"/straightWingTorsionalDivergence_AoAatUD.pdf"))

# Tip AoA vs U/U_D
plt_tipAoA = plot(xlabel="\$U/U_D\$", ylabel="Tip angle of attack [deg]", xlims=[0,1])
plot!(URange[1:iD-1]/UD,tipAoA[1:iD-1], c=:black, ls=:solid, lw=lw, label="AeroBeams")
plot!(URange[1:iD-1]/UDRef,αTipRef[1:iD-1], c=:black, ls=:dot, lw=lw, marker=:circle, ms=ms, label="Hodges and Pierce (2011)")
display(plt_tipAoA)
savefig(string(absPath,"/straightWingTorsionalDivergence_AoA.pdf"))

# Tip twist angle vs q/q_D
plt_tipTwist = plot(xlabel="\$q/q_D\$", ylabel="Tip twist angle [deg]", xlims=[0,1], ylims=[0,12], legend=:topleft)
plot!((URange[1:iD-1]/UD).^2,tipAoA[1:iD-1].-θ*180/π, c=:black, ls=:solid, lw=lw, label="AeroBeams")
plot!((URange[1:iD-1]/UDRef).^2,θTipRef[1:iD-1], c=:black, ls=:dot, lw=lw, marker=:circle, ms=ms, label="Hodges and Pierce (2011)")
display(plt_tipTwist)
savefig(string(absPath,"/straightWingTorsionalDivergence_twist.pdf"))

# Damping plot
plt_damp = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", ylims=[-50,5],legend=:topleft)
for i in eachindex(URange)
    scatter!(URange[i]*ones(length(dampingsNonOscillatory[i])), dampingsNonOscillatory[i], c=:black, ms=ms, msw=0, label=false)
end
display(plt_damp)
savefig(string(absPath,"/straightWingTorsionalDivergence_damp.pdf"))

println("Finished straightWingTorsionalDivergence.jl")