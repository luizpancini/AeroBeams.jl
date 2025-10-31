using Plots

# Run the script
include("../examples/ChengSlantedDoublePendulum.jl")

# Set paths
relPath = "/test/outputs/figures/ChengSlantedDoublePendulum"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
for (i,Λ) in enumerate(ΛRange)
    anim = plot_dynamic_deformation(problem[i],plotDistLoads=false,plotFrequency=5,fps=20,scale=1,plotUndeformed=false,plotLimits=([-L,L],[-L,0],[-L,L]),save=true,savePath=string(relPath,"/ChengSlantedDoublePendulum_deformation_Λ",round(Int,Λ*180/π),"deg.gif"),displayProgress=true)
    display(anim)
end

# Plot configurations
lw = 2
ms = 4
msw = 0
ts = 12
fs = 16
lfs = 10
Λlabels = ["\$\\Lambda = $(round(Int,Λ*180/π))^\\circ\$" for Λ in ΛRange]
colors = cgrad(:rainbow, length(ΛRange), categorical=true)
gr()

# Normalized tip u1 displacement
plt_u1 = plot(xlabel="\$t\$ [s]", ylabel="\$u_1/L\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:top)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Cheng et al. (2025)")
for (i,Λ) in enumerate(ΛRange)
    plot!(t[i], u1_tip[i]/L, lw=lw, c=colors[i], label=Λlabels[i])
    if Λ == 0
        scatter!(Λ0_r1_ref[1,:], Λ0_r1_ref[2,:].-L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif Λ == π/4
        scatter!(Λ45_r1_ref[1,:], Λ45_r1_ref[2,:].-L, ms=ms, msw=msw, c=colors[i], label=false)
    elseif Λ == π/2
        scatter!(Λ90_r1_ref[1,:], Λ90_r1_ref[2,:].-L, ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_u1)
savefig(string(absPath,"/ChengSlantedDoublePendulum_u1.pdf"))

# Normalized tip u2 displacement
plt_u2 = plot(xlabel="\$t\$ [s]", ylabel="\$u_2/L\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,Λ) in enumerate(ΛRange)
    plot!(t[i], u2_tip[i]/L, lw=lw, c=colors[i], label=false)
    if Λ == π/4
        scatter!(Λ45_r2_ref[1,:], Λ45_r2_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_u2)
savefig(string(absPath,"/ChengSlantedDoublePendulum_u2.pdf"))

# Normalized tip u3 displacement
plt_u3 = plot(palette=colors, xlabel="\$t\$ [s]", ylabel="\$u_3/L\$", tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for (i,Λ) in enumerate(ΛRange)
    plot!(t[i], u3_tip[i]/L, lw=lw, c=colors[i], label=false)
    if Λ == 0
        scatter!(Λ0_r3_ref[1,:], Λ0_r3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif Λ == π/4
        scatter!(Λ45_r3_ref[1,:], Λ45_r3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    elseif Λ == π/2
        scatter!(Λ90_r3_ref[1,:], Λ90_r3_ref[2,:], ms=ms, msw=msw, c=colors[i], label=false)
    end
end
display(plt_u3)
savefig(string(absPath,"/ChengSlantedDoublePendulum_u3.pdf"))

println("Finished ChengSlantedDoublePendulumPlotGenerator.jl")