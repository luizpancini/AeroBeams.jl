using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazySteadySweepRangePitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazySteadySweepRangePitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
xΛ = [ 0.85, 0.95, 0.95, 0.95]
yΛ = [0.65,  0.43, 0.27, 0.15]
colors = cgrad(:rainbow, length(ΛRange), categorical=true)
ls = [:solid, :dash, :dot, :dashdot]
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 3
msw = 0
gr()

# Deformed position at θ=7 deg, U=70 m/s
plt_defPos = plot(xlabel="Normalized horizontal position", ylabel="Normalized vertical position", title=string("\$\\alpha_r=7^\\circ\$",", \$U=70\$ m/s"), xlims=[0,1], ylims=[0,0.8], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topleft)
plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[1], label="No loss")
plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[2], label="Exponential loss")
plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[3], label="VLM - undeformed")
plot!([NaN], [NaN], c=:black, lw=lw, ls=ls[4], label="VLM - deformed")
scatter!([NaN], [NaN], c=:black, ms=ms, msw=msw, label="Sharpy (VLM)")
for c in 1:length(tipLossTypeConfig)
    for i in eachindex(ΛRange)
        plot!(x1_def[c,i,end,end]/L, x3_def[c,i,end,end]/L, c=colors[i], lw=lw, ls=ls[c], label=false)
        annotate!([xΛ[i]],[yΛ[i]], text("\$\\Lambda=$(round(Int,ΛRange[i]*180/π)) ^\\circ\$", 12, colors[i]))
        if i == 1
            scatter!(dispΛ0θ7U70_Sharpy[1,:], dispΛ0θ7U70_Sharpy[2,:], c=colors[i], ms=ms, msw=msw, label=false)
        elseif i == 2
            scatter!(dispΛ10θ7U70_Sharpy[1,:], dispΛ10θ7U70_Sharpy[2,:], c=colors[i], ms=ms, msw=msw, label=false)
        elseif i == 3
            scatter!(dispΛ20θ7U70_Sharpy[1,:], dispΛ20θ7U70_Sharpy[2,:], c=colors[i], ms=ms, msw=msw, label=false)
        elseif i == 4
            scatter!(dispΛ30θ7U70_Sharpy[1,:], dispΛ30θ7U70_Sharpy[2,:], c=colors[i], ms=ms, msw=msw, label=false)    
        end
    end
end
display(plt_defPos)
# savefig(string(absPath,"/sweptPazySteadySweepRangePitchRange_defPos.pdf"))

println("Finished sweptPazySteadySweepRangePitchRangePlotGenerator.jl")