using Plots, DelimitedFiles, ColorSchemes

# k2
k2 = 0.045

# Airspeeds
URange = k2 == 0.03 ? [37, 37.5, 39.5, 40, 41.5] : [20, 23, 34, 35, 36.5]

# Initialize
t = Array{Vector{Float64}}(undef,length(URange))
wingAoA = Array{Vector{Float64}}(undef,length(URange))

# Load data
relPathFig = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_dynamic"
relPathData = "/dev/cHALE/Flexible/outputs/data/cHALEwing_dynamic"
absPathFig = string(pwd(),relPathFig)
absPathData= string(pwd(),relPathData)
for (i,U) in enumerate(URange)
    Ustr = rem(U,1) ≈ 0 ? string(round(Int,U)) : string(U)
    println(Ustr)
    t[i] = vec(readdlm(string(absPathData,"/cHALEwing_dynamic_t_U",Ustr,"_k2",k2,".txt")))
    wingAoA[i] = vec(readdlm(string(absPathData,"/cHALEwing_dynamic_AoA_U",Ustr,"_k2",k2,".txt")))
end

# Plot configurations
lw = 2
ts = 10
fs = 16
lfs = 12
tsz = 10
gr()
xlims = k2 == 0.03 ? [0,30] : [0,30]
ylims = [0,2]
legPos = k2 == 0.03 ? (0.15,0.3) : (0.2,0.3)

# Wing root AoA
colors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold])
plt_AoA = plot(xlabel="Time [s]", ylabel="Normalized tip AoA", xlims=xlims, ylims=ylims, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=legPos)
for i in eachindex(URange)
    plot!(t[i], wingAoA[i]./wingAoA[i][1], c=colors[i], lw=lw, label=string("\$U_{\\infty} = ",URange[i],"\$ m/s"))
end
display(plt_AoA)
savefig(string(absPathFig,string("/cHALEwing_dynamic_AoA_k2",k2,"_comparison.pdf")))

# # Show interactive plot
# plotlyjs()
# i2plot = 1
# plt_AoA2 = plot(xlabel="Time [s]", ylabel="Normalized wing root angle of attack")
# plot!(t[i2plot], wingAoA[i2plot]./wingAoA[i2plot][1], c=:black, lw=lw, label=false)
# display(plt_AoA2)

println("Finished cHALEwing_dynamic_plot_data.jl")