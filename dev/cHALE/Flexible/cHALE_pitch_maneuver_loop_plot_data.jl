using Plots, DelimitedFiles, ColorSchemes

# k2
k2 = 0.045

# Airspeeds
URange = k2 == 0.045 ? float([24,24.5,33.5]) : float([46,47])

# Initialize
t = Array{Vector{Float64}}(undef,length(URange))
wingAoA = Array{Vector{Float64}}(undef,length(URange))

# Load data
relPathFig = "/dev/cHALE/Flexible/outputs/figures/cHALE_pitch_maneuver_loop"
relPathData = "/dev/cHALE/Flexible/outputs/data/cHALE_pitch_maneuver_loop"
absPathFig = string(pwd(),relPathFig)
absPathData= string(pwd(),relPathData)
for (i,U) in enumerate(URange)
    t[i] = vec(readdlm(string(absPathData,"/cHALE_pitch_maneuver_loop_t_k2",k2,"_U",U,".txt")))
    wingAoA[i] = vec(readdlm(string(absPathData,"/cHALE_pitch_maneuver_loop_AoA_k2",k2,"_U",U,".txt")))
end

# Plot configurations
lw = 2
ts = 10
fs = 16
lfs = 12
tsz = 10
gr()
xlims = k2 == 0.045 ? [0,121] : [0,121]
ylims = k2 == 0.045 ? [0,2] : [0.75,1.25]

# Wing root AoA
colors = palette([:royalblue, :deeppink, :darkorange], length(URange))
plt_AoA = plot(xlabel="Time [s]", ylabel="Normalized root AoA", xticks=0:30:xlims[2], xlims=xlims, ylims=ylims, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:topright)
for i in eachindex(URange)
    plot!(t[i], wingAoA[i]./wingAoA[i][1], c=colors[i], lw=lw, label="\$U_{\\infty} = $(URange[i])\\;\\mathrm{m/s}\$")
end
if k2 == 0.045
    textPos1 = [97.25, 1.25]
    annotate!(textPos1[1], textPos1[2], text("\$\\tau\\approx20\$ s", 14, colors[2]))
    quiver!([textPos1[1]-15,textPos1[1]+15], [textPos1[2],textPos1[2]], quiver=([+5,-5], [0,0]), arrow=:closed, linecolor=colors[2])
    plot!(t[3],wingAoA[3]./wingAoA[3][1], c=colors[3], lw=lw, inset=(1, bbox(0.03, 0.08, 0.47, 0.3, :bottom, :right)), xlims=[9,10], ylims=[1.055,1.06], label=nothing, subplot=2, bg_inside=nothing)
    textPos2 = [83, 0.73]
    annotate!(textPos2[1], textPos2[2], text("\$\\tau\\approx0.325\$ s", 12, colors[3]))
    quiver!([textPos2[1]-14,textPos2[1]+14], [textPos2[2]-0.1,textPos2[2]-0.1], quiver=([+5,-5], [0,0]), arrow=:closed, linecolor=colors[3])
elseif k2 == 0.015
    plot!(t[2],wingAoA[2]./wingAoA[2][1], c=colors[2], lw=lw, inset=(1, bbox(0.03, 0.08, 0.47, 0.3, :bottom, :right)), xlims=[20,21], ylims=[0.9875,0.99], label=nothing, subplot=2, bg_inside=nothing)
    textPos2 = [87, 0.9]
    annotate!(textPos2[1], textPos2[2], text("\$\\tau\\approx0.25\$ s", 12, colors[2]))
    quiver!([textPos2[1]-10,textPos2[1]+10], [textPos2[2]-0.02,textPos2[2]-0.02], quiver=([+3,-3], [0,0]), arrow=:closed, linecolor=colors[2])
end
display(plt_AoA)
savefig(string(absPathFig,string("/cHALE_pitch_maneuver_AoA_k2",k2,"_comparison.pdf")))

# # Show interactive plot
# plotlyjs()
# i2plot = 1
# plt_AoA2 = plot(xlabel="Time [s]", ylabel="Normalized wing root angle of attack")
# plot!(t[i2plot], wingAoA[i2plot]./wingAoA[i2plot][1], c=:black, lw=lw, label=false)
# display(plt_AoA2)

println("Finished cHALE_pitch_maneuver_plot_data.jl")