using Plots, DelimitedFiles, ColorSchemes

# k2
k2 = 0.030

# Airspeeds
URange = k2 == 0.045 ? [33,25,24] : [30,32]

# Initialize
t = Array{Vector{Float64}}(undef,length(URange))
wingAoA = Array{Vector{Float64}}(undef,length(URange))

# Load data
relPathFig = "/dev/cHALE/outputs/figures/cHALE_pitch_maneuver"
relPathData = "/dev/cHALE/outputs/data/cHALE_pitch_maneuver"
absPathFig = string(pwd(),relPathFig)
absPathData= string(pwd(),relPathData)
for i in eachindex(URange)
    t[i] = vec(readdlm(string(absPathData,"/cHALE_pitch_maneuver_t_U",URange[i],"_k2",k2,".txt")))
    wingAoA[i] = vec(readdlm(string(absPathData,"/cHALE_pitch_maneuver_AoA_U",URange[i],"_k2",k2,".txt")))
end

# Plot configurations
lw = 2
ts = 10
fs = 16
lfs = 12
tsz = 10
gr()
xlims = k2 == 0.045 ? [0,120] : [0,180]
ylims = k2 == 0.045 ? [0,2] : [0.5,1.5]

# Wing root AoA
colors = get(colorschemes[:rainbow], range(0, 1, length(URange)))
plt_AoA = plot(xlabel="Time [s]", ylabel="Normalized root AoA", xlims=xlims, ylims=ylims, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
for i in eachindex(URange)
    plot!(t[i], wingAoA[i]./wingAoA[i][1], c=colors[i], lw=lw, label=string("\$U_{\\infty} = ",URange[i],"\$ m/s"))
end
if k2 == 0.045
    textPos1 = [97.25, 1.25]
    annotate!(textPos1[1], textPos1[2], text("\$\\tau\\approx20\$ s", 14, :green))
    quiver!([textPos1[1]-15,textPos1[1]+15], [textPos1[2],textPos1[2]], quiver=([+5,-5], [0,0]), arrow=:closed, linecolor=colors[2])
    plot!(t[1],wingAoA[1]./wingAoA[1][1], c=colors[1], lw=lw, inset=(1, bbox(0.03, 0.08, 0.47, 0.3, :bottom, :right)), xlims=[8,10], ylims=[1.02,1.06], label=nothing, subplot=2, bg_inside=nothing)
    textPos2 = [95, 0.67]
    annotate!(textPos2[1], textPos2[2], text("\$\\tau\\approx0.37\$ s", 12, :purple))
    quiver!([textPos2[1]-10,textPos2[1]+10], [textPos2[2]-0.1,textPos2[2]-0.1], quiver=([+5,-5], [0,0]), arrow=:closed, linecolor=colors[1])
elseif k2 == 0.030
    textPos1 = [91, 1.2]
    annotate!(textPos1[1], textPos1[2], text("\$\\tau\\approx23.3\$ s", 14, :purple))
    quiver!([textPos1[1]-19,textPos1[1]+15], [textPos1[2]-0.07,textPos1[2]-0.07], quiver=([+5,-5], [0,0]), arrow=:closed, linecolor=colors[1])
    plot!(plt_AoA, inset=(1, bbox(0.03, 0.08, 0.47, 0.3, :bottom, :right)), xlims=[31.5,33.5], ylims=[1.125,1.14], subplot=2, bg_inside=nothing)
    plot!(plt_AoA[2], t[2],wingAoA[2]./wingAoA[2][1], c=colors[2], lw=lw, label=false)
    textPos2 = [32.6, 1.128]
    annotate!(plt_AoA[2], textPos2[1], textPos2[2], text("\$\\tau\\approx0.47\$ s", 12, :red))
    quiver!(plt_AoA[2], [textPos2[1]-0.6,textPos2[1]+0.48], [textPos2[2]+3e-3,textPos2[2]+3e-3], quiver=([+0.3,-0.3], [0,0]), arrow=:closed, linecolor=colors[2])
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