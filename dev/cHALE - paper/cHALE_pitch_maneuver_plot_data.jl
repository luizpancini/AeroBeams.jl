using Plots, DelimitedFiles, ColorSchemes
lw = 2
gr()

# k2
k2 = 0.030

# Airspeeds
URange = [30,32]

# Initialize
t = Array{Vector{Float64}}(undef,length(URange))
wingAoA = Array{Vector{Float64}}(undef,length(URange))
airspeed = Array{Vector{Float64}}(undef,length(URange))

# Load data
relPathFig = "/dev/outputs/figures/cHALE_pitch_maneuver"
relPathData = "/dev/outputs/data/cHALE_pitch_maneuver"
absPathFig = string(pwd(),relPathFig)
absPathData= string(pwd(),relPathData)
for i in eachindex(URange)
    t[i] = vec(readdlm(string(absPathData,"/cHALE_pitch_maneuver_t_U",URange[i],"_k2",k2,".txt")))
    wingAoA[i] = vec(readdlm(string(absPathData,"/cHALE_pitch_maneuver_AoA_U",URange[i],"_k2",k2,".txt")))
    airspeed[i] = vec(readdlm(string(absPathData,"/cHALE_pitch_maneuver_airspeed_U",URange[i],"_k2",k2,".txt")))
end

# Body AoA
colors = get(colorschemes[:rainbow], range(0, 1, length(URange)))
plt1 = plot(xlabel="Time [s]", ylabel="Normalized root wing angle of attack", xlims=[0,180],legendfontsize=12)
for i in eachindex(URange)
    plot!(t[i], wingAoA[i]./wingAoA[i][1], c=colors[i], lw=lw, label=string("\$U_{\\infty} = ",URange[i],"\$ m/s"))
end
display(plt1)
savefig(string(absPathFig,string("/cHALE_pitch_maneuver_AoA_k2",k2,"_comparison.pdf")))

# Airspeed
plt2 = plot(xlabel="Time [s]", ylabel="Normalized airspeed", xlims=[0,180])
for i in eachindex(URange)
    plot!(t[i], airspeed[i]./airspeed[i][1], c=colors[i], lw=lw, label=false)
end
display(plt2)
savefig(string(absPathFig,string("/cHALE_pitch_maneuver_airspeed_k2",k2,"_comparison.pdf")))

println("Finished cHALE_pitch_maneuver_plot_data.jl")