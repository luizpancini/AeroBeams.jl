using AeroBeams, DelimitedFiles

# Payload, airspeed, gust duration and intensity
P = 150
U = 40*.3048
τ = 360
γ = 0.1

# Seeds
seeds = 1:1:21

# Seed groups
# seedGroups = [seeds[1:3], seeds[1:9], seeds[1:21], seeds[1:27], seeds[1:33]]
seedGroups = [seeds[sort(unique(vcat(1,6,12)))], seeds[sort(unique(vcat(1,6,12,1:5)))], seeds[sort(unique(vcat(1,6,12,1:8)))], seeds[1:15], seeds[1:21]]

# Markers for frequencies of exceedance
rootAoAMarkers = vcat(0:1:90)
root_κ1Markers = vcat(0:5e-5:1e-2)
root_κ2Markers = vcat(0:1e-4:5e-2)
tip_u3Markers = vcat(0:1:50)

# Strings
Pstr = round(P)
τstr = round(τ,digits=1)
γstr = round(Int,γ*100)

# Set paths
relPath = "/dev/helios/figures/heliosVKGust"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
using Plots, ColorSchemes
gr()
ts = 10
fs = 15
lfs = 12
lw = 2
ms = 3
msw = 0
colors = cgrad(:rainbow, length(seedGroups), categorical=true)

# Initialize
COE_rootAoA = Array{Vector{Float64}}(undef,length(seeds),2)
COE_rootκ1 = Array{Vector{Float64}}(undef,length(seeds),2)
COE_rootκ2 = Array{Vector{Float64}}(undef,length(seeds),2)
COE_tipu3 = Array{Vector{Float64}}(undef,length(seeds),2)
gustDur = Array{Float64}(undef,length(seeds),2)

# Read counts of exceedance data
for (i,seed) in enumerate(seeds)
    # Sweep aero solvers
    for j in 1:2
        COE_rootAoA[i,j] = vec(readdlm(pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("COE_rootAoA_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt"))
        COE_rootκ1[i,j] = vec(readdlm(pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("COE_rootκ1_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt"))
        COE_rootκ2[i,j] = vec(readdlm(pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("COE_rootκ2_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt"))
        COE_tipu3[i,j] = vec(readdlm(pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("COE_tipu3_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt"))
        gustDur[i,j] = readdlm(pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("gustDur_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt")[1]
    end
end

# Compute frequencies of exceedance for each group [/min]
FOE_rootAoA_AF = Array{Vector{Float64}}(undef,length(seedGroups))
FOE_rootAoA_DS = Array{Vector{Float64}}(undef,length(seedGroups))
for (i,seedGroup) in enumerate(seedGroups)
    FOE_rootAoA_AF[i] = reduce(+, COE_rootAoA[seedGroup, 1]./gustDur[seedGroup,1])/(length(seedGroup)/60)
    FOE_rootAoA_DS[i] = reduce(+, COE_rootAoA[seedGroup, 2]./gustDur[seedGroup,2])/(length(seedGroup)/60)
    FOE_rootAoA_AF[i][FOE_rootAoA_AF[i] .== 0] .= NaN
    FOE_rootAoA_DS[i][FOE_rootAoA_DS[i] .== 0] .= NaN
end

# Plot frequencies of exceedance - attached flow model
plt_AoA_AF = plot(xlabel="Root AoA above trim value [deg]", ylabel="Frequency of exceedance [/min]", xlims=[0,45], ylims=[1e-3,1e2], yscale=:log10, tickfont=font(ts), guidefont=font(fs), legend=:topright, legendfontsize=lfs)
for (i,seedGroup) in enumerate(seedGroups)
    plot!(rootAoAMarkers, FOE_rootAoA_AF[i], c=colors[i], lw=lw, label=string("\$\\overline{T}_g=\$",round(sum(gustDur[seedGroup,1])/60,digits=1)," min"))
end
display(plt_AoA_AF)
savefig(plt_AoA_AF,string(absPath,"/heliosVKGustFOEconvAF_AoA_P",Pstr,"_tau",τstr,"_gamma",γstr,".pdf"))

# Plot frequencies of exceedance - dynamic stall model
plt_AoA_DS = plot(xlabel="Root AoA above trim value [deg]", ylabel="Frequency of exceedance [/min]", xlims=[0,45], ylims=[1e-3,1e2], yscale=:log10, tickfont=font(ts), guidefont=font(fs), legend=:topright, legendfontsize=lfs)
for (i,seedGroup) in enumerate(seedGroups)
    plot!(rootAoAMarkers, FOE_rootAoA_DS[i], c=colors[i], lw=lw, label=string("\$\\overline{T}_g=\$",round(sum(gustDur[seedGroup,2])/60,digits=1)," min"))
end
display(plt_AoA_DS)
savefig(plt_AoA_DS,string(absPath,"/heliosVKGustFOEconvDS_AoA_P",Pstr,"_tau",τstr,"_gamma",γstr,".pdf"))