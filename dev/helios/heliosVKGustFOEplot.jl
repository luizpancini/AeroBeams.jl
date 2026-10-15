using AeroBeams, DelimitedFiles

# Payload, gust duration and intensity
P = 150
τ = 360
γ = 0.1

# Seeds
seeds = 1:1:21

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
colors = cgrad(:rainbow, 2, categorical=true)

# Initialize
COE_rootAoA = Array{Vector{Float64}}(undef,length(seeds),2)
COE_rootκ1 = Array{Vector{Float64}}(undef,length(seeds),2)
COE_rootκ2 = Array{Vector{Float64}}(undef,length(seeds),2)
COE_tipu3 = Array{Vector{Float64}}(undef,length(seeds),2)
gustDur = Array{Float64}(undef,length(seeds),2)

# Read counts of exceedance data
for (i,seed) in enumerate(seeds)
    # Sweep aero solvers
    for (j,aeroSolver) in enumerate(aeroSolvers)
        COE_rootAoA[i,j] = vec(readdlm(pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("COE_rootAoA_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt"))
        COE_rootκ1[i,j] = vec(readdlm(pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("COE_rootκ1_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt"))
        COE_rootκ2[i,j] = vec(readdlm(pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("COE_rootκ2_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt"))
        COE_tipu3[i,j] = vec(readdlm(pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("COE_tipu3_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt"))
        gustDur[i,j] = readdlm(pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("gustDur_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".txt")[1]
    end
end

# Compute frequencies of exceedance [/min]
FOE_rootAoA_AF = reduce(+, COE_rootAoA[:, 1]./gustDur[:,1])/(length(seeds)/60)
FOE_rootAoA_DS = reduce(+, COE_rootAoA[:, 2]./gustDur[:,2])/(length(seeds)/60)
FOE_rootκ1_AF = reduce(+, COE_rootκ1[:, 1]./gustDur[:,1])/(length(seeds)/60)
FOE_rootκ1_DS = reduce(+, COE_rootκ1[:, 2]./gustDur[:,2])/(length(seeds)/60)
FOE_rootκ2_AF = reduce(+, COE_rootκ2[:, 1]./gustDur[:,1])/(length(seeds)/60)
FOE_rootκ2_DS = reduce(+, COE_rootκ2[:, 2]./gustDur[:,2])/(length(seeds)/60)
FOE_tipu3_AF = reduce(+, COE_tipu3[:, 1]./gustDur[:,1])/(length(seeds)/60)
FOE_tipu3_DS = reduce(+, COE_tipu3[:, 2]./gustDur[:,2])/(length(seeds)/60)

# Set zeros to NaNs
FOE_rootAoA_AF[FOE_rootAoA_AF .== 0] .= NaN
FOE_rootAoA_DS[FOE_rootAoA_DS .== 0] .= NaN
FOE_rootκ1_AF[FOE_rootκ1_AF .== 0] .= NaN
FOE_rootκ1_DS[FOE_rootκ1_DS .== 0] .= NaN
FOE_rootκ2_AF[FOE_rootκ2_AF .== 0] .= NaN
FOE_rootκ2_DS[FOE_rootκ2_DS .== 0] .= NaN
FOE_tipu3_AF[FOE_tipu3_AF .== 0] .= NaN
FOE_tipu3_DS[FOE_tipu3_DS .== 0] .= NaN

# Plot frequencies of exceedance
plt_AoA = plot(xlabel="Root AoA above trim value [deg]", ylabel="Frequency of exceedance [/min]", xlims=[0,45], ylims=[1e-3,1e2], yscale=:log10, tickfont=font(ts), guidefont=font(fs), legend=:topright, legendfontsize=lfs)
plot!(rootAoAMarkers, FOE_rootAoA_AF, c=colors[1], lw=lw, label="Attached flow")
plot!(rootAoAMarkers, FOE_rootAoA_DS, c=colors[2], lw=lw, label="Dynamic stall")
display(plt_AoA)
savefig(plt_AoA,string(absPath,"/heliosVKGustFOE_AoA_P",Pstr,"_tau",τstr,"_gamma",γstr,".pdf"))

plt_κ1 = plot(xlabel="Root torsional curvature above trim value [1/m]", ylabel="Frequency of exceedance [/min]", xlims=[0,0.0021], ylims=[1e-3,1e2], yscale=:log10, tickfont=font(ts), guidefont=font(fs))
plot!(root_κ1Markers, FOE_rootκ1_AF, c=colors[1], lw=lw, label=false)
plot!(root_κ1Markers, FOE_rootκ1_DS, c=colors[2], lw=lw, label=false)
display(plt_κ1)
savefig(plt_κ1,string(absPath,"/heliosVKGustFOE_k1_P",Pstr,"_tau",τstr,"_gamma",γstr,".pdf"))

plt_κ2 = plot(xlabel="Root bending curvature above trim value [1/m]", ylabel="Frequency of exceedance [/min]", xlims=[0,0.01], ylims=[1e-3,1e2], yscale=:log10, tickfont=font(ts), guidefont=font(fs))
plot!(root_κ2Markers, FOE_rootκ2_AF, c=colors[1], lw=lw, label=false)
plot!(root_κ2Markers, FOE_rootκ2_DS, c=colors[2], lw=lw, label=false)
display(plt_κ2)
savefig(plt_κ2,string(absPath,"/heliosVKGustFOE_k2_P",Pstr,"_tau",τstr,"_gamma",γstr,".pdf"))

plt_u3 = plot(xlabel="Tip OOP disp. above trime value [% semispan]", ylabel="Frequency of exceedance [/min]", xlims=[0,25], ylims=[1e-3,1e2], yscale=:log10, tickfont=font(ts), guidefont=font(fs))
plot!(tip_u3Markers, FOE_tipu3_AF, c=colors[1], lw=lw, label=false)
plot!(tip_u3Markers, FOE_tipu3_DS, c=colors[2], lw=lw, label=false)
display(plt_u3)
savefig(plt_u3,string(absPath,"/heliosVKGustFOE_u3_P",Pstr,"_tau",τstr,"_gamma",γstr,".pdf"))