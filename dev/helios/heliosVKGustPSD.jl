using AeroBeams, DelimitedFiles, Interpolations

# Payload, gust duration and intensity
P = 150
τ = 360
γ = 0.15

# Seeds
seeds = 1:1:33

# Strings
Pstr = round(P)
τstr = round(τ,digits=1)
γstr = round(Int,γ*100)

# Set paths
relPath = "/dev/helios/figures/heliosVKGustPSDs"
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
XTICKS = [1e-3,1e-2,1e-1,1e0,1e1]
colors = cgrad(:rainbow, 2, categorical=true)
labels = ["Attached flow" "Dynamic stall"]

# Initialize
timeInGust = Array{Vector{Float64}}(undef)
rootAoAInGust = Array{Vector{Float64}}(undef)
airspeedInGust = Array{Vector{Float64}}(undef)
rootκ1InGust = Array{Vector{Float64}}(undef)
rootκ2InGust = Array{Vector{Float64}}(undef)
ωPSD_rootAoA = Array{Vector{Float64}}(undef,length(seeds),2)
ωPSD_airspeed = Array{Vector{Float64}}(undef,length(seeds),2)
ωPSD_rootκ1 = Array{Vector{Float64}}(undef,length(seeds),2)
ωPSD_rootκ2 = Array{Vector{Float64}}(undef,length(seeds),2)
PSD_rootAoA = Array{Vector{Float64}}(undef,length(seeds),2)
PSD_airspeed = Array{Vector{Float64}}(undef,length(seeds),2)
PSD_rootκ1 = Array{Vector{Float64}}(undef,length(seeds),2)
PSD_rootκ2 = Array{Vector{Float64}}(undef,length(seeds),2)
avgωPSD_rootAoA = Array{Vector{Float64}}(undef,2)
avgωPSD_airspeed = Array{Vector{Float64}}(undef,2)
avgωPSD_rootκ1 = Array{Vector{Float64}}(undef,2)
avgωPSD_rootκ2 = Array{Vector{Float64}}(undef,2)
avgPSD_rootAoA = Array{Vector{Float64}}(undef,2)
avgPSD_airspeed = Array{Vector{Float64}}(undef,2)
avgPSD_rootκ1 = Array{Vector{Float64}}(undef,2)
avgPSD_rootκ2 = Array{Vector{Float64}}(undef,2)
global ωmin = 0
global ωmax = Inf

# Read arrays
for (i,seed) in enumerate(seeds)
    for (j,aeroSolver) in enumerate(aeroSolvers)
        # Load arrays
        @load pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("timeInGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".jld2" timeInGust
        @load pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("rootAoAInGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".jld2" rootAoAInGust
        @load pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("airspeedInGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".jld2" airspeedInGust
        @load pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("rootκ1InGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".jld2" rootκ1InGust
        @load pkgdir(AeroBeams)*"/dev/helios/outputs/heliosVKGust/"*string("rootκ2InGust_P",Pstr,"_tau",τstr,"_gamma",γstr,"_seed",seed,"_",j)*".jld2" rootκ2InGust
        # Compute PSDs
        rootAoAωPSD,_,rootAoAPSD = get_FFT_and_PSD(timeInGust,rootAoAInGust)
        airspeedωPSD,_,airspeedPSD = get_FFT_and_PSD(timeInGust,airspeedInGust)
        rootκ1ωPSD,_,rootκ1PSD = get_FFT_and_PSD(timeInGust,rootκ1InGust)
        rootκ2ωPSD,_,rootκ2PSD = get_FFT_and_PSD(timeInGust,rootκ2InGust)
        # Assign
        ωPSD_rootAoA[i,j] = rootAoAωPSD
        ωPSD_airspeed[i,j] = airspeedωPSD
        ωPSD_rootκ1[i,j] = rootκ1ωPSD
        ωPSD_rootκ2[i,j] = rootκ2ωPSD
        PSD_rootAoA[i,j] = rootAoAPSD
        PSD_airspeed[i,j] = airspeedPSD
        PSD_rootκ1[i,j] = rootκ1PSD
        PSD_rootκ2[i,j] = rootκ2PSD
        # Update maximum and minimum frequencies
        global ωmin = max(ωmin, rootAoAωPSD[2], airspeedωPSD[2], rootκ1ωPSD[2], rootκ2ωPSD[2])
        global ωmax = min(ωmax, rootAoAωPSD[end], airspeedωPSD[end], rootκ1ωPSD[end], rootκ2ωPSD[end])
    end
end

# Common frequency grid for interpolation
dω = ωmin
ωPSD = ωmin:dω:ωmax

# Interpolate PSDs
for (i,seed) in enumerate(seeds)
    for (j,aeroSolver) in enumerate(aeroSolvers)
        PSD_rootAoA[i,j] = linear_interpolation(ωPSD_rootAoA[i,j], PSD_rootAoA[i,j])(ωPSD)
        PSD_airspeed[i,j] = linear_interpolation(ωPSD_airspeed[i,j], PSD_airspeed[i,j])(ωPSD)
        PSD_rootκ1[i,j] = linear_interpolation(ωPSD_rootκ1[i,j], PSD_rootκ1[i,j])(ωPSD)
        PSD_rootκ2[i,j] = linear_interpolation(ωPSD_rootκ2[i,j], PSD_rootκ2[i,j])(ωPSD)
    end
end

# Average PSDs
for (j,aeroSolver) in enumerate(aeroSolvers)
    avgPSD_rootAoA[j] = reduce(+, PSD_rootAoA[:,j])/(length(seeds))
    avgPSD_airspeed[j] = reduce(+, PSD_airspeed[:,j])/(length(seeds))
    avgPSD_rootκ1[j] = reduce(+, PSD_rootκ1[:,j])/(length(seeds))
    avgPSD_rootκ2[j] = reduce(+, PSD_rootκ2[:,j])/(length(seeds))
end

# Root AoA PSD
plt_AoAPSD = plot(xlabel="Frequency [rad/s]", ylabel="Root AoA PSD", tickfont=font(ts), guidefont=font(fs), xscale=:log10, yscale=:log10, xticks=XTICKS, legend=:topright, legendfontsize=lfs)
for (j,aeroSolver) in enumerate(aeroSolvers)
    plot!(ωPSD, avgPSD_rootAoA[j], c=colors[j], lw=lw, label=labels[j])
end
display(plt_AoAPSD)
savefig(plt_AoAPSD,string(absPath,"/heliosVKGustPSD_rootAoA_P",Pstr,"_tau",τstr,"_gamma",γstr,".pdf"))

# Airspeed PSD
plt_airspeedPSD = plot(xlabel="Frequency [rad/s]", ylabel="Airspeed PSD", tickfont=font(ts), guidefont=font(fs), xscale=:log10, yscale=:log10, xticks=XTICKS)
for (j,aeroSolver) in enumerate(aeroSolvers)
    plot!(ωPSD, avgPSD_airspeed[j], c=colors[j], lw=lw, label=false)
end
display(plt_airspeedPSD)
savefig(plt_airspeedPSD,string(absPath,"/heliosVKGustPSD_U_P",Pstr,"_tau",τstr,"_gamma",γstr,".pdf"))

# Root torsional curvature PSD
plt_k1PSD = plot(xlabel="Frequency [rad/s]", ylabel="\$\\kappa_1\$ PSD", tickfont=font(ts), guidefont=font(fs), xscale=:log10, yscale=:log10, xticks=XTICKS, legend=:topright, legendfontsize=lfs)
for (j,aeroSolver) in enumerate(aeroSolvers)
    plot!(ωPSD, avgPSD_rootκ1[j], c=colors[j], lw=lw, label=labels[j])
end
display(plt_k1PSD)
savefig(plt_k1PSD,string(absPath,"/heliosVKGustPSD_k1_P",Pstr,"_tau",τstr,"_gamma",γstr,".pdf"))

# Root bending curvature PSD
plt_k2PSD = plot(xlabel="Frequency [rad/s]", ylabel="\$\\kappa_2\$ PSD", tickfont=font(ts), guidefont=font(fs), xscale=:log10, yscale=:log10, xticks=XTICKS)
for (j,aeroSolver) in enumerate(aeroSolvers)
    plot!(ωPSD, avgPSD_rootκ2[j], c=colors[j], lw=lw, label=false)
end
display(plt_k2PSD)
savefig(plt_k2PSD,string(absPath,"/heliosVKGustPSD_k2_P",Pstr,"_tau",τstr,"_gamma",γstr,".pdf"))
