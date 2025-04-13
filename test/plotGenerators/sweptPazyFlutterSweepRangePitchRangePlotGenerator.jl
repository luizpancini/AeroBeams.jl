using Plots, ColorSchemes

# Run the script
include("../examples/sweptPazyFlutterSweepRangePitchRange.jl")

# Set paths
relPath = "/test/outputs/figures/sweptPazyFlutterSweepRangePitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
configColors = cgrad([:red,:green,:yellow,:blue], length(tipLossTypeConfig), categorical=true)
ΛColors = cgrad([:red,:green,:yellow,:blue], length(ΛRange), categorical=true)
modeColors = cgrad([:blue,:red,:green], nModes, categorical=true)
configLabels = ["No loss" "Exponential loss" "VLM - undeformed" "VLM - deformed"]
ts = 10
fs = 16
lfs = 10
lw = 2
ms = 4
msw = 0
colorFeniax = :purple
colorSharpy = :olive
gr()

# Initialize plot variables
speedOnsetHump = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange))
speedOffsetHump = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange))
speedOnsetHard = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange))
speedOffsetHard = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange))
dispOnsetHump = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange))
dispOffsetHump = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange))
dispOnsetHard = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange))
dispOffsetHard = Array{Vector{Float64}}(undef,length(tipLossTypeConfig),length(ΛRange))

# Sweep configurations
for (c,hasTipCorrection,tipLossType) in zip(1:length(tipLossTypeConfig),hasTipCorrectionConfig,tipLossTypeConfig)

    # Loop angle of sweep
    for (i,Λ) in enumerate(ΛRange)

        Λstr = string(round(Int,Λ*180/pi))

        ULim = 120

        # Set sweep-dependent plot configurations
        humpN = [1 for j in eachindex(θRange)]
        hardN = [1 for j in eachindex(θRange)]
        if tipLossType == "None"
            if Λ == 0
                dispLim = 70
                xθ = [26, 32, 35, 37, 39, 40]
                yθ = [72, 67, 60, 55, 49, 43]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = zeros(Int,length(θRange))
                for (j,θ) in enumerate(θRange)
                    hardMode[j] = θ >= 6.0*π/180 ? 2 : 2
                    hardN[j] = θ >= 6.0*π/180 ? 1 : 1
                end
            elseif Λ == 10*π/180
                dispLim = 40
                xθ = [2.5,  7, 15, 21.5, 33, 37]
                yθ = [ 85, 85, 83,   80, 77, 64]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = [2 for j in eachindex(θRange)]
            elseif Λ == 20*π/180
                dispLim = 23.5
                xθ = [ 1,  3, 6.5, 10, 16, 21]
                yθ = [73, 73,  73, 72, 64, 61]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = zeros(Int,length(θRange))
                for (j,θ) in enumerate(θRange)
                    hardMode[j] = θ >= 6.51*π/180 ? 1 : 1
                end
            elseif Λ == 30*π/180
                dispLim = 15.3
                xθ = [0.7, 1.8,  4, 6.5, 10.5, 14]
                yθ = [ 74,  74, 73,  72,   67, 64]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = [1 for j in eachindex(θRange)]   
            end
        elseif tipLossType == "Exponential"
            if Λ == 0
                dispLim = 70
                xθ = [27, 36, 33, 35, 37, 39]
                yθ = [85, 82, 70, 64, 58, 52]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = zeros(Int,length(θRange))
                for (j,θ) in enumerate(θRange)
                    hardMode[j] = ifelse(θ <= 0.1*π/180,1,ifelse(θ >= 6.5*π/180,3,2))
                    hardN[j] = ifelse(θ <= 0.1*π/180,1,ifelse(θ >= 6.5*π/180,2,1))
                end
            elseif Λ == 10*π/180
                dispLim = 40
                xθ = [ 2,  6, 12, 18, 23, 25]
                yθ = [80, 79, 75, 73, 60, 51]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = [1 for j in eachindex(θRange)]
            elseif Λ == 20*π/180
                dispLim = 23.5
                xθ = [ 1,  3,  6,  9, 15, 21]
                yθ = [72, 72, 71, 70, 68, 67]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = [1 for j in eachindex(θRange)]
            elseif Λ == 30*π/180
                dispLim = 15.3
                xθ = [0.7, 1.8, 3.5,  6, 9.5, 13.5]
                yθ = [ 72,  72,  71, 70,  68,   67]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = [1 for j in eachindex(θRange)] 
            end
        elseif tipLossType == "VLM-undef"
            if Λ == 0
                dispLim = 70
                xθ = [27, 33, 32, 35, 37, 39]
                yθ = [87, 82, 70, 64, 56, 51]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = zeros(Int,length(θRange))
                for (j,θ) in enumerate(θRange)
                    hardMode[j] = θ >= 6.5*π/180 ? 2 : 2
                    hardN[j] = θ >= 6.5*π/180 ? 1 : 1
                end
            elseif Λ == 10*π/180
                dispLim = 40
                xθ = [ 2,  6, 12, 19, 27, 33]
                yθ = [85, 84, 83, 82, 75, 67]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = zeros(Int,length(θRange))
                for (j,θ) in enumerate(θRange)
                    hardMode[j] = θ >= 6*π/180 ? 2 : 1
                end
            elseif Λ == 20*π/180
                dispLim = 23.5
                xθ = [ 1,  3,  6, 9.5, 15, 20]
                yθ = [80, 80, 79,  78, 75, 72]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = [1 for j in eachindex(θRange)]
            elseif Λ == 30*π/180
                dispLim = 15
                xθ = [0.7, 1.8,  3.8,  6, 10, 13.5]
                yθ = [ 85,  85, 84, 83, 81,   78]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = [1 for j in eachindex(θRange)]   
            end    
        elseif tipLossType == "VLM-def"
            if Λ == 0
                dispLim = 70
                xθ = [30, 35, 33, 36, 37, 38]
                yθ = [90, 83, 70, 63, 57, 52]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = zeros(Int,length(θRange))
                for (j,θ) in enumerate(θRange)
                    hardMode[j] = θ >= 6.5*π/180 ? 2 : 2
                    hardN[j] = θ >= 6.5*π/180 ? 1 : 1
                end
            elseif Λ == 10*π/180
                dispLim = 40
                xθ = [ 2,  6, 12, 18, 23, 25]
                yθ = [85, 84, 80, 78, 65, 55]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = zeros(Int,length(θRange))
                for (j,θ) in enumerate(θRange)
                    hardMode[j] = θ >= 5.01*π/180 ? 1 : 1
                end
            elseif Λ == 20*π/180
                dispLim = 23.5
                xθ = [ 1,  3,  6, 9.5, 15, 20]
                yθ = [85, 85, 84,  82, 80, 75]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = [1 for j in eachindex(θRange)]
            elseif Λ == 30*π/180
                dispLim = 15
                xθ = [0.7, 1.8, 3.8,  6, 10, 13.5]
                yθ = [ 90,  90,  89, 87, 85,   82]
                humpMode = [3 for j in eachindex(θRange)]
                hardMode = [1 for j in eachindex(θRange)]
            end    
        end

        # Flutter onset/offset speeds and displacements
        speedOnsetHump[c,i] = [flutterOnsetSpeedsOfMode[c,i,j,humpMode[j]][humpN[j]] for j in eachindex(θRange)]
        speedOffsetHump[c,i] = [flutterOffsetSpeedsOfMode[c,i,j,humpMode[j]][humpN[j]] for j in eachindex(θRange)]
        speedOnsetHard[c,i] = [flutterOnsetSpeedsOfMode[c,i,j,hardMode[j]][hardN[j]] for j in eachindex(θRange)]
        speedOffsetHard[c,i] = [150 for j in eachindex(θRange)]
        dispOnsetHump[c,i] = [flutterOnsetDispOfMode[c,i,j,humpMode[j]][humpN[j]] for j in eachindex(θRange)]
        dispOffsetHump[c,i] = [flutterOffsetDispOfMode[c,i,j,humpMode[j]][humpN[j]] for j in eachindex(θRange)]
        dispOnsetHard[c,i] = [flutterOnsetDispOfMode[c,i,j,hardMode[j]][hardN[j]] for j in eachindex(θRange)]
        dispOffsetHard = dispOnsetHard

        if tipLossType == "None"
            if Λ == 0*π/180
                speedOnsetHump[c,i][1:4] .= NaN
                speedOffsetHump[c,i][1:4] .= NaN
                speedOnsetHump[c,i][end-5:end] .= NaN
                speedOffsetHump[c,i][end-5:end] .= NaN
                dispOnsetHump[c,i][1:4] .= NaN
                dispOffsetHump[c,i][1:4] .= NaN
                dispOnsetHump[c,i][end-5:end] .= NaN
                dispOffsetHump[c,i][end-5:end] .= NaN
            end
        elseif tipLossType == "Exponential" 
            if Λ == 30*π/180
                speedOffsetHump[c,i][end] = 150
            end
        end

        # V-g-f
        θ2plot = [0 1 3 5 7]*π/180
        plt_Vf = plot(ylabel="Frequency [Hz]", xlims=[0,ULim], ylims=[0,50], tickfont=font(ts), guidefont=font(fs))
        for (j,θ) in enumerate(θRange)
            if !any(y -> isapprox(θ, y; atol=1e-4), θ2plot)
                continue
            end
            for mode in 1:nModes
                plot!(URange, modeFrequencies[c,i,j,mode]/(2π), c=modeColors[mode], lw=lw, alpha=0.15+0.85*j/length(θRange), label=false)
            end
        end
        plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[0,ULim], ylims=[-0.2,0.1], tickfont=font(ts), guidefont=font(fs), legend=:topleft)
        for (j,θ) in enumerate(θRange)
            if !any(y -> isapprox(θ, y; atol=1e-4), θ2plot)
                continue
            end
            for mode in 1:nModes
                plot!(URange, modeDampingRatios[c,i,j,mode], c=modeColors[mode], lw=lw, alpha=0.15+0.85*j/length(θRange), label=false)
            end
        end
        plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
        display(plt_Vgf)
        savefig(string(absPath,"/sweptPazyFlutterSweepRangePitchRange_Vgf_",tipLossType,"_Lambda",Λstr,".pdf"))

        # Lines of flutter onset/offset speeds
        x_hump = vcat(speedOnsetHump[c,i],reverse(speedOffsetHump[c,i]))
        y_hump = vcat(θRange*180/π,reverse(θRange*180/π))
        nanIndices = findall(isnan, x_hump)
        x_hump = x_hump[.!isnan.(x_hump)]
        y_hump = y_hump[setdiff(1:length(y_hump), nanIndices)]
        x_hard = vcat(speedOnsetHard[c,i],reverse(speedOffsetHard[c,i]))
        y_hard = vcat(θRange*180/π,reverse(θRange*180/π))

        # Flutter onset and offset speeds vs root pitch angle
        plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[0,ULim], ylims=[0,7], xticks=collect(0:20:ULim), yticks=collect(0:1:7), legend=:bottomleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
        plot!(Shape(x_hump,y_hump), fillcolor=plot_color(:green, 0.25), lw=lw, label="Hump flutter region")
        plot!(Shape(x_hard,y_hard), fillcolor=plot_color(:red, 0.25), lw=lw, label="Hard flutter region")
        if i > 1
            plot!(legend=false)
        end
        display(plt1)
        savefig(string(absPath,"/sweptPazyFlutterSweepRangePitchRange_",tipLossType,"_Lambda",Λstr,"_flutterBoundaryPitch.pdf"))

        # Lines of flutter onset/offset displacements
        x_hump = vcat(dispOnsetHump[c,i],reverse(dispOffsetHump[c,i]))
        y_hump = vcat(speedOnsetHump[c,i],reverse(speedOffsetHump[c,i]))
        nanIndices = findall(isnan, x_hump)
        x_hump = x_hump[.!isnan.(x_hump)]
        y_hump = y_hump[setdiff(1:length(y_hump), nanIndices)]
        x_hard = vcat(dispOnsetHard[c,i],reverse(dispOffsetHard[c,i]))
        y_hard = vcat(speedOnsetHard[c,i],reverse(speedOffsetHard[c,i]))

        # Flutter onset and offset speeds vs tip OOP displacement for varying root pitch angle
        θ2plot = [0.5; 1; 2; 3; 5; 7]*π/180
        indθ2plot = findall(x -> any(isapprox(x, θ; atol=1e-10) for θ in θ2plot), θRange)
        θcolors = cgrad([:blue, :red], length(θ2plot), categorical=true)
        plt2 = plot(xlabel="Tip OOP displacement [% semispan]", ylabel="Airspeed [m/s]", xlims=[0,dispLim], ylims=[30,ULim], xticks=collect(0:5:dispLim), yticks=collect(30:10:ULim), legend=:bottomright, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
        plot!(Shape(x_hump,y_hump), fillcolor=plot_color(:green, 0.25), lw=lw, label="Hump flutter region")
        plot!(Shape(x_hard,y_hard), fillcolor=plot_color(:red, 0.25), lw=lw, label="Hard flutter region")
        for (n,ind) in enumerate(indθ2plot)
            if Λ == 0
                if tipLossType == "None"
                    U90ind = findfirst(x-> x > 90, URange) - 1
                    range = n == 1 ? (1:U90ind) : (1:length(URange))
                else
                    U110ind = findfirst(x-> x > 110, URange) - 1
                    range = n == 1 ? (1:U110ind) : (1:length(URange))
                end
            else
                range = 1:length(URange)
            end    
            plot!(tipOOP[c,i,ind,range]/L*100, URange[range], c=θcolors[n], ls=:dash, lw=lw, label=false)
            if i == 1
                annotate!([xθ[n]],[yθ[n]], text("\$$(round(θ2plot[n]*180/π,digits=2)) ^\\circ\$", 10, :bottom, θcolors[n]))
            end
        end
        if i==1
            plot!(flutterBoundary_dispVsU_Lambda0_Sharpy[1,:]*100,flutterBoundary_dispVsU_Lambda0_Sharpy[2,:], lw=1, ls=:dash, c=colorSharpy, marker=:circle, ms=ms, msw=msw, label="Sharpy")
            plot!(flutterBoundary_dispVsU_Lambda0_Feniax[1,:]*100,flutterBoundary_dispVsU_Lambda0_Feniax[2,:], lw=1, ls=:dash, c=colorFeniax, marker=:square, ms=ms, msw=msw, label="Feniax")
            annotate!([15],[ULim-5], text("\$\\Lambda=$(round(Int,Λ*180/π)) ^\\circ\$", 20, :black))
        elseif i==2
            plot!(flutterBoundary_dispVsU_Lambda10_Sharpy[1,:]*100,flutterBoundary_dispVsU_Lambda10_Sharpy[2,:], lw=1, ls=:dash, c=colorSharpy, marker=:circle, ms=ms, msw=msw, label=false)
            plot!(flutterBoundary_dispVsU_Lambda10_Feniax[1,:]*100,flutterBoundary_dispVsU_Lambda10_Feniax[2,:], lw=1, ls=:dash, c=colorFeniax, marker=:square, ms=ms, msw=msw, label=false)
            annotate!([37],[35], text("\$\\Lambda=$(round(Int,Λ*180/π)) ^\\circ\$", 20, :black))
        elseif i==3
            plot!(flutterBoundary_dispVsU_Lambda20_Feniax[1,:]*100,flutterBoundary_dispVsU_Lambda20_Feniax[2,:], lw=1, ls=:dash, c=colorFeniax, marker=:square, ms=ms, msw=msw, label=false)
            annotate!([18],[35], text("\$\\Lambda=$(round(Int,Λ*180/π)) ^\\circ\$", 20, :black))
        elseif i==4
            plot!(flutterBoundary_dispVsU_Lambda30_Feniax[1,:]*100,flutterBoundary_dispVsU_Lambda30_Feniax[2,:], lw=1, ls=:dash, c=colorFeniax, marker=:square, ms=ms, msw=msw, label=false)
            annotate!([13],[35], text("\$\\Lambda=$(round(Int,Λ*180/π)) ^\\circ\$", 20, :black))       
        end
        if i > 1
            plot!(legend=false)
        end
        display(plt2)
        savefig(string(absPath,"/sweptPazyFlutterSweepRangePitchRange_",tipLossType,"_Lambda",Λstr,"_flutterBoundaryDisp.pdf"))
    end
end


# Plot flutter onset and offset speeds vs root pitch angle accross configurations
fill_alpha = 0.75
# Loop angle of sweep
for (i,Λ) in enumerate(ΛRange)

    Λstr = string(round(Int,Λ*180/pi))
    
    # Initialize
    plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[0,120], ylims=[0,7], xticks=collect(0:20:120), yticks=collect(0:1:7), legend=:bottomleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
    if i == 1
        for c in eachindex(tipLossTypeConfig)
            plot!(Shape([NaN],[NaN]), fillcolor=plot_color(configColors[c], fill_alpha), lw=lw, label=configLabels[c])
        end
    end
    annotate!([20],[6], text("\$\\Lambda=$(round(Int,Λ*180/π)) ^\\circ\$", 20, :bottom, :black))

    # Sweep configurations
    for c in eachindex(tipLossTypeConfig)
        # Lines of flutter onset/offset speeds
        x_hump = vcat(speedOnsetHump[c,i],reverse(speedOffsetHump[c,i]))
        y_hump = vcat(θRange*180/π,reverse(θRange*180/π))
        nanIndices = findall(isnan, x_hump)
        x_hump = x_hump[.!isnan.(x_hump)]
        y_hump = y_hump[setdiff(1:length(y_hump), nanIndices)]
        x_hard = vcat(speedOnsetHard[c,i],reverse(speedOffsetHard[c,i]))
        y_hard = vcat(θRange*180/π,reverse(θRange*180/π))
        # Plot
        plot!(Shape(x_hump,y_hump), fillcolor=plot_color(configColors[c], fill_alpha), lw=lw, label=false)
        plot!(Shape(x_hard,y_hard), fillcolor=plot_color(configColors[c], fill_alpha), lw=lw, label=false)
    end
    display(plt1)
    savefig(string(absPath,"/sweptPazyFlutterSweepRangePitchRange_Lambda",Λstr,"_flutterBoundaryPitchAcrossConfigs.pdf"))
end

# Plot flutter onset and offset speeds vs root pitch angle accross sweep angles
fill_alpha = 0.75
for (c,tipLossType) in enumerate(tipLossTypeConfig)
    plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[0,120], ylims=[0,7], xticks=collect(0:20:120), yticks=collect(0:1:7), legend=:bottomleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs,legendtitle=configLabels[c])
    for (j,Λ) in enumerate(reverse(ΛRange))
        i = length(ΛRange)-j+1
        # Lines of flutter onset/offset speeds
        x_hump = vcat(speedOnsetHump[c,i],reverse(speedOffsetHump[c,i]))
        y_hump = vcat(θRange*180/π,reverse(θRange*180/π))
        nanIndices = findall(isnan, x_hump)
        x_hump = x_hump[.!isnan.(x_hump)]
        y_hump = y_hump[setdiff(1:length(y_hump), nanIndices)]
        x_hard = vcat(speedOnsetHard[c,i],reverse(speedOffsetHard[c,i]))
        y_hard = vcat(θRange*180/π,reverse(θRange*180/π))
        # Plot
        plot!(Shape(x_hump,y_hump), fillcolor=plot_color(ΛColors[i], fill_alpha), lw=lw, label="\$\\Lambda=$(round(Int,Λ*180/π)) ^\\circ\$")
        plot!(Shape(x_hard,y_hard), fillcolor=plot_color(ΛColors[i], fill_alpha), lw=lw, label=false)
    end
    display(plt3)
    savefig(string(absPath,"/sweptPazyFlutterSweepRangePitchRange_",tipLossType,"_flutterBoundaryPitchAcrossSweep.pdf"))
end

println("Finished sweptPazyFlutterSweepRangePitchRangePlotGenerator.jl")