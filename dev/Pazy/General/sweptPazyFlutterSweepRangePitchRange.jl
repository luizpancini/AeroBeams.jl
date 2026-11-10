using AeroBeams

# --- Configurations ---
# Flag for tip correction
hasTipCorrectionConfig = [false,true,true,true]
# Tip correction function type
tipLossTypeConfig = ["None","Exponential","VLM-undef","VLM-def"]

# Sweep angle [rad]
ΛRange = deg2rad.(vcat(0:10:30))

# Root pitch angle range [rad]
θRange = deg2rad.(unique(sort([vcat(0:0.25:1)...,vcat(1:0.5:7)...])))

# Airspeed range [m/s]
URange = unique(sort([vcat(30:0.5:60)...,vcat(60:1:120)...]))

# Flag for ad hoc corrections on sectional stiffness matrix
sweepStructuralCorrections = false

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil section
airfoil = deepcopy(FlatPlate)

# Flag for upright position
upright = false

# Gravity
g = 0

# Number of modes to track
nModes = 4

# Eigensolver method and minimum number of eigenvalues to be computed
eigenmethod = :krylovkit
nev = 40

# Fixed geometrical and discretization properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Initialize outputs
Nc = length(tipLossTypeConfig)
NΛ = length(ΛRange)
Nθ = length(θRange)
NU = length(URange)
eigensols = Array{EigenSolution}(undef,NU)
freqsOsc = Array{Matrix{Float64}}(undef,Nc,NΛ,Nθ)
dampsOsc = Array{Matrix{Float64}}(undef,Nc,NΛ,Nθ)
flutterOnsetSpeeds = Array{Vector{Float64}}(undef,Nc,NΛ,Nθ)
flutterOnsetFreqs = Array{Vector{Float64}}(undef,Nc,NΛ,Nθ)
flutterOnsetModes = Array{Vector{Int}}(undef,Nc,NΛ,Nθ)
flutterOnsetTipDisps = Array{Vector{Float64}}(undef,Nc,NΛ,Nθ)
flutterOffsetSpeeds = Array{Vector{Float64}}(undef,Nc,NΛ,Nθ)
flutterOffsetFreqs = Array{Vector{Float64}}(undef,Nc,NΛ,Nθ)
flutterOffsetModes = Array{Vector{Int}}(undef,Nc,NΛ,Nθ)
flutterOffsetTipDisps = Array{Vector{Float64}}(undef,Nc,NΛ,Nθ)
tipOOP = Array{Float64}(undef,Nc,NΛ,Nθ,NU)

# Sweep configurations
for (c,hasTipCorrection,tipLossType) in zip(1:length(tipLossTypeConfig),hasTipCorrectionConfig,tipLossTypeConfig)
    # Sweep angle of sweep
    for (i,Λ) in enumerate(ΛRange)
        # Sweep root pitch angle
        for (j,θ) in enumerate(θRange)
            # Initialize guess solution
            x0 = zeros(0)
            # Sweep airspeed
            for (k,U) in enumerate(URange)
                # Display progress
                println("Solving for configuration $c, Λ = $(round(Int,rad2deg(Λ))) deg, θ = $(round(rad2deg(θ),digits=2)) deg, U = $U m/s")
                # Model
                model,_ = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,upright=upright,sweepAngle=Λ,pitchAngle=θ,airspeed=U,gravityVector=[0,0,-g],hasTipCorrection=hasTipCorrection,tipLossType=tipLossType,sweepStructuralCorrections=sweepStructuralCorrections)
                # Create and solve problem
                problem = create_EigenProblem(model=model,eigenmethod=eigenmethod,nev=nev,x0=x0)
                solve!(problem)
                # Eigensolution
                eigensols[k] = problem.eigensolution
                # Get OOP displacement
                tipOOP[c,i,j,k] = problem.loadstepHistory.states[end][nElem].nodal.u2_b[3]
                # Set initial guess solution for next iteration
                x0 = problem.systemArrays.x
            end
            # Apply mode tracking
            freqs,damps,_ = mode_tracking(eigensols)
            # Get frequencies and dampings of oscillatory modes
            freqsOsc[c,i,j],dampsOsc[c,i,j],_ = filter_tracked_modes(freqs,damps; fmin=Hz2rads(1e0), maxTracks=nModes)
            # Flutter speeds, frequencies and tip displacements
            flutterOnsetSpeeds[c,i,j],flutterOnsetFreqs[c,i,j],flutterOnsetModes[c,i,j],flutterOnsetTipDisps[c,i,j],flutterOffsetSpeeds[c,i,j],flutterOffsetFreqs[c,i,j],flutterOffsetModes[c,i,j],flutterOffsetTipDisps[c,i,j] = find_flutter(URange,freqsOsc[c,i,j],dampsOsc[c,i,j],tipOOP[c,i,j,:])
        end
    end
end

# Sharpy and Feniax data (from AePW4 meetings)
flutterBoundary_dispVsU_Lambda0_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda0_Sharpy.txt")
flutterBoundary_dispVsU_Lambda10_Sharpy = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda10_Sharpy.txt")
flutterBoundary_dispVsU_Lambda0_Feniax = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda0_Feniax.txt")
flutterBoundary_dispVsU_Lambda10_Feniax = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda10_Feniax.txt")
flutterBoundary_dispVsU_Lambda20_Feniax = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda20_Feniax.txt")
flutterBoundary_dispVsU_Lambda30_Feniax = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/flutterBoundary_dispVsU_Lambda30_Feniax.txt")

# Test data by Avin et al.
flutterOnsetVelUp = [49, 43, 38]
flutterOnsetVelDown = [55, 48]
flutterOnsetDispUp = [23, 24.5, 26]
flutterOnsetDispDown = [28.5, 31.5]

## Plots
# Set paths
relPath = "/dev/Pazy/General/outputs/sweptPazyFlutterSweepRangePitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
configColors = cgrad([:red,:green,:yellow,:blue], length(tipLossTypeConfig), categorical=true)
ΛColors = cgrad([:red,:green,:yellow,:blue], length(ΛRange), categorical=true)
modeColors = cgrad(:rainbow, nModes, categorical=true)
configLabels = ["No loss" "Exponential loss" "VLM - undeformed" "VLM - deformed"]
ts = 10
fs = 16
lfs = 10
lw = 2
ms = 4
msw = 0
colorFeniax = :purple
colorSharpy = :gold
gr()

# Values (e.g., speeds) at which flutter onset/offset occurs for a given mode
function arrays_of(vals,onsetModes,modes,inds)
    v = similar(θRange)
    for i in eachindex(θRange)
        j = findall(onsetModes[i] .== modes[i])
        v[i] = isempty(j) ? NaN : vals[i][j[inds[i]]]
    end
    return v
end

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

        Λstr = string(round(Int,rad2deg(Λ)))

        # Set sweep-dependent plot configurations
        humpMode = [3 for j in eachindex(θRange)]
        hardMode = [2 for j in eachindex(θRange)]
        humpInd = [1 for j in eachindex(θRange)]
        hardInd = [1 for j in eachindex(θRange)]
        if tipLossType == "None"
            if Λ == 0
                dispLim = 70
                xθ = [26, 32, 35, 37, 39, 40]
                yθ = [72, 67, 60, 55, 49, 43]
                for (j,θ) in enumerate(θRange)
                    if rad2deg(θ) == 0
                        humpMode[j] = 2
                        hardMode[j] = 3
                    elseif rad2deg(θ) > 5
                        humpMode[j] = 0
                    end
                end
            elseif Λ == deg2rad(10)
                dispLim = 40
                xθ = [2.5,  7, 15, 21.5, 33, 37]
                yθ = [ 85, 85, 83,   80, 77, 64]
            elseif Λ == deg2rad(20)
                dispLim = 23.5
                xθ = [ 1,  3, 6.5, 10, 16, 21]
                yθ = [73, 73,  73, 72, 64, 61]
                hardMode = [1 for j in eachindex(θRange)]
            elseif Λ == deg2rad(30)
                dispLim = 15.3
                xθ = [0.7, 1.8,  4, 6.5, 10.5, 14]
                yθ = [ 74,  74, 73,  72,   67, 64]
                humpMode = [4 for j in eachindex(θRange)]
                for (j,θ) in enumerate(θRange)
                    if rad2deg(θ) < 3
                        hardMode[j] = 1
                    end
                end
            end
        elseif tipLossType == "Exponential"
            if Λ == deg2rad(0)
                dispLim = 70
                xθ = [27, 36, 33, 35, 37, 39]
                yθ = [85, 82, 70, 64, 58, 52]
                for (j,θ) in enumerate(θRange)
                    if rad2deg(θ) == 0
                        hardMode[j] = 1
                        humpInd[j] = 1
                    elseif rad2deg(θ) >= 6.5
                        hardMode[j] = 3
                        hardInd[j] = 2
                    end
                end
            elseif Λ == deg2rad(10)
                dispLim = 40
                xθ = [ 2,  6, 12, 18, 23, 25]
                yθ = [80, 79, 75, 73, 60, 51]
                hardMode = [1 for j in eachindex(θRange)]
            elseif Λ == deg2rad(20)
                dispLim = 23.5
                xθ = [ 1,  3,  6,  9, 15, 21]
                yθ = [72, 72, 71, 70, 68, 67]
                hardMode = [1 for j in eachindex(θRange)]
            elseif Λ == deg2rad(30)
                dispLim = 15.3
                xθ = [0.7, 1.8, 3.5,  6, 9.5, 13.5]
                yθ = [ 72,  72,  71, 70,  68,   67]
                hardMode = [1 for j in eachindex(θRange)]
            end
        elseif tipLossType == "VLM-undef"
            if Λ == deg2rad(0)
                dispLim = 70
                xθ = [27, 33, 32, 35, 37, 39]
                yθ = [87, 82, 70, 64, 56, 51]
            elseif Λ == deg2rad(10)
                dispLim = 40
                xθ = [ 2,  6, 12, 19, 27, 33]
                yθ = [85, 84, 83, 82, 75, 67]
                for (j,θ) in enumerate(θRange)
                    if rad2deg(θ) < 7
                        hardMode[j] = 1
                    end
                end
            elseif Λ == deg2rad(20)
                dispLim = 23.5
                xθ = [ 1,  3,  6, 9.5, 15, 20]
                yθ = [80, 80, 79,  78, 75, 72]
                hardMode = [1 for j in eachindex(θRange)]
            elseif Λ == deg2rad(30)
                dispLim = 15
                xθ = [0.7, 1.8,  3.8,  6, 10, 13.5]
                yθ = [ 85,  85, 84, 83, 81,   78]
                hardMode = [1 for j in eachindex(θRange)]
            end    
        elseif tipLossType == "VLM-def"
            if Λ == deg2rad(0)
                dispLim = 70
                xθ = [30, 35, 33, 36, 37, 38]
                yθ = [90, 83, 70, 63, 57, 52]
            elseif Λ == deg2rad(10)
                dispLim = 40
                xθ = [ 2,  6, 12, 18, 23, 25]
                yθ = [85, 84, 80, 78, 65, 55]
                hardMode = [1 for j in eachindex(θRange)]
            elseif Λ == deg2rad(20)
                dispLim = 23.5
                xθ = [ 1,  3,  6, 9.5, 15, 20]
                yθ = [85, 85, 84,  82, 80, 75]
                hardMode = [1 for j in eachindex(θRange)]
            elseif Λ == deg2rad(30)
                dispLim = 15
                xθ = [0.7, 1.8, 3.8,  6, 10, 13.5]
                yθ = [ 90,  90,  89, 87, 85,   82]
                hardMode = [1 for j in eachindex(θRange)]
            end    
        end

        # Flutter onset/offset speeds and displacements
        speedOnsetHump[c,i] = arrays_of(flutterOnsetSpeeds[c,i,:],flutterOnsetModes[c,i,:],humpMode,humpInd)
        speedOffsetHump[c,i] = arrays_of(flutterOffsetSpeeds[c,i,:],flutterOffsetModes[c,i,:],humpMode,humpInd)
        speedOnsetHard[c,i] = arrays_of(flutterOnsetSpeeds[c,i,:],flutterOnsetModes[c,i,:],hardMode,hardInd)
        speedOffsetHard[c,i] = [150 for j in eachindex(θRange)]
        dispOnsetHump[c,i] = arrays_of(100/L*flutterOnsetTipDisps[c,i,:],flutterOnsetModes[c,i,:],humpMode,humpInd)
        dispOffsetHump[c,i] = arrays_of(100/L*flutterOffsetTipDisps[c,i,:],flutterOffsetModes[c,i,:],humpMode,humpInd)
        dispOnsetHard[c,i] = arrays_of(100/L*flutterOnsetTipDisps[c,i,:],flutterOnsetModes[c,i,:],hardMode,hardInd)
        dispOffsetHard[c,i] = dispOnsetHard[c,i]

        if tipLossType == "None"
        elseif tipLossType == "Exponential" 
            if Λ == deg2rad(30)
                speedOffsetHump[c,i][end] = 150
            end
        end

        # V-g-f
        θ2plot = deg2rad.([0, 1, 3, 5, 7])
        plt_Vf = plot(ylabel="Frequency [Hz]", xlims=extrema(URange), ylims=[0,50], tickfont=font(ts), guidefont=font(fs), title="Configuration $c - $(configLabels[c])")
        annotate!([110],[45], text("\$\\Lambda=$(round(Int,rad2deg(Λ))) ^\\circ\$", 16, :black))   
        for (j,θ) in enumerate(θRange)
            if !any(y -> isapprox(θ, y; atol=1e-4), θ2plot)
                continue
            end
            for mode in 1:nModes
                plot!(URange, rads2Hz.(freqsOsc[c,i,j][mode,:]), c=modeColors[mode], lw=lw, alpha=0.15+0.85*j/length(θRange), label=false)
            end
        end
        plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=extrema(URange), ylims=[-0.25,0.15], tickfont=font(ts), guidefont=font(fs), legend=:topleft)
        for (j,θ) in enumerate(θRange)
            if !any(y -> isapprox(θ, y; atol=1e-4), θ2plot)
                continue
            end
            for mode in 1:nModes
                plot!(URange, dampsOsc[c,i,j][mode,:]./freqsOsc[c,i,j][mode,:], c=modeColors[mode], lw=lw, alpha=0.15+0.85*j/length(θRange), label=false)
            end
        end
        plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
        display(plt_Vgf)
        savefig(string(absPath,"/sweptPazyFlutterSweepRangePitchRange_Vgf_",tipLossType,"_Lambda",Λstr,".pdf"))

        # Lines of flutter onset/offset speeds
        x_hump = vcat(speedOnsetHump[c,i],reverse(speedOffsetHump[c,i]))
        y_hump = vcat(rad2deg.(θRange),reverse(rad2deg.(θRange)))
        nanIndices = findall(isnan, x_hump)
        x_hump = x_hump[.!isnan.(x_hump)]
        y_hump = y_hump[setdiff(1:length(y_hump), nanIndices)]
        x_hard = vcat(speedOnsetHard[c,i],reverse(speedOffsetHard[c,i]))
        y_hard = vcat(rad2deg.(θRange),reverse(rad2deg.(θRange)))

        # Flutter onset and offset speeds vs root pitch angle
        plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[0,URange[end]], ylims=[0,7], xticks=collect(0:20:URange[end]), yticks=collect(0:1:7), legend=:bottomleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, title="Configuration $c - $(configLabels[c])")
        annotate!([15],[2.5], text("\$\\Lambda=$(round(Int,rad2deg(Λ))) ^\\circ\$", 16, :black)) 
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
        θ2plot = deg2rad.([0.5, 1, 2, 3, 5, 7])
        indθ2plot = findall(x -> any(isapprox(x, θ; atol=1e-10) for θ in θ2plot), θRange)
        θcolors = cgrad([:blue, :red], length(θ2plot), categorical=true)
        plt2 = plot(xlabel="Tip OOP displacement [% semispan]", ylabel="Airspeed [m/s]", xlims=[0,dispLim], ylims=extrema(URange), xticks=collect(0:5:dispLim), yticks=collect(30:10:URange[end]), legend=:bottomright, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, title="Configuration $c - $(configLabels[c])")
        plot!(Shape(x_hump,y_hump), fillcolor=plot_color(:green, 0.25), lw=lw, label="Hump flutter region")
        plot!(Shape(x_hard,y_hard), fillcolor=plot_color(:red, 0.25), lw=lw, label="Hard flutter region")
        for (n,ind) in enumerate(indθ2plot)
            if Λ == 0
                if tipLossType == "None"
                    U90ind = findfirst(x-> x > 90, URange) - 1
                    range = n == 1 ? (1:U90ind) : (1:length(URange))
                else
                    # U110ind = findfirst(x-> x > 110, URange) - 1
                    range = 1:length(URange)
                end
            else
                range = 1:length(URange)
            end    
            plot!(tipOOP[c,i,ind,range]/L*100, URange[range], c=θcolors[n], ls=:dash, lw=lw, label=false)
            if i == 1
                annotate!([xθ[n]],[yθ[n]], text("\$$(round(rad2deg(θ2plot[n]),digits=2)) ^\\circ\$", 10, :bottom, θcolors[n]))
            end
        end
        if i==1
            plot!(flutterBoundary_dispVsU_Lambda0_Sharpy[1,:]*100,flutterBoundary_dispVsU_Lambda0_Sharpy[2,:], lw=1, ls=:dash, c=colorSharpy, marker=:circle, ms=ms, msw=msw, label="Sharpy")
            plot!(flutterBoundary_dispVsU_Lambda0_Feniax[1,:]*100,flutterBoundary_dispVsU_Lambda0_Feniax[2,:], lw=1, ls=:dash, c=colorFeniax, marker=:square, ms=ms, msw=msw, label="Feniax")
            annotate!([15],[URange[end]-5], text("\$\\Lambda=$(round(Int,rad2deg(Λ))) ^\\circ\$", 20, :black))
        elseif i==2
            plot!(flutterBoundary_dispVsU_Lambda10_Sharpy[1,:]*100,flutterBoundary_dispVsU_Lambda10_Sharpy[2,:], lw=1, ls=:dash, c=colorSharpy, marker=:circle, ms=ms, msw=msw, label=false)
            plot!(flutterBoundary_dispVsU_Lambda10_Feniax[1,:]*100,flutterBoundary_dispVsU_Lambda10_Feniax[2,:], lw=1, ls=:dash, c=colorFeniax, marker=:square, ms=ms, msw=msw, label=false)
            annotate!([37],[35], text("\$\\Lambda=$(round(Int,rad2deg(Λ))) ^\\circ\$", 20, :black))
        elseif i==3
            plot!(flutterBoundary_dispVsU_Lambda20_Feniax[1,:]*100,flutterBoundary_dispVsU_Lambda20_Feniax[2,:], lw=1, ls=:dash, c=colorFeniax, marker=:square, ms=ms, msw=msw, label=false)
            annotate!([18],[35], text("\$\\Lambda=$(round(Int,rad2deg(Λ))) ^\\circ\$", 20, :black))
        elseif i==4
            plot!(flutterBoundary_dispVsU_Lambda30_Feniax[1,:]*100,flutterBoundary_dispVsU_Lambda30_Feniax[2,:], lw=1, ls=:dash, c=colorFeniax, marker=:square, ms=ms, msw=msw, label=false)
            annotate!([13],[35], text("\$\\Lambda=$(round(Int,rad2deg(Λ))) ^\\circ\$", 20, :black))       
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

    Λstr = string(round(Int,rad2deg(Λ)))
    
    # Initialize
    plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[0,121], ylims=[0,7], xticks=collect(0:20:120), yticks=collect(0:1:7), legend=:bottomleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
    if i == 1
        for c in eachindex(tipLossTypeConfig)
            plot!(Shape([NaN],[NaN]), fillcolor=plot_color(configColors[c], fill_alpha), lw=lw, label=configLabels[c])
        end
    end
    annotate!([20],[6], text("\$\\Lambda=$(round(Int,rad2deg(Λ))) ^\\circ\$", 20, :bottom, :black))

    # Sweep configurations
    for c in eachindex(tipLossTypeConfig)
        # Lines of flutter onset/offset speeds
        x_hump = vcat(speedOnsetHump[c,i],reverse(speedOffsetHump[c,i]))
        y_hump = vcat(rad2deg.(θRange),reverse(rad2deg.(θRange)))
        nanIndices = findall(isnan, x_hump)
        x_hump = x_hump[.!isnan.(x_hump)]
        y_hump = y_hump[setdiff(1:length(y_hump), nanIndices)]
        x_hard = vcat(speedOnsetHard[c,i],reverse(speedOffsetHard[c,i]))
        y_hard = vcat(rad2deg.(θRange),reverse(rad2deg.(θRange)))
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
    plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[0,121], ylims=[0,7], xticks=collect(0:20:120), yticks=collect(0:1:7), legend=:bottomleft, tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs,legendtitle=configLabels[c])
    for (j,Λ) in enumerate(reverse(ΛRange))
        i = length(ΛRange)-j+1
        # Lines of flutter onset/offset speeds
        x_hump = vcat(speedOnsetHump[c,i],reverse(speedOffsetHump[c,i]))
        y_hump = vcat(rad2deg.(θRange),reverse(rad2deg.(θRange)))
        nanIndices = findall(isnan, x_hump)
        x_hump = x_hump[.!isnan.(x_hump)]
        y_hump = y_hump[setdiff(1:length(y_hump), nanIndices)]
        x_hard = vcat(speedOnsetHard[c,i],reverse(speedOffsetHard[c,i]))
        y_hard = vcat(rad2deg.(θRange),reverse(rad2deg.(θRange)))
        # Plot
        plot!(Shape(x_hump,y_hump), fillcolor=plot_color(ΛColors[i], fill_alpha), lw=lw, label="\$\\Lambda=$(round(Int,rad2deg(Λ))) ^\\circ\$")
        plot!(Shape(x_hard,y_hard), fillcolor=plot_color(ΛColors[i], fill_alpha), lw=lw, label=false)
    end
    display(plt3)
    savefig(string(absPath,"/sweptPazyFlutterSweepRangePitchRange_",tipLossType,"_flutterBoundaryPitchAcrossSweep.pdf"))
end

println("Finished sweptPazyFlutterSweepRangePitchRange.jl")