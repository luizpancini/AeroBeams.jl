using AeroBeams, LinearAlgebra, LinearInterpolations, DelimitedFiles

# Option for mode tracking
modeTracking = true

# Aerodynamic solver and derivatives method
aeroSolver = Indicial()
derivationMethod = AD()

# Altitude
h = 20e3

# Gravity
g = 0

# Discretization
nElem = 20

# Pitch angle
θ = 0

# System solver
σ0 = 1
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Set precurvature, tip force and airspeed ranges, and initialize outputs
kRange = collect(0:0.018:0.018)
F3Range = collect(0:2:40)
URange = vcat(1e-3,collect(20:0.5:35))
untrackedFreqs = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(kRange),length(F3Range),length(URange))
freqs = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),length(URange))
damps = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),length(URange))
tip_u3 = Array{Float64}(undef,length(kRange),length(F3Range),length(URange))
flutterSpeed = Array{Float64}(undef,length(kRange),length(F3Range))
flutterFreq = Array{Float64}(undef,length(kRange),length(F3Range))
flutterMode = Array{Int64}(undef,length(kRange),length(F3Range))
flutterTipDisp = Array{Float64}(undef,length(kRange),length(F3Range))
modeDampings = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),nModes)
modeFrequencies = Array{Vector{Float64}}(undef,length(kRange),length(F3Range),nModes)

model = Array{Model}(undef,length(kRange),length(F3Range),length(URange))

# Set number of vibration modes
nModes = 5

# Sweep wing precurvature
for (ki,k) in enumerate(kRange)
    # Sweep tip force
    for (i,F3) in enumerate(F3Range)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for k=$k, F3 = $F3 N, U = $U m/s")
            # Update model
            model[ki,i,j],_ = create_SMW(aeroSolver=aeroSolver,derivationMethod=derivationMethod,θ=θ*π/180,k2=k,nElem=nElem,altitude=h,g=g,tipF3=F3,airspeed=U)
            # Create and solve problem
            problem = create_EigenProblem(model=model[ki,i,j],systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
            solve!(problem)
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[ki,i,j] = problem.frequenciesOscillatory
            untrackedDamps[ki,i,j] = round_off!(problem.dampingsOscillatory,1e-8)
            untrackedEigenvectors[ki,i,j] = problem.eigenvectorsOscillatoryCplx
            # Tip OOP displacement
            tip_u3[ki,i,j] = problem.nodalStatesOverσ[end][nElem].u_n2[3]
        end
        # Frequencies and dampings after mode tracking
        if modeTracking
            freqs[ki,i,:],damps[ki,i,:],_ = mode_tracking_hungarian(URange,untrackedFreqs[ki,i,:],untrackedDamps[ki,i,:],untrackedEigenvectors[ki,i,:])
        else
            freqs[ki,i,:],damps[ki,i,:] = untrackedFreqs[ki,i,:],untrackedDamps[ki,i,:]
        end
        # Separate frequencies and dampings by mode
        for mode in 1:nModes
            modeFrequencies[ki,i,mode] = [freqs[ki,i,j][mode] for j in eachindex(URange)]
            modeDampings[ki,i,mode] = [damps[ki,i,j][mode] for j in eachindex(URange)]
        end
        # Flutter speeds, frequencies and tip displacements of modes at current tip force
        dampsCurrentF3 = Array{Vector{Float64}}(undef,nModes)
        freqsCurrentF3 = Array{Vector{Float64}}(undef,nModes)
        flutterSpeedOfMode = Array{Float64}(undef,nModes)
        flutterFreqOfMode = Array{Float64}(undef,nModes)
        flutterTipDispOfMode = Array{Float64}(undef,nModes)
        for mode in 1:nModes
            dampsCurrentF3[mode] = [damps[ki,i,j][mode] for j in eachindex(URange)]
            freqsCurrentF3[mode] = [freqs[ki,i,j][mode] for j in eachindex(URange)]
            indexInstability = findfirst(x->x>0,dampsCurrentF3[mode])
            if isnothing(indexInstability)
                flutterSpeedOfMode[mode] = NaN
                flutterFreqOfMode[mode] = NaN
                flutterTipDispOfMode[mode] = NaN
                continue
            end
            flutterSpeedOfMode[mode] = interpolate(dampsCurrentF3[mode][indexInstability-1:indexInstability],URange[indexInstability-1:indexInstability],0)
            flutterFreqOfMode[mode] = interpolate(dampsCurrentF3[mode][indexInstability-1:indexInstability],freqsCurrentF3[mode][indexInstability-1:indexInstability],0)
            flutterTipDispOfMode[mode] = interpolate(dampsCurrentF3[mode][indexInstability-1:indexInstability],tip_u3[ki,i,indexInstability-1:indexInstability],0)
        end
        # Set flutter speed as the greatest (for compatibility with reference), and get corresponding flutter mode and frequency, and tip displacement
        flutterSpeed[ki,i] = maximum(filter(!isnan,flutterSpeedOfMode))
        flutterMode[ki,i] = findfirst(x->x==flutterSpeed[ki,i],flutterSpeedOfMode)
        flutterFreq[ki,i] = flutterFreqOfMode[flutterMode[ki,i]]
        flutterTipDisp[ki,i] = flutterTipDispOfMode[flutterMode[ki,i]]
    end
end

# Load reference data
flutterSpeedVsTipLoadk0 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/SMW/flutterSpeedVsTipLoadk0.txt")
flutterSpeedVsTipLoadk2 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/SMW/flutterSpeedVsTipLoadk2.txt")
flutterSpeedVsTipDispk0 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/SMW/flutterSpeedVsTipDisp.txt")
freqsVsTipOOP = readdlm(pkgdir(AeroBeams)*"/test/referenceData/SMW/freqsVsTipOOP.txt")

# Adjust data for plots by padding matrices with NaN
matrices = Dict(:freqsVsTipOOP => freqsVsTipOOP)
for key in keys(matrices)
    matrices[key] .= [x == "" ? NaN64 : x for x in matrices[key]]
end
freqsVsTipOOP = matrices[:freqsVsTipOOP]

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_flutter_compare_Patil"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = cgrad(:rainbow, length(kRange), categorical=true)
modeColors = cgrad(:rainbow, nModes, categorical=true)
mshape = [:utriangle, :utriangle, :circle, :circle, :utriangle]
ts = 10
fs = 16
lfs = 12
lw = 2
ms = 4
msw = 0
L = 16
gr()

# Flutter speed vs. tip load
plt1 = plot(xlabel="Tip load [N]", ylabel="Flutter speed [m/s]", xlims=[0,40], ylims=[0,40], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend_position=:bottomright)
for (i,k) in enumerate(kRange)
    plot!(F3Range, flutterSpeed[i,:], c=colors[i], lw=lw, label="AeroBeams - \$k_2=$k\$")
end
plot!(flutterSpeedVsTipLoadk0[1,:], flutterSpeedVsTipLoadk0[2,:], c=:black, ls=:dash, lw=lw, label="Patil et al. (2001)")
plot!(flutterSpeedVsTipLoadk2[1,:], flutterSpeedVsTipLoadk2[2,:], c=:black, ls=:dash, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/cHALEwing_flutter_compare_Patil_flutterSpeedVsTipLoad.pdf"))

# Flutter speed vs. tip disp
plt2 = plot(xlabel="Tip vertical position [% semispan]", ylabel="Flutter speed [m/s]", xlims=[-20,20], ylims=[0,40], tickfont=font(ts), guidefont=font(fs))
for (i,k) in enumerate(kRange)
    plot!((flutterTipDisp[i,:] .+ model[i,1,1].elements[end].r_n2[3])/L*100, flutterSpeed[i,:], c=colors[i], lw=lw, label=false)
end
# plot!(flutterSpeedVsTipDispk0[1,:]/L*100, flutterSpeedVsTipDispk0[2,:], c=:black, ls=:dash, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/cHALEwing_flutter_compare_Patil_flutterSpeedVsTipDisp.pdf"))

# Frequency vs tip displacement of initially straight wing (at lowest airspeed)
plt3 = plot(xlabel = "Tip OOP displacement [m]", ylabel="Frequency [rad/s]", xlims=[0,2.5], ylims=[0,50], tickfont=font(ts), guidefont=font(fs), legendfontsize=10, legend_position=(0.65,0.65))
scatter!([NaN], [NaN], c=:black, shape=mshape[1], ms=ms, msw=msw, label="OOP bending modes")
scatter!([NaN], [NaN], c=:black, shape=mshape[3], ms=ms, msw=msw, label="T-IP bending modes")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dash, label="Patil et al. (2001)")
for (i,F3) in enumerate(F3Range)
    for mode in 1:nModes
        # Fix mode swap
        if mode == 4 && i > 11
            colorNow = modeColors[5]
            shapeNow = mshape[5]
        elseif mode == 5 && i > 11
            colorNow = modeColors[4]
            shapeNow = mshape[4]
        else
            colorNow = modeColors[mode] 
            shapeNow = mshape[mode]           
        end
        # Plot
        scatter!(tip_u3[1,i,:], [modeFrequencies[1,i,mode][1]], c=colorNow, shape=shapeNow, ms=ms, msw=msw, label=false)
        plot!(freqsVsTipOOP[2*mode-1,:], freqsVsTipOOP[2*mode,:], c=modeColors[mode], lw=lw, ls=:dash, label=false)
    end
end
display(plt3)
savefig(string(absPath,"/cHALEwing_flutter_compare_Patil_modalFreqVsDisp.pdf"))

println("Finished cHALEwing_flutter_compare_Patil.jl")