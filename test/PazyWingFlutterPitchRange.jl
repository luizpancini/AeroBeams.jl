using AeroBeams, LinearAlgebra, LinearInterpolations, Plots, ColorSchemes, DelimitedFiles

# Option for mode tracking
modeTracking = true

# Pazy wing
wing,L,nElem,chord,normSparPos,airfoil,surf = create_Pazy(p0=[0;-π/2;0],airfoil=flatPlate)

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
PazyWingFlutterPitchRange = create_Model(name="PazyWingFlutterPitchRange",beams=[wing],BCs=[clamp],gravityVector=[0;0;-9.807])

# Set system solver options (limit initial load factor)
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Number of modes
nModes = 5

# Set pitch angle and airspeed ranges, and initialize outputs
θRange = collect(-0.25:0.25:7)
URange = collect(0:0.5:90)
untrackedFreqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
tip_OOP = Array{Float64}(undef,length(θRange),length(URange))
flutterOnsetSpeedsOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)
flutterOnsetFreqsOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)
flutterOnsetDispOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)
flutterOffsetSpeedsOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)
flutterOffsetFreqsOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)
flutterOffsetDispOfMode = Array{Vector{Float64}}(undef,length(θRange),nModes)

# Sweep root angle
for (i,θ) in enumerate(θRange)
    # Update root angle on beam 
    wing.p0[3]=θ*π/180
    update_beam!(wing)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $θ deg, U = $U m/s")
        # Set tip loss function at current airspeed and root angle
        surf.tipLossDecayFactor = Pazy_tip_loss_factor(θ,U)
        update_beam!(wing)
        # Update velocity of basis A (and update model)
        set_motion_basis_A!(model=PazyWingFlutterPitchRange,v_A=[0;U;0])
        # Create and solve problem
        problem = create_EigenProblem(model=PazyWingFlutterPitchRange,nModes=nModes,systemSolver=NR)
        solve!(problem)
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = problem.frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(problem.dampingsOscillatory,1e-12)
        untrackedEigenvectors[i,j] = problem.eigenvectorsOscillatoryCplx
        # Get OOP displacement at midchord
        tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
        R,_ = rotation_tensor_WM(tip_p)
        Δ = R*[0; 1; 0]
        tip_twist = asind(Δ[3])
        tip_OOP[i,j] = -(problem.nodalStatesOverσ[end][nElem].u_n2[1] - chord*(1/2-normSparPos)*sind(tip_twist))
    end
    # Apply mode tracking, if applicable
    if modeTracking
        freqs[i,:],damps[i,:],_ = mode_tracking(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    else
        freqs[i,:],damps[i,:] = untrackedFreqs[i,:],untrackedDamps[i,:]
    end
    # Separate frequencies and damping ratios by mode
    modeFrequencies = Array{Vector{Float64}}(undef,nModes)
    modeDampings = Array{Vector{Float64}}(undef,nModes)
    modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
    for mode in 1:nModes
        modeFrequencies[mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[mode] = [damps[i,j][mode] for j in eachindex(URange)]
        modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
    end
    # Loop over modes: compute flutter onset and offset speeds, respective frequencies and OOP displacements
    for mode in 1:nModes
        # Find flutter onset indices
        onsetIndices = findall((modeDampings[mode][2:end] .> 0) .& (modeDampings[mode][1:end-1] .< 0)) .+ 1 
        nIndOn = length(onsetIndices)
        # Loop flutter onset indices
        flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetDisp = Vector{Float64}(undef,nIndOn),Vector{Float64}(undef,nIndOn),Vector{Float64}(undef,nIndOn)
        for (n,k) in enumerate(onsetIndices)
            flutterOnsetSpeeds[n] = interpolate(modeDampings[mode][k-1:k],URange[k-1:k],0)
            flutterOnsetFreqs[n] = interpolate(modeDampings[mode][k-1:k],modeFrequencies[mode][k-1:k],0)
            flutterOnsetDisp[n] = interpolate(modeDampings[mode][k-1:k],tip_OOP[i,k-1:k]/L*100,0)
        end
        if nIndOn == 0
            flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetDisp = [NaN],[NaN],[NaN]
        end
        # Set flutter onset variables for current mode
        flutterOnsetSpeedsOfMode[i,mode] = flutterOnsetSpeeds
        flutterOnsetFreqsOfMode[i,mode] = flutterOnsetFreqs
        flutterOnsetDispOfMode[i,mode] = flutterOnsetDisp
        # Find flutter offset indices
        offsetIndices = findall((modeDampings[mode][2:end] .< 0) .& (modeDampings[mode][1:end-1] .> 0)) .+ 1
        nIndOff = length(offsetIndices)
        # Find flutter offset variables
        flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetDisp = Vector{Float64}(undef,nIndOff),Vector{Float64}(undef,nIndOff),Vector{Float64}(undef,nIndOff)
        # Loop flutter offset indices
        for (n,k) in enumerate(offsetIndices)
            flutterOffsetSpeeds[n] = interpolate(-modeDampings[mode][k-1:k],URange[k-1:k],0)
            flutterOffsetFreqs[n] = interpolate(-modeDampings[mode][k-1:k],modeFrequencies[mode][k-1:k],0)
            flutterOffsetDisp[n] = interpolate(-modeDampings[mode][k-1:k],tip_OOP[i,k-1:k]/L*100,0)
        end
        if nIndOff == 0
            flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetDisp = [NaN],[NaN],[NaN]
        end
        # Set flutter offset variables for current mode
        flutterOffsetSpeedsOfMode[i,mode] = flutterOffsetSpeeds
        flutterOffsetFreqsOfMode[i,mode] = flutterOffsetFreqs
        flutterOffsetDispOfMode[i,mode] = flutterOffsetDisp
    end
end

# Load reference data from DRACHINSKI et al.: Flutter Tests of the Pazy Wing(2021)
rootPitchVelUp = [3; 5; 7]
rootPitchVelDown = [3; 5]
flutterOnsetVelUp = [49; 43; 38]
flutterOffsetVelUp = [58; 51; 46]
flutterOnsetVelDown = [55; 48]
flutterOffsetVelDown = [40; 36]
flutterOnsetDispUp = [23; 24.5; 26]
flutterOnsetDispDown = [28.5; 31.5]

# Plots
# ------------------------------------------------------------------------------
lw = 2
ms = 10
# Flutter onset and offset speeds vs root pitch angle
mode2plot = 3
x1 = [flutterOnsetSpeedsOfMode[i,mode2plot][1] for i in eachindex(θRange)]
x2 = [flutterOffsetSpeedsOfMode[i,mode2plot][1] for i in eachindex(θRange)]
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Root pitch angle [deg]", xlims=[0,90], ylims=[0,7.25], xticks=collect(0:15:90), yticks=collect(0:1:7), legend=:bottomleft)
plot!(Shape(vcat(x1,reverse(x2[3:end]),97,100),vcat(θRange,reverse(θRange))), fillcolor = plot_color(:red, 0.25), lw=lw, label="AeroBeams flutter region")
scatter!(flutterOnsetVelUp, rootPitchVelUp, shape=:rtriangle, mc=:red, ms=ms, msw=0, label="Test onset up")
scatter!(flutterOffsetVelUp, rootPitchVelUp, shape=:rtriangle, mc=:green, ms=ms, msw=0, label="Test offset up")
scatter!(flutterOnsetVelDown, rootPitchVelDown, shape=:ltriangle, mc=:red, ms=ms, msw=0, label="Test onset down")
scatter!(flutterOffsetVelDown, rootPitchVelDown, shape=:ltriangle, mc=:green, ms=ms, msw=0, label="Test offset down")
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/PazyWingFlutterPitchRange_1.pdf"))
# Flutter onset and offset speeds vs tip OOP displacement for varying root pitch angle
θ2plot = [0.5,1,2,3,5,7]
indθ2plot = findall(vec(any(θRange .== θ2plot', dims=2)))
x1 = [flutterOnsetDispOfMode[i,mode2plot][1] for i in eachindex(θRange)]
x2 = [flutterOffsetDispOfMode[i,mode2plot][1] for i in eachindex(θRange)]
y1 = [flutterOnsetSpeedsOfMode[i,mode2plot][1] for i in eachindex(θRange)]
y2 = [flutterOffsetSpeedsOfMode[i,mode2plot][1] for i in eachindex(θRange)]
plt2 = plot(xlabel="Tip OOP displacement [% semispan]", ylabel="Airspeed [m/s]", xlims=[0,32], ylims=[30,90], xticks=collect(0:5:30), yticks=collect(30:10:90), legend=:bottomleft)
plot!(Shape(vcat(x1,reverse(x2[3:end]),22,22,0),vcat(y1,reverse(y2[3:end]),97,100,90)), fillcolor = plot_color(:red, 0.25), lw=lw, label="AeroBeams flutter region")
for (n,ind) in enumerate(indθ2plot)
    θ = θ2plot[n]
    plot!(tip_OOP[ind,:]/L*100, URange, c=:black, ls=:dash, lw=lw, label=false)
    if θ==0.5
        xind,yind = 16,79
    elseif θ==1
        xind,yind = 20,72
    elseif θ==2
        xind,yind = 22.5,63
    elseif θ==3
        xind,yind = 24,56
    elseif θ==5
        xind,yind = 26,49
    elseif θ==7
        xind,yind = 27,42.5  
    end
    annotate!([xind],[yind], text("$(θ2plot[n]) deg", 10, :bottom))
end
scatter!(flutterOnsetDispUp, flutterOnsetVelUp, shape=:rtriangle, mc=:red, ms=ms, msw=0, label="Test onset up")
scatter!(flutterOnsetDispDown, flutterOnsetVelDown, shape=:ltriangle, mc=:red, ms=ms, msw=0, label="Test onset down")
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/PazyWingFlutterPitchRange_2.pdf"))

println("Finished PazyWingFlutterPitchRange.jl")