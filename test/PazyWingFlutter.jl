using AeroBeams, LinearAlgebra, LinearInterpolations, Plots, ColorSchemes, DelimitedFiles

# Option for mode tracking
modeTracking = true

# Pazy wing
θ = 0
wing,L,nElem,chord,normSparPos,airfoil,surf = create_Pazy(p0=[0;-π/2;θ*π/180])

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
PazyWingFlutter = create_Model(name="PazyWingFlutter",beams=[wing],BCs=[clamp],gravityVector=[0;0;-9.807])

# Set system solver options (limit initial load factor)
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Number of modes
nModes = 5

# Set airspeed range, and initialize outputs
URange = collect(0:0.5:95)
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))
tip_OOP = Array{Float64}(undef,length(URange))

# Sweep airspeed
for (i,U) in enumerate(URange)
    # Display progress
    println("Solving for U = $U m/s")
    # Set tip loss function at current airspeed and root angle
    surf.tipLossDecayFactor = Pazy_tip_loss_factor(θ,U)
    update_beam!(wing)
    # Update velocity of basis A (and update model)
    set_motion_basis_A!(model=PazyWingFlutter,v_A=[0;U;0])
    # Create and solve problem
    problem = create_EigenProblem(model=PazyWingFlutter,nModes=nModes,systemSolver=NR)
    solve!(problem)
    # Frequencies, dampings and eigenvectors
    untrackedFreqs[i] = problem.frequenciesOscillatory
    untrackedDamps[i] = round_off!(problem.dampingsOscillatory,1e-12)
    untrackedEigenvectors[i] = problem.eigenvectorsOscillatoryCplx
    # Get OOP displacement at midchord
    tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
    R,_ = rotation_tensor_WM(tip_p)
    Δ = R*[0; 1; 0]
    tip_twist = asind(Δ[3])
    tip_OOP[i] = -(problem.nodalStatesOverσ[end][nElem].u_n2[1] - chord*(1/2-normSparPos)*sind(tip_twist))
end

# Apply mode tracking, if applicable
if modeTracking
    freqs,damps,_,matchedModes = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
else
    freqs,damps = untrackedFreqs,untrackedDamps
end

# Separate frequencies and damping ratios by mode
modeFrequencies = Array{Vector{Float64}}(undef,nModes)
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeDampingRatios = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeDampingRatios[mode] = modeDampings[mode]./modeFrequencies[mode]
end

# Flutter onset and offset speeds, respective frequencies and OOP displacements
flutterOnsetSpeedsOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOnsetFreqsOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOnsetDispOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOffsetSpeedsOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOffsetFreqsOfMode = Array{Vector{Float64}}(undef,nModes)
flutterOffsetDispOfMode = Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    # Find flutter onset indices
    onsetIndices = findall((modeDampings[mode][2:end] .> 0) .& (modeDampings[mode][1:end-1] .< 0)) .+ 1 
    nIndOn = length(onsetIndices)
    # Loop flutter onset indices
    flutterOnsetSpeeds,flutterOnsetFreqs,flutterOnsetDisp = Vector{Float64}(undef,nIndOn),Vector{Float64}(undef,nIndOn),Vector{Float64}(undef,nIndOn)
    for (n,i) in enumerate(onsetIndices)
        flutterOnsetSpeeds[n] = interpolate(modeDampings[mode][i-1:i],URange[i-1:i],0)
        flutterOnsetFreqs[n] = interpolate(modeDampings[mode][i-1:i],modeFrequencies[mode][i-1:i],0)
        flutterOnsetDisp[n] = interpolate(modeDampings[mode][i-1:i],tip_OOP[i-1:i]/L*100,0)
        println("Flutter onset speed = $(flutterOnsetSpeeds[n]) m/s, flutter frequency = $(flutterOnsetFreqs[n]/(2*π)) Hz, flutter OOP disp = $(flutterOnsetDisp[n]) % semispan (mode $mode)")
    end
    # Set flutter onset variables for current mode
    flutterOnsetSpeedsOfMode[mode] = flutterOnsetSpeeds
    flutterOnsetFreqsOfMode[mode] = flutterOnsetFreqs
    flutterOnsetDispOfMode[mode] = flutterOnsetDisp
    # Find flutter offset indices
    offsetIndices = findall((modeDampings[mode][2:end] .< 0) .& (modeDampings[mode][1:end-1] .> 0)) .+ 1
    nIndOff = length(offsetIndices)
    # Find flutter offset variables
    flutterOffsetSpeeds,flutterOffsetFreqs,flutterOffsetDisp = Vector{Float64}(undef,nIndOff),Vector{Float64}(undef,nIndOff),Vector{Float64}(undef,nIndOff)
    # Loop flutter offset indices
    for (n,i) in enumerate(offsetIndices)
        flutterOffsetSpeeds[n] = interpolate(-modeDampings[mode][i-1:i],URange[i-1:i],0)
        flutterOffsetFreqs[n] = interpolate(-modeDampings[mode][i-1:i],modeFrequencies[mode][i-1:i],0)
        flutterOffsetDisp[n] = interpolate(-modeDampings[mode][i-1:i],tip_OOP[i-1:i]/L*100,0)
        println("Flutter offset speed = $(flutterOffsetSpeeds[n]) m/s, flutter frequency = $(flutterOffsetFreqs[n]/(2*π)) Hz, flutter OOP disp = $(flutterOffsetDisp[n]) % semispan (mode $mode)")
    end
    # Set flutter offset variables for current mode
    flutterOffsetSpeedsOfMode[mode] = flutterOffsetSpeeds
    flutterOffsetFreqsOfMode[mode] = flutterOffsetFreqs
    flutterOffsetDispOfMode[mode] = flutterOffsetDisp
end

# Plots
# ------------------------------------------------------------------------------
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
# V-g-f
plt11 = plot(ylabel="Frequency [Hz]")
for mode in 1:nModes
    plot!(URange, modeFrequencies[mode]/(2*π), c=modeColors[mode], lw=lw,  label=false)
end
plt12 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", ylims=[-0.15,0.05], legend=:bottomleft)
for mode in 1:nModes
    plot!(URange, modeDampingRatios[mode], c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/PazyWingFlutter_1.pdf"))
# Frequencies and dampings vs tip OOP displacement
plt21 = plot(ylabel="Frequency [Hz]")
for mode in 1:nModes
    plot!(tip_OOP/L*100, modeFrequencies[mode]/(2*π), c=modeColors[mode], lw=lw,  label=false)
end
plt22 = plot(xlabel="Tip OOP displacement [% semispan]", ylabel="Damping Ratio", ylims=[-0.15,0.05], legend=:bottomleft)
for mode in 1:nModes
    plot!(tip_OOP/L*100, modeDampingRatios[mode], c=modeColors[mode], lw=lw,  label="Mode $mode")
end
plt2 = plot(plt21,plt22, layout=(2,1))
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/PazyWingFlutter_2.pdf"))

println("Finished PazyWingFlutter.jl")