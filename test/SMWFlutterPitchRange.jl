using AeroBeams, LinearAlgebra, LinearInterpolations, Plots, ColorSchemes, DelimitedFiles

# Wing surface
chord = 1.0
normSparPos = 0.5
aeroSolver = Indicial()
derivationMethod = AD()
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=flatPlate,c=chord,normSparPos=normSparPos)

# Wing beam
L = 16
GJ,EIy,EIz = 1e4,2e4,4e6
ρA,ρIs = 0.75,0.1
ρIy = ρIs*EIy/EIz
ρIz = ρIs-ρIy
nElem = 16
∞ = 1e12
wing = create_Beam(name="wing",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs)],rotationParametrization="E321",aeroSurface=surf)

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
g = 9.80665
h = 20e3
SMWFlutterPitchRange = create_Model(name="SMWFlutterPitchRange",beams=[wing],BCs=[clamp],gravityVector=[0;0;-g],altitude=h)

# Set system solver options (limit initial load factor)
NR1 = create_NewtonRaphson(initialLoadFactor=0.5,maximumLoadFactorStep=0.25)
NR2 = create_NewtonRaphson(initialLoadFactor=0.5,maximumLoadFactorStep=0.5)

# Set number of vibration modes
nModes = 5

# Set root angle and airspeed ranges, and initialize outputs
θRange = collect(0:0.1:5.0)
URange = collect(0:0.1:35)
untrackedFreqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(θRange),length(URange))
tip_u3 = Array{Float64}(undef,length(θRange),length(URange))
flutterOnsetSpeed = Array{Float64}(undef,length(θRange),nModes)
flutterOffsetSpeed = Array{Float64}(undef,length(θRange),nModes)
flutterOnsetFreq = Array{Float64}(undef,length(θRange),nModes)
flutterOffsetFreq = Array{Float64}(undef,length(θRange),nModes)
flutterOnsetTipDisp = Array{Float64}(undef,length(θRange),nModes)
flutterOffsetTipDisp = Array{Float64}(undef,length(θRange),nModes)

# Sweep root angle
for (i,θ) in enumerate(θRange)
    # Update root angle on beam 
    wing.p0=[0;0;θ*π/180]
    update_beam!(wing)
    # Update system solver
    NR = θ <= 1.6 ? NR1 : NR2
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $θ deg, U = $U m/s")
        # Update velocity of basis A (and update model)
        set_motion_basis_A!(model=SMWFlutterPitchRange,v_A=[0;U;0])
        # Create and solve problem
        problem = create_EigenProblem(model=SMWFlutterPitchRange,systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
        solve!(problem)
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = problem.frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(problem.dampingsOscillatory,1e-12)
        untrackedEigenvectors[i,j] = problem.eigenvectorsOscillatoryCplx
        # Tip OOP displacement
        tip_u3[i,j] = problem.nodalStatesOverσ[end][nElem].u_n2[3]
    end
    # Frequencies and dampings after mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Flutter speeds, frequencies and tip displacements of modes at current root angle
    dampsCurrentθ = Array{Vector{Float64}}(undef,nModes)
    freqsCurrentθ = Array{Vector{Float64}}(undef,nModes)
    for mode in 1:nModes
        # Dampings and frequencies of mode
        dampsCurrentθ[mode] = [damps[i,j][mode] for j in eachindex(URange)]
        freqsCurrentθ[mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        # Flutter onset 
        iOnset = findfirst(x->x>0,dampsCurrentθ[mode])
        if isnothing(iOnset) || iOnset == 1
            flutterOnsetSpeed[i,mode] = NaN
            flutterOnsetFreq[i,mode] = NaN
            flutterOnsetTipDisp[i,mode] = NaN
            flutterOffsetSpeed[i,mode] = NaN
            flutterOffsetFreq[i,mode] = NaN
            flutterOffsetTipDisp[i,mode] = NaN
            continue
        end
        flutterOnsetSpeed[i,mode] = interpolate(dampsCurrentθ[mode][iOnset-1:iOnset],URange[iOnset-1:iOnset],0)
        flutterOnsetFreq[i,mode] = interpolate(dampsCurrentθ[mode][iOnset-1:iOnset],freqsCurrentθ[mode][iOnset-1:iOnset],0)
        flutterOnsetTipDisp[i,mode] = interpolate(dampsCurrentθ[mode][iOnset-1:iOnset],tip_u3[i,iOnset-1:iOnset],0)
        # Flutter offset
        iOffset = findnext(x->x<0,dampsCurrentθ[mode],iOnset)
        if isnothing(iOffset)
            flutterOffsetSpeed[i,mode] = NaN
            flutterOffsetFreq[i,mode] = NaN
            flutterOffsetTipDisp[i,mode] = NaN
            continue
        end
        flutterOffsetSpeed[i,mode] = interpolate(-dampsCurrentθ[mode][iOffset-1:iOffset],URange[iOffset-1:iOffset],0)
        flutterOffsetFreq[i,mode] = interpolate(-dampsCurrentθ[mode][iOffset-1:iOffset],freqsCurrentθ[mode][iOffset-1:iOffset],0)
        flutterOffsetTipDisp[i,mode] = interpolate(-dampsCurrentθ[mode][iOffset-1:iOffset],tip_u3[i,iOffset-1:iOffset],0)
    end
end

# Load reference data
flutterSpeedRef = readdlm(string(pwd(),"/test/referenceData/SMW/flutterSpeedVsRootAoA.txt"))
flutterFreqRef = readdlm(string(pwd(),"/test/referenceData/SMW/flutterFreqVsRootAoA.txt"))
flutterTipDispRef = readdlm(string(pwd(),"/test/referenceData/SMW/flutterTipDispVsRootAoA.txt"))

speedVsDispRootAoA0 = readdlm(string(pwd(),"/test/referenceData/SMW/speedVsDispRootAoA0_0.txt"))
speedVsDispRootAoA05 = readdlm(string(pwd(),"/test/referenceData/SMW/speedVsDispRootAoA0_5.txt"))
speedVsDispRootAoA1 = readdlm(string(pwd(),"/test/referenceData/SMW/speedVsDispRootAoA1_0.txt"))
speedVsDispRootAoA2 = readdlm(string(pwd(),"/test/referenceData/SMW/speedVsDispRootAoA2_0.txt"))

freqVsSpeedRootAoA2 = readdlm(string(pwd(),"/test/referenceData/SMW/V-f-RootAoA2.txt"))
dampVsSpeedRootAoA2 = readdlm(string(pwd(),"/test/referenceData/SMW/V-g-RootAoA2.txt"))
freqVsSpeedRootAoA3 = readdlm(string(pwd(),"/test/referenceData/SMW/V-f-RootAoA3.txt"))
dampVsSpeedRootAoA3 = readdlm(string(pwd(),"/test/referenceData/SMW/V-g-RootAoA3.txt"))
freqVsSpeedRootAoA5 = readdlm(string(pwd(),"/test/referenceData/SMW/V-f-RootAoA5.txt"))
dampVsSpeedRootAoA5 = readdlm(string(pwd(),"/test/referenceData/SMW/V-g-RootAoA5.txt"))

# Plots
# ------------------------------------------------------------------------------
modeColors = get(colorschemes[:rainbow], LinRange(0, 1, nModes))
lw = 2
ms = 3
flutterModes = [2,5]
flutterOnsetModesLabels = ["2nd OOP - onset" "2nd T-IP - onset"]
flutterOffsetModesLabels = ["2nd OOP - offset" "2nd T-IP - offset"]
# Flutter speed and frequency vs. root pitch angle
plt11 = plot(ylabel="Flutter speed [m/s]", xlims=[0,5], ylims=[0,35], legend=:bottomleft)
for (i,m) in enumerate(flutterModes)
    plot!(θRange, flutterOnsetSpeed[:,m], c=modeColors[m], ls=:solid, lw=lw, label=flutterOnsetModesLabels[i])
    plot!(θRange, flutterOffsetSpeed[:,m], c=modeColors[m], ls=:dash, lw=lw, label=flutterOffsetModesLabels[i])
end
plot!(flutterSpeedRef[1,:], flutterSpeedRef[2,:], c=:black, ls=:solid, lw=lw, label="Patil et al. (2001)")
plt12 = plot(xlabel="Root pitch angle, \$\\theta\$ [deg]", ylabel="Flutter frequency [rad/s]", xlims=[0,5], ylims=[0,35])
for (i,m) in enumerate(flutterModes)
    plot!(θRange, flutterOnsetFreq[:,m], c=modeColors[m], ls=:solid, lw=lw, label=false)
    plot!(θRange, flutterOffsetFreq[:,m], c=modeColors[m], ls=:dash, lw=lw, label=false)
end
plot!(flutterFreqRef[1,:], flutterFreqRef[2,:], c=:black, ls=:solid, lw=lw, label=false)
plt1 = plot(plt11,plt12, layout=(2,1))
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/SMWFlutterPitchRange_1.pdf"))
# Flutter tip displacement vs. root pitch angle
plt2 = plot(xlabel="Root pitch angle, \$\\theta\$ [deg]",ylabel="Flutter tip displacement [m]", xlims=[0,5], ylims=[-3,3])
for (i,m) in enumerate(flutterModes)
    plot!(θRange, flutterOnsetTipDisp[:,m], c=modeColors[m], ls=:solid, lw=lw, label=false)
    plot!(θRange, flutterOffsetTipDisp[:,m], c=modeColors[m], ls=:dash, lw=lw, label=false)
end
plot!(flutterTipDispRef[1,:], flutterTipDispRef[2,:], c=:black, ls=:solid, lw=lw, label=false)
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/SMWFlutterPitchRange_2.pdf"))
# Tip displacement vs airspeed
θplot = [0; 0.5; 1.0; 2.0]
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(θplot)))
indθ2plot = findall(vec(any(θRange .== θplot', dims=2)))
plt3 = plot(xlabel="Tip displacement [m]",ylabel="Airspeed [m/s]", xlims=[-3,3], ylims=[0,35])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="Patil et al. (2001)")
for i in eachindex(indθ2plot)
    plot!(tip_u3[indθ2plot[i],:], URange, c=colors[i], lw=lw, label="\\theta=$(θplot[i]) deg")
    if θplot[i] == 0.0
        scatter!(speedVsDispRootAoA0[1,:],speedVsDispRootAoA0[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 0.5
        scatter!(speedVsDispRootAoA05[1,:],speedVsDispRootAoA05[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 1.0
        scatter!(speedVsDispRootAoA1[1,:],speedVsDispRootAoA1[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 2.0
        scatter!(speedVsDispRootAoA2[1,:],speedVsDispRootAoA2[2,:], mc=colors[i], ms=ms, msw=0, label=false)    
    end
end
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/SMWFlutterPitchRange_3.pdf"))
# V-g-f of selected mode for selected root angles
mode2plot = 2
θplot = [2.0; 3.0; 5.0]
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(θplot)))
indθ2plot = findall(vec(any(θRange .== θplot', dims=2)))
plt41 = plot(ylabel="Frequency [rad/s]", xlims=[10,30], legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, msw=0, label="Patil et al. (2001)")
for i in eachindex(indθ2plot)
    plot!(URange, [freqs[indθ2plot[i],j][mode2plot] for j in eachindex(URange)], c=colors[i], lw=lw, label=false)
    if θplot[i] == 2.0
        scatter!(freqVsSpeedRootAoA2[1,:],freqVsSpeedRootAoA2[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 3.0
        scatter!(freqVsSpeedRootAoA3[1,:],freqVsSpeedRootAoA3[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 5.0
        scatter!(freqVsSpeedRootAoA5[1,:],freqVsSpeedRootAoA5[2,:], mc=colors[i], ms=ms, msw=0, label=false)   
    end
end
plt42 = plot(xlabel="Airspeed [m/s]", ylabel="Damping [1/s]", xlims=[10,30], ylims=[-1.,0.25], legend=:bottomleft)
for i in eachindex(indθ2plot)
    plot!(URange, [damps[indθ2plot[i],j][mode2plot] for j in eachindex(URange)], c=colors[i], lw=lw, label="\\theta=$(θplot[i]) deg")
    if θplot[i] == 2.0
        scatter!(dampVsSpeedRootAoA2[1,:],dampVsSpeedRootAoA2[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 3.0
        scatter!(dampVsSpeedRootAoA3[1,:],dampVsSpeedRootAoA3[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θplot[i] == 5.0
        scatter!(dampVsSpeedRootAoA5[1,:],dampVsSpeedRootAoA5[2,:], mc=colors[i], ms=ms, msw=0, label=false)   
    end
end
plt4 = plot(plt41,plt42, layout=(2,1))
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/SMWFlutterPitchRange_4.pdf"))

println("Finished SMWFlutterPitchRange.jl")