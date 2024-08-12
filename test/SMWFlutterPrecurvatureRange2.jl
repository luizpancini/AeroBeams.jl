using AeroBeams, LinearAlgebra, LinearInterpolations, Plots, ColorSchemes, DelimitedFiles

# Wing surface
airfoil = deepcopy(flatPlate)
chord = 1.0
normSparPos = 0.5
aeroSolver = Indicial()
derivationMethod = AD()
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=chord,normSparPos=normSparPos)

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
SMWFlutterPrecurvatureRange2 = create_Model(name="SMWFlutterPrecurvatureRange2",beams=[wing],BCs=[clamp],gravityVector=[0;0;-g],altitude=h)

# Set system solver options (limit initial load factor)
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep,displayStatus=false)

# Set number of vibration modes
nModes = 5

# Set precurvature, root angle, and airspeed ranges, and initialize outputs
kRange = collect(0:0.018:0.018)
θRange = collect(0:0.2:5)
URange = collect(0:0.5:33)
untrackedFreqs = Array{Vector{Float64}}(undef,length(kRange),length(θRange),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(kRange),length(θRange),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(kRange),length(θRange),length(URange))
freqs = Array{Vector{Float64}}(undef,length(kRange),length(θRange),length(URange))
damps = Array{Vector{Float64}}(undef,length(kRange),length(θRange),length(URange))
tip_u3 = Array{Float64}(undef,length(kRange),length(θRange),length(URange))
flutterOnsetSpeed = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]
flutterOffsetSpeed = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]
flutterOnsetFreq = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]
flutterOffsetFreq = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]
flutterOnsetTipDisp = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]
flutterOffsetTipDisp = [[[Float64[] for _ in 1:nModes] for _ in θRange] for _ in kRange]

# Sweep wing precurvature
for (ki,k) in enumerate(kRange)
    # Update precurvature on beam
    wing.k[2] = k
    # Sweep root angle
    for (i,θ) in enumerate(θRange)
        # Update tip force on model
        wing.p0[3] = θ*π/180
        update_beam!(wing)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for k=$k, θ = $θ deg, U = $U m/s")
            # Update velocity of basis A (and update model)
            set_motion_basis_A!(model=SMWFlutterPrecurvatureRange2,v_A=[0;U;0])
            # Create and solve problem
            problem = create_EigenProblem(model=SMWFlutterPrecurvatureRange2,systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
            solve!(problem)
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[ki,i,j] = problem.frequenciesOscillatory
            untrackedDamps[ki,i,j] = round_off!(problem.dampingsOscillatory,1e-8)
            untrackedEigenvectors[ki,i,j] = problem.eigenvectorsOscillatoryCplx
            # Tip OOP displacement
            tip_u3[ki,i,j] = problem.nodalStatesOverσ[end][nElem].u_n2[3]
        end
        # Frequencies and dampings after mode tracking
        freqs[ki,i,:],damps[ki,i,:],_ = mode_tracking(URange,untrackedFreqs[ki,i,:],untrackedDamps[ki,i,:],untrackedEigenvectors[ki,i,:])
        # Flutter speeds, frequencies and tip displacements of modes at current combination of curvature and root pitch
        dampsCurrentkθ = Array{Vector{Float64}}(undef,nModes)
        freqsCurrentkθ = Array{Vector{Float64}}(undef,nModes)
        for mode in 1:nModes
            dampsCurrentkθ[mode] = [damps[ki,i,j][mode] for j in eachindex(URange)]
            freqsCurrentkθ[mode] = [freqs[ki,i,j][mode] for j in eachindex(URange)]
            # Flutter onset
            iOnset = 1 .+ findall(i -> dampsCurrentkθ[mode][i] < 0 && dampsCurrentkθ[mode][i+1] > 0, 1:length(dampsCurrentkθ[mode])-1)
            if isempty(iOnset) || isempty(filter!(x->x!=1,iOnset))
                push!(flutterOnsetSpeed[ki][i][mode],NaN)
                push!(flutterOnsetFreq[ki][i][mode],NaN)
                push!(flutterOnsetTipDisp[ki][i][mode],NaN)
                push!(flutterOffsetSpeed[ki][i][mode],NaN)
                push!(flutterOffsetFreq[ki][i][mode],NaN)
                push!(flutterOffsetTipDisp[ki][i][mode],NaN)
                continue
            end
            for j in iOnset
                push!(flutterOnsetSpeed[ki][i][mode],interpolate(dampsCurrentkθ[mode][j-1:j],URange[j-1:j],0))
                push!(flutterOnsetFreq[ki][i][mode],interpolate(dampsCurrentkθ[mode][j-1:j],freqsCurrentkθ[mode][j-1:j],0))
                push!(flutterOnsetTipDisp[ki][i][mode],interpolate(dampsCurrentkθ[mode][j-1:j],tip_u3[ki,i,j-1:j],0))
            end
            # Flutter offset
            iOffset = 1 .+ findall(i -> dampsCurrentkθ[mode][i] > 0 && dampsCurrentkθ[mode][i+1] < 0, 1:length(dampsCurrentkθ[mode])-1)
            if isempty(iOffset)
                push!(flutterOffsetSpeed[ki][i][mode],NaN)
                push!(flutterOffsetFreq[ki][i][mode],NaN)
                push!(flutterOffsetTipDisp[ki][i][mode],NaN)
                continue
            end
            for j in iOffset
                push!(flutterOffsetSpeed[ki][i][mode],interpolate(-dampsCurrentkθ[mode][j-1:j],URange[j-1:j],0))
                push!(flutterOffsetFreq[ki][i][mode],interpolate(-dampsCurrentkθ[mode][j-1:j],freqsCurrentkθ[mode][j-1:j],0))
                push!(flutterOffsetTipDisp[ki][i][mode],interpolate(-dampsCurrentkθ[mode][j-1:j],tip_u3[ki,i,j-1:j],0))
            end
        end
    end
end

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(kRange)))
lw = 2
ms = 3
relPath = "/test/outputs/figures/SMWFlutterPrecurvatureRange2"
absPath = string(pwd(),relPath)
mkpath(absPath)
gr()
# Flutter speed vs. tip load
plt1 = plot(xlabel="Root angle [deg]", ylabel="Flutter speed [m/s]", xlims=[0,5], ylims=[0,40])
for (ki,k) in enumerate(kRange)
    plot!([NaN], [NaN], c=colors[ki], ls=:solid, lw=lw, label="k=$k - onset")
    plot!([NaN], [NaN], c=colors[ki], ls=:dash, lw=lw, label="k=$k - offset")
end
for m in 1:nModes
    for (ki,k) in enumerate(kRange)
        fSpeedOnsetMode = vcat([flutterOnsetSpeed[ki][j][m][1] for j in eachindex(θRange)]...)
        fSpeedOffsetMode = vcat([flutterOffsetSpeed[ki][j][m][1] for j in eachindex(θRange)]...)
        plot!(θRange, fSpeedOnsetMode, c=colors[ki], ls=:solid, lw=lw, label=false)
        plot!(θRange, fSpeedOffsetMode, c=colors[ki], ls=:dash, lw=lw, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/SMWFlutterPrecurvatureRange2_flutter.pdf"))

println("Finished SMWFlutterPrecurvatureRange2.jl")