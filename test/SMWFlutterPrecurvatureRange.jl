using AeroBeams, LinearAlgebra, LinearInterpolations, Plots, ColorSchemes, DelimitedFiles

# Wing surface
chord = 1.0
normSparPos = 0.5
aeroSolver = Indicial()
derivationMethod = AD()
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=flatPlate,c=chord,normSparPos=normSparPos)

# Wing beam
θ = π/180*0
L = 16
GJ,EIy,EIz = 1e4,2e4,4e6
ρA,ρIs = 0.75,0.1
ρIy = ρIs*EIy/EIz
ρIz = ρIs-ρIy
nElem = 16
∞ = 1e12
wing = create_Beam(name="wing",length=L,nElements=nElem,C=[isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs)],rotationParametrization="E321",p0=[0;0;θ],aeroSurface=surf)

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
g = 0
h = 20e3
SMWFlutterPrecurvatureRange = create_Model(name="SMWFlutterPrecurvatureRange",beams=[wing],BCs=[clamp],gravityVector=[0;0;-g],altitude=h)

# Set system solver options (limit initial load factor)
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Set precurvature, tip force and airspeed ranges, and initialize outputs
kRange = collect(0:0.018:0.018)
F3Range = collect(0:2:40)
URange = collect(20:0.5:35)
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


# Set number of vibration modes
nModes = 5

# Sweep wing precurvature
for (ki,k) in enumerate(kRange)
    # Update precurvature on beam
    wing.k[2] = k
    update_beam!(wing)
    # Sweep tip force
    for (i,F3) in enumerate(F3Range)
        # Update tip force on model
        tipForce = create_BC(name="tipForce",beam=wing,node=nElem+1,types=["F3A"],values=[F3])
        SMWFlutterPrecurvatureRange.BCs = [clamp,tipForce]
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            # Display progress
            println("Solving for k=$k, F3 = $F3 N, U = $U m/s")
            # Update velocity of basis A (and update model)
            set_motion_basis_A!(model=SMWFlutterPrecurvatureRange,v_A=[0;U;0])
            # Create and solve problem
            problem = create_EigenProblem(model=SMWFlutterPrecurvatureRange,systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
            solve!(problem)
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[ki,i,j] = problem.frequenciesOscillatory
            untrackedDamps[ki,i,j] = round_off!(problem.dampingsOscillatory,1e-12)
            untrackedEigenvectors[ki,i,j] = problem.eigenvectorsOscillatoryCplx
            # Tip OOP displacement
            tip_u3[ki,i,j] = problem.nodalStatesOverσ[end][nElem].u_n2[3]
        end
        # Frequencies and dampings after mode tracking
        freqs[ki,i,:],damps[ki,i,:],_ = mode_tracking(URange,untrackedFreqs[ki,i,:],untrackedDamps[ki,i,:],untrackedEigenvectors[ki,i,:])
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
flutterSpeedVsTipLoadk0 = readdlm(string(pwd(),"/test/referenceData/SMW/flutterSpeedVsTipLoadk0.txt"))
flutterSpeedVsTipLoadk2 = readdlm(string(pwd(),"/test/referenceData/SMW/flutterSpeedVsTipLoadk2.txt"))

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(kRange)))
lw = 2
ms = 3
# Flutter speed vs. tip load
plt1 = plot(xlabel="Tip load [N]", ylabel="Flutter speed [m/s]", xlims=[0,40], ylims=[0,40])
for (ki,k) in enumerate(kRange)
    plot!(F3Range, flutterSpeed[ki,:], c=colors[ki], lw=lw, label="AeroBeams - \$k\$=$k")
end
plot!(flutterSpeedVsTipLoadk0[1,:], flutterSpeedVsTipLoadk0[2,:], c=:black, ls=:dash, lw=lw, label="Patil et al. (2001)")
plot!(flutterSpeedVsTipLoadk2[1,:], flutterSpeedVsTipLoadk2[2,:], c=:black, ls=:dash, lw=lw, label=false)
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/SMWFlutterPrecurvatureRange_1.pdf"))

println("Finished SMWFlutterPrecurvatureRange.jl")