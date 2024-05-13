using AeroBeams, LinearAlgebra, LinearInterpolations, Plots, ColorSchemes, DelimitedFiles

# Wing surface
chord = 1.0
normSparPos = 0.5
aeroSolver = Indicial()
derivationMethod = AD()
surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=flatPlate,c=chord,normSparPos=normSparPos)

# Set torsion / chordwise bending coupling factor range
ΨRange = collect(-0.2:0.2:0.2)

# Wing beams
wings = Vector{Beam}(undef,length(ΨRange))
θ = π/180*0
L = 16
GJ,EIy,EIz = 1e4,2e4,4e6
ρA,ρIs = 0.75,0.1
ρIy = ρIs*EIy/EIz
ρIz = ρIs-ρIy
nElem = 16
∞ = 1e12
Ciso = isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)
for (i,Ψ) in enumerate(ΨRange)
    C = Ciso
    C[4,6] = C[6,4] = -Ψ*sqrt(GJ*EIz)
    wings[i] = create_Beam(name="beam",length=L,nElements=nElem,C=[C],I=[inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs)],rotationParametrization="E321",p0=[0;0;θ],aeroSurface=surf)
end

# BCs
clamps = Vector{BC}(undef,length(ΨRange))
for (i,Ψ) in enumerate(ΨRange)
    clamps[i] = create_BC(name="clamp",beam=wings[i],node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
end

# Models
SMWFlutterStructuralCouplingRange = Vector{Model}(undef,length(ΨRange))
g = 0
h = 20e3
for (i,Ψ) in enumerate(ΨRange)
    SMWFlutterStructuralCouplingRange[i] = create_Model(name="SMWFlutterStructuralCouplingRange",beams=[wings[i]],BCs=[clamps[i]],gravityVector=[0;0;-g],altitude=h)
end

# Set system solver options (limit initial load factor)
σ0 = 0.5
σstep = 0.5
NR = create_NewtonRaphson(initialLoadFactor=σ0,maximumLoadFactorStep=σstep)

# Set number of vibration modes
nModes = 5

# Set tip load and airspeed ranges, and initialize outputs
F3Range = collect(0:1:35)
URange = collect(20:0.25:33)
freqs = Array{Vector{Float64}}(undef,length(ΨRange),length(F3Range),length(URange))
damps = Array{Vector{Float64}}(undef,length(ΨRange),length(F3Range),length(URange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
flutterSpeed = Array{Float64}(undef,length(ΨRange),length(F3Range))

# Sweep structural coupling 
for (i,Ψ) in enumerate(ΨRange)
    # Sweep tip force
    for (j,F3) in enumerate(F3Range)
        # Update tip force on model
        tipForce = create_BC(name="tipForce",beam=wings[i],node=nElem+1,types=["F3A"],values=[F3])
        SMWFlutterStructuralCouplingRange[i].BCs = [clamps[i],tipForce]
        # Sweep airspeed
        for (k,U) in enumerate(URange)
            # Display progress
            println("Solving for Ψ=$Ψ, F=$F3 N, U = $U m/s")
            # Update velocity of basis A (and update model)
            set_motion_basis_A!(model=SMWFlutterStructuralCouplingRange[i],v_A=[0;U;0])
            # Create and solve problem
            problem = create_EigenProblem(model=SMWFlutterStructuralCouplingRange[i],systemSolver=NR,nModes=nModes,frequencyFilterLimits=[1e-3,Inf64])
            solve!(problem)
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[k] = problem.frequenciesOscillatory
            untrackedDamps[k] = round_off!(problem.dampingsOscillatory,1e-12)
            untrackedEigenvectors[k] = problem.eigenvectorsOscillatoryCplx
        end
        # Frequencies and dampings after mode tracking
        freqs[i,j,:],damps[i,j,:],_ = mode_tracking(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)
        # Flutter speeds of modes at current structural coupling and tip force
        dampsCurrentΨF3 = Array{Vector{Float64}}(undef,nModes)
        flutterSpeedOfMode = Array{Float64}(undef,nModes)
        for mode in 1:nModes
            dampsCurrentΨF3[mode] = [damps[i,j,k][mode] for k in eachindex(URange)]
            indexInstability = findfirst(x->x>0,dampsCurrentΨF3[mode])
            flutterSpeedOfMode[mode] = isnothing(indexInstability) ? NaN : interpolate(dampsCurrentΨF3[mode][indexInstability-1:indexInstability],URange[indexInstability-1:indexInstability],0)
        end
        # Set flutter speed as the greatest (for compatibility with reference)
        flutterSpeed[i,j] = maximum(filter(!isnan,flutterSpeedOfMode))
    end
end

# Load reference data
flutterSpeedVsTipLoadΨ0 = readdlm(string(pwd(),"/test/referenceData/SMW/flutterSpeedVsTipLoadPsi0_0.txt"))
flutterSpeedVsTipLoadΨp02 = readdlm(string(pwd(),"/test/referenceData/SMW/flutterSpeedVsTipLoadPsi0_2.txt"))
flutterSpeedVsTipLoadΨm02 = readdlm(string(pwd(),"/test/referenceData/SMW/flutterSpeedVsTipLoadPsi-0_2.txt"))

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(ΨRange)))
lw = 2
ms = 3
# Flutter speed vs tip load for several structural couplings
plt1 = plot(xlabel="Tip Load [N]", ylabel="Flutter speed [m/s]", xlims=[0,35], ylims=[0,35])
plot!([NaN], [NaN], c=:black, lw=lw, legend=:bottomleft, label="AeroBeams")
scatter!([NaN], [NaN], mc=:black, ms=ms, label="Patil et al. (2001)")
for (i,Ψ) in enumerate(ΨRange)
    plot!(F3Range, flutterSpeed[i,:], c=colors[i], lw = lw, label="\\Psi=$Ψ")
    if Ψ == -0.2
        scatter!(flutterSpeedVsTipLoadΨm02[1,:], flutterSpeedVsTipLoadΨm02[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif Ψ == 0.0
        scatter!(flutterSpeedVsTipLoadΨ0[1,:], flutterSpeedVsTipLoadΨ0[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif Ψ == +0.2
        scatter!(flutterSpeedVsTipLoadΨp02[1,:], flutterSpeedVsTipLoadΨp02[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    end
end
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/SMWFlutterStructuralCouplingRange_1.pdf"))

println("Finished SMWFlutterStructuralCouplingRange.jl")