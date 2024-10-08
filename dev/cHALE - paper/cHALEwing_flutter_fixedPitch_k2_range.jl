using AeroBeams

# Mode tracking option
modeTracking = true

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor for the aircraft's wing
λ = 1e0

# Altitude [ gravity is matched to altitude automatically in create_SMW() ]
h = 20e3

# Option to include induced drag
hasInducedDrag = true

# Parasite drag
cd0 = 1e-2

# Discretization
nElem = 20

# System solvers
relaxFactor = 0.5
maxIter = 100
σ0 = 0.5
σstep = 0.5
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,displayStatus=false)
NReigen = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,displayStatus=false)

# Set number of vibration modes
nModes = 6

# Set root pitch angle [deg], bending curvature and airspeed ranges
θRange = [0,5,10,15]
k2Range = range(-0.015,0.045,5)
URange = collect(20:0.5:35)

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(θRange),length(k2Range),length(URange))
freqs = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
damps = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
modeDampings = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),nModes)
modeFrequencies = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),nModes)

x1_0 = Array{Vector{Float64}}(undef,length(θRange),length(k2Range))
x3_0 = Array{Vector{Float64}}(undef,length(θRange),length(k2Range))
x1_e = Array{Vector{Float64}}(undef,length(θRange),length(k2Range))
u1_of_x1 = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
x1_def = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
x3_def = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
α = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
cn = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
ct = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
cl = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
cd = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))

eigenProblem = Array{EigenProblem}(undef,length(θRange),length(k2Range),length(URange))

# Sweep root pitch angle
for (n,θ) in enumerate(θRange)
    # Sweep bending curvature
    for (i,k2) in enumerate(k2Range)
        # Sweep airspeed
        for (j,U) in enumerate(URange)
            println("Solving for θ = $θ deg, k2 = $k2, U = $U m/s")
            # Model
            wingModel,_ = create_SMW(aeroSolver=aeroSolver,airspeed=U,nElem=nElem,stiffnessFactor=λ,∞=1e12,altitude=h,θ=θ*π/180,k2=k2,cd0=cd0,hasInducedDrag=hasInducedDrag)
            # Set initial guess solution as previous known solution
            x0Eig = j == 1 ? zeros(0) : eigenProblem[n,i,j-1].x
            # Create and solve eigen problem
            eigenProblem[n,i,j] = create_EigenProblem(model=wingModel,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],systemSolver=NReigen,x0=x0Eig)
            solve!(eigenProblem[n,i,j])
            # Frequencies, dampings and eigenvectors
            untrackedFreqs[n,i,j] = eigenProblem[n,i,j].frequenciesOscillatory
            untrackedDamps[n,i,j] = round_off!(eigenProblem[n,i,j].dampingsOscillatory,1e-8)
            untrackedEigenvectors[n,i,j] = eigenProblem[n,i,j].eigenvectorsOscillatoryCplx
            # Undeformed jig-shape properties
            if j == 1
                # Undeformed nodal positions of right wing
                x1_0[i] = vcat([vcat(wingModel.elements[e].r_n1[1],wingModel.elements[e].r_n2[1]) for e in 1:nElem]...)
                x3_0[i] = vcat([vcat(wingModel.elements[e].r_n1[3],wingModel.elements[e].r_n2[3]) for e in 1:nElem]...)
                # Undeformed elemental positions
                x1_e[i] = [wingModel.elements[e].x1 for e in 1:nElem]
            end
            # Displacements over span
            u1_of_x1[n,i,j] = vcat([vcat(eigenProblem[n,i,j].nodalStatesOverσ[end][e].u_n1[1],eigenProblem[n,i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
            u3_of_x1[n,i,j] = vcat([vcat(eigenProblem[n,i,j].nodalStatesOverσ[end][e].u_n1[3],eigenProblem[n,i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
            u1_of_x1[n,i,j] .-= u1_of_x1[n,i,j][1]
            u3_of_x1[n,i,j] .-= u3_of_x1[n,i,j][1]
            # Deformed nodal positions
            x1_def[n,i,j] = x1_0[i] .+ u1_of_x1[n,i,j]
            x3_def[n,i,j] = x3_0[i] .+ u3_of_x1[n,i,j]
            # Angle of attack and force coefficients over span
            α[n,i,j] = 180/π*[eigenProblem[n,i,j].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:nElem]
            cn[n,i,j] = [eigenProblem[n,i,j].aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:nElem]
            ct[n,i,j] = [eigenProblem[n,i,j].aeroVariablesOverσ[end][e].aeroCoefficients.ct for e in 1:nElem]
            cl[n,i,j] = @. cn[n,i,j]*cosd(α[n,i,j]) + ct[n,i,j]*sind(α[n,i,j])
            cd[n,i,j] = @. cn[n,i,j]*sind(α[n,i,j]) - ct[n,i,j]*cosd(α[n,i,j])
        end
        # Frequencies and dampings after mode tracking
        if modeTracking
            freqs[n,i,:],damps[n,i,:],_ = mode_tracking(URange,untrackedFreqs[n,i,:],untrackedDamps[n,i,:],untrackedEigenvectors[n,i,:])
        else
            freqs[n,i,:],damps[n,i,:] = untrackedFreqs[n,i,:],untrackedDamps[n,i,:]
        end
        # Separate frequencies and dampings by mode
        for mode in 1:nModes
            modeFrequencies[n,i,mode] = [freqs[n,i,j][mode] for j in eachindex(URange)]
            modeDampings[n,i,mode] = [damps[n,i,j][mode] for j in eachindex(URange)]
        end
    end
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/outputs/figures/cHALEwing_flutter_fixedPitch_k2_range"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
colors = get(colorschemes[:rainbow], range(0, 1, length(k2Range)))
lw = 2
ms = 3
msw = 0
mshape = [:circle, :star, :utriangle, :pentagon, :diamond]
labels = ["\$k_2 = $(k2) \$" for k2 in k2Range]
L = 16
gr()

# Loop root pitch angles
for (n,θ) in enumerate(θRange)

    # Root locus
    plt0 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,2], ylims=[0,50])
    scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
    for (i,k2) in enumerate(k2Range)
        scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
        for mode in 1:nModes
            scatter!(modeDampings[n,i,mode], modeFrequencies[n,i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
            scatter!([modeDampings[n,i,mode][1]], [modeFrequencies[n,i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
        end
    end
    display(plt0)
    savefig(string(absPath,string("/cHALEwing_flutter_fixedPitch",round(Int,θ),"_k2_range_rootlocus.pdf")))

    # V-g-f
    plt31 = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,50])
    for (i,k2) in enumerate(k2Range)
        for mode in 1:nModes
            scatter!(URange, modeFrequencies[n,i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        end
    end
    plt32 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.25,0.15], legend=:bottomright)
    for (i,k2) in enumerate(k2Range)
        scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
        for mode in 1:nModes
            scatter!(URange, modeDampings[n,i,mode]./modeFrequencies[n,i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        end
    end
    plt3 = plot(plt31,plt32, layout=(2,1))
    display(plt3)
    savefig(string(absPath,string("/cHALEwing_flutter_fixedPitch",round(Int,θ),"_k2_range_Vgf.pdf")))

    # Normalized deformed span at lowest airspeed
    plt4 = plot(xlabel="\$x_1/L\$", ylabel="\$x_3/L\$", xlims=[0,1], title="Airspeed = $(URange[1]) m/s")
    for (i,k2) in enumerate(k2Range)
        plot!(x1_def[n,i,1]/L, x3_def[n,i,1]/L, c=colors[i], lw=lw, label=labels[i])
    end
    display(plt4)
    savefig(string(absPath,string("/cHALEwing_flutter_fixedPitch",round(Int,θ),"_k2_range_disp.pdf")))
end

println("Finished cHALEwing_flutter_fixedPitch_k2_range.jl")