using AeroBeams

# Mode tracking option
modeTracking = true

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor for the aircraft's wing
λ = 1e0

# Altitude
h = 20e3

# No gravity (in-vacuo)
g = 0

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
θRange = [0]
k2Range = range(-0.015,0.045,5)
URange = collect(1:0.5:35)

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(θRange),length(k2Range),length(URange))
freqs = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
damps = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),length(URange))
modeDampings = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),nModes)
modeFrequencies = Array{Vector{Float64}}(undef,length(θRange),length(k2Range),nModes)
flutterOnsetSpeedOfMode = Array{Float64}(undef,length(θRange),length(k2Range),nModes)
flutterOffsetSpeedOfMode = Array{Float64}(undef,length(θRange),length(k2Range),nModes)
flutterOnsetSpeed = Array{Float64}(undef,length(θRange),length(k2Range))

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
            wingModel,_ = create_SMW(aeroSolver=aeroSolver,airspeed=U,nElem=nElem,stiffnessFactor=λ,∞=1e12,altitude=h,g=g,θ=θ*π/180,k2=k2,cd0=cd0,hasInducedDrag=hasInducedDrag)
            # Set initial guess solution as previous known solution
            x0Eig = j == 1 ? zeros(0) : eigenProblem[n,i,j-1].x
            # Create and solve eigen problem
            eigenProblem[n,i,j] = create_EigenProblem(model=wingModel,nModes=nModes,frequencyFilterLimits=[1e-1,Inf64],systemSolver=NReigen,x0=x0Eig)
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
        # Flutter speed of each mode
        for mode in 1:nModes
            iOnset = findfirst(j -> modeDampings[n,i,mode][j] < 0 && modeDampings[n,i,mode][j+1] > 0, 1:length(URange)-1)
            iOffset = findfirst(j -> modeDampings[n,i,mode][j] > 0 && modeDampings[n,i,mode][j+1] < 0, 1:length(URange)-1)
            if isnothing(iOnset)
                flutterOnsetSpeedOfMode[n,i,mode] = Inf64
            else
                flutterOnsetSpeedOfMode[n,i,mode] = interpolate(modeDampings[n,i,mode][iOnset:iOnset+1],URange[iOnset:iOnset+1],0)
            end
            if isnothing(iOffset) || isnothing(iOnset)
                flutterOffsetSpeedOfMode[n,i,mode] = Inf64
            else
                flutterOffsetSpeedOfMode[n,i,mode] = interpolate(-modeDampings[n,i,mode][iOffset:iOffset+1],URange[iOffset:iOffset+1],0)
            end
        end
        flutterOnsetSpeed[n,i] = minimum(filter(!isinf,flutterOnsetSpeedOfMode[n,i,:]),init=Inf64)
    end
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/outputs/figures/cHALEwing_invacuo_flutter_fixedPitch_k2_range"
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

# Undeformed jig-shapes
plt0 = plot(xlabel="Normalized spanwise direction", ylabel="Normalized OOP direction", xlims=[0,1], legend=:bottomleft)
for (i,k2) in enumerate(k2Range)
    plot!(x1_0[i]/16, x3_0[i]/16, c=colors[i], lw=lw, label=labels[i])
end
display(plt0)
savefig(string(absPath,"/cHALEwing_undef_k2range.pdf"))

# Loop root pitch angles
for (n,θ) in enumerate(θRange)

    # Root locus
    plt1 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,2], ylims=[0,50])
    scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
    for (i,k2) in enumerate(k2Range)
        scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
        for mode in 1:nModes
            scatter!(modeDampings[n,i,mode], modeFrequencies[n,i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
            scatter!([modeDampings[n,i,mode][1]], [modeFrequencies[n,i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
        end
    end
    if n == 1
        tsz = 8
        OOPtext = "OOP"
        TIP1text = "1st T-IP"
        TIP2text = "2nd T-IP"
        Ttext = "T"
        IPtext = "IP"
        OOPtextPos = [-8, 20]
        TIP1textPos = [1, 5]
        TIP2textPos = [0.75, 45]
        TtextPos = [0.25, 30]
        IPtextPos = [0.25, 32]
        annotate!(OOPtextPos[1], OOPtextPos[2], text(OOPtext, tsz))
        annotate!(TIP1textPos[1], TIP1textPos[2], text(TIP1text, tsz))
        annotate!(TIP2textPos[1], TIP2textPos[2], text(TIP2text, tsz))
        annotate!(TtextPos[1], TtextPos[2], text(Ttext, tsz))
        annotate!(IPtextPos[1], IPtextPos[2], text(IPtext, tsz))
        quiver!([OOPtextPos[1],OOPtextPos[1]+0.5,OOPtextPos[1]+0.5], [OOPtextPos[2]-2,OOPtextPos[2]-0.5,OOPtextPos[2]+1], quiver=([0,1.75,3.75], [-10,-3.5,17]), arrow=:closed, linecolor=:black)
    end
    display(plt1)
    savefig(string(absPath,string("/cHALEwing_invacuo_flutter_fixedPitch",round(Int,θ),"_k2_range_rootlocus.pdf")))

    # V-g-f
    plt21 = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,50])
    for (i,k2) in enumerate(k2Range)
        for mode in 1:nModes
            scatter!(URange, modeFrequencies[n,i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        end
    end
    plt22 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.25,0.15], legend=:bottomright)
    for (i,k2) in enumerate(k2Range)
        scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
        for mode in 1:nModes
            scatter!(URange, modeDampings[n,i,mode]./modeFrequencies[n,i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        end
    end
    plt2 = plot(plt21,plt22, layout=(2,1))
    display(plt2)
    savefig(string(absPath,string("/cHALEwing_invacuo_flutter_fixedPitch",round(Int,θ),"_k2_range_Vgf.pdf")))

    # Normalized deformed span at highest airspeed
    plt3 = plot(xlabel="Normalized spanwise direction", ylabel="Normalized OOP direction", xlims=[0,1])
    for (i,k2) in enumerate(k2Range)
        plot!(x1_def[n,i,end]/L, x3_def[n,i,end]/L, c=colors[i], lw=lw, label=labels[i])
    end
    display(plt3)
    savefig(string(absPath,string("/cHALEwing_invacuo_flutter_fixedPitch",round(Int,θ),"_k2_range_disp.pdf")))

    # AoA at highest airspeed
    plt4 = plot(xlabel="\$x_1/L\$", ylabel="Angle of attack [deg]", xlims=[0,1])
    for (i,k2) in enumerate(k2Range)
        plot!(x1_e[i]/L, α[n,i,end], c=colors[i], lw=lw, label=labels[i])
    end
    display(plt4)
    savefig(string(absPath,string("/cHALEwing_invacuo_flutter_fixedPitch",round(Int,θ),"_k2_range_AoA.pdf")))

    # Flutter onset and offset speeds vs k2
    plt5 = plot(xlabel="\$k_2\$", ylabel="Flutter speed [m/s]", xlims=[-0.02,0.05], ylims=[0,40], xticks=k2Range)
    scatter!([NaN], [NaN], shape=:utriangle, ms=10, msw=0, c=:red, lw=lw, label="Flutter onset")
    scatter!([NaN], [NaN], shape=:utriangle, ms=10, msw=0, c=:green, lw=lw, label="Flutter offset")
    for m = 1:nModes
        plot!(k2Range,flutterOnsetSpeedOfMode[n,:,m], marker=(:utriangle, 10), msw=0, c=:red, lw=lw, label=false)
        plot!(k2Range,flutterOffsetSpeedOfMode[n,:,m], marker=(:dtriangle, 10), msw=0, c=:green, lw=lw, label=false)
    end
    display(plt5)
    savefig(string(absPath,string("/cHALEwing_invacuo_flutter_fixedPitch",round(Int,θ),"_k2_range_speed.pdf")))

    # Flutter onset speeds vs k2
    plt6 = plot(xlabel="\$k_2\$", ylabel="Flutter speed [m/s]", xlims=[-0.02,0.05], ylims=[0,35], xticks=k2Range)
    plot!(k2Range,flutterOnsetSpeed[n,:], marker=(:circle, 10), msw=0, c=:black, lw=lw, label=false)
    display(plt6)
    savefig(string(absPath,string("/cHALEwing_invacuo_flutter_fixedPitch",round(Int,θ),"_k2_range_speedOn.pdf")))
end


println("Finished cHALEwing_invacuo_flutter_fixedPitch_k2_range.jl")