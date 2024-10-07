using AeroBeams

# Mode tracking option
modeTracking = true

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Stiffness factor for the aircraft's wing
λ = 1e0

# Altitude
h = 20e3

# Discretization
nElemWing = 20
nElemTailBoom = 10
nElemHorzStabilizer = 10
nElemVertStabilizer = 5

# System solvers
relaxFactor = 0.5
maxIter = 100
σ0 = 1.0
σstep = 0.5
NRtrim = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,displayStatus=false)
NReigen = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=σ0,maximumLoadFactorStep=σstep,displayStatus=true)

# Set number of vibration modes
nModes = 6

# Set bending curvature and airspeed ranges
k2Range = range(-0.015,0.045,5)
URange = collect(20:1:35)

# Initialize outputs
trimAoA = Array{Float64}(undef,length(k2Range),length(URange))

untrackedFreqs = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(k2Range),length(URange))
freqs = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
damps = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
modeDampings = Array{Vector{Float64}}(undef,length(k2Range),nModes)
modeFrequencies = Array{Vector{Float64}}(undef,length(k2Range),nModes)
flutterOnsetSpeedOfMode = Array{Float64}(undef,length(k2Range),nModes)
flutterOffsetSpeedOfMode = Array{Float64}(undef,length(k2Range),nModes)
flutterOnsetSpeed = Array{Float64}(undef,length(k2Range))

x1_0 = Array{Vector{Float64}}(undef,length(k2Range))
x3_0 = Array{Vector{Float64}}(undef,length(k2Range))
x1_e = Array{Vector{Float64}}(undef,length(k2Range))
u1_of_x1 = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
x1_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
x3_def = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
α = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
cn = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
ct = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
cl = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
cd = Array{Vector{Float64}}(undef,length(k2Range),length(URange))

eigenProblem = Array{EigenProblem}(undef,length(k2Range),length(URange))

# Attachment springs
μu = [1e-1; 1e-1; 1e-1; 1e-1; 1e-1]
μp = [10; 1; 10; 1; 0.1]

# ELement ranges
elemRangeRightWing = 1 + div(nElemWing,2) : nElemWing

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Set attachment springs
    spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=μu[i]*[1; 1; 1],kp=μp[i]*[1; 1; 1])
    spring2 = create_Spring(elementsIDs=[nElemTailBoom],nodesSides=[2],ku=μu[i]*[1; 1; 1],kp=μp[i]*[1; 1; 1])
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model for trim problem
        cHALEtrim,_,_,tailBoom,_ = create_conventional_HALE(aeroSolver=aeroSolver,airspeed=U,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stiffnessFactor=λ,∞=1e12,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag,altitude=h)
        # Add springs
        add_springs_to_beam!(beam=tailBoom,springs=[spring1,spring2])
        # Update model
        cHALEtrim.skipValidationMotionBasisA = true
        update_model!(cHALEtrim)
        # Set initial guess solution as previous known solution
        x0Trim = j == 1 ? zeros(0) : trimProblem.x
        # Create and trim problem
        global trimProblem = create_TrimProblem(model=cHALEtrim,systemSolver=NRtrim,x0=x0Trim)
        solve!(trimProblem)
        # Extract trim variables
        trimAoA[i,j] = trimProblem.aeroVariablesOverσ[end][div(nElemWing,2)].flowAnglesAndRates.αₑ
        println("Trim AoA = $(trimAoA[i,j]*180/π)")
        # Undeformed jig-shape properties
        if j == 1
            # Undeformed nodal positions of right wing
            x1_0[i] = vcat([vcat(cHALEtrim.elements[e].r_n1[1],cHALEtrim.elements[e].r_n2[1]) for e in elemRangeRightWing]...)
            x3_0[i] = vcat([vcat(cHALEtrim.elements[e].r_n1[3],cHALEtrim.elements[e].r_n2[3]) for e in elemRangeRightWing]...)
            # Undeformed elemental positions
            x1_e[i] = [cHALEtrim.elements[e].x1 for e in elemRangeRightWing]
        end
        # Displacements over span
        u1_of_x1[i,j] = vcat([vcat(trimProblem.nodalStatesOverσ[end][e].u_n1[1],trimProblem.nodalStatesOverσ[end][e].u_n2[1]) for e in elemRangeRightWing]...)
        u3_of_x1[i,j] = vcat([vcat(trimProblem.nodalStatesOverσ[end][e].u_n1[3],trimProblem.nodalStatesOverσ[end][e].u_n2[3]) for e in elemRangeRightWing]...)
        u1_of_x1[i,j] .-= u1_of_x1[i,j][1]
        u3_of_x1[i,j] .-= u3_of_x1[i,j][1]
        # Deformed nodal positions
        x1_def[i,j] = x1_0[i] .+ u1_of_x1[i,j]
        x3_def[i,j] = x3_0[i] .+ u3_of_x1[i,j]
        # Angle of attack and force coefficients over span
        α[i,j] = 180/π*[trimProblem.aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in elemRangeRightWing]
        cn[i,j] = [trimProblem.aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in elemRangeRightWing]
        ct[i,j] = [trimProblem.aeroVariablesOverσ[end][e].aeroCoefficients.ct for e in elemRangeRightWing]
        cl[i,j] = @. cn[i,j]*cosd(α[i,j]) + ct[i,j]*sind(α[i,j])
        cd[i,j] = @. cn[i,j]*sind(α[i,j]) - ct[i,j]*cosd(α[i,j])
        # Model for eigen problem
        wingModel,_ = create_SMW(aeroSolver=aeroSolver,airspeed=U,nElem=nElemWing,stiffnessFactor=λ,∞=1e12,altitude=h,cd0=wingCd0,k2=k2,hasInducedDrag=hasInducedDrag,θ=α[i,j][1]*π/180)
        # Set initial guess solution as previous known solution
        x0Eig = j == 1 ? zeros(0) : eigenProblem[i,j-1].x
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=wingModel,nModes=nModes,frequencyFilterLimits=[1e-2,Inf64],systemSolver=NReigen,x0=x0Eig)
        solve!(eigenProblem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = eigenProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(eigenProblem[i,j].dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = eigenProblem[i,j].eigenvectorsOscillatoryCplx
    end
    # Frequencies and dampings after mode tracking
    if modeTracking
        freqs[i,:],damps[i,:],_ = mode_tracking(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    else
        freqs[i,:],damps[i,:] = untrackedFreqs[i,:],untrackedDamps[i,:]
    end
    # Separate frequencies and dampings by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
    end
    # Flutter speed of each mode
    for mode in 1:nModes
        iOnset = findfirst(j -> modeDampings[i,mode][j] < 0 && modeDampings[i,mode][j+1] > 0, 1:length(URange)-1)
        iOffset = findfirst(j -> modeDampings[i,mode][j] > 0 && modeDampings[i,mode][j+1] < 0, 1:length(URange)-1)
        if isnothing(iOnset)
            flutterOnsetSpeedOfMode[i,mode] = Inf64
        else
            flutterOnsetSpeedOfMode[i,mode] = interpolate(modeDampings[i,mode][iOnset:iOnset+1],URange[iOnset:iOnset+1],0)
        end
        if isnothing(iOffset) || isnothing(iOnset)
            flutterOffsetSpeedOfMode[i,mode] = Inf64
        else
            flutterOffsetSpeedOfMode[i,mode] = interpolate(-modeDampings[i,mode][iOffset:iOffset+1],URange[iOffset:iOffset+1],0)
        end
    end
    flutterOnsetSpeed[i] = minimum(filter(!isinf,flutterOnsetSpeedOfMode[i,:]),init=Inf64)
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/outputs/figures/cHALEwing_flutter_matchedPitch_k2_range.jl"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Mode shapes at lowest airspeed
# for (n,k2) in enumerate(k2Range)
#     plt = plot_mode_shapes(eigenProblem[n,1],nModes=nModes,scale=5,view=(30,30),modalColorScheme=:rainbow,legendPos=:outertop,save=true,savePath=string(relPath,string("/cHALEwing_flutter_matchedPitch_k2_range_modeShapesk2_",n,".pdf")))
#     display(plt)
# end

# Plot configurations
colors = get(colorschemes[:rainbow], range(0, 1, length(k2Range)))
lw = 2
ms = 3
msw = 0
mshape = [:circle, :star, :utriangle, :pentagon, :diamond]
labels = ["\$k_2 = $(k2) \$" for k2 in k2Range]
L = 16
gr()

# Root locus
plt0 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,1], ylims=[0,50])
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        scatter!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        scatter!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
tsz = 10
OOPtext = "OOP"
TIP1text = "1st T-IP"
TIP2text = "2nd T-IP"
Ttext = "T"
IPtext = "IP"
OOPtextPos = [-7.5, 20]
TIP1textPos = [0.5, 10]
TIP2textPos = [0.25, 33]
annotate!(OOPtextPos[1], OOPtextPos[2], text(OOPtext, tsz))
annotate!(TIP1textPos[1], TIP1textPos[2], text(TIP1text, tsz))
annotate!(TIP2textPos[1], TIP2textPos[2], text(TIP2text, tsz))
quiver!([OOPtextPos[1],OOPtextPos[1]+0.5,OOPtextPos[1]+0.5], [OOPtextPos[2]-2,OOPtextPos[2]-0.5,OOPtextPos[2]+1], quiver=([0,1.75,3.5], [-10,-3.5,15]), arrow=:closed, linecolor=:black)
display(plt0)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_rootlocus.pdf"))

# Flutter onset and offset speeds vs k2
plt1 = plot(xlabel="\$k_2\$", ylabel="Flutter speed [m/s]", xlims=[-0.02,0.05], ylims=[0,40], xticks=k2Range)
scatter!([NaN], [NaN], shape=:utriangle, ms=10, msw=0, c=:red, lw=lw, label="Flutter onset")
scatter!([NaN], [NaN], shape=:utriangle, ms=10, msw=0, c=:green, lw=lw, label="Flutter offset")
for m = 1:nModes
    plot!(k2Range,flutterOnsetSpeedOfMode[:,m], marker=(:utriangle, 10), msw=0, c=:red, lw=lw, label=false)
    plot!(k2Range,flutterOffsetSpeedOfMode[:,m], marker=(:dtriangle, 10), msw=0, c=:green, lw=lw, label=false)
end
display(plt1)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_speed.pdf"))

# Flutter onset speeds vs k2
plt2 = plot(xlabel="\$k_2\$", ylabel="Flutter speed [m/s]", xlims=[-0.02,0.05], ylims=[0,35], xticks=k2Range)
plot!(k2Range,flutterOnsetSpeed, marker=(:circle, 10), msw=0, c=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_speedOn.pdf"))

# V-g-f
plt31 = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,50])
for (i,k2) in enumerate(k2Range)
    for mode in 1:nModes
        scatter!(URange, modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
plt32 = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.25,0.05], legend=:bottomright)
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        scatter!(URange, modeDampings[i,mode]./modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
    end
end
plt3 = plot(plt31,plt32, layout=(2,1))
display(plt3)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_Vgf.pdf"))

# Normalized deformed span at lowest and highest airspeeds
plt4 = plot(xlabel="Normalized spanwise direction", ylabel="Normalized OOP direction", xlims=[0,1], legend=:topleft)
plot!([NaN], [NaN], ls=:solid, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
plot!([NaN], [NaN], ls=:dash, c=:black, lw=lw, label=string("\$U_{\\infty} = ",URange[end],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    plot!(x1_def[i,1]/wingsemispan, x3_def[i,1]/wingsemispan, ls=:solid, c=colors[i], lw=lw, label=labels[i])
    plot!(x1_def[i,end]/wingsemispan, x3_def[i,end]/wingsemispan, ls=:dash, c=colors[i], lw=lw, label=false)
end
display(plt4)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_disp.pdf"))

# Angle of attack over span at lowest airspeed
plt5 = plot(xlabel="Normalized span", ylabel="Wing \$\\alpha\$ [deg]", xlims=[0,1])
for (i,k2) in enumerate(k2Range)
    plot!(x1_e[i]/L, α[i,1], c=colors[i], lw=lw, label=labels[i])
end
display(plt5)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_AoA.pdf"))

# Cl over span at lowest airspeed
plt6 = plot(xlabel="Normalized span", ylabel="Wing \$c_l\$ []", xlims=[0,1])
for (i,k2) in enumerate(k2Range)
    plot!(x1_e[i]/L, cl[i,1], c=colors[i], lw=lw, label=labels[i])
end
display(plt6)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_cl.pdf"))

# Cd over span at lowest airspeed
plt7 = plot(xlabel="Normalized span", ylabel="Wing \$c_d\$ []", xlims=[0,1])
for (i,k2) in enumerate(k2Range)
    plot!(x1_e[i]/L, cd[i,1], c=colors[i], lw=lw, label=labels[i])
end
display(plt7)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_cd.pdf"))

# Cn over span at lowest airspeed
plt8 = plot(xlabel="Normalized span", ylabel="Wing \$c_n\$ []", xlims=[0,1])
for (i,k2) in enumerate(k2Range)
    plot!(x1_e[i]/L, cn[i,1], c=colors[i], lw=lw, label=labels[i])
end
display(plt8)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_cn.pdf"))

# Ct over span at lowest airspeed
plt9 = plot(xlabel="Normalized span", ylabel="Wing \$c_t\$ []", xlims=[0,1])
for (i,k2) in enumerate(k2Range)
    plot!(x1_e[i]/L, ct[i,1], c=colors[i], lw=lw, label=labels[i])
end
display(plt9)
savefig(string(absPath,"/cHALEwing_flutter_matchedPitch_k2_range_ct.pdf"))

println("Finished cHALEwing_flutter_matchedPitch_k2_range.jl")