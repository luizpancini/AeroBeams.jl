using AeroBeams, LinearInterpolations, JLD2

# Aerodynamic solver
aeroSolver = Indicial()

# Stiffness factor
λ = 1

# Altitude
h = 20e3

# No gravity
g = 0

# Root pitch angle [rad]
θ = 0*π/180

# Option to include induced drag
hasInducedDrag = true

# Parasite drag
cd0 = 1e-2

# Discretization
nElem = 40

# System solver
maxIter = 50
σ0 = 1
NR = create_NewtonRaphson(maximumIterations=maxIter,initialLoadFactor=σ0,displayStatus=false)

# Set number of vibration modes
nModes = 5

# Set bending curvature and airspeed ranges
k2Range = range(-0.015,0.045,5)
URange = vcat(1e-3,0.5:0.5:50)

# Initialize outputs
untrackedFreqs = Array{Vector{Float64}}(undef,length(k2Range),length(URange))
untrackedDamps = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef, length(k2Range),length(URange))
freqs = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
damps = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
modeDampings = Array{Vector{Float64}}(undef, length(k2Range),nModes)
modeFrequencies = Array{Vector{Float64}}(undef, length(k2Range),nModes)

x1_0 = Array{Vector{Float64}}(undef, length(k2Range))
x3_0 = Array{Vector{Float64}}(undef, length(k2Range))
x1_n = Array{Vector{Float64}}(undef, length(k2Range))
x1_e = Array{Vector{Float64}}(undef, length(k2Range))
u1_of_x1 = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
u3_of_x1 = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
x1_def = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
x3_def = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
twist = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
α = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
cn = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
ct = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
cl = Array{Vector{Float64}}(undef, length(k2Range),length(URange))
cd = Array{Vector{Float64}}(undef, length(k2Range),length(URange))

eigenProblem = Array{EigenProblem}(undef, length(k2Range),length(URange))

# Sweep bending curvature
for (i,k2) in enumerate(k2Range)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        println("Solving for k2 = $k2, U = $U m/s")
        # Model
        model,L = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=U,nElem=nElem,altitude=h,g=g,θ=θ,k2=k2,cd0=cd0,hasInducedDrag=hasInducedDrag)
        # Set initial guess solution as previous known solution
        x0 = j == 1 ? zeros(0) : eigenProblem[i,j-1].x
        # Create and solve eigen problem
        eigenProblem[i,j] = create_EigenProblem(model=model,nModes=nModes,frequencyFilterLimits=[1e-1,Inf],systemSolver=NR,x0=x0)
        solve!(eigenProblem[i,j])
        # Frequencies, dampings and eigenvectors
        untrackedFreqs[i,j] = eigenProblem[i,j].frequenciesOscillatory
        untrackedDamps[i,j] = round_off!(eigenProblem[i,j].dampingsOscillatory,1e-8)
        untrackedEigenvectors[i,j] = eigenProblem[i,j].eigenvectorsOscillatoryCplx
        # Undeformed jig-shape properties
        if j == 1
            # Undeformed nodal positions of right wing
            x1_0[i] = vcat([vcat(model.elements[e].r_n1[1],model.elements[e].r_n2[1]) for e in 1:nElem]...)
            x3_0[i] = vcat([vcat(model.elements[e].r_n1[3],model.elements[e].r_n2[3]) for e in 1:nElem]...)
            # Nodal and elemental arclength positions
            x1_n[i] = vcat([vcat(model.elements[e].x1_n1,model.elements[e].x1_n2) for e in 1:nElem]...)
            x1_e[i] = [model.elements[e].x1 for e in 1:nElem]
        end
        # Displacements over span
        u1_of_x1[i,j] = vcat([vcat(eigenProblem[i,j].nodalStatesOverσ[end][e].u_n1[1],eigenProblem[i,j].nodalStatesOverσ[end][e].u_n2[1]) for e in 1:nElem]...)
        u3_of_x1[i,j] = vcat([vcat(eigenProblem[i,j].nodalStatesOverσ[end][e].u_n1[3],eigenProblem[i,j].nodalStatesOverσ[end][e].u_n2[3]) for e in 1:nElem]...)
        u1_of_x1[i,j] .-= u1_of_x1[i,j][1]
        u3_of_x1[i,j] .-= u3_of_x1[i,j][1]
        # Deformed nodal positions
        x1_def[i,j] = x1_0[i] .+ u1_of_x1[i,j]
        x3_def[i,j] = x3_0[i] .+ u3_of_x1[i,j]
        # Angle of twist
        p1_of_x1 = vcat([vcat(eigenProblem[i,j].nodalStatesOverσ[end][e].p_n1[1],eigenProblem[i,j].nodalStatesOverσ[end][e].p_n2[1]) for e in 1:nElem]...)
        p2_of_x1 = vcat([vcat(eigenProblem[i,j].nodalStatesOverσ[end][e].p_n1[2],eigenProblem[i,j].nodalStatesOverσ[end][e].p_n2[2]) for e in 1:nElem]...)
        p3_of_x1 = vcat([vcat(eigenProblem[i,j].nodalStatesOverσ[end][e].p_n1[3],eigenProblem[i,j].nodalStatesOverσ[end][e].p_n2[3]) for e in 1:nElem]...)
        twist[i,j] = [asind((first(rotation_tensor_WM([p1_of_x1[k],p2_of_x1[k],p3_of_x1[k]]))*AeroBeams.a2)[3]) for k in eachindex(p1_of_x1)]
        # Angle of attack and force coefficients over span
        α[i,j] = 180/π*[eigenProblem[i,j].aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ for e in 1:nElem]
        cn[i,j] = [eigenProblem[i,j].aeroVariablesOverσ[end][e].aeroCoefficients.cn for e in 1:nElem]
        ct[i,j] = [eigenProblem[i,j].aeroVariablesOverσ[end][e].aeroCoefficients.ct for e in 1:nElem]
        cl[i,j] = @. cn[i,j]*cosd(α[i,j]) + ct[i,j]*sind(α[i,j])
        cd[i,j] = @. cn[i,j]*sind(α[i,j]) - ct[i,j]*cosd(α[i,j])
    end
    # Frequencies and dampings after mode tracking
    freqs[i,:],damps[i,:],_ = mode_tracking_hungarian(URange,untrackedFreqs[i,:],untrackedDamps[i,:],untrackedEigenvectors[i,:])
    # Separate frequencies and dampings by mode
    for mode in 1:nModes
        modeFrequencies[i,mode] = [freqs[i,j][mode] for j in eachindex(URange)]
        modeDampings[i,mode] = [damps[i,j][mode] for j in eachindex(URange)]
    end
end

# Compute flutter variables
flutterOnsetMode = [fill(0,0) for _ in 1:length(k2Range)]
flutterOffsetMode = [fill(0,0) for _ in 1:length(k2Range)]
flutterOnsetSpeedOfMode = [fill(0.0,0) for _ in 1:length(k2Range), _ in 1:nModes]
flutterOnsetFreqOfMode = [fill(0.0,0) for _ in 1:length(k2Range), _ in 1:nModes]
flutterOnsetTipOOPOfMode = [fill(0.0,0) for _ in 1:length(k2Range), _ in 1:nModes]
flutterOffsetSpeedOfMode = [fill(0.0,0) for _ in 1:length(k2Range), _ in 1:nModes]
flutterOffsetFreqOfMode = [fill(0.0,0) for _ in 1:length(k2Range), _ in 1:nModes]
flutterOnsetSpeed = [Inf for _ in 1:length(k2Range)]
flutterOnsetFreq = [NaN for _ in 1:length(k2Range)]
flutterOnsetTipOOP = [NaN for _ in 1:length(k2Range)]
flutterOnsetSpeedsAll = [fill(0.0,0) for _ in 1:length(k2Range)]
flutterOnsetFreqsAll = [fill(0.0,0) for _ in 1:length(k2Range)]
flutterOffsetSpeedsAll = [fill(0.0,0) for _ in 1:length(k2Range)]
flutterOffsetFreqsAll = [fill(0.0,0) for _ in 1:length(k2Range)]
indicesFlutterOnset = [fill(0,0) for _ in 1:length(k2Range)]
for (i,k2) in enumerate(k2Range)
    # Flutter onset/offset data of each mode
    for mode in 1:nModes
        iOnset = findall(j -> modeDampings[i,mode][j] < 0 && modeDampings[i,mode][j+1] > 0, 1:length(URange)-1)
        iOffset = findall(j -> modeDampings[i,mode][j] > 0 && modeDampings[i,mode][j+1] < 0, 1:length(URange)-1)
        if modeDampings[i,mode][1]/modeFrequencies[i,mode][1] > 0
            push!(flutterOnsetMode[i],mode)
            push!(flutterOnsetSpeedOfMode[i,mode],URange[1])
            push!(flutterOnsetFreqOfMode[i,mode],modeFrequencies[i,mode][1])
            push!(flutterOnsetTipOOPOfMode[i,mode],x3_def[i,1][end]/16*100)
        end
        if !isempty(iOnset)
            for iO in iOnset
                push!(flutterOnsetMode[i],mode)
                push!(flutterOnsetSpeedOfMode[i,mode],interpolate(modeDampings[i,mode][iO:iO+1],URange[iO:iO+1],0))
                push!(flutterOnsetFreqOfMode[i,mode],interpolate(modeDampings[i,mode][iO:iO+1],modeFrequencies[i,mode][iO:iO+1],0))
                push!(flutterOnsetTipOOPOfMode[i,mode],(x3_def[i,iO][end]+x3_def[i,iO+1][end])/2/16*100)
            end
        end
        if !(isempty(iOffset) || isempty(iOnset))
            for iO in iOffset
                push!(flutterOffsetMode[i],mode)
                push!(flutterOffsetSpeedOfMode[i,mode],interpolate(-modeDampings[i,mode][iO:iO+1],URange[iO:iO+1],0))
                push!(flutterOffsetFreqOfMode[i,mode],interpolate(-modeDampings[i,mode][iO:iO+1],modeFrequencies[i,mode][iO:iO+1],0))
            end
        end
    end
    # All flutter onset/offset speeds and frequencies and tip OOP
    flutterOnsetSpeedsAll[i] = vcat(filter(!isempty,flutterOnsetSpeedOfMode[i,:])...)
    flutterOnsetFreqsAll[i] = vcat(filter(!isempty,flutterOnsetFreqOfMode[i,:])...)
    flutterOffsetSpeedsAll[i] = vcat(filter(!isempty,flutterOffsetSpeedOfMode[i,:])...)
    flutterOffsetFreqsAll[i] = vcat(filter(!isempty,flutterOffsetFreqOfMode[i,:])...)
    flutterOnsetTipOOPAll = vcat(filter(!isempty,flutterOnsetTipOOPOfMode[i,:])...)
    # Ordered indices (of URange) of flutter onset
    flutterOnsetSpeedsAllOrdered = sort(flutterOnsetSpeedsAll[i])
    for speed in flutterOnsetSpeedsAllOrdered
        push!(indicesFlutterOnset[i], -1 + findfirst(x-> x > speed, URange))
    end
    # Lowest flutter onset speed, corresponding frequency and tip OOP
    if !isempty(flutterOnsetSpeedsAll[i])
        iLowest = sortperm(flutterOnsetSpeedsAll[i])[1]
        flutterOnsetSpeed[i] = flutterOnsetSpeedsAll[i][iLowest]
        flutterOnsetFreq[i] = flutterOnsetFreqsAll[i][iLowest]
        flutterOnsetTipOOP[i] = flutterOnsetTipOOPAll[iLowest]
    end
end

using Plots, ColorSchemes

# Set paths
relPath = "/dev/cHALE/Flexible/outputs/figures/cHALEwing_g0_flutter_fixedPitch_k2_range"
absPath = string(pwd(),relPath)
mkpath(absPath)
relPathData = "/dev/cHALE/Flexible/outputs/data/cHALEwing_g0_flutter_fixedPitch_k2_range/"
absPathData = string(pwd(),relPathData)
mkpath(absPathData)

# Save flutter data
flutterOnsetSpeedWingg0 = flutterOnsetSpeed
flutterOnsetFreqWingg0 = flutterOnsetFreq
flutterOnsetTipOOPWingg0 = flutterOnsetTipOOP
@save absPathData*string("wing_g0_lambda",λ,"_flutterSpeed.jld2") flutterOnsetSpeedWingg0
@save absPathData*string("wing_g0_lambda",λ,"_flutterFreq.jld2") flutterOnsetFreqWingg0
@save absPathData*string("wing_g0_lambda",λ,"_flutterOOP.jld2") flutterOnsetTipOOPWingg0

# Plot configurations
colors = get(colorschemes[:rainbow], range(0, 1, length(k2Range)))
ts = 10
fs = 16
lfs = 10
lw = 2
ms = 3
msw = 0
tsz = 9
mshape = [:circle, :star, :utriangle, :pentagon, :diamond]
labels = ["\$k_2 = $(k2) \$" for k2 in k2Range]
L = 16
gr()

# Undeformed jig-shapes
plt_shapes = plot(xlabel="Normalized spanwise direction", ylabel="Normalized vertical direction", xlims=[0,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs, legend=:bottomleft)
for (i,k2) in enumerate(k2Range)
    plot!(x1_0[i]/L, x3_0[i]/L, c=colors[i], lw=lw, label=labels[i])
end
display(plt_shapes)
savefig(string(absPath,"/cHALEwing_undef_k2range.pdf"))

# Root locus
plt_RL = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-10,2], ylims=[0,50], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
scatter!([NaN], [NaN], c=:white, shape=:star8, ms=ms, msw=1, msα=1, msc=:black, markerstrokestyle=:solid, label=string("\$U_{\\infty} = ",URange[1],"\$ m/s"))
for (i,k2) in enumerate(k2Range)
    scatter!([NaN], [NaN], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=labels[i])
    for mode in 1:nModes
        plot!(modeDampings[i,mode], modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=msw, label=false)
        plot!([modeDampings[i,mode][1]], [modeFrequencies[i,mode][1]], c=colors[i], shape=mshape[i], ms=ms, msw=2, msα=1, msc=:black, markerstrokestyle=:solid, label=false)
    end
end
if λ == 1
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
display(plt_RL)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_rootlocus_lambda",λ,".pdf")))

# V-g-f
Vgf_k2_ind = [2,5]
plt_Vf = plot(ylabel="Frequency [rad/s]", xlims=[URange[1],URange[end]], ylims=[0,50], tickfont=font(ts), guidefont=font(12))
for (i,k2) in enumerate(k2Range)
    if !(i in Vgf_k2_ind)
        continue
    end
    for mode in 1:nModes
        plot!(URange, modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
plt_Vg = plot(xlabel="Airspeed [m/s]", ylabel="Damping Ratio", xlims=[URange[1],URange[end]], ylims=[-0.2,0.2], tickfont=font(ts), guidefont=font(12), legendfontsize=10, legend=:topleft)
plot!(plt_Vg,URange,zeros(length(URange)), c=:gray, lw=lw, ls=:dash, label=false)
for (i,k2) in enumerate(k2Range)
    if !(i in Vgf_k2_ind)
        continue
    end
    for mode in 1:nModes
        plot!(URange, modeDampings[i,mode]./modeFrequencies[i,mode], c=colors[i], shape=mshape[i], ms=ms, msw=0, label=false)
    end
end
plt_Vgf = plot(plt_Vf,plt_Vg, layout=(2,1))
display(plt_Vgf)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_Vgf_lambda",λ,".pdf")))

# Normalized deformed span at iminence of flutter onset
plt_disp = plot(xlabel="Normalized spanwise direction", ylabel="Normalized vertical direction", xlims=[0,1], ylims=[-1,1], yticks=-1:0.2:1, tickfont=font(ts), guidefont=font(fs), legendfontsize=8, legend_position=:topright)
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dot, label="Undeformed")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dash, label="At 1st flutter condition")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:solid, label="At 2nd flutter condition")
for (i,k2) in enumerate(k2Range)
    plot!(x1_0[i]/L, x3_0[i]/L, c=colors[i], lw=lw, ls=:dot, label=false)
    plot!(x1_def[i,indicesFlutterOnset[i][1]]/L, x3_def[i,indicesFlutterOnset[i][1]]/L, c=colors[i], lw=lw, ls=:dash, label=false)
    if length(indicesFlutterOnset[i]) > 1
        plot!(x1_def[i,indicesFlutterOnset[i][2]]/L, x3_def[i,indicesFlutterOnset[i][2]]/L, c=colors[i], lw=lw, ls=:solid, label=labels[i])
    end
end
display(plt_disp)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_disp_lambda",λ,".pdf")))

# AoA at iminence of flutter onset
plt_AoA = plot(xlabel="Normalized arclength", ylabel="Angle of attack [deg]", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(x1_e[i]/L, α[i,indicesFlutterOnset[i][1]], c=colors[i], lw=lw, ls=:dash, label=false)
    if length(indicesFlutterOnset[i]) > 1
        plot!(x1_e[i]/L, α[i,indicesFlutterOnset[i][2]], c=colors[i], lw=lw, ls=:solid, label=false)
    end
end
display(plt_AoA)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_AoA_lambda",λ,".pdf")))

# Angle of twist at iminence of flutter onset
plt_twist = plot(xlabel="Normalized arclength", ylabel="Angle of twist [deg]", xlims=[0,1], tickfont=font(ts), guidefont=font(fs))
for (i,k2) in enumerate(k2Range)
    plot!(x1_n[i]/L, twist[i,indicesFlutterOnset[i][1]], c=colors[i], lw=lw, ls=:dash, label=false)
    if length(indicesFlutterOnset[i]) > 1
        plot!(x1_n[i]/L, twist[i,indicesFlutterOnset[i][2]], c=colors[i], lw=lw, ls=:solid, label=false)
    end
end
display(plt_twist)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_twist_lambda",λ,".pdf")))

# Flutter onset speeds vs k2
plt_Uf = plot(xlabel="\$k_2\$ [1/m]", ylabel="Flutter speed [m/s]", xlims=[-0.02,0.05], ylims=[0,35], xticks=k2Range, tickfont=font(ts), guidefont=font(fs))
plot!(k2Range,flutterOnsetSpeed, marker=(:circle, 10), msw=msw, c=:black, lw=lw, label=false)
display(plt_Uf)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_speedOn.pdf")))

# Flutter onset speeds vs tip OOP position
plt_Uf = plot(xlabel="Tip OOP position [% semispan]", ylabel="Flutter speed [m/s]", xlims=[-40,20], ylims=[0,40], xticks=vcat(-40:10:20), tickfont=font(ts), guidefont=font(fs))
plot!(flutterOnsetTipOOP,flutterOnsetSpeed, marker=(:circle, 10), msw=msw, c=:black, lw=lw, label=false)
display(plt_Uf)
savefig(string(absPath,string("/cHALEwing_g0_flutter_fixedPitch_k2_range_speedOn_OOP_lambda",λ,".pdf")))

# # Mode shapes of straight wing at lowest airspeed
# modesPlot = plot_mode_shapes(eigenProblem[2,1],scale=5,view=(60,15),legendPos=(0.25,0.5),modalColorScheme=:rainbow,modeLabels=["1st OOP","2nd OOP","T","IP","3rd OOP"],save=true,savePath=string(relPath,"/cHALEwing_modeShapes_k0U1_lambda",λ,".pdf"))
# display(modesPlot)

# # Mode shapes of straight wing at highest airspeed
# modesPlot = plot_mode_shapes(eigenProblem[2,end],scale=5,view=(60,15),legendPos=(0.25,0.5),modalColorScheme=:rainbow,modeLabels=["1st OOP","2nd OOP","T","IP","3rd OOP"],save=true,savePath=string(relPath,"/cHALEwing_modeShapes_k0U35_lambda",λ,".pdf"))
# display(modesPlot)

# # Mode shapes of wing with k_2=0.045 at lowest airspeed
# modesPlot = plot_mode_shapes(eigenProblem[5,1],scale=5,view=(60,15),legendPos=(0.25,0.5),modalColorScheme=:rainbow,modeLabels=["1st OOP","1st T-IP","2nd OOP","3rd OOP","2nd T-IP"],save=true,savePath=string(relPath,"/cHALEwing_modeShapes_k0045U1_lambda",λ,".pdf"))
# display(modesPlot)

# # Mode shapes of wing with k_2=0.045 at highest airspeed
# modesPlot = plot_mode_shapes(eigenProblem[5,end],scale=5,view=(60,15),legendPos=(0.25,0.5),modalColorScheme=:rainbow,modeLabels=["1st OOP","1st T-IP","2nd OOP","2nd T-IP","3rd OOP"],save=true,savePath=string(relPath,"/cHALEwing_modeShapes_k0045U35_lambda",λ,".pdf"))
# display(modesPlot)

println("Finished cHALEwing_g0_flutter_fixedPitch_k2_range.jl")