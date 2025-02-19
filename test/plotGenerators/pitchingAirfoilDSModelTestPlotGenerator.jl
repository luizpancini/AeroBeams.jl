using Plots, ColorSchemes

# Run the script
include("../examples/pitchingAirfoilDSModelTest.jl")

# Set paths
relPath = "/test/outputs/figures/pitchingAirfoilDSModelTest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
nTotalAeroStates = problem.model.elements[1].aero.nTotalAeroStates
colors = get(colorschemes[:rainbow], LinRange(0, 1, nTotalAeroStates))
lw = 2
ms = 3
gr()
if source == "NASA"
    refLabel = "Exp - McAlister et al. (1982)"
elseif source == "GU"
    refLabel = "Exp. - Green & Giuni (2017)"
end

# Range of last cycle
indBeginLastCycle = argmin(abs.(t .- (t[end]-τ)))
rangeLastCycle = indBeginLastCycle:length(t)

# Aero states
χ_ = Array{Vector{Float64}}(undef,nTotalAeroStates)
for i in 1:nTotalAeroStates
    χ_[i] = [χ[tt][i] for tt in 1:length(t)]
end
plt_chi = plot(xlabel="Time [s]", ylabel="")
for i in 1:nTotalAeroStates
    plot!(t, χ_[i], c=colors[i], lw=lw, label="\$\\chi $(i)\$")
end
display(plt_chi)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_states.pdf"))

# Aero states' rates
χdot_ = Array{Vector{Float64}}(undef,nTotalAeroStates)
for i in 1:nTotalAeroStates
    χdot_[i] = [χdot[tt][i] for tt in 1:length(t)]
end
plt_chidot = plot(xlabel="Time [s]", ylabel="")
for i in 1:nTotalAeroStates
    plot!(t, χdot_[i], c=colors[i], lw=lw, label="\$\\dot{\\chi} $(i)\$")
end
display(plt_chidot)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_statesRates.pdf"))

# Pitch angle vs time
plt_alpha = plot(xlabel="Time [s]", ylabel="Pitch angle [deg]")
plot!(t, α*180/π, color=:black, lw=lw, label=false)
display(plt_alpha)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_alpha.pdf"))

# cn vs time
plt_cnt = plot(xlabel="Time [s]", ylabel="\$c_n\$")
plot!(t, cn, color=:black, lw=lw, label=false)
display(plt_cnt)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cnt.pdf"))

# cm vs time
plt_cmt = plot(xlabel="Time [s]", ylabel="\$c_m\$")
plot!(t, cm, color=:black, lw=lw, label=false)
display(plt_cmt)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cmt.pdf"))

# ct vs time
plt_ctt = plot(xlabel="Time [s]", ylabel="\$c_t\$")
plot!(t, ct, color=:black, lw=lw, label=false)
display(plt_ctt)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_ctt.pdf"))

# cl or cn vs α
if source == "NASA"
    plt_cla = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_l\$")
    plot!(α[rangeLastCycle]*180/π, cl[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
    scatter!(clRef[1,:], clRef[2,:], color=:black, ms=ms, label=refLabel)
    display(plt_cla)
    savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cla.pdf"))
elseif source == "GU"
    plt_cna = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_n\$")
    plot!(α[rangeLastCycle]*180/π, cn[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
    scatter!(cnRef[1,:], cnRef[2,:], color=:black, ms=ms, label=refLabel)
    display(plt_cna)
    savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cna.pdf"))
end

# cm vs α
plt_cma = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_m\$", legend=:bottomleft)
plot!(α[rangeLastCycle]*180/π, cm[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
scatter!(cmRef[1,:], cmRef[2,:], color=:black, ms=ms, label=refLabel)
display(plt_cma)
savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cma.pdf"))

# cd or ct vs α
if source == "NASA"
    plt_cda = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_d\$")
    plot!(α[rangeLastCycle]*180/π, cdrag[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
    scatter!(cdRef[1,:], cdRef[2,:], color=:black, ms=ms, label=refLabel)
    display(plt_cda)
    savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cda.pdf"))
elseif source == "GU"
    plt_cta = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_t\$")
    plot!(α[rangeLastCycle]*180/π, ct[rangeLastCycle], color=:black, lw=lw, label="AeroBeams")
    scatter!(ctRef[1,:], ctRef[2,:], color=:black, ms=ms, label=refLabel)
    display(plt_cta)
    savefig(string(absPath,"/pitchingAirfoilDSModelTest_",frameString,"_cta.pdf"))
end

println("Finished pitchingAirfoilDSModelTestPlotGenerator.jl")