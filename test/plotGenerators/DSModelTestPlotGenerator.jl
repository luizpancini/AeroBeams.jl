using Plots, ColorSchemes

# Run the script
include("../examples/DSModelTest.jl")

# Set paths
relPath = "/test/outputs/figures/DSModelTest"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,view=(30,30),plotLimits=[(0,L),(-L/2,L/2),(-L/2,L/2)],save=true,savePath=string(relPath,"/DSModelTest_deformation.gif"),displayProgress=true)

# Plot configurations
nTotalAeroStates = problem.model.elements[1].aero.nTotalAeroStates
colors = get(colorschemes[:rainbow], LinRange(0, 1, nTotalAeroStates))
lw = 2
ms = 3
gr()

# Pitch angle
plt1 = plot(xlabel="Time [s]", ylabel="Pitch angle [deg]")
plot!(t, α*180/π, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/DSModelTest_alpha.pdf"))

# Normal relative wind acceleration
plt2 = plot(xlabel="Time [s]", ylabel="\$\\dot{V}_3\$ [m/s^2]")
plot!(t, Vdot3, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/DSModelTest_Vdot3.pdf"))

# cn vs time
plt3 = plot(xlabel="Time [s]", ylabel="\$c_n\$")
plot!(t, cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/DSModelTest_cnt.pdf"))

# cm vs time
plt4 = plot(xlabel="Time [s]", ylabel="\$c_m\$")
plot!(t, cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/DSModelTest_cmt.pdf"))

# ct vs time
plt5 = plot(xlabel="Time [s]", ylabel="\$c_t\$")
plot!(t, ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/DSModelTest_ctt.pdf"))

# cl vs α
plt6 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_l\$")
plot!(α*180/π, cl, color=:black, lw=lw, label="AeroBeams")
scatter!(clRef[1,:], clRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
display(plt6)
savefig(string(absPath,"/DSModelTest_cla.pdf"))

# cm vs α
plt7 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_m\$")
plot!(α*180/π, cm, color=:black, lw=lw, label="AeroBeams")
scatter!(cmRef[1,:], cmRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
display(plt7)
savefig(string(absPath,"/DSModelTest_cma.pdf"))

# cd vs α
plt8 = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$c_d\$")
plot!(α*180/π, cd, color=:black, lw=lw, label="AeroBeams")
scatter!(cdRef[1,:], cdRef[2,:], color=:black, ms=ms, label="Exp. McAlister et al (1982)")
display(plt8)
savefig(string(absPath,"/DSModelTest_cda.pdf"))

# Aero states at 3/4-span
χ_ = Array{Vector{Float64}}(undef,nTotalAeroStates)
for i in 1:nTotalAeroStates
    χ_[i] = [χ[tt][i] for tt in 1:length(t)]
end
plt9 = plot(xlabel="Time [s]", ylabel="")
for i in 1:nTotalAeroStates
    plot!(t, χ_[i], c=colors[i], lw=lw, label="\$\\chi $(i)\$")
end
display(plt9)
savefig(string(absPath,"/DSModelTest_states.pdf"))

# Aero states' rates at 3/4-span
χdot_ = Array{Vector{Float64}}(undef,nTotalAeroStates)
for i in 1:nTotalAeroStates
    χdot_[i] = [χdot[tt][i] for tt in 1:length(t)]
end
plt10 = plot(xlabel="Time [s]", ylabel="")
for i in 1:nTotalAeroStates
    plot!(t, χdot_[i], c=colors[i], lw=lw, label="\$\\dot{\\chi} $(i)\$")
end
display(plt10)
savefig(string(absPath,"/DSModelTest_statesRates.pdf"))

println("Finished DSModelTestPlotGenerator.jl")