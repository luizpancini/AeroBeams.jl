using Plots

# Run the script
include("../examples/PazyWingOMCGust.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingOMCGust"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,plotLimits=([-L/2,L/2],[-L/2,L/2],[0,L]),save=true,savePath=string(relPath,"/PazyWingOMCGust_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Tip displacement
plt1 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, tipOOP/L*100, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/PazyWingOMCGust_disp.pdf"))

# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingOMCGust_AoA.pdf"))

# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/PazyWingOMCGust_cn.pdf"))

# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/PazyWingOMCGust_cm.pdf"))

# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/PazyWingOMCGust_ct.pdf"))

# Aero states at 3/4-span
nAeroStates = problem.model.elements[1].aero.nTotalAeroStates
colors = get(colorschemes[:rainbow], LinRange(0, 1, nAeroStates))
tqsχ_ = Array{Vector{Float64}}(undef,nAeroStates)
for i in 1:nAeroStates
    tqsχ_[i] = [tqsχ[tt][i] for tt in 1:length(t)]
end
plt6 = plot(xlabel="Time [s]", ylabel="")
for i in 1:nAeroStates
    plot!(t, tqsχ_[i], c=colors[i], lw=lw, label="\$\\chi $(i)\$")
end
display(plt6)
savefig(string(absPath,"/PazyWingOMCGust_states.pdf"))

println("Finished PazyWingOMCGustPlotGenerator.jl")