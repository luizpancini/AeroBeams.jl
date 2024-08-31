using Plots

# Run the script
include("../examples/PazyWingDARPAGust.jl")

# Set paths
relPath = "/test/outputs/figures/PazyWingDARPAGust"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Animation
plot_dynamic_deformation(problem,refBasis="A",plotFrequency=10,plotLimits=[(-L/2,L/2),(-L/2,L/2),(0,L)],save=true,savePath=string(relPath,"/PazyWingDARPAGust_deformation.gif"),displayProgress=true)

# Plot configurations
lw = 2
gr()

# Tip displacement
plt1 = plot(xlabel="Time [s]", ylabel="Tip OOP disp. [% semispan]")
plot!(t, tipOOP/L*100, color=:black, lw=lw, label=false)
display(plt1)
savefig(string(absPath,"/PazyWingDARPAGust_disp.pdf"))

# Tip AoA
plt2 = plot(xlabel="Time [s]", ylabel="Tip angle of attack [deg]")
plot!(t, tipAoA*180/π, color=:black, lw=lw, label=false)
display(plt2)
savefig(string(absPath,"/PazyWingDARPAGust_AoA.pdf"))

# 3/4-span cn
plt3 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_n\$")
plot!(t, tqSpan_cn, color=:black, lw=lw, label=false)
display(plt3)
savefig(string(absPath,"/PazyWingDARPAGust_cn.pdf"))

# 3/4-span cm
plt4 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_m\$")
plot!(t, tqSpan_cm, color=:black, lw=lw, label=false)
display(plt4)
savefig(string(absPath,"/PazyWingDARPAGust_cm.pdf"))

# 3/4-span ct
plt5 = plot(xlabel="Time [s]", ylabel="3/4-span \$c_t\$")
plot!(t, tqSpan_ct, color=:black, lw=lw, label=false)
display(plt5)
savefig(string(absPath,"/PazyWingDARPAGust_ct.pdf"))

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
savefig(string(absPath,"/PazyWingDARPAGust_states.pdf"))

# Gust velocity profile
s = collect(0:0.001:L)
time = collect(t0:0.01:t0+τ)
UGust(x,t) = t0<=t<=t0+τ ? (gust.UGustInertial([0;U*t;x],t))[1] : 0
X = [x for _ = time for x = s]
T = [t for t = time for _ = s]
Xnorm = X/L
Tnorm = @. (T-t0)/τ
Znorm = UGust.(X,T)/abs(Ug)
plt7 = surface(ylabel="Normalized time", xlabel="Normalized span", zlabel="Normalized gust velocity")
surface!(Xnorm,Tnorm,Znorm,xlims=[0,1],ylims=[0,1],zlims=[-1,1])
display(plt7)
savefig(string(absPath,"/PazyWingDARPAGust_gust.pdf"))

println("Finished PazyWingDARPAGustPlotGenerator.jl")