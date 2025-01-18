using LinearAlgebra, Plots

# Set unbounded rotation angle range
ϕRange = vcat(-4π:0.01:4π)

# Arrays for the generating function value of the unscaled, scaled and extended set of parameters
pu = Vector{Float64}(undef,length(ϕRange))
ps = Vector{Float64}(undef,length(ϕRange))
pe = Vector{Float64}(undef,length(ϕRange))

# Loop rotation angle
for (i,ϕ) in enumerate(ϕRange)
    # Unscaled parameters
    pu[i] = 4*tan(ϕ/4)
    # Scaled parameters (apply equation only if abs(pu[i]) > 4)
    ps[i] = abs(pu[i]) > 4 ? -4*cot(ϕ/4) : pu[i]
    # Extended parameters
    q = round(ϕ/(2π))
    pe[i] = ps[i] + 8*q # or pe[i] = 4*tan((ϕ-2π*q)/4) + 8*q
end

# Set paths
relPath = "/dev/misc/outputs/figures/WienerMilenkovicParameters"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
lw = 2
gr()

# Plot
plt = plot(xlabel="\\phi/2\\pi",ylabel="Generating function value", xlims=[ϕRange[1],ϕRange[end]+1e-2]/(2π), ylims=[-16,16], xticks=-2:0.5:2, yticks=-20:4:20, legend=:top)
plot!(ϕRange/(2π),pu, lw=lw, c=:red, label="Unscaled")
plot!(ϕRange/(2π),ps, lw=lw, c=:black, label="Scaled")
plot!(ϕRange/(2π),pe, lw=lw, c=:green, label="Extended")
display(plt)
savefig(string(absPath,"/WMParameters.pdf"))