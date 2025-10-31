using Plots, Measures

# Set paths
relPath = "/dev/helios/figures/heliosElevatorProfile"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot configurations
gr()
fs = 30

# Fictitious data
t = 0:1e-3:10
trimδ = 5
Δδ = 10
tδinit = 1
tδramp = 1
tδpeak = tδinit+tδramp
tδfinal = tδpeak+tδramp

# Flap elevator deflection
δ = t -> ifelse(
    t <= tδinit, 
    trimδ,
    ifelse(
        t <= tδpeak, 
        trimδ + Δδ * ((t-tδinit) / (tδpeak-tδinit)),
        ifelse(
            t <= tδfinal, 
            trimδ + Δδ - Δδ * ((t-tδpeak) / (tδfinal-tδpeak)),
            trimδ
        )
    )
)

# Plot
plt_δ = plot(xlabel="\$ t \$", ylabel="\$\\delta_f\$", xlims=extrema(t), ylims=[0,trimδ+Δδ], xticks=false, yticks=false, guidefont=font(fs), grid=false, left_margin=5mm)
plot!(t, δ.(t), c=:black, lw=2, label=false)
display(plt_δ)
savefig(plt_δ,string(absPath,"/heliosElevatorProfile.pdf"))