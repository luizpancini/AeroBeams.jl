using AeroBeams, JLD2, Plots, ColorSchemes

# Stiffness factor range to compare
λ = 1

# Bending pre-curvature range to compare
k2Range = -0.015:0.015:0.045

# Set paths
relPathFig = "/dev/cHALE/Flexible/outputs/figures/cHALE_compare_stability_disp_plots_k2Range"
absPathFig = string(pwd(),relPathFig)
mkpath(absPathFig)
relPathData = "/dev/cHALE/Flexible/outputs/data/cHALE_compare_stability/"
absPathData = string(pwd(),relPathData)
mkpath(absPathData)

# Plot configurations
colors = palette([:royalblue, :blueviolet, :deeppink, :darkorange, :gold])
ts = 10
fs = 16
lfs = 9
tsz = 10
lw = 2
ms = 3
msw = 0
L = 16
gr()

# Initialize plot
plt_u3 = plot(xlabel="Normalized spanwise direction", ylabel="Normalized OOP direction", xlims=[0,1], ylims=[-0.4,1], tickfont=font(ts), guidefont=font(fs), legendfontsize=lfs)
plot!([NaN], [NaN], c=:gray, lw=lw, ls=:dashdot, label="Undeformed")
plot!([NaN], [NaN], c=:gray, lw=lw, ls=:solid, shape=:circle, ms=2, msw=0, label="At flutter - aircraft")
plot!([NaN], [NaN], c=:gray, lw=lw, ls=:dash, label="At flutter - wing")

# Sweep bending pre-curvature
for (i,k2) in enumerate(k2Range)

    # Load data
    @load absPathData*string("x1_0_lambda",λ,"_k2",k2)*".jld2" x1_0
    @load absPathData*string("x3_0_lambda",λ,"_k2",k2)*".jld2" x3_0
    @load absPathData*string("URange_lambda",λ,"_k2",k2)*".jld2" URange
    @load absPathData*string("x1_def_lambda",λ,"_k2",k2)*".jld2" x1_def
    @load absPathData*string("x3_def_lambda",λ,"_k2",k2)*".jld2" x3_def
    @load absPathData*string("freqsAircraft_lambda",λ,"_k2",k2)*".jld2" freqsAircraft
    @load absPathData*string("dampsAircraft_lambda",λ,"_k2",k2)*".jld2" dampsAircraft
    @load absPathData*string("flutterSpeedAircraftAll_lambda",λ,"_k2",k2)*".jld2" flutterSpeedAircraftAll
    @load absPathData*string("flutterFreqAircraftAll_lambda",λ,"_k2",k2)*".jld2" flutterFreqAircraftAll
    @load absPathData*string("fluttertipOOPAircraftAll_lambda",λ,"_k2",k2)*".jld2" fluttertipOOPAircraftAll
    @load absPathData*string("flutterSpeedIndexAircraft_lambda",λ,"_k2",k2)*".jld2" flutterSpeedIndexAircraft
    @load absPathData*string("flutterModeAircraft_lambda",λ,"_k2",k2)*".jld2" flutterModeAircraft
    @load absPathData*string("flutterSpeedAircraft_lambda",λ,"_k2",k2)*".jld2" flutterSpeedAircraft
    @load absPathData*string("flutterFreqAircraft_lambda",λ,"_k2",k2)*".jld2" flutterFreqAircraft
    @load absPathData*string("flutterTipOOPAircraft_lambda",λ,"_k2",k2)*".jld2" flutterTipOOPAircraft
    @load absPathData*string("freqsWing_lambda",λ,"_k2",k2)*".jld2" freqsWing
    @load absPathData*string("dampsWing_lambda",λ,"_k2",k2)*".jld2" dampsWing
    @load absPathData*string("flutterSpeedWingAll_lambda",λ,"_k2",k2)*".jld2" flutterSpeedWingAll
    @load absPathData*string("flutterFreqWingAll_lambda",λ,"_k2",k2)*".jld2" flutterFreqWingAll
    @load absPathData*string("fluttertipOOPWingAll_lambda",λ,"_k2",k2)*".jld2" fluttertipOOPWingAll
    @load absPathData*string("flutterSpeedIndexWing_lambda",λ,"_k2",k2)*".jld2" flutterSpeedIndexWing
    @load absPathData*string("flutterModeWing_lambda",λ,"_k2",k2)*".jld2" flutterModeWing
    @load absPathData*string("flutterSpeedWing_lambda",λ,"_k2",k2)*".jld2" flutterSpeedWing
    @load absPathData*string("flutterFreqWing_lambda",λ,"_k2",k2)*".jld2" flutterFreqWing
    @load absPathData*string("flutterTipOOPWing_lambda",λ,"_k2",k2)*".jld2" flutterTipOOPWing

    # Plot undeformed and deformed shape at flutter
    plot!(x1_0/L, x3_0/L, c=colors[i], lw=lw, ls=:dashdot, label=false)
    plot!(x1_def[flutterSpeedIndexAircraft]/L, x3_def[flutterSpeedIndexAircraft]/L, c=colors[i], lw=lw, ls=:solid, shape=:circle, ms=2, msw=0, label=string("\$k_2=\$",k2))
    plot!(x1_def[flutterSpeedIndexWing]/L, x3_def[flutterSpeedIndexWing]/L, c=colors[i], lw=lw, ls=:dash, label=false)

end
display(plt_u3)
savefig(string(absPathFig,"/cHALE_compare_stability_disp_plots_k2Range_lambda_",λ,".pdf"))
