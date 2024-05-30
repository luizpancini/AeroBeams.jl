using Plots
gr()
histogram2d(randn(10000), randn(10000), nbins = 20, color=:rainbow,
colorbar_title="AIRSPEED")