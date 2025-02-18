using Plots

# Set example parameters
α1₀ = 15*π/180
fb = 0.7
f₀ = 0.01
S1 = 4*π/180
S2 = 2*π/180
fb_smooth = (S2-S1*f₀)/(S1+S2)
δ = 5e-3

# Angle of attack range
α = vcat(-20:0.1:20)*π/180

# Separation point expressions
f(α) = abs(α) <= α1₀ ? 1-(1-fb)*exp((abs(α)-α1₀)/S1) : f₀+(fb-f₀)*exp((α1₀-abs(α))/S2)
f_smooth(α) = abs(α) <= α1₀ ? 1-(1-fb_smooth)*exp((abs(α)-α1₀)/S1) : f₀+(fb_smooth-f₀)*exp((α1₀-abs(α))/S2)
f_low(α) = 1-(1-fb)*exp((abs(α)-α1₀)/S1)
f_high(α) = f₀+(fb-f₀)*exp((α1₀-abs(α))/S2)
f_tanh(α) = f_low(α) + (f_high(α)-f_low(α)) * (1+tanh((abs(α)-α1₀)/δ))/2

# Plot
lw = 2

plt = plot(xlabel="\$\\alpha\$ [deg]", ylabel="\$f\$", xlims=[-20,20], ylims=[0,1], legend_position=:bottom)
plot!(α*180/π, f.(α), c=:black, lw=lw, label="Standard")
plot!(α*180/π, f_smooth.(α), c=:blue, lw=lw, label="Smooth")
plot!(α*180/π, f_tanh.(α), c=:red, lw=lw, label="tanh()")
display(plt)