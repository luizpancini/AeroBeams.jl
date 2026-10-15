using Plots, Polynomials, Interpolations, DelimitedFiles, ColorSchemes

# Span points of data collection
x_data = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
nSpanPoints = length(x_data)

# VLM data for normalized lift slope of undeformed wing
y_Λ0  = [5.59, 5.56, 5.46, 5.20, 4.56, 1.40]/(2π)/cosd(0)^2
y_Λ10 = [5.19, 5.29, 5.29, 5.11, 4.55, 1.41]/(2π)/cosd(10)^2
y_Λ20 = [4.52, 4.73, 4.80, 4.71, 4.26, 1.33]/(2π)/cosd(20)^2
y_Λ30 = [3.69, 3.94, 4.07, 4.04, 3.72, 1.18]/(2π)/cosd(30)^2
y = [y_Λ0,y_Λ10,y_Λ20,y_Λ30]

# Define ranges
ΛRange = [0,10,20,30]*π/180

# Order of polynomial
order = 8

# Define coefficient storage
polyCoeffs = zeros(order+1, length(ΛRange))

# Sweep conditions
for (i, Λ) in enumerate(ΛRange)
    # Lift slope at current condition
    y_data = y[i]
    # Fit polynomial
    p = fit(x_data, y_data, order)
    # Store coefficients
    for n=1:order+1
        polyCoeffs[n,i] = p.coeffs[n]
    end
end

# Save the coefficients
mkpath(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy")
writedlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/polyCoeffsUndef.txt", polyCoeffs)

# Construct polynomial at query point
function interpolate_polynomial(polyCoeffs::Matrix{Float64}, ΛRange::Vector{Float64}, Λq::Float64)
    order = size(polyCoeffs,1)-1
    q = zeros(order+1)
    for i=1:order+1
        itp = Interpolations.interpolate((ΛRange,), polyCoeffs[i,:], Gridded(Linear()))
        q[i] = itp(Λq)
    end
    return Polynomial(q)
end

# Interpolation polynomials
p = Array{Polynomial}(undef,length(ΛRange))
for (i, Λ) in enumerate(ΛRange)
    p[i] = interpolate_polynomial(polyCoeffs, ΛRange, Λ)
end

# Domain
x = 0:0.01:1

# Plot
colors = cgrad(:rainbow, length(ΛRange), categorical=true)
plt = plot(xlims=[0,1], ylims=[0,1], xlabel="Normalized span", ylabel="Normalized lift slope, \$C_{n_\\alpha}/(2\\pi)\$")
scatter!([NaN], [NaN], c=:black, label="Data points")
plot!([NaN], [NaN], lw=2, c=:black, label="Polynomial fit")
for (i, Λ) in enumerate(ΛRange)
    scatter!(x_data, y[i], c=colors[i], label=false)
    plot!(x,p[i].(x), c=colors[i], lw=2, label="\$\\Lambda=$(round(Int,Λ*180/π))^\\circ\$")
end
display(plt)