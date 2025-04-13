using AeroBeams, Plots, DelimitedFiles, Interpolations, ColorSchemes

function VLM_tip_loss(type,Λ,θ,U)

    # Sweep angle, root pitch angle and airspeed ranges of the reference data
    ΛRange = [0,10,20,30]*π/180
    θRange = type == "VLM-undef" ? nothing : [0,3,5,7]*π/180
    URange = type == "VLM-undef" ? nothing : vcat(0:10:70,120)

    # Polynomial order
    order = 8

    # Load polynomial coefficients at reference data points
    polyCoeffs = type == "VLM-undef" ? zeros(order+1, length(ΛRange)) : zeros(order+1, length(ΛRange), length(θRange), length(URange))
    if type == "VLM-undef"
        polyCoeffs = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/polyCoeffsUndef.txt")
    elseif type == "VLM-def"
        for j=1:length(θRange)
            for k=1:length(URange)
                polyCoeffs[:,:,j,k] = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/polyCoeffsDef_"*string(j)*"_"*string(k)*".txt")
            end
        end
    end

    # Interpolate polynomial coefficients
    q = zeros(9)
    for i=1:9
        if type == "VLM-undef"
            itp = Interpolations.interpolate((ΛRange,), polyCoeffs[i,:], Gridded(Linear()))
            q[i] = itp(max(min(Λ,ΛRange[end]),ΛRange[1]))*sqrt(1-(U/343)^2) # Correct with compressibility factor
        elseif type == "VLM-def"
            itp = Interpolations.interpolate((ΛRange, θRange, URange), polyCoeffs[i,:,:,:], Gridded(Linear()))
            q[i] = itp(max(min(Λ,ΛRange[end]),ΛRange[1]), max(min(θ,θRange[end]),θRange[1]), min(U,URange[end]))
        end
    end

    # Tip loss function
    return s -> q[1]+q[2]*s+q[3]*s^2+q[4]*s^3+q[5]*s^4+q[6]*s^5+q[7]*s^6+q[8]*s^7+q[9]*s^8
end

# Set case
case = 2

# Airspeed range
URange = vcat(10:10:60)

# Sweep angle
Λ = case == 1 ? 0*π/180 : 20*π/180

# Pitch angle
θ = 5*π/180

# Domain
s = 0:0.01:1

# Load reference data (from AePW4 meetings)
clalapha_Lambda0_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda0_aoa5_NastranTechnion.txt")
clalapha_Lambda20_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda20_aoa5_NastranTechnion.txt")

ylim = 0
refData = 0
if Λ == 0 && θ == 5*π/180
    refData = clalapha_Lambda0_aoa5_NastranTechnion
    ylim = 7.5
elseif Λ == 20*π/180 && θ == 5*π/180
    refData = clalapha_Lambda20_aoa5_NastranTechnion
    ylim = 6
end

# Set paths
relPath = "/dev/misc/outputs/figures/tipLossComparisonVLM_Nastran"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot normal force distribution at the specified condition
colors = cgrad(:rainbow, length(URange), categorical=true)
plt = plot(xlims=[0,1], ylims=[0,ylim], xlabel="Normalized span", ylabel="\$C_{n_\\alpha}\$ [1/rad]", title=string("\$\\Lambda=",round(Int,Λ*180/pi),"^\\circ\$",", \$\\alpha_r=",round(Int,θ*180/pi),"^\\circ\$"), yticks=0:1:ylim, tickfont=font(10), guidefont=font(16), legendfontsize=11, legend=:bottomleft)
plot!([NaN], [NaN], c=:black, lw=2, label="Present VLM")
scatter!([NaN], [NaN], c=:black, ms=4, msw=0, label="Nastran (Technion)")
for (i,U) in enumerate(URange)
    β = sqrt(1-(U/343)^2)
    W = VLM_tip_loss("VLM-def",Λ,θ,U)
    cnαVLMdef = 2π/β * W.(s)
    plot!(s, cnαVLMdef, c=colors[i], lw=2, label="\$U=$U\$ m/s")
    j = Int(U/10)
    scatter!(refData[2*j-1,:],refData[2*j,:], c=colors[i], ms=4, msw=0, label=false)
end
if case != 1
    plot!(legend=false)
end
display(plt)
savefig(string(absPath,"/cnalpha_compare_VLM_Nastran_case",case,".pdf"))