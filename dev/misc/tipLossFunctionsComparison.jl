using AeroBeams, Plots, DelimitedFiles, Interpolations, ColorSchemes

function exponential_loss(θ,U)

    # Bound inputs
    θ = min(7*π/180,max(θ,0))
    U = min(60,U)

    # Coefficients as a function of root angle of attack and airspeed
    θRange = [0; 1; 2; 3; 4; 5; 6; 7]*π/180
    τ₀Range = [6.58; 6.29; 6.00; 5.92; 5.91; 6.08; 6.43; 6.88]
    τ₁Range = 1e-2*[0; 3.33; 5.19; 6.01; 6.49; 5.98; 4.43; 2.30]
    τ₂Range = -1e-4*[0; 5.56; 8.59; 10.6; 12.3; 12.8; 11.9; 10.2]
    τ₀ = LinearInterpolations.interpolate(θRange,τ₀Range,θ)
    τ₁ = LinearInterpolations.interpolate(θRange,τ₁Range,θ)
    τ₂ = LinearInterpolations.interpolate(θRange,τ₂Range,θ)    
    τ = τ₀ + τ₁*U + τ₂*U^2

    # Tip loss function
    return s -> 1-exp(-τ*(1-s))
end

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

# Sweep angle
Λ = 0*π/180

# Pitch angle
θ = 0*π/180

# Airspeed
U = 100

# Compressibility factor
M = U/343
β = sqrt(1-M^2)

# Domain
s = 0:0.01:1

# No tip loss
WConst = s -> 1

# Exponential loss
WExp = exponential_loss(θ,U)

# VLM - undeformed shape
WVLMundef = VLM_tip_loss("VLM-undef",Λ,θ,U)

# VLM - deformed shape
WVLMdef = VLM_tip_loss("VLM-def",Λ,θ,U)

# Load reference data (from AePW4 meetings)
clalapha_Lambda0_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda0_aoa5_NastranTechnion.txt")
clalapha_Lambda20_aoa5_NastranTechnion = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/clalapha_Lambda20_aoa5_NastranTechnion.txt")

# Set paths
relPath = "/dev/misc/outputs/figures/tipLossFunctionsComparison"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Plot normal force distribution at the specified condition
colors = cgrad(:rainbow, 4, categorical=true)
plt = plot(xlims=[0,1], ylims=[0,1.25], xlabel="Normalized span", ylabel="\$ W(s) = C_{n_\\alpha}/(2\\pi/\\beta)\$", title=string("\$\\Lambda=",round(Int,Λ*180/pi),"^\\circ\$",", \$\\alpha_r=",round(Int,θ*180/pi),"^\\circ\$",", \$U=",U," \$ m/s"), tickfont=font(10), guidefont=font(16), legendfontsize=12)
plot!(s, WConst.(s), c=colors[1], lw=2, ls=:solid, label="No loss")
plot!(s, WExp.(s), c=colors[2], lw=2, ls=:dash, label="Exponential loss")
plot!(s, WVLMundef.(s), c=colors[3], lw=2, ls=:dot, label="VLM - undeformed")
plot!(s, WVLMdef.(s), c=colors[4], lw=2, ls=:dashdot, label="VLM - deformed")
if Λ == 0 && θ == 5*π/180
    i = Int(U/10)
    scatter!(clalapha_Lambda0_aoa5_NastranTechnion[2*i-1,:],clalapha_Lambda0_aoa5_NastranTechnion[2*i,:]/(2π/β), c=:black, ms=4, msw=0, label="Nastran (Technion)")
elseif Λ == 20*π/180 && θ == 5*π/180
    i = Int(U/10)
    scatter!(clalapha_Lambda20_aoa5_NastranTechnion[2*i-1,:],clalapha_Lambda20_aoa5_NastranTechnion[2*i,:]/(2π/β), c=:black, ms=4, msw=0, label="Nastran (Technion)")
end
display(plt)
savefig(string(absPath,"/cnalpha.pdf"))
