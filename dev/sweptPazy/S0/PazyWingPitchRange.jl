using AeroBeams, DelimitedFiles, Plots, ColorSchemes

# Aerodynamic solver
aeroSolver = Inflow()

# Airfoil section
airfoil = deepcopy(NACA0018)

# Flag for upright position
upright = true

# Gravity
g = 9.80665

# Flag for small angles approximation
smallAngles = false

# Fixed geometrical and discretization properties
nElem,L,chord,normSparPos = geometrical_properties_Pazy()

# Tip mass (for test 1, 0.01 kg, 40 mm behind the trailing-edge)
tipMass = 0.01
ηtipMass = [0; -chord*(1-normSparPos)-0.04; 0]

# Root angle (in degrees) and airspeed ranges
θRange = [3, 5, 7]
URange = collect(0:1:60)

# Initialize outputs
problem = Array{SteadyProblem}(undef,length(θRange),length(URange))
tip_OOP = Array{Float64}(undef,length(θRange),length(URange))
tip_IP = Array{Float64}(undef,length(θRange),length(URange))
tip_twist = Array{Float64}(undef,length(θRange),length(URange))
tip_AoA = Array{Float64}(undef,length(θRange),length(URange))

# Sweep root angle
for (i,θ) in enumerate(θRange)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $θ deg, U = $U m/s")
        # Update model
        PazyWingPitchRange,_ = create_Pazy(aeroSolver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,upright=upright,θ=θ*π/180,airspeed=U,g=g,tipMass=tipMass,ηtipMass=ηtipMass,smallAngles=smallAngles)
        # Create and solve problem
        problem[i,j] = create_SteadyProblem(model=PazyWingPitchRange)
        solve!(problem[i,j])
        # Get tip twist, AoA, IP and OOP displacement at beam reference line
        tip_p = problem[i,j].nodalStatesOverσ[end][nElem].p_n2_b
        R,_ = rotation_tensor_WM(tip_p)
        Δ = R*[0; 1; 0]
        tip_twist[i,j] = asind(Δ[3])
        tip_OOP[i,j] = problem[i,j].nodalStatesOverσ[end][nElem].u_n2_b[3]
        tip_IP[i,j] = -problem[i,j].nodalStatesOverσ[end][nElem].u_n2_b[2]
        tip_AoA[i,j] = problem[i,j].model.elements[end].aero.flowAnglesAndRates.αₑ*180/π
    end
end

# Load reference data
tip_u3VsU_rootPitch5_Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3VsU_rootPitch5_Exp.txt")
tip_u3VsU_rootPitch5_UMNAST = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3VsU_rootPitch5_UMNAST.txt")
tip_u3VsU_rootPitch7_Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3VsU_rootPitch7_Exp.txt")
tip_u3VsU_rootPitch7_UMNAST = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3VsU_rootPitch7_UMNAST.txt")
tip_thetaVsU_rootPitch5_Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_thetaVsU_rootPitch5_Exp.txt")
tip_thetaVsU_rootPitch5_UMNAST = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_thetaVsU_rootPitch5_UMNAST.txt")
tip_thetaVsU_rootPitch7_Exp = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_thetaVsU_rootPitch7_Exp.txt")
tip_thetaVsU_rootPitch7_UMNAST = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_thetaVsU_rootPitch7_UMNAST.txt")
tip_u3Vsq_rootPitch3_ExpTest1 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3Vsq_rootPitch3_ExpTest1.txt")
tip_u3Vsq_rootPitch5_ExpTest1 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3Vsq_rootPitch5_ExpTest1.txt")
tip_u3Vsq_rootPitch7_ExpTest1 = readdlm(pkgdir(AeroBeams)*"/test/referenceData/Pazy/tip_u3Vsq_rootPitch7_ExpTest1.txt")
ρTest2 = 2*1050/43^2 # Estimated air density for test 2, from Table 13 of Avin et al.'s paper

# Set paths
relPath = "/dev/sweptPazy/S0/outputs/PazyWingPitchRange"
absPath = string(pwd(),relPath)
mkpath(absPath)

# Deformed state of last problem
deformationPlot = plot_steady_deformation(problem,view=(45,30),save=true,savePath="/dev/sweptPazy/S0/outputs/PazyWingPitchRange/PazyWingPitchRange_deformation.pdf")
display(deformationPlot)

# Plot configurations
colors = cgrad(:rainbow, length(θRange), categorical=true)
ts = 10
fs = 16
lfs = 9
lw = 2
ms = 3
msw = 0
gr()

# Tip OOP displacement vs. airspeed for root several pitch angles 
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP disp. [% semispan]", xlims=[0,60], ylims=[0,50], tickfont=font(ts), guidefont=font(fs), legend=:topleft, legendfontsize=lfs)
plot!([NaN], [NaN], c=:black, lw=lw, ls=:solid, label="AeroBeams")
plot!([NaN], [NaN], c=:black, lw=lw, ls=:dashdot, label="UM/NAST (Exp. loss) - Riso & Cesnik (2022,2023)")
scatter!([NaN], [NaN], c=:black, ms=ms, label="Exp. (test 1) - Avin et al. (2022)")
scatter!([NaN], [NaN], c=:black, shape=:diamond, ms=ms, label="Exp. (test 2) - Avin et al. (2022)")
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_OOP[i,:]/L*100, c=colors[i], lw=lw, ls=:solid, label=string("\$\\alpha_r=",round(Int,θ),"^\\circ\$"))
    if θ==3
        scatter!(sqrt.(2*tip_u3Vsq_rootPitch3_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch3_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    elseif θ==5
        plot!(tip_u3VsU_rootPitch5_UMNAST[1,:], tip_u3VsU_rootPitch5_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_u3VsU_rootPitch5_Exp[1,:], tip_u3VsU_rootPitch5_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        scatter!(sqrt.(2*tip_u3Vsq_rootPitch5_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch5_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    elseif θ==7
        plot!(tip_u3VsU_rootPitch7_UMNAST[1,:], tip_u3VsU_rootPitch7_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_u3VsU_rootPitch7_Exp[1,:], tip_u3VsU_rootPitch7_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
        scatter!(sqrt.(2*tip_u3Vsq_rootPitch7_ExpTest1[1,:]/ρTest2), tip_u3Vsq_rootPitch7_ExpTest1[2,:]*100/L, mc=colors[i], shape=:diamond, ms=ms, msw=msw, label=false)
    end
end
display(plt1)
savefig(string(absPath,"/PazyWingPitchRange_tipOOP.pdf"))

# Tip twist vs. airspeed for root several pitch angles 
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_twist[i,:], c=colors[i], lw=lw, label=false)
    if θ==5
        plot!(tip_thetaVsU_rootPitch5_UMNAST[1,:], tip_thetaVsU_rootPitch5_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_thetaVsU_rootPitch5_Exp[1,:], tip_thetaVsU_rootPitch5_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
    elseif θ==7
        plot!(tip_thetaVsU_rootPitch7_UMNAST[1,:], tip_thetaVsU_rootPitch7_UMNAST[2,:], lw=lw, ls=:dashdot, c=colors[i], label=false)
        scatter!(tip_thetaVsU_rootPitch7_Exp[1,:], tip_thetaVsU_rootPitch7_Exp[2,:], mc=colors[i], ms=ms, msw=msw, label=false)
    end
end
display(plt2)
savefig(string(absPath,"/PazyWingPitchRange_tipTwist.pdf"))

# Tip AoA vs. airspeed for root several pitch angles 
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Tip AoA [deg]", xlims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_AoA[i,:], c=colors[i], lw=lw, label=false)
end
display(plt3)
savefig(string(absPath,"/PazyWingPitchRange_tipAoA.pdf"))

# Tip in-plane displacement vs. airspeed for root several pitch angles 
plt4 = plot(xlabel="Airspeed [m/s]", ylabel="Tip IP disp. [% semispan]", xlims=[0,60], tickfont=font(ts), guidefont=font(fs))
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_IP[i,:]/L*100, c=colors[i], lw=lw, label=false)
end
display(plt4)
savefig(string(absPath,"/PazyWingPitchRange_tipIP.pdf"))

# Lift coefficient over span
plot_steady_outputs(problem,outputs=["cl","κ2","F3"],colorScheme=:grays,lw=lw,save=true,saveFolder=string(relPath,"/"))

println("Finished PazyWingPitchRange.jl")