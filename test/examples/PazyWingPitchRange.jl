using AeroBeams, LinearAlgebra, Plots, ColorSchemes, DelimitedFiles

# Aerodynamic solver
aeroSolver = Indicial()

# Airfoil
airfoil = deepcopy(flatPlate)

# Derivation method
derivationMethod = AD()

# Pazy wing
wing,L,nElem,chord,normSparPos,airfoil,surf = create_Pazy(aeroSolver=aeroSolver,airfoil=airfoil,derivationMethod=derivationMethod,p0=[0;-π/2;0])

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
PazyWingPitchRange = create_Model(name="PazyWingPitchRange",beams=[wing],BCs=[clamp],gravityVector=[0;0;-9.80665])

# Set NR system solver 
displayStatus = false
NR = create_NewtonRaphson(displayStatus=displayStatus)

# Set root angle and airspeed ranges, and initialize outputs
θRange = [3;5;7]
URange = collect(0:1:60)
tip_OOP = Array{Float64}(undef,length(θRange),length(URange))
tip_IP = Array{Float64}(undef,length(θRange),length(URange))
tip_twist = Array{Float64}(undef,length(θRange),length(URange))
tip_AoA = Array{Float64}(undef,length(θRange),length(URange))

# Sweep root angle
for (i,θ) in enumerate(θRange)
    # Update root angle on beam 
    wing.p0[3] = θ*π/180
    update_beam!(wing)
    # Sweep airspeed
    for (j,U) in enumerate(URange)
        # Display progress
        println("Solving for θ = $θ deg, U = $U m/s")
        # Set tip loss function at current airspeed and root angle
        surf.tipLossDecayFactor = tip_loss_factor_Pazy(θ,U)
        update_beam!(wing)
        # Update velocity of basis A (and update model)
        set_motion_basis_A!(model=PazyWingPitchRange,v_A=[0;U;0])
        # Create and solve problem
        global problem = create_SteadyProblem(model=PazyWingPitchRange,systemSolver=NR)
        solve!(problem)
        # @profview solve!(problem)
        # Get tip twist, AoA and OOP displacement at midchord
        tip_p = problem.nodalStatesOverσ[end][nElem].p_n2_b
        R,_ = rotation_tensor_WM(tip_p)
        Δ = R*[0; 1; 0]
        tip_twist[i,j] = asind(Δ[3])
        tip_OOP[i,j] = -(problem.nodalStatesOverσ[end][nElem].u_n2[1] - chord*(1/2-normSparPos)*sind(tip_twist[i,j]))
        tip_IP[i,j] = -problem.nodalStatesOverσ[end][nElem].u_n2[2]
        tip_AoA[i,j] = problem.model.elements[end].aero.flowAnglesAndRates.αₑ*180/π
    end
end

# Load reference data
tip_u3VsU_rootPitch5 = readdlm(string(pwd(),"/test/referenceData/Pazy/tip_u3VsU_rootPitch5.txt"))
tip_u3VsU_rootPitch7 = readdlm(string(pwd(),"/test/referenceData/Pazy/tip_u3VsU_rootPitch7.txt"))

# Plots
# ------------------------------------------------------------------------------
colors = get(colorschemes[:rainbow], LinRange(0, 1, length(θRange)))
lw = 2
ms = 3
# Deformed state of last problem
deformationPlot = plot_steady_deformation(problem,view=(45,30),save=true,savePath="/test/outputs/figures/PazyWingPitchRange/PazyWingPitchRange_deformation.pdf")
display(deformationPlot)
# Tip midchord OOP displacement vs. airspeed for root several pitch angles 
gr()
plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Tip OOP displacement [% semispan]", xlims=[0,60], ylims=[0,50])
plot!([NaN], [NaN], c=:black, lw=lw, label="AeroBeams")
scatter!([NaN], [NaN], c=:black, ms=ms, label="Exp.")
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_OOP[i,:]/L*100, c=colors[i], lw=lw, label="θ = $θ deg")
    if θ==5
        scatter!(tip_u3VsU_rootPitch5[1,:], tip_u3VsU_rootPitch5[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    elseif θ==7
        scatter!(tip_u3VsU_rootPitch7[1,:], tip_u3VsU_rootPitch7[2,:], mc=colors[i], ms=ms, msw=0, label=false)
    end
end
display(plt1)
savefig(string(pwd(),"/test/outputs/figures/PazyWingPitchRange/PazyWingPitchRange_tipOOP.pdf"))
# Tip twist vs. airspeed for root several pitch angles 
plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Tip twist [deg]", xlims=[0,60])
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_twist[i,:], c=colors[i], lw=lw, label="θ = $θ deg")
end
display(plt2)
savefig(string(pwd(),"/test/outputs/figures/PazyWingPitchRange/PazyWingPitchRange_tipTwist.pdf"))
# Tip AoA vs. airspeed for root several pitch angles 
plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Tip AoA [deg]", xlims=[0,60])
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_AoA[i,:], c=colors[i], lw=lw, label="θ = $θ deg")
end
display(plt3)
savefig(string(pwd(),"/test/outputs/figures/PazyWingPitchRange/PazyWingPitchRange_tipAoA.pdf"))
# Tip in-plane displacement vs. airspeed for root several pitch angles 
plt4 = plot(xlabel="Airspeed [m/s]", ylabel="Tip IP displacement [% semispan]", xlims=[0,60])
for (i,θ) in enumerate(θRange)
    plot!(URange, tip_IP[i,:]/L*100, c=colors[i], lw=lw, label="θ = $θ deg")
end
display(plt4)
savefig(string(pwd(),"/test/outputs/figures/PazyWingPitchRange/PazyWingPitchRange_tipIP.pdf"))
# Lift coefficient
plot_steady_outputs(problem,outputs=["cl"],colorScheme=:grays,lw=lw,save=true,saveFolder=string(relPath,"/"))

println("Finished PazyWingPitchRange.jl")