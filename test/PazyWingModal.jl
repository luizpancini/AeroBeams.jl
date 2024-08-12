using AeroBeams, LinearAlgebra

# Pazy wing
wing,L,_ = create_Pazy()

# BCs
clamp = create_BC(name="clamp",beam=wing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

# Model
PazyWingModal = create_Model(name="PazyWingModal",beams=[wing],BCs=[clamp],gravityVector=[0;0;-9.80665],units=create_UnitsSystem(frequency="Hz"))

# Number of modes
nModes = 5

# Wing position range (horizontal and vertical configurations)
θRange = [0,-π/2]

# Initialize outputs
problem = Array{EigenProblem}(undef,length(θRange))
freqs = Array{Vector{Float64}}(undef,length(θRange))

# Loop configurations
for (i,θ) in enumerate(θRange)
    # Set rotation and update model
    wing.p0 = [0; θ; 0]
    update_beam!(wing)
    update_model!(PazyWingModal)
    # Create and solve problem
    problem[i] = create_EigenProblem(model=PazyWingModal,nModes=nModes)
    solve!(problem[i])
    # Get frequencies in Hz
    freqs[i] = problem[i].frequenciesOscillatory/(2π)
    # Plot mode shapes
    mkpath(string(pwd(),"/test/outputs/figures/PazyWingModal"))
    savePath = string("/test/outputs/figures/PazyWingModal/PazyWingModal_",i,".pdf")
    modesPlot = plot_mode_shapes(problem[i],scale=0.1,view=(30,30),frequencyLabel="frequency",save=true,savePath=savePath)
    display(modesPlot)
end

# Load reference data
modalFreqsGVT = [4.39 29.80 41.00 82.50 NaN; 4.26 28.50 42.00 81.50 60.70]
modalFreqsUMNAST = [4.19 28.49 41.58 83.03 100.07; 4.19 28.49 41.88 83.06 105.89]

# Compare solutions
ϵ_rel_horz_GVT = freqs[1] ./ modalFreqsGVT[1,:] .- 1
ϵ_rel_vert_GVT = freqs[2] ./ modalFreqsGVT[2,:] .- 1
ϵ_rel_horz_UMNAST = freqs[1] ./ modalFreqsUMNAST[1,:] .- 1
ϵ_rel_vert_UMNAST = freqs[2] ./ modalFreqsUMNAST[2,:] .- 1

println("Relative erros for horizontal wing:")
println("GVT: $ϵ_rel_horz_GVT")
println("UMNAST: $ϵ_rel_horz_UMNAST")
println("Relative erros for vertical wing:")
println("GVT: $ϵ_rel_vert_GVT")
println("UMNAST: $ϵ_rel_vert_UMNAST")

println("Finished PazyWingModal.jl")