# Run the script
include("../examples/twoStoryFrame.jl")

# Reference frequencies (in Hz) by PETYT - Introduction to Finite Element Vibration Analysis - [2nd Ed.] (2010)
refFreqs = [11.8; 34.1]

# Display relative errors
ϵ_rel = freqs[[1,4]]/(2π)./refFreqs .- 1.0
println("Relative errors: $ϵ_rel")

# Plot mode shapes
relPath = "/test/outputs/figures/twoStoryFrame"
absPath = string(pwd(),relPath)
mkpath(absPath)

modesPlot = plot_mode_shapes(problem,scale=1,view=(45,30),legendPos=(0.3,0.1),frequencyLabel="frequency",save=true,savePath=string(relPath,"/twoStoryFrame_modeShapes.pdf"))
display(modesPlot)

println("Finished twoStoryFramePlotGenerator.jl")