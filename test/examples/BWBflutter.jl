# # Flutter of a Blended-Wing-Body
# This example illustrates how to set up a flutter analysis of an aircraft in free flight. For that we take a Blended-Wing-Body (BWB) vehicle, a swept flying-wing with trailing-edge control surfaces. This aircraft model was described in [Weihua Su's PhD thesis](https://www.proquest.com/docview/304572531?pq-origsite=gscholar&fromopenview=true&sourcetype=Dissertations%20&%20Theses):

#md # ![](../assets/BWB.png)
#md # *BWB model geometry* by [Su](https://www.proquest.com/docview/304572531?pq-origsite=gscholar&fromopenview=true&sourcetype=Dissertations%20&%20Theses)

#md # ![](../assets/BWBbodyprops.png)
#md # *BWB's body properties* by [Su](https://www.proquest.com/docview/304572531?pq-origsite=gscholar&fromopenview=true&sourcetype=Dissertations%20&%20Theses)

#md # ![](../assets/BWBwingprops.png)
#md # *BWB's wing properties* by [Su](https://www.proquest.com/docview/304572531?pq-origsite=gscholar&fromopenview=true&sourcetype=Dissertations%20&%20Theses)

#md # !!! tip
#md #     The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/BWBflutter.jl).

# ### Problem setup
# Let's begin by setting up the variables of our problem.
using AeroBeams, DelimitedFiles, LinearInterpolations

## Aerodynamic solver
aeroSolver = Indicial()

## Flight altitude
h = 0e3

## Airspeed range
URange = collect(30:5:160)

## Number of vibration modes
nModes = 8

## Pre-allocate memory and initialize output arrays
trimProblem = Array{TrimProblem}(undef,length(URange))
trimProblemSpringed = Array{TrimProblem}(undef,length(URange))
eigenProblem = Array{EigenProblem}(undef,length(URange))
trimAoA = Array{Float64}(undef,length(URange))
trimThrust = Array{Float64}(undef,length(URange))
trimδ = Array{Float64}(undef,length(URange))
untrackedFreqs = Array{Vector{Float64}}(undef,length(URange))
untrackedDamps = Array{Vector{Float64}}(undef,length(URange))
untrackedEigenvectors = Array{Matrix{ComplexF64}}(undef,length(URange))
freqs = Array{Vector{Float64}}(undef,length(URange))
damps = Array{Vector{Float64}}(undef,length(URange))
#md nothing #hide

# For the trim problem, we set a Newton-Raphson solver for the system of equations, with the adequate relaxation factor for trim problems (`relaxFactor = 0.5`), and an increased number of maximum iterations (`maxIter = 50`, the default is 20).
## System solver
relaxFactor = 0.5
maxIter = 50
NR = create_NewtonRaphson(ρ=relaxFactor,maximumIterations=maxIter)
#md nothing #hide

# Next, we address an important step to be taken when performing flutter analyses in free flight with AeroBeams: to attach the model to light springs in displacement and rotation. This step is necessary for the solver to find the flight dynamic (rigid-body) modes of the vehicle, by introducing some sensitivity of the finite element states to those degrees-of-freedom. In the present case, we will attach two springs to the vehicle, one at each of the nodes where the transition from the body (fuselage) to the wing begins. An appropriate value for the stiffness of the springs is specified by the variable `μ`. The vectors `ku` and `kp` denote the stiffness values in the three orthogonal directions.
## Attachment springs
μ = 1e-1
ku = kp = μ*[1; 1; 1]
spring1 = create_Spring(elementsIDs=[1],nodesSides=[1],ku=ku,kp=kp)
spring2 = create_Spring(elementsIDs=[3],nodesSides=[2],ku=ku,kp=kp)
#md nothing #hide

# ### Problem solution
# At this point we can sweep the airspeed vector to find the solution at each value.
## Sweep airspeed range
for (i,U) in enumerate(URange)
    println("Solving for U = $U m/s") #src
    ## The first step of the solution is to trim the aircraft at that flight condition (combination of altitude and airspeed). We leverage the built-in function in AeroBeams to create our model for the trim problem.
    BWBtrim = first(create_BWB(aeroSolver=aeroSolver,altitude=h,airspeed=U,δElevIsTrimVariable=true,thrustIsTrimVariable=true))

    ## Let's set an initial guess for the trim problem as the previously converged solution, if available
    x0Trim = i == 1 ? zeros(0) : trimProblem[i-1].x

    ## Now we create and solve the trim problem.
    trimProblem[i] = create_TrimProblem(model=BWBtrim,systemSolver=NR,x0=x0Trim)
    solve!(trimProblem[i])

    ## We extract the trim variables at the current airspeed and set them into our pre-allocated arrays. The trimmed angle of attack at the root, `trimAoA[i]`, is not necessary for the flutter analyses, it is merely an output of interest.
    trimAoA[i] = trimProblem[i].aeroVariablesOverσ[end][BWBtrim.beams[3].elementRange[1]].flowAnglesAndRates.αₑ
    trimThrust[i] = trimProblem[i].x[end-1]*BWBtrim.forceScaling 
    trimδ[i] = trimProblem[i].x[end]
    println("Trim AoA = $(round(trimAoA[i]*180/π,digits=2)) deg, trim thrust = $(round(trimThrust[i],digits=2)) N, trim δ = $(round(trimδ[i]*180/π,digits=2)) deg") #src

    ## Now let's create a new model, add springs to it, and update
    BWBtrimSpringed = first(create_BWB(aeroSolver=aeroSolver,altitude=h,airspeed=U,δElevIsTrimVariable=true,thrustIsTrimVariable=true))
    add_springs_to_beam!(beam=BWBtrimSpringed.beams[2],springs=[spring1])
    add_springs_to_beam!(beam=BWBtrimSpringed.beams[3],springs=[spring2])
    update_model!(BWBtrimSpringed)

    ## Let's set an initial guess for the springed trim problem as the previously converged solution, if available
    x0TrimSpringed = i == 1 ? zeros(0) : trimProblemSpringed[i-1].x

    ## Now we create and solve the trim problem for the springed model.
    trimProblemSpringed[i] = create_TrimProblem(model=BWBtrimSpringed,systemSolver=NR,x0=x0TrimSpringed)
    solve!(trimProblemSpringed[i])

    ## Let's make sure that the trim outputs of the springed model are very close to the one without springs, that is, the springs did not affect the trim solution. Notice that all ratios are very close to 1.
    trimAoASpringed = trimProblemSpringed[i].aeroVariablesOverσ[end][BWBtrimSpringed.beams[3].elementRange[1]].flowAnglesAndRates.αₑ
    trimThrustSpringed = trimProblemSpringed[i].x[end-1]*BWBtrimSpringed.forceScaling 
    trimδSpringed = trimProblemSpringed[i].x[end]
    println("Trim outputs' ratios springed/nominal: AoA = $(round(trimAoASpringed/trimAoA[i],sigdigits=4)), thrust = $(round(trimThrustSpringed/trimThrust[i],sigdigits=4)), δ = $(round(trimδSpringed/trimδ[i],sigdigits=4))")

    ## All the variables needed for the stability analysis are now in place. We create the model for the eigenproblem, using the trim variables found previously with the springed model in order to solve for the stability around that exact state. 
    BWBeigen = first(create_BWB(aeroSolver=aeroSolver,altitude=h,airspeed=U,δElev=trimδSpringed,thrust=trimThrustSpringed))

    ## Now we create and solve the eigenproblem. Notice that by using `solve_eigen!()`, we skip the step of finding the steady state of the problem, making use of the known trim solution (with the keyword argument `refTrimProblem`). We apply a filter to find only modes whose frequencies are greater than 1 rad/s through the keyword argument `frequencyFilterLimits`
    eigenProblem[i] = create_EigenProblem(model=BWBeigen,nModes=nModes,frequencyFilterLimits=[1,Inf],refTrimProblem=trimProblemSpringed[i])
    solve_eigen!(eigenProblem[i])

    ## The final step in the loop is extracting the frequencies, dampings and eigenvectors of the solution
    untrackedFreqs[i] = eigenProblem[i].frequenciesOscillatory
    untrackedDamps[i] = eigenProblem[i].dampingsOscillatory
    untrackedEigenvectors[i] = eigenProblem[i].eigenvectorsOscillatoryCplx
end
#md nothing #hide

# ### Post-processing
# We can use the AeroBeams built-in function [`mode_tracking_hungarian`](@ref mode_tracking_hungarian) to enhance the chances of correctly tracking the frequencies and dampings of each mode
## Mode tracking
freqs,damps,_,matchedModes = mode_tracking_hungarian(URange,untrackedFreqs,untrackedDamps,untrackedEigenvectors)

## Separate frequencies and damping ratios by mode
modeDampings = Array{Vector{Float64}}(undef,nModes)
modeFrequencies =  Array{Vector{Float64}}(undef,nModes)
for mode in 1:nModes
    modeDampings[mode] = [damps[i][mode] for i in eachindex(URange)]
    modeFrequencies[mode] = [freqs[i][mode] for i in eachindex(URange)]
end
#md nothing #hide

# The flutter onset speeds and frequencies are computed for all modes
## Flutter onset speed of each mode
flutterOnsetSpeedOfMode = fill(NaN, nModes)
flutterOnsetFreqOfMode = fill(NaN, nModes)
for mode in 1:nModes
    iOnset = findfirst(j -> modeDampings[mode][j] < 0 && modeDampings[mode][j+1] > 0, 1:length(URange)-1)
    flutterOnsetSpeedOfMode[mode] = isnothing(iOnset) ? Inf : interpolate(modeDampings[mode][iOnset:iOnset+1],URange[iOnset:iOnset+1],0)
    flutterOnsetFreqOfMode[mode] = isnothing(iOnset) ? Inf : interpolate(modeDampings[mode][iOnset:iOnset+1],modeFrequencies[mode][iOnset:iOnset+1],0)
end
flutterModes = findall(!isinf, flutterOnsetSpeedOfMode)
flutterOnsetSpeedsAll = filter(!isinf,flutterOnsetSpeedOfMode)
flutterOnsetFreqsAll = filter(!isinf,flutterOnsetFreqOfMode)

flutterMode = flutterModes[argmin(flutterOnsetSpeedsAll)]
flutterOnsetSpeed = minimum(flutterOnsetSpeedsAll)
flutterSpeedInd = findfirst(x-> x>flutterOnsetSpeed, URange)
for (mode,speed,freq) in zip(flutterModes,flutterOnsetSpeedsAll,flutterOnsetFreqsAll)
    println("Mode $mode: flutter speed = $(round(speed,sigdigits=4)) m/s, flutter frequency = $(round(freq,digits=2)) rad/s")
end

# We can load the reference solution found with the University of Michigan's Nonlinear Aeroelastic Simulation Tool (UM/NAST) in its 2024 version (not the version in Su's thesis).
## Load reference data
trimAoARef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/BWB/trimAoA.txt")
trimThrustRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/BWB/trimThrust.txt")
trimδRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/BWB/trimDelta.txt")
freqsRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/BWB/freqs.txt")
dampsRef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/BWB/damps.txt")
#md nothing #hide

#md # We are ready to plot the results. The following plots show the trim root angle of attack, motor thrust and elevator deflection as functions of the airspeed. The correlation with the reference solution is very good, except for the thrust. This difference arises because of the way in which the airfoil tangential force is computed when a flap deflection is present: either as proportional to the total normal force, which includes the flap-induced component (in AeroBeams), or proportional to the effective angle of attack (in UM/NAST).
#md using Suppressor #hide
#md using Plots, ColorSchemes
#md gr()
#md ENV["GKSwstype"] = "100" #hide
#md nothing #hide

## Root AoA
#md plt1 = plot(xlabel="Airspeed [m/s]", ylabel="Trim root AoA [deg]", xlims=[0,160], xticks=vcat(30:10:160))
#md plot!(URange, trimAoA*180/π, c=:black, lw=2, label="AeroBeams")
#md scatter!(trimAoARef[1,:],trimAoARef[2,:], c=:black, ms=4, label="UM/NAST")
#md savefig("BWBflutter_AoA.svg") #hide

## Thrust
#md plt2 = plot(xlabel="Airspeed [m/s]", ylabel="Trim thrust [N]", xlims=[30,160], xticks=vcat(0:10:160), legend=:bottomright)
#md plot!(URange, trimThrust, c=:black, lw=2, label="AeroBeams")
#md scatter!(trimThrustRef[1,:],trimThrustRef[2,:], c=:black, ms=4, label="UM/NAST")
#md savefig("BWBflutter_thrust.svg") #hide

## Elevator deflection
#md plt3 = plot(xlabel="Airspeed [m/s]", ylabel="Trim elevator deflection [deg]", xlims=[30,160], xticks=vcat(0:10:160), legend=:bottomright)
#md plot!(URange, trimδ*180/π, c=:black, lw=2, label="AeroBeams")
#md scatter!(trimδRef[1,:],trimδRef[2,:], c=:black, ms=4, label="UM/NAST")
#md savefig("BWBflutter_delta.svg") #hide
#md nothing #hide

#md # ![](BWBflutter_AoA.svg)
#md # ![](BWBflutter_thrust.svg)
#md # ![](BWBflutter_delta.svg)

#md # The stability results can be visualized through the following root locus and V-g-f (frequency and damping evolution) plots. The lowest flutter speed is computed by AeroBeams at approximately 120.7 m/s, at a frequency of 47.4 rad/s. Conversely, the flutter speed and frequency predicted by UM/NAST are around 137.7 m/s and 47.5 rad/s, yielding differences of 12.3% and 0.2%, on speed and frequency, respectively. Despite the flutter speed disagreement, there is qualitatively good agreement among the models for the behavior of the roots.
## Colormap
#md modeColors = cgrad(:rainbow, nModes, categorical=true)

## Root locus
#md plt4 = plot(xlabel="Damping [1/s]", ylabel="Frequency [rad/s]", xlims=[-30,5],ylims=[0,125])
#md scatter!([NaN],[NaN], c=:black, shape=:circle, ms=4, msw=0, label="AeroBeams")
#md scatter!([NaN],[NaN], c=:black, shape=:utriangle, ms=4, msw=0, label="UM/NAST")
#md for mode in 1:nModes
#md scatter!(dampsRef[mode+1,:], 2π*freqsRef[mode+1,:], c=:black, shape=:utriangle, ms=4, msw=0, label=false)
#md scatter!(modeDampings[mode], modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=4, msw=0, label=false)
#md end
#md savefig("BWBflutter_rootlocus.svg") #hide

## V-g-f
#md plt51 = plot(ylabel="Frequency [rad/s]", xlims=[30,160], xticks=0:10:160, ylims=[0,125])
#md for mode in 1:nModes
#md plot!(URange, modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=4, msw=0, label=false)
#md end
#md plt52 = plot(xlabel="Airspeed [m/s]", ylabel="Damping ratio", xlims=[30,160], ylims=[-0.3,0.2], xticks=30:10:160, yticks=-0.3:0.1:0.2)
#md for mode in 1:nModes
#md plot!(URange, modeDampings[mode]./modeFrequencies[mode], c=modeColors[mode], shape=:circle, ms=4, msw=0, label=false)
#md end
#md plt5 = plot(plt51,plt52, layout=(2,1))
#md savefig("BWBflutter_Vgf.svg") #hide
#md nothing #hide

#md # ![](BWBflutter_rootlocus.svg)
#md # ![](BWBflutter_Vgf.svg)

#md # Finally, we may visualize the mode shapes of the last eigenproblem (at highest airspeed), making use of the [`plot_mode_shapes`](@ref plot_mode_shapes) function with the appropriate inputs. Modes 1 and 2 seem to respectively be lateral-directional and longitudinal flight dynamic modes, whereas the others are structural.
## Plot mode shapes  
#md @suppress begin #hide
#md modesPlot = plot_mode_shapes(eigenProblem[end],scale=2,view=(45,30),element2centralize=11,legendPos=:top,modalColorScheme=:rainbow)
#md end #hide
#md savefig("BWBflutter_modeShapes.svg") #hide
#md nothing #hide
#md # ![](BWBflutter_modeShapes.svg)

println("Finished BWBflutter.jl") #src