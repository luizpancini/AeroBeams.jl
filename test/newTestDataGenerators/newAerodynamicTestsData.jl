# New reference data for aerodynamic tests

# Dynamic analysis of an airfoils with harmonic flap deflection profile
include("../examples/flapOscillation.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flapOscillation"))
writedlm("test/newTestDataGenerators/flapOscillation/cn.txt", cn)
writedlm("test/newTestDataGenerators/flapOscillation/cm.txt", cm)

# Dynamic analysis of two airfoils with linked harmonic flap deflection profiles
include("../examples/flapOscillationLinked.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/flapOscillationLinked"))
writedlm("test/newTestDataGenerators/flapOscillationLinked/cnMaster.txt", cnMaster)
writedlm("test/newTestDataGenerators/flapOscillationLinked/cmMaster.txt", cmMaster)
writedlm("test/newTestDataGenerators/flapOscillationLinked/cnSlave.txt", cnSlave)
writedlm("test/newTestDataGenerators/flapOscillationLinked/cmSlave.txt", cmSlave)

# One-minus-cosine gust response of an airfoil section at several pitch angles
include("../examples/OMCgustTests.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/OMCgustTests"))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
        for (k,testCase) in enumerate(tests)
            if typeof(aeroSolver) == QuasiSteady
                aeroSolverName = "QS"
            elseif typeof(aeroSolver) == Indicial
                aeroSolverName = "Indicial"
            elseif typeof(aeroSolver) == Inflow
                aeroSolverName = "Inflow"
            elseif typeof(aeroSolver) == BLi
                aeroSolverName = "BLi"
            elseif typeof(aeroSolver) == BLo
                aeroSolverName = "BLo"
            end
            gustSolverName = gustLoadsSolver.indicialFunctionName
            relPath = string("/test/newTestDataGenerators/OMCgustTests/OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase)
            absPath = string(pwd(),relPath)
            mkpath(absPath)
            writedlm(string("test/newTestDataGenerators/OMCgustTests/OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/tau",i,j,k,".txt"), τ[i,j,k])
            writedlm(string("test/newTestDataGenerators/OMCgustTests/OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dcl",i,j,k,".txt"), Δcl[i,j,k])
            writedlm(string("test/newTestDataGenerators/OMCgustTests/OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dclRef",i,j,k,".txt"), ΔclRef[i,j,k])
        end
    end
end


# Dynamic analysis of a harmonically pitching airfoil, using the dynamic stall model
include("../examples/pitchingAirfoilDSModelTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/pitchingAirfoilDSModelTest"))
writedlm("test/newTestDataGenerators/pitchingAirfoilDSModelTest/"*frameString*"_cn.txt", cn)
writedlm("test/newTestDataGenerators/pitchingAirfoilDSModelTest/"*frameString*"_cm.txt", cm)
writedlm("test/newTestDataGenerators/pitchingAirfoilDSModelTest/"*frameString*"_ct.txt", ct)

# Dynamic analysis of a harmonically plunging airfoil, using the dynamic stall model
include("../examples/plungingAirfoilDSModelTest.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/plungingAirfoilDSModelTest"))
writedlm("test/newTestDataGenerators/plungingAirfoilDSModelTest/cn.txt", cn)
writedlm("test/newTestDataGenerators/plungingAirfoilDSModelTest/cm.txt", cm)

# Sharp-edged gust response of an airfoil section at several pitch angles
include("../examples/SEgustTests.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/SEgustTests"))
for (i,aeroSolver) in enumerate(aeroSolvers)
    for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
        for (k,testCase) in enumerate(tests)
            if typeof(aeroSolver) == QuasiSteady
                aeroSolverName = "QS"
            elseif typeof(aeroSolver) == Indicial
                aeroSolverName = "Indicial"
            elseif typeof(aeroSolver) == Inflow
                aeroSolverName = "Inflow"
            elseif typeof(aeroSolver) == BLi
                aeroSolverName = "BLi"
            elseif typeof(aeroSolver) == BLo
                aeroSolverName = "BLo"
            end
            gustSolverName = gustLoadsSolver.indicialFunctionName
            relPath = string("/test/newTestDataGenerators/SEgustTests/SEgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase)
            absPath = string(pwd(),relPath)
            mkpath(absPath)
            writedlm(string("test/newTestDataGenerators/SEgustTests/SEgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/tau",i,j,k,".txt"), τ[i,j,k])
            writedlm(string("test/newTestDataGenerators/SEgustTests/SEgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dcl",i,j,k,".txt"), Δcl[i,j,k])
            writedlm(string("test/newTestDataGenerators/SEgustTests/SEgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dclRef",i,j,k,".txt"), ΔclRef[i,j,k])
        end
    end
end

# Dynamic analysis of an airfoil section sinusoidally surging (facing a time-varying freestream)
include("../examples/timeVaryingFreestream.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/timeVaryingFreestream"))
writedlm("test/newTestDataGenerators/timeVaryingFreestream/cn.txt", cn[end])
writedlm("test/newTestDataGenerators/timeVaryingFreestream/cm.txt", cm[end])

# Dynamic analysis of an airfoil section sinusoidally pitching and surging (facing a time-varying freestream)
include("../examples/timeVaryingFreestreamAndPitch.jl")
mkpath(string(pwd(),"/test/newTestDataGenerators/timeVaryingFreestreamAndPitch"))
writedlm("test/newTestDataGenerators/timeVaryingFreestreamAndPitch/cn.txt", cn[end])
writedlm("test/newTestDataGenerators/timeVaryingFreestreamAndPitch/cm.txt", cm[end])