# Aerodynamic problems

@testset "Dynamic analysis of a harmonically pitching airfoil, using the dynamic stall model" begin
    include("examples/DSModelTest.jl")
    # Self-comparison
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "DSModelTest", "cn.txt"))
    cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "DSModelTest", "cm.txt"))
    ct_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "DSModelTest", "ct.txt"))
    @test cn ≈ cn_ atol=SELFatol
    @test cm ≈ cm_ atol=SELFatol
    @test ct ≈ ct_ atol=SELFatol
end

@testset "Dynamic analysis of an airfoils with harmonic flap deflection profile" begin
    include("examples/flapOscillation.jl")
    # Self-comparison
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillation", "cn.txt"))
    cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillation", "cm.txt"))
    @test cn ≈ cn_ atol=SELFatol
    @test cm ≈ cm_ atol=SELFatol
end

@testset "Dynamic analysis of two airfoils with linked harmonic flap deflection profiles" begin
    include("examples/flapOscillationLinked.jl")
    # Self-comparison
    cnMaster_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillationLinked", "cnMaster.txt"))
    cmMaster_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillationLinked", "cmMaster.txt"))
    cnSlave_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillationLinked", "cnSlave.txt"))
    cmSlave_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "flapOscillationLinked", "cmSlave.txt"))
    @test cnMaster ≈ cnMaster_ atol=SELFatol
    @test cmMaster ≈ cmMaster_ atol=SELFatol
    @test cnSlave ≈ cnSlave_ atol=SELFatol
    @test cmSlave ≈ cmSlave_ atol=SELFatol
end

@testset "One-minus-cosine gust response of an airfoil section at several pitch angles" begin
    include("examples/OMCgustTests.jl")
    # Self-comparison
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
            for (k,testCase) in enumerate(1:6)
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
                Δcl_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "OMCgustTests", string("OMCgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dcl",i,j,k,".txt")))
                @test Δcl[i,j,k] ≈ Δcl_ atol=SELFatol
            end
        end
    end
end

@testset "Sharp-edged gust response of an airfoil section at several pitch angles" begin
    include("examples/SEgustTests.jl")
    # Self-comparison
    for (i,aeroSolver) in enumerate(aeroSolvers)
        for (j,gustLoadsSolver) in enumerate(gustLoadsSolvers)
            for (k,testCase) in enumerate(1:3)
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
                Δcl_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "SEgustTests", string("SEgustTests_",aeroSolverName,"_",gustSolverName,"_test",testCase,"/dcl",i,j,k,".txt")))
                @test Δcl[i,j,k] ≈ Δcl_ atol=SELFatol
            end
        end
    end
end

@testset "Dynamic analysis of a harmonically plunging airfoil, using the dynamic stall model" begin
    include("examples/plungingAirfoil.jl")
    # Self-comparison
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "plungingAirfoil", "cn.txt"))
    cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "plungingAirfoil", "cm.txt"))
    @test cn ≈ cn_ atol=SELFatol
    @test cm ≈ cm_ atol=SELFatol
end

@testset "Dynamic analysis of an airfoil section sinusoidally surging (facing a time-varying freestream)" begin
    include("examples/timeVaryingFreestream.jl")
    # Loop normalized airspeed amplitude
    for i in eachindex(λᵤRange)
        # Reference comparison
        @test cn[i][end]/(2π*θ) ≈ cnCFD[i][2,end] atol=0.3
        @test cm[i][end] ≈ cmCFD[i][2,end] atol=1e-2
        # Analytical comparison
        @test Vdot2[i][end] ≈ Vdot2Analytical[i][end] rtol=1e-2
        @test Vdot3[i][end] ≈ Vdot3Analytical[i][end] rtol=1e-2
    end
    # Self-comparison (for the last λᵤ)
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "timeVaryingFreestream", "cn.txt"))
    cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "timeVaryingFreestream", "cm.txt"))
    @test cn[end] ≈ cn_ atol=SELFatol
    @test cm[end] ≈ cm_ atol=SELFatol
end

@testset "Dynamic analysis of an airfoil section sinusoidally pitching and surging (facing a time-varying freestream)" begin
    include("examples/timeVaryingFreestreamAndPitch.jl")
    # Loop normalized airspeed amplitude
    for i in eachindex(λᵤRange)
        # index at 3/4 of the cycle
        indexTqcycle = rangeLastCycle[i][round(Int,length(rangeLastCycle[i])*3/4)]
        # Reference comparison
        @test cn[i][end]/cn_qs(t[i][end]) ≈ cnCFD[i][2,end] atol=0.5
        @test cm[i][end] ≈ cmCFD[i][2,end] atol=1e-2
        # Analytical comparison
        @test Vdot2[i][end] ≈ Vdot2Analytical[i][end] rtol=1e-2
        @test Vdot3[i][end] ≈ Vdot3Analytical[i][end] rtol=1e-2
        @test Ω1[i][end] ≈ Ω1Analytical(t[i][end]) rtol=1e-2
        @test Ωdot1[i][indexTqcycle] ≈ Ωdot1Analytical(t[i][indexTqcycle]) rtol=0.1
    end
    # Self-comparison (for the last λᵤ)
    cn_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "timeVaryingFreestreamAndPitch", "cn.txt"))
    cm_ = readdlm(joinpath(@__DIR__, "newTestDataGenerators", "timeVaryingFreestreamAndPitch", "cm.txt"))
    @test cn[end] ≈ cn_ atol=SELFatol
    @test cm[end] ≈ cm_ atol=SELFatol
end
