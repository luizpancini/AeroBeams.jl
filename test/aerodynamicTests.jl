@testset "Dynamic analysis of an airfoils with harmonic flap deflection profile" begin
    include_example("examples/aerodynamic/flapOscillation.jl")
    @test matches_test_data(cn, "flapOscillation", "cn.txt")
    @test matches_test_data(cm, "flapOscillation", "cm.txt")
end

@testset "Dynamic analysis of two airfoils with linked harmonic flap deflection profiles" begin
    include_example("examples/aerodynamic/flapOscillationLinked.jl")
    @test matches_test_data(cnMaster, "flapOscillationLinked", "cnMaster.txt")
    @test matches_test_data(cmMaster, "flapOscillationLinked", "cmMaster.txt")
    @test matches_test_data(cnSlave, "flapOscillationLinked", "cnSlave.txt")
    @test matches_test_data(cmSlave, "flapOscillationLinked", "cmSlave.txt")
end

!is_ci && @testset "One-minus-cosine gust response of an airfoil section at several pitch angles" begin
    include_example("examples/aerodynamic/OMCgustTests.jl")
    @test matches_test_data(Δcl, "OMCgustTests", "dcl.txt")
end

@testset "One-minus-cosine gust response of an airfoil section at several pitch angles" begin
    include_example("examples/aerodynamic/OMCgustTests2.jl")
    @test matches_test_data(Δcl, "OMCgustTests2", "dcl.txt")
end

@testset "Dynamic analysis of a harmonically pitching airfoil, using the dynamic stall model" begin
    include_example("examples/aerodynamic/pitchingAirfoilDSModelTest.jl")
    @test matches_test_data(cn, "pitchingAirfoilDSModelTest", frameString*"_cn.txt")
    @test matches_test_data(cm, "pitchingAirfoilDSModelTest", frameString*"_cm.txt")
    @test matches_test_data(ct, "pitchingAirfoilDSModelTest", frameString*"_ct.txt")
end

@testset "Dynamic analysis of a harmonically plunging airfoil, using the dynamic stall model" begin
    include_example("examples/aerodynamic/plungingAirfoilDSModelTest.jl")
    @test matches_test_data(cn, "plungingAirfoilDSModelTest", "cn.txt")
    @test matches_test_data(cm, "plungingAirfoilDSModelTest", "cm.txt")
end

!is_ci && @testset "Sharp-edged gust response of an airfoil section at several pitch angles" begin
    include_example("examples/aerodynamic/SEgustTests.jl")
    @test matches_test_data(Δcl, "SEgustTests", "dcl.txt")
end

@testset "Dynamic analysis of an airfoil section sinusoidally surging (facing a time-varying freestream)" begin
    include_example("examples/aerodynamic/timeVaryingFreestream.jl")
    @test matches_test_data(cn[end], "timeVaryingFreestream", "cn.txt")
    @test matches_test_data(cm[end], "timeVaryingFreestream", "cm.txt")
    # The relative wind accelerations match the analytical ones from the first time step on, without oscillations (the airspeed does not need to start smoothly)
    for i in eachindex(λᵤRange)
        @test maximum(abs, Vdot2[i] - Vdot2Analytical[i]) < 1e-3 * maximum(abs, Vdot2Analytical[i])
        @test maximum(abs, Vdot3[i] - Vdot3Analytical[i]) < 1e-3 * maximum(abs, Vdot3Analytical[i])
    end
end

@testset "Dynamic analysis of an airfoil section sinusoidally pitching and surging (facing a time-varying freestream)" begin
    include_example("examples/aerodynamic/timeVaryingFreestreamAndPitch.jl")
    @test matches_test_data(cn[end], "timeVaryingFreestreamAndPitch", "cn.txt")
    @test matches_test_data(cm[end], "timeVaryingFreestreamAndPitch", "cm.txt")
end