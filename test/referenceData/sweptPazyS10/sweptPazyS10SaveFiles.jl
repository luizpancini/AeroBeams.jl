using AeroBeams, JLD2, MAT

# File IDs
IDs = [4023; 4024; 4032; 4033; 4034; 4035; 4303; 4306; 4307; 5208; 5209]

# Loop
for ID in IDs
    # Load .mat file
    matfile = matopen(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/ID_"*string(ID)*"_flutter_data.mat")
    data = read(matfile)
    close(matfile)
    # Save as .jld2
    @save pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS10/ID_"*string(ID)*"_flutter_data.jld2" data
end

println("Finished sweptPazyS10SaveFiles.jl")