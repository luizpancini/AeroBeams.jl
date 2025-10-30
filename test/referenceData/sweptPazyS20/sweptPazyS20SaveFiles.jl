using AeroBeams, JLD2, MAT

# File IDs
IDs = vcat(3020,3031,3034:3038,3201:3208)

# Loop
for ID in IDs
    # Load .mat file
    matfile = matopen(pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/ID_"*string(ID)*"_flutter_data.mat")
    data = read(matfile)
    close(matfile)
    # Save as .jld2
    @save pkgdir(AeroBeams)*"/test/referenceData/sweptPazyS20/ID_"*string(ID)*"_flutter_data.jld2" data
end

println("Finished sweptPazyS20SaveFiles.jl")