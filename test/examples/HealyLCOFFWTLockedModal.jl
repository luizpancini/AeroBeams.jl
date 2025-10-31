using AeroBeams, PrettyTables

# Flags
wingtipTab = true
wingtipLocked = true
verticalMount = true

# Discretization
nElementsInner = 16
nElementsOuter = 4
nElementsFFWT = 8

# Number of modes
nModes = 6

# Model
model = create_HealyLCOFFWT(verticalMount=verticalMount,wingtipLocked=wingtipLocked,wingtipTab=wingtipTab,nElementsInner=nElementsInner,nElementsOuter=nElementsOuter,nElementsFFWT=nElementsFFWT)

# Create and solve problem
problem = create_EigenProblem(model=model,nModes=nModes)
solve!(problem)

# Get frequencies in Hz
freqs = problem.frequenciesOscillatory/(2π)

# Reference data (Tables 5.3 and 5.4 of Healy's thesis)
freqsGVT = [1.62; 8.81; 9.31; 22.37; 30.53; 30.74]
freqsNastran = freqsGVT .* [1.03; 0.9; 0.8; 0.92; 0.97; 0.63]

# Mode labels
modeLabels = ["OOP1","OOP2","IP1","OOP3","IP2","T"]

# Relative errors
ϵ_rel_GVT = freqs ./ freqsGVT .- 1
ϵ_rel_Nastran = freqs ./ freqsNastran .- 1

# Build table
data = hcat(modeLabels, freqs, freqsGVT, freqsNastran, ϵ_rel_GVT, ϵ_rel_Nastran)
pretty_table(
    data,
    header = ["Mode shape", "AeroBeams [Hz]", "GVT [Hz]", "Nastran [Hz]", "Rel. Err. GVT", "Rel. Err. Nastran"],
    formatters = ft_round(2), # round to 2 decimals
    alignment = :c
)

println("Finished HealyLCOFFWTLockedModal.jl")