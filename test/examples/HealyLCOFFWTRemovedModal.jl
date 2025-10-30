using AeroBeams, PrettyTables

# Flags
wingtipRemoved = true
verticalMount = true

# Discretization
nElementsInner = 16
nElementsOuter = 4
nElementsFFWT = 8

# Number of modes
nModes = 9

# Model
model = create_HealyLCOFFWT(verticalMount=verticalMount,wingtipRemoved=wingtipRemoved,nElementsInner=nElementsInner,nElementsOuter=nElementsOuter,nElementsFFWT=nElementsFFWT)

# Create and solve problem
problem = create_EigenProblem(model=model,nModes=nModes)
solve!(problem)

# Get frequencies in Hz
freqs = problem.frequenciesOscillatory/(2π)

# Reference data (Fig. 5.13 of Healy's thesis)
freqsNastran = [2.6; 15.8; 15.9; 36.2; 57.2; 64.8; 99.4; 167.8; 112.7]

# Mode labels
modeLabels = ["OOP1","IP1","OOP2","OOP3","T1","OOP4","IP2","T2","OOP5"]

# Relative errors
ϵ_rel_Nastran = freqs ./ freqsNastran .- 1

# Build table
data = hcat(modeLabels, freqs, freqsNastran, ϵ_rel_Nastran)
pretty_table(
    data,
    header = ["Mode shape", "AeroBeams [Hz]", "Nastran [Hz]", "Rel. Err. Nastran"],
    formatters = ft_round(2), # round to 2 decimals
    alignment = :c
)

println("Finished HealyLCOFFWTRemovedModal.jl")