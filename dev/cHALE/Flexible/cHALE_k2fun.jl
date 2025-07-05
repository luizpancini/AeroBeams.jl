using AeroBeams, LinearInterpolations

# Aerodynamic solver
aeroSolver = Indicial()

# Options for stabilizers
stabilizersAero = true
includeVS = true
wingCd0 = stabsCd0 = 1e-2

# Option to include induced drag
hasInducedDrag = true

# Stiffness factor
λ = 1

# Altitude
h = 20e3

# Discretization
nElemWing = 20
nElemTailBoom = 2
nElemHorzStabilizer = 2
nElemVertStabilizer = 1

# Linearly spanwise-varying bending curvature
k2root = 0.03
k2 = x1 -> k2root*(1-(x1/16)^(1/2))

# Wing model
wingModel,_ = create_SMW(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=1,nElem=div(nElemWing,2),altitude=h,cd0=wingCd0,k2=k2,hasInducedDrag=hasInducedDrag,θ=0)

# Plot undeformed wing model
pltWing = plot_undeformed_assembly(wingModel; view=(0,0))
display(pltWing)

# Aircraft model
cHALE,_ = create_conventional_HALE(aeroSolver=aeroSolver,stiffnessFactor=λ,airspeed=1,nElemWing=nElemWing,nElemTailBoom=nElemTailBoom,nElemHorzStabilizer=nElemHorzStabilizer,nElemVertStabilizer=nElemVertStabilizer,stabilizersAero=stabilizersAero,includeVS=includeVS,wingCd0=wingCd0,stabsCd0=stabsCd0,δElevIsTrimVariable=stabilizersAero,thrustIsTrimVariable=true,k2=k2,hasInducedDrag=hasInducedDrag,altitude=h)

# Plot undeformed wing model
pltAircraft = plot_undeformed_assembly(cHALE; view=(0,0))
display(pltAircraft)

println("Finished cHALE_k2fun.jl")