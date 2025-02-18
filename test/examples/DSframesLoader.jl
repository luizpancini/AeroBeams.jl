# Loads dynamic stall data from a frame of the NASA (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19830003778.pdf)
function NASA_frames_loader(frame::Int)

    # Atmosphere
    altitude = 0
    atmosphere = standard_atmosphere(altitude)

    # Frame data lookup (airfoil, mean angle, angle amplitude, airfoil semichord, reduced frequency, Mach number)
    frames = Dict(
        10_022 => (airfoil=deepcopy(NACA0012), a₀=12π/180, a₁=9.9π/180, b=0.305, k=0.098, Ma=0.301)
    )

    # Retrieve frame data
    frame_data = get(frames, frame, nothing)
    isnothing(frame_data) && error("Frame $frame is unavailable")

    # Extract variables
    airfoil,a₀,a₁,b,k,Ma = frame_data

    # Airspeed
    U = Ma*atmosphere.a

    # Update airfoil parameters
    update_Airfoil_params!(airfoil,Ma=Ma,U=U,b=b)

    return airfoil,a₀,a₁,b,k,Ma,U
end

# Loads dynamic stall data from a frame of the University of Glasgow (DOI:10.5525/gla.researchdata.464)
function GU_frames_loader(frame::Int)

    # Frame data lookup (airfoil, mean angle, angle amplitude, airfoil semichord, reduced frequency, Mach number, airspeed)
    frames = Dict(
        02_000_101 => (airfoil=deepcopy(NACA23012A), a₀=0.24817, a₁=0.26967, b=0.275, k=0.001, Ma=0.11931, U=41.570),
        02_010_151 => (airfoil=deepcopy(NACA23012A), a₀=0.18100, a₁=0.10401, b=0.275, k=0.12489, Ma=0.11635, U=40.605),
        02_010_351 => (airfoil=deepcopy(NACA23012A), a₀=0.17899, a₁=0.17741, b=0.275, k=0.174, Ma=0.11694, U=40.815)
    )

    # Retrieve frame data
    frame_data = get(frames, frame, nothing)
    isnothing(frame_data) && error("Frame $frame is unavailable")

    # Extract variables
    airfoil,a₀,a₁,b,k,Ma,U = frame_data

    # Update airfoil parameters
    update_Airfoil_params!(airfoil,Ma=Ma,U=U,b=b)

    return airfoil,a₀,a₁,b,k,Ma,U
end