## CONSTANTS
# ------------------------------------------------------------------------------
# Cache for swept Pazy wing VLM tip loss data
const _VLM_cache = let

    # Common sweep angle range
    ΛRange = deg2rad.([0,10,20,30])

    # -------- VLM-undef --------
    polyCoeffsUndef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/polyCoeffsUndef.txt")
    itpsUndef = [LinearInterpolations.Interpolate(ΛRange, polyCoeffsUndef[i,:]) for i in 1:9]

    # -------- VLM-def ---------
    θRange = deg2rad.([0,3,5,7])
    URange = vcat(0:10:70, 120)
    polyCoeffsDef = zeros(9, length(ΛRange), length(θRange), length(URange))
    for j in eachindex(θRange), k in eachindex(URange)
        polyCoeffsDef[:,:,j,k] = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/polyCoeffsDef_$(j)_$(k).txt")
    end
    itpsDef = [LinearInterpolations.Interpolate((ΛRange, θRange, URange), polyCoeffsDef[i,:,:,:]) for i in 1:9]

    (ΛRange = ΛRange,
    itpsUndef = itpsUndef,
    θRange = θRange,
    URange = URange,
    itpsDef = itpsDef,)

end

## PRIVATE API
# ------------------------------------------------------------------------------
# Chordwise position of the reference axis (spar) of the straight Pazy wing in the original AePW3 data (https://github.com/UM-A2SRL/AePW3-LDWG), a fraction of the chord
const normSparPos_Pazy_AePW3 = 0.44096

# Chordwise position of the reference axis (spar) of the Pazy wing, a fraction of the chord: the shear center, which depends on the sweep angle (the S20 value is provisional)
function shear_center_Pazy(Λ)

    if isapprox(rad2deg(Λ), 0; atol=1e-8)
        return 0.4653
    elseif isapprox(rad2deg(Λ), 10; atol=1e-8)
        return 0.4470
    else
        return 0.4400
    end

end

# Checks the mesh of the Pazy wing: :nominal (the standard beam mesh, with the nodes at the root, at the ribs and at the tip of the spar) or :refined (with a node at the middle of each element of the :nominal mesh, which has approximately twice as many elements)
function validate_mesh_Pazy(mesh)

    @assert mesh in (:nominal, :refined) "mesh must be :nominal or :refined"

    return mesh

end

# Chord of the Pazy wing [m], measured along the chordwise (x) direction of the MSC Nastran frame
chord_Pazy(Λ::Real) = isapprox(rad2deg(Λ), 0; atol=1e-8) ? 0.0988504 : 0.0950

# Path of the CSV file with the reference data of the Pazy wing with sweep angle Λ (0, 10 or 20 degrees) for the given mesh (:nominal or :refined), in test/referenceData/Pazy, or in test/referenceData/PazyNoRod for the wing without the tip rod (the straight wing, and only the inertia of the swept wing). `kind` is :inertia (S<Λ>_inertia_matrices_<mesh>mesh.csv for the straight wing, and S<Λ>_<config>_inertia_matrices_<mesh>mesh.csv for the swept one, with the tip mass configuration `config`, which the straight wing does not have, so it is ignored; S<Λ>_norod_… without the tip rod) or :stiffness, S<Λ>_stiffness_matrices_<mesh>mesh_<tip>.csv (S<Λ>_norod_stiffness_matrices_<mesh>mesh.csv without the tip rod, which has no `tip`), the least squares ones, selected by the Nastran model with the rigid (`tip = :rigidtip`, the tip section of the wing is rigid, as the stiffness matrices require) or flexible (`:flextip`) tip rod. With `rootTorsionCorr = true` the stiffness data is that with the correction of the root torsional stiffness (and bending-torsion coupling) of the straight wing with the tip rod, S0_stiffness_matrices_<mesh>mesh_<tip>_corr.csv. Throws a descriptive error if the data is not available
function Pazy_data_file(kind::Symbol, Λ, mesh::Symbol; config::Symbol=:baseline, withWingtipRod::Bool=true, tip::Symbol=:rigidtip, rootTorsionCorr::Bool=false)

    @assert kind in (:inertia, :stiffness) "kind must be :inertia or :stiffness"
    validate_mesh_Pazy(mesh)
    ΛDeg = round(Int, rad2deg(Λ))
    isStraight = ΛDeg == 0
    @assert withWingtipRod || isStraight || kind == :inertia "the stiffness matrices of the wing without the tip rod are only available for the straight wing"
    @assert config in (:baseline, :LE, :TE) "config must be :baseline, :LE or :TE"
    @assert tip in (:rigidtip, :flextip) "tip must be :rigidtip or :flextip"
    @assert !(rootTorsionCorr && kind == :stiffness) || (isStraight && withWingtipRod) "the root torsion correction of the stiffness matrices is only available for the straight wing with the tip rod"
    prefix = "S$(ΛDeg)" * (withWingtipRod ? "" : "_norod")
    # Name of the file without the extension, with the mesh in it: S<Λ>_inertia_matrices_<mesh>mesh, or S<Λ>_stiffness_matrices_<mesh>mesh_<tip> (with the suffix _corr for the root torsion correction)
    stem = if kind == :inertia
        prefix * (isStraight ? "" : "_" * String(config)) * "_inertia_matrices_" * String(mesh) * "mesh"
    else
        prefix * "_stiffness_matrices_" * String(mesh) * "mesh" * (withWingtipRod ? "_" * String(tip) : "") * (rootTorsionCorr ? "_corr" : "")
    end
    dir = joinpath(pkgdir(AeroBeams), "test", "referenceData", withWingtipRod ? "Pazy" : "PazyNoRod")
    file = joinpath(dir, stem * ".csv")

    if !isfile(file)
        # List the meshes with data (the same name with another mesh)
        pattern = Regex("^" * replace(stem, String(mesh) * "mesh" => "(nominal|refined)mesh") * "\\.csv\$")
        matches = filter(!isnothing, match.(pattern, readdir(dir)))
        available = sort([Symbol(m[1]) for m in matches])
        availableStr = isempty(available) ? "none" : join(["mesh=:$m" for m in available], "; ")
        throw(ArgumentError("no $(kind == :inertia ? "inertia" : "stiffness") matrices for the Pazy wing with sweep angle $(ΛDeg)° and mesh=:$mesh (expected file $file). Available data: $availableStr. Generate the data with the BeamROM scripts and add it to $dir"))
    end

    return file

end

# Returns the nodal coordinates and chord-normalized position of the reference axis (shear center of the spar) of the Pazy wing and swept variations, in the MSC Nastran default frame (x chordwise from the leading edge, y spanwise unswept, y=0 at the root of the spar).
function reference_axis_Pazy(Λ; mesh)

    validate_mesh_Pazy(mesh)

    # Reference axis x, y and z coordinates
    file = Pazy_data_file(:inertia, Λ, mesh)
    x = Float64.(readdlm(file, ','; header=true)[1][:, 2])
    y = Float64.(readdlm(file, ','; header=true)[1][:, 3])
    z = Float64.(readdlm(file, ','; header=true)[1][:, 4])

    # Reference axis position (shear center) for the sweep angle
    ra = shear_center_Pazy(Λ)

    return x,y,z,ra

end

# Sectional inertia matrix for a beam element given its mass, center of mass offset, inertia tensor, and length
sectional_inertia_elem(m, η, J, Δℓ) = 1/Δℓ * [m*I3          -m*tilde(η); 
                                              m*tilde(η) J-m*tilde(η)^2]

# Transforms the nodal mass data of a swept wing from the MSC Nastran basic frame (x chordwise positive aft, y spanwise along the unswept direction, z vertical) to the local beam frame (x1 along the swept reference axis, x2 toward the leading edge, x3 vertical), and sets the elemental sectional inertia matrices
function inertia_matrices_from_Nastran_data(; Λ,nElem,L,nodalPositions,m,m_x,m_y,m_z,Ixx,Iyy,Izz,Ixy,Ixz,Iyz)

    # Rotation tensor from the Nastran basis to the undeformed beam basis b
    R = [ sin(Λ) cos(Λ) 0; 
         -cos(Λ) sin(Λ) 0; 
              0      0  1]

    # Center of mass offsets and inertia tensors about the center of mass, resolved in the local beam frame
    η_cm = [R*[m_x[n]; m_y[n]; m_z[n]] for n in eachindex(m)]
    J = [R*[Ixx[n] Ixy[n] Ixz[n]; 
            Ixy[n] Iyy[n] Iyz[n]; 
            Ixz[n] Iyz[n] Izz[n]]*R' for n in eachindex(m)]

    # Set sectional inertia matrices - distribute nodal masses and inertias to elements
    I = [zeros(6, 6) for _ in 1:nElem]
    for k in 1:nElem
        Δℓ = L*(nodalPositions[k+1]-nodalPositions[k])
        wStart = k == 1 ? 1.0 : 0.5
        wEnd = k == nElem ? 1.0 : 0.5
        I[k] += sectional_inertia_elem(wStart * m[k], [-Δℓ / 2, 0, 0] + η_cm[k], wStart * J[k], Δℓ)
        I[k] += sectional_inertia_elem(wEnd * m[k+1], [Δℓ / 2, 0, 0] + η_cm[k+1], wEnd * J[k+1], Δℓ)
    end

    return I

end                                              

# Sectional inertia matrices of the Pazy wing with sweep angle Λ, read from `file` (default: the inertia matrices of the given mesh and tip mass configuration `config` (only for the swept wing), see `Pazy_data_file`). The nodal data is in the MSC Nastran basic frame. `applyTipInertiaCorrection` applies an ad hoc factor on the tip node's inertia about the Nastran y axis for the :LE and :TE configurations to approximate GVT torsional mode data. This correction accounts for the flexibility of the wingtip rod and the breaking beam of the flutter protection system.
function inertia_matrices_Pazy(Λ; mesh::Symbol=:nominal, config::Symbol=:baseline, applyTipInertiaCorrection::Bool=false, withWingtipRod::Bool=true, file::String=Pazy_data_file(:inertia, Λ, mesh; config=config, withWingtipRod=withWingtipRod))

    # Geometrical properties
    nElem,L,_,_,nodalPositions = geometrical_properties_Pazy(mesh=mesh, Λ=Λ)

    # Nodal data of the equivalent beam mass model (Nastran frame): mass, center of mass offsets from the reference axis and inertia tensor about the center of mass
    data, header = readdlm(file, ','; header=true)
    column(name) = Float64.(data[:, findfirst(==(name), vec(header))])
    @assert size(data,1) == nElem+1 "the beam mass model data of the Pazy wing with sweep angle $(round(Int, rad2deg(Λ)))° for mesh=:$mesh has $(size(data,1)) nodes, but $(nElem+1) are expected"

    # Ad hoc factor on the tip node's torsional inertia for the LE and TE configurations (to approximate wingtip rod flexibility)
    Iyy = column("Jyy")
    if applyTipInertiaCorrection && config != :baseline && !(Λ ≈ 0)
        Iyy[end] *= 1.4
    end

    return inertia_matrices_from_Nastran_data(Λ=Λ,nElem=nElem,L=L,nodalPositions=nodalPositions,m=column("m"),m_x=column("dx"),m_y=column("dy"),m_z=column("dz"),Ixx=column("Jxx"),Iyy=Iyy,Izz=column("Jzz"),Ixy=column("Jxy"),Ixz=column("Jxz"),Iyz=column("Jyz"))

end

# Sectional stiffness matrices of the Pazy wing with sweep angle Λ, read from the least squares data of the nodal segments (see `Pazy_data_file`), with the rigid or flexible tip rod of the Nastran model `tip` (`:rigidtip` or `:flextip`). The matrices of the segments are assigned to the elements they contain. The coordinates of the ends of the segments are along the reference axis, from the root. With `rootTorsionCorr = true` (straight wing with the tip rod only) the matrices have the correction of the torsional stiffness (S44, and the S45 = S54 coupling for the nominal mesh) of the root, which compensates the extra stiffness at the clamp of the Nastran model and the flexibility of the real one
function stiffness_matrices_Pazy(Λ; mesh::Symbol=:nominal, withWingtipRod::Bool=true, tip::Symbol=:rigidtip, rootTorsionCorr::Bool=false)

    # Geometrical properties
    nElem,L,_,_,nodalPositions = geometrical_properties_Pazy(mesh=mesh, Λ=Λ)

    # Stiffness matrix of each segment (or element), given by its components row by row (S11, S12, ..., S66) and the coordinates of its ends from the root
    file = Pazy_data_file(:stiffness, Λ, mesh; withWingtipRod=withWingtipRod, tip=tip, rootTorsionCorr=rootTorsionCorr)
    data, header = readdlm(file, ','; header=true)
    columns = [findfirst(==("S$i$j"), vec(header)) for i in 1:6, j in 1:6]
    y1, y2 = Float64.(data[:, 2]), Float64.(data[:, 3])

    # The segments must span the beam from the root to the tip, and each element belongs to the segment that contains its middle
    desc = "the stiffness matrices data of the Pazy wing with sweep angle $(round(Int, rad2deg(Λ)))° for mesh=:$mesh ($file)"
    @assert isapprox(y1[1], 0; atol=1e-6*L) && isapprox(y2[end], L; atol=1e-6*L) "$desc do not span the beam, from 0 to $L m"
    yMid = L*(nodalPositions[1:end-1] + nodalPositions[2:end])/2
    segmentOfElement = [findfirst(i -> y1[i] <= y && y <= y2[i], eachindex(y1)) for y in yMid]
    @assert !any(isnothing, segmentOfElement) "$desc do not match the nodes of the beam"
    @assert size(data,1) == nElem "$desc has $(size(data,1)) rows, but $nElem are expected"

    return [Float64.(data[segmentOfElement[k], columns]) for k in 1:nElem]

end

# Data of the wingtip rod of the Pazy wing (Technion Nastran models nastran/Pazy_Technion/S0/S0.bdf and S10/<config>/S10_<config>.bdf), to model it as a beam. The rod is a circular bar along the chordwise direction. Positions are along the chordwise axis of the Nastran frame (x positive aft), relative to the tip node of the beam (reference axis). The slight tilt of the rod in the Nastran models (about 1e-3 rad) is neglected. The fields are:
# - the rod section and material (`d`, `A`, `I`, `J`, `E`, `G`, `ρ`, `Ks`) and its ends `xStart` and `xEnd`
# - `offset`: the position of the attachment point of the rod, relative to the same node, in the Nastran frame (with origin at the wing root shear center)
# - `tied`: ends of the parts of the rod that are modeled as rigid, due to being glued to the wing along the chord, for the bending in the plane of the wing (`inPlane`, displacements along the span) and out of it (`outOfPlane`, vertical displacements); the attachment point is inside both
# - `calibrationMass`, `xCalibration`: calibration masses and their position along the rod
# - `hinge`: for the S10 LE and TE configurations, the breaking beam (CBEAM property 24 in the Nastran models, a cone of PLA at the end of the rod, is modeled as a hinge with springs), the fixed adapter (CONM2 205, a ring at the end of the breaking beam) and the FPS weights (CONM2 433, with the offset `arm` from the end of the breaking beam, and the transverse inertia `fpsInertia` about its center of mass)
function wingtip_rod_Pazy(Λ,config)

    if isapprox(rad2deg(Λ), 0; atol=1e-8)
        # S0: the rod goes from x = -0.0841 to 0.1841 m, at y = 0.5569833 to 0.5566833 m and z = -3.167e-4 to -1.667e-5 m, and the tip node of the spar is at x = 0.045995091120000005, y = 0.5519937 m, z = 9.2695e-8 m. The grids of the rod from x = 0.00033333 to 0.0996667 m (1202 to 1222) are tied by RBE2 to the shell grids of the tip section (the chord is from x = 0 to 0.0988502 m), which is not rigid: it has the spar, the leading and trailing edge stringers (ending at the tip section, which has no rib) and a very thin skin. The stringers run along the span, so they restrain the rod much more against bending in the plane of the wing (axially) than out of it (by their bending), and the rigid part of the rod, which represents this restraint, is longer in the plane: its limits were calibrated with the natural frequencies and mode shapes of the Nastran model with the flexible rod (S0_modal_flextip.op2), as the effective clamps of the cantilevers of the rod beyond the chord (in the plane, x = 7.0 to 85.4 mm, as given by the frequencies of the modes of the rod alone, and out of it, x = 12 to 73 mm, within the spar)
        ℓ = 0.2682                  # length of the rod
        c = chord_Pazy(Λ)           # chord of the wing
        sc = shear_center_Pazy(Λ)   # shear center of the center wing stations
        xSpar = sc*c                # x-coordinate of the tip node of the spar at the reference axis (shear center)
        d = 0.01                    # diameter of the rod
        s = (0.0841 + 0.05)/ℓ       # normalized distance from the LE to the attachment point (coinciding with the midpoint) of the rod
        return (d = d, A = π*d^2/4, I = π*d^4/64, J = π*d^4/32,
                E = 1.1e9, G = 3.94548e8, ρ = 930.0, Ks = 0.886177,
                xStart = -xSpar - 0.0841, xEnd = -xSpar + 0.1841,
                offset = [0.0, 
                          0.5569833 + s*(0.5566833 - 0.5569833) - 0.5519937, 
                          -3.167e-4 + s*(-1.667e-5 + 3.167e-4) - 9.2695e-8],
                tied = (inPlane = (0.0070 - xSpar, 0.0854 - xSpar),
                        outOfPlane = (0.012 - xSpar, 0.073 - xSpar)),
                calibrationMass = 0.015 + 0.0011, xCalibration = -xSpar + 0.05,
                hinge = nothing)
    else
        # S10: the rod goes from x = 0.05289 to 0.23278 m, at y = 0.5516 to 0.55139 m, and the tip node is at x = 0.138871796846351, y = 0.54673 m. The rod crosses the plane of the tip node (the attachment point) at y = 0.5515 m, and it is made of nylon. The grids of the rod from x = 0.09734 to 0.19305 m are tied by RBE2 to the shell grids of the tip section (the chord is from x = 0.0964 to 0.1914 m), which is not rigid: in the modes of the Nastran model with the flexible rod (S10_<config>_modal_flextip.op2) the rod bends from about x = 0.17 m, so the rigid part is from 10 mm from the leading edge (0.1064 m) to the trailing edge of the spar (0.17337 m), as in the S0, which are the limits that best reproduce the natural frequencies of the Nastran models
        xSpar = 0.138871796846351
        d = 0.0095
        yAttachment = 0.5516 + (xSpar - 0.05289)/(0.23278 - 0.05289)*(0.55139 - 0.5516)
        hinge = if config == :LE
            (side = :LE, bbLength = 0.05289 - 0.052, bbDiameter = 0.0017, arm = 0.0213, adapterMass = 0.0035, fpsMass = 0.0087, fpsInertia = 3.3e-7)
        elseif config == :TE
            (side = :TE, bbLength = 0.233 - 0.23278, bbDiameter = 0.0015, arm = 0.0223, adapterMass = 0.0035, fpsMass = 0.0087, fpsInertia = 3.3e-7)
        else
            nothing
        end
        return (d = d, A = π*d^2/4, I = π*d^4/64, J = π*d^4/32,
                E = 1.4e9, G = 5.0215e8, ρ = 930.0, Ks = 0.886179,
                xStart = 0.05289 - xSpar, xEnd = 0.23278 - xSpar,
                offset = [0.0, yAttachment - 0.54673, 0.0],
                tied = (inPlane = (0.1064 - xSpar, 0.17337 - xSpar), outOfPlane = (0.1064 - xSpar, 0.17337 - xSpar)),
                calibrationMass = 0.0, xCalibration = 0.0,
                hinge = hinge)
    end

end

# Creates the beams that model the wingtip rod of the Pazy wing attached to its tip, and their connections: a rigid link from the tip node of the wing `mainWing` (with `nElem` elements) to the attachment point of the rod, and the rod, connected at the attachment point to the end of the link. The rod is uniform, with `nElementsRod` elements: the part glued to the wing is rigid, with one element between each pair of nodes at the ends of its rigid parts (`tied` in `wingtip_rod_Pazy`, one for each plane of bending, so that in the part where only one of them is rigid the other bending is flexible) and at the attachment point, and the others are distributed beyond it in proportion to their lengths. In the S10 LE and TE configurations the rod is extended beyond its end with one element, up to the center of mass of the FPS weights, which is connected to the rod by a hinge (the breaking beam), restrained by springs of the stiffness of the breaking beam in torsion and bending, with the masses of the adapter and of the FPS weights as in the Nastran model
function create_wingtip_rod_beams_Pazy(mainWing,nElem,Λ,config,nElementsRod)

    # Get the wingtip rod model for the given configuration
    rod = wingtip_rod_Pazy(Λ,config)

    # Rotation tensor of the wing's beam basis from basis A, and rotation tensor from the Nastran frame to the wing's beam basis (as in `inertia_matrices_from_Nastran_data`)
    R0wing = mainWing.geometry.R0
    RNastran = [sin(Λ) cos(Λ) 0; -cos(Λ) sin(Λ) 0; 0 0 1]

    # Attachment point relative to the tip node, resolved in basis A
    offsetA = R0wing * (RNastran * rod.offset)

    # Rigid link from the tip node to the attachment point (centerline) of the rod: its x1 axis is along the offset, and x3 is as close as possible to that of the wing
    linkLength = norm(offsetA)
    x1link = offsetA / linkLength
    x2link = normalize(cross(R0wing[:,3], x1link))
    x3link = cross(x1link, x2link)
    linkStiffness = 1e12
    link = create_Beam(
        name="wingtipRodLink",
        length=linkLength,
        nElements=1,
        S=[isotropic_stiffness_matrix(∞=linkStiffness)],
        I=[inertia_matrix()],
        p0=rotation_parameters_WM([x1link x2link x3link])
        )

    # Positions of the nodes of the rod along the chordwise axis, relative to the tip node of the wing, and the attachment node. The part of the rod glued to the wing, from the start of the first rigid part to the end of the last one (of the two planes of bending), has an element between each pair of nodes at the ends of the rigid parts and at the attachment point, and the parts beyond it, in cantilever, are divided in proportion to their lengths
    xAttachment = rod.offset[1]
    @assert all(lim -> lim[1] < xAttachment < lim[2], rod.tied) "the attachment point of the rod must be inside its rigid parts"
    xRigid = sort(unique(vcat(collect(rod.tied.inPlane), collect(rod.tied.outOfPlane), xAttachment)))
    nRigid = length(xRigid) - 1
    nFree = nElementsRod - nRigid
    @assert nFree >= 2 "nElementsWingtipRod must be at least $(nRigid+2) for this rod"
    lengthStart = xRigid[1] - rod.xStart
    lengthEnd = rod.xEnd - xRigid[end]
    nStart = clamp(round(Int, nFree*lengthStart/(lengthStart+lengthEnd)), 1, nFree-1)
    nEnd = nFree - nStart
    xNodes = vcat(range(rod.xStart, xRigid[1], length=nStart+1)[1:end-1], xRigid, range(xRigid[end], rod.xEnd, length=nEnd+1)[2:end])
    attachmentNode = findfirst(x -> isapprox(x, xAttachment; atol=1e-9), xNodes)

    # Extension of the rod beyond its end up to the center of mass of the FPS weights, and the node of the hinge
    hinge = rod.hinge
    hingeNode = nothing
    if !isnothing(hinge)
        extension = hinge.bbLength + hinge.arm
        if hinge.side == :LE
            pushfirst!(xNodes, rod.xStart - extension)
            attachmentNode += 1
            hingeNode = 2
        else
            push!(xNodes, rod.xEnd + extension)
            hingeNode = length(xNodes) - 1
        end
    end
    nElements = length(xNodes) - 1
    rodLength = xNodes[end] - xNodes[1]
    normalizedNodalPositions = (xNodes .- xNodes[1]) / rodLength
    normalizedNodalPositions[end] = 1.0

    # Sectional matrices: rod (nylon), rigid in the part tied to the wing. The rod has its x2 axis along the span, so the bending in the plane of the wing is that of the shear GAy and bending EIz, and the one out of it, of GAz and EIy. Where only one of the planes is tied, the rod is rigid in it (and in axial and torsional deformation), and flexible in the other
    Srod = isotropic_stiffness_matrix(E=rod.E,G=rod.G,A=rod.A,Iy=rod.I,Iz=rod.I,J=rod.J,Ksy=rod.Ks,Ksz=rod.Ks)
    Srigid = isotropic_stiffness_matrix(∞=linkStiffness)
    SrigidInPlane = isotropic_stiffness_matrix(∞=linkStiffness,GAz=rod.G*rod.A*rod.Ks,EIy=rod.E*rod.I)
    SrigidOutOfPlane = isotropic_stiffness_matrix(∞=linkStiffness,GAy=rod.G*rod.A*rod.Ks,EIz=rod.E*rod.I)
    Irod = inertia_matrix(ρ=rod.ρ,A=rod.A,Iy=rod.I,Iz=rod.I)
    isTied(e,lim) = xNodes[e] >= lim[1] - 1e-9 && xNodes[e+1] <= lim[2] + 1e-9
    function elementStiffness(e)
        inPlane = isTied(e,rod.tied.inPlane)
        outOfPlane = isTied(e,rod.tied.outOfPlane)
        return inPlane && outOfPlane ? Srigid : inPlane ? SrigidInPlane : outOfPlane ? SrigidOutOfPlane : Srod
    end

    # Rod: x1 axis along the chord, positive aft (the rod's first node is at its leading edge end), x2 along the span and x3 vertical, as in the Nastran frame
    R0rod = R0wing * RNastran
    wingtipRod = create_Beam(
        name="wingtipRod",
        length=rodLength,
        firstNodeOffset=(xNodes[1]-xAttachment) * R0rod[:,1],
        nElements=nElements,
        normalizedNodalPositions=normalizedNodalPositions,
        S=[elementStiffness(e) for e in 1:nElements],
        I=fill(Irod,nElements),
        p0=rotation_parameters_WM(R0rod),
        hingedNodes=isnothing(hingeNode) ? Int[] : [hingeNode],
        hingedNodesDoF=isnothing(hingeNode) ? Vector{Vector{Bool}}() : [trues(3)]
        )

    # Calibration masses, at their position on the rod
    if rod.calibrationMass > 0
        calibrationElement = findfirst(e -> xNodes[e] <= rod.xCalibration <= xNodes[e+1], 1:nElements)
        add_point_inertia_to_beam!(
            beam=wingtipRod,
            elementID=calibrationElement,
            pointInertia=create_PointInertia(mass=rod.calibrationMass, η=[rod.xCalibration - (xNodes[calibrationElement]+xNodes[calibrationElement+1])/2, 0, 0])
            )
    end

    # Hinge of the breaking beam, with the adapter and the FPS weights on the element after it
    if !isnothing(hinge)
        # Material (PLA) of the breaking beam of the Sx LE and TE configurations
        breakingBeamE_Pazy = 2.3e9
        breakingBeamG_Pazy = 8.5185e8
        # Springs around the hinge, with the stiffnesses of the breaking beam (a cone of the PLA, whose length is that of the hinge)
        bbI = π*hinge.bbDiameter^4/64
        bbJ = π*hinge.bbDiameter^4/32
        kTorsion = breakingBeamG_Pazy * bbJ / hinge.bbLength
        kBending = breakingBeamE_Pazy * bbI / hinge.bbLength
        for (k,p) in zip((kTorsion,kBending,kBending), ([0,0,0],[π/2,0,0],[0,π/2,0]))
            add_spring_to_beams!(beams=[wingtipRod,wingtipRod],spring=create_Spring(kp=k),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=p)
        end
        # Element after the hinge, whose far end is the center of mass of the FPS weights (a node), and its midpoint
        extensionElement = hinge.side == :LE ? 1 : nElements
        xFar = hinge.side == :LE ? xNodes[1] : xNodes[end]
        xMid = (xNodes[extensionElement] + xNodes[extensionElement+1]) / 2
        extensionLength = abs(xNodes[extensionElement+1] - xNodes[extensionElement])
        # Adapter, at the end of the breaking beam
        xAdapter = hinge.side == :LE ? rod.xStart - hinge.bbLength : rod.xEnd + hinge.bbLength
        add_point_inertia_to_beam!(
            beam=wingtipRod,
            elementID=extensionElement,
            pointInertia=create_PointInertia(mass=hinge.adapterMass, η=[xAdapter-xMid, 0, 0])
            )
        # FPS weights: the element has the distributed mass and inertia of the rod, so a concentrated mass is added to complete the mass of the FPS weights and keep their center of mass, and the inertia of the FPS weights about it
        elementMass = rod.ρ * rod.A * extensionLength
        concentratedMass = hinge.fpsMass - elementMass
        @assert concentratedMass > 0 "the mass of the extension of the rod is larger than that of the FPS weights"
        xConcentrated = (hinge.fpsMass*xFar - elementMass*xMid) / concentratedMass
        elementInertia = elementMass*extensionLength^2/12 + elementMass*(xMid-xFar)^2 + concentratedMass*(xConcentrated-xFar)^2
        concentratedInertia = max(hinge.fpsInertia - elementInertia, 0)
        add_point_inertia_to_beam!(
            beam=wingtipRod,
            elementID=extensionElement,
            pointInertia=create_PointInertia(mass=concentratedMass, η=[xConcentrated-xMid, 0, 0], Iyy=concentratedInertia, Izz=concentratedInertia, Ixx=0)
            )
    end

    # Set connections (the rod is connected to the end of the link at the attachment node)
    connections = [
        (link => mainWing, 1 => nElem+1),
        (wingtipRod => link, attachmentNode => 2)
    ]

    return [mainWing,link,wingtipRod],connections

end

# Path of the CSV file with the original structural data of the straight Pazy wing (https://github.com/UM-A2SRL/AePW3-LDWG), for the :nominal mesh (test/referenceData/Pazy/S0_AePW3_<kind>_matrices.csv, with kind :inertia or :stiffness). `kind` is :inertia (nodal data in the MSC Nastran frame, same format as the beam mass model) or :stiffness
function Pazy_S0_AePW3_data_file(kind::Symbol)

    @assert kind in (:inertia, :stiffness) "kind must be :inertia or :stiffness"
    file = joinpath(pkgdir(AeroBeams), "test", "referenceData", "Pazy", "S0_AePW3_" * (kind == :inertia ? "inertia" : "stiffness") * "_matrices.csv")
    isfile(file) || throw(ArgumentError("original AePW3 $kind data of the straight Pazy wing not found (expected file $file)"))

    return file

end

# Sectional stiffness matrices of the straight Pazy wing from the original AePW3 data (with skin), one for each element of the :nominal mesh. The shear stiffnesses are set to 3e6 N, and the couplings not present in the data are zero
function stiffness_matrices_PazyS0_AePW3()

    # Geometrical properties
    nElem,L,_,_,nodalPositions = geometrical_properties_Pazy(mesh=:nominal, Λ=0)

    # Stiffness matrix of each element, given by its components row by row (S11, S12, ..., S66) and the coordinates of its ends from the root
    file = Pazy_S0_AePW3_data_file(:stiffness)
    data, header = readdlm(file, ','; header=true)
    @assert size(data,1) == nElem "the original AePW3 stiffness data has $(size(data,1)) rows, but $nElem are expected"
    @assert isapprox(Float64.(data[:, 2]), L*nodalPositions[1:end-1]; atol=1e-6*L) && isapprox(Float64.(data[:, 3]), L*nodalPositions[2:end]; atol=1e-6*L) "the original AePW3 stiffness data does not match the nodes of the :nominal mesh"
    columns = [findfirst(==("S$i$j"), vec(header)) for i in 1:6, j in 1:6]

    return [Float64.(data[k, columns]) for k in 1:nElem]

end

# Computes the tip loss function for the Pazy wing's tip correction function using a fixed exponential decay model
function tip_loss_fixed_exponential(type)

    # Extract tip loss decay factor
    m = match(r"FixedExponential-([\d.]+)", type)
    τ = m !== nothing ? parse(Float64, m[1]) : nothing
    @assert !isnothing(τ)
    @assert τ > 0

    # Set tip loss function
    return s -> 1 - exp(-τ*(1-s))

end

# Computes the tip loss function for the Pazy wing's tip correction function using a root pitch and airspeed dependent exponential decay model
function tip_loss_exponential(θ,U)

    # Bound inputs
    θ = clamp(θ, 0, 7π/180)
    U = clamp(U, 0, 60)

    # Coefficients as a function of root pitch angle and airspeed
    θRange = deg2rad.([0; 1; 2; 3; 4; 5; 6; 7])
    τ₀Range = [6.58; 6.29; 6.00; 5.92; 5.91; 6.08; 6.43; 6.88]
    τ₁Range = 1e-2*[0; 3.33; 5.19; 6.01; 6.49; 5.98; 4.43; 2.30]
    τ₂Range = -1e-4*[0; 5.56; 8.59; 10.6; 12.3; 12.8; 11.9; 10.2]
    τ₀ = LinearInterpolations.interpolate(θRange,τ₀Range,θ)
    τ₁ = LinearInterpolations.interpolate(θRange,τ₁Range,θ)
    τ₂ = LinearInterpolations.interpolate(θRange,τ₂Range,θ)    
    τ = τ₀ + τ₁*U + τ₂*U^2

    # Set tip loss function
    return s -> 1 - exp(-τ*(1-s))

end

# Computes the tip loss function for the Pazy wing's tip correction function using VLM data
function tip_loss_vlm(type,θ,U,Λ)

    Λc = clamp(Λ, _VLM_cache.ΛRange[1], _VLM_cache.ΛRange[end])
    q = zeros(9)

    if type == "VLM-undef"
        for i in 1:9
            q[i] = _VLM_cache.itpsUndef[i](Λc)
        end
        # Compressibility correction
        scale = sqrt(1 - (U/343)^2)
        q .*= scale
    elseif type == "VLM-def"
        θc = clamp(θ, _VLM_cache.θRange[1], _VLM_cache.θRange[end])
        Uc = min(U, _VLM_cache.URange[end])
        for i in 1:9
            q[i] = _VLM_cache.itpsDef[i]((Λc, θc, Uc))
        end
    else
        error("Unknown VLM type: $type")
    end

    # Return polynomial (Horner evaluation)
    return s -> begin
        p = q[9]
        @inbounds for k=8:-1:1
            p = p*s + q[k]
        end
        p
    end
    
end


## PUBLIC API
# ------------------------------------------------------------------------------
"""
    typical_section_data(name::String)

Returns the data of the typical section of given name

# Arguments
- `name::String`: name of the typical section

# Outputs
- `a`: position of elastic axis (semichords aft of midchord)
- `e`: position of mass axis (semichords aft of midchord)
- `μ`: mass ratio
- `rα²`: squared radius of gyration in pitch 
- `σ`: ratio of plunge to pitch frequency in vacuo
- `ωα`: pitch frequency in vacuo [rad/s]
- `c`: chord [m]
"""
function typical_section_data(name::String)

    if name == "HP-1" # Defined in problem 5.5 of Hodges & Pierce
        a = -1/5
        e = -1/10
        μ = 20
        rα² = 6/25
        σ = 2/5
        ωα = 30 # A numerical value was not given
        c = 1   # A numerical value was not given
    elseif name == "HP-2" # Defined in problem 5.9 of Hodges & Pierce
        a = -1/3
        e = -1/10
        μ = 50
        rα² = 4/25
        σ = 2/5
        ωα = 30 # A numerical value was not given
        c = 1   # A numerical value was not given
    elseif name == "Fung" # Defined by Fung (Sec. 6.9) - An Introduction to Theory of Aeroelasticity
        a = -0.15
        e = 0.1
        μ = 76
        rα² = 0.6229^2
        ωα = 64.1
        ωh = 55.9
        σ = ωh/ωα
        c = 0.127*2
    end

    return a,e,μ,rα²,σ,ωα,c

end
export typical_section_data


"""
    geometrical_properties_Pazy(; mesh, Λ, useAePW3Matrices)

Returns the fixed geometrical (and discretization) properties of the Pazy wing (sweep angles of 0, 10 and 20 degrees), computed from the reference axis, which is at the shear center (see `reference_axis_Pazy`), whose length is the total length of the beam

# Keyword Arguments
- `mesh::Symbol`: beam mesh, `:nominal` (nodes at the root, at the ribs and at the tip of the spar) or `:refined` (with a node at the middle of each element of the nominal mesh, approximately twice as many elements) (default: `:nominal`)
- `Λ::Real`: sweep angle [rad] (default: 0)
- `useAePW3Matrices::Bool`: flag to return the normalized spar position of the original AePW3 sectional matrices (0.44096), instead of that of the shear center. Only available for the straight wing (default: `false`)

# Outputs
- `nElem::Int`: number of elements
- `L::Float64`: total length
- `chord::Float64`: chord
- `normSparPos::Float64`: normalized spar position, equal to the chordwise position of the shear center (0.4653 for sweep 0°, 0.447 for 10° and 0.440 for 20°), or to that of the original AePW3 reference axis (0.44096) if `useAePW3Matrices` is `true`
- `nodalPositions::Vector{Float64}`: normalized nodal positions
"""
function geometrical_properties_Pazy(; Λ=0, mesh=:nominal, useAePW3Matrices=false)

    validate_mesh_Pazy(mesh)
    @assert any(isapprox.(rad2deg(Λ), [0,10,20]; atol=1e-8)) "sweepAngle must be 0, 10 or 20 degrees"
    @assert !useAePW3Matrices || Λ ≈ 0 "the original AePW3 matrices are only available for the straight wing"

    # Chord [m]
    chord = chord_Pazy(Λ)

    # Reference axis coordinates (chordwise from the leading edge, spanwise)
    x,y,_,ra = reference_axis_Pazy(Λ; mesh=mesh)

    # Length of each segment and total length [m]
    segmentLengths = sqrt.(diff(x).^2 .+ diff(y).^2)
    L = sum(segmentLengths)

    # Normalized spar position (the original AePW3 data refers to a reference axis at a different chordwise position than the shear center of the condensed model)
    normSparPos = useAePW3Matrices ? normSparPos_Pazy_AePW3 : ra

    # Normalized nodal positions
    nodalPositions = [0.0; cumsum(segmentLengths)]/L
    nodalPositions[end] = 1.0 # avoid rounding error

    # Number of elements
    nElem = length(nodalPositions) - 1

    return nElem,L,chord,normSparPos,nodalPositions

end
export geometrical_properties_Pazy


"""
    get_sweptPazy_tip_inertias_from_baseline(; Λ, configs=(:baseline, :LE, :TE), withFPSMass=true, withAdapterMass=true, withCalibrationMass=true, sweptTipMassInertiaCorrection=false) -> (tipMasses, ηtipMasses, tipMassesInertias, tipMassesElementIDs)

Get tip masses and inertias for the swept Pazy wing LE and TE configurations from the baseline configuration.

# Keyword arguments
- `Λ::Float64`: sweep angle [rad]
- `configs::Tuple`: tip mass configurations to consider (default: `(:baseline, :LE, :TE)`)
- `withFPSMass::Bool`: include flight protection system (FPS) mass (default: `true`)
- `withAdapterMass::Bool`: include adapter mass (default: `true`)
- `withCalibrationMass::Bool`: include calibration mass (default: `true`)
- `sweptTipMassInertiaCorrection::Bool`: apply correction to tip mass inertias (default: `false`). This correction is roughly equivalent on setting the adapter and FPS masses 20% further from the reference axis.

# Returns
- `tipMasses::Vector{Vector{Float64}}`: tip masses for each configuration
- `ηtipMasses::Vector{Vector{Vector{Float64}}}`: tip mass positions for each configuration
- `tipMassesInertias::Vector{Vector{Matrix{Float64}}}`: tip mass inertias for each configuration
- `tipMassesElementIDs::Vector{Vector{Int}}`: IDs of the elements to which the tip masses of each configuration must be attached (pass it to `create_Pazy`). Following the NASTRAN beam mass model, each mass belongs to the element ending at the beam node nearest to it in the unswept spanwise coordinate. The positions `ηtipMasses` remain relative to the tip node.
"""
function get_sweptPazy_tip_inertias_from_baseline(; Λ, configs=(:baseline, :LE, :TE), withFPSMass=true, withAdapterMass=true, withCalibrationMass=true, sweptTipMassInertiaCorrection=false)

    # Validate
    Λdeg = rad2deg(Λ)
    Λis10deg = Λdeg ≈ 10
    Λis20deg = Λdeg ≈ 20
    @assert Λis10deg || Λis20deg "Sweep angle Λ must be 10 or 20 degrees (in rad)"
    for config in configs
        @assert config in (:baseline, :LE, :TE) "Invalid configuration: $config"
    end

    # Set masses, positions and inertias (in the MSC NASTRAN frame) for the tip masses in the LE and TE configurations
    if Λis10deg
        fpsMassLE = fpsMassTE = 8.7e-3
        adapterMassLE = adapterMassTE = 3.5e-3
        calibrationMassLE = calibrationMassTE = 5e-3
        fpsMassPosLE = [-106.99e-3, 4.87e-3, 0]
        fpsMassPosTE = [117.61e-3, 4.66e-3, 0]
        adapterMassPosLE = [-85.69e-3, 4.87e-3, 0]
        adapterMassPosTE = [95.31e-3, 4.66e-3, 0]
        calibrationMassPosLE = [0.31e-3, -23.54e-3, 0]
        calibrationMassPosTE = [-26.45e-3, -23.54e-3, 0]
    else 
        fpsMassLE = fpsMassTE = 11.3e-3
        adapterMassLE = adapterMassTE = 3.1e-3
        calibrationMassLE = 7e-3
        calibrationMassTE = 3e-3
        fpsMassPosLE = [-108.54e-3, 4.65e-3, 0]
        fpsMassPosTE = [110.35e-3, 4.45e-3, 0]
        adapterMassPosLE = [-93.93e-3, 4.65e-3, 0]
        adapterMassPosTE = [90.88e-3, 4.45e-3, 0]
        calibrationMassPosLE = [-83.97e-3, 4.65e-3, 0]
        calibrationMassPosTE = [-83.97e-3, 4.65e-3, 0]
    end
    fpsMassInertia = sweptTipMassInertiaCorrection ? 
                        diagm([0, 3e-5, 3.3e-7]) : 
                        diagm([0, 3.3e-7, 3.3e-7])
    adapterMassInertia = sweptTipMassInertiaCorrection ? 
                         diagm([0, 1.7e-5, 0]) : 
                         zeros(3,3)
    fpsMassInertiaLE = fpsMassInertia
    fpsMassInertiaTE = fpsMassInertia
    adapterMassInertiaLE = adapterMassInertia
    adapterMassInertiaTE = adapterMassInertia
    calibrationMassInertiaLE = zeros(3,3)
    calibrationMassInertiaTE = zeros(3,3)

    # Rotation tensor from the Nastran frame (x chordwise positive aft, y spanwise unswept, z vertical) to the local beam frame (basis b)
    R = [sin(Λ) cos(Λ) 0; -cos(Λ) sin(Λ) 0; 0 0 1]

    # Element that carries a tip mass at position `posNastran` (relative to the tip beam node, in the Nastran frame): the one ending at the beam node nearest to the mass in the unswept y coordinate, as in the NASTRAN beam mass model
    nElem = geometrical_properties_Pazy(Λ=Λ)[1]
    _,yRef,_ = reference_axis_Pazy(Λ; mesh=:nominal)
    host_element(posNastran) = clamp(argmin(abs.(yRef .- (yRef[end] + posNastran[2]))) - 1, 1, nElem)

    # Initialize tip mass and inertia containers for each configuration
    tipMasses = Vector{Vector{Float64}}(undef,length(configs))
    ηtipMasses = Vector{Vector{Vector{Float64}}}(undef,length(configs))
    tipMassesInertias = Vector{Vector{Matrix{Float64}}}(undef,length(configs))
    tipMassesElementIDs = Vector{Vector{Int}}(undef,length(configs))

    # Loop over each configuration and populate the tip mass, position, and inertia arrays, rotating the positions and inertias to the local beam frame
    for (i,config) in enumerate(configs)
        m = Vector{Float64}()
        η = Vector{Vector{Float64}}()
        I = Vector{Matrix{Float64}}()
        elementIDs = Vector{Int}()
        if config == :baseline
            push!(m, 0.0)
            push!(η, zeros(3))
            push!(I, zeros(3,3))
            push!(elementIDs, nElem)
        elseif config == :LE
            if withFPSMass
                push!(m, fpsMassLE)
                push!(η, R*fpsMassPosLE)
                push!(I, R*fpsMassInertiaLE*R')
                push!(elementIDs, host_element(fpsMassPosLE))
            end
            if withAdapterMass
                push!(m, adapterMassLE)
                push!(η, R*adapterMassPosLE)
                push!(I, R*adapterMassInertiaLE*R')
                push!(elementIDs, host_element(adapterMassPosLE))
            end
            if withCalibrationMass
                push!(m, calibrationMassLE)
                push!(η, R*calibrationMassPosLE)
                push!(I, R*calibrationMassInertiaLE*R')
                push!(elementIDs, host_element(calibrationMassPosLE))
            end
        elseif config == :TE
            if withFPSMass
                push!(m, fpsMassTE)
                push!(η, R*fpsMassPosTE)
                push!(I, R*fpsMassInertiaTE*R')
                push!(elementIDs, host_element(fpsMassPosTE))
            end
            if withAdapterMass
                push!(m, adapterMassTE)
                push!(η, R*adapterMassPosTE)
                push!(I, R*adapterMassInertiaTE*R')
                push!(elementIDs, host_element(adapterMassPosTE))
            end
            if withCalibrationMass
                push!(m, calibrationMassTE)
                push!(η, R*calibrationMassPosTE)
                push!(I, R*calibrationMassInertiaTE*R')
                push!(elementIDs, host_element(calibrationMassPosTE))
            end
        end
        tipMasses[i] = m
        ηtipMasses[i] = η
        tipMassesInertias[i] = I
        tipMassesElementIDs[i] = elementIDs
    end

    return tipMasses, ηtipMasses, tipMassesInertias, tipMassesElementIDs

end
export get_sweptPazy_tip_inertias_from_baseline


"""
    tip_loss_function_Pazy(type::String,θ::Real,U::Real,Λ::Real)

Computes the tip loss function for the Pazy wing's tip correction function

# Arguments
- `type::String`: tip loss function type
- `θ::Real`: root pitch angle [rad]
- `U::Real`: airspeed [m/s]
- `Λ::Real`: sweep angle [rad]
"""
function tip_loss_function_Pazy(type::String,θ::Real,U::Real,Λ::Real)

    # Fixed exponential decay tip loss
    if occursin("FixedExponential", type)
        return tip_loss_fixed_exponential(type)
    # Exponential decay tip loss (root pitch and airspeed dependent)
    elseif type == "Exponential"
        return tip_loss_exponential(θ,U)
    # VLM: interpolate 8th order polynomial   
    elseif startswith(type, "VLM")
        return tip_loss_vlm(type, θ, U, Λ)
    # No tip loss function    
    elseif type == "None"
        return s -> 1.0
    else
        error("Unknown tip loss function type: $type")
    end

end
export tip_loss_function_Pazy


"""
    tip_loss_function_HealyLCOFFWT(θ::Real,U::Real)

Computes the tip loss function for Healy's wing with flared folding wingtip used in LCO studies

# Arguments
- `θ::Real`: root pitch angle [rad]
- `U::Real`: airspeed [m/s]
"""
function tip_loss_function_HealyLCOFFWT(θ::Real,U::Real)

    # Bound inputs
    θ = clamp(rad2deg(θ), 2.5, 10)
    U = min(30,U)

    # Coefficients as a function of root pitch angle
    θRange  = [2.5; 5.0; 7.5; 10]
    τ₀Range = [15; 13; 12; 11]
    τ₁Range = 1/30*[5; 7; 10; 10]
    τ₀ = LinearInterpolations.interpolate(θRange,τ₀Range,θ)
    τ₁ = LinearInterpolations.interpolate(θRange,τ₁Range,θ)

    # Set tip loss function
    return τ₀ + τ₁*U

end
export tip_loss_function_HealyLCOFFWT


"""
    create_Pazy(; kwargs...)

Creates the Pazy wing

Returns the wing model and geometrical properties

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `upright::Bool`: flag to set the wing in the upright position
- `airfoil::Airfoil`: airfoil section
- `pitchAngle::Real`: pitch angle at the root [rad]
- `sweepAngle::Real`: sweep angle [rad]
- `hasTipCorrection::Bool`: flag for aerodynamic tip correction 
- `altitude::Real`: altitude
- `gravityVector::Real`: acceleration of gravity vector (resolved in the inertial basis I)
- `airspeed::Union{Real,<:Function}`: airspeed
- `smallAngles::Bool`: flag for small angles approximation
- `tipLossType::String`: tip loss function type
- `tipLossFunctionIsAirspeedDependent::Bool`: flag to update the tip loss function with the instantaneous airspeed along the problem
- `mesh::Symbol`: beam mesh, can be `:nominal` (nodes at the root, at the ribs and at the tip of the spar) or `:refined` (approximately twice as many elements, with a node at the middle of each element of the nominal mesh)
- `useAePW3Matrices::Bool`: flag to use the original sectional stiffness and inertia matrices of the straight wing (with skin, from https://github.com/UM-A2SRL/AePW3-LDWG, with shear stiffnesses of 3e6 N, and with the reference axis at 0.44096 of the chord) instead of those obtained from the condensation of the finite element model. Only available for the straight wing with `mesh=:nominal` (default: `false`)
- `withWingtipRod::Bool`: flag to include the wingtip rod (its inertia is lumped in the tip node of the straight wing's beam mass model). If `false`, the wing has no rod, and its inertia and stiffness matrices are read from `test/referenceData/PazyNoRod` (least squares stiffness matrices only). Only available for the straight wing (default: `true`)
- `rootTorsionCorr::Bool`: flag to use the sectional stiffness matrices with the correction of the root torsional stiffness (files `S0_stiffness_matrices_<mesh>mesh_<tip>_corr.csv` in `test/referenceData/Pazy`, instead of those of the condensation of the finite element model, `S0_stiffness_matrices_<mesh>mesh_<tip>.csv`). The rigid clamp of the Nastran model restrains the warping and stiffens the root in torsion, and the real clamp is more flexible. The correction reduces S44 of the root element (and also changes the bending-torsion coupling S45 = S54 for the nominal mesh), calibrated so that the flutter onset speeds of the hump mode agree with the onset (up) and onset (down) speeds measured by Drachinski et al. It reduces the frequency of the first torsion mode by about 1% (default: `true`). Only available for the straight wing with the tip rod (S0), and ignored for other configurations and for `useAePW3Matrices=true`
- `beamWingtipRod::Bool`: flag to model the wingtip rod as a beam attached to the tip of the wing (`true`), instead of lumping its inertia in the tip node (`false`, the conventional model). Available for the straight wing (S0) and for the wing with 10 degrees of sweep (S10), and requires `withWingtipRod=true`. The wing beam then has the inertia data of the wing without the rod (`test/referenceData/PazyNoRod`), and the stiffness data of the Nastran model with the flexible tip rod (`flextip`) for the straight wing, and of the one with the rigid tip rod for the S10, and the model has two more beams, in this order: a rigid link, from the tip node of the wing to the attachment point of the rod, and the rod (flexible, with the properties of the Nastran models), connected to the end of the link. The part of the rod glued to the wing along the chord is rigid (see `wingtip_rod_Pazy` for its limits, which for the straight wing are not the same for the bending in the plane of the wing and out of it), and the parts beyond it are in cantilever. For the straight wing, the rod carries the calibration masses. For the S10, the calibration mass is part of the wing; in the `:LE` and `:TE` configurations of `sweptTipMassConfig`, the rod is extended with an element up to the center of mass of the FPS weights, which is connected to the rod by a hinge (the breaking beam) with springs of the stiffness of the breaking beam in torsion and bending, and has the masses of the fixed adapter and of the FPS weights (default: `false`)
- `nElementsWingtipRod::Int`: number of elements of the wingtip rod beam (at least 6 for the straight wing and 4 for the S10), without the element of the FPS weights of the S10 `:LE` and `:TE` configurations. Those in the part of the rod glued to the wing are fixed (4 for the straight wing, which has different rigid parts in the two planes of bending, and 2 for the S10), and the others are distributed over the parts beyond it in proportion to their lengths (default: 12)
- `tipMasses::Vector{<:Real}`: masses of point inertias added to the tip of the wing
- `ηtipMasses::Vector{<:Real}`: position vectors of the tip masses relative to the tip of the spar, resolved in the local basis
- `tipMassesInertias::Vector{<:Matrix{<:Real}}`: mass moment of inertia matrices of the tip masses
- `tipMassesElementIDs::Union{Nothing,Vector{<:Integer}}`: IDs of the elements to which the tip masses are attached (default: `nothing`, which attaches all tip masses to the last element). The positions `ηtipMasses` are always relative to the tip of the spar, regardless of the element to which each mass is attached
- `additionalBCs::Vector{BC}`: additional BCs (beyond the clamp)
- `flightDirection::Vector{<:Real}`: flight direction vector, resolved in basis A
- `tipLossPitchAngle::Real`: pitch angle at the root used to evaluate the tip loss function [rad] (default: `pitchAngle`). Useful when the angle of attack is set through `flightDirection` instead of `pitchAngle`
- `gust::Union{Nothing,Gust}`: gust
- `hinged::Bool`: flag for hinged wing
- `hingeNode::Union{Nothing,Int}`: hinge node
- `foldAngle::Union{Nothing,Real}`: fold angle [rad] (for hinge)
- `flareAngle::Union{Nothing,Real}`: flare angle [rad] (for hinge)
- `kSpring::Union{Nothing,Real}`: stiffness of the spring around the hinge
- `sweptTipMassConfig::Symbol`: configuration of the swept tip mass, can be `:baseline`, `:TE` or `:LE` (with `beamWingtipRod=true`, this sets the rod and the masses at its end)
- `sweptTipMassInertiaCorrection::Bool`: flag to apply an ad hoc factor on the tip node's torsional inertia for the `:LE` and `:TE` swept tip mass configurations (approximates the GVT torsional frequencies, which are lowered by the tip rod flexibility not captured by the stiffened-tip NASTRAN model from which the swept wing mass data was derived)
"""
function create_Pazy(; 
    aeroSolver::AeroSolver=Indicial(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    updateAirfoilParameters::Bool=true,
    upright::Bool=false,
    airfoil::Airfoil=deepcopy(NACA0018),
    pitchAngle::Real=0,
    sweepAngle::Real=0,
    hasTipCorrection::Bool=true,
    altitude::Real=0,
    gravityVector::AbstractVector{<:Real}=[0,0,-standard_atmosphere(altitude).g],
    airspeed::Union{Real,<:Function}=0,
    smallAngles::Bool=false,
    tipLossType::String="Exponential",
    tipLossFunctionIsAirspeedDependent::Bool=false,
    mesh::Symbol=:nominal,
    useAePW3Matrices::Bool=false,
    withWingtipRod::Bool=true,
    rootTorsionCorr::Bool=true,
    beamWingtipRod::Bool=false,
    nElementsWingtipRod::Int=12,
    tipMasses::AbstractVector{<:Real}=zeros(0),
    ηtipMasses::AbstractVector{<:AbstractVector{<:Real}}=Vector{Vector{Float64}}(),
    tipMassesInertias::AbstractVector{<:AbstractMatrix{<:Real}}=Vector{Matrix{Float64}}(),
    tipMassesElementIDs::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    additionalBCs::Vector{BC}=Vector{BC}(),
    flightDirection::AbstractVector{<:Real}=a2,
    tipLossPitchAngle::Real=pitchAngle,
    gust::Union{Nothing,Gust}=nothing,
    hinged::Bool=false,
    hingeNode::Union{Nothing,Int} = hinged ? 14 : nothing,
    foldAngle::Union{Missing,Real}=missing,
    flareAngle::Union{Nothing,Real} = hinged ? 0 : nothing,
    kSpring::Union{Nothing,Real} = hinged ? 1e-5 : nothing,
    sweptTipMassConfig::Symbol=:baseline,
    sweptTipMassInertiaCorrection::Bool=false,
    )

    # Validate
    @assert altitude >= 0 "altitude must be non-negative"
    @assert length(gravityVector) == 3 "gravityVector must be a 3-element vector"
    @assert any(isapprox.(rad2deg(sweepAngle), [0,10,20]; atol=1e-8)) "sweepAngle must be 0, 10 or 20 degrees"
    validate_mesh_Pazy(mesh)
    @assert length(tipMasses) == length(ηtipMasses) == length(tipMassesInertias) "tipMasses, ηtipMasses, and tipMassesInertias must have the same length"
    @assert all(tipMasses .>= 0) "tip masses must be non-negative"
    for ηtipMass in ηtipMasses
        @assert length(ηtipMass) == 3 "each ηtipMass must be a 3-element vector"
    end
    for tipMassInertia in tipMassesInertias
        @assert size(tipMassInertia) == (3,3) "each tipMassInertia must be a 3x3 matrix"
    end
    @assert isnothing(tipMassesElementIDs) || length(tipMassesElementIDs) == length(tipMasses) "tipMassesElementIDs must have the same length as tipMasses"
    if !isnothing(tipMassesElementIDs)
        nElem = first(geometrical_properties_Pazy(mesh=mesh))
        @assert all(1 .<= tipMassesElementIDs .<= nElem) "tip mass element IDs must be between 1 and $nElem"
    end
    if hasTipCorrection
        @assert tipLossType in ["Exponential","VLM-undef","VLM-def"] || occursin("FixedExponential", tipLossType) "tipLossType must be one of \"Exponential\", \"VLM-undef\", \"VLM-def\", or \"FixedExponential-<value>\""
    end
    if hinged
        @assert 2 <= hingeNode <= 15 "hingeNode must be between 2 and 15"
        if !ismissing(foldAngle)
            @assert -π < foldAngle <= π "set foldAngle between -π and π (rad)"
        end
        @assert 0 <= flareAngle < π/4 "set flareAngle between 0 and π/4 (rad)"
        @assert kSpring >= 0 "kSpring must be non-negative"
    end
    @assert sweptTipMassConfig in [:baseline, :TE, :LE] "sweptTipMassConfig must be either :baseline, :TE or :LE"
    if useAePW3Matrices
        @assert sweepAngle ≈ 0 && mesh == :nominal "the original AePW3 matrices are only available for the straight wing with mesh=:nominal"
    end
    if beamWingtipRod
        @assert withWingtipRod "beamWingtipRod=true requires withWingtipRod=true"
        @assert (sweepAngle ≈ 0 || rad2deg(sweepAngle) ≈ 10) && !useAePW3Matrices "the wingtip rod modeled as a beam is only available for the straight wing and the wing with 10 degrees of sweep, and not with the original AePW3 matrices"
        @assert !sweptTipMassInertiaCorrection "sweptTipMassInertiaCorrection cannot be used with beamWingtipRod=true, as the rod flexibility is modeled"
        @assert nElementsWingtipRod >= 4 "nElementsWingtipRod must be at least 4"
    end
    if !withWingtipRod
        @assert sweepAngle ≈ 0 && !useAePW3Matrices "the wing without the tip rod is only available for the straight wing, and not with the original AePW3 matrices"
    end

    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,airspeed(0)/atmosphere.a)

    # Tip loss function, as a function of the airspeed
    tipLossFunction = U -> tip_loss_function_Pazy(tipLossType,tipLossPitchAngle,U,sweepAngle)

    # Fixed geometrical and discretization properties
    nElem,L,chord,normSparPos,nodalPositions = geometrical_properties_Pazy(mesh=mesh,Λ=sweepAngle,useAePW3Matrices=useAePW3Matrices)

    # Stiffness matrices
    S = if useAePW3Matrices
        stiffness_matrices_PazyS0_AePW3()
    else
        stiffness_matrices_Pazy(sweepAngle; mesh=mesh, withWingtipRod=withWingtipRod, tip=(beamWingtipRod && sweepAngle ≈ 0) ? :flextip : :rigidtip, rootTorsionCorr=rootTorsionCorr && sweepAngle ≈ 0 && withWingtipRod)
    end

    # Inertia matrices
    I = if useAePW3Matrices
        inertia_matrices_Pazy(sweepAngle; mesh=:nominal, file=Pazy_S0_AePW3_data_file(:inertia))
    else
        inertia_matrices_Pazy(sweepAngle; mesh=mesh, config=sweptTipMassConfig, applyTipInertiaCorrection=sweptTipMassInertiaCorrection, withWingtipRod=withWingtipRod && !beamWingtipRod)
    end

    # Rotation parametrization and rotation parameters from basis A to basis b(0)
    rotationParametrization = upright ? 
                              :E231 : 
                              :E321
    p0 = upright ? 
         [-π/2, -sweepAngle, pitchAngle] : 
         [-sweepAngle, 0, pitchAngle]

    # Aerodynamic surface
    surf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=airfoil,
        c=chord,
        normSparPos=normSparPos,
        hasTipCorrection=hasTipCorrection,
        tipLossFunction=tipLossFunction,
        updateAirfoilParameters=updateAirfoilParameters,
        Λ=sweepAngle,
        smallAngles=smallAngles,
        tipLossFunctionIsAirspeedDependent=hasTipCorrection,
        updateTipLossWithAirspeed=tipLossFunctionIsAirspeedDependent
        )

    # Wing beam
    mainWing = create_Beam(
        name="wingSpar",
        length=L,
        nElements=nElem,
        normalizedNodalPositions=nodalPositions,
        S=S,
        I=I,
        p0=transform_to_WM(p0; from=rotationParametrization),
        surface=surf,
        hingedNodes= hinged ? [hingeNode] : Int[],
        hingedNodesDoF= hinged ? [trues(3)] : Vector{Vector{Bool}}()
        )

    # Tip masses (attached to the last element unless specified otherwise)
    elementIDsTipMasses = isnothing(tipMassesElementIDs) ? 
                          fill(nElem,length(tipMasses)) : 
                          tipMassesElementIDs
    for (m, η, I, elementID) in zip(tipMasses, ηtipMasses, tipMassesInertias, elementIDsTipMasses)
        # Position of the tip mass relative to the element's midpoint: distance from the tip to the element's end plus half of its length
        ΔL = L*(nodalPositions[elementID+1]-nodalPositions[elementID])
        tipToElementMidpoint = L*(nodalPositions[end]-nodalPositions[elementID+1]) + ΔL/2
        tipInertia = create_PointInertia(
            mass=m,
            η=η+[tipToElementMidpoint, 0, 0],
            inertiaMatrix=I
            )
        add_point_inertia_to_beam!(beam=mainWing, elementID=elementID, pointInertia=tipInertia)
    end

    # Wingtip rod modeled as beams (rigid link and rod, with their connections to the wing)
    beams = [mainWing]
    connections = Tuple{Pair{Beam,Beam},Pair{Int,Int}}[]
    if beamWingtipRod
        beams,connections = create_wingtip_rod_beams_Pazy(mainWing,nElem,sweepAngle,sweptTipMassConfig,nElementsWingtipRod)
    end

    # Update beam of additional BCs, if applicable
    for BC in additionalBCs
        reset_BC_beam!(BC,mainWing)
    end
    BCs = additionalBCs

    # Add root clamp to BCs
    clamp = create_BC(
        name="clamp",
        beam=mainWing,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0)
        )
    push!(BCs,clamp)

    # Set basis A velocity vector
    v_A(t) = airspeed(t) * flightDirection / norm(flightDirection)

    # Hinge axis constraint and spring around hinge
    if hinged
        # Hinge axis (defined in the local, undeformed beam basis)
        localHingeAxis = rotation_tensor_E321([-flareAngle; 0; 0]) * a2
        # Prescribed hinge rotation value
        pHValue = !isnothing(foldAngle) ? 4*tan(foldAngle/4) : missing
        # Set hinge axis constraint
        hingeAxisConstraint = create_HingeAxisConstraint(
                                  beam=mainWing,
                                  localHingeAxis=localHingeAxis,
                                  pHValue=pHValue
                                  )
        # Springs around hinge
        if kSpring > 0
            spring = create_Spring(kp=kSpring)
            add_spring_to_beams!(beams=[mainWing,mainWing],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1])
            add_spring_to_beams!(beams=[mainWing,mainWing],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[π/2,0,0])
            add_spring_to_beams!(beams=[mainWing,mainWing],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[0,π/2,0])
        end
    end
    hingeAxisConstraints = hinged ? 
                           [hingeAxisConstraint] : 
                           HingeAxisConstraint[]

    # Model
    pazy = create_Model(
        name="Pazy",
        beams=beams,
        connections=connections,
        BCs=BCs,
        gravityVector=gravityVector,
        v_A=v_A,
        atmosphere=atmosphere,
        gust=gust,
        hingeAxisConstraints=hingeAxisConstraints,
        units=create_UnitsSystem(frequency="Hz")
        )

    return pazy,nElem,L,chord,normSparPos

end
export create_Pazy


"""
    create_SMW(; kwargs...)

Creates the 16-meter wing (SMW) of the conventional HALE aircraft described by Patil, Hodges and Cesnik in: Nonlinear Aeroelasticity and Flight Dynamics of HALE (2001)

Returns the wing model and its span

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `flapSolver::FlapAeroSolver`: aerodynamic solver for flap loads
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `pitchAngle::Real`: pitch angle [rad]
- `k1::Real`: twisting curvature
- `k2::Real`: flapwise bending curvature
- `nElem::Int`: number of elements for discretization
- `altitude::Real`: altitude
- `airspeed::Union{Real,<:Function}`: airspeed
- `g::Real`: acceleration of gravity
- `stiffnessFactor::Real`: stiffness factor for beam structural properties
- `∞::Real`: value of rigid structural properties
- `Ψ::Real`: torsion-IP-bending stiffness coupling factor
- `tipF3::Union{Real,<:Function}`: tip dead transverse force applied at the tip
- `cd0::Real`: parasite drag coefficient for the wing
- `cnδ::Real`: cn vs δ slope for the wing
- `cmδ::Real`: cm vs δ slope for the wing
- `cdδ::Real`: cd vs δ slope for the wing
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters
- `hasInducedDrag::Bool`: flag to include induced drag
- `δAil::Union{Nothing,Real,<:Function}`: aileron deflection
- `additionalBCs::Vector{BC}`: additional BCs (beyond the clamp and tip load)
"""
function create_SMW(; 
    aeroSolver::AeroSolver=Indicial(),
    flapSolver::FlapAeroSolver=TableLookup(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(cHALEairfoil),
    pitchAngle::Real=0,
    k1::Real=0,
    k2::Real=0,
    nElem::Int=32,
    altitude::Real=0,
    airspeed::Union{Real,<:Function}=0,
    g::Real=standard_atmosphere(altitude).g,
    stiffnessFactor::Real=1,
    torsionStiffnessFactor::Real=1,
    ∞::Real=1e12,
    Ψ::Real=0,
    tipF3::Union{Real,<:Function}=0,
    cd0::Real=0,
    cnδ::Real=2.5,
    cmδ::Real=-0.35,
    cdδ::Real=0.15,
    updateAirfoilParameters::Bool=false,
    hasInducedDrag::Bool=false,
    δAil::Union{Nothing,Real,<:Function}=nothing,
    additionalBCs::Vector{BC}=Vector{BC}()
    )

    # Validate
    @assert -π/2 < pitchAngle < π/2 "set pitchAngle between -π/2 and π/2 (rad)"
    @assert nElem > 1 "nElem must be greater than 1"
    @assert altitude >= 0 "altitude must be non-negative"
    @assert stiffnessFactor > 0 "stiffnessFactor must be positive"
    @assert torsionStiffnessFactor > 0 "torsionStiffnessFactor must be positive"
    if δAil isa Real
        δAilConst = δAil
        δAil = t -> δAilConst
    end

    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Create new airfoil with modified database
    flapSiteID = 75
    db = modified_database(
        airfoil.name;
        attached = (
            cd₀ = cd0,
        ),
        flap = (
            flapSiteID => create_FlapParameters(cdδ,cmδ,cnδ),
        )
    )
    airfoil = create_Airfoil(database=db,Ma=airspeed(0)/atmosphere.a)

    # Wing surface
    chord = 1
    normSparPos = 0.5
    normFlapPos = flapSiteID/100
    ailSize = 0.25
    normFlapSpan = [1-ailSize, 1]

    wingSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        flapSolver=flapSolver,
        airfoil=airfoil,
        c=chord,
        normSparPos=normSparPos,
        normFlapPos=normFlapPos,
        normFlapSpan=normFlapSpan,
        δ=δAil,
        updateAirfoilParameters=updateAirfoilParameters,
        hasInducedDrag=hasInducedDrag,
        hasSymmetricCounterpart=true
        )

    # Wing properties
    L = 16
    GJ,EIy,EIz = 1e4, 2e4, 4e6
    GJ *= torsionStiffnessFactor
    GJ_EIz = -Ψ * sqrt(GJ * EIz)
    GJ,EIy,EIz,GJ_EIz = multiply_scalars!(stiffnessFactor, GJ,EIy,EIz,GJ_EIz)
    ρA,ρIx = 0.75,0.1
    ρIy,ρIz = (EIy/EIz)*ρIx, (1-EIy/EIz)*ρIx

    # Sectional stiffness and inertia matrices
    S = anisotropic_stiffness_matrix(
        ∞=∞,
        GJ=GJ,
        EIy=EIy,
        EIz=EIz,
        c46=GJ_EIz,
        )
    I = inertia_matrix(
        ρA=ρA,
        ρIy=ρIy,
        ρIz=ρIz,
        )

    # Wing beam
    beam = create_Beam(
        name="rightWing",
        length=L,
        nElements=nElem,
        S=[S],
        I=[I],
        p0=transform_to_WM([0; 0; pitchAngle],from=:E321),
        k=[k1; k2; 0],
        surface=wingSurf,
        )
    
    ## Update beam of additional BCs, if applicable
    [reset_BC_beam!(BC,beam) for BC in additionalBCs]
    BCs = additionalBCs

    # Add root clamp and tip load to BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0, :u2A=>0, :u3A=>0, :p1A=>0, :p2A=>0, :p3A=>0)
        )
    tipLoad = create_BC(
        name="tipF3",
        beam=beam,
        node=nElem+1,
        conditions=Dict(:F3A=>tipF3)
        )
    push!(BCs,clamp,tipLoad)

    # Wing model
    wing = create_Model(
        name="SMW",
        beams=[beam],
        BCs=BCs,
        gravityVector=[0; 0; -g],
        v_A=t->[0; airspeed(t); 0],
        atmosphere=atmosphere,
        )

    return wing,L

end
export create_SMW


"""
    create_TDWing(; kwargs...)

Creates the wing described by Tang and Dowell in: Experimental and Theoretical Study on Aeroelastic Response of High-Aspect-Ratio Wings (2001)

Returns the wing model

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `pitchAngle::Real`: pitch angle
- `nElem::Int`: number of elements for discretization
- `altitude::Real`: altitude
- `airspeed::Union{Real,<:Function}`: airspeed
- `g::Real`: acceleration of gravity
- `∞::Real`: value of rigid structural properties
"""
function create_TDWing(; 
    aeroSolver::AeroSolver=Indicial(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(FlatPlate),
    updateAirfoilParameters::Bool=true,
    pitchAngle::Real=0,
    nElem::Int=20,
    altitude::Real=0,
    airspeed::Union{Real,<:Function}=0,
    g::Real=standard_atmosphere(altitude).g,
    ∞::Real=1e6
    )

    # Validate
    @assert -π/2 <= pitchAngle <= π/2 "set pitchAngle between -π/2 and π/2 (rad)"
    @assert nElem > 1 "nElem must be greater than 1"
    @assert altitude >= 0 "altitude must be non-negative"
    
    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,airspeed(0)/atmosphere.a)

    # Wing surface
    chord = 0.0508
    normSparPos = 0.5
    wingSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=airfoil,
        c=chord,
        normSparPos=normSparPos,
        updateAirfoilParameters=updateAirfoilParameters
        )

    # Wing properties
    L = 0.4508
    GJ,EIy,EIz = 0.9539, 0.4186, 18.44
    ρA,ρIx,ρIy = 0.2351, 0.2056e-4, 1e-6
    ρIz = ρIx - ρIy
    e3 = 1e-2 * chord

    # Sectional stiffness and inertia matrices
    S = isotropic_stiffness_matrix(
        ∞=∞,
        GJ=GJ,
        EIy=EIy,
        EIz=EIz
        )
    I = inertia_matrix(
        ρA=ρA,
        ρIy=ρIy,
        ρIz=ρIz,
        e3=e3
        )

    # Wing beam
    beam = create_Beam(
        name="wingBeam",
        length=L,
        nElements=nElem,
        S=[S],
        I=[I],
        p0=transform_to_WM([0; 0; pitchAngle], from=:E321),
        surface=wingSurf,
        )

    # Wing's tip store
    tipMass = 0.0417
    tipIyy = 0.3783e-5
    tipIzz = 0.9753e-4
    tipStore = create_PointInertia(
        mass=tipMass,
        η=[L/nElem/2;0;0],
        Iyy=tipIyy,
        Izz=tipIzz
        )
    add_point_inertia_to_beam!(beam=beam,elementID=nElem,pointInertia=tipStore)

    # BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0)
        )

    # Wing model
    wing = create_Model(
        name="TDWing",
        beams=[beam],
        BCs=[clamp],
        gravityVector=[0 ;0; -g],
        v_A=t->[0; airspeed(t); 0],
        atmosphere=atmosphere
        )

    return wing

end
export create_TDWing


"""
    create_Helios(; kwargs...)

Creates a model based on the flying-wing aircraft described by Patil and Hodges in: Flight Dynamics of Highly Flexible Flying Wings (2006)

# Keyword arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `flapSolver::FlapAeroSolver`: flap aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `altitude::Real`: altitude
- `g::Real`: local acceleration of gravity
- `wingAirfoil::Airfoil`: airfoil section of the wing
- `podAirfoil::Airfoil`: airfoil section of the pods
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `beamPods::Bool`: flag to include pods
- `stiffnessFactor::Real`: stiffness factor for the wing structure
- `∞::Real`: value for rigid structural properties
- `nElemStraightSemispan::Int`: number of elements in the straight section of the semispan
- `nElemDihedralSemispan::Int`: number of elements in the dihedral section of the semispan
- `nElemPod::Int`: number of elements in the pods
- `payloadPounds::Real`: payload, in pounds
- `airspeed::Real`: local initial/trim airspeed
- `δIsTrimVariable::Bool`: flag for flap deflection being a trim variable
- `thrustIsTrimVariable::Bool`: flag for motors' thrust being a trim variable
- `δ::Union{Nothing,Real,<:Function}`: flap deflection
- `thrust::Union{Real,<:Function}`: motors' thrust
- `reducedChord::Bool`: flag to employ a reduced (7 ft) chord
- `payloadOnWing::Bool`: flag to set the payload on the wing's reference line
- `wingRootAoA::Real`: root angle of attack for clamped wing model [rad]
- `withAttachmentSprings::Bool`: flag to include attachment springs for eigen-analyses
- `μ::Real`: stiffness of the attachment springs
- `gust::Union{Nothing,Gust}`: gust
"""
function create_Helios(; 
    aeroSolver::AeroSolver=Indicial(),
    flapSolver::FlapAeroSolver=TableLookup(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    altitude::Real=0,
    g::Real=standard_atmosphere(altitude).g,
    wingAirfoil::Airfoil=deepcopy(HeliosWingAirfoil),
    podAirfoil::Airfoil=deepcopy(HeliosPodAirfoil),
    updateAirfoilParameters::Bool=true,
    beamPods::Bool=false,
    stiffnessFactor::Real=1.0,
    ∞::Real=1e12,
    nElemStraightSemispan::Int=10,
    nElemDihedralSemispan::Int=5,
    nElemPod::Int=2,
    payloadPounds::Real=0,
    airspeed::Real=0,
    δIsTrimVariable::Bool=false,
    thrustIsTrimVariable::Bool=false,
    δ::Union{Nothing,Real,<:Function}=nothing,
    thrust::Union{Real,<:Function}=0,
    reducedChord::Bool=false,
    payloadOnWing::Bool=false,
    wingRootAoA::Real=0,
    withAttachmentSprings::Bool=false,
    μ::Real=1e-1,
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate
    @assert ∞ > 1e10 "rigid structural properties must be greater than 1e10"
    @assert stiffnessFactor > 0 "stiffnessFactor must be greater than 0"
    @assert payloadPounds >= 0 "payloadPounds must be non-negative"
    @assert airspeed >= 0 "airspeed must be non-negative"
    @assert μ >= 0 "spring stiffness μ must be non-negative"
    @assert nElemStraightSemispan > 0 "nElemStraightSemispan must be positive"
    @assert iseven(nElemStraightSemispan) "nElemStraightSemispan must be even"
    @assert nElemDihedralSemispan > 0 "nElemDihedralSemispan must be positive"
    @assert nElemPod > 0 "nElemPod must be positive"
    δIsInput = !isnothing(δ)
    if δIsInput
        @assert !δIsTrimVariable "cannot set δIsTrimVariable to true when δ is provided as input"
    end
    if δ isa Real
        δconst = δ
        δ = t -> δconst
    end
    if thrust isa Real
        thrustConst = thrust
        thrust = t -> thrustConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(wingAirfoil,airspeed/atmosphere.a)
    update_all_airfoil_parameters!(podAirfoil,airspeed/atmosphere.a)

    # Wing surface
    wChord = reducedChord ? 
             feet2meters(7) : 
             feet2meters(8)
    wNormSparPos = 0.25
    wNormFlapPos = 0.75
    wNormFlapSpan = [0,1]

    wingSurf = create_AeroSurface(
            ppSolver=aeroSolver,
            gustSolver=gustSolver,
            derivationMethod=derivationMethod,
            flapSolver=flapSolver,
            airfoil=wingAirfoil,
            c=wChord,
            normSparPos=wNormSparPos,
            normFlapPos=wNormFlapPos,
            normFlapSpan=wNormFlapSpan,
            δ = δIsInput ? δ : nothing,
            δIsTrimVariable=δIsTrimVariable,
            updateAirfoilParameters=updateAirfoilParameters
            ) 
    
    # Wing properties
    wingSection = feet2meters(40)
    Γ = deg2rad(10)
    EA,GJ,EIy,EIz = ∞, 1.6530e5, 1.0331e6, 1.2398e7
    wρA,wρIy,wρIz = 8.929, 0.691, 3.456
    EA,GJ,EIy,EIz = multiply_scalars!(stiffnessFactor, EA,GJ,EIy,EIz)

    # Wing sectional stiffness and inertia matrices
    Swing = isotropic_stiffness_matrix(
        ∞=∞,
        EA=EA,
        GJ=GJ,
        EIy=EIy,
        EIz=EIz
        )
    Iwing = inertia_matrix(
        ρA=wρA,
        ρIy=wρIy,
        ρIz=wρIz
        )

    # Wing beams
    leftWingDihedral = create_Beam(
        name="leftWingDihedral",
        length=wingSection,
        nElements=nElemDihedralSemispan,
        S=[Swing],
        I=[Iwing],
        surface=deepcopy(wingSurf),
        p0=transform_to_WM([0;Γ;0],from=:E321)
        )

    leftWingStraight = create_Beam(
        name="leftWingStraight",
        length=2*wingSection,
        nElements=nElemStraightSemispan,
        S=[Swing],
        I=[Iwing],
        surface=deepcopy(wingSurf)
        )

    rightWingStraight = create_Beam(
        name="rightWingStraight",
        length=2*wingSection,
        nElements=nElemStraightSemispan,
        S=[Swing],
        I=[Iwing],
        surface=deepcopy(wingSurf)
        )

    rightWingDihedral = create_Beam(
        name="rightWingDihedral",
        length=wingSection,
        nElements=nElemDihedralSemispan,
        S=[Swing],
        I=[Iwing],
        surface=deepcopy(wingSurf),
        p0=transform_to_WM([0;-Γ;0],from=:E321)
        )

    # Link elevators
    elevatorLink = create_FlapLink(
        masterBeam=rightWingStraight,
        slaveBeams=[leftWingDihedral,leftWingStraight,rightWingDihedral]
        )

    # Pods' point inertias
    podLength = feet2meters(6)
    podMassVertOffset = payloadOnWing ? 
                        0 : 
                        -podLength/2
    sidePodMass = pounds2kilograms(50)
    centerPodMass = pounds2kilograms(60)
    payloadMass = pounds2kilograms(payloadPounds)

    leftSidePod = create_PointInertia(
        η=[-(2*wingSection/nElemStraightSemispan)/2, 0, podMassVertOffset],
        mass=sidePodMass
        )

    centerPod = create_PointInertia(
        η=[(2*wingSection/nElemStraightSemispan)/2, 0, podMassVertOffset],
        mass=centerPodMass+payloadMass
        )

    rightSidePod = create_PointInertia(
        η=[(2*wingSection/nElemStraightSemispan)/2, 0, podMassVertOffset],
        mass=sidePodMass
        )

    add_point_inertia_to_beam!(beam=leftWingStraight,pointInertia=leftSidePod, elementID=1)
    add_point_inertia_to_beam!(beam=leftWingStraight,pointInertia=centerPod, elementID=nElemStraightSemispan)
    add_point_inertia_to_beam!(beam=rightWingStraight,pointInertia=rightSidePod,elementID=nElemStraightSemispan)

    # Pod surface
    podSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=podAirfoil,
        c=wChord,
        normSparPos=wNormSparPos,
        updateAirfoilParameters=false
        )

    # Pod properties
    pρA,pρIy,pρIz = 1e-3, wρIy, wρIz

    # Pod sectional stiffness and inertia matrices
    Spod = isotropic_stiffness_matrix(∞=∞)
    Ipod = inertia_matrix(
        ρA=pρA,
        ρIy=pρIy,
        ρIz=pρIz
        )

    # Pod beams
    leftPod = create_Beam(
        name="leftPod",
        length=podLength,
        nElements=nElemPod,
        S=[Spod],
        I=[Ipod],
        surface=deepcopy(podSurf),
        p0=transform_to_WM([0,π/2,0],from=:E321),
        )

    centerPod = create_Beam(
        name="centerPod",
        length=podLength,
        nElements=nElemPod,
        S=[Spod],
        I=[Ipod],
        surface=deepcopy(podSurf),
        p0=transform_to_WM([0,π/2,0],from=:E321),
        )

    rightPod = create_Beam(
        name="rightPod",
        length=podLength,
        nElements=nElemPod,
        S=[Spod],
        I=[Ipod],
        surface=deepcopy(podSurf),
        p0=transform_to_WM([0,π/2,0],from=:E321),
        )

    # Beam connections for aircraft model
    connections = beamPods ?
    [
        (leftPod => leftWingStraight, 1 => 1),
        (centerPod => rightWingStraight, 1 => 1),
        (rightPod => rightWingStraight, 1 => nElemStraightSemispan+1),
    ] : 
    Tuple{Pair{Beam,Beam},Pair{Int,Int}}[]

    # Attachment springs
    if withAttachmentSprings
        
        spring = create_Spring(ku=μ,kp=μ)

        orientations = ([0,0,0], [π/2,0,0], [0,π/2,0])

        for p in orientations
            add_spring_to_beam!(beam=rightWingStraight, spring=deepcopy(spring), nodeID=1, p=p)
        end
        
    end

    # Propellers thrust force
    thrustLeftOut = create_BC(
        name="thrustLeftOut",
        beam=leftWingStraight,
        node=1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    thrustLeftIn = create_BC(
        name="thrustLeftIn",
        beam=leftWingStraight,
        node=div(nElemStraightSemispan,2)+1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    thrustCenter = create_BC(
        name="thrustCenter",
        beam=leftWingStraight,
        node=nElemStraightSemispan+1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    thrustRightIn = create_BC(
        name="thrustRightIn",
        beam=rightWingStraight,
        node=div(nElemStraightSemispan,2)+1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    thrustRightOut = create_BC(
        name="thrustRightOut",
        beam=rightWingStraight,
        node=nElemStraightSemispan+1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    aircraftBCs = [thrustLeftOut,thrustLeftIn,thrustCenter,thrustRightIn,thrustRightOut]

    # Link propellers' thrust, if applicable
    thrustsLink = thrustIsTrimVariable ? 
        [create_TrimLoadsLink(
            masterBC=thrustCenter,
            slaveBCs=[thrustLeftOut,thrustLeftIn,thrustRightIn,thrustRightOut]
            )
        ] : 
        Vector{TrimLoadsLink}()

    # Beams of aircraft model
    aircraftBeams = beamPods ? 
        [leftWingDihedral,leftWingStraight,rightWingStraight,rightWingDihedral,leftPod,centerPod,rightPod] : 
        [leftWingDihedral,leftWingStraight,rightWingStraight,rightWingDihedral]

    # Position of the first node of the aircraft model (left wingtip)
    firstNodePosition = [-wingSection*(2+cos(Γ)); 0; wingSection*sin(Γ)]

    # Generate beam copies for wing model (has to be done before aircraft model creation)
    wingStraight = deepcopy(rightWingStraight)
    wingDihedral = deepcopy(rightWingDihedral)
    wingPod = deepcopy(rightPod)

    # Set root angle of attack for beams of wing model
    for beam in [wingStraight,wingDihedral,wingPod]
        update_beam_orientation!(beam=beam, p0=beam.geometry.p0+transform_to_WM([0, 0, wingRootAoA],from=:E321))
    end        

    # Aircraft model (with initial position such that aircraft center is coincident with the origin of frame A)
    helios = create_Model(
        name="Helios",
        beams=aircraftBeams,
        connections=connections,
        BCs=aircraftBCs,
        firstNodePosition=firstNodePosition,
        gravityVector=[0;0;-g],
        v_A=[0; airspeed; 0],
        flapLinks=[elevatorLink],
        trimLoadsLinks=thrustsLink,
        atmosphere=atmosphere,
        gust=gust
        )

    # Set midpsan element for the aircraft model
    midSpanElem = nElemDihedralSemispan + nElemStraightSemispan        

    ### --- Wing model --- ###
    # --------------------------------------------------------------------------
    # Beams of wing model
    wingBeams = beamPods ? 
        [wingStraight,wingDihedral,wingPod] : 
        [wingStraight,wingDihedral]

    # Clamp for wing model
    clamp = create_BC(
        name="clamp",
        beam=wingStraight,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0),
        )

    # Link elevators
    elevatorLinkWing = create_FlapLink(
        masterBeam=wingStraight,
        slaveBeams=[wingDihedral]
        )

    # Beam connections for wing model
    wingConnections = beamPods ?
    [
        (wingPod => wingStraight, 1 => nElemStraightSemispan+1),
    ] : 
    Tuple{Pair{Beam,Beam},Pair{Int,Int}}[]        

    # Wing model
    wing = create_Model(
        name="HeliosWing",
        beams=wingBeams,
        BCs=[clamp],
        connections=wingConnections,
        gravityVector=[0; 0; -g],
        v_A=[0; airspeed; 0],
        flapLinks=[elevatorLinkWing],
        atmosphere=atmosphere,
        gust=gust
        )

    return helios,midSpanElem,wing,leftWingStraight,rightWingStraight,leftWingDihedral,rightWingDihedral,leftPod,rightPod,centerPod

end
export create_Helios


"""
    create_cHALE(; kwargs...)

Creates a model based on the conventional HALE aircraft described by Patil, Hodges and Cesnik in: Nonlinear Aeroelasticity and Flight Dynamics of HALE (2001)

# Keyword arguments
- `altitude::Real`: altitude
- `aeroSolver::AeroSolver`: aerodynamic solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `flapSolver::FlapAeroSolver`: aerodynamic solver for flap loads
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `stabilizersAero::Bool`: flag for stabilizers with aerodynamic surfaces
- `includeVS::Bool`: flag to include a vertical stabilizer in the model
- `wAirfoil::Airfoil`: wing airfoil section
- `sAirfoil::Airfoil`: stabilizers airfoil section
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `nElemWing::Int`: number of elements of the full wing
- `nElemHorzStabilizer::Int`: number of elements of the horizontal stabilizer
- `nElemTailBoom::Int`: number of elements of the tail boom
- `nElemVertStabilizer::Int`: number of elements of the vertical stabilizer
- `∞::Real=1e12`: value of rigid structural properties
- `stiffnessFactor::Real`: stiffness factor for the wing structure
- `k1::Real`: undeformed wing torsional curvature
- `k2::Real`: undeformed wing flapwise bending curvature
- `airspeed::Real`: local initial/trim airspeed
- `δElevIsTrimVariable::Bool`: flag for elevator deflection being a trim variable
- `δAilIsTrimVariable::Bool`: flag for aileron deflection being a trim variable
- `δRuddIsTrimVariable::Bool`: flag for rudder deflection being a trim variable    
- `thrustIsTrimVariable::Bool`: flag for motors' thrust being a trim variable
- `δElev::Union{Nothing,Real,<:Function}`: elevator deflection
- `δAil::Union{Nothing,Real,<:Function}`: aileron deflection
- `δRudd::Union{Nothing,Real,<:Function}`: rudder deflection
- `thrust::Union{Real,<:Function}`: motors' thrust
- `g::Real`: local acceleration of gravity
- `wingCd0::Real`: parasite drag coefficient for the wing
- `wingcnδ::Real`: cn vs δ slope for the wing
- `wingcmδ::Real`: cm vs δ slope for the wing
- `wingcdδ::Real`: cd vs δ slope for the wing
- `stabsCd0::Real`: parasite drag coefficient for the stabilizers
- `stabscnδ::Real`: cn vs δ slope for the stabilizers
- `stabscmδ::Real`: cm vs δ slope for the stabilizers
- `stabscdδ::Real`: cd vs δ slope for the stabilizers
- `hasInducedDrag::Bool`: flag to include induced drag
- `rootClamped::Bool`: flag for being clamped at wing root
- `withAttachmentSprings::Bool`: flag to include attachment springs for eigen-analyses
- `μ::Real`: stiffness coefficient for the attachment springs
- `gust::Union{Nothing,Gust}`: gust input for the model
"""
function create_cHALE(; 
    altitude::Real=20e3,
    aeroSolver::AeroSolver=Indicial(),
    derivationMethod::DerivationMethod=AD(),
    flapSolver::FlapAeroSolver=TableLookup(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    stabilizersAero::Bool=true,
    includeVS::Bool=true,
    wAirfoil::Airfoil=deepcopy(cHALEairfoil),
    sAirfoil::Airfoil=deepcopy(cHALEairfoil),
    updateAirfoilParameters::Bool=false,
    nElemWing::Int=20,
    nElemHorzStabilizer::Int=10,
    nElemTailBoom::Int=10,
    nElemVertStabilizer::Int=5,
    ∞::Real=1e12,
    stiffnessFactor::Real=1,
    torsionStiffnessFactor::Real=1,
    k1::Real=0,
    k2::Real=0,
    airspeed::Real=0,
    δElevIsTrimVariable::Bool=false,
    δAilIsTrimVariable::Bool=false,
    δRuddIsTrimVariable::Bool=false,
    thrustIsTrimVariable::Bool=false,
    δElev::Union{Nothing,Real,<:Function}=nothing,
    δAil::Union{Nothing,Real,<:Function}=nothing,
    δRudd::Union{Nothing,Real,<:Function}=nothing,
    thrust::Union{Real,<:Function}=0,
    g::Real=standard_atmosphere(altitude).g,
    wingCd0::Real=0,
    wingcnδ::Real=2.5,
    wingcmδ::Real=-0.35,
    wingcdδ::Real=0.15,
    stabsCd0::Real=0,
    stabscnδ::Real=2.5,
    stabscmδ::Real=-0.35,
    stabscdδ::Real=0.15,
    hasInducedDrag::Bool=false,
    rootClamped::Bool=false,
    withAttachmentSprings::Bool=false,
    μ::Real=1e-1,
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate
    @assert iseven(nElemWing) "nElemWing must be even"
    @assert iseven(nElemHorzStabilizer) "nElemHorzStabilizer must be even"
    @assert ∞ > 1e8 "rigid structural properties must be greater than 1e8"
    @assert stiffnessFactor > 0 "stiffnessFactor must be greater than 0"
    @assert wingCd0 >= 0 "wingCd0 must be non-negative"
    @assert stabsCd0 >= 0 "stabsCd0 must be non-negative"
    @assert airspeed >= 0 "airspeed must be non-negative"
    @assert g >= 0 "g must be non-negative"
    @assert μ >= 0 "spring stiffness μ must be non-negative"
    δElevIsInput = !isnothing(δElev)
    δAilIsInput = !isnothing(δAil)
    δRuddIsInput = !isnothing(δRudd)
    if δElevIsInput
        @assert !δElevIsTrimVariable "cannot set δElevIsTrimVariable to true when δElev is provided as input"
    end
    if δElev isa Real
        δElevConst = δElev
        δElev = t -> δElevConst
    end
    if δAilIsInput
        @assert !δAilIsTrimVariable "cannot set δAilIsTrimVariable to true when δAil is provided as input"
    end
    if δAil isa Real
        δAilConst = δAil
        δAil = t -> δAilConst
    end
    if δRuddIsInput
        @assert !δRuddIsTrimVariable "cannot set δRuddIsTrimVariable to true when δRudd is provided as input"
    end
    if δRudd isa Real
        δRuddConst = δRudd
        δRudd = t -> δRuddConst
    end
    if thrust isa Real
        thrustConst = thrust
        thrust = t -> thrustConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Set modified airfoils
    wFlapSiteID = 75
    dbWingAirfoil = modified_database(
        wAirfoil.name;
        attached = (
            cd₀ = wingCd0,
        ),
        flap = (
            wFlapSiteID => create_FlapParameters(wingcdδ,wingcmδ,wingcnδ),
        )
    )
    wAirfoil = create_Airfoil(database=dbWingAirfoil,Ma=airspeed/atmosphere.a)

    sFlapSiteID = 75
    dbStabAirfoil = modified_database(
        sAirfoil.name;
        attached = (
            cd₀ = stabsCd0,
        ),
        flap = (
            sFlapSiteID => create_FlapParameters(stabscdδ,stabscmδ,stabscnδ),
        )
    )
    sAirfoil = create_Airfoil(database=dbStabAirfoil,Ma=airspeed/atmosphere.a)

    # Wing surface
    wChord = 1
    wNormSparPos = 0.5
    wNormFlapPos = wFlapSiteID/100
    wFlapSize = 0.25
    wLNormFlapSpan = [0, wFlapSize]
    wRNormFlapSpan = [1-wFlapSize, 1]

    wingSurfLeft = create_AeroSurface(
            ppSolver=aeroSolver,
            gustSolver=gustSolver,
            derivationMethod=derivationMethod,
            flapSolver=flapSolver,
            airfoil=wAirfoil,
            c=wChord,
            normSparPos=wNormSparPos,
            normFlapPos=wNormFlapPos,
            normFlapSpan=wLNormFlapSpan,
            δ = δAilIsInput ? δAil : nothing,
            δIsTrimVariable=δAilIsTrimVariable,
            updateAirfoilParameters=updateAirfoilParameters,
            hasInducedDrag=hasInducedDrag,
            hasSymmetricCounterpart=true
            )

    wingSurfRight = create_AeroSurface(
            ppSolver=aeroSolver,
            gustSolver=gustSolver,
            derivationMethod=derivationMethod,
            flapSolver=flapSolver,
            airfoil=wAirfoil,
            c=wChord,
            normSparPos=wNormSparPos,
            normFlapPos=wNormFlapPos,
            normFlapSpan=wRNormFlapSpan,
            δ = δAilIsInput ? δAil : nothing,
            δIsTrimVariable=δAilIsTrimVariable,
            updateAirfoilParameters=updateAirfoilParameters,
            hasInducedDrag=hasInducedDrag,
            hasSymmetricCounterpart=true
            )

    # Wing properties
    Lw = 16
    wGJ,wEIy,wEIz = 1e4, 2e4, 4e6
    wGJ,wEIy,wEIz = multiply_scalars!(stiffnessFactor, wGJ,wEIy,wEIz)
    wGJ *= torsionStiffnessFactor
    wρA,wρIx = 0.75,0.1
    wρIy,wρIz = (wEIy/wEIz)*wρIx, (1-wEIy/wEIz)*wρIx

    # Wing sectional stiffness and inertia matrices
    Swing = isotropic_stiffness_matrix(
        ∞=∞,
        GJ=wGJ,
        EIy=wEIy,
        EIz=wEIz
        )
    Iwing = inertia_matrix(
        ρA=wρA,
        ρIy=wρIy,
        ρIz=wρIz
        )

    # Position of left wingtip node (first node of the assembly)
    ρ = 1/k2
    θ = Lw/ρ
    x = ρ*sin(θ)
    z = ρ*(1-cos(θ))
    firstNodePosition = iszero(k2) ? 
                        [-Lw; 0; 0] : 
                        [-x; 0; -z]

    # Initial angle of twist
    r = 1/k1
    ψ = Lw/r

    # Wing beams
    leftWing = create_Beam(
        name="leftWing",
        length=Lw,
        nElements=div(nElemWing,2),
        S=[Swing],
        I=[Iwing],
        k=[-k1; k2; 0],
        p0=transform_to_WM([0; -θ; ψ], from=:E321),
        surface=wingSurfLeft,
        )

    rightWing = create_Beam(
        name="rightWing",
        length=Lw,
        nElements=div(nElemWing,2),
        S=[Swing],
        I=[Iwing],
        k=[k1; k2; 0],
        surface=wingSurfRight,
        )

    # Link wing ailerons
    aileronLink = create_FlapLink(
        masterBeam=rightWing,
        slaveBeams=[leftWing],
        δMultipliers=[-1]
        )

    # Payload (@ wing root)
    payloadMass = 50
    payloadInertia = 200
    payload = create_PointInertia(
        mass=payloadMass,
        η=[-Lw/div(nElemWing,2)/2; 0; 0],
        Iyy=payloadInertia,
        Izz=payloadInertia,
        Ixx=payloadInertia
        )
    add_point_inertia_to_beam!(beam=rightWing,elementID=1,pointInertia=payload)

    # Tailboom properties
    Lt = 10
    tρA,tρIy,tρIz = 0.08, wρIy/10, wρIz/10
    Stail = isotropic_stiffness_matrix(∞=∞)
    Itail = inertia_matrix(ρA=tρA,ρIy=tρIy,ρIz=tρIz)

    # Tailboom beam
    tailBoom = create_Beam(
        name="tailBoom",
        length=Lt,
        nElements=nElemTailBoom,
        S=[Stail],
        I=[Itail],
        p0=transform_to_WM([-π/2; 0; 0], from=:E321),
        )
    
    # Horizontal stabilizer surface
    hChord = 0.5
    hNormSparPos = 0.5
    hNormFlapPos = sFlapSiteID/100
    hNormFlapSpan = [0, 1]

    hsSurf = stabilizersAero ? 
             create_AeroSurface(
                ppSolver=aeroSolver,
                gustSolver=gustSolver,
                derivationMethod=derivationMethod,
                flapSolver=flapSolver,
                airfoil=sAirfoil,
                c=hChord,
                normSparPos=hNormSparPos,
                normFlapPos=hNormFlapPos,
                normFlapSpan=hNormFlapSpan,
                δ = δElevIsInput ? δElev : nothing,
                δIsTrimVariable=δElevIsTrimVariable,
                updateAirfoilParameters=updateAirfoilParameters,
                hasInducedDrag=hasInducedDrag,
                hasSymmetricCounterpart=false
            ) : 
            nothing

    # Horizontal stabilizer properties
    Lh = 5
    hρA,hρIy,hρIz = 0.08, wρIy/10, wρIz/10
    Shorz = isotropic_stiffness_matrix(∞=∞)
    Ihorz = inertia_matrix(ρA=hρA,ρIy=hρIy,ρIz=hρIz)

    # Horizontal stabilizer beam
    horzStabilizer = create_Beam(
        name="horzStabilizer",
        length=Lh,
        firstNodeOffset=[-Lh/2; 0; 0],
        nElements=nElemHorzStabilizer,
        S=[Shorz],
        I=[Ihorz],
        surface=hsSurf
        )

    # Vertical stabilizer surface
    vChord = 0.5
    vNormSparPos = 0.5
    vNormFlapPos = 0.75
    vNormFlapSpan = [0, 1]

    vsSurf = stabilizersAero ? 
             create_AeroSurface(
                ppSolver=aeroSolver,
                gustSolver=gustSolver,
                derivationMethod=derivationMethod,
                flapSolver=flapSolver,
                airfoil=sAirfoil,
                c=vChord,
                normSparPos=vNormSparPos,
                normFlapPos=vNormFlapPos,
                normFlapSpan=vNormFlapSpan,
                δ = δRuddIsInput ? δRudd : nothing,
                δIsTrimVariable=δRuddIsTrimVariable,
                updateAirfoilParameters=updateAirfoilParameters,
                hasInducedDrag=hasInducedDrag,
                hasSymmetricCounterpart=false
            ) : 
            nothing

    # Vertical stabilizer properties
    Lv = 2.5
    vρA,vρIy,vρIz = 0.08, hρIy, hρIz
    Svert = isotropic_stiffness_matrix(∞=∞)
    Ivert = inertia_matrix(ρA=vρA,ρIy=vρIy,ρIz=vρIz)

    # Vertical stabilizer beam
    vertStabilizer = create_Beam(
        name="vertStabilizer",
        length=Lv,
        nElements=nElemVertStabilizer,
        S=[Svert],
        I=[Ivert],
        p0=transform_to_WM([0; -π/2; 0], from=:E321),
        surface=vsSurf
        )

    # Attachment springs
    if withAttachmentSprings
        
        spring = create_Spring(ku=μ,kp=μ)

        orientations = ([0,0,0], [π/2,0,0], [0,π/2,0])

        for nodeID in (1, nElemTailBoom+1), p in orientations
            add_spring_to_beam!(beam=tailBoom, spring=deepcopy(spring), nodeID=nodeID, p=p)
        end
        
    end

    # Propeller thrust force
    propThrust = create_BC(
        name="propThrust",
        beam=rightWing,
        node=1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    # Clamp wing root (if applicable)
    clamp = create_BC(
        name="clamp",
        beam=rightWing,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0)
        )

    # Beams of the model
    beams = includeVS ? 
            [leftWing,rightWing,tailBoom,horzStabilizer,vertStabilizer] : 
            [leftWing,rightWing,tailBoom,horzStabilizer]

    # Beam connections
    connections = includeVS ?
    [
        (rightWing => leftWing, 1 => div(nElemWing,2)+1),
        (tailBoom => rightWing, 1 => 1),
        (horzStabilizer => tailBoom, div(nElemHorzStabilizer,2)+1 => nElemTailBoom+1),
        (vertStabilizer => tailBoom, 1 => nElemTailBoom+1)
    ] : 
    [
        (rightWing => leftWing, 1 => div(nElemWing,2)+1),
        (tailBoom => rightWing, 1 => 1),
        (horzStabilizer => tailBoom, div(nElemHorzStabilizer,2)+1 => nElemTailBoom+1),
    ]

    # BCs of the model
    BCs = rootClamped ? 
          [clamp,propThrust] : 
          [propThrust]

    # Aircraft model
    cHALE = create_Model(
        name="cHALE",
        beams=beams,
        BCs=BCs,
        firstNodePosition=firstNodePosition,
        gravityVector=[0; 0; -g],
        connections=connections,
        v_A=[0; airspeed; 0],
        flapLinks=[aileronLink],
        atmosphere=atmosphere,
        gust=gust,
        )

    return cHALE,leftWing,rightWing,tailBoom,horzStabilizer,vertStabilizer

end
export create_cHALE


"""
    create_BWB(; kwargs...)

Creates a model based on the blended-wing-body described by Weihua Su's PhD thesis and in https://arc.aiaa.org/doi/10.2514/1.47317

# Keyword arguments
- `altitude::Real`: altitude
- `aeroSolver::AeroSolver`: aerodynamic solver
- `flapSolver::FlapAeroSolver`: flap-induced loads solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil of the wing
- `∞::Real=1e12`: value of rigid structural properties
- `stiffnessFactor::Real`: stiffness factor for the wing structure
- airspeed::Real = local initial/trim airspeed
- `δElevIsTrimVariable::Bool`: flag for elevator deflection being a trim variable
- `thrustIsTrimVariable::Bool`: flag for motors' thrust being a trim variable
- `δElev::Union{Nothing,Real,<:Function}`: elevator deflection
- `thrust::Union{Real,<:Function}`: motors' thrust
- `g::Real`: local acceleration of gravity
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `hasTipCorrection::Bool`: flag to employ aerodynamic tip correction
- `tipLossDecayFactor::Real`: tip loss decay factor
- `withAttachmentSprings::Bool`: flag to include attachment springs for eigen-analyses
- `μ::Real`: attachment spring constant (effective only if `withAttachmentSprings` is true)
- `wingRootAoA::Real`: wing root pitch angle for wing model [rad]
"""
function create_BWB(; 
    altitude::Real=0,
    aeroSolver::AeroSolver=Indicial(),
    flapSolver::FlapAeroSolver=TableLookup(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(BWBAirfoil),
    ∞::Real=1e12,
    stiffnessFactor::Real=1,
    airspeed::Real=0,
    δElevIsTrimVariable::Bool=false,
    thrustIsTrimVariable::Bool=false,
    δElev::Union{Nothing,Real,<:Function}=nothing,
    thrust::Union{Real,<:Function}=0,
    g::Real=standard_atmosphere(altitude).g,
    updateAirfoilParameters::Bool=true,
    hasTipCorrection::Bool=false,
    tipLossDecayFactor::Real=40,
    withAttachmentSprings::Bool=false,
    μ::Real=1e-1,
    wingRootAoA::Real=0,
    )

    # Validate
    @assert ∞ > 1e10 "rigid structural properties must be greater than 1e10"
    @assert stiffnessFactor > 0 "stiffnessFactor must be greater than 0"
    @assert airspeed >= 0 "airspeed must be non-negative"
    @assert μ >= 0 "attachment spring constant spring stiffness μ must be non-negative"
    δElevIsInput = !isnothing(δElev)
    if δElevIsInput
        @assert !δElevIsTrimVariable "δElevIsTrimVariable cannot be true when δElev is provided as input"
    end
    if δElev isa Real
        δElevConst = δElev
        δElev = t -> δElevConst
    end
    if thrust isa Real
        thrustConst = thrust
        thrust = t -> thrustConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Keypoints
    kp1 = [0; -1.126236+1.126236; 0]
    kp2 = [0.889; -0.997712+1.126236; 0]
    kp3 = [3.248152; -2.359914+1.126236; 0]

    # Length of fuselage section
    fusLength = sqrt( (kp2[1]-kp1[1])^2 + (kp2[2]-kp1[2])^2 )

    # Angle of sweep of the fuselage section 
    fΛ = -acos( (kp2[1]-kp1[1]) / fusLength )             

    # Length of wing section
    wingSemispan = sqrt( (kp3[1]-kp2[1])^2 + (kp3[2]-kp2[2])^2 )  
    
    # Angle of sweep of the wing section 
    wΛ = acos( (kp3[1]-kp2[1]) / wingSemispan )
    
    # Chord and spar position
    rootChord = 1.38557
    tipChord = 0.54864
    rootSparPos = 0.6438
    tipSparPos = 0.4560

    # Airfoil
    update_all_airfoil_parameters!(airfoil,airspeed/atmosphere.a)

    # Wing surfaces
    wChord = tipChord
    wNormSparPos = tipSparPos
    wNormFlapPos = 0.75
    normFlapSpanL = [1/4; 1]
    normFlapSpanR = [0; 3/4]

    leftWingSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        flapSolver=flapSolver,
        airfoil=airfoil,
        Λ=-wΛ,
        c=wChord,
        normSparPos=wNormSparPos,
        normFlapPos=wNormFlapPos,
        normFlapSpan=normFlapSpanL,
        δ = δElevIsInput ? δElev : nothing,
        δIsTrimVariable=δElevIsTrimVariable,
        updateAirfoilParameters=updateAirfoilParameters,
        hasTipCorrection=hasTipCorrection,
        tipLossDecayFactor=-tipLossDecayFactor
        )

    rightWingSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        flapSolver=flapSolver,
        airfoil=airfoil,
        Λ=wΛ,
        c=wChord,
        normSparPos=wNormSparPos,
        normFlapPos=wNormFlapPos,
        normFlapSpan=normFlapSpanR,
        δ = δElevIsInput ? δElev : nothing,
        δIsTrimVariable=δElevIsTrimVariable,
        updateAirfoilParameters=updateAirfoilParameters,
        hasTipCorrection=hasTipCorrection,
        tipLossDecayFactor=tipLossDecayFactor
        )

    # Wing properties
    nElemWing = 8
    wEA,wGJ,wEIy,wEIz = 155_000_000, 11_000, 11_700, 130_000
    wGJ,wEIy,wEIz = multiply_scalars!(stiffnessFactor, wGJ,wEIy,wEIz)
    wρA,wρIy,wρIz = 6.2, 5e-4, 4.62e-3

    # Wing sectional stiffness and inertia matrices
    Swing = isotropic_stiffness_matrix(
        ∞=∞,
        EA=wEA,
        GJ=wGJ,
        EIy=wEIy,
        EIz=wEIz
        )
    Iwing = inertia_matrix(
        ρA=wρA,
        ρIy=wρIy,
        ρIz=wρIz
        )

    # Wing beams
    leftWing = create_Beam(
        name="leftWing",
        length=wingSemispan,
        nElements=nElemWing,
        S=[Swing],
        I=[Iwing],
        surface=leftWingSurf,
        p0=transform_to_WM([wΛ;0;0],from=:E321),
        )

    rightWing = create_Beam(
        name="rightWing",
        length=wingSemispan,
        nElements=nElemWing,
        S=[Swing],
        I=[Iwing],
        surface=rightWingSurf,
        p0=transform_to_WM([-wΛ;0;0],from=:E321),
        )    

    # Concentrated wing inertias
    wConMass = 2.0
    wLTipConInertia = create_PointInertia(mass=wConMass,η=[-wingSemispan/nElemWing/2;0;0])
    wRTipConInertia = create_PointInertia(mass=wConMass,η=[wingSemispan/nElemWing/2;0;0])
    wLConInertia = create_PointInertia(mass=wConMass,η=[wingSemispan/nElemWing/2;0;0])
    wRConInertia = create_PointInertia(mass=wConMass,η=[-wingSemispan/nElemWing/2;0;0])

    # Add concentrated wing inertias to beams
    add_point_inertia_to_beam!(beam=leftWing,elementID=1,pointInertia=wLTipConInertia)
    add_point_inertia_to_beam!(beam=rightWing,elementID=nElemWing,pointInertia=wRTipConInertia)
    for e=1:nElemWing
        add_point_inertia_to_beam!(beam=leftWing,elementID=e,pointInertia=wLConInertia)
        add_point_inertia_to_beam!(beam=rightWing,elementID=e,pointInertia=wRConInertia)
    end

    # Link wing elevons
    elevonLink = create_FlapLink(
        masterBeam=rightWing,
        slaveBeams=[leftWing],
        δMultipliers=[1]
        )

    # Fuselage surfaces
    leftFusChord = x1 -> tipChord + (rootChord-tipChord) * x1 / fusLength
    rightFusChord = x1 -> rootChord + (tipChord-rootChord) * x1 / fusLength
    leftFusSparPos = x1 -> tipSparPos + (rootSparPos-tipSparPos) * x1 / fusLength
    rightFusSparPos = x1 -> rootSparPos + (tipSparPos-rootSparPos) * x1 / fusLength

    leftFusSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        derivationMethod=derivationMethod,
        flapSolver=flapSolver,
        airfoil=airfoil,
        Λ=-fΛ,
        c=leftFusChord,
        normSparPos=leftFusSparPos,
        updateAirfoilParameters=updateAirfoilParameters
        )

    rightFusSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        derivationMethod=derivationMethod,
        flapSolver=flapSolver,
        airfoil=airfoil,
        Λ=fΛ,
        c=rightFusChord,
        normSparPos=rightFusSparPos,
        updateAirfoilParameters=updateAirfoilParameters
        )

    # Fuselage properties
    nElemFus = 3
    fEA,fGJ,fEIy,fEIz = 169_000_000, 2_250_000, 750_000, 35_000_000
    fρA,fρIy,fρIz = 50, 0.7, 22.0

    # Fuselage sectional stiffness and inertia matrices
    Sfus = isotropic_stiffness_matrix(
        ∞=∞,
        EA=fEA,
        GJ=fGJ,
        EIy=fEIy,
        EIz=fEIz
        )
    Ifus = inertia_matrix(
        ρA=fρA,
        ρIy=fρIy,
        ρIz=fρIz
        )

    # Fuselage beams
    leftFus = create_Beam(
        name="leftFus",
        length=fusLength,
        nElements=nElemFus,
        S=[Sfus],
        I=[Ifus],
        surface=leftFusSurf,
        p0=transform_to_WM([fΛ;0;0],from=:E321),
        )

    rightFus = create_Beam(
        name="rightFus",
        length=fusLength,
        nElements=nElemFus,
        S=[Sfus],
        I=[Ifus],
        surface=rightFusSurf,
        p0=transform_to_WM([-fΛ;0;0],from=:E321),
        )        

    # Concentrated fuselage inertias
    fConMass = 40
    fConMassOffset = 0.891968
    fLeftConInertia = create_PointInertia(
        mass=fConMass,
        η=[fusLength/nElemFus/2; fConMassOffset; 0]
        )
    fRightConInertia = create_PointInertia(
        mass=fConMass,
        η=[-fusLength/nElemFus/2; fConMassOffset; 0]
        )

    # Add concentrated fuselage inertias to the beams
    add_point_inertia_to_beam!(beam=leftFus,elementID=nElemFus,pointInertia=fLeftConInertia)
    add_point_inertia_to_beam!(beam=rightFus,elementID=1,pointInertia=fRightConInertia)

    # Propellers thrust force
    thrustLeft = create_BC(
        name="thrustLeft",
        beam=leftFus,
        node=1,
        conditions=Dict(:Ff2A=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    thrustRight = create_BC(
        name="thrustRight",
        beam=rightFus,
        node=nElemFus+1,
        conditions=Dict(:Ff2A=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    # Link propellers' thrust, if applicable
    thrustsLink = thrustIsTrimVariable ? 
                  [create_TrimLoadsLink(
                   masterBC=thrustRight,
                   slaveBCs=[thrustLeft]
                   )] :
                  Vector{TrimLoadsLink}()

    # Attachment springs
    if withAttachmentSprings
        
        spring = create_Spring(ku=μ,kp=μ)

        orientations = ([0,0,0], [π/2,0,0], [0,π/2,0])

        for (beam, nodeID) in ((leftFus, 1), (rightFus, nElemFus+1)), p in orientations
            add_spring_to_beam!(beam=beam, spring=deepcopy(spring), nodeID=nodeID, p=p)
        end

    end

    # Initial position for left wingtip
    firstNodePosition = [-kp3[1]; kp3[2]; kp3[3]]

    # Generate beam copies for wing model (has to be done before aircraft model creation)
    fus = deepcopy(rightFus)
    wing = deepcopy(rightWing)

    # Set root angle of attack for beams of wing model
    for beam in [fus,wing]
        update_beam_orientation!(beam=beam, p0=beam.geometry.p0+transform_to_WM([0, 0, wingRootAoA],from=:E321))
    end

    # BCs for wing model (clamp + thrust)
    clamp = create_BC(
        name="clamp",
        beam=fus,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0),
        )
    wingThrust = create_BC(
        name="thrust",
        beam=fus,
        node=nElemFus+1,
        conditions=Dict(:Ff2A=>thrust)
        )

    # Aircraft model (with initial position such that aircraft center is coincident with the origin of frame A)
    BWB = create_Model(
        name="BWB",
        beams=[leftWing,leftFus,rightFus,rightWing],
        BCs=[thrustLeft,thrustRight],
        firstNodePosition=firstNodePosition,
        v_A=[0; airspeed; 0],
        gravityVector=[0; 0; -g],
        flapLinks=[elevonLink],
        trimLoadsLinks=thrustsLink,
        atmosphere=atmosphere,
        )

    # Right fuselage + wing model
    rightWingModel = create_Model(
        name="rightWingModel",
        beams=[fus,wing],
        BCs=[clamp,wingThrust],
        v_A=[0; airspeed; 0],
        gravityVector=[0; 0; -g],
        atmosphere=atmosphere,
        )

    return BWB,rightWingModel

end
export create_BWB


"""
    create_HealyBaselineFFWT(; kwargs...)

Creates a version of Healy's baseline wing with flared folding wingtip (FFWT). See Healy's PhD thesis, chapter 3.

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `hasTipCorrection::Bool`: flag for aerodynamic tip correction
- `tipLossDecayFactor::Real`: tip loss decay factor, in case hasTipCorrection is true
- `hingeConfiguration::Symbol`: hinge hingeConfiguration (:free or :locked)
- `flareAngle::Real`: flare angle [rad]
- `pitchAngle::Real`: root pitch angle [rad]
- `foldAngle::Union{Real,Nothing}`: fold angle [rad]
- `kSpring::Real`: stiffness of the rotational spring around the hinge
- `altitude::Real`: altitude
- `airspeed::Union{Real,<:Function}`: local airspeed
- `g::Real`: local acceleration of gravity
- `nElementsInner::Int`: number of elements for discretization of inner part of the wing
- `nElementsFFWT::Int`: number of elements for discretization of the wingtip
- `gust::Union{Nothing,Gust}=nothing`: gust
"""
function create_HealyBaselineFFWT(;
    aeroSolver::AeroSolver=Indicial(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(FlatPlate),
    updateAirfoilParameters::Bool=true,
    hasTipCorrection::Bool=false,
    tipLossDecayFactor::Real=12,
    hingeConfiguration::Symbol,
    flareAngle::Real=deg2rad(15),
    pitchAngle::Real=0,
    foldAngle::Union{Real,Nothing}=nothing,
    kSpring::Real=1.6e-4,
    altitude::Real=0,
    airspeed::Union{Real,<:Function}=0,
    g::Real=standard_atmosphere(altitude).g,
    nElementsInner::Int=16,
    nElementsFFWT::Int=4,
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate
    @assert hingeConfiguration in [:free,:locked] " set 'hingeConfiguration' as :free or :locked'"
    if !isnothing(foldAngle)
        @assert -π < foldAngle <= π "set foldAngle between -π and π (rad) "
    end
    @assert altitude >= 0 "altitude must be non-negative"
    @assert tipLossDecayFactor >= 0 "tipLossDecayFactor must be non-negative"
    @assert g >= 0 "g must be non-negative"
    @assert kSpring >= 0 "kSpring must be non-negative"
    @assert nElementsInner >= 8 "set at least 8 elements for inner wing"
    @assert nElementsFFWT >= 2 "set at least 2 elements for wingtip"

    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,airspeed(0)/atmosphere.a)

    # Wing (semi-)span, chord and normalized spar position
    L = 1
    chord = 0.12
    normSparPos = 0.25

    # Hinge normalized position over (semi-)span at the beam reference line (See Fig. 3.1 of Healy's thesis)
    hingePos = 0.8 + 22e-3 * tan(flareAngle) / L

    # Normalized nodal positions
    normalizedNodalPositions = unique(
        vcat(
             LinRange(0, hingePos, nElementsInner+1),
             LinRange(hingePos, 1, nElementsFFWT+1)
             )
        )
    
    # Total number of elements
    nElem = nElementsInner + nElementsFFWT    

    # Spar geometric and material properties
    b = 30e-3
    h = 4e-3
    E = 193e9
    G = 76e9
    ρ = 8e3
    A = b*h
    Iy = b*h^3/12
    Iz = h*b^3/12
    J = Iy+Iz

    # Stiffness correction factors (assumed)
    KOOP = 1.115
    KIP = 1.18
    KT = 0.070

    # Central spar stiffness matrix
    S1 = isotropic_stiffness_matrix(
        EA=E*A,
        GJ=G*J*KT,
        EIy=E*Iy*KOOP,
        EIz=E*Iz*KIP,
        ∞=1e8
        )

    # Wingtip stiffness matrix (assume 50% of central spar's stiffness)
    S2 = 0.5*S1

    # Central spar inertia matrix
    I1 = inertia_matrix(
        ρA=ρ*A,
        ρIy=ρ*Iy,
        ρIz=ρ*Iz
        )

    # Wingtip baseline inertia matrix (assume 10% of central spar's mass per unit length)
    I2 = inertia_matrix(
        ρA=0.1*ρ*A
        )

    # Aerodynamic surface
    surf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=airfoil,
        c=chord,
        normSparPos=normSparPos,
        hasTipCorrection=hasTipCorrection,
        tipLossDecayFactor=tipLossDecayFactor,
        updateAirfoilParameters=updateAirfoilParameters
        )

    # Hinge node
    hingeNode = nElementsInner + 1        

    # Arrays of hinged node and DOFs
    hingedNodes = hingeConfiguration == :free ? 
                  [hingeNode] : 
                  Vector{Int}()
    hingedNodesDoF = hingeConfiguration == :free ? 
                     [trues(3)] : 
                     [falses(3)]

    # Beam
    beam = create_Beam(
        name="beam",
        length=L,
        nElements=nElem,
        normalizedNodalPositions=normalizedNodalPositions,
        p0=transform_to_WM([0; 0; pitchAngle],from=:E321),
        S=vcat(fill(S1,nElementsInner),fill(S2,nElementsFFWT)),
        I=vcat(fill(I1,nElementsInner),fill(I2,nElementsFFWT)),
        hingedNodes=hingedNodes,
        hingedNodesDoF=hingedNodesDoF,
        surface=surf
        )

    # Point inertias on wing (see Table 3.2 of Healy's thesis)
    m_PI = 1e-3*[75; 75; 75; 75; 75; 56; 167]
    x1_PI = 1e-3*[70; 210; 350; 490; 630; 767; 887]
    x2_PI = -1e-3*[21; 21; 21; 21; 21; 17; 22]
    Ixx_PI = 1e-6*[73; 73; 73; 73; 73; 32; 122]
    Iyy_PI = 1e-6*[82; 82; 82; 82; 82; 26; 942]
    Izz_PI = 1e-6*[151; 151; 151; 151; 151; 56; 1057]
    for (m,x1,x2,Ixx,Iyy,Izz) in zip(m_PI,x1_PI,x2_PI,Ixx_PI,Iyy_PI,Izz_PI)
        elem = findfirst(x -> x > x1, normalizedNodalPositions * L) - 1
        pointInertia = create_PointInertia(
            η=[x1 - normalizedNodalPositions[elem] * L, x2, 0],
            mass=m,
            Ixx=Ixx,
            Iyy=Iyy,
            Izz=Izz
            )
        add_point_inertia_to_beam!(beam=beam,pointInertia=pointInertia,elementID=elem)
    end      

    # Hinge axis (defined in the local, undeformed beam basis)
    localHingeAxis = rotation_tensor_E321([-flareAngle, 0, 0]) * a2

    # Hinge axis constraint
    hingeAxisConstraints = Vector{HingeAxisConstraint}()
    if hingeConfiguration == :free
        push!(hingeAxisConstraints, create_HingeAxisConstraint(
            beam=beam,
            localHingeAxis=localHingeAxis,
            pHValue = !isnothing(foldAngle) ? 4*tan(foldAngle/4) : missing
            )
        )
    end

    # Springs around hinge
    if kSpring > 0 && hingeConfiguration == :free
        spring = create_Spring(kp=kSpring)
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[π/2,0,0])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[0,π/2,0])
    end

    # BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0),
        )

    # Beam model
    HealyBaselineFFWT = create_Model(
        name="HealyBaselineFFWT",
        beams=[beam],
        BCs=[clamp],
        v_A=t->[0; airspeed(t); 0],
        gravityVector=[0; 0; -g],
        hingeAxisConstraints=hingeAxisConstraints,
        gust=gust,
        atmosphere=atmosphere,
        units=create_UnitsSystem(frequency="Hz")
        )

    return HealyBaselineFFWT

end
export create_HealyBaselineFFWT


"""
    create_HealySideslipFFWT(; kwargs...)

Creates a version of Healy's wing with flared folding wingtip (FFWT) used for sideslip studies. See Healy's PhD thesis, chapter 7.

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `hasTipCorrection::Bool`: flag for aerodynamic tip correction
- `tipLossDecayFactor::Real`: tip loss decay factor, in case hasTipCorrection is true
- `wingtipConfiguration::Symbol`: configuration of the wingtip, either :free, :locked or :removed
- `pitchAngle::Real`: root pitch angle [rad]
- `wingtipTwist::Real`: wingtip twist angle [rad]
- `foldAngle::Union{Real,Nothing}`: fold angle [rad]
- `flareAngle::Real`: flare angle [rad]
- `kSpring::Real`: stiffness of the spring around the hinge
- `altitude::Real`: altitude
- `g::Real`: local acceleration of gravity
- `airspeed::Union{Real,<:Function}`: local airspeed
- `flightDirection::Vector{<:Real}`: flight direction vector, resolved on basis A
- `nElementsInner::Int`: number of elements for discretization of inner part of the wing
- `nElementsFFWT::Int`: number of elements for discretization of the wingtip
- `gust::Union{Nothing,Gust}=nothing`: gust
"""
function create_HealySideslipFFWT(; 
    aeroSolver::AeroSolver=Indicial(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(FlatPlate),
    updateAirfoilParameters::Bool=false,
    hasTipCorrection::Bool=true,
    tipLossDecayFactor::Real=Inf,
    wingtipConfiguration::Symbol=:free,
    pitchAngle::Real,
    wingtipTwist::Real=0,
    foldAngle::Union{Real,Nothing}=nothing,
    flareAngle::Real,
    kSpring::Real=1e-5,
    altitude::Real=0,
    g::Real=standard_atmosphere(altitude).g,
    airspeed::Union{Real,<:Function}=0,
    flightDirection::Vector{<:Real}=[0,1,0],
    nElementsInner::Int=15,
    nElementsFFWT::Int=5,
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate
    @assert wingtipConfiguration in [:free,:locked,:removed] "set `wingtipConfiguration` as :free, :locked or :removed"
    if !isnothing(foldAngle)
        @assert -π < foldAngle <= π "set foldAngle between -π and π (rad) "
    end
    @assert 0 <= flareAngle < π/4 "set flareAngle between 0 and π/4 (rad)"
    @assert kSpring >= 0 "set kSpring greater than or equal to 0"
    @assert g >= 0 "gravity must be non-negative"
    @assert altitude >= 0 "altitude must be non-negative"
    @assert length(flightDirection) == 3 "flightDirection must be a 3-element vector"
    @assert tipLossDecayFactor >= 0 "tipLossDecayFactor must be non-negative"
    @assert nElementsInner >= 8 "set at least 8 elements for inner wing"
    @assert nElementsFFWT >= 2 "set at least 2 elements for wingtip"

    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,airspeed(0)/atmosphere.a)

    # Wing (semi-)span [m]
    L = wingtipConfiguration == :removed ? 
        0.365 : 
        0.5

    # Chord [m]
    chord = 0.078

    # Normalized spar position (visually assumed from Fig. 7.3)
    normSparPos = 0.25

    # Total number of elements
    nElem = wingtipConfiguration == :removed ? 
            nElementsInner : 
            nElementsInner + nElementsFFWT

    # Hinge normalized position over (semi-)span (See Fig. 7.3 of Healy's thesis)
    hingePos = 365/500

    # Normalized nodal positions
    normalizedNodalPositions = wingtipConfiguration == :removed ? 
                               LinRange(0,L,nElementsInner+1) : 
                               unique(
                                vcat(
                                    LinRange(0, hingePos, nElementsInner+1),
                                    LinRange(hingePos, 1, nElementsFFWT+1)
                                    )
                                )

    # Spar geometric and material properties
    b = 19e-3
    h = 4.8e-3
    E = 70e9
    G = 27e9
    ρ = 2.7e3
    A = b*h
    Iy = b*h^3/12
    Iz = h*b^3/12
    J = Iy+Iz

    # Stiffness correction factors
    KOOP = 1.15
    KIP = 1
    KT = 0.2

    # Baseline stiffness matrix
    S1 = isotropic_stiffness_matrix(
        EA=E*A,
        GJ=G*J*KT,
        EIy=E*Iy*KOOP,
        EIz=E*Iz*KIP,
        ∞=1e8
        )

    # Wingtip baseline stiffness matrix (assume 50% of central spar's stiffness)
    S2 = 0.5*S1

    # Central spar inertia matrix
    I1 = inertia_matrix(
        ρA=ρ*A,
        ρIy=ρ*Iy,
        ρIz=ρ*Iz
        )

    # Wingtip baseline inertia matrix (assume 10% of central spar's inertia)
    I2 = inertia_matrix(
        ρA=0.1*ρ*A,
        ρIy=0.1*ρ*Iy,
        ρIz=0.1*ρ*Iz
        )

    # Array of sectional stiffness matrices
    S = wingtipConfiguration == :removed ? 
        fill(S1,nElementsInner) : 
        vcat(fill(S1,nElementsInner), fill(S2,nElementsFFWT))

    # Array of sectional inertia matrices
    I = wingtipConfiguration == :removed ? 
        fill(I1,nElementsInner) : 
        vcat(fill(I1,nElementsInner), fill(I2,nElementsFFWT))

    # Twist angle distribution over wingspan (given in Fig. 7.32)
    φ_of_x1 = wingtipConfiguration == :free ? 
              x1 -> ifelse(x1<=(0.365+0.06), 
                           0.0, 
                           ifelse(x1>=0.365+0.085, 
                                  wingtipTwist, 
                                  wingtipTwist*(x1-0.365-0.06)/0.025
                                  )
                          ) : 
              x1 -> 0.0

    # Aerodynamic surface
    surf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=airfoil,
        c=chord,
        φ=φ_of_x1,
        normSparPos=normSparPos,
        hasTipCorrection=hasTipCorrection,
        tipLossDecayFactor=tipLossDecayFactor,
        updateAirfoilParameters=updateAirfoilParameters
        )

    # Hinge node
    hingeNode = nElementsInner + 1    

    # Hinged nodes and DOFs arrays
    hingedNodes = wingtipConfiguration == :free ? 
                  [hingeNode] : 
                  Vector{Int}()
    hingedNodesDoF = wingtipConfiguration == :free ? 
                     [trues(3)] : 
                     [falses(3)]

    # Wing beam
    beam = create_Beam(
        name="beam",
        length=L,
        nElements=nElem,
        normalizedNodalPositions=normalizedNodalPositions,
        S=S,
        I=I,
        p0=transform_to_WM([0, 0, pitchAngle],from=:E321),
        surface=surf,
        hingedNodes=hingedNodes,
        hingedNodesDoF=hingedNodesDoF,
        )

    # Point inertias on wing (see Table 7.2 of Healy's thesis)
    m_PI = 1e-3*[20.3; 29.6; 29.6; 23.0; 61.0]
    x1_PI = 1e-3*[17; 155; 255; 328; 424]
    x2_PI = 1e-3*[-12.3; -12.3; -12.3; -12.4; -13.9]
    Ixx_PI = 1e-6*[9.1; 13.4; 13.4; 6.5; 18.9]
    Iyy_PI = 1e-6*[8.4; 24.8; 24.8; 2.8; 116]
    Izz_PI = Ixx_PI + Iyy_PI
    for (m,x1,x2,Ixx,Iyy,Izz) in zip(m_PI,x1_PI,x2_PI,Ixx_PI,Iyy_PI,Izz_PI)
        elem = findfirst(x -> x > x1, normalizedNodalPositions * L) - 1
        pointInertia = create_PointInertia(
            η=[x1 - normalizedNodalPositions[elem] * L, x2, 0],
            mass=m,
            Ixx=Ixx,
            Iyy=Iyy,
            Izz=Izz
            )
        add_point_inertia_to_beam!(beam=beam,pointInertia=pointInertia,elementID=elem)
    end

    # Hinge axis (defined in the local, undeformed beam basis)
    localHingeAxis = rotation_tensor_E321([-flareAngle, 0, 0]) * a2

    # Hinge axis constraint
    hingeAxisConstraints = Vector{HingeAxisConstraint}()
    if wingtipConfiguration == :free
        push!(hingeAxisConstraints, create_HingeAxisConstraint(
            beam=beam,
            localHingeAxis=localHingeAxis,
            pHValue = !isnothing(foldAngle) ? 4*tan(foldAngle/4) : missing
            )
        )
    end

    # Springs around hinge
    if kSpring > 0
        spring = create_Spring(kp=kSpring)
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[π/2,0,0])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[0,π/2,0])
    end

    # BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0),
        )

    # Basis A velocity vector
    v_A(t) = airspeed(t) * flightDirection / norm(flightDirection)

    # Wing model
    HealySideslipFFWT = create_Model(
        name="HealySideslipFFWT",
        beams=[beam],
        BCs=[clamp],
        gravityVector=[0; 0; -g],
        v_A=v_A,
        hingeAxisConstraints=hingeAxisConstraints,
        gust=gust,
        units=create_UnitsSystem(frequency="Hz")
        )

    return HealySideslipFFWT

end
export create_HealySideslipFFWT


"""
    create_HealyLCOFFWT(; kwargs...)

Creates a version of Healy's wing with flared folding wingtip (FFWT) used for flutter and LCOs studies. See Healy's PhD thesis, chapter 7 and CHEUNG et al. (2020) - Analyzing the Dynamic Behavior of a High Aspect Ratio Wing Incorporating a Folding Wingtip [SCITECH].

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `hasTipCorrection::Bool`: flag for aerodynamic tip correction
- `tipLossDecayFactor::Union{Real,Missing}`: tip loss decay factor, in case hasTipCorrection is true
- `verticalMount::Bool`: flag for wing mounted vertically
- `wingtipRemoved::Bool`: flag for removed wingtip
- `wingtipLocked::Bool`: flag for locked wingtip
- `wingtipTab::Bool`: flag for wingtip tab
- `pitchAngle::Real`: root pitch angle [rad]
- `foldAngle::Union{Real,Nothing}`: fold angle [rad]
- `flareAngle::Real`: flare angle [rad]
- `kSpring::Real`: stiffness of the spring around the hinge
- `altitude::Real`: altitude
- `g::Real`: local acceleration of gravity
- `airspeed::Union{Real,<:Function}`: local airspeed
- `flightDirection::Vector{<:Real}`: flight direction vector, resolved on basis A
- `nElementsInner::Int`: number of elements for discretization of inner part of the wing
- `nElementsOuter::Int`: number of elements for discretization of outer part of the wing
- `nElementsFFWT::Int`: number of elements for discretization of the wingtip
- `gust::Union{Nothing,Gust}=nothing`: gust
"""
function create_HealyLCOFFWT(; 
    aeroSolver::AeroSolver=Indicial(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(FlatPlate),
    updateAirfoilParameters::Bool=false,
    hasTipCorrection::Bool=true,
    tipLossDecayFactor::Union{Real,Missing}=missing,
    verticalMount::Bool=false,
    wingtipRemoved::Bool=false,
    wingtipLocked::Bool=false,
    wingtipTab::Bool=true,
    pitchAngle::Real=0,
    foldAngle::Union{Real,Nothing}=nothing,
    flareAngle::Real=deg2rad(10),
    kSpring::Real=1e-5,
    altitude::Real=0,
    g::Real=standard_atmosphere(altitude).g,
    airspeed::Union{Real,<:Function}=0,
    flightDirection::Vector{<:Real}=[0,1,0],
    nElementsInner::Int=8,
    nElementsOuter::Int=2,
    nElementsFFWT::Int=4,
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate
    isnothing(foldAngle) || @assert -π < foldAngle ≤ π "set foldAngle between -π and π (rad)"
    @assert 0 <= flareAngle < π/4 "set flareAngle between 0 and π/4 (rad)"
    @assert kSpring >= 0 "set kSpring greater than or equal to 0"
    @assert g >= 0 "gravity must be non-negative"
    @assert altitude >= 0 "altitude must be non-negative"
    @assert length(flightDirection) == 3 "flightDirection must be a 3-element vector"
    ismissing(tipLossDecayFactor) || @assert tipLossDecayFactor ≥ 0 "tipLossDecayFactor must be non-negative"

    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Flags
    wingtipLocked = wingtipLocked || wingtipRemoved
    wingtipTab = wingtipTab && !wingtipRemoved

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,airspeed(0)/atmosphere.a)

    # Compute tip loss decay factor, if applicable
    tipLossDecayFactor = ismissing(tipLossDecayFactor) ? 
                         tip_loss_function_HealyLCOFFWT(pitchAngle,airspeed(0)) : 
                         tipLossDecayFactor

    # Wing (semi-)span [m]
    L = wingtipRemoved ? 
        1 : 
        1.345

    # Chord [m]
    chord = 0.15

    # Chord-normalized spar position
    normSparPos = 0.25    

    # Spar transition normalized position
    sparTransitionPos = 0.875/L

    # Hinge normalized position over (semi-)span
    hingePos = 1.0/L

    # Normalized tab span and position over the chord
    tabSpan = [L-0.16; L]/L
    tabPos = 0.75

    # Total number of elements
    nElem = wingtipRemoved ? 
            nElementsInner + nElementsOuter : 
            nElementsInner + nElementsOuter + nElementsFFWT

    # Normalized nodal positions
    normalizedNodalPositions = wingtipRemoved ? 
        unique(
            vcat(
                LinRange(0, sparTransitionPos, nElementsInner+1),
                LinRange(sparTransitionPos, hingePos, nElementsOuter+1)
                )
            ) : 
        unique(
            vcat(
                LinRange(0, sparTransitionPos, nElementsInner+1),
                LinRange(sparTransitionPos, hingePos, nElementsOuter+1),
                LinRange(hingePos, 1, nElementsFFWT+1)
                )
            )

    # Hinge node
    hingeNode = nElementsInner + nElementsOuter + 1

    # Inner spar geometric and material properties
    b1 = 30e-3
    h1 = 5e-3
    E1 = 193e9
    G1 = 76e9
    ρ1 = 8e3
    A1 = b1*h1
    Iy1 = b1*h1^3/12
    Iz1 = h1*b1^3/12
    J1 = Iy1+Iz1

    # Stiffness correction factors for inner spar
    KOOP1 = 1.1
    KIP1 = 1.1
    KT1 = 0.17

    # Outer spar geometric and material properties
    b2 = 19.3e-3
    h2 = 13.3e-3
    E2 = 1.65e9
    G2 = 0.7e9
    ρ2 = 3.3e3
    A2 = b2*h2
    Iy2 = b2*h2^3/12
    Iz2 = h2*b2^3/12
    J2 = Iy2+Iz2

    # Stiffness correction factors for outer spar
    KOOP2 = 2
    KIP2 = 2
    KT2 = 0.5

    # Inner spar stiffness matrix
    S1 = isotropic_stiffness_matrix(
        EA=E1*A1,
        GJ=G1*J1*KT1,
        EIy=E1*Iy1*KOOP1,
        EIz=E1*Iz1*KIP1,
        ∞=1e8
        )

    # Outer spar stiffness matrix
    S2 = isotropic_stiffness_matrix(
        EA=E2*A2,
        GJ=G2*J2*KT2,
        EIy=E2*Iy2*KOOP2,
        EIz=E2*Iz2*KIP2,
        ∞=1e8
        )

    # Wingtip stiffness matrix (assume fraction of inner spar)
    S3 = 0.5*S1

    # Inner spar inertia matrix
    I1 = inertia_matrix(
        ρA=ρ1*A1,
        ρIy=ρ1*Iy1,
        ρIz=ρ1*Iz1
        )

    # Outer spar inertia matrix
    I2 = inertia_matrix(
        ρA=ρ2*A2,
        ρIy=ρ2*Iy2,
        ρIz=ρ2*Iz2
        )

    # Baseline wingtip inertia matrix (assume small fraction of inner spar)
    γ = 0.1
    I3 = inertia_matrix(
        ρA=γ*ρ1*A1,
        ρIy=γ*ρ1*Iy1,
        ρIz=γ*ρ1*Iz1,
        )

    # Array of sectional stiffness and inertia matrices
    S = wingtipRemoved ? 
        vcat(
            fill(S1,nElementsInner),
            fill(S2,nElementsOuter)
            ) : 
        vcat(
            fill(S1,nElementsInner),
            fill(S2,nElementsOuter),
            fill(S3,nElementsFFWT)
            )
    I = wingtipRemoved ? 
        vcat(
            fill(I1,nElementsInner),
            fill(I2,nElementsOuter)
            ) : 
        vcat(
            fill(I1,nElementsInner),
            fill(I2,nElementsOuter),
            fill(I3,nElementsFFWT)
            )

    # Aerodynamic surface
    surf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=airfoil,
        c=chord,
        normSparPos=normSparPos,
        normFlapPos = wingtipTab ? tabPos : nothing,
        normFlapSpan = wingtipTab ? tabSpan : nothing,
        hasTipCorrection=hasTipCorrection,
        tipLossDecayFactor=tipLossDecayFactor,
        updateAirfoilParameters=updateAirfoilParameters
        )  

    # Hinged nodes and DOFs arrays
    hingedNodes = wingtipLocked ? 
                  Vector{Int}() : 
                  [hingeNode]
    hingedNodesDoF = wingtipLocked ? 
                     [falses(3)] : 
                     [trues(3)]

    # Rotation parameters from basis A to basis b
    p0 = verticalMount ? 
         [0, -π/2, pitchAngle] : 
         [0, 0, pitchAngle]

    # Wing beam
    beam = create_Beam(
        name="beam",
        length=L,
        nElements=nElem,
        normalizedNodalPositions=normalizedNodalPositions,
        S=S,
        I=I,
        p0=transform_to_WM(p0,from=:E321),
        surface=surf,
        hingedNodes=hingedNodes,
        hingedNodesDoF=hingedNodesDoF,
        )        

    # Inner pannels point inertias
    mInner = 142e-3*ones(8)
    x1Inner = 1e-3 * (LinRange(87.5,822.5,8) .- 35)
    x2Inner = -30.5e-3*ones(8)
    IxxInner = 171e-6*ones(8)
    IyyInner = 85e-6*ones(8)

    # Outer (hinge) pannel point inertia
    mOuter = 353e-3
    x1Outer = 1e-3 * (952 - 35)
    x2Outer = -23.2e-3
    IxxOuter = 269e-6
    IyyOuter = 355e-6

    # Wingtip point inertias
    mFWT = wingtipTab ? 563e-3 : 520e-3
    x1FWT = wingtipTab ? 1e-3 * (1159.6 - 35) : 1e-3 * (1179.9 - 35)
    x2FWT = wingtipTab ? -47.5e-3 : -23.8e-3
    IxxFWT = wingtipTab ? 615e-6 : 680e-6
    IyyFWT = wingtipTab ? 2162e-6 : 3858e-6

    # Point inertias on wing (see Table 5.2 of Healy's thesis)
    m_PI   = wingtipRemoved ? vcat(mInner, mOuter) : vcat(mInner, mOuter, mFWT)
    x1_PI  = wingtipRemoved ? vcat(x1Inner, x1Outer) : vcat(x1Inner, x1Outer, x1FWT)
    x2_PI  = wingtipRemoved ? vcat(x2Inner, x2Outer) : vcat(x2Inner, x2Outer, x2FWT)
    Ixx_PI = wingtipRemoved ? vcat(IxxInner, IxxOuter) : vcat(IxxInner, IxxOuter, IxxFWT)
    Iyy_PI = wingtipRemoved ? vcat(IyyInner, IyyOuter) : vcat(IyyInner, IyyOuter, IyyFWT)
    Izz_PI = Ixx_PI + Iyy_PI
    for (m,x1,x2,Ixx,Iyy,Izz) in zip(m_PI,x1_PI,x2_PI,Ixx_PI,Iyy_PI,Izz_PI)
        elem = findfirst(x -> x > x1, normalizedNodalPositions * L) - 1
        pointInertia = create_PointInertia(
            η=[x1 - normalizedNodalPositions[elem] * L, x2, 0],
            mass=m,
            Ixx=Ixx,
            Iyy=Iyy,
            Izz=Izz
            )
        add_point_inertia_to_beam!(beam=beam,pointInertia=pointInertia,elementID=elem)
    end

    # Hinge axis (defined in the local, undeformed beam basis)
    localHingeAxis = rotation_tensor_E321([-flareAngle, 0, 0]) * a2

    # Hinge axis constraint
    hingeAxisConstraints = Vector{HingeAxisConstraint}()
    if !wingtipLocked
        push!(hingeAxisConstraints, create_HingeAxisConstraint(
            beam=beam,
            localHingeAxis=localHingeAxis,
            pHValue = !isnothing(foldAngle) ? 4*tan(foldAngle/4) : missing
            )
        )
    end

    # Springs around hinge
    if kSpring > 0 && !wingtipLocked
        spring = create_Spring(kp=kSpring)
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[π/2,0,0])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[0,π/2,0])
    end

    # BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0),
        )

    # Basis A velocity vector
    v_A(t) = airspeed(t) * flightDirection / norm(flightDirection)

    # Wing model
    HealyLCOFFWT = create_Model(
        name="HealyLCOFFWT",
        beams=[beam],
        BCs=[clamp],
        gravityVector=[0, 0, -g],
        v_A=v_A,
        hingeAxisConstraints=hingeAxisConstraints,
        gust=gust,
        units=create_UnitsSystem(frequency="Hz")
        )

    return HealyLCOFFWT

end
export create_HealyLCOFFWT