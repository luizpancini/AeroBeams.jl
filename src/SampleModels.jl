using DelimitedFiles

# Returns the normalized nodal positions of the Pazy wing
function nodal_positions_Pazy()
    return [0.0; 0.06956521653730675; 0.13913043671201064; 0.208695655068016; 0.2782608734240213; 0.34782609178002666; 0.41739131195473056; 0.4869565303107358; 0.5565217486667412; 0.626086968841445; 0.6956521871974504; 0.7652174055534556; 0.8347826239094611; 0.9043478440841649; 0.9652173080712125; 1.0]
end


# Returns the sectional stiffness matrices of the Pazy wing
function stiffness_matrices_Pazy(GAy::Real,GAz::Real)

    @assert GAy > 0
    @assert GAz > 0

    # Load elemental sectional stiffness values
    # Axial stiffness
    EA = [9.79449259e06 9.66669055e06 9.64764696e06 9.64877085e06 9.64996612e06 9.65062048e06 9.65094702e06 9.65092085e06 9.65060878e06 9.64997519e06 9.64887183e06 9.64715800e06 9.65344706e06 9.69921416e06 1.00196328e07]
    # Torsional stiffness
    GJ = [7.58259714e00 7.58259714e00 7.58259714e00 7.58259714e00 6.51183018e00 6.51183018e00 6.51183018e00 6.51183018e00 6.51583594e00 6.51583594e00 6.51583594e00 6.51583594e00 7.11892627e00 7.11892627e00 1.73730728e01]
    # Out-of-plane bending stiffness
    EIy = [5.24743501e00 4.50448461e00 4.47822861e00 4.47438433e00 4.47481272e00 4.47492882e00 4.47492903e00 4.47492905e00 4.47492889e00 4.47488638e00 4.47471764e00 4.47618354e00 4.48916524e00 4.54819972e00 4.77037999e00]
    # In-plane bending stiffness
    EIz = [3.31757932e03 3.27862894e03 3.28287728e03 3.28653821e03 3.28881887e03 3.29012279e03 3.29071820e03 3.29071311e03 3.29015111e03 3.28892576e03 3.28685980e03 3.28351136e03 3.27904566e03 3.27164430e03 3.37372291e03]
    # Axial-torsion coupling stiffness
    c_EA_GJ = [-5.69828967e-01 -5.69828967e-01 -5.69828967e-01 -5.69828967e-01 -2.44730557e-01 -2.44730557e-01 -2.44730557e-01 -2.44730557e-01 5.55131815e-01 5.55131815e-01 5.55131815e-01 5.55131815e-01 1.14679646e00 1.14679646e00 6.11339950e00]
    # Axial-OOP bending coupling stiffness
    c_EA_EIy = [-1.37141817e00 -1.85955625e00 -2.03745676e00 -2.36861943e00 -2.34621942e00 -2.33093397e00 -2.32194657e00 -2.31978268e00 -2.32431249e00 -2.33594699e00 -2.35794054e00 -2.42282990e00 -2.63402308e00 -2.79279201e00 -2.28331802e00]
    # Axial-IP bending coupling stiffness
    c_EA_EIz = [5.44855583e04 5.35354102e04 5.32071297e04 5.31204588e04 5.30688769e04 5.30373657e04 5.30224268e04 5.30209363e04 5.30331840e04 5.30603697e04 5.31061305e04 5.31755122e04 5.33603082e04 5.41016226e04 5.44519115e04]
    # Torsion-OOP bending coupling stiffness
    c_GJ_EIy = [9.33080027e-02 9.33080027e-02 9.33080027e-02 9.33080027e-02 4.22659591e-05 4.22659591e-05 4.22659591e-05 4.22659591e-05 -1.02012319e-03 -1.02012319e-03 -1.02012319e-03 -1.02012319e-03 -8.09780725e-02 -8.09780725e-02 -1.39455813e00]
    # Torsion-IP bending coupling stiffness
    c_GJ_EIz = [1.52918906e-02 1.52918906e-02 1.52918906e-02 1.52918906e-02 3.94787982e-03 3.94787982e-03 3.94787982e-03 3.94787982e-03 -8.69600350e-03 -8.69600350e-03 -8.69600350e-03 -8.69600350e-03 1.01277971e-02 1.01277971e-02 4.24981724e-01]
    # OOP-IP bending coupling stiffness
    c_EIy_EIz = [-1.17141160e-01 -1.12442858e-01 -1.15192612e-01 -1.18950224e-01 -1.20649484e-01 -1.21389557e-01 -1.21682056e-01 -1.21652571e-01 -1.21317479e-01 -1.20584487e-01 -1.19291993e-01 -1.17184960e-01 -1.14807975e-01 -1.16130851e-01 -1.14249441e-01]

    # Set matrices
    S = [[    EA[i]  0    0  c_EA_GJ[i]  c_EA_EIy[i]  c_EA_EIz[i];
                  0 GAy   0           0            0            0;
                  0  0  GAz           0            0            0;
         c_EA_GJ[i]  0    0       GJ[i]  c_GJ_EIy[i]  c_GJ_EIz[i];
        c_EA_EIy[i]  0    0 c_GJ_EIy[i]       EIy[i] c_EIy_EIz[i];
        c_EA_EIz[i]  0    0 c_GJ_EIz[i] c_EIy_EIz[i]       EIz[i]] for i in 1:15]

    return S
end


# Returns the sectional inertia matrices of the Pazy wing
function inertia_matrices_Pazy()

    # Length and nodal position
    L = 0.549843728 
    nodalPositions = nodal_positions_Pazy()

    # Load nodal inertia values
    # Point mass
    m = [1.91186108e-02 2.06788141e-02 2.06828340e-02 2.06828274e-02 2.06828571e-02 2.06828792e-02 2.06828498e-02 2.06828527e-02 2.06828616e-02 2.06828621e-02 2.06828626e-02 2.06828630e-02 2.06828635e-02 2.06828635e-02 2.54479524e-02 4.31558925e-02]
    # Center of mass offset from reference axis in x1 axis direction
    m_x1 = [6.44284658e-03 -1.87501083e-03 -1.88019018e-03 -1.88021112e-03 -1.88022963e-03 -1.88021127e-03 -1.88020545e-03 -1.88022219e-03 -1.88022845e-03 -1.88022769e-03 -1.88022693e-03 -1.88022617e-03 -1.88022541e-03 -1.88022469e-03 1.34295741e-03 3.28782957e-03]
    # Center of mass offset from reference axis in x2 axis direction
    m_x2 = [6.09115952e-06 -9.03114429e-04 -9.13241152e-04 -9.13181107e-04 -9.13136269e-04 -9.13139373e-04 -9.13148460e-04 -9.13133755e-04 -9.13138940e-04 -9.13138146e-04 -9.13137351e-04 -9.13136557e-04 -9.13135762e-04 -9.13135762e-04 8.11101073e-04 -5.09272488e-03]
    # Center of mass offset  from reference axis in x3 axis direction
    m_x3 = [-2.71488305e-05 -1.48960084e-05 -1.49441628e-05 -1.50172024e-05 -1.50375210e-05 -1.50375096e-05 -1.50376952e-05 -1.50378716e-05 -1.50380388e-05 -1.50382121e-05 -1.50383854e-05 -1.50385588e-05 -1.50387321e-05 -1.50387321e-05 -3.00392326e-05 -1.43641715e-04]
    # x1 axis mass moment of inertia
    Ixx = [1.21416697e-05 1.09721424e-05 1.09788396e-05 1.09788493e-05 1.09788664e-05 1.09788786e-05 1.09788644e-05 1.09788687e-05 1.09788717e-05 1.09788723e-05 1.09788729e-05 1.09788735e-05 1.09788741e-05 1.09788741e-05 1.45515348e-05 1.22200220e-04]
    # x1-x3 axes mass product of inertia
    Ixy = [1.29099112e-08 -4.29590762e-08 -4.66241072e-08 -4.66239396e-08 -4.66186558e-08 -4.66208462e-08 -4.66181019e-08 -4.66178877e-08 -4.66200515e-08 -4.66200824e-08 -4.66201133e-08 -4.66201442e-08 -4.66201750e-08 -4.66201750e-08 -2.55761995e-07 3.06994479e-07]
    # x1-x3 axes mass product of inertia
    Ixz = [2.61326836e-10 -4.05020470e-10 -4.09113627e-10 -4.09779617e-10 -4.05575499e-10 -4.05572069e-10 -4.05576866e-10 -4.05587362e-10 -4.05596285e-10 -4.05603025e-10 -4.05609766e-10 -4.05616506e-10 -4.05623247e-10 -4.05623246e-10 -2.22503423e-09 -7.38220857e-09]
    # x2 axis mass moment of inertia
    Iyy = [6.63182204e-07 2.08434944e-06 2.08700528e-06 2.08700318e-06 2.08700972e-06 2.08701867e-06 2.08700826e-06 2.08700778e-06 2.08701031e-06 2.08701035e-06 2.08701038e-06 2.08701041e-06 2.08701047e-06 2.08701047e-06 1.68654882e-06 8.76965373e-07]
    # x2-x3 axes mass product of inertia
    Iyz = [9.44646195e-09 2.45046373e-09 2.45271020e-09 2.45762425e-09 2.45962684e-09 2.45962929e-09 2.45978322e-09 2.45991134e-09 2.46005086e-09 2.46019166e-09 2.46033249e-09 2.46047329e-09 2.46061128e-09 2.46061778e-09 1.47231285e-08 1.29956112e-07]
    # x3-axis mass moment of inertia
    Izz = [1.20937988e-05 1.28286585e-05 1.28380090e-05 1.28380097e-05 1.28380332e-05 1.28380544e-05 1.28380298e-05 1.28380335e-05 1.28380389e-05 1.28380395e-05 1.28380401e-05 1.28380406e-05 1.28380412e-05 1.28380412e-05 1.55828514e-05 1.22614163e-04]

    # Set matrices
    I = Vector{Matrix{Float64}}(undef,15)
    for n=2:16
        Δℓ = L*(nodalPositions[n]-nodalPositions[n-1])
        η = [Δℓ/2+m_x1[n]; m_x2[n]; m_x3[n]]
        inertiaMatrix = [Ixx[n] Ixy[n] Ixz[n]; Ixy[n] Iyy[n] Iyz[n]; Ixz[n] Iyz[n] Izz[n]]
        I[n-1] = 1/Δℓ*[       m[n]*I3                -m[n]*tilde(η);
                        m[n]*tilde(η) inertiaMatrix-m[n]*tilde(η)^2]
        if n==2
            η = [-Δℓ/2+m_x1[1]; m_x2[1]; m_x3[1]]
            inertiaMatrix = [Ixx[1] Ixy[1] Ixz[1]; Ixy[1] Iyy[1] Iyz[1]; Ixz[1] Iyz[1] Izz[1]]
            I[1] += 1/Δℓ*[        m[1]*I3                -m[1]*tilde(η);
                            m[1]*tilde(η) inertiaMatrix-m[1]*tilde(η)^2]
        end
    end

    return I
end


"""
    tip_loss_factor_Pazy(θ::Real,U::Real)

Computes the tip loss factor for the Pazy wing's tip correction function

# Arguments
- `θ::Real` = root pitch angle, in degrees
- `U::Real` = airspeed
"""
function tip_loss_factor_Pazy(θ::Real,U::Real)

    # Bound inputs
    θ = min(7,max(θ,0))
    U = min(60,max(U,0))

    # Coefficients as a function of root angle of attack
    θRange = [0; 1; 2; 3; 4; 5; 6; 7]
    τ₀Range = [6.58; 6.29; 6.00; 5.92; 5.91; 6.08; 6.43; 6.88]
    τ₁Range = 1e-2*[0; 3.33; 5.19; 6.01; 6.49; 5.98; 4.43; 2.30]
    τ₂Range = -1e-4*[0; 5.56; 8.59; 10.6; 12.3; 12.8; 11.9; 10.2]
    τ₀ = interpolate(θRange,τ₀Range,θ)
    τ₁ = interpolate(θRange,τ₁Range,θ)
    τ₂ = interpolate(θRange,τ₂Range,θ)
    
    return τ₀ + τ₁*U + τ₂*U^2
end
export tip_loss_factor_Pazy


"""
    geometrical_properties_Pazy()

Returns the fixed geometrical (and discretization) properties of the Pazy wing

"""
function geometrical_properties_Pazy()

    # Number of elements
    nElem = 15

    # Length
    L = 0.549843728

    # Chord
    chord = 0.0989

    # Normalized spar position
    normSparPos = 0.44096

    return nElem,L,chord,normSparPos
end
export geometrical_properties_Pazy


"""
    typical_section_data(name::String)

Returns the data of the typical section of given name

# Arguments
- `name::String` = name of the typical section

# Outputs
- `a` = position of elastic axis (semichords aft of midchord)
- `e` = position of mass axis (semichords aft of midchord)
- `μ` = mass ratio
- `rα²` = squared radius of gyration in pitch 
- `σ` = ratio of plunge to pitch frequency in vacuo
- `ωα` = pitch frequency in vacuo [rad/s]
- `c` = chord [m]
"""
function typical_section_data(name::String)

    if name == "HP-1" # Defined in problem 5.5 of Hodges & Pierce
        a = -1/5
        e = -1/10
        μ = 20
        rα² = 6/25
        σ = 2/5
        ωα = 50 # A numerical value was not given
        c = 1   # A numerical value was not given
    elseif name == "HP-2" # Defined in problem 5.9 of Hodges & Pierce
        a = -1/3
        e = -1/10
        μ = 50
        rα² = 4/25
        σ = 2/5
        ωα = 40 # A numerical value was not given
        c = 1   # A numerical value was not given
    elseif name == "Fung" # Defined by Fung - An Introduction to Theory of Aeroelasticity
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
    create_Pazy(; kwargs...)

Creates the Pazy wing

Returns the wing model and geometrical properties

# Arguments
- `aeroSolver::AeroSolver` = aerodynamic solver
- `gustLoadsSolver::GustAeroSolver` = indicial gust loads solver
- `derivationMethod::DerivationMethod` = method for aerodynamic derivatives
- `updateAirfoilParameters::Bool` = flag to update airfoil parameters with airspeed
- `upright::Bool` = flag to set the wing in the upright position
- `airfoil::Airfoil` = airfoil section
- `θ::Real` = pitch angle at the root [rad]
- `Λ::Real` = sweep angle [rad]
- `withTipCorrection::Bool` = flag for aerodynamic tip correction 
- `GAy::Real` = shear stiffness in the x2 direction
- `GAz::Real` = shear stiffness in the x3 direction
- `altitude::Real` = altitude
- `g::Real` = acceleration of gravity
- `airspeed::Real` = airspeed
- `tipMass::Real` = mass of a point inertia added to the tip of the wing
- `ξtipMass::Vector{<:Real}` = position vector of the tip mass relative to the tip of the spar, resolved in the local basis
- `tipMassInertia::Matrix{<:Real}` = mass moment of inertia matrix of the tip mass
- `additionalBCs::Vector{BC}` = additional BCs (beyond the clamp)
- `gust::Union{Nothing,Gust}` = gust
"""
function create_Pazy(; aeroSolver::AeroSolver=Indicial(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),derivationMethod::DerivationMethod=AD(),updateAirfoilParameters::Bool=true,upright::Bool=false,airfoil::Airfoil=deepcopy(NACA0018),θ::Real=0,Λ::Real=0,withTipCorrection::Bool=true,GAy::Real=1e16,GAz::Real=GAy,altitude::Real=0,g::Real=-9.80665,airspeed::Real=0,tipMass::Real=0,ξtipMass::Vector{<:Real}=zeros(3),tipMassInertia::Matrix{<:Real}=zeros(3,3),additionalBCs::Vector{BC}=Vector{BC}(),gust::Union{Nothing,Gust}=nothing)

    # Validate
    @assert altitude >= 0
    @assert airspeed >= 0
    @assert GAy >= 1e6
    @assert GAz >= 1e6
    @assert tipMass >= 0
    @assert length(ξtipMass) == 3

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Fixed geometrical and discretization properties
    nElem,L,chord,normSparPos = geometrical_properties_Pazy()

    # Normalized nodal positions
    nodalPositions = nodal_positions_Pazy()

    # Stiffness matrices
    S = stiffness_matrices_Pazy(GAy,GAz)

    # Inertia matrices
    I = inertia_matrices_Pazy()

    # Tip loss factor
    τ = tip_loss_factor_Pazy(θ*180/pi,airspeed)

    # Rotation parameters from basis A to basis b
    p0 = upright ? [-Λ; -π/2; θ] : [-Λ; 0; θ]

    # Update airfoil parameters
    update_Airfoil_params!(airfoil,Ma=airspeed/atmosphere.a,U=airspeed,b=chord/2)

    # Aerodynamic surface
    surf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=chord,normSparPos=normSparPos,hasTipCorrection=withTipCorrection,tipLossDecayFactor=τ,updateAirfoilParameters=updateAirfoilParameters)

    # Tip mass
    tipInertia = PointInertia(elementID=nElem,η=[L/nElem/2+ξtipMass[1];ξtipMass[2];ξtipMass[3]],mass=tipMass,inertiaMatrix=tipMassInertia)

    # Wing beam
    beam = create_Beam(name="wingBeam",length=L,nElements=nElem,normalizedNodalPositions=nodalPositions,S=S,I=I,rotationParametrization="E321",p0=p0,aeroSurface=surf,pointInertias=[tipInertia])

    # Update beam of additional BCs, if applicable
    for BC in additionalBCs
        BC.beam = beam
    end
    BCs = additionalBCs

    # Add root clamp to BCs
    clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
    push!(BCs,clamp)

    # Model
    pazy = create_Model(name="Pazy",beams=[beam],BCs=BCs,gravityVector=[0;0;g],v_A=[0;airspeed;0],gust=gust,units=create_UnitsSystem(frequency="Hz"))

    return pazy,nElem,L,chord,normSparPos
end
export create_Pazy


"""
    create_PazyFFWT(; kwargs...)

Creates a version of the Pazy wing with flared folding wingtip (FFWT)

# Arguments
- `p0::Vector{<:Real}` = initial rotation parameters in the Euler 3-2-1 sequence (yaw-pitch-roll)
- `airfoil::Airfoil` = airfoil section
- `aeroSolver::AeroSolver` = aerodynamic solver
- `gustLoadsSolver::GustAeroSolver` = indicial gust loads solver
- `derivationMethod::DerivationMethod` = method for aerodynamic derivatives
- `withTipCorrection::Bool` = flag for aerodynamic tip correction 
- `GAy::Real` = shear stiffness in the x2 direction
- `GAz::Real` = shear stiffness in the x3 direction
- `hingeNode::Int64` = hinge node
- `foldAngle::Union{Real,Nothing}` = fold angle [rad]
- `flareAngle::Real` = flare angle [rad]
- `kSpring::Real` = stiffness of the hinge
- `g::Real` = local acceleration of gravity
- `altitude::Real` = altitude
- `airspeed::Real` = local airspeed
"""
function create_PazyFFWT(; p0::Vector{<:Real}=zeros(3),airfoil::Airfoil=deepcopy(NACA0018),aeroSolver::AeroSolver=Indicial(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),derivationMethod::DerivationMethod=AD(),withTipCorrection::Bool=true,GAy::Real=1e16,GAz::Real=GAy,hingeNode::Int64=14,foldAngle::Union{Real,Nothing}=nothing,flareAngle::Real=0,kSpring::Real=1e-4,g::Real=9.80665,altitude::Real=0,airspeed::Real=0)

    # Validate
    @assert 2 <= hingeNode <= 15 "hingeNode must be between 2 and 15"
    if !isnothing(foldAngle)
        @assert -π < foldAngle <= π "set foldAngle between -π and π (rad) "
    end
    @assert 0 <= flareAngle < π/4 "set flareAngle between 0 and π/4 (rad)"
    @assert kSpring >= 0
    @assert g >= 0
    @assert airspeed >= 0

    # Flags
    foldAngleIsInput = !isnothing(foldAngle)

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Total length 
    L = 0.549843728

    # Number of elements
    nElem = 15

    # Normalized nodal positions
    nodalPositions = nodal_positions_Pazy()

    # Stiffness matrices
    S = stiffness_matrices_Pazy(GAy,GAz)

    # Inertia matrices
    I = inertia_matrices_Pazy()

    # Chord
    chord = 0.0989

    # Normalized spar position
    normSparPos = 0.44096

    # Elements inboard and outboard of the hinge
    inboardElem = hingeNode-1
    outboardElem = hingeNode

    # Update airfoil parameters
    update_Airfoil_params!(airfoil,Ma=airspeed/atmosphere.a,U=airspeed,b=chord/2)

    # Tip loss factor
    τ = tip_loss_factor_Pazy(p0[3]*180/π,airspeed)

    # Aerodynamic surface
    surf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=chord,normSparPos=normSparPos,hasTipCorrection=withTipCorrection,tipLossDecayFactor=τ)

    # Beams 
    beam = create_Beam(name="beam",length=L,nElements=nElem,normalizedNodalPositions=nodalPositions,S=S,I=I,rotationParametrization="E321",p0=p0,aeroSurface=surf,hingedNodes=[hingeNode],hingedNodesDoF=[trues(3)])

    # Hinge axis (defined in the local, undeformed beam basis)
    localHingeAxis = rotation_tensor_E321([-flareAngle; 0; 0]) * a2

    # Set the reference DOF as the greatest component of hinge axis
    refDOF = argmax(abs.(localHingeAxis))

    # Guess value for fold (applicable only when fold angle is unknown)
    foldGuessValue = foldAngleIsInput ? nothing : 0

    # Set hinge axis constraint
    hingeAxisConstraint = create_HingeAxisConstraint(beam=beam,masterElementLocalID=inboardElem,slaveElementLocalID=outboardElem,localHingeAxis=localHingeAxis,loadBalanceLocalNode=hingeNode+1,foldGuessValue=foldGuessValue)

    # Set fold angle constraint, if applicable
    foldAngleConstraint = foldAngleIsInput ? create_RotationConstraint(beam=beam,masterElementLocalID=inboardElem,slaveElementLocalID=outboardElem,masterDOF=refDOF,slaveDOF=refDOF,value=4*tan(foldAngle/4)*hingeAxisConstraint.initialHingeAxis[refDOF],loadBalanceLocalNode=hingeNode+1) : Vector{RotationConstraint}()
    
    rotationConstraints = foldAngleIsInput ? [foldAngleConstraint] : Vector{RotationConstraint}()

    # OOP bending spring around hinge
    if kSpring > 0
        spring = create_Spring(basis="A",elementsIDs=[inboardElem,outboardElem],nodesSides=[1,2],kTwist=kSpring,kIPBending=kSpring,kOOPBending=kSpring)
        add_spring_to_beams!(beams=[beam,beam],spring=spring)
    end

    # BCs
    clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

    # Wing model
    pazyFFWT = create_Model(name="pazyFFWT",beams=[beam],BCs=[clamp],gravityVector=[0;0;-g],v_A=[0;airspeed;0],rotationConstraints=rotationConstraints,hingeAxisConstraints=[hingeAxisConstraint],units=create_UnitsSystem(frequency="Hz"))

    return pazyFFWT
end
export create_PazyFFWT


"""
    create_SMW(; kwargs...)

Creates the 16 meters wing (SMW) of the conventional HALE aircraft described by Patil, Hodges and Cesnik in: Nonlinear Aeroelasticity and Flight Dynamics of HALE (2001)

Returns the wing model and its span

# Arguments
- `aeroSolver::AeroSolver` = aerodynamic solver
- `flapLoadsSolver::FlapAeroSolver` = aerodynamic solver for flap loads
- `gustLoadsSolver::GustAeroSolver` = indicial gust loads solver
- `derivationMethod::DerivationMethod` = method for aerodynamic derivatives
- `airfoil::Airfoil` = airfoil section
- `θ::Real` = pitch angle [rad]
- `k1::Real` = twisting curvature
- `k2::Real` = flapwise bending curvature
- `nElem::Int64` = number of elements for discretization
- `altitude::Real` = altitude
- `airspeed::Real` = airspeed
- `g::Real` = acceleration of gravity
- `stiffnessFactor::Real` = stiffness factor for beam structural properties
- `∞::Real` = value of rigid structural properties
- `tipF3::Real` = tip dead transverse force applied at the tip
- `cd0::Real` = parasite drag coefficient for the wing
- `cnδ::Real` = cn vs δ slope for the wing
- `cmδ::Real` = cm vs δ slope for the wing
- `cdδ::Real` = cd vs δ slope for the wing
- `hasInducedDrag::Bool` = flag to include induced drag
- `δAil::Union{Nothing,Real,<:Function}` = aileron deflection
- `additionalBCs::Vector{BC}` = additional BCs (beyond the clamp and tip load)
"""
function create_SMW(; aeroSolver::AeroSolver=Indicial(),flapLoadsSolver::FlapAeroSolver=TableLookup(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),derivationMethod::DerivationMethod=AD(),airfoil::Airfoil=deepcopy(flatPlate),θ::Real=0,k1::Real=0,k2::Real=0,nElem::Int64=32,altitude::Real=0,airspeed::Real=0,g::Real=standard_atmosphere(altitude).g,stiffnessFactor::Real=1,∞::Real=1e12,tipF3::Real=0,cd0::Real=0,cnδ::Real=2.5,cmδ::Real=-0.35,cdδ::Real=0.15,hasInducedDrag::Bool=false,δAil::Union{Nothing,Real,<:Function}=nothing,additionalBCs::Vector{BC}=Vector{BC}())

    # Validate
    @assert -π/2 < θ < π/2
    @assert nElem > 1
    @assert altitude >= 0
    @assert airspeed >= 0
    if δAil isa Real
        δAilConst = deepcopy(δAil)
        δAil = t -> δAilConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Wing surface
    chord = 1
    normSparPos = 0.5
    normFlapPos = 0.75
    ailSize = 0.25
    normFlapSpan = [1-ailSize,1]
    update_Airfoil_params!(airfoil,Ma=airspeed/atmosphere.a,U=airspeed,b=chord/2)

    wingSurf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=airfoil,c=chord,normSparPos=normSparPos,normFlapPos=normFlapPos,normFlapSpan=normFlapSpan,δ=δAil,updateAirfoilParameters=false,hasInducedDrag=hasInducedDrag,hasSymmetricCounterpart=true)

    # Override wing airfoil parameters
    wingSurf.airfoil.attachedFlowParameters.cd₀ = cd0
    wingSurf.airfoil.parametersBLi.cd₀ = cd0
    wingSurf.airfoil.parametersBLo.cd₀ = cd0
    wingSurf.airfoil.flapParameters.cnδ = cnδ
    wingSurf.airfoil.flapParameters.cmδ = cmδ
    wingSurf.airfoil.flapParameters.cdδ = cdδ

    # Wing properties
    L = 16
    GJ,EIy,EIz = 1e4,2e4,4e6
    GJ,EIy,EIz = multiply_inplace!(stiffnessFactor, GJ,EIy,EIz)
    ρA,ρIs = 0.75,0.1
    ρIy,ρIz = (EIy/EIz)*ρIs,(1-EIy/EIz)*ρIs
    S = isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)
    I = inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs)

    # Wing beam
    beam = create_Beam(name="rightWing",length=L,nElements=nElem,S=[S],I=[I],aeroSurface=wingSurf,rotationParametrization="E321",p0=[0;0;θ],k=[k1;k2;0])
    
    # Update beam of additional BCs, if applicable
    for BC in additionalBCs
        BC.beam = beam
    end
    BCs = additionalBCs

    # Add root clamp and tip load to BCs
    clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])
    push!(BCs,clamp)

    # Add tip load to BCs, if applicable
    if isempty(additionalBCs)
        tipLoad = create_BC(name="tipLoad",beam=beam,node=nElem+1,types=["F3A"],values=[tipF3])
        push!(BCs,tipLoad)
    end

    # Wing model
    wing = create_Model(name="SMW",beams=[beam],BCs=BCs,altitude=altitude,gravityVector=[0;0;-g],v_A=[0;airspeed;0])

    return wing,L
end
export create_SMW


"""
    create_TDWing(; kwargs...)

Creates the wing described by Tang and Dowell in: Experimental and Theoretical Study on Aeroelastic Response of High-Aspect-Ratio Wings (2001)

Returns the wing model

# Arguments
- `aeroSolver::AeroSolver` = aerodynamic solver
- `gustLoadsSolver::GustAeroSolver` = indicial gust loads solver
- `derivationMethod::DerivationMethod` = method for aerodynamic derivatives
- `updateAirfoilParameters::Bool` = flag to update airfoil parameters with airspeed
- `airfoil::Airfoil` = airfoil section
- `θ::Real` = pitch angle
- `nElem::Int64` = number of elements for discretization
- `altitude::Real` = altitude
- `airspeed::Real` = airspeed
- `g::Real` = acceleration of gravity
- `∞::Real` = value of rigid structural properties
"""
function create_TDWing(; aeroSolver::AeroSolver=Indicial(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),derivationMethod::DerivationMethod=AD(),updateAirfoilParameters::Bool=true,airfoil::Airfoil=deepcopy(flatPlate),θ::Real=0,nElem::Int64=20,altitude::Real=0,airspeed::Real=0,g::Real=9.80665,∞::Real=1e6)

    # Validate
    @assert -π/2 <= θ <= π/2
    @assert nElem > 1
    @assert altitude >= 0
    @assert airspeed >= 0

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Wing surface
    chord = 0.0508
    normSparPos = 0.5
    update_Airfoil_params!(airfoil,Ma=airspeed/atmosphere.a,U=airspeed,b=chord/2)
    wingSurf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=chord,normSparPos=normSparPos,updateAirfoilParameters=updateAirfoilParameters)

    # Wing properties
    L = 0.4508
    GJ,EIy,EIz = 0.9539,0.4186,18.44
    ρA,ρIs,ρIy = 0.2351,0.2056e-4,1e-6
    ρIz = ρIy*EIz/EIy
    e3 = 1e-2*chord
    S = isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)
    I = inertia_matrix(ρA=ρA,ρIy=ρIy,ρIz=ρIz,ρIs=ρIs,e3=e3)

    # Wing beam
    beam = create_Beam(name="wingBeam",length=L,nElements=nElem,S=[S],I=[I],aeroSurface=wingSurf,rotationParametrization="E321",p0=[0;0;θ])

    # Wing's tip store
    tipMass = 0.0417
    tipIyy = 0.3783e-5
    tipIzz = 0.9753e-4
    tipStore = PointInertia(elementID=nElem,η=[L/nElem/2;0;0],mass=tipMass,Iyy=tipIyy,Izz=tipIzz)
    add_point_inertias_to_beam!(beam,inertias=[tipStore])

    # BCs
    clamp = create_BC(name="clamp",beam=beam,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

    # Wing model
    wing = create_Model(name="TDWing",beams=[beam],BCs=[clamp],gravityVector=[0;0;-g],v_A=[0;airspeed;0])

    return wing
end
export create_TDWing


"""
    create_Helios(; kwargs...)

Creates a model based on the flying-wing aircraft described by Patil and Hodges in: Flight Dynamics of Highly Flexible Flying Wings (2006)

# Keyword arguments
- `altitude::Real` = altitude
- `aeroSolver::AeroSolver` = aerodynamic solver
- `gustLoadsSolver::GustAeroSolver` = indicial gust loads solver
- `derivationMethod::DerivationMethod` = method for aerodynamic derivatives
- `g::Real` = local acceleration of gravity
- `wingAirfoil::Airfoil` = airfoil section of the wing
- `podAirfoil::Airfoil` = airfoil section of the pods
- `beamPods::Bool` = flag to include pods
- `stiffnessFactor::Real` = stiffness factor for the wing structure
- `∞::Real` = value for rigid structural properties
- `nElemStraightSemispan::Int64` = number of elements in the straight section of the semispan
- `nElemDihedralSemispan::Int64` = number of elements in the dihedral section of the semispan
- `nElemPod::Int64` = number of elements in the pods
- `payloadPounds::Real` = payload, in pounds
- `airspeed::Real` = local initial/trim airspeed
- `δIsTrimVariable::Bool` = flag for flap deflection being a trim variable
- `thrustIsTrimVariable::Bool` = flag for motors' thrust being a trim variable
- `δ::Union{Nothing,Real,<:Function}` = flap deflection
- `thrust::Union{Real,<:Function}` = motors' thrust
- `reducedChord::Bool` = flag to employ a reduced (7 ft) chord
- `payloadOnWing::Bool` = flag to set the payload on the wing's reference line
"""
function create_Helios(; altitude::Real=0,aeroSolver::AeroSolver=Indicial(),derivationMethod::DerivationMethod=AD(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),g::Real=-9.80665,wingAirfoil::Airfoil=deepcopy(HeliosWingAirfoil),podAirfoil::Airfoil=HeliosPodAirfoil,beamPods::Bool=false,stiffnessFactor::Real=1.0,∞::Real=1e12,nElemStraightSemispan::Int64=10,nElemDihedralSemispan::Int64=5,nElemPod::Int64=1,payloadPounds::Real=0,airspeed::Real=0,δIsTrimVariable::Bool=false,thrustIsTrimVariable::Bool=false,δ::Union{Nothing,Real,<:Function}=nothing,thrust::Union{Real,<:Function}=0,reducedChord::Bool=false,payloadOnWing::Bool=false)

    # Validate
    @assert ∞ > 1e8
    @assert stiffnessFactor > 0
    @assert payloadPounds >= 0
    @assert airspeed >= 0
    @assert iseven(nElemStraightSemispan)
    δIsInput = !isnothing(δ)
    if δIsInput
        @assert !δIsTrimVariable
    end
    if δ isa Real
        δconst = deepcopy(δ)
        δ = t -> δconst
    end
    if thrust isa Real
        thrustConst = deepcopy(thrust)
        thrust = t -> thrustConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Wing surface
    wChord = reducedChord ? 7*0.3048 : 8*0.3048
    wNormSparPos = 0.25
    wNormFlapPos = 0.75
    wNormFlapSpan = [0,1]
    update_Airfoil_params!(wingAirfoil,Ma=airspeed/atmosphere.a,U=airspeed,b=wChord/2)
    wingSurf = δIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=wingAirfoil,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=wNormFlapSpan,δ=δ,updateAirfoilParameters=false) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=wingAirfoil,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=wNormFlapSpan,δIsTrimVariable=δIsTrimVariable,updateAirfoilParameters=false)

    # Wing properties
    wingSection = 40*0.3048
    Γ = 10*π/180
    GJ,EIy,EIz = 1.6530e5,1.0331e6,1.2398e7
    wρA,wρIy,wρIz = 8.929,0.691,3.456
    GJ,EIy,EIz = multiply_inplace!(stiffnessFactor, GJ,EIy,EIz)
    Swing = isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)
    Iwing = inertia_matrix(ρA=wρA,ρIy=wρIy,ρIz=wρIz)

    # Wing beams
    leftWingDihedral = create_Beam(name="leftWingDihedral",length=wingSection,nElements=nElemDihedralSemispan,S=[Swing],I=[Iwing],aeroSurface=deepcopy(wingSurf),rotationParametrization="E321",p0=[0;Γ;0])

    leftWingStraight = create_Beam(name="leftWingStraight",length=2*wingSection,nElements=nElemStraightSemispan,S=[Swing],I=[Iwing],aeroSurface=deepcopy(wingSurf))

    rightWingStraight = create_Beam(name="rightWingStraight",length=2*wingSection,nElements=nElemStraightSemispan,S=[Swing],I=[Iwing],aeroSurface=deepcopy(wingSurf))

    rightWingDihedral = create_Beam(name="rightWingDihedral",length=wingSection,nElements=nElemDihedralSemispan,S=[Swing],I=[Iwing],aeroSurface=deepcopy(wingSurf),rotationParametrization="E321",p0=[0;-Γ;0])

    # Link elevators
    elevatorLink = create_FlapLink(masterBeam=rightWingStraight,slaveBeams=[leftWingDihedral,leftWingStraight,rightWingDihedral])

    # Pods' point inertias
    podLength = 6*0.3048
    podMassVertOffset = payloadOnWing ? 0 : -podLength/2
    sidePodMass = 50*0.453592
    centerPodMass = 60*0.453592
    payloadMass = payloadPounds*0.453592

    leftSidePod = PointInertia(elementID=1,η=[-(2*wingSection/nElemStraightSemispan)/2;0;podMassVertOffset],mass=sidePodMass)

    centerPod = PointInertia(elementID=nElemStraightSemispan,η=[(2*wingSection/nElemStraightSemispan)/2;0;podMassVertOffset],mass=centerPodMass+payloadMass)

    rightSidePod = PointInertia(elementID=nElemStraightSemispan,η=[(2*wingSection/nElemStraightSemispan)/2;0;podMassVertOffset],mass=sidePodMass)

    add_point_inertias_to_beam!(leftWingStraight,inertias=[leftSidePod,centerPod])

    add_point_inertias_to_beam!(rightWingStraight,inertias=[rightSidePod])

    # Pod surface
    update_Airfoil_params!(podAirfoil,Ma=airspeed/atmosphere.a,U=airspeed,b=wChord/2)
    podSurf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=podAirfoil,c=wChord,normSparPos=wNormSparPos,updateAirfoilParameters=false)

    # Pod properties
    pρA,pρIy,pρIz = 0,wρIy,wρIz
    Spod = isotropic_stiffness_matrix(∞=∞)
    Ipod = inertia_matrix(ρA=pρA,ρIy=pρIy,ρIz=pρIz)

    # Pod beams
    leftPod = create_Beam(name="leftPod",length=podLength,nElements=nElemPod,S=[Spod],I=[Ipod],aeroSurface=deepcopy(podSurf),rotationParametrization="E321",p0=[0;π/2;0],connectedBeams=[leftWingStraight],connectedNodesThis=[1],connectedNodesOther=[1])

    centerPod = create_Beam(name="centerPod",length=podLength,nElements=nElemPod,S=[Spod],I=[Ipod],aeroSurface=deepcopy(podSurf),rotationParametrization="E321",p0=[0;π/2;0],connectedBeams=[rightWingStraight],connectedNodesThis=[1],connectedNodesOther=[1])

    rightPod = create_Beam(name="rightPod",length=podLength,nElements=nElemPod,S=[Spod],I=[Ipod],aeroSurface=deepcopy(podSurf),rotationParametrization="E321",p0=[0;π/2;0],connectedBeams=[rightWingStraight],connectedNodesThis=[1],connectedNodesOther=[nElemStraightSemispan+1])

    # Generate copies for wing model (has to be done before aircraft model creation)
    wingStraight = deepcopy(rightWingStraight)
    wingDihedral = deepcopy(rightWingDihedral)
    wingPod = deepcopy(rightPod)

    # Propellers thrust force
    thrustValue = thrustIsTrimVariable ? t -> 0 : thrust

    thrustLeftOut = create_BC(name="thrustLeftOut",beam=leftWingStraight,node=1,types=["Ff2b"],values=[t->thrustValue(t)],toBeTrimmed=[thrustIsTrimVariable])

    thrustLeftIn = create_BC(name="thrustLeftIn",beam=leftWingStraight,node=div(nElemStraightSemispan,2)+1,types=["Ff2b"],values=[t->thrustValue(t)],toBeTrimmed=[thrustIsTrimVariable])

    thrustCenter = create_BC(name="thrustCenter",beam=leftWingStraight,node=nElemStraightSemispan+1,types=["Ff2b"],values=[t->thrustValue(t)],toBeTrimmed=[thrustIsTrimVariable])

    thrustRightIn = create_BC(name="thrustRightIn",beam=rightWingStraight,node=div(nElemStraightSemispan,2)+1,types=["Ff2b"],values=[t->thrustValue(t)],toBeTrimmed=[thrustIsTrimVariable])

    thrustRightOut = create_BC(name="thrustRightOut",beam=rightWingStraight,node=nElemStraightSemispan+1,types=["Ff2b"],values=[t->thrustValue(t)],toBeTrimmed=[thrustIsTrimVariable])

    aircraftBCs = [thrustLeftOut,thrustLeftIn,thrustCenter,thrustRightIn,thrustRightOut]

    # Link propellers' thrust, if applicable
    thrustsLink = thrustIsTrimVariable ? [create_TrimLoadsLink(masterBC=thrustCenter,slaveBCs=[thrustLeftOut,thrustLeftIn,thrustRightIn,thrustRightOut])] : Vector{TrimLoadsLink}()

    # Beams of aircraft model
    aircraftBeams = beamPods ? [leftWingDihedral,leftWingStraight,rightWingStraight,rightWingDihedral,leftPod,centerPod,rightPod] : [leftWingDihedral,leftWingStraight,rightWingStraight,rightWingDihedral]

    # Aircraft model (with initial position such that aircraft center is coincident with the origin of frame A)
    helios = create_Model(name="Helios",beams=aircraftBeams,BCs=aircraftBCs,initialPosition=[-wingSection*(2+cos(Γ));0;wingSection*sin(Γ)],altitude=altitude,gravityVector=[0;0;g],v_A=[0;airspeed;0],flapLinks=[elevatorLink],trimLoadsLinks=thrustsLink)

    # Set midpsan element for the aircraft model
    midSpanElem = nElemDihedralSemispan + nElemStraightSemispan

    # Beams of wing model
    wingBeams = beamPods ? [wingStraight,wingDihedral,wingPod] : [wingStraight,wingDihedral]

    # Clamp for wing model
    clamp = create_BC(name="clamp",beam=wingStraight,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

    # Link elevators
    elevatorLinkWing = create_FlapLink(masterBeam=wingStraight,slaveBeams=[wingDihedral])

    # Wing model
    wing = create_Model(name="HeliosWing",beams=wingBeams,BCs=[clamp],altitude=altitude,gravityVector=[0;0;g],v_A=[0;airspeed;0],flapLinks=[elevatorLinkWing])

    return helios,midSpanElem,wing,leftWingStraight,rightWingStraight,leftWingDihedral,rightWingDihedral,leftPod,rightPod,centerPod

end
export create_Helios


"""
    create_conventional_HALE(; kwargs...)

Creates a model based on the conventional HALE aircraft described by Patil, Hodges and Cesnik in: Nonlinear Aeroelasticity and Flight Dynamics of HALE (2001)

# Keyword arguments
- `altitude::Real` = altitude
- `aeroSolver::AeroSolver` = aerodynamic solver
- `derivationMethod::DerivationMethod` = method for aerodynamic derivatives
- `flapLoadsSolver::FlapAeroSolver` = aerodynamic solver for flap loads
- `gustLoadsSolver::GustAeroSolver` = indicial gust loads solver
- `stabilizersAero::Bool` = flag for stabilizers with aerodynamic surfaces
- `includeVS::Bool` = flag to include a vertical stabilizer in the model
- `wAirfoil::Airfoil` = wing airfoil section
- `sAirfoil::Airfoil` = stabilizers airfoil section
- `nElemWing::Int64` = number of elements of the full wing
- `nElemHorzStabilizer::Int64` = number of elements of the horizontal stabilizer
- `nElemTailBoom::Int64` = number of elements of the tail boom
- `nElemVertStabilizer::Int64` = number of elements of the vertical stabilizer
- `∞::Real=1e12` = value of rigid structural properties
- `stiffnessFactor::Real` = stiffness factor for the wing structure
- `k1::Real` = undeformed wing torsional curvature
- `k2::Real` = undeformed wing flapwise bending curvature
- `airspeed::Real` = local initial/trim airspeed
- `δElevIsTrimVariable::Bool` = flag for elevator deflection being a trim variable
- `δAilIsTrimVariable::Bool` = flag for aileron deflection being a trim variable
- `δRuddIsTrimVariable::Bool` = flag for rudder deflection being a trim variable    
- `thrustIsTrimVariable::Bool` = flag for motors' thrust being a trim variable
- `δElev::Union{Nothing,Real,<:Function}` = elevator deflection
- `δAil::Union{Nothing,Real,<:Function}` = aileron deflection
- `δRudd::Union{Nothing,Real,<:Function}` = rudder deflection
- `thrust::Union{Real,<:Function}` = motors' thrust
- `g::Real` = local acceleration of gravity
- `wingCd0::Real` = parasite drag coefficient for the wing
- `wingcnδ::Real` = cn vs δ slope for the wing
- `wingcmδ::Real` = cm vs δ slope for the wing
- `wingcdδ::Real` = cd vs δ slope for the wing
- `stabsCd0::Real` = parasite drag coefficient for the stabilizers
- `stabscnδ::Real` = cn vs δ slope for the stabilizers
- `stabscmδ::Real` = cm vs δ slope for the stabilizers
- `stabscdδ::Real` = cd vs δ slope for the stabilizers
- `hasInducedDrag::Bool` = flag to include induced drag
"""
function create_conventional_HALE(; altitude::Real=20e3,aeroSolver::AeroSolver=Indicial(),derivationMethod::DerivationMethod=AD(),flapLoadsSolver::FlapAeroSolver=TableLookup(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),stabilizersAero::Bool=true,includeVS::Bool=true,wAirfoil::Airfoil=deepcopy(flatPlate),sAirfoil::Airfoil=deepcopy(flatPlate),nElemWing::Int64=20,nElemHorzStabilizer::Int64=10,nElemTailBoom::Int64=10,nElemVertStabilizer::Int64=5,∞::Real=1e12,stiffnessFactor::Real=1,k1::Real=0,k2::Real=0,airspeed::Real=0,δElevIsTrimVariable::Bool=false,δAilIsTrimVariable::Bool=false,δRuddIsTrimVariable::Bool=false,thrustIsTrimVariable::Bool=false,δElev::Union{Nothing,Real,<:Function}=nothing,δAil::Union{Nothing,Real,<:Function}=nothing,δRudd::Union{Nothing,Real,<:Function}=nothing,thrust::Union{Real,<:Function}=0,g::Real=-9.80665,wingCd0::Real=0,wingcnδ::Real=2.5,wingcmδ::Real=-0.35,wingcdδ::Real=0.15,stabsCd0::Real=0,stabscnδ::Real=2.5,stabscmδ::Real=-0.35,stabscdδ::Real=0.15,hasInducedDrag::Bool=false)

    # Validate
    @assert iseven(nElemWing)
    @assert iseven(nElemHorzStabilizer)
    @assert ∞ > 1e8
    @assert stiffnessFactor > 0
    @assert wingCd0 >= 0
    @assert stabsCd0 >= 0
    @assert airspeed >= 0
    δElevIsInput = !isnothing(δElev)
    δAilIsInput = !isnothing(δAil)
    δRuddIsInput = !isnothing(δRudd)
    if δElevIsInput
        @assert !δElevIsTrimVariable
    end
    if δElev isa Real
        δElevConst = deepcopy(δElev)
        δElev = t -> δElevConst
    end
    if δAilIsInput
        @assert !δAilIsTrimVariable
    end
    if δAil isa Real
        δAilConst = deepcopy(δAil)
        δAil = t -> δAilConst
    end
    if δRuddIsInput
        @assert !δRuddIsTrimVariable
    end
    if δRudd isa Real
        δRuddConst = deepcopy(δRudd)
        δRudd = t -> δRuddConst
    end
    if thrust isa Real
        thrustConst = deepcopy(thrust)
        thrust = t -> thrustConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Wing surface
    wChord = 1
    wNormSparPos = 0.5
    wNormFlapPos = 0.75
    wFlapSize = 0.25
    wLNormFlapSpan = [0,wFlapSize]
    wRNormFlapSpan = [1-wFlapSize,1]
    update_Airfoil_params!(wAirfoil,Ma=airspeed/atmosphere.a,U=airspeed,b=wChord/2)

    wingSurfLeft = δAilIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=wAirfoil,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=wLNormFlapSpan,δ=δAil,updateAirfoilParameters=false,hasInducedDrag=hasInducedDrag,hasSymmetricCounterpart=true) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=wAirfoil,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=wLNormFlapSpan,δIsTrimVariable=δAilIsTrimVariable,updateAirfoilParameters=false,hasInducedDrag=hasInducedDrag,hasSymmetricCounterpart=true)

    wingSurfRight = δAilIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=wAirfoil,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=wRNormFlapSpan,δ=δAil,updateAirfoilParameters=false,hasInducedDrag=hasInducedDrag,hasSymmetricCounterpart=true) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=wAirfoil,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=wRNormFlapSpan,δIsTrimVariable=δAilIsTrimVariable,updateAirfoilParameters=false,hasInducedDrag=hasInducedDrag,hasSymmetricCounterpart=true)

    # Override wing airfoil parameters
    wingSurfLeft.airfoil.attachedFlowParameters.cd₀ = wingCd0
    wingSurfLeft.airfoil.parametersBLi.cd₀ = wingCd0
    wingSurfLeft.airfoil.parametersBLo.cd₀ = wingCd0
    wingSurfLeft.airfoil.flapParameters.cnδ = wingcnδ
    wingSurfLeft.airfoil.flapParameters.cmδ = wingcmδ
    wingSurfLeft.airfoil.flapParameters.cdδ = wingcdδ

    wingSurfRight.airfoil.attachedFlowParameters.cd₀ = wingCd0
    wingSurfRight.airfoil.parametersBLi.cd₀ = wingCd0
    wingSurfRight.airfoil.parametersBLo.cd₀ = wingCd0
    wingSurfRight.airfoil.flapParameters.cnδ = wingcnδ
    wingSurfRight.airfoil.flapParameters.cmδ = wingcmδ
    wingSurfRight.airfoil.flapParameters.cdδ = wingcdδ

    # Wing properties
    Lw = 16
    wGJ,wEIy,wEIz = 1e4,2e4,4e6
    wGJ,wEIy,wEIz = multiply_inplace!(stiffnessFactor, wGJ,wEIy,wEIz)
    wρA,wρIs = 0.75,0.1
    wρIy,wρIz = (wEIy/wEIz)*wρIs,(1-wEIy/wEIz)*wρIs
    Swing = isotropic_stiffness_matrix(∞=∞,GJ=wGJ,EIy=wEIy,EIz=wEIz)
    Iwing = inertia_matrix(ρA=wρA,ρIy=wρIy,ρIz=wρIz,ρIs=wρIs)

    # Initial position for first node of left wing
    ρ = 1/k2
    θ = Lw/ρ
    x = ρ*sin(θ)
    z = ρ*(1-cos(θ))
    initialPosition = k2 == 0 ? [-Lw; 0; 0] : [-x; 0; -z]

    # Initial angle of twist
    r = 1/k1
    ψ = Lw/r

    # Wing beams
    leftWing = create_Beam(name="leftWing",length=Lw,nElements=div(nElemWing,2),S=[Swing],I=[Iwing],aeroSurface=wingSurfLeft,k=[-k1;k2;0],rotationParametrization="E321",p0=[0;-θ;ψ])

    rightWing = create_Beam(name="rightWing",length=Lw,nElements=div(nElemWing,2),S=[Swing],I=[Iwing],aeroSurface=wingSurfRight,k=[k1;k2;0],connectedBeams=[leftWing],connectedNodesThis=[1],connectedNodesOther=[div(nElemWing,2)+1])

    # Link wing ailerons
    aileronLink = create_FlapLink(masterBeam=rightWing,slaveBeams=[leftWing],δMultipliers=[-1])

    # Payload
    payloadMass = 50
    payloadInertia = 200
    payload = PointInertia(elementID=1,η=[-Lw/div(nElemWing,2)/2;0;0],mass=payloadMass,Iyy=payloadInertia,Izz=payloadInertia,Ixx=payloadInertia)
    add_point_inertias_to_beam!(rightWing,inertias=[payload])

    # Tail boom
    Lt = 10
    tρA,tρIy,tρIz = 0.08,wρIy/10,wρIz/10
    tailBoom = create_Beam(name="tailBoom",length=Lt,nElements=nElemTailBoom,S=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=tρA,ρIy=tρIy,ρIz=tρIz)],rotationParametrization="E321",p0=[-π/2;0;0],connectedBeams=[rightWing],connectedNodesThis=[1],connectedNodesOther=[1])

    # Horizontal stabilizer surface
    hChord = 0.5
    hNormSparPos = 0.5
    hNormFlapPos = 0.75
    hNormFlapSpan = [0,1]
    update_Airfoil_params!(sAirfoil,Ma=airspeed/atmosphere.a,U=airspeed,b=hChord/2)
    hsSurf = δElevIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=sAirfoil,c=hChord,normSparPos=hNormSparPos,normFlapPos=hNormFlapPos,normFlapSpan=hNormFlapSpan,δ=δElev,updateAirfoilParameters=false,hasInducedDrag=hasInducedDrag,hasSymmetricCounterpart=false) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=sAirfoil,c=hChord,normSparPos=hNormSparPos,normFlapPos=hNormFlapPos,normFlapSpan=hNormFlapSpan,δIsTrimVariable=δElevIsTrimVariable,updateAirfoilParameters=false,hasInducedDrag=hasInducedDrag,hasSymmetricCounterpart=false)

    # Override horizontal stabilizer airfoil parameters
    hsSurf.airfoil.attachedFlowParameters.cd₀ = stabsCd0
    hsSurf.airfoil.parametersBLi.cd₀ = stabsCd0
    hsSurf.airfoil.parametersBLo.cd₀ = stabsCd0
    hsSurf.airfoil.flapParameters.cnδ = stabscnδ
    hsSurf.airfoil.flapParameters.cmδ = stabscmδ
    hsSurf.airfoil.flapParameters.cdδ = stabscdδ

    # Horizontal stabilizer beam
    Lh = 5
    hρA,hρIy,hρIz = 0.08,wρIy/10,wρIz/10
    horzStabilizer = create_Beam(name="horzStabilizer",length=Lh,initialPosition=[-Lh/2;0;0],nElements=nElemHorzStabilizer,S=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=hρA,ρIy=hρIy,ρIz=hρIz)],connectedBeams=[tailBoom],connectedNodesThis=[div(nElemHorzStabilizer,2)+1],connectedNodesOther=[nElemTailBoom+1])
    if stabilizersAero
        horzStabilizer.aeroSurface = hsSurf
        update_beam!(horzStabilizer)
    end

    # Vertical stabilizer surface
    vChord = 0.5
    vNormSparPos = 0.5
    vNormFlapPos = 0.75
    vNormFlapSpan = [0,1]
    update_Airfoil_params!(sAirfoil,Ma=airspeed/atmosphere.a,U=airspeed,b=vChord/2)
    vsSurf = δRuddIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=sAirfoil,c=vChord,normSparPos=vNormSparPos,normFlapPos=vNormFlapPos,normFlapSpan=vNormFlapSpan,δ=δRudd,updateAirfoilParameters=false,hasInducedDrag=hasInducedDrag,hasSymmetricCounterpart=false) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=sAirfoil,c=vChord,normSparPos=vNormSparPos,normFlapPos=vNormFlapPos,normFlapSpan=vNormFlapSpan,δIsTrimVariable=δRuddIsTrimVariable,updateAirfoilParameters=false,hasInducedDrag=hasInducedDrag,hasSymmetricCounterpart=false)

    # Override vertical stabilizer airfoil parameters
    vsSurf.airfoil.attachedFlowParameters.cd₀ = stabsCd0
    vsSurf.airfoil.parametersBLi.cd₀ = stabsCd0
    vsSurf.airfoil.parametersBLo.cd₀ = stabsCd0
    vsSurf.airfoil.flapParameters.cnδ = stabscnδ
    vsSurf.airfoil.flapParameters.cmδ = stabscmδ
    vsSurf.airfoil.flapParameters.cdδ = stabscdδ

    # Vertical stabilizer beam
    Lv = 2.5
    vρA,vρIy,vρIz = 0.08,hρIy,hρIz
    vertStabilizer = create_Beam(name="vertStabilizer",length=Lv,nElements=nElemVertStabilizer,S=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=vρA,ρIy=vρIy,ρIz=vρIz)],rotationParametrization="E321",p0=[0;-π/2;0],connectedBeams=[tailBoom],connectedNodesThis=[1],connectedNodesOther=[nElemTailBoom+1])
    if stabilizersAero
        vertStabilizer.aeroSurface = vsSurf
        update_beam!(vertStabilizer)
    end

    # Propeller thrust force 
    thrustValue = thrustIsTrimVariable ? t -> 0 : thrust

    propThrust = create_BC(name="propThrust",beam=rightWing,node=1,types=["Ff2b"],values=[t->thrustValue(t)],toBeTrimmed=[thrustIsTrimVariable])

    # Beams of the model
    beams = includeVS ? [leftWing,rightWing,tailBoom,horzStabilizer,vertStabilizer] : [leftWing,rightWing,tailBoom,horzStabilizer]

    # Aircraft model
    conventionalHALE = create_Model(name="conventionalHALE",beams=beams,BCs=[propThrust],initialPosition=initialPosition,v_A=[0;airspeed;0],altitude=altitude,gravityVector=[0;0;g],flapLinks=[aileronLink])

    return conventionalHALE,leftWing,rightWing,tailBoom,horzStabilizer,vertStabilizer
end
export create_conventional_HALE


"""
    create_BWB(; kwargs...)

Creates a model based on the blended-wing-body described by Weihua Su's PhD thesis

# Keyword arguments
- `altitude::Real` = altitude
- `aeroSolver::AeroSolver` = aerodynamic solver
- `gustLoadsSolver::GustAeroSolver` = indicial gust loads solver
- `derivationMethod::DerivationMethod` = method for aerodynamic derivatives 
- `∞::Real=1e12` = value of rigid structural properties
- `stiffnessFactor::Real` = stiffness factor for the wing structure
- airspeed::Real = local initial/trim airspeed
- `δElevIsTrimVariable::Bool` = flag for elevator deflection being a trim variable
- `thrustIsTrimVariable::Bool` = flag for motors' thrust being a trim variable
- `δElev::Union{Nothing,Real,<:Function}` = elevator deflection
- `thrust::Union{Real,<:Function}` = motors' thrust
- `g::Real` = local acceleration of gravity
- `updateAirfoilParameters::Bool` = flag to update airfoil parameters with airspeed
- `hasTipCorrection::Bool` = flag to employ aerodynamic tip correction
- `tipLossDecayFactor::Real` = tip loss decay factor
"""
function create_BWB(; altitude::Real=0,aeroSolver::AeroSolver=Indicial(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),derivationMethod::DerivationMethod=AD(),∞::Real=1e12,stiffnessFactor::Real=1,airspeed::Real=0,δElevIsTrimVariable::Bool=false,thrustIsTrimVariable::Bool=false,δElev::Union{Nothing,Real,<:Function}=nothing,thrust::Union{Real,<:Function}=0,g::Real=-9.80665,updateAirfoilParameters::Bool=false,hasTipCorrection::Bool=false,tipLossDecayFactor::Real=40)

    # Validate
    @assert ∞ > 1e8
    @assert stiffnessFactor > 0
    @assert airspeed >= 0
    δElevIsInput = !isnothing(δElev)
    if δElevIsInput
        @assert !δElevIsTrimVariable
    end
    if δElev isa Real
        δElevConst = deepcopy(δElev)
        δElev = t -> δElevConst
    end
    if thrust isa Real
        thrustConst = deepcopy(thrust)
        thrust = t -> thrustConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Keypoints
    kp1 = [0; -1.126236+1.126236; 0]
    kp2 = [0.889; -0.997712+1.126236; 0]
    kp3 = [3.248152; -2.359914+1.126236; 0]

    # Length of fuselage section
    fusLength = sqrt((kp2[1]-kp1[1])^2+(kp2[2]-kp1[2])^2)

    # Angle of sweep of the fuselage section 
    fΛ = acos((kp2[1]-kp1[1])/fusLength)             

    # Length of wing section
    wingSemispan = sqrt((kp3[1]-kp2[1])^2+(kp3[2]-kp2[2])^2)  
    
    # Angle of sweep of the wing section 
    wΛ = acos((kp3[1]-kp2[1])/wingSemispan)
    
    # Chord and spar position
    rootChord = 1.38557
    tipChord = 0.54864
    rootSparPos = 0.6438
    tipSparPos = 0.4560

    # Airfoil
    airfoil = deepcopy(BWBAirfoil)
    if updateAirfoilParameters
        update_Airfoil_params!(airfoil,Ma=airspeed/atmosphere.a,U=airspeed,b=tipChord/2)
    end

    # Wing surfaces
    wChord = tipChord
    wNormSparPos = tipSparPos
    wNormFlapPos = 0.75

    leftWingSurf = δElevIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=airfoil,Λ=wΛ,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=[1/4; 1],δ=δElev,updateAirfoilParameters=updateAirfoilParameters,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=-tipLossDecayFactor) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=airfoil,Λ=wΛ,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=[1/4; 1],δIsTrimVariable=δElevIsTrimVariable,updateAirfoilParameters=updateAirfoilParameters,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=-tipLossDecayFactor)

    rightWingSurf = δElevIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=airfoil,Λ=-wΛ,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=[0; 3/4],δ=δElev,updateAirfoilParameters=updateAirfoilParameters,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=airfoil,Λ=-wΛ,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=[0; 3/4],δIsTrimVariable=δElevIsTrimVariable,updateAirfoilParameters=updateAirfoilParameters,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor)

    # Wing properties
    nElemWing = 8
    wEA,wGJ,wEIy,wEIz = 155_000_000,11_000,11_700,130_000
    wGJ,wEIy,wEIz = multiply_inplace!(stiffnessFactor, wGJ,wEIy,wEIz)
    wρA,wρIy,wρIz = 6.2,0.0005,0.00462
    Swing = isotropic_stiffness_matrix(∞=∞,EA=wEA,GJ=wGJ,EIy=wEIy,EIz=wEIz)
    Iwing = inertia_matrix(ρA=wρA,ρIy=wρIy,ρIz=wρIz)

    # Concentrated wing inertias
    wConMass = 2
    wLConInertias = Vector{PointInertia}()
    wRConInertias = Vector{PointInertia}()
    push!(wLConInertias,PointInertia(elementID=1,η=[-wingSemispan/nElemWing/2;0;0],mass=wConMass))
    push!(wRConInertias,PointInertia(elementID=nElemWing,η=[wingSemispan/nElemWing/2;0;0],mass=wConMass))
    for e=1:8
        push!(wLConInertias,PointInertia(elementID=e,η=[wingSemispan/nElemWing/2;0;0],mass=wConMass))
        push!(wRConInertias,PointInertia(elementID=e,η=[-wingSemispan/nElemWing/2;0;0],mass=wConMass))
    end

    # Wing beams
    leftWing = create_Beam(name="leftWing",length=wingSemispan,nElements=nElemWing,S=[Swing],I=[Iwing],aeroSurface=leftWingSurf,rotationParametrization="E321",p0=[wΛ;0;0],pointInertias=wLConInertias)

    rightWing = create_Beam(name="rightWing",length=wingSemispan,nElements=nElemWing,S=[Swing],I=[Iwing],aeroSurface=rightWingSurf,rotationParametrization="E321",p0=[-wΛ;0;0],pointInertias=wRConInertias)

    # Link wing elevons
    elevonLink = create_FlapLink(masterBeam=rightWing,slaveBeams=[leftWing],δMultipliers=[1])

    # Fuselage surfaces
    leftFusChord = x1 -> tipChord + (rootChord-tipChord)*x1/fusLength
    rightFusChord = x1 -> rootChord + (tipChord-rootChord)*x1/fusLength
    leftFusSparPos = x1 -> tipSparPos + (rootSparPos-tipSparPos)*x1/fusLength
    rightFusSparPos = x1 -> rootSparPos + (tipSparPos-rootSparPos)*x1/fusLength

    leftFusSurf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=airfoil,Λ=-fΛ,c=leftFusChord,normSparPos=leftFusSparPos,updateAirfoilParameters=updateAirfoilParameters)

    rightFusSurf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=airfoil,Λ=fΛ,c=rightFusChord,normSparPos=rightFusSparPos,updateAirfoilParameters=updateAirfoilParameters)

    # Fuselage properties
    nElemFus = 3
    fEA,fGJ,fEIy,fEIz = 169_000_000,2_250_000,750_000,35_000_000
    fρA,fρIy,fρIz = 50,0.7,22.0
    Sfus = isotropic_stiffness_matrix(∞=∞,EA=fEA,GJ=fGJ,EIy=fEIy,EIz=fEIz)
    Ifus = inertia_matrix(ρA=fρA,ρIy=fρIy,ρIz=fρIz)

    # Concentrated fuselage inertias
    fConMass = 40
    fConMassOffset = 0.891968
    fLeftConInertia = PointInertia(elementID=nElemFus,η=[fusLength/nElemFus/2;fConMassOffset;0],mass=fConMass)
    fRightConInertia = PointInertia(elementID=1,η=[-fusLength/nElemFus/2;fConMassOffset;0],mass=fConMass)

    # Fuselage beams
    leftFus = create_Beam(name="leftFus",length=fusLength,nElements=nElemFus,S=[Sfus],I=[Ifus],aeroSurface=leftFusSurf,rotationParametrization="E321",p0=[-fΛ;0;0],pointInertias=[fLeftConInertia])

    rightFus = create_Beam(name="rightFus",length=fusLength,nElements=nElemFus,S=[Sfus],I=[Ifus],aeroSurface=rightFusSurf,rotationParametrization="E321",p0=[fΛ;0;0],pointInertias=[fRightConInertia])

    # Propellers thrust force
    thrustValue = thrustIsTrimVariable ? t -> 0 : thrust

    thrustLeft = create_BC(name="thrustLeft",beam=leftFus,node=1,types=["Ff2A"],values=[t->thrustValue(t)],toBeTrimmed=[thrustIsTrimVariable])

    thrustRight = create_BC(name="thrustRight",beam=rightFus,node=nElemFus+1,types=["Ff2A"],values=[t->thrustValue(t)],toBeTrimmed=[thrustIsTrimVariable])

    # Link propellers' thrust
    thrustsLink = thrustIsTrimVariable ? [create_TrimLoadsLink(masterBC=thrustRight,slaveBCs=[thrustLeft])] : Vector{TrimLoadsLink}()

    # Initial position for left wingtip
    initialPosition = [-kp3[1]; kp3[2]; kp3[3]]

    # Aircraft model (with initial position such that aircraft center is coincident with the origin of frame A)
    BWB = create_Model(name="BWB",beams=[leftWing,leftFus,rightFus,rightWing],BCs=[thrustLeft,thrustRight],initialPosition=initialPosition,v_A=[0;airspeed;0],altitude=altitude,gravityVector=[0;0;g],flapLinks=[elevonLink],trimLoadsLinks=thrustsLink)

    return BWB
end
export create_BWB