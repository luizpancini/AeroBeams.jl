# Gets the normalized nodal positions of the Pazy wing
function Pazy_nodal_positions()
    return [0.0; 0.06956521653730675; 0.13913043671201064; 0.208695655068016; 0.2782608734240213; 0.34782609178002666; 0.41739131195473056; 0.4869565303107358; 0.5565217486667412; 0.626086968841445; 0.6956521871974504; 0.7652174055534556; 0.8347826239094611; 0.9043478440841649; 0.9652173080712125; 1.0]
end


# Gets the sectional stiffness matrices of the Pazy wing
function Pazy_stiffness_matrices(GAy::Number,GAz::Number)

    @assert GAy > 0
    @assert GAz > 0

    # Load elemental sectional stiffness values
    EA = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/EA.txt")))
    GJ = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/GJ.txt")))
    EIy = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/EIy.txt")))
    EIz = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/EIz.txt")))
    c_EA_GJ = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/c_EA_GJ.txt")))
    c_EA_EIy = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/c_EA_EIy.txt")))
    c_EA_EIz = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/c_EA_EIz.txt")))
    c_GJ_EIy = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/c_GJ_EIy.txt")))
    c_GJ_EIz = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/c_GJ_EIz.txt")))
    c_EIy_EIz = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/c_EIy_EIz.txt")))

    # Set matrices
    C = [[    EA[i]  0    0  c_EA_GJ[i]  c_EA_EIy[i]  c_EA_EIz[i];
                  0 GAy   0           0            0            0;
                  0  0  GAz           0            0            0;
         c_EA_GJ[i]  0    0       GJ[i]  c_GJ_EIy[i]  c_GJ_EIz[i];
        c_EA_EIy[i]  0    0 c_GJ_EIy[i]       EIy[i] c_EIy_EIz[i];
        c_EA_EIz[i]  0    0 c_GJ_EIz[i] c_EIy_EIz[i]       EIz[i]] for i in 1:15]

    return C
end


# Gets the sectional inertia matrices of the Pazy wing
function Pazy_inertia_matrices()

    # Length and nodal position
    L = 0.549843728 
    nodalPositions = Pazy_nodal_positions()

    # Load nodal inertia values
    m = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/m.txt")))
    m_x1 = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/m_x1.txt")))
    m_x2 = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/m_x2.txt")))
    m_x3 = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/m_x3.txt")))
    Ixx = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/Ixx.txt")))
    Ixy = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/Ixy.txt")))
    Ixz = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/Ixz.txt")))
    Iyy = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/Iyy.txt")))
    Iyz = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/Iyz.txt")))
    Izz = vec(readdlm(string(pwd(),"/test/referenceData/Pazy/Izz.txt")))

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
    create_Pazy(; kwargs...)

Creates the Pazy wing

Returns the beam and geometrical properties

# Arguments
- `p0::Vector{<:Number}` = initial rotation parameters
- `airfoil::Airfoil=NACA0018` = airfoil section
- `aeroSolver::AeroSolver=Indicial()` = aerodynamic solver
- `gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner")` = indicial gust loads solver
- `derivationMethod::DerivationMethod=AD()` = method for aerodynamic derivatives
- `withTipCorrection::Bool=true` = flag for aerodynamic tip correction 
- `GAy::Number` = shear stiffness in the x2 direction
- `GAz::Number` = shear stiffness in the x3 direction
"""
function create_Pazy(; p0::Vector{<:Number}=zeros(3),airfoil::Airfoil=deepcopy(NACA0018),aeroSolver::AeroSolver=Indicial(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),derivationMethod::DerivationMethod=AD(),withTipCorrection::Bool=true,GAy::Number=1e16,GAz::Number=GAy)

    # Length
    L = 0.549843728

    # Number of elements
    nElem = 15

    # Normalized nodal positions
    nodalPositions = Pazy_nodal_positions()

    # Stiffness matrices
    C = Pazy_stiffness_matrices(GAy,GAz)

    # Inertia matrices
    I = Pazy_inertia_matrices()

    # Chord
    chord = 0.0989

    # Normalized spar position
    normSparPos = 0.44096

    # Aerodynamic surface
    surf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=chord,normSparPos=normSparPos,hasTipCorrection=withTipCorrection)

    # Wing 
    wing = create_Beam(name="pazy",length=L,nElements=nElem,normalizedNodalPositions=nodalPositions,C=C,I=I,rotationParametrization="E321",p0=p0,aeroSurface=surf)

    return wing,L,nElem,chord,normSparPos,airfoil,surf
end
export create_Pazy


"""
    create_PazyFFWT(; kwargs...)

Creates a version of the Pazy wing with flared folding wingtip (FFWT)

# Arguments
- `p0::Vector{<:Number}` = initial rotation parameters
- `airfoil::Airfoil=NACA0018` = airfoil section
- `aeroSolver::AeroSolver=Indicial()` = aerodynamic solver
- `gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner")` = indicial gust loads solver
- `derivationMethod::DerivationMethod=AD()` = method for aerodynamic derivatives
- `withTipCorrection::Bool=true` = flag for aerodynamic tip correction 
- `GAy::Number` = shear stiffness in the x2 direction
- `GAz::Number` = shear stiffness in the x3 direction
- `hingeNode::Int64` = hinge node
- `hingeAngle::Number` = hinge (fold) angle
- `flareAngle::Number` = flare angle
- `kSpring::Number` = stiffness of the hinge
- `g::Number=0` = local acceleration of gravity
- `airspeed::Number` = local airspeed
"""
function create_PazyFFWT(; p0::Vector{<:Number}=zeros(3),airfoil::Airfoil=deepcopy(NACA0018),aeroSolver::AeroSolver=Indicial(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),derivationMethod::DerivationMethod=AD(),withTipCorrection::Bool=false,GAy::Number=1e16,GAz::Number=GAy,hingeNode::Int64=14,hingeAngle::Number=0,flareAngle::Number=10,kSpring::Number=1e6,g::Number=0,airspeed::Number)

    @assert 2 <= hingeNode <= 15 "hingeNode must be between 2 and 15"
    @assert -3π/4 <= hingeAngle <= 3π/4 "set hingeAngle between -3π/4 and 3π/4 "
    @assert 0 <= flareAngle <= 20 "set flareAngle between 0 and 20 (degrees)"
    @assert kSpring >= 0
    @assert g >= 0
    @assert airspeed >= 0

    # Total length 
    L = 0.549843728

    # Number of elements in each part of the wing
    nElem = 15
    nElem1 = hingeNode-1
    nElem2 = nElem-nElem1

    # Normalized nodal positions in each part of the wing
    nodalPositions = Pazy_nodal_positions()
    nodalPositions1 = nodalPositions[1:hingeNode]/nodalPositions[hingeNode]
    nodalPositions2 = (nodalPositions[hingeNode:end].-nodalPositions[hingeNode])./(nodalPositions[end]-nodalPositions[hingeNode])

    # Length in each part of the wing
    L1 = L*nodalPositions[hingeNode]
    L2 = L*(nodalPositions[end]-nodalPositions[hingeNode])

    # Stiffness matrices
    C = Pazy_stiffness_matrices(GAy,GAz)
    C1 = C[1:nElem1]
    C2 = C[nElem1+1:end]

    # Inertia matrices
    I = Pazy_inertia_matrices()
    I1 = I[1:nElem1]
    I2 = I[nElem1+1:end]

    # Chord
    chord = 0.0989

    # Normalized spar position in each part of the wing
    normSparPos = 0.44096
    normSparPos1 = normSparPos
    normSparPos2 = x1 -> normSparPos + x1*tand(flareAngle)/chord
    @assert 0 < normSparPos2(L2) < 1 "flareAngle is too large for the specified hingeNode"

    # Aerodynamic surfaces
    surf1 = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=chord,normSparPos=normSparPos1,hasTipCorrection=withTipCorrection)

    surf2 = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=chord,Λ=-flareAngle*π/180,normSparPos=normSparPos2,hasTipCorrection=withTipCorrection)

    # Beams 
    mainWing = create_Beam(name="mainWing",length=L1,nElements=nElem1,normalizedNodalPositions=nodalPositions1,C=C1,I=I1,rotationParametrization="E321",p0=p0,aeroSurface=surf1,hingedNodes=[nElem1+1],hingedNodesDoF=[[false,true,false]])

    wingTip = create_Beam(name="wingTip",length=L2,nElements=nElem2,normalizedNodalPositions=nodalPositions2,C=C2,I=I2,rotationParametrization="E321",p0=p0+[-flareAngle*π/180;0;0],aeroSurface=surf2)

    # Relative rotation constraint
    rotationConstraint = create_RotationConstraint(masterBeam=mainWing,slaveBeam=wingTip,masterElementLocalID=nElem1,slaveElementLocalID=1,DOF=2,value=4*tan(hingeAngle/4))

    # Spring around hinge
    spring = create_Spring(basis="A",elementsIDs=[nElem1,1],nodesSides=[1,2],kIPBending=kSpring)
    add_spring_to_beams!(beams=[mainWing,wingTip],spring=spring)

    # BCs
    clamp = create_BC(name="clamp",beam=mainWing,node=1,types=["u1A","u2A","u3A","p1A","p2A","p3A"],values=[0,0,0,0,0,0])

    # Wing model
    pazyFFWT = create_Model(name="pazyFFWT",beams=[mainWing,wingTip],BCs=[clamp],gravityVector=[0;0;-g],v_A=[0;airspeed;0],rotationConstraints=[rotationConstraint],units=create_UnitsSystem(frequency="Hz"))

    return pazyFFWT,hingeNode
end
export create_PazyFFWT


"""
    Pazy_tip_loss_factor(αᵣ::Number,U::Number)

Computes the tip loss factor for the Pazy wing's tip correction function

# Arguments
- `αᵣ::Number` = root pitch angle, in degrees
- `U::Number` = airspeed
"""
function Pazy_tip_loss_factor(αᵣ::Number,U::Number)

    # Bound inputs
    αᵣ = min(7,max(αᵣ,0))
    U = min(60,max(U,0))

    # Coefficients as a function of root angle of attack
    αᵣRange = [0; 1; 2; 3; 4; 5; 6; 7]
    τ₀Range = [6.58; 6.29; 6.00; 5.92; 5.91; 6.08; 6.43; 6.88]
    τ₁Range = 1e-2*[0; 3.33; 5.19; 6.01; 6.49; 5.98; 4.43; 2.30]
    τ₂Range = -1e-4*[0; 5.56; 8.59; 10.6; 12.3; 12.8; 11.9; 10.2]
    τ₀ = interpolate(αᵣRange,τ₀Range,αᵣ)
    τ₁ = interpolate(αᵣRange,τ₁Range,αᵣ)
    τ₂ = interpolate(αᵣRange,τ₂Range,αᵣ)
    
    return τ₀ + τ₁*U + τ₂*U^2
end
export Pazy_tip_loss_factor


"""
    create_Helios(; kwargs...)

Creates a model based on the flying-wing aircraft described by Patil and Hodges in: Flight Dynamics of Highly Flexible Flying Wings (2006)

# Keyword arguments
- `altitude::Number` = altitude
- `aeroSolver::AeroSolver` = aerodynamic solver
- `gustLoadsSolver::GustAeroSolver` = indicial gust loads solver
- `derivationMethod::DerivationMethod` = method for aerodynamic derivatives
- `g::Number` = local acceleration of gravity
- `wingAirfoil::Airfoil` = airfoil section of the wing
- `podAirfoil::Airfoil` = airfoil section of the pods
- `beamPods::Bool` = flag to include pods
- `stiffnessFactor::Number` = stiffness factor for the wing structure
- `∞::Number` = value for rigid structural properties
- `nElemStraightSemispan::Int64` = number of elements in the straight section of the semispan
- `nElemDihedralSemispan::Int64` = number of elements in the dihedral section of the semispan
- `nElemPod::Int64` = number of elements in the pods
- `payloadPounds::Number` = payload, in pounds
- `airspeed::Number` = local initial/trim airspeed
- `δIsTrimVariable::Bool` = flag for flap deflection being a trim variable
- `thrustIsTrimVariable::Bool` = flag for motors' thrust being a trim variable
- `δ::Union{Nothing,Number,<:Function}` = flap deflection
- `thrust::Union{Number,<:Function}` = motors' thrust
- `reducedChord::Bool` = flag to employ a reduced (7 ft) chord
- `payloadOnWing::Bool` = flag to set the payload on the wing's reference line
"""
function create_Helios(; altitude::Number=0,aeroSolver::AeroSolver=Indicial(),derivationMethod::DerivationMethod=AD(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),g::Number=-9.80665,wingAirfoil::Airfoil=deepcopy(HeliosWingAirfoil),podAirfoil::Airfoil=HeliosPodAirfoil,beamPods::Bool=false,stiffnessFactor::Number=1.0,∞::Number=1e12,nElemStraightSemispan::Int64=10,nElemDihedralSemispan::Int64=5,nElemPod::Int64=1,payloadPounds::Number=0,airspeed::Number=0,δIsTrimVariable::Bool=false,thrustIsTrimVariable::Bool=false,δ::Union{Nothing,Number,<:Function}=nothing,thrust::Union{Number,<:Function}=0,reducedChord::Bool=false,payloadOnWing::Bool=false)

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
    if δ isa Number
        δconst = deepcopy(δ)
        δ = t -> δconst
    end
    if thrust isa Number
        thrustConst = deepcopy(thrust)
        thrust = t -> thrustConst
    end

    # Wing surface
    wChord = reducedChord ? 7*0.3048 : 8*0.3048
    wNormSparPos = 0.25
    wNormFlapPos = 0.75
    wNormFlapSpan = [0,1]
    wingSurf = δIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=wingAirfoil,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=wNormFlapSpan,δ=δ,updateAirfoilParameters=false) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=wingAirfoil,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=wNormFlapSpan,δIsTrimVariable=δIsTrimVariable,updateAirfoilParameters=false)

    # Wing properties
    wingSection = 40*0.3048
    Γ = 10*π/180
    GJ,EIy,EIz = 1.6530e5,1.0331e6,1.2398e7
    wρA,wρIy,wρIz = 8.929,0.691,3.456
    GJ,EIy,EIz = multiply_inplace!(stiffnessFactor, GJ,EIy,EIz)
    wingC = isotropic_stiffness_matrix(∞=∞,GJ=GJ,EIy=EIy,EIz=EIz)
    wingI = inertia_matrix(ρA=wρA,ρIy=wρIy,ρIz=wρIz)

    # Wing beams
    leftWingDihedral = create_Beam(name="leftWingDihedral",length=wingSection,nElements=nElemDihedralSemispan,C=[wingC],I=[wingI],aeroSurface=deepcopy(wingSurf),rotationParametrization="E321",p0=[0;Γ;0])

    leftWingStraight = create_Beam(name="leftWingStraight",length=2*wingSection,nElements=nElemStraightSemispan,C=[wingC],I=[wingI],aeroSurface=deepcopy(wingSurf))

    rightWingStraight = create_Beam(name="rightWingStraight",length=2*wingSection,nElements=nElemStraightSemispan,C=[wingC],I=[wingI],aeroSurface=deepcopy(wingSurf))

    rightWingDihedral = create_Beam(name="rightWingDihedral",length=wingSection,nElements=nElemDihedralSemispan,C=[wingC],I=[wingI],aeroSurface=deepcopy(wingSurf),rotationParametrization="E321",p0=[0;-Γ;0])

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
    podSurf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,airfoil=podAirfoil,c=wChord,normSparPos=wNormSparPos,updateAirfoilParameters=false)

    # Pod properties
    pρA,pρIy,pρIz = 0,wρIy,wρIz
    podC = isotropic_stiffness_matrix(∞=∞)
    podI = inertia_matrix(ρA=pρA,ρIy=pρIy,ρIz=pρIz)

    # Pod beams
    leftPod = create_Beam(name="leftPod",length=podLength,nElements=nElemPod,C=[podC],I=[podI],aeroSurface=deepcopy(podSurf),rotationParametrization="E321",p0=[0;π/2;0],connectedBeams=[leftWingStraight],connectedNodesThis=[1],connectedNodesOther=[1])

    centerPod = create_Beam(name="centerPod",length=podLength,nElements=nElemPod,C=[podC],I=[podI],aeroSurface=deepcopy(podSurf),rotationParametrization="E321",p0=[0;π/2;0],connectedBeams=[rightWingStraight],connectedNodesThis=[1],connectedNodesOther=[1])

    rightPod = create_Beam(name="rightPod",length=podLength,nElements=nElemPod,C=[podC],I=[podI],aeroSurface=deepcopy(podSurf),rotationParametrization="E321",p0=[0;π/2;0],connectedBeams=[rightWingStraight],connectedNodesThis=[1],connectedNodesOther=[nElemStraightSemispan+1])

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
- `altitude::Number` = altitude
- `aeroSolver::AeroSolver` = aerodynamic solver
- `derivationMethod::DerivationMethod` = method for aerodynamic derivatives
- `flapLoadsSolver::FlapAeroSolver` = aerodynamic solver for flap loads
- `gustLoadsSolver::GustAeroSolver` = indicial gust loads solver
- `stabilizersAero::Bool` = flag for stabilizers with aerodynamic surfaces
- `includeVS::Bool` = flag to include a vertical stabilizer in the model   
- `nElemWing::Int64` = number of elements of the full wing
- `nElemHorzStabilizer::Int64` = number of elements of the horizontal stabilizer
- `nElemTailBoom::Int64` = number of elements of the tail boom
- `nElemVertStabilizer::Int64` = number of elements of the vertical stabilizer
- `∞::Number=1e12` = value of rigid structural properties
- `stiffnessFactor::Number` = stiffness factor for the wing structure
- `k1::Number` = undeformed wing torsional curvature
- `k2::Number` = undeformed wing flapwise bending curvature
- `airspeed::Number` = local initial/trim airspeed
- `δElevIsTrimVariable::Bool` = flag for elevator deflection being a trim variable
- `thrustIsTrimVariable::Bool` = flag for motors' thrust being a trim variable
- `δElev::Union{Nothing,Number,<:Function}` = elevator deflection
- `thrust::Union{Number,<:Function}` = motors' thrust
- `g::Number` = local acceleration of gravity
- `wingCd0::Number` = parisite drag coefficient for the wing
- `wingcnδ::Number` = cn vs δ slope for the wing
- `wingcmδ::Number` = cm vs δ slope for the wing
- `wingcdδ::Number` = cd vs δ slope for the wing
- `stabsCd0::Number` = parisite drag coefficient for the stabilizers
- `stabscnδ::Number` = cn vs δ slope for the stabilizers
- `stabscmδ::Number` = cm vs δ slope for the stabilizers
- `stabscdδ::Number` = cd vs δ slope for the stabilizers
"""
function create_conventional_HALE(; altitude::Number=20e3,aeroSolver::AeroSolver=Indicial(),derivationMethod::DerivationMethod=AD(),flapLoadsSolver::FlapAeroSolver=TableLookup(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),stabilizersAero::Bool=true,includeVS::Bool=true,nElemWing::Int64=20,nElemHorzStabilizer::Int64=10,nElemTailBoom::Int64=10,nElemVertStabilizer::Int64=5,∞::Number=1e12,stiffnessFactor::Number=1,k1::Number=0,k2::Number=0,airspeed::Number=0,δElevIsTrimVariable::Bool=false,thrustIsTrimVariable::Bool=false,δElev::Union{Nothing,Number,<:Function}=nothing,thrust::Union{Number,<:Function}=0,g::Number=-9.80665,wingCd0::Number=0,wingcnδ::Number=2.5,wingcmδ::Number=-0.35,wingcdδ::Number=0.15,stabsCd0::Number=0,stabscnδ::Number=2.5,stabscmδ::Number=-0.35,stabscdδ::Number=0.15)

    # Validate
    @assert iseven(nElemWing)
    @assert iseven(nElemHorzStabilizer)
    @assert ∞ > 1e8
    @assert stiffnessFactor > 0
    @assert wingCd0 >= 0
    @assert stabsCd0 >= 0
    @assert airspeed >= 0
    δElevIsInput = !isnothing(δElev)
    if δElevIsInput
        @assert !δElevIsTrimVariable
    end
    if δElev isa Number
        δElevConst = deepcopy(δElev)
        δElev = t -> δElevConst
    end
    if thrust isa Number
        thrustConst = deepcopy(thrust)
        thrust = t -> thrustConst
    end

    # Wing airfoil
    wAirfoil = deepcopy(flatPlate)

    # Wing surface
    wChord = 1
    wNormSparPos = 0.5
    wNormFlapPos = 0.75
    wNormFlapSpan = [0.75,1]
    wingSurf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=wAirfoil,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=wNormFlapSpan,updateAirfoilParameters=false)
    wingSurf.airfoil.attachedFlowParameters.cd₀ = wingCd0
    wingSurf.airfoil.attachedFlowParameters.cnδ = wingcnδ
    wingSurf.airfoil.attachedFlowParameters.cmδ = wingcmδ
    wingSurf.airfoil.attachedFlowParameters.cdδ = wingcdδ

    # Wing properties
    Lw = 16
    wGJ,wEIy,wEIz = 1e4,2e4,4e6
    wGJ,wEIy,wEIz = multiply_inplace!(stiffnessFactor, wGJ,wEIy,wEIz)
    wρA,wρIs = 0.75,0.1
    wρIy,wρIz = (wEIy/wEIz)*wρIs,(1-wEIy/wEIz)*wρIs
    Cwing = isotropic_stiffness_matrix(∞=∞,GJ=wGJ,EIy=wEIy,EIz=wEIz)
    Iwing = inertia_matrix(ρA=wρA,ρIy=wρIy,ρIz=wρIz,ρIs=wρIs)

    # Initial position for first node of left wing
    ρ = 1/k2
    θ = Lw/ρ
    x = ρ*sin(θ)
    y = ρ*(1-cos(θ))
    initialPosition = k2 == 0 ? [-Lw; 0; 0] : [-x; 0; -y]

    # Initial angle of twist
    r = 1/k1
    ψ = Lw/r

    # Wing beams
    leftWing = create_Beam(name="leftWing",length=Lw,nElements=div(nElemWing,2),C=[Cwing],I=[Iwing],aeroSurface=deepcopy(wingSurf),k=[-k1;k2;0],rotationParametrization="E321",p0=[0;-θ;-ψ])

    rightWing = create_Beam(name="rightWing",length=Lw,nElements=div(nElemWing,2),C=[Cwing],I=[Iwing],aeroSurface=deepcopy(wingSurf),k=[k1;k2;0],connectedBeams=[leftWing],connectedNodesThis=[1],connectedNodesOther=[div(nElemWing,2)+1])

    # Link wing ailerons
    aileronLink = create_FlapLink(masterBeam=rightWing,slaveBeams=[leftWing],δMultipliers=[-1])

    # Tail boom
    Lt = 10
    tρA,tρIy,tρIz = 0.08,wρIy/10,wρIz/10
    tailBoom = create_Beam(name="tailBoom",length=Lt,nElements=nElemTailBoom,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=tρA,ρIy=tρIy,ρIz=tρIz)],rotationParametrization="E321",p0=[-π/2;0;0],connectedBeams=[rightWing],connectedNodesThis=[1],connectedNodesOther=[1])

    # Payload (at wing's spar position)
    payloadMass = 50
    payloadInertia = 200
    payload = PointInertia(elementID=1,η=[-Lw/div(nElemWing,2)/2;0;0],mass=payloadMass,Iyy=payloadInertia,Izz=payloadInertia,Ixx=payloadInertia)
    add_point_inertias_to_beam!(rightWing,inertias=[payload])

    # Stabilizers airfoil
    sAirfoil = deepcopy(flatPlate)

    # Horizontal stabilizer surface
    hChord = 0.5
    hNormSparPos = 0.5
    hNormFlapPos = 0.75
    hNormFlapSpan = [0,1]
    hsSurf = δElevIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=sAirfoil,c=hChord,normSparPos=hNormSparPos,normFlapPos=hNormFlapPos,normFlapSpan=hNormFlapSpan,δ=δElev,updateAirfoilParameters=false) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=sAirfoil,c=hChord,normSparPos=hNormSparPos,normFlapPos=hNormFlapPos,normFlapSpan=hNormFlapSpan,δIsTrimVariable=δElevIsTrimVariable,updateAirfoilParameters=false)
    hsSurf.airfoil.attachedFlowParameters.cd₀ = stabsCd0
    hsSurf.airfoil.attachedFlowParameters.cnδ = stabscnδ
    hsSurf.airfoil.attachedFlowParameters.cmδ = stabscmδ
    hsSurf.airfoil.attachedFlowParameters.cdδ = stabscdδ

    # Horizontal stabilizer beam
    Lh = 5
    hρA,hρIy,hρIz = 0.08,wρIy/10,wρIz/10
    horzStabilizer = create_Beam(name="horzStabilizer",length=Lh,initialPosition=[-Lh/2;0;0],nElements=nElemHorzStabilizer,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=hρA,ρIy=hρIy,ρIz=hρIz)],connectedBeams=[tailBoom],connectedNodesThis=[div(nElemHorzStabilizer,2)+1],connectedNodesOther=[nElemTailBoom+1])
    if stabilizersAero
        horzStabilizer.aeroSurface = hsSurf
        update_beam!(horzStabilizer)
    end

    # Vertical stabilizer surface
    vChord = 0.5
    vNormSparPos = 0.5
    vNormFlapPos = 0.75
    vNormFlapSpan = [0,1]
    vsSurf = create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=flapLoadsSolver,airfoil=sAirfoil,c=vChord,normSparPos=vNormSparPos,normFlapPos=vNormFlapPos,normFlapSpan=vNormFlapSpan,updateAirfoilParameters=false)
    vsSurf.airfoil.attachedFlowParameters.cd₀ = stabsCd0
    vsSurf.airfoil.attachedFlowParameters.cnδ = stabscnδ
    vsSurf.airfoil.attachedFlowParameters.cmδ = stabscmδ
    vsSurf.airfoil.attachedFlowParameters.cdδ = stabscdδ

    # Vertical stabilizer beam
    Lv = 2.5
    vρA,vρIy,vρIz = 0.08,hρIy,hρIz
    vertStabilizer = create_Beam(name="vertStabilizer",length=Lv,nElements=nElemVertStabilizer,C=[isotropic_stiffness_matrix(∞=∞)],I=[inertia_matrix(ρA=vρA,ρIy=vρIy,ρIz=vρIz)],rotationParametrization="E321",p0=[0;-π/2;0],connectedBeams=[tailBoom],connectedNodesThis=[1],connectedNodesOther=[nElemTailBoom+1])
    if stabilizersAero
        vertStabilizer.aeroSurface = vsSurf
        update_beam!(vertStabilizer)
    end

    # Propeller thrust force 
    thrustValue = thrustIsTrimVariable ? t -> 0 : thrust

    propThrust = create_BC(name="propThrust",beam=rightWing,node=1,types=["Ff2b"],values=[t->thrustValue(t)],toBeTrimmed=[thrustIsTrimVariable])

    # Beams of the model
    beams = includeVS ? [leftWing,rightWing,tailBoom,horzStabilizer,vertStabilizer] : [leftWing,rightWing,tailBoom,horzStabilizer]

    # Aircraft model (with initial position such that aircraft center is coincident with the origin of frame A)
    conventional_HALE = create_Model(name="conventionalHALE",beams=beams,BCs=[propThrust],initialPosition=initialPosition,v_A=[0;airspeed;0],altitude=altitude,gravityVector=[0;0;g],flapLinks=[aileronLink])

    return conventional_HALE,leftWing,rightWing,tailBoom,horzStabilizer,vertStabilizer
end
export create_conventional_HALE


"""
    create_BWB(; kwargs...)

Creates a model based on the blended-wing-body described by Weihua Su's PhD thesis

# Keyword arguments
- `altitude::Number` = altitude
- `aeroSolver::AeroSolver` = aerodynamic solver
- `gustLoadsSolver::GustAeroSolver` = indicial gust loads solver
- `derivationMethod::DerivationMethod` = method for aerodynamic derivatives 
- `∞::Number=1e12` = value of rigid structural properties
- `stiffnessFactor::Number` = stiffness factor for the wing structure
- airspeed::Number = local initial/trim airspeed
- `δElevIsTrimVariable::Bool` = flag for elevator deflection being a trim variable
- `thrustIsTrimVariable::Bool` = flag for motors' thrust being a trim variable
- `δElev::Union{Nothing,Number,<:Function}` = elevator deflection
- `thrust::Union{Number,<:Function}` = motors' thrust
- `g::Number` = local acceleration of gravity
- `updateAirfoilParameters::Bool` = flag to update airfoil parameters with airspeed
- `hasTipCorrection::Bool` = flag to employ aerodynamic tip correction
- `tipLossDecayFactor::Number` = tip loss decay factor
"""
function create_BWB(; altitude::Number=0,aeroSolver::AeroSolver=Indicial(),gustLoadsSolver::GustAeroSolver=IndicialGust("Kussner"),derivationMethod::DerivationMethod=AD(),∞::Number=1e12,stiffnessFactor::Number=1,airspeed::Number=0,δElevIsTrimVariable::Bool=false,thrustIsTrimVariable::Bool=false,δElev::Union{Nothing,Number,<:Function}=nothing,thrust::Union{Number,<:Function}=0,g::Number=-9.80665,updateAirfoilParameters::Bool=false,hasTipCorrection::Bool=false,tipLossDecayFactor::Number=40)

    # Validate
    @assert ∞ > 1e8
    @assert stiffnessFactor > 0
    @assert airspeed >= 0
    δElevIsInput = !isnothing(δElev)
    if δElevIsInput
        @assert !δElevIsTrimVariable
    end
    if δElev isa Number
        δElevConst = deepcopy(δElev)
        δElev = t -> δElevConst
    end
    if thrust isa Number
        thrustConst = deepcopy(thrust)
        thrust = t -> thrustConst
    end

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

    # Wing surfaces
    wChord = tipChord
    wNormSparPos = tipSparPos
    wNormFlapPos = 0.75

    leftWingSurf = δElevIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=deepcopy(BWBAirfoil),Λ=wΛ,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=[1/4; 1],δ=δElev,updateAirfoilParameters=updateAirfoilParameters,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=-tipLossDecayFactor) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=deepcopy(BWBAirfoil),Λ=wΛ,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=[1/4; 1],δIsTrimVariable=δElevIsTrimVariable,updateAirfoilParameters=updateAirfoilParameters,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=-tipLossDecayFactor)

    rightWingSurf = δElevIsInput ? create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=deepcopy(BWBAirfoil),Λ=-wΛ,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=[0; 3/4],δ=δElev,updateAirfoilParameters=updateAirfoilParameters,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor) : create_AeroSurface(solver=aeroSolver,gustLoadsSolver=gustLoadsSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=deepcopy(BWBAirfoil),Λ=-wΛ,c=wChord,normSparPos=wNormSparPos,normFlapPos=wNormFlapPos,normFlapSpan=[0; 3/4],δIsTrimVariable=δElevIsTrimVariable,updateAirfoilParameters=updateAirfoilParameters,hasTipCorrection=hasTipCorrection,tipLossDecayFactor=tipLossDecayFactor)

    # Wing properties
    nElemWing = 8
    wEA,wGJ,wEIy,wEIz = 155_000_000,11_000,11_700,130_000
    wGJ,wEIy,wEIz = multiply_inplace!(stiffnessFactor, wGJ,wEIy,wEIz)
    wρA,wρIy,wρIz = 6.2,0.0005,0.00462
    Cwing = isotropic_stiffness_matrix(∞=∞,EA=wEA,GJ=wGJ,EIy=wEIy,EIz=wEIz)
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
    leftWing = create_Beam(name="leftWing",length=wingSemispan,nElements=nElemWing,C=[Cwing],I=[Iwing],aeroSurface=leftWingSurf,rotationParametrization="E321",p0=[wΛ;0;0],pointInertias=wLConInertias)

    rightWing = create_Beam(name="rightWing",length=wingSemispan,nElements=nElemWing,C=[Cwing],I=[Iwing],aeroSurface=rightWingSurf,rotationParametrization="E321",p0=[-wΛ;0;0],pointInertias=wRConInertias)

    # Link wing elevons
    elevonLink = create_FlapLink(masterBeam=rightWing,slaveBeams=[leftWing],δMultipliers=[1])

    # Fuselage surfaces
    leftFusChord = x1 -> tipChord + (rootChord-tipChord)*x1/fusLength
    rightFusChord = x1 -> rootChord + (tipChord-rootChord)*x1/fusLength
    leftFusSparPos = x1 -> tipSparPos + (rootSparPos-tipSparPos)*x1/fusLength
    rightFusSparPos = x1 -> rootSparPos + (tipSparPos-rootSparPos)*x1/fusLength

    leftFusSurf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=deepcopy(BWBAirfoil),Λ=-fΛ,c=leftFusChord,normSparPos=leftFusSparPos,updateAirfoilParameters=updateAirfoilParameters)

    rightFusSurf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,flapLoadsSolver=TableLookup(),airfoil=deepcopy(BWBAirfoil),Λ=fΛ,c=rightFusChord,normSparPos=rightFusSparPos,updateAirfoilParameters=updateAirfoilParameters)

    # Fuselage properties
    nElemFus = 3
    fEA,fGJ,fEIy,fEIz = 169_000_000,2_250_000,750_000,35_000_000
    fρA,fρIy,fρIz = 50,0.7,22.0
    Cfus = isotropic_stiffness_matrix(∞=∞,EA=fEA,GJ=fGJ,EIy=fEIy,EIz=fEIz)
    Ifus = inertia_matrix(ρA=fρA,ρIy=fρIy,ρIz=fρIz)

    # Concentrated fuselage inertias
    fConMass = 40
    fConMassOffset = 0.891968
    fLeftConInertia = PointInertia(elementID=3,η=[fusLength/nElemFus/2;fConMassOffset;0],mass=fConMass)
    fRightConInertia = PointInertia(elementID=1,η=[-fusLength/nElemFus/2;fConMassOffset;0],mass=fConMass)

    # Fuselage beams
    leftFus = create_Beam(name="leftFus",length=fusLength,nElements=nElemFus,C=[Cfus],I=[Ifus],aeroSurface=leftFusSurf,rotationParametrization="E321",p0=[-fΛ;0;0],pointInertias=[fLeftConInertia])

    rightFus = create_Beam(name="rightFus",length=fusLength,nElements=nElemFus,C=[Cfus],I=[Ifus],aeroSurface=rightFusSurf,rotationParametrization="E321",p0=[fΛ;0;0],pointInertias=[fRightConInertia])

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