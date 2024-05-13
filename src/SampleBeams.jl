"""
Pazy_nodal_positions()

Gets the normalized nodal positions of the Pazy wing

# Arguments
"""
function Pazy_nodal_positions()
    return [0.0; 0.06956521653730675; 0.13913043671201064; 0.208695655068016; 0.2782608734240213; 0.34782609178002666; 0.41739131195473056; 0.4869565303107358; 0.5565217486667412; 0.626086968841445; 0.6956521871974504; 0.7652174055534556; 0.8347826239094611; 0.9043478440841649; 0.9652173080712125; 1.0]
end


"""
Pazy_stiffness_matrices(GAy::Number,GAz::Number)

Gets the sectional stiffness matrices of the Pazy wing

# Arguments
- GAy::Number
- GAz::Number
"""
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

"""
Pazy_inertia_matrices()

Gets the sectional inertia matrices of the Pazy wing

"""
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
create_Pazy(; p0::Vector{<:Number}=zeros(3),airfoil::Airfoil=NACA0018,aeroSolver::AeroSolver=Indicial(),derivationMethod::DerivationMethod=AD(),withTipCorrection::Bool=true,GAy::Number=1e16,GAz::Number=GAy)

Creates the Pazy wing

# Arguments
- p0::Vector{<:Number}
- airfoil::Airfoil=NACA0018
- aeroSolver::AeroSolver=Indicial()
- derivationMethod::DerivationMethod=AD()
- withTipCorrection::Bool=true
- GAy::Number=1e16
- GAz::Number=GAy
"""
function create_Pazy(; p0::Vector{<:Number}=zeros(3),airfoil::Airfoil=NACA0018,aeroSolver::AeroSolver=Indicial(),derivationMethod::DerivationMethod=AD(),withTipCorrection::Bool=true,GAy::Number=1e16,GAz::Number=GAy)

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
    surf = create_AeroSurface(solver=aeroSolver,derivationMethod=derivationMethod,airfoil=airfoil,c=chord,normSparPos=normSparPos,hasTipCorrection=withTipCorrection)

    # Wing 
    wing = create_Beam(name="pazy",length=L,nElements=nElem,normalizedNodalPositions=nodalPositions,C=C,I=I,rotationParametrization="E321",p0=p0,aeroSurface=surf)

    return wing,L,nElem,chord,normSparPos,airfoil,surf
end
export create_Pazy


"""
Pazy_tip_loss_factor(αᵣ::Number,U::Number)

Computes the tip loss factor for the Pazy wing's tip correction function

# Arguments
- αᵣ::Number
- U::Number
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