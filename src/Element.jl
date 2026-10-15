## ELEMENT COMPOSITE TYPES
# ------------------------------------------------------------------------------
# Element geometry data composite type
# Fields
# - `Δℓ`: element (arc)length
# - `x1`: position of element's midpoint along the beam arclength
# - `x1_n1`: position of node 1 along the beam arclength
# - `x1_n2`: position of node 2 along the beam arclength
# - `x1_norm`: normalized position of element's midpoint along the beam arclength
# - `x1_n1_norm`: normalized position of node 1 along the beam arclength
# - `x1_n2_norm`: normalized position of node 2 along the beam arclength
# - `k`: element's curvature vector per unit length
# - `R0`: element's midpoint rotation matrix from basis A to basis b
# - `R0T`: element's midpoint rotation matrix from basis b to basis A
# - `R0_n1`: node 1 rotation matrix from basis A to basis b
# - `R0_n2`: node 2 rotation matrix from basis A to basis b
# - `R0T_n1`: node 1 rotation matrix from basis b to basis A
# - `R0T_n2`: node 2 rotation matrix from basis b to basis A
# - `hingedNode1Mat`: boolean matrix indicating which DOFs of node 1 are hinged
# - `notHingedNode1Mat`: boolean matrix indicating which DOFs of node 1 are not hinged
# - `notHingedNode2Mat`: boolean matrix indicating which DOFs of node 2 are not hinged
#
struct ElementGeometry

    Δℓ::Float64
    x1::Float64
    x1_n1::Float64
    x1_n2::Float64
    x1_norm::Float64
    x1_n1_norm::Float64
    x1_n2_norm::Float64

    k::SVector{3,Float64}
    R0::SMatrix{3,3,Float64,9}
    R0T::SMatrix{3,3,Float64,9}
    R0_n1::SMatrix{3,3,Float64,9}
    R0_n2::SMatrix{3,3,Float64,9}
    R0T_n1::SMatrix{3,3,Float64,9}
    R0T_n2::SMatrix{3,3,Float64,9}

    hingedNode1Mat::SMatrix{3,3,Bool,9}
    notHingedNode1Mat::SMatrix{3,3,Bool,9}
    notHingedNode2Mat::SMatrix{3,3,Bool,9}

end


# Element sectional data composite type
# Fields
# - `C`: sectional compliance matrix
# - `I`: sectional inertia matrix
# - `C11`: sectional compliance matrix block 11
# - `C12`: sectional compliance matrix block 12
# - `C21`: sectional compliance matrix block 21
# - `C22`: sectional compliance matrix block 22
# - `I11`: sectional inertia matrix block 11
# - `I12`: sectional inertia matrix block 12
# - `I21`: sectional inertia matrix block 21
# - `I22`: sectional inertia matrix block 22
# - `ηtilde`: skew-symmetric matrix of the sectional mass center position vector
# - `μ`: sectional mass per unit length
#
struct ElementSectional

    C::SMatrix{6,6,Float64,36}
    I::SMatrix{6,6,Float64,36}

    C11::SMatrix{3,3,Float64,9}
    C12::SMatrix{3,3,Float64,9}
    C21::SMatrix{3,3,Float64,9}
    C22::SMatrix{3,3,Float64,9}

    I11::SMatrix{3,3,Float64,9}
    I12::SMatrix{3,3,Float64,9}
    I21::SMatrix{3,3,Float64,9}
    I22::SMatrix{3,3,Float64,9}

    ηtilde::SMatrix{3,3,Float64,9}
    μ::Float64

end


# Element distributed loads data composite type
# Fields
# - `f_A_of_ζt`: distributed dead forces defined in basis A as a function of ζ and time
# - `m_A_of_ζt`: distributed dead moments defined in basis A as a function of ζ and time
# - `f_b_of_ζt`: distributed dead forces defined in basis b as a function of ζ and time
# - `m_b_of_ζt`: distributed dead moments defined in basis b as a function of ζ and time
# - `ff_A_of_ζt`: distributed follower forces defined in basis A as a function of ζ and time
# - `mf_A_of_ζt`: distributed follower moments defined in basis A as a function of ζ and time
# - `ff_b_of_ζt`: distributed follower forces defined in basis b as a function of ζ and time
# - `mf_b_of_ζt`: distributed follower moments defined in basis b as a function of ζ and time
# - `has_f_A`: flag indicating if there are non-zero distributed dead forces defined in basis A
# - `has_m_A`: flag indicating if there are non-zero distributed dead moments defined in basis A
# - `has_f_b`: flag indicating if there are non-zero distributed dead forces defined in basis b
# - `has_m_b`: flag indicating if there are non-zero distributed dead moments defined in basis b
# - `has_ff_A`: flag indicating if there are non-zero distributed follower forces defined in basis A
# - `has_mf_A`: flag indicating if there are non-zero distributed follower moments defined in basis A
# - `has_ff_b`: flag indicating if there are non-zero distributed follower forces defined in basis b
# - `has_mf_b`: flag indicating if there are non-zero distributed follower moments defined in basis b
# - `f_A`: distributed dead forces' nodal values defined in basis A over time
# - `m_A`: distributed dead moments' nodal values defined in basis A over time
# - `f_b`: distributed dead forces' nodal values defined in basis b over time
# - `m_b`: distributed dead moments' nodal values defined in basis b over time
# - `ff_A`: distributed follower forces' nodal values defined in basis A over time
# - `mf_A`: distributed follower moments' nodal values defined in basis A over time
# - `ff_b`: distributed follower forces' nodal values defined in basis b over time
# - `mf_b`: distributed follower moments' nodal values defined in basis b over time
#
mutable struct ElementDistributedLoads

    f_A_of_ζt::Function
    m_A_of_ζt::Function
    f_b_of_ζt::Function
    m_b_of_ζt::Function
    ff_A_of_ζt::Function
    mf_A_of_ζt::Function
    ff_b_of_ζt::Function
    mf_b_of_ζt::Function

    has_f_A::Bool
    has_m_A::Bool
    has_f_b::Bool
    has_m_b::Bool
    has_ff_A::Bool
    has_mf_A::Bool
    has_ff_b::Bool
    has_mf_b::Bool

    f_A::Array{Float64,3}
    m_A::Array{Float64,3}
    f_b::Array{Float64,3}
    m_b::Array{Float64,3}
    ff_A::Array{Float64,3}
    mf_A::Array{Float64,3}
    ff_b::Array{Float64,3}
    mf_b::Array{Float64,3}

end
ElementDistributedLoads() = ElementDistributedLoads(
    ntuple(_ -> (ζ,t) -> z3, 8)...,
    ntuple(_ -> false, 8)...,
    ntuple(_ -> Array{Float64,3}(undef,2,3,0), 8)...,
)


# Element assembly data composite type
# Fields
# - `parentBeamID`: parent beam's global ID on assembly
# - `localID`: element's local ID (on parent beam)
# - `globalID`: element's global ID (on beam assembly)
# - `nodesLocalIDs`: local IDs (on parent beam) of the element's nodes
# - `nodesGlobalIDs`: global IDs (on beam assembly) of the element's nodes
# - `r`: element's midpoint position vector w.r.t. origin of assembly (basis A)
# - `r_n1`: node 1 position vector w.r.t. origin of assembly (basis A)
# - `r_n2`: node 2 position vector w.r.t. origin of assembly (basis A)
#
struct ElementAssemblyData

    parentBeamID::Int
    localID::Int
    globalID::Int
    nodesLocalIDs::SVector{2,Int}
    nodesGlobalIDs::SVector{2,Int}

    r::SVector{3,Float64}
    r_n1::SVector{3,Float64}
    r_n2::SVector{3,Float64}

end
ElementAssemblyData(localID::Int,nodesLocalIDs::AbstractVector{Int}) =
    ElementAssemblyData(
        0, localID, 0, SVector{2,Int}(nodesLocalIDs), zeros(SVector{2,Int}),
        ntuple(_ -> zeros(SVector{3,Float64}), 3)...,
    )   

ElementAssemblyData(element; parentBeamID::Int,globalID::Int,nodesGlobalIDs::AbstractVector{Int},r::AbstractVector{Float64},r_n1::AbstractVector{Float64},r_n2::AbstractVector{Float64}) =
    ElementAssemblyData(
        parentBeamID, element.assemblyData.localID, globalID, element.assemblyData.nodesLocalIDs, SVector{2,Int}(nodesGlobalIDs),
        SVector{3,Float64}(r), SVector{3,Float64}(r_n1), SVector{3,Float64}(r_n2),
    )


# Element assembly indices composite type
# Fields
# - `eqs_Fu1`: global indices of the first node translational equilibrium equations
# - `eqs_Fu2`: global indices of the second node translational equilibrium equations
# - `eqs_Fp1`: global indices of the first node rotational equilibrium equations
# - `eqs_Fp2`: global indices of the second node rotational equilibrium equations
# - `eqs_FF1`: global indices of the first node translational compatibility equations
# - `eqs_FF2`: global indices of the second node translational compatibility equations
# - `eqs_FM1`: global indices of the first node rotational compatibility equations
# - `eqs_FM2`: global indices of the second node rotational compatibility equations
# - `eqs_FV`: global indices of the midpoint linear velocity-displacement compatibility equations
# - `eqs_FΩ`: global indices of the midpoint angular velocity-displacement compatibility equations
# - `eqs_Fχ`: global indices of the aerodynamic equations
# - `DOF_u`: global DOFs of the translational displacements
# - `DOF_p`: global DOFs of the rotational displacements
# - `DOF_F`: global DOFs of the sectional forces
# - `DOF_M`: global DOFs of the sectional moments
# - `DOF_V`: global DOFs of the sectional velocities
# - `DOF_Ω`: global DOFs of the sectional angular velocities
# - `DOF_χ`: global DOFs of the aerodynamic states
# - `isSpecialNode1`: flag indicating if node 1 is a special node
# - `isSpecialNode2`: flag indicating if node 2 is a special node
# - `eqsNode1Set`: flag indicating if the equations indices of node 1 have been set during the assembly process
# - `eqsNode2Set`: flag indicating if the equations indices of node 2 have been set during the assembly process
#
struct ElementAssemblyIndices{N}

    eqs_Fu1::SVector{3,Int} 
    eqs_Fu2::SVector{3,Int}
    eqs_Fp1::SVector{3,Int}
    eqs_Fp2::SVector{3,Int}
    eqs_FF1::SVector{3,Int}
    eqs_FF2::SVector{3,Int}
    eqs_FM1::SVector{3,Int}
    eqs_FM2::SVector{3,Int}
    eqs_FV::SVector{3,Int}
    eqs_FΩ::SVector{3,Int}
    eqs_Fχ::SVector{N,Int}

    DOF_u::SVector{3,Int}
    DOF_p::SVector{3,Int}
    DOF_F::SVector{3,Int}
    DOF_M::SVector{3,Int}
    DOF_V::SVector{3,Int}
    DOF_Ω::SVector{3,Int}
    DOF_χ::SVector{N,Int}

    isSpecialNode1::Bool
    isSpecialNode2::Bool
    eqsNode1Set::Bool
    eqsNode2Set::Bool

end
ElementAssemblyIndices{N}() where N = ElementAssemblyIndices{N}(
    ntuple(_ -> zeros(SVector{3,Int}), 10)..., zeros(SVector{N,Int}),
    ntuple(_ -> zeros(SVector{3,Int}), 6)..., zeros(SVector{N,Int}),
    ntuple(_ -> false, 4)...,
)


# Jacobian indices maps composite type
# Fields (global indices of the Jacobian block mappings)
# - `Fu1_p`: first node translational equilibrium equations and rotational displacements DOF
# - `Fu1_F`: first node translational equilibrium equations and sectional forces DOF
# - `Fu1_V`: first node translational equilibrium equations and sectional velocities DOF
# - `Fu1_Ω`: first node translational equilibrium equations and sectional angular velocities DOF
# - `Fu1_χ`: first node translational equilibrium equations and aerodynamic DOF
# - `Fu1_δ`: first node translational equilibrium equations and trim flap deflections DOF
# - `Fp1_p`: first node rotational equilibrium equations and rotational displacements DOF
# - `Fp1_F`: first node rotational equilibrium equations and sectional forces DOF
# - `Fp1_M`: first node rotational equilibrium equations and sectional moments DOF
# - `Fp1_V`: first node rotational equilibrium equations and sectional velocities DOF
# - `Fp1_Ω`: first node rotational equilibrium equations and sectional angular velocities DOF
# - `Fp1_χ`: first node rotational equilibrium equations and aerodynamic DOF
# - `Fp1_δ`: first node rotational equilibrium equations and trim flap deflections DOF
# - `FF1_u`: first node translational compatibility equations and translational displacements DOF
# - `FF1_p`: first node translational compatibility equations and rotational displacements DOF
# - `FF1_F`: first node translational compatibility equations and sectional forces DOF
# - `FF1_M`: first node translational compatibility equations and sectional moments DOF
# - `FM1_p`: first node rotational compatibility equations and rotational displacements DOF
# - `FM1_F`: first node rotational compatibility equations and sectional forces DOF
# - `FM1_M`: first node rotational compatibility equations and sectional moments DOF
# - `Fu2_p`: second node translational equilibrium equations and rotational displacements DOF
# - `Fu2_F`: second node translational equilibrium equations and sectional forces DOF
# - `Fu2_V`: second node translational equilibrium equations and sectional velocities DOF
# - `Fu2_Ω`: second node translational equilibrium equations and sectional angular velocities DOF
# - `Fu2_χ`: second node translational equilibrium equations and aerodynamic DOF
# - `Fu2_δ`: second node translational equilibrium equations and trim flap deflections DOF
# - `Fp2_p`: second node rotational equilibrium equations and rotational displacements DOF
# - `Fp2_F`: second node rotational equilibrium equations and sectional forces DOF
# - `Fp2_M`: second node rotational equilibrium equations and sectional moments DOF
# - `Fp2_V`: second node rotational equilibrium equations and sectional velocities DOF
# - `Fp2_Ω`: second node rotational equilibrium equations and sectional angular velocities DOF
# - `Fp2_χ`: second node rotational equilibrium equations and aerodynamic DOF
# - `Fp2_δ`: second node rotational equilibrium equations and trim flap deflections DOF
# - `FF2_u`: second node translational compatibility equations and translational displacements DOF
# - `FF2_p`: second node translational compatibility equations and rotational displacements DOF
# - `FF2_F`: second node translational compatibility equations and sectional forces DOF
# - `FF2_M`: second node translational compatibility equations and sectional moments DOF
# - `FM2_p`: second node rotational compatibility equations and rotational displacements DOF
# - `FM2_F`: second node rotational compatibility equations and sectional forces DOF
# - `FM2_M`: second node rotational compatibility equations and sectional moments DOF
# - `FV_u`: midpoint linear velocity-displacement compatibility equations and translational displacements DOF
# - `FV_p`: midpoint linear velocity-displacement compatibility equations and rotational displacements DOF
# - `FV_V`: midpoint linear velocity-displacement compatibility equations and sectional velocities DOF
# - `FΩ_p`: midpoint angular velocity-rotation compatibility equations and rotational displacements DOF
# - `FΩ_Ω`: midpoint angular velocity-rotation compatibility equations and sectional angular velocities DOF
# - `Fχ_p`: aerodynamic equations and rotational DOF
# - `Fχ_V`: aerodynamic equations and sectional velocities DOF
# - `Fχ_Ω`: aerodynamic equations and sectional angular velocities DOF
# - `Fχ_χ`: aerodynamic equations and aerodynamic forces DOF
#
struct JacobianIndicesMaps{N3,Nδ3,NN}

    Fu1_p::SVector{9,Int}
    Fu1_F::SVector{9,Int}
    Fu1_V::SVector{9,Int}
    Fu1_Ω::SVector{9,Int}
    Fu1_χ::SVector{N3,Int}
    Fu1_δ::SVector{Nδ3,Int}

    Fp1_p::SVector{9,Int}
    Fp1_F::SVector{9,Int}
    Fp1_M::SVector{9,Int}
    Fp1_V::SVector{9,Int}
    Fp1_Ω::SVector{9,Int}
    Fp1_χ::SVector{N3,Int}
    Fp1_δ::SVector{Nδ3,Int}

    FF1_u::SVector{9,Int}
    FF1_p::SVector{9,Int}
    FF1_F::SVector{9,Int}
    FF1_M::SVector{9,Int}

    FM1_p::SVector{9,Int}
    FM1_F::SVector{9,Int}
    FM1_M::SVector{9,Int}

    Fu2_p::SVector{9,Int}
    Fu2_F::SVector{9,Int}
    Fu2_V::SVector{9,Int}
    Fu2_Ω::SVector{9,Int}
    Fu2_χ::SVector{N3,Int}
    Fu2_δ::SVector{Nδ3,Int}

    Fp2_p::SVector{9,Int}
    Fp2_F::SVector{9,Int}
    Fp2_M::SVector{9,Int}
    Fp2_V::SVector{9,Int}
    Fp2_Ω::SVector{9,Int}
    Fp2_χ::SVector{N3,Int}
    Fp2_δ::SVector{Nδ3,Int}

    FF2_u::SVector{9,Int}
    FF2_p::SVector{9,Int}
    FF2_F::SVector{9,Int}
    FF2_M::SVector{9,Int}

    FM2_p::SVector{9,Int}
    FM2_F::SVector{9,Int}
    FM2_M::SVector{9,Int}

    FV_u::SVector{9,Int}
    FV_p::SVector{9,Int}
    FV_V::SVector{9,Int}

    FΩ_p::SVector{9,Int}
    FΩ_Ω::SVector{9,Int}

    Fχ_p::SVector{N3,Int}
    Fχ_V::SVector{N3,Int}
    Fχ_Ω::SVector{N3,Int}
    Fχ_χ::SVector{NN,Int}

    function JacobianIndicesMaps{N3,Nδ3,NN}(Fu1_p,Fu1_F,Fu1_V,Fu1_Ω,Fu1_χ,Fu1_δ,Fp1_p,Fp1_F,Fp1_M,Fp1_V,Fp1_Ω,Fp1_χ,Fp1_δ,FF1_u,FF1_p,FF1_F,FF1_M,FM1_p,FM1_F,FM1_M,Fu2_p,Fu2_F,Fu2_V,Fu2_Ω,Fu2_χ,Fu2_δ,Fp2_p,Fp2_F,Fp2_M,Fp2_V,Fp2_Ω,Fp2_χ,Fp2_δ,FF2_u,FF2_p,FF2_F,FF2_M,FM2_p,FM2_F,FM2_M,FV_u,FV_p,FV_V,FΩ_p,FΩ_Ω,Fχ_p,Fχ_V,Fχ_Ω,Fχ_χ) where {N3,Nδ3,NN}
        new{N3,Nδ3,NN}(
            Fu1_p,Fu1_F,Fu1_V,Fu1_Ω,Fu1_χ,Fu1_δ,
            Fp1_p,Fp1_F,Fp1_M,Fp1_V,Fp1_Ω,Fp1_χ,Fp1_δ,
            FF1_u,FF1_p,FF1_F,FF1_M,
            FM1_p,FM1_F,FM1_M,
            Fu2_p,Fu2_F,Fu2_V,Fu2_Ω,Fu2_χ,Fu2_δ,
            Fp2_p,Fp2_F,Fp2_M,Fp2_V,Fp2_Ω,Fp2_χ,Fp2_δ,
            FF2_u,FF2_p,FF2_F,FF2_M,
            FM2_p,FM2_F,FM2_M,
            FV_u,FV_p,FV_V,
            FΩ_p,FΩ_Ω,
            Fχ_p,Fχ_V,Fχ_Ω,Fχ_χ
        )
    end

end
JacobianIndicesMaps{N,Nδ}() where {N,Nδ} = 
    JacobianIndicesMaps{3*N,3*Nδ,N*N}(
        ntuple(_ -> zeros(SVector{9, Int}), 4)..., zeros(SVector{3*N, Int}), zeros(SVector{3*Nδ, Int}),
        ntuple(_ -> zeros(SVector{9, Int}), 5)..., zeros(SVector{3*N, Int}), zeros(SVector{3*Nδ, Int}),
        ntuple(_ -> zeros(SVector{9, Int}), 4)...,
        ntuple(_ -> zeros(SVector{9, Int}), 3)...,
        ntuple(_ -> zeros(SVector{9, Int}), 4)..., zeros(SVector{3*N, Int}), zeros(SVector{3*Nδ, Int}),
        ntuple(_ -> zeros(SVector{9, Int}), 5)..., zeros(SVector{3*N, Int}), zeros(SVector{3*Nδ, Int}),
        ntuple(_ -> zeros(SVector{9, Int}), 4)...,
        ntuple(_ -> zeros(SVector{9, Int}), 3)...,
        ntuple(_ -> zeros(SVector{9, Int}), 3)...,
        ntuple(_ -> zeros(SVector{9, Int}), 2)...,
        zeros(SVector{3*N, Int}), zeros(SVector{3*N, Int}), zeros(SVector{3*N, Int}), zeros(SVector{N*N, Int})
    )


# Inertia indices maps composite type
# Fields (global indices of the inertia matrix mappings)    
# - `Fu1_pdot`: first node translational equilibrium equations and rotational displacements DOF
# - `Fu1_Vdot`: first node translational equilibrium equations and sectional velocities DOF
# - `Fu1_Ωdot`: first node translational equilibrium equations and sectional angular velocities DOF
# - `Fp1_pdot`: first node rotational equilibrium equations and rotational displacements DOF
# - `Fp1_Vdot`: first node rotational equilibrium equations and sectional velocities DOF
# - `Fp1_Ωdot`: first node rotational equilibrium equations and sectional angular velocities DOF
# - `Fu2_pdot`: second node translational equilibrium equations and rotational displacements DOF
# - `Fu2_Vdot`: second node translational equilibrium equations and sectional velocities DOF
# - `Fu2_Ωdot`: second node translational equilibrium equations and sectional angular velocities DOF
# - `Fp2_pdot`: second node rotational equilibrium equations and rotational displacements DOF
# - `Fp2_Vdot`: second node rotational equilibrium equations and sectional velocities DOF
# - `Fp2_Ωdot`: second node rotational equilibrium equations and sectional angular velocities DOF
# - `FV_udot`: midpoint linear velocity-displacement compatibility equations and translational displacements DOF
# - `FΩ_pdot`: midpoint angular velocity-rotation compatibility equations and rotational displacements DOF
# - `Fχ_Vdot`: aerodynamic equations and sectional velocities DOF
# - `Fχ_Ωdot`: aerodynamic equations and sectional angular velocities DOF
# - `Fχ_χdot`: aerodynamic equations and aerodynamic states DOF
#
struct InertiaIndicesMaps{N3,NN}

    Fu1_pdot::SVector{9,Int}
    Fu1_Vdot::SVector{9,Int}
    Fu1_Ωdot::SVector{9,Int}

    Fp1_pdot::SVector{9,Int}
    Fp1_Vdot::SVector{9,Int}
    Fp1_Ωdot::SVector{9,Int}

    Fu2_pdot::SVector{9,Int}
    Fu2_Vdot::SVector{9,Int}
    Fu2_Ωdot::SVector{9,Int}

    Fp2_pdot::SVector{9,Int}
    Fp2_Vdot::SVector{9,Int}
    Fp2_Ωdot::SVector{9,Int}

    FV_udot::SVector{9,Int}
    FΩ_pdot::SVector{9,Int}

    Fχ_Vdot::SVector{N3,Int}
    Fχ_Ωdot::SVector{N3,Int}
    Fχ_χdot::SVector{NN,Int}

    function InertiaIndicesMaps{N3,NN}(Fu1_pdot,Fu1_Vdot,Fu1_Ωdot,Fp1_pdot,Fp1_Vdot,Fp1_Ωdot,Fu2_pdot,Fu2_Vdot,Fu2_Ωdot,Fp2_pdot,Fp2_Vdot,Fp2_Ωdot,FV_udot,FΩ_pdot,Fχ_Vdot,Fχ_Ωdot,Fχ_χdot) where {N3,NN}
        new{N3,NN}(
            Fu1_pdot,Fu1_Vdot,Fu1_Ωdot,
            Fp1_pdot,Fp1_Vdot,Fp1_Ωdot,
            Fu2_pdot,Fu2_Vdot,Fu2_Ωdot,
            Fp2_pdot,Fp2_Vdot,Fp2_Ωdot,
            FV_udot,FΩ_pdot,
            Fχ_Vdot,Fχ_Ωdot,Fχ_χdot
        )
    end
    
end
InertiaIndicesMaps{N}() where {N} = 
    InertiaIndicesMaps{3*N,N*N}(
        ntuple(_ -> zeros(SVector{9, Int}), 3)...,
        ntuple(_ -> zeros(SVector{9, Int}), 3)...,
        ntuple(_ -> zeros(SVector{9, Int}), 3)...,
        ntuple(_ -> zeros(SVector{9, Int}), 3)...,
        ntuple(_ -> zeros(SVector{9, Int}), 2)...,
        zeros(SVector{3*N, Int}), zeros(SVector{3*N, Int}),
        zeros(SVector{N*N, Int})
    )


# MidpointStates composite type
# Fields
# - `u`: displacements resolved in basis A
# - `p`: rotation parameters resolved in basis A
# - `F`: sectional forces resolved in basis B
# - `M`: sectional moments resolved in basis B
# - `V`: sectional velocities resolved in basis B
# - `Ω`: sectional angular velocities resolved in basis B
# - `γ`: force strains resolved in basis b
# - `κ`: moment strains resolved in basis b
# - `P`: linear momenta resolved in basis B
# - `H`: angular momenta resolved in basis B
# - `χ`: aerodynamic states
#
struct MidpointStates{N}
    
    u::MVector{3,Float64}
    p::MVector{3,Float64}
    F::MVector{3,Float64}
    M::MVector{3,Float64}
    V::MVector{3,Float64}
    Ω::MVector{3,Float64}

    γ::MVector{3,Float64}
    κ::MVector{3,Float64}
    P::MVector{3,Float64}
    H::MVector{3,Float64}

    χ::MVector{N,Float64}

end
MidpointStates{N}(surface) where {N} = MidpointStates{N}(
    ntuple(_ -> zeros(MVector{3,Float64}), 6)...,
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
    initial_χ(N,surface),
)
MidpointStates(element,N::Int) = begin
    s = element.state.midpoint
    MidpointStates{N}(
        s.u, s.p, s.F, s.M, s.V, s.Ω,
        s.γ, s.κ, s.P, s.H,
        initial_χ(N,element.aero.surface)
    )
end
MidpointNaNStates(element,n::Float64) = MidpointStates{length(element.state.midpoint.χ)}(
    ntuple(_ -> fill(n, MVector{3,Float64}), 6)...,
    ntuple(_ -> fill(n, MVector{3,Float64}), 4)...,
    fill(n, MVector{length(element.state.midpoint.χ)})
)
Base.copy(x::MidpointStates{N}) where {N} = copy_fields(x)


# NodalStates composite type
# Fields
# - `u1`: displacements of node 1 resolved in basis A
# - `u2`: displacements of node 2 resolved in basis A
# - `p1`: rotation parameters of node 1 resolved in basis A
# - `p2`: rotation parameters of node 2 resolved in basis A
# - `u1_b`: displacements of node 1 resolved in basis b
# - `u2_b`: displacements of node 2 resolved in basis b
# - `p1_b`: rotation parameters of node 1 resolved in basis b
# - `p2_b`: rotation parameters of node 2 resolved in basis b
# - `F1`: sectional forces of node 1 resolved in basis B
# - `F2`: sectional forces of node 2 resolved in basis B
# - `M1`: sectional moments of node 1 resolved in basis B
# - `M2`: sectional moments of node 2 resolved in basis B
# - `θ`: array with rotation angles of nodes 1 and 2
#
struct NodalStates
    
    u1::MVector{3,Float64}
    u2::MVector{3,Float64}
    p1::MVector{3,Float64}
    p2::MVector{3,Float64}

    u1_b::MVector{3,Float64}
    u2_b::MVector{3,Float64}
    p1_b::MVector{3,Float64}
    p2_b::MVector{3,Float64}

    F1::MVector{3,Float64}
    F2::MVector{3,Float64}
    M1::MVector{3,Float64}
    M2::MVector{3,Float64}

    θ::MVector{2,Float64}

end
NodalStates() = NodalStates(
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
    zeros(MVector{2,Float64}),
)
NodalNaNStates() = NodalStates(
    ntuple(_ -> fill(NaN, MVector{3,Float64}), 4)...,
    ntuple(_ -> fill(NaN, MVector{3,Float64}), 4)...,
    ntuple(_ -> fill(NaN, MVector{3,Float64}), 4)...,
    fill(NaN, MVector{2,Float64})
)
Base.copy(x::NodalStates) = copy_fields(x)


# StatesRates composite type
# Fields
# - `udot`: displacements time derivative
# - `pdot`: rotation parameters time derivative
# - `Vdot`: sectional velocities time derivative
# - `Ωdot`: sectional angular velocities time derivative
# - `Pdot`: linear momenta time derivative
# - `Hdot`: angular momenta time derivative
# - `χdot`: aerodynamic states time derivative
#
struct StatesRates{N}
    
    udot::MVector{3,Float64}
    pdot::MVector{3,Float64}
    Vdot::MVector{3,Float64}
    Ωdot::MVector{3,Float64}

    Pdot::MVector{3,Float64}
    Hdot::MVector{3,Float64}

    χdot::MVector{N,Float64}
    
end
StatesRates{N}() where {N} = StatesRates{N}(
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
    ntuple(_ -> zeros(MVector{3,Float64}), 2)...,
    zeros(MVector{N,Float64}),
)
StatesRates(element,N::Int) = begin
    s = element.state.rates
    StatesRates{N}(
        s.udot, s.pdot, s.Vdot, s.Ωdot,
        s.Pdot, s.Hdot,
        zeros(MVector{N,Float64})
    )
end
StatesRatesNaN(element,n::Float64) = StatesRates{length(element.state.rates.χdot)}(
    ntuple(_ -> fill(n, MVector{3,Float64}), 4)...,
    ntuple(_ -> fill(n, MVector{3,Float64}), 2)...,
    fill(n, MVector{length(element.state.rates.χdot)})
)
Base.copy(x::StatesRates{N}) where {N} = copy_fields(x)


# EquivStatesRates composite type
# Fields
# - `udotEquiv`: equivalent displacements time derivative at begin of time step
# - `pdotEquiv`: equivalent rotation parameters time derivative at begin of time step
# - `VdotEquiv`: equivalent sectional velocities time derivative at begin of time step
# - `ΩdotEquiv`: equivalent sectional angular velocities time derivative at begin of time step
# - `χdotEquiv`: equivalent aerodynamic states time derivative at begin of time step
#
struct EquivStatesRates{N}
    
    udotEquiv::MVector{3,Float64}
    pdotEquiv::MVector{3,Float64}
    VdotEquiv::MVector{3,Float64}
    ΩdotEquiv::MVector{3,Float64}

    χdotEquiv::MVector{N,Float64}

end
EquivStatesRates{N}() where {N} = EquivStatesRates{N}(
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
    zeros(MVector{N,Float64}),
)


# Kinematics of basis b composite type
# Fields
# - `v`: element's b basis linear velocity vector
# - `ω`: element's b basis angular velocity vector
# - `vdot`: element's b basis linear acceleration vector
# - `ωdot`: element's b basis angular acceleration vector
#
struct KinematicsBasisb

    v::MVector{3,Float64}
    ω::MVector{3,Float64}
    vdot::MVector{3,Float64}
    ωdot::MVector{3,Float64}

end
KinematicsBasisb() = KinematicsBasisb(
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
)


# RotationMatrices composite type
# Fields
# - `R`: rotation matrix from basis b to basis B
# - `HT`: transpose of tangent tensor from basis B to basis b
# - `HTinv`: inverse transpose of tangent tensor from basis B to basis b
# - `RR0`: rotation matrix from basis A to basis B
# - `RR0T`: transpose of rotation matrix from basis A to basis B
# - `RdotR0`: time derivative of rotation matrix from basis A to basis B
# - `R_p1`: derivative of R w.r.t. first rotation parameter
# - `R_p2`: derivative of R w.r.t. second rotation parameter
# - `R_p3`: derivative of R w.r.t. third rotation parameter
# - `RT_p1`: derivative of R' w.r.t. first rotation parameter
# - `RT_p2`: derivative of R' w.r.t. second rotation parameter
# - `RT_p3`: derivative of R' w.r.t. third rotation parameter
# - `HT_p1`: derivative of HT w.r.t. first rotation parameter
# - `HT_p2`: derivative of HT w.r.t. second rotation parameter
# - `HT_p3`: derivative of HT w.r.t. third rotation parameter
# - `HTinv_p1`: derivative of HTinv w.r.t. first rotation parameter
# - `HTinv_p2`: derivative of HTinv w.r.t. second rotation parameter
# - `HTinv_p3`: derivative of HTinv w.r.t. third rotation parameter
# - `Rdot_p1`: derivative of RdotR0 w.r.t. first rotation parameter
# - `Rdot_p2`: derivative of RdotR0 w.r.t. second rotation parameter
# - `Rdot_p3`: derivative of RdotR0 w.r.t. third rotation parameter
#
struct RotationMatrices
    
    R::MMatrix{3,3,Float64,9}
    RR0::MMatrix{3,3,Float64,9}
    RR0T::MMatrix{3,3,Float64,9}
    RdotR0::MMatrix{3,3,Float64,9}

    HT::MMatrix{3,3,Float64,9}
    HTinv::MMatrix{3,3,Float64,9}

    R_p1::MMatrix{3,3,Float64,9}
    R_p2::MMatrix{3,3,Float64,9}
    R_p3::MMatrix{3,3,Float64,9}
    RT_p1::MMatrix{3,3,Float64,9}
    RT_p2::MMatrix{3,3,Float64,9}
    RT_p3::MMatrix{3,3,Float64,9}

    HT_p1::MMatrix{3,3,Float64,9}
    HT_p2::MMatrix{3,3,Float64,9}
    HT_p3::MMatrix{3,3,Float64,9}
    HTinv_p1::MMatrix{3,3,Float64,9}
    HTinv_p2::MMatrix{3,3,Float64,9}
    HTinv_p3::MMatrix{3,3,Float64,9}

    Rdot_p1::MMatrix{3,3,Float64,9}
    Rdot_p2::MMatrix{3,3,Float64,9}
    Rdot_p3::MMatrix{3,3,Float64,9}

end
RotationMatrices(R0) = RotationMatrices(
    MMatrix{3,3,Float64,9}(I3), MMatrix{3,3,Float64,9}(R0), MMatrix{3,3,Float64,9}(R0'), zeros(MMatrix{3,3,Float64,9}),
    ntuple(_ -> MMatrix{3,3,Float64,9}(I3), 2)...,
    ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 6)...,
    ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 6)...,
    ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 3)...,
)


# LoadsResultants composite type
# All resolved in basis A
# Fields
# - `f1`: first node total (structural + aero) forces resultant
# - `f2`: second node total (structural + aero) forces resultant
# - `m1`: first node total (structural + aero) moments resultant
# - `m2`: second node total (structural + aero) moments resultant
#
struct LoadsResultants
    
    f1::MVector{3,Float64}
    f2::MVector{3,Float64}
    m1::MVector{3,Float64}
    m2::MVector{3,Float64}

end
LoadsResultants() = LoadsResultants(
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
)


# StructuralLoadsResultants composite type
# All resolved in basis A
# Fields
# - `f1_ext`: first node external forces (dead + follower) resultant
# - `f2_ext`: second node external forces (dead + follower) resultant
# - `m1_ext`: first node external moments (dead + follower) resultant
# - `m2_ext`: second node external moments (dead + follower) resultant
# - `ff1_A`: first node external follower forces resultant defined in basis A
# - `ff2_A`: second node external follower forces resultant defined in basis A
# - `mf1_A`: first node external follower moments resultant defined in basis A
# - `mf2_A`: second node external follower moments resultant defined in basis A
# - `ff1_b`: first node external follower forces resultant defined in basis b
# - `ff2_b`: second node external follower forces resultant defined in basis b
# - `mf1_b`: first node external follower moments resultant defined in basis b
# - `mf2_b`: second node external follower moments resultant defined in basis b
# - `f_g`: gravitational forces resultant (same for both nodes)
# - `m_g`: gravitational moments resultant (same for both nodes)
#
struct StructuralLoadsResultants

    f1_ext::MVector{3,Float64}
    f2_ext::MVector{3,Float64}
    m1_ext::MVector{3,Float64}
    m2_ext::MVector{3,Float64}

    ff1_A::MVector{3,Float64}
    ff2_A::MVector{3,Float64}
    mf1_A::MVector{3,Float64}
    mf2_A::MVector{3,Float64}

    ff1_b::MVector{3,Float64}
    ff2_b::MVector{3,Float64}
    mf1_b::MVector{3,Float64}
    mf2_b::MVector{3,Float64}

    f_g::MVector{3,Float64}
    m_g::MVector{3,Float64}

end
StructuralLoadsResultants() = StructuralLoadsResultants(
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
    ntuple(_ -> zeros(MVector{3,Float64}), 2)...,
)


# LoadsJacobians composite type
# Fields
# - `f1_p`: derivatives of first node total (structural + aerodynamic) forces resultant w.r.t. rotation
# - `f2_p`: derivatives of second node total (structural + aerodynamic) forces resultant w.r.t. rotation
# - `m1_p`: derivatives of first node total (structural + aerodynamic) moments resultant w.r.t. rotation
# - `m2_p`: derivatives of second node total (structural + aerodynamic) moments resultant w.r.t. rotation
#
struct LoadsJacobians

    f1_p::MMatrix{3,3,Float64,9}
    f2_p::MMatrix{3,3,Float64,9}
    m1_p::MMatrix{3,3,Float64,9}
    m2_p::MMatrix{3,3,Float64,9}

end
LoadsJacobians() = LoadsJacobians(
    ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 4)...,
)


# StructuralLoadsJacobians composite type
# Fields
# - `f1_ext_p`: derivatives of first node external forces resultant w.r.t. rotation parameters
# - `f2_ext_p`: derivatives of second node external forces resultant w.r.t. rotation parameters
# - `m1_ext_p`: derivatives of first node external moments resultant w.r.t. rotation parameters
# - `m2_ext_p`: derivatives of second node external moments resultant w.r.t. rotation parameters
# - `mg_p`: derivatives of gravitational moments resultant w.r.t. rotation parameters
#
struct StructuralLoadsJacobians

    f1_ext_p::MMatrix{3,3,Float64,9}
    f2_ext_p::MMatrix{3,3,Float64,9}
    m1_ext_p::MMatrix{3,3,Float64,9}
    m2_ext_p::MMatrix{3,3,Float64,9}

    mg_p::MMatrix{3,3,Float64,9}

end
StructuralLoadsJacobians() = StructuralLoadsJacobians(
    ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 4)...,
    zeros(MMatrix{3,3,Float64,9}),
)


# AeroLoadsJacobians composite type
# Fields
# - `f1χ_p`: aerodynamic force derivatives w.r.t. rotation parameters
# - `f2χ_p`: aerodynamic force derivatives w.r.t. rotation parameters
# - `m1χ_p`: aerodynamic moment derivatives w.r.t. rotation parameters
# - `m2χ_p`: aerodynamic moment derivatives w.r.t. rotation parameters
# - `f1χ_V`: aerodynamic force derivatives w.r.t. velocity
# - `f2χ_V`: aerodynamic force derivatives w.r.t. velocity
# - `m1χ_V`: aerodynamic moment derivatives w.r.t. velocity
# - `m2χ_V`: aerodynamic moment derivatives w.r.t. velocity
# - `f1χ_Ω`: aerodynamic force derivatives w.r.t. angular velocity
# - `f2χ_Ω`: aerodynamic force derivatives w.r.t. angular velocity
# - `m1χ_Ω`: aerodynamic moment derivatives w.r.t. angular velocity
# - `m2χ_Ω`: aerodynamic moment derivatives w.r.t. angular velocity
# - `f1χ_χ`: aerodynamic force derivatives w.r.t. aerodynamic states
# - `f2χ_χ`: aerodynamic force derivatives w.r.t. aerodynamic states
# - `m1χ_χ`: aerodynamic moment derivatives w.r.t. aerodynamic states
# - `m2χ_χ`: aerodynamic moment derivatives w.r.t. aerodynamic states
# - `f1χ_δ`: aerodynamic force derivatives w.r.t. control surface deflection
# - `f2χ_δ`: aerodynamic force derivatives w.r.t. control surface deflection
# - `m1χ_δ`: aerodynamic moment derivatives w.r.t. control surface deflection
# - `m2χ_δ`: aerodynamic moment derivatives w.r.t. control surface deflection

# - `f1χ_Vdot`: aerodynamic force derivatives w.r.t. velocity rate
# - `f2χ_Vdot`: aerodynamic force derivatives w.r.t. velocity rate
# - `m1χ_Vdot`: aerodynamic moment derivatives w.r.t. velocity rate
# - `m2χ_Vdot`: aerodynamic moment derivatives w.r.t. velocity rate
# - `f1χ_Ωdot`: aerodynamic force derivatives w.r.t. angular velocity rate
# - `f2χ_Ωdot`: aerodynamic force derivatives w.r.t. angular velocity rate
# - `m1χ_Ωdot`: aerodynamic moment derivatives w.r.t. angular velocity rate
# - `m2χ_Ωdot`: aerodynamic moment derivatives w.r.t. angular velocity rate
#
struct AeroLoadsJacobians{N,Nδ,N3,Nδ3}

    f1χ_p::MMatrix{3,3,Float64,9}
    f2χ_p::MMatrix{3,3,Float64,9}
    m1χ_p::MMatrix{3,3,Float64,9}
    m2χ_p::MMatrix{3,3,Float64,9}

    f1χ_V::MMatrix{3,3,Float64,9}
    f2χ_V::MMatrix{3,3,Float64,9}
    m1χ_V::MMatrix{3,3,Float64,9}
    m2χ_V::MMatrix{3,3,Float64,9}

    f1χ_Ω::MMatrix{3,3,Float64,9}
    f2χ_Ω::MMatrix{3,3,Float64,9}
    m1χ_Ω::MMatrix{3,3,Float64,9}
    m2χ_Ω::MMatrix{3,3,Float64,9}

    f1χ_χ::MMatrix{3,N,Float64,N3}
    f2χ_χ::MMatrix{3,N,Float64,N3}
    m1χ_χ::MMatrix{3,N,Float64,N3}
    m2χ_χ::MMatrix{3,N,Float64,N3}

    f1χ_δ::MMatrix{3,Nδ,Float64,Nδ3}
    f2χ_δ::MMatrix{3,Nδ,Float64,Nδ3}
    m1χ_δ::MMatrix{3,Nδ,Float64,Nδ3}
    m2χ_δ::MMatrix{3,Nδ,Float64,Nδ3}

    f1χ_Vdot::MMatrix{3,3,Float64,9}
    f2χ_Vdot::MMatrix{3,3,Float64,9}
    m1χ_Vdot::MMatrix{3,3,Float64,9}
    m2χ_Vdot::MMatrix{3,3,Float64,9}

    f1χ_Ωdot::MMatrix{3,3,Float64,9}
    f2χ_Ωdot::MMatrix{3,3,Float64,9}
    m1χ_Ωdot::MMatrix{3,3,Float64,9}
    m2χ_Ωdot::MMatrix{3,3,Float64,9}

    function AeroLoadsJacobians{N,Nδ,N3,Nδ3}(f1χ_p,f2χ_p,m1χ_p,m2χ_p,f1χ_V,f2χ_V,m1χ_V,m2χ_V,f1χ_Ω,f2χ_Ω,m1χ_Ω,m2χ_Ω,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ,f1χ_Vdot,f2χ_Vdot,m1χ_Vdot,m2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Ωdot,m2χ_Ωdot) where {N,Nδ,N3,Nδ3}
        new{N,Nδ,N3,Nδ3}(f1χ_p,f2χ_p,m1χ_p,m2χ_p,f1χ_V,f2χ_V,m1χ_V,m2χ_V,f1χ_Ω,f2χ_Ω,m1χ_Ω,m2χ_Ω,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,f1χ_δ,f2χ_δ,m1χ_δ,m2χ_δ,f1χ_Vdot,f2χ_Vdot,m1χ_Vdot,m2χ_Vdot,f1χ_Ωdot,f2χ_Ωdot,m1χ_Ωdot,m2χ_Ωdot)
    end

end
AeroLoadsJacobians{N,Nδ}() where {N,Nδ} =
    AeroLoadsJacobians{N,Nδ,3*N,3*Nδ}(
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 4)...,
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 4)...,
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 4)...,
        ntuple(_ -> zeros(MMatrix{3,N,Float64,3*N}), 4)...,
        ntuple(_ -> zeros(MMatrix{3,Nδ,Float64,3*Nδ}), 4)...,
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 4)...,
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 4)...,
    )


# Residuals composite type
# Fields
# - `F_u1`: first node translational equilibrium equations' residuals
# - `F_u2`: second node translational equilibrium equations' residuals
# - `F_p1`: first node rotational equilibrium equations' residuals
# - `F_p2`: second node rotational equilibrium equations' residuals
# - `F_F1`: first node translational compatibility equations' residuals
# - `F_F2`: second node translational compatibility equations' residuals
# - `F_M1`: first node rotational compatibility equations' residuals
# - `F_M2`: second node rotational compatibility equations' residuals
# - `F_V`: midpoint linear velocity-displacement compatibility equations' residuals
# - `F_Ω`: midpoint angular velocity-displacement compatibility equations' residuals
# - `F_χ`: aerodynamic equations' residuals
#
struct Residuals{N}

    F_u1::MVector{3,Float64}
    F_u2::MVector{3,Float64}

    F_p1::MVector{3,Float64}
    F_p2::MVector{3,Float64}

    F_F1::MVector{3,Float64}
    F_F2::MVector{3,Float64}

    F_M1::MVector{3,Float64}
    F_M2::MVector{3,Float64}

    F_V::MVector{3,Float64}

    F_Ω::MVector{3,Float64}

    F_χ::MVector{N,Float64}

end
Residuals{N}() where {N} =
    Residuals{N}(
        ntuple(_ -> zeros(MVector{3,Float64}), 2)...,
        ntuple(_ -> zeros(MVector{3,Float64}), 2)...,
        ntuple(_ -> zeros(MVector{3,Float64}), 2)...,
        ntuple(_ -> zeros(MVector{3,Float64}), 2)...,
        zeros(MVector{3,Float64}),
        zeros(MVector{3,Float64}),
        zeros(MVector{N,Float64})
    )


# Jacobians composite type
# Fields
# - `F_u1_p`: ∂F_u1/∂p
# - `F_u2_p`: ∂F_u2/∂p
# - `F_u1_F`: ∂F_u1/∂F
# - `F_u2_F`: ∂F_u2/∂F
# - `F_u1_V`: ∂F_u1/∂V
# - `F_u2_V`: ∂F_u2/∂V
# - `F_u1_Ω`: ∂F_u1/∂Ω
# - `F_u2_Ω`: ∂F_u2/∂Ω
# - `F_u1_χ`: ∂F_u1/∂χ
# - `F_u2_χ`: ∂F_u2/∂χ
# - `F_u1_δ`: ∂F_u1/∂δ
# - `F_u2_δ`: ∂F_u2/∂δ
# - `F_p1_p`: ∂F_p1/∂p
# - `F_p2_p`: ∂F_p2/∂p
# - `F_p1_F`: ∂F_p1/∂F
# - `F_p2_F`: ∂F_p2/∂F
# - `F_p1_M`: ∂F_p1/∂M
# - `F_p2_M`: ∂F_p2/∂M
# - `F_p1_V`: ∂F_p1/∂V
# - `F_p2_V`: ∂F_p2/∂V
# - `F_p1_Ω`: ∂F_p1/∂Ω
# - `F_p2_Ω`: ∂F_p2/∂Ω
# - `F_p1_χ`: ∂F_p1/∂χ
# - `F_p2_χ`: ∂F_p2/∂χ
# - `F_p1_δ`: ∂F_p1/∂δ
# - `F_p2_δ`: ∂F_p2/∂δ
# - `F_F1_u`: ∂F_F1/∂u
# - `F_F2_u`: ∂F_F2/∂u
# - `F_F1_p`: ∂F_F1/∂p
# - `F_F2_p`: ∂F_F2/∂p
# - `F_F1_F`: ∂F_F1/∂F
# - `F_F2_F`: ∂F_F2/∂F
# - `F_F1_M`: ∂F_F1/∂M
# - `F_F2_M`: ∂F_F2/∂M
# - `F_M1_p`: ∂F_M1/∂p
# - `F_M2_p`: ∂F_M2/∂p
# - `F_M1_F`: ∂F_M1/∂F
# - `F_M2_F`: ∂F_M2/∂F
# - `F_M1_M`: ∂F_M1/∂M
# - `F_M2_M`: ∂F_M2/∂M
# - `F_V_u`: ∂F_V/∂u
# - `F_V_p`: ∂F_V/∂p
# - `F_V_V`: ∂F_V/∂V
# - `F_Ω_p`: ∂F_Ω/∂p
# - `F_Ω_Ω`: ∂F_Ω/∂Ω
# - `F_χ_p`: ∂F_χ/∂p
# - `F_χ_V`: ∂F_χ/∂V
# - `F_χ_Ω`: ∂F_χ/∂Ω
# - `F_χ_χ`: ∂F_χ/∂χ
# - `F_u1_pdot`: ∂F_u1/∂pdot
# - `F_u2_pdot`: ∂F_u2/∂pdot
# - `F_u1_Vdot`: ∂F_u1/∂Vdot
# - `F_u2_Vdot`: ∂F_u2/∂Vdot
# - `F_u1_Ωdot`: ∂F_u1/∂Ωdot
# - `F_u2_Ωdot`: ∂F_u2/∂Ωdot
# - `F_p1_pdot`: ∂F_p1/∂pdot
# - `F_p2_pdot`: ∂F_p2/∂pdot
# - `F_p1_Vdot`: ∂F_p1/∂Vdot
# - `F_p2_Vdot`: ∂F_p2/∂Vdot
# - `F_p1_Ωdot`: ∂F_p1/∂Ωdot
# - `F_p2_Ωdot`: ∂F_p2/∂Ωdot
# - `F_V_udot`: ∂F_V/∂udot
# - `F_Ω_pdot`: ∂F_Ω/∂pdot
# - `F_χ_Vdot`: ∂F_χ/∂Vdot
# - `F_χ_Ωdot`: ∂F_χ/∂Ωdot
# - `F_χ_χdot`: ∂F_χ/∂χdot
#    
struct Jacobians{N,Nδ,N3,Nδ3,NN}

    F_u1_p::MMatrix{3,3,Float64,9}
    F_u2_p::MMatrix{3,3,Float64,9}
    F_u1_F::MMatrix{3,3,Float64,9}
    F_u2_F::MMatrix{3,3,Float64,9}
    F_u1_V::MMatrix{3,3,Float64,9}
    F_u2_V::MMatrix{3,3,Float64,9}
    F_u1_Ω::MMatrix{3,3,Float64,9}
    F_u2_Ω::MMatrix{3,3,Float64,9}
    F_u1_χ::MMatrix{3,N,Float64,N3}
    F_u2_χ::MMatrix{3,N,Float64,N3}
    F_u1_δ::MMatrix{3,Nδ,Float64,Nδ3}
    F_u2_δ::MMatrix{3,Nδ,Float64,Nδ3}

    F_p1_p::MMatrix{3,3,Float64,9}
    F_p2_p::MMatrix{3,3,Float64,9}
    F_p1_F::MMatrix{3,3,Float64,9}
    F_p2_F::MMatrix{3,3,Float64,9}
    F_p1_M::MMatrix{3,3,Float64,9}
    F_p2_M::MMatrix{3,3,Float64,9}
    F_p1_V::MMatrix{3,3,Float64,9}
    F_p2_V::MMatrix{3,3,Float64,9}
    F_p1_Ω::MMatrix{3,3,Float64,9}
    F_p2_Ω::MMatrix{3,3,Float64,9}
    F_p1_χ::MMatrix{3,N,Float64,N3}
    F_p2_χ::MMatrix{3,N,Float64,N3}
    F_p1_δ::MMatrix{3,Nδ,Float64,Nδ3}
    F_p2_δ::MMatrix{3,Nδ,Float64,Nδ3}

    F_F1_u::SMatrix{3,3,Float64,9}
    F_F2_u::SMatrix{3,3,Float64,9}
    F_F1_p::MMatrix{3,3,Float64,9}
    F_F2_p::MMatrix{3,3,Float64,9}
    F_F1_F::MMatrix{3,3,Float64,9}
    F_F2_F::MMatrix{3,3,Float64,9}
    F_F1_M::MMatrix{3,3,Float64,9}
    F_F2_M::MMatrix{3,3,Float64,9}

    F_M1_p::MMatrix{3,3,Float64,9}
    F_M2_p::MMatrix{3,3,Float64,9}
    F_M1_F::MMatrix{3,3,Float64,9}
    F_M2_F::MMatrix{3,3,Float64,9}
    F_M1_M::MMatrix{3,3,Float64,9}
    F_M2_M::MMatrix{3,3,Float64,9}

    F_V_u::MMatrix{3,3,Float64,9}
    F_V_p::MMatrix{3,3,Float64,9}
    F_V_V::MMatrix{3,3,Float64,9}

    F_Ω_p::MMatrix{3,3,Float64,9}
    F_Ω_Ω::SMatrix{3,3,Float64,9}

    F_χ_p::MMatrix{N,3,Float64,N3}
    F_χ_V::MMatrix{N,3,Float64,N3}
    F_χ_Ω::MMatrix{N,3,Float64,N3}
    F_χ_χ::MMatrix{N,N,Float64,NN}

    F_u1_pdot::MMatrix{3,3,Float64,9}
    F_u2_pdot::MMatrix{3,3,Float64,9}
    F_u1_Vdot::MMatrix{3,3,Float64,9}
    F_u2_Vdot::MMatrix{3,3,Float64,9}
    F_u1_Ωdot::MMatrix{3,3,Float64,9}
    F_u2_Ωdot::MMatrix{3,3,Float64,9}

    F_p1_pdot::MMatrix{3,3,Float64,9}
    F_p2_pdot::MMatrix{3,3,Float64,9}
    F_p1_Vdot::MMatrix{3,3,Float64,9}
    F_p2_Vdot::MMatrix{3,3,Float64,9}
    F_p1_Ωdot::MMatrix{3,3,Float64,9}
    F_p2_Ωdot::MMatrix{3,3,Float64,9}

    F_V_udot::SMatrix{3,3,Float64,9}

    F_Ω_pdot::MMatrix{3,3,Float64,9}

    F_χ_Vdot::MMatrix{N,3,Float64,N3}
    F_χ_Ωdot::MMatrix{N,3,Float64,N3}
    F_χ_χdot::SMatrix{N,N,Float64,NN}

    function Jacobians{N,Nδ,N3,Nδ3,NN}(F_u1_p,F_u2_p,F_u1_F,F_u2_F,F_u1_V,F_u2_V,F_u1_Ω,F_u2_Ω,F_u1_χ,F_u2_χ,F_u1_δ,F_u2_δ,F_p1_p,F_p2_p,F_p1_F,F_p2_F,F_p1_M,F_p2_M,F_p1_V,F_p2_V,F_p1_Ω,F_p2_Ω,F_p1_χ,F_p2_χ,F_p1_δ,F_p2_δ,F_F1_u,F_F2_u,F_F1_p,F_F2_p,F_F1_F,F_F2_F,F_F1_M,F_F2_M,F_M1_p,F_M2_p,F_M1_F,F_M2_F,F_M1_M,F_M2_M,F_V_u,F_V_p,F_V_V,F_Ω_p,F_Ω_Ω,F_χ_p,F_χ_V,F_χ_Ω,F_χ_χ,F_u1_pdot,F_u2_pdot,F_u1_Vdot,F_u2_Vdot,F_u1_Ωdot,F_u2_Ωdot,F_p1_pdot,F_p2_pdot,F_p1_Vdot,F_p2_Vdot,F_p1_Ωdot,F_p2_Ωdot,F_V_udot,F_Ω_pdot,F_χ_Vdot,F_χ_Ωdot,F_χ_χdot) where {N,Nδ,N3,Nδ3,NN}
        new{N,Nδ,N3,Nδ3,NN}(
            F_u1_p, F_u2_p, F_u1_F, F_u2_F, F_u1_V, F_u2_V, F_u1_Ω, F_u2_Ω, F_u1_χ, F_u2_χ, F_u1_δ, F_u2_δ,
            F_p1_p, F_p2_p, F_p1_F, F_p2_F, F_p1_M, F_p2_M, F_p1_V, F_p2_V, F_p1_Ω, F_p2_Ω, F_p1_χ, F_p2_χ, F_p1_δ, F_p2_δ,
            F_F1_u, F_F2_u, F_F1_p, F_F2_p, F_F1_F, F_F2_F, F_F1_M, F_F2_M,
            F_M1_p, F_M2_p, F_M1_F, F_M2_F, F_M1_M, F_M2_M,
            F_V_u, F_V_p, F_V_V,
            F_Ω_p, F_Ω_Ω,
            F_χ_p, F_χ_V, F_χ_Ω, F_χ_χ,
            F_u1_pdot, F_u2_pdot, F_u1_Vdot, F_u2_Vdot, F_u1_Ωdot, F_u2_Ωdot,
            F_p1_pdot, F_p2_pdot, F_p1_Vdot, F_p2_Vdot, F_p1_Ωdot, F_p2_Ωdot,
            F_V_udot,
            F_Ω_pdot,
            F_χ_Vdot, F_χ_Ωdot, F_χ_χdot
        )
    end

end
Jacobians{N,Nδ}(surface,geometry) where {N,Nδ} =
    Jacobians{N,Nδ,3*N,3*Nδ,N*N}(
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 8)..., ntuple(_ -> zeros(MMatrix{3,N,Float64,3*N}), 2)..., ntuple(_ -> zeros(MMatrix{3,Nδ,Float64,3*Nδ}), 2)...,
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 10)..., ntuple(_ -> zeros(MMatrix{3,N,Float64,3*N}), 2)..., ntuple(_ -> zeros(MMatrix{3,Nδ,Float64,3*Nδ}), 2)...,
        SMatrix{3,3,Float64,9}(I), SMatrix{3,3,Float64,9}(-I), ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 6)...,
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 6)...,
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 3)...,
        zeros(MMatrix{3,3,Float64,9}), SMatrix{3,3,Float64,9}(I),
        ntuple(_ -> zeros(MMatrix{N,3,Float64,3*N}), 3)..., initial_F_χ_χ(N,surface),
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 6)...,
        ntuple(_ -> zeros(MMatrix{3,3,Float64,9}), 6)...,
        SMatrix{3,3,Float64,9}(-I),
        zeros(MMatrix{3,3,Float64,9}),
        initial_F_χ_Vdot(N,surface), initial_F_χ_Ωdot(N,surface,geometry), SMatrix{N,N,Float64,N*N}(I),
    )


# Element state composite type
# Fields
# - `midpoint`: midpoint states
# - `nodal`: nodal states
# - `rates`: midpoint states' rates
#
struct ElementState{N}
    
    midpoint::MidpointStates{N}
    nodal::NodalStates
    rates::StatesRates{N}

end
ElementState{N}(surface) where {N} = ElementState(
    MidpointStates{N}(surface),
    NodalStates(),
    StatesRates{N}(),
)
ElementState(element,N::Int) = 
    ElementState{N}(
        MidpointStates(element,N),
        element.state.nodal,
        StatesRates(element,N)
    )
ElementNaNState(element,n::Float64) = ElementState(
    MidpointNaNStates(element,n),
    NodalNaNStates(),
    StatesRatesNaN(element,n),
)
Base.copy(x::ElementState{N}) where {N} = copy_fields(x)


# Element state cache composite type
# Fields
# - `equivStatesRates`: equivalent states' rates at begin of time step
# - `kinematicsBasisb`: velocities and accelerations of basis b
# - `rotationMatrices`: rotation matrices
# - `loadsResultants`: total loads resultants
# - `structuralLoadsResultants`: structural loads resultants
# - `loadsJacobians`: total loads jacobians
# - `structuralLoadsJacobians`: structural loads jacobians
# - `aeroLoadsJacobians`: aerodynamic loads jacobians
# - `residuals`: equations' residuals
# - `jacobians`: jacobians (Jacobian and inertia matrices)
#
struct ElementStateCache{N,Nδ,N3,Nδ3,NN}

    equivStatesRates::EquivStatesRates{N}
    kinematicsBasisb::KinematicsBasisb
    rotationMatrices::RotationMatrices
    loadsResultants::LoadsResultants
    structuralLoadsResultants::StructuralLoadsResultants
    loadsJacobians::LoadsJacobians
    structuralLoadsJacobians::StructuralLoadsJacobians
    aeroLoadsJacobians::AeroLoadsJacobians{N,Nδ,N3,Nδ3}
    residuals::Residuals{N}
    jacobians::Jacobians{N,Nδ,N3,Nδ3,NN}

end
ElementStateCache{N,Nδ}(R0,surface,geometry) where {N,Nδ} = 
    ElementStateCache{N,Nδ,3*N,3*Nδ,N*N}(
        EquivStatesRates{N}(),
        KinematicsBasisb(),
        RotationMatrices(R0),
        LoadsResultants(),
        StructuralLoadsResultants(),
        LoadsJacobians(),
        StructuralLoadsJacobians(),
        AeroLoadsJacobians{N,Nδ}(),
        Residuals{N}(),
        Jacobians{N,Nδ}(surface,geometry),
    )

ElementStateCache(element,N::Int,Nδ::Int) =
    ElementStateCache{N,Nδ,3*N,3*Nδ,N*N}(
        EquivStatesRates{N}(),
        element.stateCache.kinematicsBasisb,
        element.stateCache.rotationMatrices,
        element.stateCache.loadsResultants,
        element.stateCache.structuralLoadsResultants,
        element.stateCache.loadsJacobians,
        element.stateCache.structuralLoadsJacobians,
        AeroLoadsJacobians{N,Nδ}(),
        Residuals{N}(),
        Jacobians{N,Nδ}(element.aero.surface,element.aero.geometry),
    )


# ModalMidpointStates composite type
# Same as MidpointStates, but without the aerodynamic states
# Fields
# - `u`: displacements resolved in basis A
# - `p`: rotation parameters resolved in basis A
# - `F`: sectional forces resolved in basis B
# - `M`: sectional moments resolved in basis B
# - `V`: sectional velocities resolved in basis B
# - `Ω`: sectional angular velocities resolved in basis B
# - `γ`: force strains resolved in basis b
# - `κ`: moment strains resolved in basis b
# - `P`: linear momenta resolved in basis B
# - `H`: angular momenta resolved in basis B
#
struct ModalMidpointStates

    u::MVector{3,Float64}
    p::MVector{3,Float64}
    F::MVector{3,Float64}
    M::MVector{3,Float64}
    V::MVector{3,Float64}
    Ω::MVector{3,Float64}

    γ::MVector{3,Float64}
    κ::MVector{3,Float64}
    P::MVector{3,Float64}
    H::MVector{3,Float64}

end
ModalMidpointStates() = ModalMidpointStates(
    ntuple(_ -> zeros(MVector{3,Float64}), 6)...,
    ntuple(_ -> zeros(MVector{3,Float64}), 4)...,
)
Base.copy(x::ModalMidpointStates) = copy_fields(x)


# ModalState composite type
# Fields
# - `midpoint`: midpoint states
# - `nodal`: nodal states
#
struct ModalState

    midpoint::ModalMidpointStates
    nodal::NodalStates

end
ModalState() = 
    ModalState(
        ModalMidpointStates(),
        NodalStates(),
    )
Base.copy(x::ModalState) = copy_fields(x)


# Element composite type
# Fields
# - `geometry`: element geometry data
# - `sectional`: element sectional (geometric/material) data
# - `aero`: aerodynamic properties associated with the element
# - `distLoads`: element distributed loads data
# - `assemblyData`: element assembly data
# - `assemblyIndices`: element assembly indices
# - `state`: element state data
# - `stateCache`: element state cache data
#
@with_kw mutable struct Element{N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,N3,Nc3,Ni3,Nn3,Nf3,Ng3,Nδ3,NcN,NiN,NnN,NfN,NgN,NN3,NNN,L,P,S1,S2,S3,F}

    # Deduced fields (from Beam definition)
    geometry::ElementGeometry
    sectional::ElementSectional
    aero::AeroProperties{N,Nc,Ni,Nn,Nf,Ng,Nδ,NN,N3,Nc3,Ni3,Nn3,Nf3,Ng3,NcN,NiN,NnN,NfN,NgN,NN3,NNN,L,P,S1,S2,S3,F}
    distLoads::ElementDistributedLoads = ElementDistributedLoads()

    # Assembly-dependent fields (to be set by Model and Problem)
    assemblyData::ElementAssemblyData
    assemblyIndices::ElementAssemblyIndices{N} = ElementAssemblyIndices{N}()
    jacobianMaps::JacobianIndicesMaps{N3,Nδ3,NN} = JacobianIndicesMaps{N,Nδ}()
    inertiaMaps::InertiaIndicesMaps{N3,NN} = InertiaIndicesMaps{N}()

    # State fields
    state::ElementState{N} = ElementState{N}(aero.surface)
    stateCache::ElementStateCache{N,Nδ,N3,Nδ3,NN} = ElementStateCache{N,Nδ}(geometry.R0,aero.surface,aero.geometry)

end


## PRIVATE API
# ------------------------------------------------------------------------------
# Element constructor
function Element(nElementsCreated,beamLength,kBeam,beamR0,beamS,beamI,normalizedNodalPositions,hingedNodes,hingedNodesDoF,surface)

    # ID on beam (local)
    localID = nElementsCreated + 1

    # Nodes' IDs on beam (local)
    nodesLocalIDs = [localID, localID+1]

    # Get element geometry data
    geometry = element_geometry(nodesLocalIDs,beamLength,kBeam,beamR0,normalizedNodalPositions,hingedNodes,hingedNodesDoF)

    # Get element sectional data
    sectional = element_sectional(beamS,beamI,localID)

    # Get aerodynamic properties and associated sizes and types
    aero,N,Nc,Ni,Nn,Nf,Ng,Nδ,S1,S2,S3,F = element_aero_properties(surface,geometry)

    # Initialize element assembly data
    assemblyData = ElementAssemblyData(localID,nodesLocalIDs)
    
    return Element{N,Nc,Ni,Nn,Nf,Ng,Nδ,N*N,3*N,3*Nc,3*Ni,3*Nn,3*Nf,3*Ng,3*Nδ,Nc*N,Ni*N,Nn*N,Nf*N,Ng*N,3*N*N,N*N*N,9+N+Nδ,12+Nc+Ni+Nn+Nf+Ng+N,S1,S2,S3,F}(
        geometry=geometry,
        sectional=sectional,
        aero=aero,
        assemblyData=assemblyData,
    )

end


# Element geometry
function element_geometry(nodesLocalIDs,beamLength,kBeam,beamR0,normalizedNodalPositions,hingedNodes,hingedNodesDoF)

    # Nodal arclength positions
    x1_n1,x1_n2 = beamLength * normalizedNodalPositions[nodesLocalIDs]

    # Length
    Δℓ = x1_n2 - x1_n1

    # Midpoint arclength position
    x1 = (x1_n1 + x1_n2)/2

    # Midpoint arclength position normalized by parent beam length
    x1_norm = x1/beamLength

    # Nodal arclength coordinates normalized by parent beam length
    x1_n1_norm,x1_n2_norm = normalizedNodalPositions[nodesLocalIDs]

    # Curvature vector
    k = kBeam

    # Midpoint rotation tensor from basis A to basis b and its transpose
    R0 = rotation_tensor_from_curvature(beamR0,k,x1)
    R0T = R0'

    # Nodal rotation tensors
    R0_n1 = rotation_tensor_from_curvature(beamR0,k,x1_n1)
    R0_n2 = rotation_tensor_from_curvature(beamR0,k,x1_n2)
    R0T_n1 = R0_n1'
    R0T_n2 = R0_n2'

    # Hinged nodes bit matrices
    hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat = get_hinged_nodes_matrices(hingedNodes,hingedNodesDoF,nodesLocalIDs)

    # Create element geometry data
    geometry = ElementGeometry(Δℓ,x1,x1_n1,x1_n2,x1_norm,x1_n1_norm,x1_n2_norm,k,R0,R0T,R0_n1,R0_n2,R0T_n1,R0T_n2,hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat)

    return geometry

end


# Gets the bit matrices resulting from hinged and not hinged nodal DoFs times the identity matrix
function get_hinged_nodes_matrices(hingedNodes,hingedNodesDoF,nodesLocalIDs)

    # Create 3x3 boolean identity matrix
    I3bit = Matrix(true*LinearAlgebra.I,3,3)

    # Initialize outputs
    hingedNode1Mat,hingedNode2Mat,notHingedNode1Mat = falses(3,3),falses(3,3),I3bit

    # Initialize flags of hinged DoFs on element's nodes
    hingedNode1DoF,hingedNode2DoF = falses(3),falses(3)

    # Loop over beam's hinged nodes
    for (hingedNode,hingedNodeDoF) in zip(hingedNodes,hingedNodesDoF)
        # Skip if the element does not contain the hinged node
        !(hingedNode in nodesLocalIDs) && continue
        # Flag of hinged DoFs on element's first node
        hingedNode1DoF = nodesLocalIDs[1] == hingedNode ? 
                         hingedNodeDoF : 
                         falses(3)
        # Flag of hinged DoFs on element's second node
        hingedNode2DoF = nodesLocalIDs[2] == hingedNode ? 
                         hingedNodeDoF : 
                         falses(3)
        # Skip once the hinged node is found                 
        break
    end

    # Hinged node bit matrices
    hingedNode1Mat = (hingedNode1DoF*hingedNode1DoF') .* I3bit
    hingedNode2Mat = (hingedNode2DoF*hingedNode2DoF') .* I3bit
    notHingedNode1Mat = .!hingedNode1Mat .* I3bit
    notHingedNode2Mat = .!hingedNode2Mat .* I3bit

    # Convert to correct type
    hingedNode1Mat = SMatrix{3,3,Bool,9}(hingedNode1Mat)
    notHingedNode1Mat = SMatrix{3,3,Bool,9}(notHingedNode1Mat)
    notHingedNode2Mat = SMatrix{3,3,Bool,9}(notHingedNode2Mat)

    return hingedNode1Mat,notHingedNode1Mat,notHingedNode2Mat

end


# Element sectional properties
function element_sectional(beamS,beamI,localID)

    # Sectional compliance (force-strain) matrix and submatrices
    C = length(beamS) == 1 ?
        inv(beamS[1]) :
        inv(beamS[localID])
    C11 = C[1:3,1:3]
    C12 = C[1:3,4:6]
    C21 = C[4:6,1:3]
    C22 = C[4:6,4:6]
    
    # Sectional inertia (momentum-velocity) matrix and submatrices
    I = length(beamI) == 1 ? 
        beamI[1] :
        beamI[localID]
    I11 = I[1:3,1:3]
    I12 = I[1:3,4:6]
    I21 = I[4:6,1:3]
    I22 = I[4:6,4:6]

    # Mass per unit length and sectional centroid position's skew-symmetric matrix with respect to reference point, resolved in basis b
    μ = I[1,1]
    ηtilde = μ > 0 ?
             I21/μ :
             Z3

    # Create element sectional data
    sectional = ElementSectional(C,I,C11,C12,C21,C22,I11,I12,I21,I22,ηtilde,μ)             

    return sectional

end


# Gets aerodynamic properties of the element and numbers of aerodynamic states
function element_aero_properties(surface,geometry)

    @unpack R0,x1,x1_norm,x1_n1_norm,x1_n2_norm = geometry

    return AeroProperties(surface,R0,x1,x1_norm,x1_n1_norm,x1_n2_norm)

end


# Update the distributed loads in the element's local coordinate
function update_element_distributed_loads!(element,loadType,loadFun)

    @unpack x1_n1,Δℓ = element.geometry

    # Function for arclength in terms of local coordinate ζ
    x1(ζ) = x1_n1 + Δℓ*ζ

    # New distributed loads function of the given type
    newLoadFun(ζ,t) = loadFun(x1(ζ),t)

    # Load type symbol
    loadTypeSymbol = Symbol(chopsuffix(String(loadType), "_x1t") * "_ζt")

    # Symbol for the flag indicating if the element has non-zero distributed loads of the given type
    hasLoadSymbol = Symbol("has_" * chopsuffix(String(loadType), "_of_x1t"))

    # Update fields of the element's distributed loads data
    setfield!(element.distLoads, loadTypeSymbol, newLoadFun)
    setfield!(element.distLoads, hasLoadSymbol, true)
    
end


# Adds the point inertia's properties to the element's sectional data
function add_point_inertia_to_element!(element,pointInertia)

    @unpack Δℓ = element.geometry
    @unpack I,C,C11,C12,C21,C22 = element.sectional
    @unpack mass,η,inertiaMatrix = pointInertia

    # Add to sectional inertia matrix 
    I += 1/Δℓ * [      mass*I3                -mass*tilde(η);
                 mass*tilde(η) inertiaMatrix-mass*tilde(η)^2]

    # Inertia submatrices
    I11 = I[1:3,1:3]
    I12 = I[1:3,4:6]
    I21 = I[4:6,1:3]
    I22 = I[4:6,4:6]

    # Mass per unit length and skew-symmetric matrix of sectional position of the element
    μ = I[1,1]
    ηtilde = μ > 0 ? 
             I21/μ :
             Z3

    # Update element sectional data
    setfield!(element, :sectional, ElementSectional(C,I,C11,C12,C21,C22,I11,I12,I21,I22,ηtilde,μ))

end


# Rebuilds an element with new type parameters, after changing the number of aerodynamic states (dependent on the presence of a gust and of flap trim states)
function rebuild_element(element,N,Nδ,withGust)

    # Recreate fields that depend on N, Nδ and gust presence
    aero,Nc,Ni,Nn,Nf,Ng,Nδ,S1,S2,S3,F = AeroProperties(element,N,Nδ,withGust)
    assemblyIndices                   = ElementAssemblyIndices{N}()
    jacobianMaps                      = JacobianIndicesMaps{N,Nδ}()
    inertiaMaps                       = InertiaIndicesMaps{N}()
    state                             = ElementState{N}(aero.surface)
    stateCache                        = ElementStateCache{N,Nδ}(element.geometry.R0,aero.surface,aero.geometry)

    # Reconstruct the element with new type parameters
    return Element{N,Nc,Ni,Nn,Nf,Ng,Nδ,N*N,3*N,3*Nc,3*Ni,3*Nn,3*Nf,3*Ng,3*Nδ,Nc*N,Ni*N,Nn*N,Nf*N,Ng*N,3*N*N,N*N*N,9+N+Nδ,12+Nc+Ni+Nn+Nf+Ng+N,S1,S2,S3,F}(
        geometry        = element.geometry,
        sectional       = element.sectional,
        aero            = aero,
        distLoads       = element.distLoads,
        assemblyData    = element.assemblyData,
        assemblyIndices = assemblyIndices,
        jacobianMaps    = jacobianMaps,
        inertiaMaps     = inertiaMaps,
        state           = state,
        stateCache      = stateCache,
    )

end