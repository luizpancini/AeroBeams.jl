"""
@with_kw mutable struct AttachedFlowParameters

    AttachedFlowParameters composite type

# Fields
- α₀N::Number
- ϵₙ::Number
- ϵₘ::Number
- cd₀::Number
- cdδ::Number
- cm₀::Number
- cmα::Number
- cmδ::Number
- cnα::Number
- cnδ::Number
"""
@with_kw mutable struct AttachedFlowParameters

    α₀N::Number
    ϵₙ::Number
    ϵₘ::Number
    cd₀::Number
    cdδ::Number
    cm₀::Number
    cmα::Number
    cmδ::Number
    cnα::Number
    cnδ::Number

    function AttachedFlowParameters(name::String; Re::Number=0,Ma::Number=0,flapSiteID::Int64=100)

        # Airfoil parameters' tables (as functions of Mach)
        if name in ["flatPlate","NACA0002","NACA0006"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.01,min(0.5,Ma))
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [ 0.01;   0.5]
            α₀NRng = π/180*[  0.0;   0.0]
            ϵₙRng =        [  0.7;   0.7]
            ϵₘRng =        [ 0.96;  0.96]
            cd₀Rng =  1e-2*[  0.0;   0.0]
            cm₀Rng =  1e-3*[  0.0;   0.0]
            cmαRng =       [  0.0;   0.0]
            cnαRng =  2π/β*[  1.0;   1.0]
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            elseif flapSiteID == 75
                cdδRng =  0.25*ones(length(MaRng))
                cmδRng = -0.25*ones(length(MaRng))
                cnδRng =  1.0*ones(length(MaRng))
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["NACA0012"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.035,min(0.3,Ma))
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [0.035; 0.072; 0.110; 0.185; 0.215; 0.25; 0.28;  0.3]
            α₀NRng = π/180*[  0.0;   0.0;   0.0;   0.0;   0.0;  0.0;  0.0;  0.0]
            ϵₙRng =        [  0.7;   0.7;   0.7;   0.7;   0.7;  0.7;  0.7;  0.7]
            ϵₘRng =        [ 0.96;  0.96;  0.96;  0.96;  0.96; 0.96; 0.96; 0.96]
            cd₀Rng =  1e-2*[  1.2;   1.2;   0.8;   0.5;   0.5;  0.5;  0.5;  0.5]
            cm₀Rng =  1e-3*[    0;   -14;    -5;    -5;    -5;   -5;   -5;   -5]
            cmαRng =       [  0.0;   0.0;   0.0;   0.0;   0.0;  0.0;  0.0;  0.0]
            cnαRng =  2π/β*[1.049; 0.972; 0.998; 0.995; 0.980; 1.00; 1.00; 1.00]
            if flapSiteID == 100
                cdδRng =   0.0*ones(length(MaRng))
                cmδRng = -0.25*ones(length(MaRng))
                cnδRng =   1.0*ones(length(MaRng))
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["NACA0018"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.062,min(0.15,Ma))
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [0.062; 0.080; 0.120; 0.150]
            α₀NRng = π/180*[  0.0;   0.0;   0.0;   0.0]
            ϵₙRng =        [  0.7;   0.7;   0.7;   0.7]
            ϵₘRng =        [ 0.96;  0.96;  0.96;  0.96]
            cd₀Rng =  1e-2*[  1.0;   0.8;   0.5;   0.6] 
            cm₀Rng =  1e-3*[ -1.0;   2.0;   2.0;   1.0]
            cmαRng =       [  0.0;   0.0;   0.0;   0.0]
            cnαRng =  2π/β*[  1.0;   1.0;   1.0;   1.0]
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            elseif flapSiteID == 75
                cdδRng =   0.0*ones(length(MaRng))
                cmδRng = -0.25*ones(length(MaRng))
                cnδRng =   1.0*ones(length(MaRng))
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["NACA23012A"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.001,min(0.3,Ma))
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  =       [  0.001;     0.3]
            α₀NRng = π/180*[    1.2;     1.2]
            ϵₙRng =        [    0.7;     0.7]
            ϵₘRng =        [   0.96;    0.96]
            cd₀Rng =  1e-2*[    0.3;     0.3] 
            cm₀Rng =  1e-2*[    5.0;     5.0]
            cmαRng =       [-0.0573; -0.0573]
            cnαRng =    2π*[   1.08;    1.08]
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            elseif flapSiteID == 75
                cdδRng =   0.0*ones(length(MaRng))
                cmδRng = -0.25*ones(length(MaRng))
                cnδRng =   1.0*ones(length(MaRng))
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["HeliosWingAirfoil"]
            # Bound Mach 
            Ma = max(0.001,min(0.3,Ma))
            # Mach-dependent parameters
            MaRng  =       [0.001; 0.3]
            α₀NRng = π/180*[  0.0; 0.0]
            ϵₙRng =        [  1.0; 1.0]
            ϵₘRng =        [  1.0; 1.0]
            cd₀Rng =  1e-2*[  1.0; 1.0]  
            cm₀Rng =  1e-2*[  2.5; 2.5]
            cmαRng =       [  0.0; 0.0]
            cnαRng =    2π*[  1.0; 1.0] 
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            elseif flapSiteID == 75
                cdδRng =   0.0*ones(length(MaRng))
                cmδRng = -0.25*ones(length(MaRng))
                cnδRng =   1.0*ones(length(MaRng))
            else
                error("Unavailable flap site ID")
            end     
        elseif name in ["HeliosPodAirfoil"]
            # Bound Mach 
            Ma = max(0.001,min(0.3,Ma))
            # Mach-dependent parameters
            MaRng  =       [0.001; 0.3]
            α₀NRng = π/180*[  0.0; 0.0]
            ϵₙRng =        [  1.0; 1.0]
            ϵₘRng =        [  1.0; 1.0]
            cd₀Rng =  1e-2*[  2.0; 2.0]  
            cm₀Rng =  1e-3*[  0.0; 0.0]
            cmαRng =       [  0.0; 0.0]
            cnαRng =       [  5.0; 5.0] 
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["BWBAirfoil"] 
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.01,min(0.3,Ma))
            β = sqrt(1-Ma^2)    
            # Mach-dependent parameters
            MaRng  =       [0.01; 0.3]
            α₀NRng = π/180*[ 0.0; 0.0]
            ϵₙRng =        [ 1.0; 1.0]
            ϵₘRng =        [ 1.0; 1.0]
            cd₀Rng =  1e-2*[ 1.0; 1.0]
            cm₀Rng =  1e-1*[ 1.0; 1.0]
            cmαRng =       [ 0.0; 0.0]
            cnαRng =  2π/β*[ 1.0; 1.0]
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            elseif flapSiteID == 75
                cdδRng =  0.0123*ones(length(MaRng))
                cmδRng = -0.4610*ones(length(MaRng))
                cnδRng =  3.5180*ones(length(MaRng))  
            else
                error("Unavailable flap site ID")
            end
        else
            error("Airfoil not listed")
        end

        # Interpolated values
        α₀N = interpolate(MaRng,α₀NRng,Ma)
        ϵₙ  = interpolate(MaRng,ϵₙRng,Ma)
        ϵₘ  = interpolate(MaRng,ϵₘRng,Ma)
        cd₀ = interpolate(MaRng,cd₀Rng,Ma)
        cdδ = interpolate(MaRng,cdδRng,Ma)
        cm₀ = interpolate(MaRng,cm₀Rng,Ma)
        cmα = interpolate(MaRng,cmαRng,Ma)
        cmδ = interpolate(MaRng,cmδRng,Ma)
        cnα = interpolate(MaRng,cnαRng,Ma)
        cnδ = interpolate(MaRng,cnδRng,Ma)

        return new(α₀N,ϵₙ,ϵₘ,cd₀,cdδ,cm₀,cmα,cmδ,cnα,cnδ)
    end

end


"""
@with_kw mutable struct SeparatedFlowParameters

    SeparatedFlowParameters composite type

# Fields
- α₀N::Number
- αds₀::Number
- αₛₛ::Number
- α1₀N::Number
- α1₀M::Number
- α1₀T::Number
- βσ1N::Number
- βσ1T::Number
- βσ2N::Number
- βS2Nlpr::Number
- βS2Tlpr::Number
- βS1Nu::Number
- βS1Mu::Number
- βS1Tu::Number
- βS1Nd::Number
- βS1Md::Number
- βS1Td::Number
- βS2Nu::Number
- βS2Mu::Number
- βS2Tu::Number
- βS2Nd::Number
- βS2Md::Number
- βS2Td::Number
- γLS::Number
- δα₀::Number
- δα₁::Number
- ϵₙ::Number
- ϵₘ::Number
- η::Number
- κ₀::Number
- κ₁::Number
- κ₂::Number
- κ₃::Number
- λ₁::Number
- λ₂::Number
- μv₂::Number
- ν₁::Number
- ν₂::Number
- ν₃::Number
- ν₄::Number
- ν₅::Number
- χu::Number
- χd::Number
- ξ::Number
- ζₐ::Number
- cd₀::Number
- cdδ::Number
- cm₀::Number
- cmδ::Number
- cnα::Number
- cnδ::Number
- dt::Number
- dm::Number
- E₀::Number
- E₁::Number
- f₀N::Number
- f₀M::Number
- f₀T::Number
- fbN::Number
- fbM::Number
- fbT::Number
- gᵥ::Number
- gᵥ₂::Number
- K₀::Number
- K₁::Number
- K₂::Number
- r₀::Number
- S1N::Number
- S1M::Number
- S1T::Number
- S2N::Number
- S2M::Number
- S2T::Number
- Ta::Number
- Tf::Number
- Tv::Number
- Tv₂::Number
- Vn₁::Number
- Vn₂::Number
- Vn₃::Number
- Vm::Number
- Vt::Number
- ztd::Number
- ztu::Number
- zm::Number
- λbWDiag::Matrix{Float64}
"""
@with_kw mutable struct SeparatedFlowParameters

    α₀N::Number
    αds₀::Number
    αₛₛ::Number
    α1₀N::Number
    α1₀M::Number
    α1₀T::Number
    βσ1N::Number
    βσ1T::Number
    βσ2N::Number
    βS2Nlpr::Number
    βS2Tlpr::Number
    βS1Nu::Number
    βS1Mu::Number
    βS1Tu::Number
    βS1Nd::Number
    βS1Md::Number
    βS1Td::Number
    βS2Nu::Number
    βS2Mu::Number
    βS2Tu::Number
    βS2Nd::Number
    βS2Md::Number
    βS2Td::Number
    γLS::Number
    δα₀::Number
    δα₁::Number
    ϵₙ::Number
    ϵₘ::Number
    η::Number
    κ₀::Number
    κ₁::Number
    κ₂::Number
    κ₃::Number
    λ₁::Number
    λ₂::Number
    μv₂::Number
    ν₁::Number
    ν₂::Number
    ν₃::Number
    ν₄::Number
    ν₅::Number
    χu::Number
    χd::Number
    ξ::Number
    ζₐ::Number
    cd₀::Number
    cdδ::Number
    cm₀::Number
    cmδ::Number
    cnα::Number
    cnδ::Number
    dt::Number
    dm::Number
    E₀::Number
    E₁::Number
    f₀N::Number
    f₀M::Number
    f₀T::Number
    fbN::Number
    fbM::Number
    fbT::Number
    gᵥ::Number
    gᵥ₂::Number
    K₀::Number
    K₁::Number
    K₂::Number
    r₀::Number
    S1N::Number
    S1M::Number
    S1T::Number
    S2N::Number
    S2M::Number
    S2T::Number
    Ta::Number
    Tf::Number
    Tv::Number
    Tv₂::Number
    Vn₁::Number
    Vn₂::Number
    Vn₃::Number
    Vm::Number
    Vt::Number
    ztd::Number
    ztu::Number
    zm::Number
    λbWDiag::Matrix{Float64}
    
    function SeparatedFlowParameters(name::String; Re::Number=0,Ma::Number=0,flapSiteID::Int64=100,U::Number=0,b::Number=0)

        # Airfoil parameters' tables (as functions of Mach)
        if name in ["NACA0012"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.035,min(0.3,Ma))
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [0.035; 0.072; 0.110; 0.185; 0.215; 0.25; 0.28;  0.3]
            α₀NRng = π/180*[  0.0;   0.0;   0.0;   0.0;   0.0;  0.0;  0.0;  0.0]
            αds₀Rng =π/180*[ 14.5;  17.6;  18.9;  18.8;  18.9; 17.5; 15.5; 15.0]
            αₛₛRng = π/180*[ 12.0;  13.5;  14.8;  16.1;  16.2; 15.9; 14.5; 13.6]
            α1₀NRng =π/180*[ 12.3;  13.7;  15.0;  16.3;  16.4; 15.9; 14.6; 13.8]
            α1₀MRng = π/180*[12.0;  13.7;  14.8;  16.2;  15.9;  5.7; 14.6; 13.8]
            α1₀TRng = π/180*[12.0;  13.6;  14.9;  16.0;  16.3; 16.1; 14.6; 13.6]
            βσ1NRng =       [ 3.0;  0.02;   0.7;   1.1;   0.7; 0.97;  0.0;  2.8]
            βσ1TRng =       [ 0.3;   0.3;   0.3;  -0.5;   0.7; -0.5;  1.5;  0.0]
            βσ2NRng =       [ 1.1;   1.8;   1.0;   1.0;   2.7;  2.1;  3.1;  2.0]
            βS2NlprRng =    [0.40;  0.47;  0.67;  0.78;  1.00; 0.36; 0.54; 0.38]
            βS2TlprRng =    [0.50;  0.50;  0.99;  0.83;  0.60; 0.99; 1.00; 0.88]
            βS1NuRng =      [1.09;  0.31;  1.62;  0.79;  1.80; 1.97;-0.63;-0.30]
            βS1MuRng =      [1.93; -0.80; -0.66;  1.74;  2.00; 0.09; 0.94; 1.78]
            βS1TuRng =      [1.97; -0.44;  0.17; -0.65;  0.62; 0.56; 1.98;-0.79]
            βS1NdRng =     [-0.63;  1.04;  0.51;  1.45; -0.27; 1.90; 1.98; 0.73]
            βS1MdRng =      [0.43;  1.10; -0.15;  0.47;  1.08; 0.59; 0.66; 0.27]
            βS1TdRng =      [1.98;  0.39; -0.10;  1.94;  1.93; 1.95; 1.91; 1.41]
            βS2NuRng =     [-0.04; -0.29; -0.54; -0.52; -0.36; 0.12;-0.63;-0.03]
            βS2MuRng =     [-0.80;  1.15; -0.79; -0.78; -0.67;-0.78; 2.00;-0.61]
            βS2TuRng =     [-0.70; -0.38; -0.71; -0.74; -0.80;-0.75;-0.80;-0.50]
            βS2NdRng =     [-0.80; -0.29;  0.14;  0.65;  1.93; 0.32; 1.84; 0.30]
            βS2MdRng =     [-0.61; -0.52; -0.48;  1.08;  0.49; 2.40; 2.47; 1.83]
            βS2TdRng =      [2.07;  1.31;  1.32; -0.32; -0.50; 1.96; 1.99; 1.34]
            γLSRng =        [0.20;  0.38;  0.35;  0.33;  0.66; 0.83; 0.90; 0.95]
            δα₀Rng=   π/180*[2.50;  0.50;  1.80;  2.80;  2.60; 1.80; 0.80; 0.80]
            δα₁Rng=   π/180*[0.95;  3.14;  0.50;  1.50;  1.80; 1.00;  3.5; 2.7]
            ϵₙRng =         [0.70;  0.70;  0.70;  0.70;  0.70; 0.70; 0.70; 0.70]
            ϵₘRng =         [0.96;  0.96;  0.96;  0.96;  0.96; 0.96; 0.96; 0.96]
            ηRng =         [0.058; 0.123; 0.050; 0.080; 0.120;0.090;0.100;0.100]
            κ₀Rng =         [2.26;   3.0;  2.46;  2.15;  2.12; 2.24; 2.35; 2.25]
            κ₁Rng =         [0.00;  0.07;  0.20;  0.03;  0.16; 0.01; 0.30; 0.00]
            κ₂Rng =         [0.03;  0.15;  0.18;  0.11;  0.13; 0.11; 0.10; 0.04]
            κ₃Rng =         [10.0;  4.75;  3.97;  1.79;  0.76; 1.37; 1.93; 1.00]
            λ₁Rng =           [10;    20;    25;    44;    20;   15;   10;    0]
            λ₂Rng =          [7.7;   4.0;   7.8;   7.8;   6.5;  7.7;  8.0;  8.0]
            μv₂Rng =        [0.00;  1.32;  1.09;  1.13;  0.80; 0.52; 0.50; 0.50]
            ν₁Rng =         [2.00;  1.80;  1.80;  1.77;  1.06; 0.43; 0.96; 0.57]
            ν₂Rng =         [0.30;  0.36;  1.00;  1.35;  1.36; 1.23; 1.23; 1.03]
            ν₃Rng =         [0.00;  0.00;  0.00;  0.00;  0.00; 0.00; 0.00; 0.00]
            ν₄Rng =         [1.57;  0.83;  1.03;  2.56;  1.50; 1.31;  1.0;  2.5]
            ν₅Rng =         [1.00;  1.00;  1.00;  1.00;  1.00; 1.00; 1.00; 1.00]
            χuRng =        [-0.00; -0.00; -0.00; -0.00; -0.00;-0.00;-0.00;-0.00]
            χdRng =        [-0.35; -0.35; -0.35;  0.23;  0.25; 0.15; 0.15; 0.11]
            ξRng =          [0.00;  0.05;  0.08;  0.24;  0.25; 0.37; 0.40; 0.44]
            ζₐRng =         [0.75;  0.31;  0.57;  0.47;  0.41; 0.62; 0.70; 0.75]
            cd₀Rng =       [0.012; 0.012; 0.008; 0.005; 0.005;0.005;0.005;0.005]
            cm₀Rng =   [-0.000;-0.014;-0.005;-0.005;-0.005;-0.005;-0.005;-0.005]
            cnαRng =  2π/β*[1.049; 0.972; 0.998; 0.995; 0.980; 1.00; 1.00; 1.00]
            dtRng =   π/180*[-1.9;   1.0;   0.4;  -0.7;   0.2; -0.8; -1.2; -0.7]
            dmRng =   π/180*[ 1.4;   1.9;   2.0;   1.9;   2.0;  2.3;  2.5;  0.9]
            E₀Rng =         [2.69;  1.22;  2.46;  1.84;  2.52; 2.00;  1.0; 2.57]
            E₁Rng =        [0.060; 0.091; 0.060; 0.076; 0.072;0.050;0.060;0.050]
            f₀NRng =       [0.020; 0.010; 0.020; 0.014; 0.020;0.016;0.010;0.003]
            f₀MRng =       [0.000; 0.000; 0.012; 0.010; 0.016;0.016;0.000;0.005]
            f₀TRng =       [0.000; 0.011; 0.020; 0.000; 0.002;0.020;0.020;0.010]
            fbNRng =        [0.71;  0.67;  0.67;  0.69;  0.68; 0.72; 0.70; 0.74]
            fbMRng =        [0.65;  0.67;  0.79;  0.71;  0.77; 0.78; 0.71; 0.70]
            fbTRng =        [0.75;  0.76;  0.76;  0.84;  0.80; 0.80; 0.81;0.825]
            gᵥRng =         [0.10; -0.05; -0.23; -0.22; -0.26;-0.27;-0.30;-0.30]
            gᵥ₂Rng =       [-1.00; -1.00; -1.00; -1.00; -1.00;-1.00;-1.00;-1.00]
            K₀Rng =        [0.015; 0.010; 0.000; 0.000; 0.008;0.010;0.008;0.000]
            K₁Rng =    [-0.141;-0.128;-0.107;-0.161;-0.131;-0.103;-0.105;-0.108]
            K₂Rng =        [0.010; 0.064; 0.041; 0.064; 0.052;0.042;0.040;0.045]
            r₀Rng =    1e-2*[1.58;  1.53;  1.58;  1.58;  1.20; 1.20; 1.20; 1.13]
            S1NRng =  π/180*[4.53;  5.00;  4.18;  3.31;  2.24; 3.80; 3.50; 2.83]
            S1MRng =  π/180*[4.47;  2.15;  4.00;  3.03;  4.00; 2.27; 2.57; 3.88]
            S1TRng =  π/180*[2.49;  3.37;  3.90;  3.08;  2.00; 3.02; 2.00; 2.93]
            S2NRng =  π/180*[1.00;  1.17;  1.02;  1.31;  1.08; 1.00; 1.24; 1.68]
            S2MRng =  π/180*[0.78;  1.02;  1.50;  2.50;  2.10; 1.27; 1.37; 1.40]
            S2TRng =  π/180*[0.97;  1.08;  1.47;  2.16;  1.00; 1.18; 1.50; 1.95]
            TaRng =         [3.66;  5.00;  5.97;  4.67;  3.05; 2.00; 2.66; 1.06]
            TfRng =         [3.02;  3.23;  3.00;  3.76;  3.00; 2.73; 2.80; 3.57]
            TvRng =         [3.34;  4.10;  3.72;  4.37;  4.46; 4.35; 4.30; 4.70]
            Tv₂Rng =        [4.06;  4.50;  5.19;  4.64;  4.50; 4.50; 4.55; 6.25]
            VmRng =         [0.33;  0.35;  0.39;  0.37;  0.37; 0.37; 0.37; 0.37]
            VtRng =         [0.10;  0.10;  0.10;  0.10;  0.00; 0.00; 0.00; 0.00]
            Vn₁Rng =        [1.08;  1.40;  1.57;  1.40;  1.16; 1.05; 1.00; 0.60]
            Vn₂Rng =        [1.21;  0.81;  1.00;  0.90;  0.75; 0.65; 0.51; 0.25]
            Vn₃Rng =        [0.60;  0.40;  0.50;  0.45;  0.38; 0.32; 0.25; 0.12]
            ztdRng =        [0.79;  1.18;  1.25;  1.13;  1.26; 1.23; 1.00; 1.00]
            ztuRng =        [1.00;  1.00;  1.00;  1.00;  1.00; 1.00; 1.00; 1.00]
            zmRng =         [2.31;  2.02;  2.50;  1.68;  1.47; 1.49; 1.48; 1.25]
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            elseif flapSiteID == 75
                cdδRng =   0.0*ones(length(MaRng))
                cmδRng = -0.25*ones(length(MaRng))
                cnδRng =   1.0*ones(length(MaRng))
            else
                error("Unavailable flap site ID")
            end
            # Fixed parameters
            λbWDiag = Diagonal([2.5; 0.8])
        elseif name in ["NACA0012-GU","NACA0015","NACA0015-s"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.078,min(0.155,Ma))
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [0.078;   0.117;   0.155]
            α₀NRng =     π/180*[0.6;     0.6;     0.5]    
            αds₀Rng =    π/180*[16.4;    18.8;    19.3]   
            αₛₛRng =     π/180*[14.0;    15.1;    14.9]   
            α1₀NRng =    π/180*[14.3;    15.1;    14.8]
            α1₀MRng =    π/180*[13.7;    14.9;    15.0]
            α1₀TRng =    π/180*[14.0;    15.0;    15.0]
            βσ1NRng =          [0.75;    0.25;    1.32]
            βσ1TRng =          [2.96;    2.08;    1.50]        
            βσ2NRng =          [0.50;    4.54;    4.17]         
            βS2NlprRng =        [0.20;    0.00;    0.20]   
            βS2TlprRng =        [0.62;    1.99;    0.38]    
            βS1NuRng =          [1.54;    1.10;    4.56]       
            βS1MuRng =          [3.12;    -0.68;   4.79]        
            βS1TuRng =          [0.17;    2.67;    4.99]       
            βS1NdRng =          [4.57;    4.93;    2.71]           
            βS1MdRng =          [3.51;    2.48;    4.99]      
            βS1TdRng =          [0.36;    3.77;    4.42]
            βS2NuRng =          [2.59;    4.80;    1.07]
            βS2MuRng =          [1.83;    -0.05;   -0.57]
            βS2TuRng =          [0.08;    -0.76;   -0.09]                 
            βS2NdRng =          [0.59;    -0.20;    2.20] 
            βS2MdRng =          [2.97;    2.44;    4.68]         
            βS2TdRng =          [3.89;    4.62;    4.98]       
            γLSRng =            [0.81;    0.81;    0.47]
            δα₀Rng= π/180*[1.43;    1.20;    1.32]    
            δα₁Rng= π/180*[2.97;    1.90;    1.59]      
            ϵₙRng =              [0.70;    0.70;    0.70]  
            ϵₘRng =              [0.96;    0.96;    0.96]   
            ηRng =                 [0.145;   0.003;   0.005] 
            κ₀Rng =             [1.84;    1.90;    1.82]
            κ₁Rng =             [0.72;    1.88;    0.93]
            κ₂Rng =             [0.20;    0.10;    0.20]
            κ₃Rng =             [4.18;    2.56;    3.68]
            λ₁Rng =            [5;       5;       5]
            λ₂Rng =            [0.5;     0.5;     0.5]
            μv₂Rng =               [2.43;    0.95;    1.34]
            ν₁Rng =                [1.70;    0.83;    1.64]    
            ν₂Rng =                [1.99;    2.52;    0.78]
            ν₃Rng =                [2.00;    2.00;    2.00]
            ν₄Rng =                [0.54;    0.82;    0.55]
            ν₅Rng =                [1.00;    1.00;    1.00]
            χuRng =               [-0.10;   0.08;    -0.10]
            χdRng =               [-0.06;   0.28;    -0.24]
            ξRng =                  [0.15;    0.02;    0.60]
            ζₐRng =              [0.27;    0.19;    0.41]
            cd₀Rng =           1e-2*[1.2;     0.8;     0.4]  
            cm₀Rng =           1e-3*[-1.0;    -1.2;    -2.0]
            cnαRng = 2*π/β*[0.988;   0.980;   1.016]  
            dtRng =         π/180*[2.2;     2.8;     0.6]    
            dmRng =         π/180*[2.9;     1.3;     1.8] 
            E₀Rng =                  [2.62;    0.85;    2.70]  
            E₁Rng =                  [0.055;   0.055;   0.053]
            f₀NRng =                [0.050;   0.042;   0.050]
            f₀MRng =                [0.027;   0.037;   0.027]
            f₀TRng =                [0.005;   0.008;   0.000]
            fbNRng =                [0.68;    0.70;    0.74]
            fbMRng =                [0.69;    0.77;    0.75]
            fbTRng =                [0.78;    0.80;    0.82]
            gᵥRng =                 [-0.23;   -0.30;   -0.33] 
            gᵥ₂Rng =                [-0.30;   -0.39;   0.50]
            K₀Rng =             1e-3*[1.6;     1.0;     1.0]
            K₁Rng =                  [-0.140;  -0.141;  -0.154]
            K₂Rng =                  [0.046;   0.049;   0.047]   
            r₀Rng =             1e-2*[1.37;    1.50;    1.30]  
            S1NRng =         π/180*[5.00;    3.14;    2.56]
            S1MRng =         π/180*[2.84;    2.25;    2.22]
            S1TRng =         π/180*[3.70;    2.50;    1.80]
            S2NRng =         π/180*[0.76;    0.76;    0.60]
            S2MRng =         π/180*[0.61;    0.78;    0.60]
            S2TRng =         π/180*[0.60;    0.67;    0.68]
            TaRng =                  [5.69;    4.82;    5.42]    
            TfRng =                  [3.14;    4.01;    3.64]    
            TvRng =                  [4.28;    4.72;    4.18] 
            Tv₂Rng =                 [4.91;    4.19;    3.66]
            VmRng =                  [0.35;    0.36;    0.38]
            VtRng =                  [0.16;    0.24;    0.23]
            Vn₁Rng =                 [2.49;    1.74;    1.10]   
            Vn₂Rng =                 [0.48;    0.46;    0.60]
            Vn₃Rng =                 [0.32;    1.02;    0.58] 
            ztdRng =               [1.00;    0.55;    1.35]
            ztuRng =               [0.00;    0.00;    0.00]
            zmRng =                [1.59;    1.65;    1.84]
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            elseif flapSiteID == 75
                cdδRng =   0.0*ones(length(MaRng))
                cmδRng = -0.25*ones(length(MaRng))
                cnδRng =   1.0*ones(length(MaRng))
            else
                error("Unavailable flap site ID")
            end
            # Fixed parameters
            λbWDiag = Diagonal([2.5; 0.8])
        elseif name in ["NACA0018"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.062,min(0.15,Ma))
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [0.062; 0.080; 0.120; 0.150]
            α₀NRng =     π/180*[0.0;    0.0;     0.0;     0.0]    
            αds₀Rng =    π/180*[16.8;   18.0;    18.1;    18.2]   
            αₛₛRng =     π/180*[15.0;   14.8;    14.7;    14.7]   
            α1₀NRng =    π/180*[14.5;   14.4;    14.7;    14.8]
            α1₀MRng =    π/180*[14.5;   14.5;    14.8;    14.8]
            α1₀TRng =    π/180*[14.9;   14.8;    14.8;    14.8]
            βσ1NRng =          [2.98;   1.20;    0.90;    1.77]
            βσ1TRng =          [1.44;   2.42;    1.01;    2.88]        
            βσ2NRng =          [1.78;   2.37;    3.02;    4.83]         
            βS2NlprRng =        [1.92;   1.91;    1.15;    0.57]   
            βS2TlprRng =        [2.00;   1.81;    2.00;    1.97]    
            βS1NuRng =          [0.86;   1.83;    1.55;    1.90]       
            βS1MuRng =          [-0.71;  1.07;    0.99;    -0.32]        
            βS1TuRng =          [-0.23;  1.69;    -0.13;   -0.62]       
            βS1NdRng =          [3.13;   0.88;    2.61;    0.98]           
            βS1MdRng =          [-0.31;  -0.10;   0.17;    0.05]      
            βS1TdRng =          [4.50;   5.00;    3.96;    4.12]
            βS2NuRng =          [1.30;   1.44;    2.16;    0.40]
            βS2MuRng =          [4.83;   2.09;    3.70;    3.90]
            βS2TuRng =          [-0.27;  -0.80;   -0.67;   -0.62]                 
            βS2NdRng =          [-0.36;  0.07;    0.63;    1.21] 
            βS2MdRng =          [2.09;   3.83;    4.65;    4.86]         
            βS2TdRng =          [4.88;   2.36;    3.69;    4.20]       
            γLSRng =            [1.00;   1.00;    1.00;    1.00]
            δα₀Rng= π/180*[0.00;   0.00;    0.00;    0.00]    
            δα₁Rng= π/180*[2.37;   1.20;    2.15;    1.86]      
            ϵₙRng =              [0.70;   0.70;    0.70;    0.70]  
            ϵₘRng =              [0.96;   0.96;    0.96;    0.96]   
            ηRng =                 [0.010;  0.030;   0.047;   0.048] 
            κ₀Rng =             [1.95;   1.88;    1.98;    1.98]
            κ₁Rng =             [1.04;   2.00;    2.86;    2.84]
            κ₂Rng =             [0.26;   0.22;    0.12;    0.13]
            κ₃Rng =             [2.77;   1.66;    0.46;    1.30]
            λ₁Rng =            [46;     68;      61;      70]
            λ₂Rng =            [0.51;   0.78;    0.51;    0.51]
            μv₂Rng =               [0.20;   0.00;    0.05;    0.99]
            ν₁Rng =                [1.97;   0.93;    1.53;    0.81]    
            ν₂Rng =                [1.89;   1.66;    1.80;    2.39]
            ν₃Rng =                [1.47;   1.34;    2.25;    2.78]
            ν₄Rng =                [1.43;   1.00;    1.56;    1.00]
            ν₅Rng =                [3.00;   2.08;    3.00;    3.00]
            χuRng =               [-0.02;  0.16;    -0.18;   0.05]
            χdRng =               [0.00;   0.19;    0.28;    0.30]
            ξRng =                  [0.71;   0.08;    0.04;    0.20]
            ζₐRng =              [0.46;   0.52;    0.52;    0.42]
            cd₀Rng =           1e-2*[1.0;    0.8;     0.5;     0.6]  
            cm₀Rng =           1e-3*[-1.0;   2.0;     2.0;     1.0]
            cnαRng = 2*π/β*[1.0;    1.0;     1.0;     1.0]  
            dtRng =         π/180*[0.0;    1.4;     1.8;     1.8]    
            dmRng =         π/180*[-3.4;   1.4;     1.0;     1.3] 
            E₀Rng =                  [2.15;   1.91;    1.54;    1.88]  
            E₁Rng =                  [0.080;  0.071;   0.070;   0.070]
            f₀NRng =                [0.029;  0.028;   0.020;   0.020]
            f₀MRng =                [0.013;  0.016;   0.011;   0.011]
            f₀TRng =                [0.001;  0.009;   0.006;   0.001]
            fbNRng =                [0.51;   0.54;    0.58;    0.63]
            fbMRng =                [0.60;   0.62;    0.60;    0.58]
            fbTRng =                [0.77;   0.80;    0.82;    0.82]
            gᵥRng =                 [0.24;   0.30;    0.21;    0.16]
            gᵥ₂Rng =                [0.10;   0.19;    0.15;    0.10]
            K₀Rng =                  [0.009;  0.010;   0.008;   0.008]
            K₁Rng =                  [-0.193; -0.190;  -0.154;  -0.159]
            K₂Rng =                  [0.063;  0.051;   0.047;   0.047]   
            r₀Rng =             1e-2*[1.38;   1.36;    1.42;    1.25]  
            S1NRng =         π/180*[6.95;   6.00;    4.86;    4.11]
            S1MRng =         π/180*[3.16;   4.15;    3.25;    3.77]
            S1TRng =         π/180*[5.30;   2.94;    2.90;    2.88]
            S2NRng =         π/180*[1.29;   1.27;    1.50;    1.41]
            S2MRng =         π/180*[6.47;   6.45;    6.49;    6.49]
            S2TRng =         π/180*[4.23;   5.41;    5.42;    5.56]
            TaRng =                  [4.12;   4.75;    3.68;    4.00]    
            TfRng =                  [3.76;   3.50;    3.80;    3.91]    
            TvRng =                  [4.00;   3.93;    5.26;    4.94] 
            Tv₂Rng =                 [3.88;   4.50;    3.82;    3.66]
            VmRng =                  [0.40;   0.60;    0.51;    0.50]
            VtRng =                  [0.35;   0.27;    0.29;    0.30]
            Vn₁Rng =                 [0.79;   0.53;    0.71;    0.80]   
            Vn₂Rng =                 [0.19;   0.43;    0.59;    0.55]
            Vn₃Rng =                 [0.08;   0.06;    0.20;    0.23]
            ztdRng =               [0.35;   0.50;    0.30;    0.30]
            ztuRng =               [0.00;   0.00;    0.00;    0.00]
            zmRng =                [2.30;   3.84;    2.83;    2.97]
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            elseif flapSiteID == 75
                cdδRng =   0.0*ones(length(MaRng))
                cmδRng = -0.25*ones(length(MaRng))
                cnδRng =   1.0*ones(length(MaRng))
            else
                error("Unavailable flap site ID")
            end
            # Fixed parameters
            λbWDiag = Diagonal([1.0; 1.0])
        elseif name in ["NACA23012A"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.001,min(0.3,Ma))
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  =       [  0.001;     0.3]
            α₀NRng =     π/180*[1.2;      1.2]
            αds₀Rng =    π/180*[17.7;     17.7]   
            αₛₛRng =     π/180*[14.4;     14.4]   
            α1₀NRng =    π/180*[14.4;     14.4]
            α1₀MRng =    π/180*[13.7;     13.7]
            α1₀TRng =    π/180*[14.1;     14.1]
            βσ1NRng =          [0.17;     0.17]
            βσ1TRng =          [0.29;     0.29]        
            βσ2NRng =          [2.47;     2.47]         
            βS2NlprRng =        [0.02;     0.02]   
            βS2TlprRng =        [1.64;     1.64]    
            βS1NuRng =          [1.43;     1.43]       
            βS1MuRng =          [2.31;     2.31]        
            βS1TuRng =          [0.97;     0.97]       
            βS1NdRng =          [5.82;     5.82]           
            βS1MdRng =          [5.67;     5.67]      
            βS1TdRng =          [4.02;     4.02]
            βS2NuRng =          [0.58;     0.58]         
            βS2MuRng =          [2.98;     2.98]          
            βS2TuRng =          [-0.80;    -0.80]                 
            βS2NdRng =          [0.77;     0.77]          
            βS2MdRng =          [1.15;     1.15]         
            βS2TdRng =          [4.50;     4.50]       
            γLSRng =            [0.73;     0.73]
            δα₀Rng= π/180*[0.00;     0.00]    
            δα₁Rng= π/180*[3.47;     3.47]      
            ϵₙRng =              [0.70;     0.70]  
            ϵₘRng =              [0.96;     0.96]   
            ηRng =                 [0.002;    0.002] 
            κ₀Rng =             [1.85;     1.85]
            κ₁Rng =             [1.51;     1.51]
            κ₂Rng =             [0.04;     0.04]
            κ₃Rng =             [0.35;     0.35]
            λ₁Rng =            [25.4;     25.4]
            λ₂Rng =            [0.62;     0.62]
            μv₂Rng =               [1.28;     1.28]
            ν₁Rng =                [1.64;     1.64]    
            ν₂Rng =                [0.55;     0.55]
            ν₃Rng =                [1.41;     1.41]
            ν₄Rng =                [1.74;     1.74]
            ν₅Rng =                [0.65;     0.65]
            χuRng =               [-0.17;    -0.17]
            χdRng =               [0.14;     0.14]
            ξRng =                  [0.00;     0.00]
            ζₐRng =              [0.20;     0.20]
            cd₀Rng =                [0.003;    0.003]  
            cm₀Rng =                [0.050;    0.050]
            cnαRng =      2*π*[1.08;     1.08]  
            dtRng =         π/180*[1.8;      1.8]    
            dmRng =         π/180*[-1.9;     -1.9] 
            E₀Rng =                  [1.76;     1.76]  
            E₁Rng =                  [0.057;    0.057]
            f₀NRng =                [0.050;    0.050]
            f₀MRng =                [0.033;    0.033]
            f₀TRng =                [0.020;    0.020]
            fbNRng =                [0.70;     0.70]
            fbMRng =                [0.70;     0.70]
            fbTRng =                [0.80;     0.80]
            gᵥRng =                 [0.18;     0.18] 
            gᵥ₂Rng =                [-0.48;    -0.48]
            K₀Rng =                  [-0.005;   -0.005]  
            K₁Rng =                  [-0.17;    -0.17] 
            K₂Rng =                  [0.015;    0.015]   
            r₀Rng =             1e-2*[1.62;     1.62]  
            S1NRng =         π/180*[2.50;     2.50]
            S1MRng =         π/180*[1.89;     1.89]
            S1TRng =         π/180*[1.27;     1.27]
            S2NRng =         π/180*[2.05;     2.05]
            S2MRng =         π/180*[4.70;     4.70]
            S2TRng =         π/180*[4.74;     4.74]
            TaRng =                  [3.74;     3.74]    
            TfRng =                  [3.75;     3.75]    
            TvRng =                  [4.99;     4.99] 
            Tv₂Rng =                 [4.67;     4.67]
            VmRng =                  [0.43;     0.43]
            VtRng =                  [0.33;     0.33]
            Vn₁Rng =                 [0.98;     0.98]   
            Vn₂Rng =                 [0.54;     0.54]
            Vn₃Rng =                 [0.13;     0.13]
            ztdRng =               [1.38;     1.38]
            ztuRng =               [0.10;     0.10]
            zmRng =                [1.72;     1.72]
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            elseif flapSiteID == 75
                cdδRng =   0.0*ones(length(MaRng))
                cmδRng = -0.25*ones(length(MaRng))
                cnδRng =   1.0*ones(length(MaRng))
            else
                error("Unavailable flap site ID")
            end
            # Fixed parameters
            λbWDiag = Diagonal([1.0; 1.0])
        elseif name in ["flatPlate","NACA0002","NACA0006","HeliosWingAirfoil","HeliosPodAirfoil","BWBAirfoil"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.001,min(0.3,Ma))
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  =       [  0.001;     0.3]
            α₀NRng =     π/180*[0.0;     0.0]
            αds₀Rng =   π/180*[14.5;    15.0]   
            αₛₛRng =     π/180*[12.0;    13.6]   
            α1₀NRng =   π/180*[12.3;    13.8]
            α1₀MRng =   π/180*[12.3;    13.8]
            α1₀TRng =   π/180*[12.3;    13.8]
            βσ1NRng =         [ 0.0;     0.0]
            βσ1TRng =         [ 0.0;     0.0]        
            βσ2NRng =         [ 0.0;     0.0]         
            βS2NlprRng =      [ 0.0;     0.0]   
            βS2TlprRng =      [ 0.0;     0.0]    
            βS1NuRng =        [ 0.0;     0.0]       
            βS1MuRng =        [ 0.0;     0.0]        
            βS1TuRng =        [ 0.0;     0.0]       
            βS1NdRng =        [ 0.0;     0.0]           
            βS1MdRng =        [ 0.0;     0.0]      
            βS1TdRng =        [ 0.0;     0.0]
            βS2NuRng =        [ 0.0;     0.0]         
            βS2MuRng =        [ 0.0;     0.0]          
            βS2TuRng =        [ 0.0;     0.0]                 
            βS2NdRng =        [ 0.0;     0.0]          
            βS2MdRng =        [ 0.0;     0.0]         
            βS2TdRng =        [ 0.0;     0.0]       
            γLSRng =          [ 1.0;     1.0]
            δα₀Rng=      π/180*[0.0;     0.0]    
            δα₁Rng=      π/180*[0.0;     0.0]      
            ϵₙRng =           [0.70;    0.70]  
            ϵₘRng =           [0.96;    0.96]   
            ηRng =             [0.0;     0.0] 
            κ₀Rng =            [2.0;     2.0]
            κ₁Rng =             [0.0;     0.0]
            κ₂Rng =             [0.0;     0.0]
            κ₃Rng =             [0.0;     0.0]
            λ₁Rng =            [0.0;     0.0]
            λ₂Rng =            [1.0;     1.0]
            μv₂Rng =               [0.0;     0.0]
            ν₁Rng =                [0.0;     0.0]    
            ν₂Rng =                [0.0;     0.0]
            ν₃Rng =                [0.0;     0.0]
            ν₄Rng =                [1.0;     1.0]
            ν₅Rng =                [1.0;     1.0]
            χuRng =               [0.0;     0.0]
            χdRng =               [0.0;     0.0]
            ξRng =                [0.0;     0.0]
            ζₐRng =              [0.0;     0.0]
            cd₀Rng =                [0.0;     0.0]  
            cm₀Rng =                [0.0;     0.0]
            cnαRng =      2π/β*[1.0;     1.0]  
            dtRng =         π/180*[0.0;     0.0]    
            dmRng =         π/180*[0.0;     0.0] 
            E₀Rng =                  [0.0;     0.0]  
            E₁Rng =                  [0.0;     0.0]
            f₀NRng =                [0.0;     0.0]
            f₀MRng =                [0.0;     0.0]
            f₀TRng =                [0.0;     0.0]
            fbNRng =                [0.70;     0.70]
            fbMRng =                [0.70;     0.70]
            fbTRng =                [0.70;     0.70]
            gᵥRng =                 [0.0;     0.0] 
            gᵥ₂Rng =                [0.0;     0.0]
            K₀Rng =                  [0.0;     0.0]  
            K₁Rng =                  [0.0;     0.0] 
            K₂Rng =                  [0.0;     0.0]   
            r₀Rng =             1e-2*[1.0;     1.0]  
            S1NRng =         π/180*[2.0;     2.0]
            S1MRng =         π/180*[2.0;     2.0]
            S1TRng =         π/180*[2.0;     2.0]
            S2NRng =         π/180*[2.0;     2.0]
            S2MRng =         π/180*[2.0;     2.0]
            S2TRng =         π/180*[2.0;     2.0]
            TaRng =                  [4.0;     4.0]    
            TfRng =                  [3.0;     3.0]    
            TvRng =                  [4.0;     4.0] 
            Tv₂Rng =                 [4.0;     4.0]
            VmRng =                  [0.4;     0.4]
            VtRng =                  [0.0;     0.0]
            Vn₁Rng =                 [0.5;     0.5]   
            Vn₂Rng =                 [0.0;     0.0]
            Vn₃Rng =                 [0.0;     0.0]
            ztdRng =               [1.0;     1.0]
            ztuRng =               [0.0;     0.0]
            zmRng =                [1.0;     1.0]
            if flapSiteID == 100
                cdδRng = 0.0*ones(length(MaRng))
                cmδRng = 0.0*ones(length(MaRng))
                cnδRng = 0.0*ones(length(MaRng))
            elseif flapSiteID == 75
                cdδRng =   0.0*ones(length(MaRng))
                cmδRng = -0.25*ones(length(MaRng))
                cnδRng =   1.0*ones(length(MaRng))    
            else
                error("Unavailable flap site ID")
            end
            # Fixed parameters
            λbWDiag = Diagonal([1.0; 1.0])
        else
            error("Airfoil not listed")
        end

        # Interpolated values
        α₀N = interpolate(MaRng,α₀NRng,Ma)
        αds₀ = interpolate(MaRng,αds₀Rng,Ma)
        αₛₛ = interpolate(MaRng,αₛₛRng,Ma)
        α1₀N = interpolate(MaRng,α1₀NRng,Ma)
        α1₀M = interpolate(MaRng,α1₀MRng,Ma)
        α1₀T = interpolate(MaRng,α1₀TRng,Ma)
        βσ1N = interpolate(MaRng,βσ1NRng,Ma)
        βσ1T = interpolate(MaRng,βσ1TRng,Ma)
        βσ2N = interpolate(MaRng,βσ2NRng,Ma)
        βS2Nlpr = interpolate(MaRng,βS2NlprRng,Ma)
        βS2Tlpr = interpolate(MaRng,βS2TlprRng,Ma)
        βS1Nu = interpolate(MaRng,βS1NuRng,Ma)
        βS1Mu = interpolate(MaRng,βS1MuRng,Ma)
        βS1Tu = interpolate(MaRng,βS1TuRng,Ma)
        βS1Nd = interpolate(MaRng,βS1NdRng,Ma)
        βS1Md = interpolate(MaRng,βS1MdRng,Ma)
        βS1Td = interpolate(MaRng,βS1TdRng,Ma)
        βS2Nu = interpolate(MaRng,βS2NuRng,Ma)
        βS2Mu = interpolate(MaRng,βS2MuRng,Ma)
        βS2Tu = interpolate(MaRng,βS2TuRng,Ma)
        βS2Nd = interpolate(MaRng,βS2NdRng,Ma)
        βS2Md = interpolate(MaRng,βS2MdRng,Ma)
        βS2Td = interpolate(MaRng,βS2TdRng,Ma)
        γLS = interpolate(MaRng,γLSRng,Ma)
        δα₀ = interpolate(MaRng,δα₀Rng,Ma)
        δα₁ = interpolate(MaRng,δα₁Rng,Ma)
        ϵₙ = interpolate(MaRng,ϵₙRng,Ma)
        ϵₘ = interpolate(MaRng,ϵₘRng,Ma)
        η = interpolate(MaRng,ηRng,Ma)
        κ₀ = interpolate(MaRng,κ₀Rng,Ma)
        κ₁ = interpolate(MaRng,κ₁Rng,Ma)
        κ₂ = interpolate(MaRng,κ₂Rng,Ma)
        κ₃ = interpolate(MaRng,κ₃Rng,Ma)
        λ₁ = interpolate(MaRng,λ₁Rng,Ma)
        λ₂ = interpolate(MaRng,λ₂Rng,Ma)
        μv₂ = interpolate(MaRng,μv₂Rng,Ma)
        ν₁ = interpolate(MaRng,ν₁Rng,Ma)
        ν₂ = interpolate(MaRng,ν₂Rng,Ma)
        ν₃ = interpolate(MaRng,ν₃Rng,Ma)
        ν₄ = interpolate(MaRng,ν₄Rng,Ma)
        ν₅ = interpolate(MaRng,ν₅Rng,Ma)
        χu = interpolate(MaRng,χuRng,Ma)
        χd = interpolate(MaRng,χdRng,Ma)
        ξ = interpolate(MaRng,ξRng,Ma)
        ζₐ = interpolate(MaRng,ζₐRng,Ma)
        cd₀ = interpolate(MaRng,cd₀Rng,Ma)
        cdδ = interpolate(MaRng,cdδRng,Ma)
        cm₀ = interpolate(MaRng,cm₀Rng,Ma)
        cmδ = interpolate(MaRng,cmδRng,Ma)
        cnα = interpolate(MaRng,cnαRng,Ma)
        cnδ = interpolate(MaRng,cnδRng,Ma)
        dt = interpolate(MaRng,dtRng,Ma)
        dm = interpolate(MaRng,dmRng,Ma)
        E₀ = interpolate(MaRng,E₀Rng,Ma)
        E₁ = interpolate(MaRng,E₁Rng,Ma)
        f₀N = interpolate(MaRng,f₀NRng,Ma)
        f₀M = interpolate(MaRng,f₀MRng,Ma)
        f₀T = interpolate(MaRng,f₀TRng,Ma)
        fbN = interpolate(MaRng,fbNRng,Ma)
        fbM = interpolate(MaRng,fbMRng,Ma)
        fbT = interpolate(MaRng,fbTRng,Ma)
        gᵥ = interpolate(MaRng,gᵥRng,Ma)
        gᵥ₂ = interpolate(MaRng,gᵥ₂Rng,Ma)
        K₀ = interpolate(MaRng,K₀Rng,Ma)
        K₁ = interpolate(MaRng,K₁Rng,Ma)
        K₂ = interpolate(MaRng,K₂Rng,Ma)
        r₀ = interpolate(MaRng,r₀Rng,Ma)
        S1N = interpolate(MaRng,S1NRng,Ma)
        S1M = interpolate(MaRng,S1MRng,Ma)
        S1T = interpolate(MaRng,S1TRng,Ma)
        S2N = interpolate(MaRng,S2NRng,Ma)
        S2M = interpolate(MaRng,S2MRng,Ma)
        S2T = interpolate(MaRng,S2TRng,Ma)
        Ta = interpolate(MaRng,TaRng,Ma)
        Tf = interpolate(MaRng,TfRng,Ma)
        Tv = interpolate(MaRng,TvRng,Ma)
        Tv₂ = interpolate(MaRng,Tv₂Rng,Ma)
        Vn₁ = interpolate(MaRng,Vn₁Rng,Ma)
        Vn₂ = interpolate(MaRng,Vn₂Rng,Ma)
        Vn₃ = interpolate(MaRng,Vn₃Rng,Ma)
        Vm = interpolate(MaRng,VmRng,Ma)
        Vt = interpolate(MaRng,VtRng,Ma)
        ztd = interpolate(MaRng,ztdRng,Ma)
        ztu = interpolate(MaRng,ztuRng,Ma)
        zm = interpolate(MaRng,zmRng,Ma)

        # Dimensionalize time delay constants
        if U > 0
            Ta *= b/U
            Tf *= b/U
            Tv *= b/U
            Tv₂ *= b/U
        end

        return new(α₀N,αds₀,αₛₛ,α1₀N,α1₀M,α1₀T,βσ1N,βσ1T,βσ2N,βS2Nlpr,βS2Tlpr,βS1Nu,βS1Mu,βS1Tu,βS1Nd,βS1Md,βS1Td,βS2Nu,βS2Mu,βS2Tu,βS2Nd,βS2Md,βS2Td,γLS,δα₀,δα₁,ϵₙ,ϵₘ,η,κ₀,κ₁,κ₂,κ₃,λ₁,λ₂,μv₂,ν₁,ν₂,ν₃,ν₄,ν₅,χu,χd,ξ,ζₐ,cd₀,cdδ,cm₀,cmδ,cnα,cnδ,dt,dm,E₀,E₁,f₀N,f₀M,f₀T,fbN,fbM,fbT,gᵥ,gᵥ₂,K₀,K₁,K₂,r₀,S1N,S1M,S1T,S2N,S2M,S2T,Ta,Tf,Tv,Tv₂,Vn₁,Vn₂,Vn₃,Vm,Vt,ztd,ztu,zm,λbWDiag)
    end

end


"""
@with_kw mutable struct Airfoil

    Airfoil composite type

# Fields
- name::String
- 
"""
@with_kw mutable struct Airfoil

    name::String
    attachedFlowParameters::AttachedFlowParameters
    separatedFlowParameters::SeparatedFlowParameters

end
export Airfoil


"""
create_Airfoil(;name::String,Re::Number=0,Ma::Number=0)

Initializes the airfoil with the predefined name 

# Arguments
- name::String
- Re::Number
- Ma::Number
"""
function create_Airfoil(;name::String,Re::Number=0,Ma::Number=0)

    attachedFlowParameters = AttachedFlowParameters(name,Re=Re,Ma=Ma)
    separatedFlowParameters = SeparatedFlowParameters(name,Re=Re,Ma=Ma)

    return Airfoil(name,attachedFlowParameters,separatedFlowParameters)
end
export create_Airfoil


"""
create_flapped_Airfoil(;name::String,flapSiteID::Int64,Re::Number=0,Ma::Number=0)

Initializes the airfoil with the predefined name and flap site ID 

# Arguments
- name::String
- flapSiteID::Int64
- Re::Number
- Ma::Number
"""
function create_flapped_Airfoil(;name::String,flapSiteID::Int64,Re::Number=0,Ma::Number=0)

    attachedFlowParameters = AttachedFlowParameters(name,Re=Re,Ma=Ma,flapSiteID=flapSiteID)
    separatedFlowParameters = SeparatedFlowParameters(name,Re=Re,Ma=Ma,flapSiteID=flapSiteID)
    
    return Airfoil(name,attachedFlowParameters,separatedFlowParameters)
end

# Sample airfoils
flatPlate = create_Airfoil(name="flatPlate")
NACA0002 = create_Airfoil(name="NACA0002")
NACA0006 = create_Airfoil(name="NACA0006")
NACA0012 = create_Airfoil(name="NACA0012")
NACA0018 = create_Airfoil(name="NACA0018")
NACA23012A = create_Airfoil(name="NACA23012A")
HeliosWingAirfoil = create_Airfoil(name="HeliosWingAirfoil")
HeliosPodAirfoil = create_Airfoil(name="HeliosPodAirfoil")
BWBAirfoil = create_Airfoil(name="BWBAirfoil")
export flatPlate, NACA0002, NACA0006, NACA0012, NACA0018, NACA23012A, HeliosWingAirfoil, HeliosPodAirfoil, BWBAirfoil