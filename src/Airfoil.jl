#
# @with_kw mutable struct AttachedFlowParameters

#     AttachedFlowParameters composite type

# Fields
# - `α₀N::Real`
# - `ϵₙ::Real`
# - `ϵₘ::Real`
# - `η::Real`
# - `cd₀::Real`
# - `cm₀::Real`
# - `cnα::Real`
# - `x_ac::Real`
#
@with_kw mutable struct AttachedFlowParameters

    α₀N::Real
    ϵₙ::Real
    ϵₘ::Real
    η::Real
    cd₀::Real
    cm₀::Real
    cnα::Real
    x_ac::Real

    function AttachedFlowParameters(name::String; Re::Real=0,Ma::Real=0)

        # Validate
        @assert all(>=(0), [Ma, Re])

        # Airfoil parameters' tables (as functions of Mach)
        if name in ["flatPlate","NACA0002","NACA0006"]
            # Bound Mach and corresponding compressibility factor
            Ma = min(0.8,Ma)
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [    0;   0.8]
            α₀NRng = π/180*[  0.0;   0.0]
            ϵₙRng =        [  0.7;   0.7]
            ϵₘRng =        [ 0.96;  0.96]
            ηRng =         [  1.0;   1.0]
            cd₀Rng =  1e-2*[  0.0;   0.0]
            cm₀Rng =  1e-3*[  0.0;   0.0]
            cnαRng =  2π/β*[  1.0;   1.0]
            x_acRng =      [ 0.25;  0.25]
        elseif name in ["NACA0012"]
            # Bound Mach and corresponding compressibility factor
            Ma = min(0.3,Ma)
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [0; 0.035; 0.072; 0.110; 0.185; 0.215; 0.25; 0.28;  0.3]
            α₀NRng = π/180*[0.0;  0.0;   0.0;   0.0;   0.0;   0.0;  0.0;  0.0;  0.0]
            ϵₙRng =        [0.7;  0.7;   0.7;   0.7;   0.7;   0.7;  0.7;  0.7;  0.7]
            ϵₘRng =        [0.96; 0.96;  0.96;  0.96;  0.96;  0.96; 0.96; 0.96; 0.96]
            ηRng =         [0.95; 0.95;  0.95;  0.95;  0.95;  0.95; 0.95; 0.95; 0.95]
            cd₀Rng =  1e-2*[1.2;  1.2;   1.2;   0.8;   0.5;   0.5;  0.5;  0.5;  0.5]
            cm₀Rng =  1e-3*[0;    0;   -14;    -5;    -5;    -5;   -5;   -5;   -5]
            cnαRng =  2π/β*[1.049; 1.049; 0.972; 0.998; 0.995; 0.980; 1.00; 1.00; 1.00]
            x_acRng =      [ 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25; 0.25]
        elseif name in ["NACA0015","NACA0018"]
            # Bound Mach and corresponding compressibility factor
            Ma = min(0.8,Ma)
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [    0;  0.8]
            α₀NRng = π/180*[  0.0;  0.0]
            ϵₙRng =        [  0.7;  0.7]
            ϵₘRng =        [ 0.96; 0.96]
            ηRng =         [ 0.95; 0.95]
            cd₀Rng =  1e-2*[  1.0;  0.6] 
            cm₀Rng =  1e-3*[  1.0;  1.0]
            cnαRng =  2π/β*[  1.0;  1.0]
            x_acRng =      [ 0.25; 0.25]
        elseif name in ["VERTOL23010"]
            # Bound Mach and corresponding compressibility factor
            Ma = min(0.8,Ma)
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  =        [    0;  0.8]
            α₀NRng = -π/180*[  0.7;  0.7]
            ϵₙRng =         [  0.7;  0.7]
            ϵₘRng =         [ 0.96; 0.96]
            ηRng =          [ 0.95; 0.95]
            cd₀Rng =   1e-2*[  1.0;  1.0] 
            cm₀Rng =   1e-2*[  0.0;  0.0]
            cnαRng =   2π/β*[  1.0;  1.0]
            x_acRng =       [ 0.25; 0.25]  
        elseif name in ["NACA23012A"]
            # Bound Mach and corresponding compressibility factor
            Ma = min(0.8,Ma)
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  =       [     0;   0.8]
            α₀NRng = -π/180*[   1.2;   1.2]
            ϵₙRng =        [   0.7;   0.7]
            ϵₘRng =        [  0.96;  0.96]
            ηRng =         [     1;     1]
            cd₀Rng =  1e-2*[   0.3;   0.3] 
            cm₀Rng =  1e-2*[   5.0;   5.0]
            cnαRng =  2π*[  1.08;  1.08]
            x_acRng =      [ 0.259; 0.259]
        elseif name in ["HeliosWingAirfoil"]
            # Bound Mach
            Ma = min(0.8,Ma)
            # Mach-dependent parameters
            MaRng  =       [   0; 0.8]
            α₀NRng = π/180*[ 0.0; 0.0]
            ϵₙRng =        [ 1.0; 1.0]
            ϵₘRng =        [ 1.0; 1.0]
            ηRng =         [ 1.0; 1.0]
            cd₀Rng =  1e-2*[ 1.0; 1.0]  
            cm₀Rng =  1e-2*[ 2.5; 2.5]
            cnαRng =    2π*[ 1.0; 1.0]
            x_acRng =      [ 0.25;  0.25]   
        elseif name in ["HeliosPodAirfoil"]
            # Bound Mach
            Ma = min(0.8,Ma)
            # Mach-dependent parameters
            MaRng  =       [   0;   0.8]
            α₀NRng = π/180*[ 0.0;   0.0]
            ϵₙRng =        [ 1.0;   1.0]
            ϵₘRng =        [ 1.0;   1.0]
            ηRng =         [ 1.0;   1.0]
            cd₀Rng =  1e-2*[ 2.0;   2.0]  
            cm₀Rng =  1e-3*[ 0.0;   0.0]
            cnαRng =       [ 5.0;   5.0]
            x_acRng =      [ 0.25; 0.25]
        elseif name in ["BWBAirfoil"] 
            # Bound Mach and corresponding compressibility factor
            Ma = min(0.8,Ma)
            β = sqrt(1-Ma^2)    
            # Mach-dependent parameters
            MaRng  =       [   0;   0.8]
            α₀NRng = π/180*[ 0.0;   0.0]
            ϵₙRng =        [ 1.0;   1.0]
            ϵₘRng =        [ 1.0;   1.0]
            ηRng =         [ 1.0;   1.0]
            cd₀Rng =  1e-2*[ 1.0;   1.0]
            cm₀Rng =  1e-1*[ 1.0;   1.0]
            cnαRng =  2π/β*[ 1.0;   1.0]
            x_acRng =      [ 0.25; 0.25]
        elseif name in ["cHALEairfoil"]
            # Bound Mach and corresponding compressibility factor
            Ma = min(0.8,Ma)
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [    0;   0.8]
            α₀NRng = π/180*[  0.0;   0.0]
            ϵₙRng =        [  0.7;   0.7]
            ϵₘRng =        [ 0.96;  0.96]
            ηRng =         [ 0.95;  0.95]
            cd₀Rng =  1e-2*[  0.0;   0.0]
            cm₀Rng =  1e-3*[  0.0;   0.0]
            cnαRng =  2π/β*[  1.0;   1.0]
            x_acRng =      [ 0.25;  0.25]
        else
            error("Airfoil not listed")
        end

        # Interpolated values
        α₀N = LinearInterpolations.interpolate(MaRng,α₀NRng,Ma)
        ϵₙ  = LinearInterpolations.interpolate(MaRng,ϵₙRng,Ma)
        ϵₘ  = LinearInterpolations.interpolate(MaRng,ϵₘRng,Ma)
        η   = LinearInterpolations.interpolate(MaRng,ηRng,Ma)
        cd₀ = LinearInterpolations.interpolate(MaRng,cd₀Rng,Ma)
        cm₀ = LinearInterpolations.interpolate(MaRng,cm₀Rng,Ma)
        cnα = LinearInterpolations.interpolate(MaRng,cnαRng,Ma)
        x_ac = LinearInterpolations.interpolate(MaRng,x_acRng,Ma)

        return new(α₀N,ϵₙ,ϵₘ,η,cd₀,cm₀,cnα,x_ac)
    end

end


#
# @with_kw mutable struct BLiParameters

#     BLiParameters composite type

# Fields
# - `α₀N::Real`
# - `αds₀::Real`
# - `αₛₛ::Real`
# - `α1₀N::Real`
# - `α1₀M::Real`
# - `α1₀T::Real`
# - `βσ1N::Real`
# - `βσ1T::Real`
# - `βσ2N::Real`
# - `βS2Nlpr::Real`
# - `βS2Tlpr::Real`
# - `βS1Nu::Real`
# - `βS1Mu::Real`
# - `βS1Tu::Real`
# - `βS1Nd::Real`
# - `βS1Md::Real`
# - `βS1Td::Real`
# - `βS2Nu::Real`
# - `βS2Mu::Real`
# - `βS2Tu::Real`
# - `βS2Nd::Real`
# - `βS2Md::Real`
# - `βS2Td::Real`
# - `γLS::Real`
# - `δα₀::Real`
# - `δα₁::Real`
# - `ϵₙ::Real`
# - `ϵₘ::Real`
# - `η::Real`
# - `κ₀::Real`
# - `κ₁::Real`
# - `κ₂::Real`
# - `κ₃::Real`
# - `λ₁::Real`
# - `λ₂::Real`
# - `μv₂::Real`
# - `ν₁::Real`
# - `ν₂::Real`
# - `ν₃::Real`
# - `ν₄::Real`
# - `ν₅::Real`
# - `χu::Real`
# - `χd::Real`
# - `ξ::Real`
# - `ζₐ::Real`
# - `cd₀::Real`
# - `cm₀::Real`
# - `cnα::Real`
# - `dt::Real`
# - `dm::Real`
# - `E₀::Real`
# - `E₁::Real`
# - `f₀N::Real`
# - `f₀M::Real`
# - `f₀T::Real`
# - `fbN::Real`
# - `fbM::Real`
# - `fbT::Real`
# - `gᵥ::Real`
# - `gᵥ₂::Real`
# - `K₀::Real`
# - `K₁::Real`
# - `K₂::Real`
# - `r₀::Real`
# - `S1N::Real`
# - `S1M::Real`
# - `S1T::Real`
# - `S2N::Real`
# - `S2M::Real`
# - `S2T::Real`
# - `Ta::Real`
# - `Tf::Real`
# - `Tg::Real`
# - `Tv::Real`
# - `Tv₂::Real`
# - `Vn₁::Real`
# - `Vn₂::Real`
# - `Vn₃::Real`
# - `Vm::Real`
# - `Vt::Real`
# - `ztd::Real`
# - `ztu::Real`
# - `zm::Real`
# - `γbC::Vector{Float64}`
# - `γbCMat::Matrix{Float64}`
#
@with_kw mutable struct BLiParameters

    α₀N::Real
    αds₀::Real
    αₛₛ::Real
    α1₀N::Real
    α1₀M::Real
    α1₀T::Real
    βσ1N::Real
    βσ1T::Real
    βσ2N::Real
    βS2Nlpr::Real
    βS2Tlpr::Real
    βS1Nu::Real
    βS1Mu::Real
    βS1Tu::Real
    βS1Nd::Real
    βS1Md::Real
    βS1Td::Real
    βS2Nu::Real
    βS2Mu::Real
    βS2Tu::Real
    βS2Nd::Real
    βS2Md::Real
    βS2Td::Real
    γLS::Real
    δα₀::Real
    δα₁::Real
    ϵₙ::Real
    ϵₘ::Real
    η::Real
    κ₀::Real
    κ₁::Real
    κ₂::Real
    κ₃::Real
    λ₁::Real
    λ₂::Real
    μv₂::Real
    ν₁::Real
    ν₂::Real
    ν₃::Real
    ν₄::Real
    ν₅::Real
    χu::Real
    χd::Real
    ξ::Real
    ζₐ::Real
    cd₀::Real
    cm₀::Real
    cnα::Real
    dt::Real
    dm::Real
    E₀::Real
    E₁::Real
    f₀N::Real
    f₀M::Real
    f₀T::Real
    fbN::Real
    fbM::Real
    fbT::Real
    gᵥ::Real
    gᵥ₂::Real
    K₀::Real
    K₁::Real
    K₂::Real
    r₀::Real
    S1N::Real
    S1M::Real
    S1T::Real
    S2N::Real
    S2M::Real
    S2T::Real
    Ta::Real
    Tf::Real
    Tg::Real
    Tv::Real
    Tv₂::Real
    Vn₁::Real
    Vn₂::Real
    Vn₃::Real
    Vm::Real
    Vt::Real
    ztd::Real
    ztu::Real
    zm::Real
    γbC::Vector{Float64}
    γbCMat::Matrix{Float64}
    
    function BLiParameters(name::String; Re::Real=0,Ma::Real=0,U::Real=0,b::Real=0)

        # Validate
        @assert all(>=(0), [Ma, Re, U, b])

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
            TgRng =         5*TaRng
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
            # Fixed parameters
            γbC = [2.5; 0.8]
            γbCMat = Diagonal(γbC)
        elseif name in ["NACA0012-GU","NACA0015","NACA0015-s"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.078,min(0.155,Ma))
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [0.078;   0.117;   0.155]
            α₀NRng =     -π/180*[0.0;     0.0;     0.0]    
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
            TgRng =             5*TaRng    
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
            # Fixed parameters
            γbC = [2.5; 0.8]
            γbCMat = Diagonal(γbC)
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
            TgRng =         5*TaRng    
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
            # Fixed parameters
            γbC = [1.0; 1.0]
            γbCMat = Diagonal(γbC)
        elseif name in ["NACA23012A"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.001,min(0.3,Ma))
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  =       [  0.001;     0.3]
            α₀NRng =     -π/180*[1.2;      1.2]
            αds₀Rng =    π/180*[17.7;     17.7]   
            αₛₛRng =     π/180*[14.4;     14.4]   
            α1₀NRng =    π/180*[14.4;     14.4]
            α1₀MRng =    π/180*[13.7;     13.7]
            α1₀TRng =    π/180*[15.5;     15.5]
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
            cnαRng =      2π*[1.08;     1.08]  
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
            K₀Rng =                  [-0.009;   -0.009]  
            K₁Rng =                  [-0.171;    -0.171] 
            K₂Rng =                  [0.015;    0.015]   
            r₀Rng =             1e-2*[1.62;     1.62]  
            S1NRng =         π/180*[1.5;     1.5]
            S1MRng =         π/180*[1.89;     1.89]
            S1TRng =         π/180*[1.27;     1.27]
            S2NRng =         π/180*[2.05;     2.05]
            S2MRng =         π/180*[4.70;     4.70]
            S2TRng =         π/180*[4.74;     4.74]
            TaRng =                  [3.74;     3.74]    
            TfRng =                  [3.75;     3.75]
            TgRng =         5*TaRng    
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
            # Fixed parameters
            γbC = [1.0; 1.0]
            γbCMat = Diagonal(γbC)
        elseif name in ["HeliosPodAirfoil"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.001,min(0.3,Ma))
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  =       [  0.001;     0.3]
            α₀NRng =     π/180*[0.0;     0.0]
            αds₀Rng =   π/180*[30.0;    30.0]   
            αₛₛRng =     π/180*[25.0;    25.0]   
            α1₀NRng =   π/180*[25.0;    25.0]
            α1₀MRng =   π/180*[25.0;    25.0]
            α1₀TRng =   π/180*[25.0;    25.0]
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
            cd₀Rng =                [0.02;     0.02]  
            cm₀Rng =                [0.0;     0.0]
            cnαRng =      [5.0;     5.0]  
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
            TgRng =         5*TaRng    
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
            # Fixed parameters
            γbC = [1.0; 1.0]
            γbCMat = Diagonal(γbC)    
        elseif name in ["VERTOL23010"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.2,min(0.6,Ma))
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng  =       [0.2;   0.4;   0.6]
            α₀NRng =     -π/180*[0.7;     0.7;      0.7]    
            αds₀Rng =    π/180*[20.0;    17.0;     20.0]   
            αₛₛRng =     π/180*[14.5;    13.5;     17.0]   
            α1₀NRng =    π/180*[14.5;    14.0;     8.0]
            α1₀MRng =    π/180*[14.5;    14.0;     8.0]
            α1₀TRng =    π/180*[14.5;    14.0;     8.0]
            βσ1NRng =          [0.00;    0.00;    0.00]
            βσ1TRng =          [0.00;    0.00;    0.00]        
            βσ2NRng =          [0.00;    0.00;    0.00]         
            βS2NlprRng =        [0.00;    0.00;    0.00]   
            βS2TlprRng =        [0.00;    0.00;    0.00]    
            βS1NuRng =          [0.00;    0.00;    0.00]       
            βS1MuRng =          [0.00;    0.00;    0.00]        
            βS1TuRng =          [0.00;    0.00;    0.00]       
            βS1NdRng =          [0.00;    0.00;    0.00]           
            βS1MdRng =          [0.00;    0.00;    0.00]      
            βS1TdRng =          [0.00;    0.00;    0.00]
            βS2NuRng =          [0.00;    0.00;    0.00]
            βS2MuRng =          [0.00;    0.00;    0.00]
            βS2TuRng =          [0.00;    0.00;    0.00]                 
            βS2NdRng =          [0.00;    0.00;    0.00] 
            βS2MdRng =          [0.00;    0.00;    0.00]         
            βS2TdRng =          [0.00;    0.00;    0.00]       
            γLSRng =            [14.5;    14.0;     8.0]
            δα₀Rng= π/180*[1.5;     0.9;      2.0]    
            δα₁Rng= π/180*[1.5;     5.0;      2.0]      
            ϵₙRng =              [0.70;    0.70;    0.70]  
            ϵₘRng =              [0.96;    0.96;    0.96]   
            ηRng =                 [0.10;    0.35;     0.10] 
            κ₀Rng =             [3.37;    3.37;    3.37]
            κ₁Rng =             [0.45;    0.45;    0.45]
            κ₂Rng =             [0.20;    0.14;     0.20]
            κ₃Rng =             [2.00;    2.00;    2.00]
            λ₁Rng =            [50;       50;       50]
            λ₂Rng =            [0.7;     0.7;     0.7]
            μv₂Rng =               [1.00;    1.00;    1.00]
            ν₁Rng =                [0.8;     0.85;     1.5]    
            ν₂Rng =                [1.0;     1.45;     1.0]
            ν₃Rng =                [0.25;    0.25;     0.25]
            ν₄Rng =                [0.50;    0.50;    0.50]
            ν₅Rng =                [1.00;    1.00;    1.00]
            χuRng =               [0.00;   0.00;    0.00]
            χdRng =               [0.00;   0.00;    0.00]
            ξRng =                  [0.00;    0.00;    0.00]
            ζₐRng =              [0.00;    0.00;    0.00]
            cd₀Rng =           1e-2*[1.0;     1.0;     1.0]  
            cm₀Rng =           1e-2*[0.0;    -0.0;    -0.0]
            cnαRng =        2π/β*[1.0;   1.0;   1.0]  
            dtRng =         π/180*[0.0;     0.0;      0.0]    
            dmRng =         π/180*[0.0;     0.0;      0.0] 
            E₀Rng =                  [1.0;     1.0;     1.0]  
            E₁Rng =                  [0.30;    0.30;     0.30]
            f₀NRng =                [0.000;   0.000;   0.000]
            f₀MRng =                [0.000;   0.000;   0.000]
            f₀TRng =                [0.000;   0.000;   0.000]
            fbNRng =                [0.80;    0.80;    0.80]
            fbMRng =                [0.80;    0.80;    0.80]
            fbTRng =                [0.80;    0.80;    0.80]
            gᵥRng =                 [1.5;     0.5;      1.5] 
            gᵥ₂Rng =                [1.00;    1.00;    1.00]
            K₀Rng =             1e-3*[0.0;     0.0;     0.0]
            K₁Rng =                  [-0.100; -0.123;   -0.080]
            K₂Rng =                  [0.020;   0.032;   -0.020]   
            r₀Rng =             1e-2*[1.00;    2.64;    1.60]  
            S1NRng =         π/180*[1.5;     1.67;      3.0]
            S1MRng =         π/180*[1.5;     1.67;      3.0]
            S1TRng =         π/180*[1.5;     1.67;      3.0]
            S2NRng =         π/180*[2.0;     1.80;      4.5]
            S2MRng =         π/180*[2.0;     1.30;      4.5]
            S2TRng =         π/180*[2.0;     1.30;      4.5]
            TaRng =                  [2.0;     2.5;      1.5]    
            TfRng =                  [3.0;     3.5;      3.0]
            TgRng =             5*TaRng    
            TvRng =                  [5.0;     6.0;      5.0] 
            Tv₂Rng =                 [5.0;     6.0;      5.0]
            VmRng =                  [0.35;    0.35;     0.35]
            VtRng =                  [0.16;    0.24;    0.23]
            Vn₁Rng =                 [1.5;     1.00;     0.50]   
            Vn₂Rng =                 [1.0;     0.45;     0.10]
            Vn₃Rng =                 [0.00;    0.00;    0.00] 
            ztdRng =               [1.25;    1.5;      1.2]
            ztuRng =               [0.00;    0.00;    0.00]
            zmRng =                [0.5;     3.0;      0.5]
            # Fixed parameters
            γbC = [2.5; 0.8]
            γbCMat = Diagonal(γbC)
        elseif name in ["flatPlate","NACA0002","NACA0006","HeliosWingAirfoil","BWBAirfoil","cHALEairfoil"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.001,min(0.8,Ma))
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  =       [  0.001;     0.8]
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
            TgRng =         5*TaRng    
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
            # Fixed parameters
            γbC = [1.0; 1.0]
            γbCMat = Diagonal(γbC)
        else
            error("Airfoil not listed")
        end

        # Set effective angles
        αds₀Rng += α₀NRng   
        αₛₛRng += α₀NRng   
        α1₀NRng += α₀NRng
        α1₀MRng += α₀NRng
        α1₀TRng += α₀NRng

        # Interpolated values
        α₀N = LinearInterpolations.interpolate(MaRng,α₀NRng,Ma)
        αds₀ = LinearInterpolations.interpolate(MaRng,αds₀Rng,Ma)
        αₛₛ = LinearInterpolations.interpolate(MaRng,αₛₛRng,Ma)
        α1₀N = LinearInterpolations.interpolate(MaRng,α1₀NRng,Ma)
        α1₀M = LinearInterpolations.interpolate(MaRng,α1₀MRng,Ma)
        α1₀T = LinearInterpolations.interpolate(MaRng,α1₀TRng,Ma)
        βσ1N = LinearInterpolations.interpolate(MaRng,βσ1NRng,Ma)
        βσ1T = LinearInterpolations.interpolate(MaRng,βσ1TRng,Ma)
        βσ2N = LinearInterpolations.interpolate(MaRng,βσ2NRng,Ma)
        βS2Nlpr = LinearInterpolations.interpolate(MaRng,βS2NlprRng,Ma)
        βS2Tlpr = LinearInterpolations.interpolate(MaRng,βS2TlprRng,Ma)
        βS1Nu = LinearInterpolations.interpolate(MaRng,βS1NuRng,Ma)
        βS1Mu = LinearInterpolations.interpolate(MaRng,βS1MuRng,Ma)
        βS1Tu = LinearInterpolations.interpolate(MaRng,βS1TuRng,Ma)
        βS1Nd = LinearInterpolations.interpolate(MaRng,βS1NdRng,Ma)
        βS1Md = LinearInterpolations.interpolate(MaRng,βS1MdRng,Ma)
        βS1Td = LinearInterpolations.interpolate(MaRng,βS1TdRng,Ma)
        βS2Nu = LinearInterpolations.interpolate(MaRng,βS2NuRng,Ma)
        βS2Mu = LinearInterpolations.interpolate(MaRng,βS2MuRng,Ma)
        βS2Tu = LinearInterpolations.interpolate(MaRng,βS2TuRng,Ma)
        βS2Nd = LinearInterpolations.interpolate(MaRng,βS2NdRng,Ma)
        βS2Md = LinearInterpolations.interpolate(MaRng,βS2MdRng,Ma)
        βS2Td = LinearInterpolations.interpolate(MaRng,βS2TdRng,Ma)
        γLS = LinearInterpolations.interpolate(MaRng,γLSRng,Ma)
        δα₀ = LinearInterpolations.interpolate(MaRng,δα₀Rng,Ma)
        δα₁ = LinearInterpolations.interpolate(MaRng,δα₁Rng,Ma)
        ϵₙ = LinearInterpolations.interpolate(MaRng,ϵₙRng,Ma)
        ϵₘ = LinearInterpolations.interpolate(MaRng,ϵₘRng,Ma)
        η = LinearInterpolations.interpolate(MaRng,ηRng,Ma)
        κ₀ = LinearInterpolations.interpolate(MaRng,κ₀Rng,Ma)
        κ₁ = LinearInterpolations.interpolate(MaRng,κ₁Rng,Ma)
        κ₂ = LinearInterpolations.interpolate(MaRng,κ₂Rng,Ma)
        κ₃ = LinearInterpolations.interpolate(MaRng,κ₃Rng,Ma)
        λ₁ = LinearInterpolations.interpolate(MaRng,λ₁Rng,Ma)
        λ₂ = LinearInterpolations.interpolate(MaRng,λ₂Rng,Ma)
        μv₂ = LinearInterpolations.interpolate(MaRng,μv₂Rng,Ma)
        ν₁ = LinearInterpolations.interpolate(MaRng,ν₁Rng,Ma)
        ν₂ = LinearInterpolations.interpolate(MaRng,ν₂Rng,Ma)
        ν₃ = LinearInterpolations.interpolate(MaRng,ν₃Rng,Ma)
        ν₄ = LinearInterpolations.interpolate(MaRng,ν₄Rng,Ma)
        ν₅ = LinearInterpolations.interpolate(MaRng,ν₅Rng,Ma)
        χu = LinearInterpolations.interpolate(MaRng,χuRng,Ma)
        χd = LinearInterpolations.interpolate(MaRng,χdRng,Ma)
        ξ = LinearInterpolations.interpolate(MaRng,ξRng,Ma)
        ζₐ = LinearInterpolations.interpolate(MaRng,ζₐRng,Ma)
        cd₀ = LinearInterpolations.interpolate(MaRng,cd₀Rng,Ma)
        cm₀ = LinearInterpolations.interpolate(MaRng,cm₀Rng,Ma)
        cnα = LinearInterpolations.interpolate(MaRng,cnαRng,Ma)
        dt = LinearInterpolations.interpolate(MaRng,dtRng,Ma)
        dm = LinearInterpolations.interpolate(MaRng,dmRng,Ma)
        E₀ = LinearInterpolations.interpolate(MaRng,E₀Rng,Ma)
        E₁ = LinearInterpolations.interpolate(MaRng,E₁Rng,Ma)
        f₀N = LinearInterpolations.interpolate(MaRng,f₀NRng,Ma)
        f₀M = LinearInterpolations.interpolate(MaRng,f₀MRng,Ma)
        f₀T = LinearInterpolations.interpolate(MaRng,f₀TRng,Ma)
        fbN = LinearInterpolations.interpolate(MaRng,fbNRng,Ma)
        fbM = LinearInterpolations.interpolate(MaRng,fbMRng,Ma)
        fbT = LinearInterpolations.interpolate(MaRng,fbTRng,Ma)
        gᵥ = LinearInterpolations.interpolate(MaRng,gᵥRng,Ma)
        gᵥ₂ = LinearInterpolations.interpolate(MaRng,gᵥ₂Rng,Ma)
        K₀ = LinearInterpolations.interpolate(MaRng,K₀Rng,Ma)
        K₁ = LinearInterpolations.interpolate(MaRng,K₁Rng,Ma)
        K₂ = LinearInterpolations.interpolate(MaRng,K₂Rng,Ma)
        r₀ = LinearInterpolations.interpolate(MaRng,r₀Rng,Ma)
        S1N = LinearInterpolations.interpolate(MaRng,S1NRng,Ma)
        S1M = LinearInterpolations.interpolate(MaRng,S1MRng,Ma)
        S1T = LinearInterpolations.interpolate(MaRng,S1TRng,Ma)
        S2N = LinearInterpolations.interpolate(MaRng,S2NRng,Ma)
        S2M = LinearInterpolations.interpolate(MaRng,S2MRng,Ma)
        S2T = LinearInterpolations.interpolate(MaRng,S2TRng,Ma)
        Ta = LinearInterpolations.interpolate(MaRng,TaRng,Ma)
        Tf = LinearInterpolations.interpolate(MaRng,TfRng,Ma)
        Tg = LinearInterpolations.interpolate(MaRng,TgRng,Ma)
        Tv = LinearInterpolations.interpolate(MaRng,TvRng,Ma)
        Tv₂ = LinearInterpolations.interpolate(MaRng,Tv₂Rng,Ma)
        Vn₁ = LinearInterpolations.interpolate(MaRng,Vn₁Rng,Ma)
        Vn₂ = LinearInterpolations.interpolate(MaRng,Vn₂Rng,Ma)
        Vn₃ = LinearInterpolations.interpolate(MaRng,Vn₃Rng,Ma)
        Vm = LinearInterpolations.interpolate(MaRng,VmRng,Ma)
        Vt = LinearInterpolations.interpolate(MaRng,VtRng,Ma)
        ztd = LinearInterpolations.interpolate(MaRng,ztdRng,Ma)
        ztu = LinearInterpolations.interpolate(MaRng,ztuRng,Ma)
        zm = LinearInterpolations.interpolate(MaRng,zmRng,Ma)

        # Dimensionalize time delay constants
        if U > 0 && b > 0
            Ta *= b/U
            Tf *= b/U
            Tg *= b/U
            Tv *= b/U
            Tv₂ *= b/U
        end

        return new(α₀N,αds₀,αₛₛ,α1₀N,α1₀M,α1₀T,βσ1N,βσ1T,βσ2N,βS2Nlpr,βS2Tlpr,βS1Nu,βS1Mu,βS1Tu,βS1Nd,βS1Md,βS1Td,βS2Nu,βS2Mu,βS2Tu,βS2Nd,βS2Md,βS2Td,γLS,δα₀,δα₁,ϵₙ,ϵₘ,η,κ₀,κ₁,κ₂,κ₃,λ₁,λ₂,μv₂,ν₁,ν₂,ν₃,ν₄,ν₅,χu,χd,ξ,ζₐ,cd₀,cm₀,cnα,dt,dm,E₀,E₁,f₀N,f₀M,f₀T,fbN,fbM,fbT,gᵥ,gᵥ₂,K₀,K₁,K₂,r₀,S1N,S1M,S1T,S2N,S2M,S2T,Ta,Tf,Tg,Tv,Tv₂,Vn₁,Vn₂,Vn₃,Vm,Vt,ztd,ztu,zm,γbC,γbCMat)
    end

end


# @with_kw mutable struct BLoParameters

#     BLoParameters composite type

# Fields
# - `α₀N::Real`
# - `α1₀::Real`
# - `δα::Real`
# - `ϵₙ::Real`
# - `ϵₘ::Real`
# - `η::Real`
# - `cd₀::Real`
# - `cm₀::Real`
# - `cn₁::Real`
# - `cnα::Real`
# - `Df::Real`
# - `E₀::Real`
# - `f₀::Real`
# - `fb::Real`
# - `K₀::Real`
# - `K₁::Real`
# - `K₂::Real`
# - `S1::Real`
# - `S2::Real`
# - `Tf₀::Real`
# - `Tp::Real`
# - `Tv₀::Real`
# - `TvL::Real`
# - `γbC::Vector{Float64}`
# - `γbCMat::Matrix{Float64}`
#
@with_kw mutable struct BLoParameters

    α₀N::Real
    α1₀::Real
    δα::Real
    ϵₙ::Real
    ϵₘ::Real
    η::Real
    cd₀::Real
    cm₀::Real
    cn₁::Real
    cnα::Real
    Df::Real
    E₀::Real
    f₀::Real
    fb::Real
    K₀::Real
    K₁::Real
    K₂::Real
    S1::Real
    S2::Real
    Tf₀::Real
    Tp::Real
    Tv₀::Real
    TvL::Real
    γbC::Vector{Float64}
    γbCMat::Matrix{Float64}
    
    function BLoParameters(name::String; Re::Real=0,Ma::Real=0,U::Real=0,b::Real=0)

        # Validate
        @assert all(>=(0), [Ma, Re, U, b])

        # Airfoil parameters' tables (as functions of Mach)
        if name in ["NACA0012"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.035,min(0.8,Ma))
            β = sqrt(1-Ma^2)
            # Mach-dependent parameters
            MaRng = [0.035; 0.072; 0.110; 0.185; 0.215; 0.25; 0.3; 0.4; 0.5; 0.6; 0.7; 0.75; 0.8]
            α₀NRng = π/180*[-0.0; -0.0; -0.0; -0.0; -0.0; -0.0; -0.0; -0.0; -0.0; -0.0; -0.0; -0.0; -0.0]
            α1₀Rng = π/180*[12.4; 13.8; 15.2; 16.3; 16.5; 16.0; 13.7; 12.5; 10.5; 8.5; 5.6; 3.5; 0.7]
            δαRng = π/180*[2.0; 1.0; 2.5; 2.5; 2.5; 2.5; 0.5; 2.0; 1.45; 1.0; 0.8; 0.2; 0.1]
            ϵₙRng = [0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70] 
            ϵₘRng = [0.96; 0.96; 0.96; 0.96; 0.96; 0.96; 0.96; 0.96; 0.96; 0.96; 0.96; 0.96; 0.96]
            cd₀Rng = [0.010; 0.010; 0.010; 0.010; 0.010; 0.010; 0.010; 0.008; 0.0077; 0.0078; 0.0078; 0.0079; 0.0114]
            cm₀Rng = [-0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037]
            ηRng = [0.95; 0.95; 0.95; 0.95; 0.95; 0.95; 0.95; 0.95; 0.95; 0.95; 0.95; 0.95; 0.95]
            cn₁Rng = [1.25; 1.45; 1.60; 1.80; 1.80; 1.85; 1.60; 1.2; 1.05; 0.92; 0.68; 0.5; 0.18]
            cnαRng = 180/π*[0.105; 0.108; 0.108; 0.110; 0.113; 0.115; 0.116; 0.113; 0.117; 0.127; 0.154; 0.175; 0.216]
            DfRng = [8.0; 8.0; 8.0; 8.0; 8.0; 8.0; 8.0; 7.75; 6.2; 6.0; 5.9; 5.5; 4.0]
            E₀Rng = [0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00]
            f₀Rng = [0.02; 0.02; 0.02; 0.02; 0.02; 0.02; 0.02; 0.02; 0.02; 0.02; 0.02; 0.02; 0.02]
            fbRng = [0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70]
            K₀Rng = [0.0025; 0.0025; 0.0025; 0.0025; 0.0025; 0.0025; 0.0025; 0.006; 0.02; 0.038; 0.03; 0.001; -0.01]
            K₁Rng = [-0.120; -0.120; -0.120; -0.120; -0.120; -0.120; -0.120; -0.135; -0.125; -0.12; -0.09; -0.13; 0.02]
            K₂Rng = [0.04; 0.04; 0.04; 0.04; 0.04; 0.04; 0.04; 0.05; 0.04; 0.04; 0.15; -0.02; -0.01]
            S1Rng = π/180*[3.0; 3.0; 3.0; 3.5; 3.5; 3.5; 3.0; 3.25; 3.5; 4.0; 4.5; 3.5; 0.70]
            S2Rng = π/180*[1.5; 1.5; 1.5; 2.0; 2.0; 2.0; 1.5; 1.6; 1.2; 0.7; 0.5; 0.8; 0.18]
            Tf₀Rng = [3.0; 3.0; 3.0; 3.0; 3.0; 3.0; 3.0; 2.5; 2.2; 2.0; 2.0; 2.0; 2.0]
            TpRng = [1.7; 1.7; 1.7; 1.7; 1.7; 1.7; 1.7; 1.8; 2.0; 2.5; 3.0; 3.3; 4.3]
            Tv₀Rng = [6.0; 6.0; 6.0; 6.0; 6.0; 6.0; 6.0; 6.0; 6.0; 6.0; 6.0; 6.0; 4.0]
            TvLRng = [4.0; 5.0; 5.0; 5.0; 5.0; 5.0; 5.0; 9.0; 9.0; 9.0; 9.0; 9.0; 9.0]
            # Fixed parameters
            γbC = [2.5; 0.8]
            γbCMat = Diagonal(γbC)
        elseif name in ["NACA23012A"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.001,min(0.3,Ma))
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  = [0.001; 0.3]
            α₀NRng = -π/180*[1.2; 1.2]
            α1₀Rng = π/180*[13.5; 13.5]
            δαRng = π/180*[1.0; 1.0]
            ϵₙRng = [0.70; 0.70]
            ϵₘRng = [0.96; 0.96]
            cd₀Rng = [0.003; 0.003]
            cm₀Rng = [0.05; 0.05]
            ηRng = [1.0; 1.0]
            cn₁Rng = [1.25; 1.85]
            cnαRng = 2π*[1.08; 1.08]
            DfRng = [8.0; 8.0]
            E₀Rng = [0.05; 0.05]
            f₀Rng = [0.02; 0.02]
            fbRng = [0.70; 0.70]
            K₀Rng = [-0.009; -0.009]
            K₁Rng = [-0.17; -0.17]
            K₂Rng = [0.015; 0.015]
            S1Rng = π/180*[2.5; 2.5]
            S2Rng = π/180*[2.0; 2.0]
            Tf₀Rng = [3.75; 3.75]
            TpRng = [1.7; 1.7]
            Tv₀Rng = [5.0; 5.0]
            TvLRng = [5.0; 5.0]
            # Fixed parameters
            γbC = [1.0; 1.0]
            γbCMat = Diagonal(γbC)
        elseif name in ["HeliosPodAirfoil"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.001,min(0.8,Ma))
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  = [0.001; 0.8]
            α₀NRng = π/180*[0.0; 0.0]
            α1₀Rng = π/180*[25.0; 25.0]
            δαRng = π/180*[0.0; 0.0]
            ϵₙRng = [0.70; 0.70] 
            ϵₘRng = [0.96; 0.96]
            cd₀Rng = [0.02; 0.02]
            cm₀Rng = [0.0; 0.0]
            ηRng = [1.0; 1.0]
            cn₁Rng = [2.0; 2.0]
            cnαRng = [5.0; 5.0]
            DfRng = [8.0; 8.0]
            E₀Rng = [0.0; 0.0]
            f₀Rng = [0.02; 0.02]
            fbRng = [0.70; 0.70]
            K₀Rng = [0.0; 0.0]
            K₁Rng = [0.0; 0.0]
            K₂Rng = [0.0; 0.0]
            S1Rng = π/180*[2.0; 2.0]
            S2Rng = π/180*[2.0; 2.0]
            Tf₀Rng = [3.0; 3.0]
            TpRng = [1.7; 1.7]
            Tv₀Rng = [6.0; 6.0]
            TvLRng = [5.0; 5.0]
            # Fixed parameters
            γbC = [1.0; 1.0]
            γbCMat = Diagonal(γbC)    
        elseif name in ["flatPlate","NACA0002","NACA0006","NACA0012-GU","NACA0015","NACA0015-s","NACA0018","VERTOL23010","HeliosWingAirfoil","BWBAirfoil","cHALEairfoil"]
            # Bound Mach and corresponding compressibility factor
            Ma = max(0.001,min(0.8,Ma))
            β = sqrt(1-Ma^2) 
            # Mach-dependent parameters
            MaRng  = [0.001; 0.8]
            α₀NRng = π/180*[0.0; 0.0]
            α1₀Rng = π/180*[12.4; 16.0]
            δαRng = π/180*[2.0; 2.5]
            ϵₙRng = [0.70; 0.70] 
            ϵₘRng = [0.96; 0.96]
            cd₀Rng = [0.010; 0.010]
            cm₀Rng = [0.0; 0.0]
            ηRng = [0.95; 0.95]
            cn₁Rng = [1.25; 1.85]
            cnαRng = 2π/β*[1.0; 1.0]
            DfRng = [8.0; 8.0]
            E₀Rng = [0.0; 0.0]
            f₀Rng = [0.02; 0.02]
            fbRng = [0.70; 0.70]
            K₀Rng = [0.0; 0.0]
            K₁Rng = [-0.120; -0.120]
            K₂Rng = [0.04; 0.04]
            S1Rng = π/180*[3.0; 3.0]
            S2Rng = π/180*[1.5; 1.5]
            Tf₀Rng = [3.0; 3.0]
            TpRng = [1.7; 1.7]
            Tv₀Rng = [6.0; 6.0]
            TvLRng = [5.0; 5.0]
            # Fixed parameters
            γbC = [1.0; 1.0]
            γbCMat = Diagonal(γbC)
        else
            error("Airfoil not listed")
        end

        # Set effective angles
        α1₀Rng += α₀NRng

        # Interpolated values
        α₀N = LinearInterpolations.interpolate(MaRng,α₀NRng,Ma)
        α1₀ = LinearInterpolations.interpolate(MaRng,α1₀Rng,Ma)
        δα = LinearInterpolations.interpolate(MaRng,δαRng,Ma)
        ϵₙ = LinearInterpolations.interpolate(MaRng,ϵₙRng,Ma)
        ϵₘ = LinearInterpolations.interpolate(MaRng,ϵₘRng,Ma)
        η = LinearInterpolations.interpolate(MaRng,ηRng,Ma)
        cd₀ = LinearInterpolations.interpolate(MaRng,cd₀Rng,Ma)
        cm₀ = LinearInterpolations.interpolate(MaRng,cm₀Rng,Ma)
        cn₁ = LinearInterpolations.interpolate(MaRng,cn₁Rng,Ma)
        cnα = LinearInterpolations.interpolate(MaRng,cnαRng,Ma)
        Df = LinearInterpolations.interpolate(MaRng,DfRng,Ma)
        E₀ = LinearInterpolations.interpolate(MaRng,E₀Rng,Ma)
        f₀ = LinearInterpolations.interpolate(MaRng,f₀Rng,Ma)
        fb = LinearInterpolations.interpolate(MaRng,fbRng,Ma)
        K₀ = LinearInterpolations.interpolate(MaRng,K₀Rng,Ma)
        K₁ = LinearInterpolations.interpolate(MaRng,K₁Rng,Ma)
        K₂ = LinearInterpolations.interpolate(MaRng,K₂Rng,Ma)
        S1 = LinearInterpolations.interpolate(MaRng,S1Rng,Ma)
        S2 = LinearInterpolations.interpolate(MaRng,S2Rng,Ma)
        Tf₀ = LinearInterpolations.interpolate(MaRng,Tf₀Rng,Ma)
        Tp = LinearInterpolations.interpolate(MaRng,TpRng,Ma)
        Tv₀ = LinearInterpolations.interpolate(MaRng,Tv₀Rng,Ma)
        TvL = LinearInterpolations.interpolate(MaRng,TvLRng,Ma)

        # Dimensionalize time delay constants
        if U > 0 && b > 0
            Tf₀ *= b/U
            Tp *= b/U
            Tv₀ *= b/U
            TvL *= b/U
        end

        return new(α₀N,α1₀,δα,ϵₙ,ϵₘ,η,cd₀,cm₀,cn₁,cnα,Df,E₀,f₀,fb,K₀,K₁,K₂,S1,S2,Tf₀,Tp,Tv₀,TvL,γbC,γbCMat)
    end

end


# @with_kw mutable struct FlapParameters

#     FlapParameters composite type

# Fields
# - `cdδ::Real`
# - `cmδ::Real`
# - `cnδ::Real`
# 
@with_kw mutable struct FlapParameters

    cdδ::Real
    cmδ::Real
    cnδ::Real

    function FlapParameters(name::String; flapSiteID::Int64)

        # Validate
        @assert 0 < flapSiteID <= 100

        # Airfoil flap parameters
        if name in ["flatPlate","NACA0002","NACA0006"]
            if flapSiteID == 100
                cdδ = 0.0
                cmδ = 0.0
                cnδ = 0.0
            elseif flapSiteID == 75
                cdδ =  0.25
                cmδ = -0.25
                cnδ =  1.0
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["NACA0012"]
            if flapSiteID == 100
                cdδ = 0.0
                cmδ = 0.0
                cnδ = 0.0
            elseif flapSiteID == 75
                cdδ =  0.25
                cmδ = -0.25
                cnδ =  1.0
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["NACA0015","NACA0018"]
            if flapSiteID == 100
                cdδ = 0.0
                cmδ = 0.0
                cnδ = 0.0
            elseif flapSiteID == 75
                cdδ =  0.25
                cmδ = -0.25
                cnδ =  1.0
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["NACA23012A"]
            if flapSiteID == 100
                cdδ = 0.0
                cmδ = 0.0
                cnδ = 0.0
            elseif flapSiteID == 75
                cdδ =  0.0
                cmδ = -0.25
                cnδ =  1.0
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["HeliosWingAirfoil"]
            if flapSiteID == 100
                cdδ = 0.0
                cmδ = 0.0
                cnδ = 0.0
            elseif flapSiteID == 75
                cdδ =  0.0
                cmδ = -0.25
                cnδ =  1.0
            else
                error("Unavailable flap site ID")
            end    
        elseif name in ["HeliosPodAirfoil"]
            if flapSiteID == 100
                cdδ = 0.0
                cmδ = 0.0
                cnδ = 0.0
            elseif flapSiteID == 75
                cdδ =  0.0
                cmδ = -0.25
                cnδ =  1.0
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["BWBAirfoil"] 
            if flapSiteID == 100
                cdδ = 0.0
                cmδ = 0.0
                cnδ = 0.0
            elseif flapSiteID == 75
                cdδ =  0.0123
                cmδ = -0.4610
                cnδ =  3.5180
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["VERTOL23010"]
            if flapSiteID == 100
                cdδ = 0.0
                cmδ = 0.0
                cnδ = 0.0
            elseif flapSiteID == 75
                cdδ =  0.0
                cmδ = -0.25
                cnδ =  1.0
            else
                error("Unavailable flap site ID")
            end
        elseif name in ["cHALEairfoil"]
            if flapSiteID == 100
                cdδ = 0.0
                cmδ = 0.0
                cnδ = 0.0
            elseif flapSiteID == 75
                cdδ =  0.15
                cmδ = -0.35
                cnδ =  2.5
            else
                error("Unavailable flap site ID")
            end  
        else
            error("Airfoil not listed")
        end

        return new(cdδ,cmδ,cnδ)
    end
end


# Gets the coordinates of the airfoil
function get_airfoil_coordinates(name::String)

    if name in ["flatPlate","NACA0002","NACA0006","cHALEairfoil"]
        coords = [1.000000	-0.000000
                0.998459	0.000019
                0.993844	0.000074
                0.986185	0.000166
                0.975528	0.000292
                0.961940	0.000450
                0.945503	0.000638
                0.926320	0.000852
                0.904508	0.001089
                0.880203	0.001346
                0.853553	0.001620
                0.824724	0.001906
                0.793893	0.002200
                0.761249	0.002500
                0.726995	0.002801
                0.691342	0.003099
                0.654508	0.003391
                0.616723	0.003671
                0.578217	0.003937
                0.539230	0.004183
                0.500000	0.004405
                0.460770	0.004599
                0.421783	0.004759
                0.383277	0.004882
                0.345492	0.004963
                0.308658	0.004999
                0.273005	0.004987
                0.238751	0.004924
                0.206107	0.004809
                0.175276	0.004642
                0.146447	0.004424
                0.119797	0.004154
                0.095492	0.003837
                0.073680	0.003475
                0.054497	0.003072
                0.038060	0.002632
                0.024472	0.002158
                0.013815	0.001654
                0.006156	0.001125
                0.001541	0.000573
                0.000000	0.000000
                0.001541	-0.000573
                0.006156	-0.001125
                0.013815	-0.001654
                0.024472	-0.002158
                0.038060	-0.002632
                0.054497	-0.003072
                0.073680	-0.003475
                0.095492	-0.003837
                0.119797	-0.004154
                0.146447	-0.004424
                0.175276	-0.004642
                0.206107	-0.004809
                0.238751	-0.004924
                0.273005	-0.004987
                0.308658	-0.004999
                0.345492	-0.004963
                0.383277	-0.004882
                0.421783	-0.004759
                0.460770	-0.004599
                0.500000	-0.004405
                0.539230	-0.004183
                0.578217	-0.003937
                0.616723	-0.003671
                0.654508	-0.003391
                0.691342	-0.003099
                0.726995	-0.002801
                0.761249	-0.002500
                0.793893	-0.002200
                0.824724	-0.001906
                0.853553	-0.001620
                0.880203	-0.001346
                0.904508	-0.001089
                0.926320	-0.000852
                0.945503	-0.000638
                0.961940	-0.000450
                0.975528	-0.000292
                0.986185	-0.000166
                0.993844	-0.000074
                0.998459	-0.000019
                1.000000	0.000000]
    elseif name in ["NACA0012","NACA0012-GU"]
        coords = [1.000000	0
                0.999416	0.001342
                0.997666	0.001587
                0.994753	0.001994
                0.990685	0.002560
                0.985471	0.003280
                0.979123	0.004152
                0.971656	0.005169
                0.963087	0.006324
                0.953437	0.007611
                0.942728	0.009022
                0.930985	0.010549
                0.918235	0.012182
                0.904509	0.013914
                0.889837	0.015735
                0.874255	0.017635
                0.857800	0.019605
                0.840508	0.021635
                0.822421	0.023714
                0.803581	0.025834
                0.784032	0.027983
                0.763820	0.030152
                0.742992	0.032329
                0.721596	0.034506
                0.699682	0.036670
                0.677303	0.038811
                0.654509	0.040917
                0.631354	0.042978
                0.607892	0.044980
                0.584179	0.046912
                0.560268	0.048762
                0.536217	0.050516
                0.512082	0.052162
                0.487918	0.053687
                0.463783	0.055077
                0.439732	0.056320
                0.415822	0.057403
                0.392108	0.058314
                0.368646	0.059042
                0.345492	0.059575
                0.322698	0.059903
                0.300318	0.060017
                0.278404	0.059910
                0.257008	0.059576
                0.236180	0.059008
                0.215968	0.058205
                0.196419	0.057164
                0.177579	0.055886
                0.159492	0.054372
                0.142201	0.052625
                0.125745	0.050651
                0.110163	0.048457
                0.095492	0.046049
                0.081765	0.043437
                0.069015	0.040631
                0.057272	0.037641
                0.046563	0.034479
                0.036913	0.031156
                0.028344	0.027683
                0.020877	0.024071
                0.014529	0.020330
                0.009315	0.016471
                0.005247	0.012501
                0.002334	0.008429
                0.000584	0.004260
                0.000000	0.000000
                0.000584	-0.004260
                0.002334	-0.008429
                0.005247	-0.012501
                0.009315	-0.016471
                0.014529	-0.020330
                0.020877	-0.024071
                0.028344	-0.027683
                0.036913	-0.031156
                0.046563	-0.034479
                0.057272	-0.037641
                0.069015	-0.040631
                0.081765	-0.043437
                0.095492	-0.046049
                0.110163	-0.048457
                0.125745	-0.050651
                0.142201	-0.052625
                0.159492	-0.054372
                0.177579	-0.055886
                0.196419	-0.057164
                0.215968	-0.058205
                0.236180	-0.059008
                0.257008	-0.059576
                0.278404	-0.059910
                0.300318	-0.060017
                0.322698	-0.059903
                0.345492	-0.059575
                0.368646	-0.059042
                0.392108	-0.058314
                0.415822	-0.057403
                0.439732	-0.056320
                0.463783	-0.055077
                0.487918	-0.053687
                0.512082	-0.052162
                0.536217	-0.050516
                0.560268	-0.048762
                0.584179	-0.046912
                0.607892	-0.044980
                0.631354	-0.042978
                0.654509	-0.040917
                0.677303	-0.038811
                0.699682	-0.036670
                0.721596	-0.034506
                0.742992	-0.032329
                0.763820	-0.030152
                0.784032	-0.027983
                0.803581	-0.025834
                0.822421	-0.023714
                0.840508	-0.021635
                0.857800	-0.019605
                0.874255	-0.017635
                0.889837	-0.015735
                0.904509	-0.013914
                0.918235	-0.012182
                0.930985	-0.010549
                0.942728	-0.009022
                0.953437	-0.007611
                0.963087	-0.006324
                0.971656	-0.005169
                0.979123	-0.004152
                0.985471	-0.003280
                0.990685	-0.002560
                0.994753	-0.001994
                0.997666	-0.001587
                0.999416	-0.001342
                1.000000	0]
    elseif name in ["NACA0015","NACA0015-s"]
        coords = [1.000000	0.001580
                0.950000	0.010080
                0.900000	0.018100
                0.800000	0.032790
                0.700000	0.045800
                0.600000	0.057040
                0.500000	0.066170
                0.400000	0.072540
                0.300000	0.075020
                0.250000	0.074270
                0.200000	0.071720
                0.150000	0.066820
                0.100000	0.058530
                0.075000	0.052500
                0.050000	0.044430
                0.025000	0.032680
                0.012500	0.023670
                0.000000	0.000000
                0.012500	-0.023670
                0.025000	-0.032680
                0.050000	-0.044430
                0.075000	-0.052500
                0.100000	-0.058530
                0.150000	-0.066820
                0.200000	-0.071720
                0.250000	-0.074270
                0.300000	-0.075020
                0.400000	-0.072540
                0.500000	-0.066170
                0.600000	-0.057040
                0.700000	-0.045800
                0.800000	-0.032790
                0.900000	-0.018100
                0.950000	-0.010080
                1.000000	-0.001580]
    elseif name in ["NACA0018","HeliosPodAirfoil"]
        coords = [1.000000	0.000000
                0.990000	0.001890
                0.950000	0.012100
                0.900000	0.021720
                0.800000	0.039350
                0.700000	0.054960
                0.600000	0.068450
                0.500000	0.079410
                0.400000	0.087050
                0.300000	0.090030
                0.250000	0.089120
                0.200000	0.086060
                0.150000	0.080180
                0.100000	0.070240
                0.075000	0.063000
                0.050000	0.053320
                0.025000	0.039220
                0.012500	0.028410
                0.000000	0.000000
                0.012500	-0.028410
                0.025000	-0.039220
                0.050000	-0.053320
                0.075000	-0.063000
                0.100000	-0.070240
                0.150000	-0.080180
                0.200000	-0.086060
                0.250000	-0.089120
                0.300000	-0.090030
                0.400000	-0.087050
                0.500000	-0.079410
                0.600000	-0.068450
                0.700000	-0.054960
                0.800000	-0.039350
                0.900000	-0.021720
                0.950000	-0.012100
                0.990000	-0.001890
                1.000000	0.000000]
    elseif name in ["NACA23012A","HeliosWingAirfoil","BWBAirfoil"]
        coords = [0            0
                -4.4000e-04   8.0200e-03
                3.3700e-03   1.6940e-02
                1.1660e-02   2.6570e-02
                2.4540e-02   3.6510e-02
                4.2070e-02   4.6260e-02
                6.4130e-02   5.5230e-02
                9.0480e-02   6.2860e-02
                1.2069e-01   6.8760e-02
                1.5421e-01   7.2760e-02
                1.9042e-01   7.5030e-02
                2.2902e-01   7.6030e-02
                2.7060e-01   7.5970e-02
                3.1507e-01   7.4790e-02
                3.6224e-01   7.2410e-02
                4.1195e-01   6.8720e-02
                4.6399e-01   6.3650e-02
                5.1816e-01   5.7250e-02
                5.7424e-01   4.9640e-02
                6.3202e-01   4.1030e-02
                6.9125e-01   3.1690e-02
                7.5169e-01   2.2020e-02
                8.1310e-01   1.2570e-02
                8.7521e-01   4.2200e-03
                9.3773e-01  -1.2500e-03
                1.0003e+00   5.1000e-04
                1.0003e+00  -5.0000e-04
                9.3768e-01  -1.7050e-02
                8.7515e-01  -2.5870e-02
                8.1306e-01  -3.1470e-02
                7.5169e-01  -3.5440e-02
                6.9128e-01  -3.8430e-02
                6.3209e-01  -4.0770e-02
                5.7436e-01  -4.2640e-02
                5.1831e-01  -4.4010e-02
                4.6418e-01  -4.4980e-02
                4.1216e-01  -4.5470e-02
                3.6247e-01  -4.5400e-02
                3.1530e-01  -4.4710e-02
                2.7083e-01  -4.3330e-02
                2.2925e-01  -4.1230e-02
                1.9077e-01  -3.8380e-02
                1.5631e-01  -3.5080e-02
                1.2588e-01  -3.1800e-02
                9.9100e-02  -2.8740e-02
                7.5640e-02  -2.5880e-02
                5.5290e-02  -2.3080e-02
                3.7910e-02  -2.0080e-02
                2.3540e-02  -1.6580e-02
                1.2290e-02  -1.2260e-02
                4.3600e-03  -6.8100e-03
                        0            0]
    elseif name in ["VERTOL23010"]
        coords = [1.000000  0.000000
                0.949997  0.005670
                0.899993  0.011540
                0.849990  0.017010
                0.799986  0.022980
                0.749983  0.028150
                0.699980  0.033720
                0.649976  0.039490
                0.599973  0.044660
                0.549970  0.049330
                0.499968  0.053700
                0.449965  0.057470
                0.399964  0.060540
                0.349962  0.063010
                0.299961  0.064980
                0.249960  0.065950
                0.199960  0.066220
                0.149961  0.064190
                0.124963  0.061275
                0.099966  0.057060
                0.074969  0.051445
                0.049974  0.042830
                0.029980  0.033118
                0.019984  0.027012
                0.009988  0.019206
                0.000000  0.000000
                0.010009 -0.015094
                0.020010 -0.017588
                0.030011 -0.018982
                0.050012 -0.020570
                0.075013 -0.022655
                0.100015 -0.024740
                0.125016 -0.026825
                0.150017 -0.029110
                0.200019 -0.032680
                0.250021 -0.035350
                0.300022 -0.036820
                0.350022 -0.037090
                0.400022 -0.036760
                0.450021 -0.035530
                0.500020 -0.033800
                0.550019 -0.031470
                0.600017 -0.028740
                0.650015 -0.025410
                0.700013 -0.022180
                0.750011 -0.018550
                0.800009 -0.014920
                0.850007 -0.011190
                0.900004 -0.007560
                0.950002 -0.003830
                1.000000  0.000000]                        
    else
        error("Airfoil coordinates not listed")
    end

    return coords
end


#

    # Airfoil composite type

#
@with_kw mutable struct Airfoil

    name::String
    coordinates::Matrix{Float64}
    attachedFlowParameters::AttachedFlowParameters
    parametersBLi::BLiParameters
    parametersBLo::BLoParameters
    flapParameters::FlapParameters

end
export Airfoil


"""
    create_Airfoil(; kwargs...)

Airfoil constructor

# Arguments
- `name::String`: name of the airfoil
- `Re::Real`: Reynolds number
- `Ma::Real`: Mach number
- `U::Real`: relative airspeed
- `b::Real`: semichord
"""
function create_Airfoil(; name::String,Re::Real=0,Ma::Real=0,U::Real=0,b::Real=0)

    coordinates = get_airfoil_coordinates(name)
    attachedFlowParameters = AttachedFlowParameters(name,Re=Re,Ma=Ma)
    parametersBLi = BLiParameters(name,Re=Re,Ma=Ma,U=U,b=b)
    parametersBLo = BLoParameters(name,Re=Re,Ma=Ma,U=U,b=b)
    flapParameters = FlapParameters(name,flapSiteID=100)

    return Airfoil(name,coordinates,attachedFlowParameters,parametersBLi,parametersBLo,flapParameters)
end
export create_Airfoil


"""
    create_flapped_Airfoil(; kwargs...)

Airfoil constructor (with a trailing-edge flap) 

# Arguments
- `name::String`: name of the airfoil
- `flapSiteID::Int64`: flap site ID
- `Re::Real`: Reynolds number
- `Ma::Real`: Mach number
- `U::Real`: relative airspeed
- `b::Real`: semichord
"""
function create_flapped_Airfoil(; name::String,flapSiteID::Int64,Re::Real=0,Ma::Real=0,U::Real=0,b::Real=0)

    coordinates = get_airfoil_coordinates(name)
    attachedFlowParameters = AttachedFlowParameters(name,Re=Re,Ma=Ma)
    parametersBLi = BLiParameters(name,Re=Re,Ma=Ma,U=U,b=b)
    parametersBLo = BLoParameters(name,Re=Re,Ma=Ma,U=U,b=b)
    flapParameters = FlapParameters(name,flapSiteID=flapSiteID)
    
    return Airfoil(name,coordinates,attachedFlowParameters,parametersBLi,parametersBLo,flapParameters)
end
export create_flapped_Airfoil


"""
    update_Airfoil_params!(airfoil::Airfoil; kwargs...)

Updates the airfoil parameters

# Arguments
- `airfoil::Airfoil`

# Keyword arguments
- `Re::Real`: Reynolds number
- `Ma::Real`: Mach number
- `U::Real`: reference relative airspeed
- `b::Real`: reference semichord
"""
function update_Airfoil_params!(airfoil::Airfoil; Re::Real=0,Ma::Real=0,U::Real,b::Real)

    airfoil.attachedFlowParameters = AttachedFlowParameters(airfoil.name,Re=Re,Ma=Ma)
    airfoil.parametersBLi = BLiParameters(airfoil.name,Re=Re,Ma=Ma,U=U,b=b)
    airfoil.parametersBLo = BLoParameters(airfoil.name,Re=Re,Ma=Ma,U=U,b=b)

end
export update_Airfoil_params!


# Sample airfoils
flatPlate = create_Airfoil(name="flatPlate")
NACA0002 = create_Airfoil(name="NACA0002")
NACA0006 = create_Airfoil(name="NACA0006")
NACA0012 = create_Airfoil(name="NACA0012")
NACA0015 = create_Airfoil(name="NACA0015")
NACA0018 = create_Airfoil(name="NACA0018")
VERTOL23010 = create_Airfoil(name="VERTOL23010")
NACA23012A = create_Airfoil(name="NACA23012A")
HeliosWingAirfoil = create_Airfoil(name="HeliosWingAirfoil")
HeliosPodAirfoil = create_Airfoil(name="HeliosPodAirfoil")
BWBAirfoil = create_Airfoil(name="BWBAirfoil")
cHALEairfoil = create_Airfoil(name="cHALEairfoil")
export flatPlate, NACA0002, NACA0006, NACA0012, NACA0015, NACA0018, NACA23012A, VERTOL23010, HeliosWingAirfoil, HeliosPodAirfoil, BWBAirfoil, cHALEairfoil