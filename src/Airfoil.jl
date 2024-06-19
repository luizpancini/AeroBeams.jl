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
@with_kw mutable struct Airfoil

    Airfoil composite type

# Fields
- name::String
- 
"""
@with_kw mutable struct Airfoil

    name::String
    attachedFlowParameters::AttachedFlowParameters

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

    return Airfoil(name,attachedFlowParameters)
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

    return Airfoil(name,attachedFlowParameters)
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