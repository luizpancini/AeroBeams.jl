# # Creating a HALE aircraft model
# This example illustrates how to create a conventional aircraft model, the high-altitude long-endurance (HALE) aircraft described by [Patil, Hodges and Cesnik](https://doi.org/10.2514/2.2738):

#md # ![](../assets/cHALE.png)
#md # *HALE model geometry* by [Patil, Hodges and Cesnik](https://doi.org/10.2514/2.2738)

#md # ![](../assets/cHALEdata.png)
#md # *HALE model data* by [Patil, Hodges and Cesnik](https://doi.org/10.2514/2.2738)

#md # !!! tip
#md #     The code for this example is available [here](https://github.com/luizpancini/AeroBeams.jl/blob/main/test/examples/steadyAeroelastic/cHALEmodel.jl).

# The contents of the present example are wrapped into the function [`create_cHALE`](@ref create_cHALE), in which you may check all the keyword arguments that we'll define below.

## Define input arguments
using AeroBeams

altitude=20e3
aeroSolver=Indicial()
derivationMethod=AD()
flapSolver=TableLookup()
gustSolver=IndicialGust(:Kussner)
stabilizersAero=true
includeVS=true
nElemWing=20
nElemHorzStabilizer=10
nElemTailBoom=10
nElemVertStabilizer=5
∞=1e12
stiffnessFactor=1
torsionStiffnessFactor=1
k1=0
k2=0
airspeed=0
δElevIsTrimVariable=false
δAilIsTrimVariable=false
δRuddIsTrimVariable=false
thrustIsTrimVariable=false
δElev=nothing
δAil=nothing
δRudd=nothing
thrust=0
g=standard_atmosphere(altitude).g
wingCd0=0
wingcnδ=2.5
wingcmδ=-0.35
wingcdδ=0.15
stabsCd0=0
stabscnδ=2.5
stabscmδ=-0.35
stabscdδ=0.15
updateAirfoilParameters=false
hasInducedDrag=true
#md nothing #hide

# ### Validation
# The first step is the validation of the inputs, for which reasonable values are necessary.
## Validate
@assert iseven(nElemWing) "nElemWing must be even"
@assert iseven(nElemHorzStabilizer) "nElemHorzStabilizer must be even"
@assert ∞ > 1e8 "rigid structural properties must be greater than 1e8"
@assert stiffnessFactor > 0 "stiffnessFactor must be greater than 0"
@assert wingCd0 >= 0 "wingCd0 must be non-negative"
@assert stabsCd0 >= 0 "stabsCd0 must be non-negative"
@assert airspeed >= 0 "airspeed must be non-negative"
@assert g >= 0 "g must be non-negative"
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
#md nothing #hide

# ### Atmosphere
# Let's initialize an instance of the International Standard Atmosphere (ISA) at the specified altitude using the function [`standard_atmosphere`](@ref standard_atmosphere).
## Atmosphere
atmosphere = standard_atmosphere(altitude)
#md nothing #hide

# ### Airfoils
# The airfoils of the wing and stabilizers are created with the function [`create_Airfoil`](@ref create_Airfoil). However, the aerodynamic parameters of the airfoils are modified with respect to a reference flat plate airfoil by first creating a new instance of the airfoils' databases with the function [`modified_database`](@ref modified_database).
## Airfoils
wFlapSiteID = 75
dbWingAirfoil = modified_database(
    "FlatPlate";
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
    "FlatPlate";
    attached = (
        cd₀ = stabsCd0,
    ),
    flap = (
        sFlapSiteID => create_FlapParameters(stabscdδ,stabscmδ,stabscnδ),
    )
)
sAirfoil = create_Airfoil(database=dbStabAirfoil,Ma=airspeed/atmosphere.a)
#md nothing #hide

# ### Wing aerodynamic surface
# In AeroBeams, each beam can have an aerodynamic surface attached to it. The wing surface is defined by an airfoil, `wAirfoil`, the normalized position of the beam's reference line (spar) on the chord, `wNormSparPos`, the normalized position of the flap (aileron) hinge on the chord, `wNormFlapPos`, and the normalized endpoints of the flap (aileron) span over the beam's arclength, `wNormFlapSpan`. To create an aerodynamic surface, we use the function [`create_AeroSurface`](@ref create_AeroSurface).
## Wing surface
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
#md nothing #hide

# ### Wing beam (spar)
# Next, we'll define the structural properties of the beam. Notice that we have the option of specifying a stiffness factor, `stiffnessFactor`, for the relevant properties of the wing. We then use the built-in functions [`isotropic_stiffness_matrix`](@ref isotropic_stiffness_matrix) and [`inertia_matrix`](@ref inertia_matrix) to create the sectional stiffness (force-strain relations) and inertia (velocity-displacement relation) matrices.
## Beam properties
Lw = 16
wGJ,wEIy,wEIz = 1e4, 2e4, 4e6
wGJ,wEIy,wEIz = AeroBeams.multiply_scalars!(stiffnessFactor, wGJ,wEIy,wEIz)
wGJ *= torsionStiffnessFactor
wρA,wρIs = 0.75,0.1
wρIy,wρIz = (wEIy/wEIz)*wρIs,(1-wEIy/wEIz)*wρIs
Swing = isotropic_stiffness_matrix(∞=∞,GJ=wGJ,EIy=wEIy,EIz=wEIz)
Iwing = inertia_matrix(ρA=wρA,ρIy=wρIy,ρIz=wρIz)
#md nothing #hide

# In order to position the center of the wing at the origin of the frame A, we need the position of the first node of the assembly, at the left wingtip. That `firstNodePosition` vector will be a function of the undeformed wing curvature `k2` (which is zero in the original model). Similarly, the initial angle of twist of the wing, `ψ`, is a function of the undeformed wing twist curvature, `k1` (also zero in the original model).
## Initial position for first node of left wing
ρ = 1/k2
θ = Lw/ρ
x = ρ*sin(θ)
z = ρ*(1-cos(θ))
firstNodePosition = iszero(k2) ? 
                    [-Lw; 0; 0] : 
                    [-x; 0; -z]

## Initial angle of twist
r = 1/k1
ψ = Lw/r
#md nothing #hide

# To create the beams of the left and right wings, we use the function [`create_Beam`](@ref create_Beam). The rotation parameters from the underformed beam basis (b) to the global basis (A), namely `p0`, are set for the left wing such that the left and right wings become symmetric about the x1-plane of basis A. We also create an antissymetric link between the ailerons of the left and right wings, using the function [`create_FlapLink`](@ref create_FlapLink).
## Wing beams
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

## Link wing ailerons
aileronLink = create_FlapLink(masterBeam=rightWing,slaveBeams=[leftWing],δMultipliers=[-1])
#md nothing #hide

# ### Payload
# The payload is modeled as a point inertia at the center of the wing spar. It is created with the constructor [`create_PointInertia`](@ref create_PointInertia), and then added to the center of the wing through the function [`add_point_inertia_to_beam!`](@ref add_point_inertia_to_beam!).
## Payload
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
#md nothing #hide

# ### Tail boom
# The tail boom is assumed to be a rigid beam, devoid of aerodynamic surfaces, and rigidly connected to wing spar. Some of its sectional inertia matrix properties are assumed with reasonable values, since they are not given in the reference. 
## Tail boom
Lt = 10
tρA,tρIy,tρIz = 0.08,wρIy/10,wρIz/10
Stail = isotropic_stiffness_matrix(∞=∞)
Itail = inertia_matrix(ρA=tρA,ρIy=tρIy,ρIz=tρIz)
tailBoom = create_Beam(
        name="tailBoom",
        length=Lt,
        nElements=nElemTailBoom,
        S=[Stail],
        I=[Itail],
        p0=transform_to_WM([-π/2; 0; 0], from=:E321),
        )
#md nothing #hide

# ### Horizontal stabilizer
# The horizontal stabilizer's aerodynamic surface is constructed similarly to the wing's, though it depends on whether the elevator deflecton is a trim variable or not (through the flags `δElevIsTrimVariable` and `δElevIsInput`).
## Horizontal stabilizer surface
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
#md nothing #hide

# The horizontal stabilizer beam (spar) is also created analogously to the wing, albeit now we input the `firstNodeOffset` keyword to [`create_Beam`](@ref create_Beam) in order to the define the position of the stabilizer's first node at the left tip.
## Horizontal stabilizer beam
Lh = 5
hρA,hρIy,hρIz = 0.08,wρIy/10,wρIz/10
Shorz = isotropic_stiffness_matrix(∞=∞)
Ihorz = inertia_matrix(ρA=hρA,ρIy=hρIy,ρIz=hρIz)
horzStabilizer = create_Beam(
        name="horzStabilizer",
        length=Lh,
        firstNodeOffset=[-Lh/2; 0; 0],
        nElements=nElemHorzStabilizer,
        S=[Shorz],
        I=[Ihorz],
        surface=hsSurf
        )
#md nothing #hide

# ### Vertical stabilizer
# The construction of the vertical stabilizer is again similar to that of the wing and horizontal counterpart.
## Vertical stabilizer surface
vChord = 0.5
vNormSparPos = 0.5
vNormFlapPos = sFlapSiteID/100
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
#md nothing #hide

## Vertical stabilizer beam
Lv = 2.5
vρA,vρIy,vρIz = 0.08,hρIy,hρIz
Svert = isotropic_stiffness_matrix(∞=∞)
Ivert = inertia_matrix(ρA=vρA,ρIy=vρIy,ρIz=vρIz)
vertStabilizer = create_Beam(
        name="vertStabilizer",
        length=Lv,
        nElements=nElemVertStabilizer,
        S=[Svert],
        I=[Ivert],
        p0=transform_to_WM([0; -π/2; 0], from=:E321),
        surface=vsSurf
        )
#md nothing #hide

# ### Thrust force
# The thrust force of the vehicle is modeled as a follower concentrated force acting at the wing's center node (node 1 of `rightWing`), initially pointing in the x2-direction of the undeformed beam basis. Whether it is to take a time-dependent user input value or to be a trim variable depends on the user's choice. It is created using the function [`create_BC`](@ref create_BC).
## Propeller thrust force 
propThrust = create_BC(
        name="propThrust",
        beam=rightWing,
        node=1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )
#md nothing #hide

# ### Model assembly
# The beams composing the model depend on whether the vertical stabilizer is to be included, through the flag `includeVS`. The connections between those beams have to be set through the `connections` dictionary, in which each key-value pair defines the beams and the corresponding nodes to be connected. The model is created with the appropriate inputs using the function [`create_Model`](@ref create_Model). Notice that the keyword argument `firstNodePosition` defines the position of the first node of the first beam of the model, namely, the node at the left wingtip. Also, the velocity of the body-attached global frame A with respect to an inertial frame is set through `v_A`.
## Beams of the model
beams = includeVS ? 
        [leftWing,rightWing,tailBoom,horzStabilizer,vertStabilizer] : 
        [leftWing,rightWing,tailBoom,horzStabilizer]

## Beam connections
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

## Aircraft model
cHALE = create_Model(
        name="cHALE",
        beams=beams,
        BCs=[propThrust],
        firstNodePosition=firstNodePosition,
        gravityVector=[0; 0; -g],
        connections=connections,
        v_A=[0; airspeed; 0],
        flapLinks=[aileronLink],
        atmosphere=atmosphere,
        )
#md nothing #hide

#md # Finally, we may visually check the undeformed assembly through the function [`plot_undeformed_assembly`](@ref plot_undeformed_assembly):
## Visualization
#md plot_undeformed_assembly(cHALE, view=(30,30))
#md nothing #hide
#md savefig("cHALE.svg") #hide
#md nothing #hide
#md # ![](cHALE.svg)
#md nothing #hide

## Plots
p1 = plot_undeformed_assembly(cHALE, view=(30,30))
display(p1)

println("Finished cHALEmodel.jl") #src