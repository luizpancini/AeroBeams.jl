"""
@with_kw mutable struct Model

Model composite type

# Fields
- name::String = name of the model
- units::UnitsSystem = unit system (length, force, angle and frequency) 
- beams::Vector{Beam} = beams that compose the model's assembly
- initialPosition::Vector{<:Number} = initial position of the first node of the first beam of the model, defined in basis A 
- gravityVector::Vector{<:Number} = gravity vector, defined in the I (inertial) frame
- BCs::Dict{Int64,Int64} = boundary conditions applied to the model
- nElementsTotal::Int64 = total number of elements of the beam's assembly
- nNodesTotal::Int64 = total number of nodes of the beam's assembly
- elementNodes::Vector{Vector{Int64}} = a nElementsTotal x 2 matrix in which the row represents the element, and the columns are the corresponding element's nodes
- p_A0::Vector{Float64} = initial Euler 3-2-1 rotation parameters that bring basis I to basis A
- R_A::Matrix{Float64}= 1.0I(3) = initial rotation tensor that brings basis I to basis A, resolved in basis A
- nTrimVariables::Int64 = number of trim variables
- trimLoadsLinks::Vector{TrimLoadsLink} = trim links among loads
# Notes
 - The default is an empty model
"""
@with_kw mutable struct Model

    # Primary (inputs for the definition of the model)
    # ------------------------------------------------
    name::String
    units::UnitsSystem
    beams::Vector{Beam}
    initialPosition::Vector{<:Number}
    gravityVector::Vector{<:Number}
    BCs::Vector{BC}
    p_A0::Vector{Float64}
    u_A::Union{Vector{<:Number},<:Function,Nothing}
    v_A::Union{Vector{<:Number},<:Function,Nothing}
    ω_A::Union{Vector{<:Number},<:Function,Nothing}
    vdot_A::Union{Vector{<:Number},<:Function,Nothing}
    ωdot_A::Union{Vector{<:Number},<:Function,Nothing}
    altitude::Union{Nothing,Number}
    atmosphere::Union{Nothing,Atmosphere}
    gust::Union{Nothing,Gust}
    trimLoadsLinks::Vector{TrimLoadsLink}
    flapLinks::Vector{FlapLink}

    # Secondary (outputs from model creation)
    # ---------------------------------------
    elements::Vector{Element} = Vector{Element}()
    nElementsTotal::Int64 = 0
    nNodesTotal::Int64 = 0
    elementNodes::Vector{Vector{Int64}} = Vector{Vector{Int64}}()
    r_n::Vector{Vector{Float64}} = Vector{Vector{Float64}}()
    BCedNodes::Vector{Int64} = Vector{Int64}()
    specialNodes::Vector{SpecialNode} = Vector{SpecialNode}()
    specialNodesGlobalIDs::Vector{Int64} = Vector{Int64}()
    systemOrder::Int64 = 0
    forceScaling::Float64 = 1.0
    R_A::Matrix{Float64} = I3
    R_AT::Matrix{Float64} = I3
    skipValidationMotionBasisA::Bool = false
    nTrimVariables::Int64 = 0
    mass::Number = 0
    centerOfMass::Vector{<:Number} = zeros(3)
    I::Vector{<:Number} = zeros(3)

end
export Model


# Constructor
function create_Model(;name::String="",units::UnitsSystem=UnitsSystem(),beams::Vector{Beam},initialPosition::Vector{<:Number}=zeros(3),gravityVector::Vector{<:Number}=zeros(3),BCs::Vector{BC}=Vector{BC}(),p_A0::Vector{Float64}=zeros(3),u_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,v_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,ω_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,vdot_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,ωdot_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,altitude::Union{Nothing,Number}=nothing,atmosphere::Union{Nothing,Atmosphere}=nothing,gust::Union{Nothing,Gust}=nothing,trimLoadsLinks::Vector{TrimLoadsLink}=Vector{TrimLoadsLink}(),flapLinks::Vector{FlapLink}=Vector{FlapLink}()) 
    
    # Initialize 
    self = Model(name=name,units=units,beams=beams,initialPosition=initialPosition,gravityVector=gravityVector,BCs=BCs,p_A0=p_A0,u_A=u_A,v_A=v_A,ω_A=ω_A,vdot_A=vdot_A,ωdot_A=ωdot_A,altitude=altitude,atmosphere=atmosphere,gust=gust,trimLoadsLinks=trimLoadsLinks,flapLinks=flapLinks)

    # Update  
    update_model!(self)

    return self
end
export create_Model


"""
update_model!(model::Model)

Updates the model with its current settings

# Arguments
- model::Model 
"""
function update_model!(model::Model)

    @unpack beams,BCs = model

    # Validate inputs
    validate_model!(model)

    # Assemble the beams into model  
    assemble_model!(model,beams)

    # Compute inertia properties of the assembly
    inertia_properties!(model)

    # Set atmosphere, if applicable
    set_atmosphere!(model)

    # Update number of gust states 
    update_number_gust_states!(model) 

    # Update the nodes' global IDs of the loads trim links
    update_loads_trim_links_global_ids!(model)

    # Update the global IDs of the springs' nodes
    update_spring_nodes_global_ids!(model::Model)

    # Update linked flap deflections
    update_linked_flap_deflections!(model)

    # Set BCs on model 
    set_BCs!(model,BCs)

    # Get rotation tensors from basis A
    initialize_basis_A_rotation!(model)

    # Get special nodes
    get_special_nodes!(model)

    # Update initial conditions
    update_initial_conditions!(model)

    # Get system indices
    get_system_indices!(model)

    return model

end
export update_model!


"""
validate_model!(model::Model)

Validates the inputs to the model

# Arguments
- model::Model
"""
function validate_model!(model::Model)

    @unpack units,initialPosition,gravityVector,p_A0,altitude = model

    # Validate unit system
    validate_unit_system(units)

    # Validate initial position of the first beam of the model
    @assert length(initialPosition) == 3

    # Validate gravity vector
    @assert length(gravityVector) == 3

    # Validate initial rotation of basis A
    @assert length(p_A0) == 3

    # Validate basis A displacements, velocities and accelerations
    validate_and_update_motion_basis_A!(model)

    # Validate altitude
    if !isnothing(altitude) 
        @assert altitude >= 0
    end

end


"""
validate_and_update_motion_basis_A!(model::Model)

Validates and updates the motion variables of basis A

# Arguments
- model::Model
"""
function validate_and_update_motion_basis_A!(model::Model)

    @unpack u_A,v_A,ω_A,vdot_A,ωdot_A,skipValidationMotionBasisA = model

    # Check flag
    if skipValidationMotionBasisA
        return
    end

    # Initialize TFs
    uIsInput,vIsInput,ωIsInput,vdotIsInput,ωdotIsInput = false,false,false,false,false

    # Assert variables' lengths and update TFs 
    if u_A isa Function && isequal(u_A(1),(t -> zeros(3))(1))
        uIsInput = false
    elseif !isnothing(u_A)
        if u_A isa Function
            @assert length(u_A(0)) == 3
        else
            @assert length(u_A) == 3
        end
        uIsInput = true
    else
        u_A = t -> zeros(3)
    end
    if v_A isa Function && isequal(v_A(1),(t -> zeros(3))(1))
        vIsInput = false 
    elseif !isnothing(v_A)
        if v_A isa Function
            @assert length(v_A(0)) == 3
        else
            @assert length(v_A) == 3
        end
        vIsInput = true   
    else
        v_A = t -> zeros(3)    
    end
    if ω_A isa Function && isequal(ω_A(1),(t -> zeros(3))(1))
        ωIsInput = false 
    elseif !isnothing(ω_A)
        if ω_A isa Function
            @assert length(ω_A(0)) == 3
        else
            @assert length(ω_A) == 3
        end
        ωIsInput = true   
    else
        ω_A = t -> zeros(3)    
    end
    if vdot_A isa Function && isequal(vdot_A(1),(t -> zeros(3))(1))
        vdotIsInput = false
    elseif !isnothing(vdot_A)
        if vdot_A isa Function
            @assert length(vdot_A(0)) == 3
        else
            @assert length(vdot_A) == 3
        end
        vdotIsInput = true    
    else
        vdot_A = t -> zeros(3)    
    end
    if ωdot_A isa Function && isequal(ωdot_A(1),(t -> zeros(3))(1))
        ωdotIsInput = false
    elseif !isnothing(ωdot_A)
        if ωdot_A isa Function
            @assert length(ωdot_A(0)) == 3
        else
            @assert length(ωdot_A) == 3
        end
        ωdotIsInput = true    
    else
        ωdot_A = t -> zeros(3)    
    end

    # Check consistency
    @assert !(uIsInput && vIsInput) "Both u_A and v_A cannot be input at the same time"
    @assert !(uIsInput && vdotIsInput) "Both u_A and vdot_A cannot be input at the same time"
    @assert !(vIsInput && vdotIsInput) "Both v_A and vdot_A cannot be input at the same time"
    @assert !(ωIsInput && ωdotIsInput) "Both ω_A and ωdot_A cannot be input at the same time"

    # Set motion variables as functions of time, if they aren't already
    if !(u_A isa Function)
        u_A_const = deepcopy(u_A)
        u_A = t -> u_A_const
    end
    if !(v_A isa Function)
        v_A_const = deepcopy(v_A)
        v_A = t -> v_A_const
    end
    if !(ω_A isa Function)
        ω_A_const = deepcopy(ω_A)
        ω_A = t -> ω_A_const
    end
    if !(vdot_A isa Function)
        vdot_A_const = deepcopy(vdot_A)
        vdot_A = t -> vdot_A_const
    end
    if !(ωdot_A isa Function)
        ωdot_A_const = deepcopy(ωdot_A)
        ωdot_A = t -> ωdot_A_const
    end

    # If u_A(t) was input, get v_A(t) and vdot_A(t)
    if uIsInput
        v_A = t -> ForwardDiff.derivative(u_A, t)
        vdot_A = t -> ForwardDiff.derivative(v_A, t)
    end

    # If v_A(t) was input, get u_A(t) and vdot_A(t)
    if vIsInput
        u_A = t -> quadgk(v_A, 0, t)[1]
        vdot_A = t -> ForwardDiff.derivative(v_A, t)
    end

    # If vdot_A(t) was input, get u_A(t) and v_A(t)
    if vdotIsInput
        v_A = t -> quadgk(vdot_A, 0, t)[1]
        u_A = t -> quadgk(v_A, 0, t)[1]
    end

    # If ω_A(t) was input, get ωdot_A(t)
    if ωIsInput
        ωdot_A = t -> ForwardDiff.derivative(ω_A, t)
    end

    # If ωdot_A(t) was input, get ω_A(t)
    if ωdotIsInput
        ω_A = t -> quadgk(ωdot_A, 0, t)[1]
    end

    @pack! model = u_A,v_A,ω_A,vdot_A,ωdot_A

end


"""
assemble_model!(model::Model,beams::Vector{Beam})

Loads the assembly of beams into the model

# Arguments
- model::Model = model at hand
- beams::Vector{Beam} = beams to load into the model
# Notes
"""
function assemble_model!(model::Model,beams::Vector{Beam})

    # Reset to a model empty of beams and elements
    model.beams = Vector{Beam}()
    model.elements = Vector{Element}()

    # Reset assembly variables
    nElementsTotal = 0
    nNodesTotal = 0
    elementNodes = Vector{Vector{Int64}}()
    r_n = Vector{Vector{Float64}}()
    forceScaling = 1.0
    S = Vector{Matrix{Float64}}()
    nodeCount = 0

    # Check input
    if isempty(beams)
        @pack! model = beams,nElementsTotal,nNodesTotal,elementNodes,r_n,forceScaling
        return 
    end

    # Loop over beams
    for (ID,beam) in enumerate(beams)

        @unpack nElements,elements,connectedBeams,connectedNodesThis,connectedNodesOther = beam

        # Update beam ID on assembly
        @pack! beam = ID

        # Adjust connectedBeams in case it is empty
        if isnothing(connectedBeams)
            connectedBeams = Vector{Beam}()
        end

        # Initialize range of global IDs of the elements of the current beam
        elementRange = Vector(1:nElements)

        # Update the element range for the current beam in the assembly
        elementRangeIncrement = 0
        for id=2:ID
            elementRangeIncrement += beams[id-1].nElements
        end
        elementRange .+= elementRangeIncrement
        @pack! beam = elementRange

        # Add elements to model and update their global IDs on the assembly
        for element in elements
            push!(model.elements,element)
            element.globalID = element.localID + elementRangeIncrement
        end

        # Initialize elements' nodes' global ID array for simply connected beam
        for e in elementRange
            push!(elementNodes, [e, e+1])
        end

        # Update elements' nodes' global ID array in case the beam is not simply connected
        if ID > 1 && elementNodes[elementRange][1][2] > nodeCount+1
            diff = elementNodes[elementRange][1][2] - (nodeCount+1)
            elementNodes[elementRange] .= [x .- [diff,diff] for x in elementNodes[elementRange]]
        end
        for (connectedBeam,connectedNodeThis,connectedNodeOther) in zip( connectedBeams,connectedNodesThis,connectedNodesOther)
            # Local IDs of the elements of the current beam to which the connected node belongs
            currentBeamConnectedElementsLocalID = ifelse(connectedNodeThis==1,[1],ifelse(connectedNodeThis==nElements+1,[nElements],[connectedNodeThis-1;connectedNodeThis]))
            # Global IDs of the above
            currentBeamConnectedElementsGlobalID = elementRange[currentBeamConnectedElementsLocalID]
            # Sides of the connected node on the current beam's connected elements
            currentBeamNodeSides = ifelse(connectedNodeThis==1,[1],ifelse(connectedNodeThis==nElements+1,[2],[2;1]))
            # Local ID of the first element of the connected beam to which the connected node belongs
            connectedBeamConnectedElementLocalID = connectedNodeOther == 1 ? connectedNodeOther : connectedNodeOther-1
            # Side of the node on the connected beam
            connectedBeamNodeSide = connectedNodeOther == connectedBeamConnectedElementLocalID ? 1 : 2
            # Global ID of the node on the connected beam
            connectedBeamsConnectedNodeGlobalID = elementNodes[connectedBeam.elementRange[connectedBeamConnectedElementLocalID]][connectedBeamNodeSide]
            # In case there is more than one connected element in the current beam, increase the indices of the nodes up to the connected node
            if length(currentBeamConnectedElementsGlobalID) > 1
                elementNodes[elementRange[1:currentBeamConnectedElementsLocalID[1]]] .= [x .+ [1,1] for x in elementNodes[elementRange[1:currentBeamConnectedElementsLocalID[1]]]]
            end
            # Loop current beam's connected elements
            for (currentBeamConnectedElementGlobalID,currentBeamNodeSide) in zip(currentBeamConnectedElementsGlobalID,currentBeamNodeSides)
                # Set connected element's nodes' global ID
                elementNodes[currentBeamConnectedElementGlobalID][currentBeamNodeSide] = connectedBeamsConnectedNodeGlobalID
            end
        end

        # Update global ID of the elements' nodes
        for (e, element) in zip(elementRange, elements)
            element.nodesGlobalID = elementNodes[e]
        end

        # Update the range of global IDs of the nodes of the current beam in the assembly
        nodeRange = unique(vcat(elementNodes[elementRange]...))
        @pack! beam = nodeRange

        # Get current node count
        nodeCount = maximum(nodeRange)

        # Update nodal coordinates
        for element in elements
            # Reset element's variables in the assembly
            element.r_n1 = position_vector_from_curvature(beam.R0,beam.k,element.x1_n1)
            element.r_n2 = position_vector_from_curvature(beam.R0,beam.k,element.x1_n2)
            element.r = position_vector_from_curvature(beam.R0,beam.k,element.x1)
            # Add the initial position of the beam to element's nodal and midpoint coordinates
            element.r_n1 .+= beam.initialPosition 
            element.r_n2 .+= beam.initialPosition 
            element.r .+= beam.initialPosition
            # Add the initial position of the model to element's nodal and midpoint coordinates, if in the first beam
            if ID == 1
                element.r_n1 .+= model.initialPosition 
                element.r_n2 .+= model.initialPosition 
                element.r .+= model.initialPosition 
            end
            # Add coordinates of the beam's first node starting at the second beam
            if ID > 1
                firstNodeOfBeam = minimum(nodeRange)
                coordinatesOfFirstNode = r_n[firstNodeOfBeam]
                element.r_n1 .+= coordinatesOfFirstNode 
                element.r_n2 .+= coordinatesOfFirstNode
                element.r .+= coordinatesOfFirstNode
            end
            # Push current element's nodal coordinates
            push!(r_n,element.r_n1,element.r_n2)  
            # Take out repeated nodal coordinates
            round_off!.(r_n,1e-8)
            r_n = unique(r_n)
        end

        # Reset and update nodal coordinates on beam level
        beam.r_n = Vector{Vector{Float64}}()
        for element in elements
            push!(beam.r_n,element.r_n1,element.r_n2)
        end

        # Update array of assembly's sectional compliance matrices
        for element in elements
            push!(S,element.S)
        end

        # Update total number of elements of the model
        nElementsTotal += nElements

    end

    # Update total number of nodes
    nNodesTotal = maximum(vcat(elementNodes...))

    # Get force scaling
    forceScaling = force_scaling(S)

    @pack! model = beams,nElementsTotal,nNodesTotal,elementNodes,r_n,forceScaling

end


"""
inertia_properties!(model::Model)

Computes the inertia properties of the undeformed assembly

# Arguments
- model::Model
"""
function inertia_properties!(model::Model)

    @unpack elements = model

    # Mass and mass moment of inertia of the elements (including attached point inertias)
    elementsMass = [element.I[1,1]*element.Δℓ for element in elements] 
    elementsI = [abs.(element.Δℓ*element.R0*[element.I[4,4]; element.I[5,5]; element.I[6,6]]) for element in elements]

    # Mass times position of point inertias attached to the elements, resolved in basis A
    elementsMassη = [-element.Δℓ*element.R0*[element.I[5,3]; element.I[6,1]; element.I[4,2]] for element in elements]

    # Position of the elements
    elements_x1 = [element.r[1] for element in elements]
    elements_x2 = [element.r[2] for element in elements]
    elements_x3 = [element.r[3] for element in elements]

    # Mass and mass moments of inertia of the model
    mass = sum(elementsMass)
    I = sum(elementsI)

    # Mass times position of all attached point inertias
    massη = sum(elementsMassη)
    
    # Position of the center of mass, resolved in basis A
    centerOfMass = [(massη[1]+dot(elements_x1,elementsMass))/mass; (massη[2]+dot(elements_x2,elementsMass))/mass; (massη[3]+dot(elements_x3,elementsMass))/mass]

    @pack! model = mass,I,centerOfMass

end


"""
set_atmosphere!(model::Model)

Sets the atmosphere data

# Arguments
- model::Model
"""
function set_atmosphere!(model::Model)

    @unpack altitude,atmosphere = model

    # Altitude was input: set corresponding International Standard Atmosphere (ISA)   
    if !isnothing(altitude)
        atmosphere = standard_atmosphere(altitude)
    # Neither altidude nor atmosphere was input, but there are aerodynamic surfaces: set altitude as zero and corresponding ISA  
    elseif any(x -> x !== nothing, [beam.aeroSurface for beam in model.beams]) && isnothing(atmosphere)
        altitude = 0
        atmosphere = standard_atmosphere(altitude)
    end

    @pack! model = altitude,atmosphere

end


"""
update_number_gust_states!(model::Model)

Updates the number of gust states in every element with aerodynamic surface

# Arguments
- model::Model
"""
function update_number_gust_states!(model::Model) 

    @unpack gust,elements = model

    # Skip if there is no active gust
    if isnothing(gust)
        return 
    end

    # Loop over elements
    for element in elements
        @unpack aero = element
        # Skip if element does not have aerodynamic surface
        if isnothing(aero)
            continue
        end
        @unpack solver,gustLoadsSolver,nTotalAeroStates,pitchPlungeStatesRange,airfoil,c,normSparPos = aero
        # Update gust states range
        nGustStates = gustLoadsSolver.nStates
        gustStatesRange = nTotalAeroStates+1:nTotalAeroStates+nGustStates
        # Update number of total aerodynamic states
        nTotalAeroStates += nGustStates
        # Resize arrays
        A = zeros(nTotalAeroStates,nTotalAeroStates)
        B = zeros(nTotalAeroStates)
        f1χ_χ = zeros(3,nTotalAeroStates)
        f2χ_χ = zeros(3,nTotalAeroStates)
        m1χ_χ = zeros(3,nTotalAeroStates)
        m2χ_χ = zeros(3,nTotalAeroStates)
        F_χ_V = zeros(nTotalAeroStates,3)
        F_χ_Ω = zeros(nTotalAeroStates,3)
        F_χ_χ = initial_F_χ_χ(solver,nTotalAeroStates)
        F_χ_δ = zeros(nTotalAeroStates)
        F_χ_Vdot = initial_F_χ_Vdot(solver,nTotalAeroStates,pitchPlungeStatesRange,airfoil.attachedFlowParameters.cnα)
        F_χ_Ωdot = initial_F_χ_Ωdot(solver,nTotalAeroStates,pitchPlungeStatesRange,c,normSparPos,airfoil.attachedFlowParameters.cnα)
        F_χ_χdot = Matrix(1.0*LinearAlgebra.I,nTotalAeroStates,nTotalAeroStates)
        # Pack data
        @pack! nTotalAeroStates,nGustStates,gustStatesRange,A,B,f1χ_χ,f2χ_χ,m1χ_χ,m2χ_χ,F_χ_V,F_χ_Ω,F_χ_χ,F_χ_δ,F_χ_Vdot,F_χ_Ωdot,F_χ_χdot = element.aero
    end

end


"""
update_loads_trim_links_global_ids!(model::Model)

Updates the nodes' global IDs of the loads trim links

# Arguments
- model::Model
"""
function update_loads_trim_links_global_ids!(model::Model)

    @unpack trimLoadsLinks = model

    # Loop trim load links
    for loadsLink in trimLoadsLinks
        # Unpack trim loads link data
        @unpack masterBeam,masterNodeLocalID,slaveBeams,slaveNodesLocalIDs = loadsLink
        # Update master node global IDs 
        masterNodeGlobalID = model.beams[masterBeam.ID].nodeRange[masterNodeLocalID]
        # Update slave node global IDs 
        slaveNodesGlobalIDs = Vector{Int64}()
        for (slaveBeam,slaveNodeLocalID) in zip(slaveBeams,slaveNodesLocalIDs)
            slaveNodeGlobalID = model.beams[slaveBeam.ID].nodeRange[slaveNodeLocalID]
            push!(slaveNodesGlobalIDs,slaveNodeGlobalID)
        end
        # Pack trim loads link data
        @pack! loadsLink = masterNodeGlobalID,slaveNodesGlobalIDs
    end

    @pack! model = trimLoadsLinks
end


"""
update_spring_nodes_global_ids!(model::Model)

Updates the nodes' global IDs of the attached springs

# Arguments
- model::Model
"""
function update_spring_nodes_global_ids!(model::Model)

    @unpack beams = model

    # Loop beams 
    for beam in beams
        @unpack springs,elements = beam
        # Loop springs
        for spring in springs
            @unpack elementID,localNode = spring
            # Global ID of the spring node
            nodeGlobalID = elements[elementID].nodesGlobalID[localNode]
            # Rotation tensor from basis A to basis b
            R0 = localNode == 1 ? elements[elementID].R0_n1 : elements[elementID].R0_n2
            # Pack
            @pack! spring = nodeGlobalID,R0
        end
    end

end


"""
update_linked_flap_deflections!(model::Model)

Updates the linked flap deflections for the slave beams

# Arguments
- model::Model
"""
function update_linked_flap_deflections!(model::Model)

    @unpack flapLinks = model

    # Loop flap links
    for flapLink in flapLinks
        @unpack masterBeam,slaveBeams,δMultipliers = flapLink
        # Loop slave beams
        for (i,slaveBeam) in enumerate(slaveBeams)
            # Flapped elements of the beam
            flappedElements = slaveBeam.elements[[element.aero.flapped for element in slaveBeam.elements]]
            # Update TF for trim flap deflection
            [setfield!(element.aero, :δIsTrimVariable, masterBeam.elements[1].aero.δIsTrimVariable) for element in flappedElements]
            # Get flap deflection and rates of elements equal to the values of the master beam times the slave's deflection multiplier
            δSlave = t -> masterBeam.elements[1].aero.δ(t)*δMultipliers[i]
            δdotSlave = t -> masterBeam.elements[1].aero.δdot(t)*δMultipliers[i]
            δddotSlave = t -> masterBeam.elements[1].aero.δddot(t)*δMultipliers[i]
            # Set flap deflection, its rates and deflection multipliers of the slave elements 
            [setfield!(element.aero, :δ, δSlave) for element in flappedElements]
            [setfield!(element.aero, :δdot, δdotSlave) for element in flappedElements]
            [setfield!(element.aero, :δddot, δddotSlave) for element in flappedElements]
            [setfield!(element.aero, :δNow, element.aero.δ(0)) for element in flappedElements]
            [setfield!(element.aero, :δdotNow, element.aero.δdot(0)) for element in flappedElements]
            [setfield!(element.aero, :δddotNow, element.aero.δddot(0)) for element in flappedElements]
            [setfield!(element.aero, :δMultiplier, δMultipliers[i]) for element in flappedElements]
        end
        # Loop beams
        for beam in vcat(masterBeam,slaveBeams...)
            # Update flag
            beam.aeroSurface.hasIndependentFlap = false
        end
    end

end


"""
initialize_basis_A_rotation!(model::Model)

Initializes the rotation tensor from basis I (fixed, inertial) to basis A, and its transpose

# Arguments
- model::Model
"""
function initialize_basis_A_rotation!(model::Model)

    @unpack p_A0 = model

    R_A,_ = rotation_tensor_WM(p_A0)
    R_AT = Matrix(R_A')

    @pack! model = R_A,R_AT

end


"""
update_initial_conditions!(model::Model)

Updates the initial condition states on all elements of the assembly

# Arguments
- model::Model
"""
function update_initial_conditions!(model::Model)

    @unpack beams,specialNodes = model

    # Loop over beams
    for beam in beams

        @unpack u0_of_x1,p0_of_x1,udot0_of_x1,pdot0_of_x1,uprime0_of_x1,pprime0_of_x1 = beam

        # Loop over elements
        for element in beam.elements

            # Unpack
            @unpack states,statesRates,compStates,nodalStates,Δℓ,x1,x1_n1,x1_n2,R0,R0T,R0T_n1,R0T_n2,S,I,k = element

            # b basis' generalized velocities at element's midpoint, resolved in basis A
            v,ω = element_velocities_basis_b!(model,element)

            # b basis' generalized accelerations at element's midpoint, resolved in basis A
            vdot,ωdot = element_accelerations_basis_b!(model,element)

            # Generalized displacements, resolved in basis A
            u = R0*u0_of_x1(x1)
            p = R0*p0_of_x1(x1)

            # Derivatives of generalized displacements, resolved in basis A
            u_prime = R0*uprime0_of_x1(x1)
            p_prime = R0*pprime0_of_x1(x1)

            # Rotation and tangent tensors and derivatives, resolved in basis A
            R, = rotation_tensor_WM(p)
            HT = tangent_operator_transpose_WM(p)
            RR0 = R*R0
            RT = Matrix(R')
            RR0T = Matrix(RR0')

            # Strains, resolved in basis b
            γ = RR0T*(R0*a1+u_prime) - a1
            κ = R0T*HT*p_prime
            round_off!(γ)

            # Sectional forces, resolved in basis B
            sectionalForces = inv(S)*[γ; κ]
            F = sectionalForces[1:3]
            M = sectionalForces[4:6]

            # Generalized displacements' rates, resolved in basis A
            udot = R0*udot0_of_x1(x1)
            pdot = R0*pdot0_of_x1(x1)

            # Generalized sectional velocities, resolved in basis B
            V = RR0T*(v+cross(ω,u)+udot)
            Ω = R0T*(HT*pdot+RT*ω)
            
            # Sectional linear and angular momenta, resolved in basis B 
            P = I[1:3,:]*[V; Ω]
            H = I[4:6,:]*[V; Ω]

            # Midpoint deformed curvature vector, resolved in basis B
            K = k+κ

            # Midpoint generalized forces' derivatives with respect to x1 (neglecting the unknown momenta rates)
            F_prime = cross(Ω,P)-cross(K,F)
            M_prime = cross(Ω,H)+cross(V,P)-cross(a1+γ,F)-cross(K,M)

            # Nodal displacements and rotations, resolved in basis A
            u_n1 = u0_of_x1(x1_n1)
            u_n2 = u0_of_x1(x1_n2)
            p_n1 = p0_of_x1(x1_n1)
            p_n2 = p0_of_x1(x1_n2)

            # Nodal displacements and rotation parameters' vectors, resolved in basis b
            u_n1_b = R0T_n1*u_n1
            u_n2_b = R0T_n2*u_n2
            p_n1_b = R0T_n1*p_n1
            p_n2_b = R0T_n2*p_n2

            # Nodal rotation angles
            θ_n1 = rotation_angle(p_n1)
            θ_n2 = rotation_angle(p_n2)

            # Nodal forces and moments (neglecting distributed loads' resultants)
            F_n1 = F - (Δℓ/2*F_prime)
            F_n2 = F + (Δℓ/2*F_prime)
            M_n1 = M - (Δℓ/2*M_prime)
            M_n2 = M + (Δℓ/2*M_prime)

            # Pack variables onto element
            @pack! states = u,p,F,M,V,Ω
            @pack! compStates = γ,κ,P,H
            @pack! statesRates = udot,pdot
            @pack! nodalStates = u_n1,u_n2,p_n1,p_n2,u_n1_b,u_n2_b,p_n1_b,p_n2_b,θ_n1,θ_n2,F_n1,F_n2,M_n1,M_n2
            @pack! element = states,statesRates,compStates,nodalStates,v,ω,vdot,ωdot,R,HT,RR0,RR0T

        end
    end

    # Loop over special nodes
    for specialNode in specialNodes
        @unpack connectedElements,ζonElements = specialNode
        u = ζonElements[1] == -1 ? connectedElements[1].nodalStates.u_n1 : connectedElements[1].nodalStates.u_n2
        p = ζonElements[1] == -1 ? connectedElements[1].nodalStates.p_n1 : connectedElements[1].nodalStates.p_n2
        F = ζonElements[1] == -1 ? connectedElements[1].nodalStates.F_n1 : connectedElements[1].nodalStates.F_n2
        M = ζonElements[1] == -1 ? connectedElements[1].nodalStates.M_n1 : connectedElements[1].nodalStates.M_n2
        @pack! specialNode = u,p,F,M
    end

end


"""
set_BCs!(model::Model)

Sets the BCs onto the model

# Arguments
- model::Model
"""
function set_BCs!(model::Model,BCs::Vector{BC})

    # Set BCs on model
    model.BCs = BCs

    # Loop BCs
    for BC in BCs
        # Set global node ID on BC
        BC.globalNodeID = BC.beam.nodeRange[BC.node]
        # Add that BC's node to the BC'ed nodes list
        push!(model.BCedNodes,BC.globalNodeID)
    end

    # Remove duplicated nodes from list
    unique!(model.BCedNodes)
end


"""
get_special_nodes!(model::Model)

Gets the special nodes in the system of equations: connection, boundary, and BC'ed nodes

# Arguments
- model::Model
"""
function get_special_nodes!(model::Model)

    @unpack elements,nNodesTotal,elementNodes,BCs,BCedNodes,beams = model

    # Initialize array of special nodes
    specialNodes = Vector{SpecialNode}()

    # List of all nodes
    nodesList = 1:nNodesTotal

    # Initialize TF for nodes being special
    special = falses(nNodesTotal)

    # BC'ed nodes (any essential or non-zero natural BC) are special 
    special[BCedNodes] .= true
    
    # Nodes with attached springs and their properties
    springNodes = [spring.nodeGlobalID for beam in beams for spring in beam.springs]
    springNodesku = [spring.ku for beam in beams for spring in beam.springs]
    springNodeskp = [spring.kp for beam in beams for spring in beam.springs]
    springNodesR0 = [spring.R0 for beam in beams for spring in beam.springs]

    # Nodes with attached springs are special
    special[springNodes] .= true

    # Loop over nodes
    for node in nodesList  

        # Elements connected to this node and their IDs
        connectedElementsGlobalIDs = findall(x -> node in x, elementNodes)
        connectedElements = elements[connectedElementsGlobalIDs]

        # Number of elements that are connected to the node
        nConnectedElements = length(connectedElements)

        # Find in which side and ζ position of the elements that node is 
        sideOnElements = Vector{Int64}(undef, nConnectedElements)
        ζonElements = Vector{Int64}(undef, nConnectedElements)
        for (i, e) in enumerate(connectedElementsGlobalIDs)
            sideOnElements[i] = elementNodes[e][1] == node ? 1 : 2
            ζonElements[i] = elementNodes[e][1] == node ? -1 : 1
        end

        # Get the local and global IDs of the node
        localID = connectedElements[1].nodesLocalID[sideOnElements[1]]
        globalID = connectedElements[1].nodesGlobalID[sideOnElements[1]]

        # If this is either a boundary (connected to only 1 element) or a connection node (connected to 3 or more elements), set as special 
        if nConnectedElements != 2                         
            special[node] = true                
        end

        # Add stiffness constants resolved in basis A for current beam, if applicable
        R0_ku = zeros(3)
        R0_kp = zeros(3)
        ids = findall(x -> x==node, springNodes)
        for id in ids
            R0_ku += springNodesR0[id] * springNodesku[id]
            R0_kp += springNodesR0[id] * springNodeskp[id]
        end

        # Add special node to array
        if special[node]
            # Find if that node is BC'ed and get the BCs
            if node in BCedNodes
                nodesBCs = Vector{BC}()
                for BC in BCs
                    if node == BC.globalNodeID
                        push!(nodesBCs,BC)
                    end
                end
                push!(specialNodes,SpecialNode(localID,globalID,connectedElementsGlobalIDs,connectedElements,ζonElements,R0_ku,R0_kp,nodesBCs))
            else
                push!(specialNodes,SpecialNode(localID,globalID,connectedElementsGlobalIDs,connectedElements,ζonElements,R0_ku,R0_kp))
            end
        end
    end

    # Global IDs of the special nodes
    specialNodesGlobalIDs = nodesList[special]

    @pack! model = specialNodes,specialNodesGlobalIDs
end


"""
get_system_indices!(model::Model)

Gets the indices (for equations and DOFs) of the system of equations

# Arguments
- model::Model
"""
function get_system_indices!(model::Model)

    @unpack beams,elements,nNodesTotal,specialNodes,specialNodesGlobalIDs,BCedNodes,nTrimVariables,trimLoadsLinks,flapLinks = model
    
    ## Initialize flags
    #---------------------------------------------------------------------------
    # Flag indicating whether indices of a node's equations/states indices have been assigned
    assigned = falses(nNodesTotal)
    # Array with element on which node's equations have been assigned (first column contains the elements, and second column which local node it is in that element)
    e_assigned::Vector{Vector{Int64}} = Vector{Vector{Int64}}()
    # Indices for next equations/states
    i_equations = 1                   
    i_states = 1
    
    ## Set indices of elemental and standard nodes' structural equations/states
    #---------------------------------------------------------------------------
    # Loop over elements
    for (e, element) in enumerate(elements)
        # Unpack element data
        @unpack nodesGlobalID,eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,eqs_FF1_sep,eqs_FM1_sep,eqs_FF2_sep,eqs_FM2_sep,isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set = element      
        # Get element's first node's global ID
        n = nodesGlobalID[1]
        # Check if the node's equations/states have been assigned indices
        if !assigned[n]
            # Update assignment flag and pointer to element
            assigned[n] = true
            push!(e_assigned,[e, 1])
            # Indices for equilibrium and compatibility equations of element's first node
            eqs_Fu1 = i_equations+0:i_equations+2
            eqs_Fp1 = i_equations+3:i_equations+5
            eqs_FF1 = i_equations+6:i_equations+8
            eqs_FM1 = i_equations+9:i_equations+11
            # Update equations' index
            i_equations += 12
            # Check whether this is a special node - has nodal states
            if n in specialNodesGlobalIDs
                # Get corresponding special node
                s = findfirst(x -> n in x, specialNodesGlobalIDs)
                specialNode = specialNodes[s]
                # Indices for nodal states
                specialNode.DOF_uF = i_states+0:i_states+2
                specialNode.DOF_pM = i_states+3:i_states+5
                # Update states' index
                i_states += 6
                # Set flag
                isSpecialNode1 = true
            end
        else   
            # Update flag for equations having already been set
            eqsNode1Set = true
            # Indices for equilibrium and compatibility equations of element's first node
            if e_assigned[n][2] == 1
                eqs_Fu1 = elements[e_assigned[n][1]].eqs_Fu1
                eqs_Fp1 = elements[e_assigned[n][1]].eqs_Fp1
                eqs_FF1 = elements[e_assigned[n][1]].eqs_FF1
                eqs_FM1 = elements[e_assigned[n][1]].eqs_FM1
            elseif e_assigned[n][2] == 2
                eqs_Fu1 = elements[e_assigned[n][1]].eqs_Fu2
                eqs_Fp1 = elements[e_assigned[n][1]].eqs_Fp2
                eqs_FF1 = elements[e_assigned[n][1]].eqs_FF2
                eqs_FM1 = elements[e_assigned[n][1]].eqs_FM2
            end
            # Check whether this is a special node - has separate compatibility equations
            if n in specialNodesGlobalIDs       
                # Indices for separate compatibility equations  
                eqs_FF1_sep = i_equations+0:i_equations+2
                eqs_FM1_sep = i_equations+3:i_equations+5
                # Update equations' index 
                i_equations += 6
                # Set flag
                isSpecialNode1 = true
            end
        end

        # Indices for structural elemental states (u, p, F, M, V, Ω)
        DOF_u = i_states+0:i_states+2
        DOF_p = i_states+3:i_states+5
        DOF_F = i_states+6:i_states+8
        DOF_M = i_states+9:i_states+11
        DOF_V = i_states+12:i_states+14
        DOF_Ω = i_states+15:i_states+17
        # Update states' index
        i_states += 18
        # Indices for element's generalized velocity-displacement equations
        eqs_FV = i_equations+0:i_equations+2
        eqs_FΩ = i_equations+3:i_equations+5
        # Update equations' index
        i_equations += 6
            
        # Get element's second node's global ID
        n = nodesGlobalID[2]
        # Check if the node's equations/states have been assigned indices
        if !assigned[n]
            # Update assignment flag and pointer to element
            assigned[n] = true
            push!(e_assigned,[e, 2])
            # Indices for equilibrium and compatibility equations of element's second node
            eqs_Fu2 = i_equations+0:i_equations+2
            eqs_Fp2 = i_equations+3:i_equations+5
            eqs_FF2 = i_equations+6:i_equations+8
            eqs_FM2 = i_equations+9:i_equations+11
            # Update equations' index
            i_equations += 12
            # Check whether this is a special node - has nodal states
            if n in specialNodesGlobalIDs
                # Get corresponding special node
                s = findfirst(x -> n in x, specialNodesGlobalIDs)
                specialNode = specialNodes[s]
                # Indices for nodal states 
                specialNode.DOF_uF = i_states+0:i_states+2
                specialNode.DOF_pM = i_states+3:i_states+5
                # Update states' index
                i_states += 6
                # Set flag
                isSpecialNode2 = true
            end
        else   
            # Update flag for equations having already been set
            eqsNode2Set = true
            # Indices for equilibrium and compatibility equations of element's second node
            if e_assigned[n][2] == 1
                eqs_Fu2 = elements[e_assigned[n][1]].eqs_Fu1
                eqs_Fp2 = elements[e_assigned[n][1]].eqs_Fp1
                eqs_FF2 = elements[e_assigned[n][1]].eqs_FF1
                eqs_FM2 = elements[e_assigned[n][1]].eqs_FM1
            elseif e_assigned[n][2] == 2
                eqs_Fu2 = elements[e_assigned[n][1]].eqs_Fu2
                eqs_Fp2 = elements[e_assigned[n][1]].eqs_Fp2
                eqs_FF2 = elements[e_assigned[n][1]].eqs_FF2
                eqs_FM2 = elements[e_assigned[n][1]].eqs_FM2
            end
            # Check whether this is a special node - has separate compatibility equations
            if n in specialNodesGlobalIDs   
                # Indices for separate compatibility equations
                eqs_FF2_sep = i_equations+0:i_equations+2
                eqs_FM2_sep = i_equations+3:i_equations+5
                # Update equations' index 
                i_equations += 6
                # Set flag
                isSpecialNode2 = true
            end
        end
        # Pack element data
        @pack! element = eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FV,eqs_FΩ,DOF_u,DOF_p,DOF_F,DOF_M,DOF_V,DOF_Ω,eqs_FF1_sep,eqs_FM1_sep,eqs_FF2_sep,eqs_FM2_sep,isSpecialNode1,isSpecialNode2,eqsNode1Set,eqsNode2Set
    end

    ## Set indices of elemental aerodynamic equations/states
    #---------------------------------------------------------------------------
    # Loop over elements
    for element in elements
        # Skip if there are no aerodynamic loads
        if isnothing(element.aero)
            continue
        end
        # Unpack element data
        @unpack nTotalAeroStates = element.aero
        # Set states and equations
        DOF_χ = i_states+0:i_states+nTotalAeroStates-1
        eqs_Fχ = i_equations+0:i_equations+nTotalAeroStates-1
        # Update indices
        i_equations += nTotalAeroStates
        i_states += nTotalAeroStates
        # Pack element data
        @pack! element = DOF_χ,eqs_Fχ
    end

    # Number of states/equations (system's order)
    systemOrder = i_equations - 1

    ## Set indices for special nodes' equations
    # --------------------------------------------------------------------------
    # Loop over special nodes
    for specialNode in specialNodes
        @unpack connectedElements,ζonElements,eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,eqs_FF_sep,eqs_FM_sep = specialNode
        # Loop over connected elements
        for (e, element) in enumerate(connectedElements)
            @unpack eqs_Fu1,eqs_Fu2,eqs_Fp1,eqs_Fp2,eqs_FF1,eqs_FF2,eqs_FM1,eqs_FM2,eqs_FF1_sep,eqs_FF2_sep,eqs_FM1_sep,eqs_FM2_sep = element
            # Check in which side of the element the node is and set indices for nodal equations accordingly
            if ζonElements[e] == -1
                eqs_Fu[e] = eqs_Fu1
                eqs_Fp[e] = eqs_Fp1
                eqs_FF[e] = eqs_FF1
                eqs_FM[e] = eqs_FM1
                eqs_FF_sep[e] = eqs_FF1_sep
                eqs_FM_sep[e] = eqs_FM1_sep
            elseif ζonElements[e] == 1
                eqs_Fu[e] = eqs_Fu2
                eqs_Fp[e] = eqs_Fp2
                eqs_FF[e] = eqs_FF2
                eqs_FM[e] = eqs_FM2
                eqs_FF_sep[e] = eqs_FF2_sep
                eqs_FM_sep[e] = eqs_FM2_sep
            end
        end
        # Pack nodal data
        @pack! specialNode = eqs_Fu,eqs_Fp,eqs_FF,eqs_FM,eqs_FF_sep,eqs_FM_sep
    end

    # Initialize number of independent trim loads and trim flap deflections
    nTrimLoads,nTrimδ = 0,0

    # Set nodal trim loads indices for master nodes
    for (specialNode,specialNodeGlobalID) in zip(specialNodes,specialNodesGlobalIDs)
        # Initialize trim loads state indices
        DOF_trimLoads = zeros(Int64,6)
        @pack! specialNode = DOF_trimLoads
        # Skip non-BC'ed nodes
        if !(specialNodeGlobalID in BCedNodes)
            continue
        end
        # Treat cases with and without trim loads links separately
        if isempty(trimLoadsLinks)
            # Loop over nodes' BCs
            for BC in specialNode.BCs
                @unpack isTrim = BC
                # Skip if there are no trim loads
                if !any(isTrim)
                    continue
                end
                # Loop load indices
                for i=1:6
                    # Update nodal trim loads count and set states indices 
                    if isTrim[i] 
                        nTrimLoads += 1
                        DOF_trimLoads[i] = systemOrder+nTrimLoads
                    end
                end
            end
        else
            # Loop trim loads links
            for link in trimLoadsLinks
                @unpack masterNodeGlobalID,masterBC = link
                @unpack isTrim = masterBC
                # Skip if not a master node 
                if specialNodeGlobalID != masterNodeGlobalID
                    continue
                end
                # Loop load indices
                for i=1:6
                    # Update nodal trim loads count and set states indices 
                    if isTrim[i] 
                        nTrimLoads += 1
                        DOF_trimLoads[i] = systemOrder+nTrimLoads
                    end
                end
            end
        end
        # Pack nodal data
        @pack! specialNode = DOF_trimLoads 
    end

    # Set nodal trim loads indices for slave nodes
    for (specialNode,specialNodeGlobalID) in zip(specialNodes,specialNodesGlobalIDs)
        # Skip non-BC'ed nodes
        if !(specialNodeGlobalID in BCedNodes)
            continue
        end
        @unpack DOF_trimLoads = specialNode
        # Skip if DOF_trimLoads has been set (this is a master node)
        if any(!iszero(DOF_trimLoads))
            continue
        end
        # This is a slave node: loop trim loads links
        for link in trimLoadsLinks
            @unpack masterNodeGlobalID,slaveNodesGlobalIDs,slaveBCs = link
            # Skip if the node is not in the current link
            if !(specialNodeGlobalID in slaveNodesGlobalIDs)
                continue
            end
            # Get the master node's trim loads DoF
            masterSpecialNodeID = findfirst(x -> x.globalID == masterNodeGlobalID, specialNodes)
            DOF_trimLoadsMaster = specialNodes[masterSpecialNodeID].DOF_trimLoads
            # Loop slave BCs 
            for slaveBC in slaveBCs
                @unpack isTrim = slaveBC
                # Loop load indices
                for i=1:6
                    # Update nodal trim state indices as those of the master node
                    if isTrim[i] 
                        DOF_trimLoads[i] = DOF_trimLoadsMaster[i]
                    end
                end
            end
        end
        # Pack nodal data
        @pack! specialNode = DOF_trimLoads 
    end

    # Set elemental trim flap deflections indices for master and independent beams
    for beam in beams
        # Skip if element does not have a trim flap deflection
        if isnothing(beam.aeroSurface) || !beam.aeroSurface.δIsTrimVariable
            continue
        end
        # Treat cases with and without flap links separately
        if isempty(flapLinks) || beam.aeroSurface.hasIndependentFlap
            # Increment trim flap deflections count
            nTrimδ += 1
            # Flapped elements of the beam
            flappedElements = beam.elements[[element.aero.flapped for element in beam.elements]]
            # Set the same DOF_δ for all flapped elements of the beam
            [setfield!(element, :DOF_δ, systemOrder+nTrimLoads+nTrimδ) for element in flappedElements]
        else
            # Loop flap links
            for flapLink in flapLinks
                # Skip if not a master beam
                if beam !== flapLink.masterBeam
                    continue
                end
                # Increment trim flap deflections count
                nTrimδ += 1
                # Flapped elements of the beam
                flappedElements = beam.elements[[element.aero.flapped for element in beam.elements]]
                # Set the same DOF_δ for all flapped elements of the beam
                [setfield!(element, :DOF_δ, systemOrder+nTrimLoads+nTrimδ) for element in flappedElements]
            end
        end
    end

    # Set elemental trim flap deflections indices for slave beams
    for beam in beams
        # Skip if element does not have a trim flap deflection
        if isnothing(beam.aeroSurface) || !beam.aeroSurface.δIsTrimVariable
            continue
        end
        # Loop flap links
        for flapLink in flapLinks
            # Skip if not a slave beam
            if !(beam in flapLink.slaveBeams)
                continue
            end
            # Flapped elements of the beam
            flappedElements = beam.elements[[element.aero.flapped for element in beam.elements]]
            # Set DOF_δ for flapped element of the beam as that of the master beam
            [setfield!(element, :DOF_δ, flapLink.masterBeam.elements[1].DOF_δ) for element in flappedElements]
        end
    end

    # Update number of trim variables
    nTrimVariables = nTrimLoads+nTrimδ
    
    @pack! model = systemOrder,elements,specialNodes,nTrimVariables

end


"""
set_motion_basis_A!(; model::Model,u_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,v_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,ω_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,vdot_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,ωdot_A::Union{Vector{<:Number},<:Function,Nothing}=nothing)

Sets the motion of basis A into the model

# Arguments
- model::Model
- u_A::Union{Vector{<:Number},<:Function,Nothing}=nothing
- v_A::Union{Vector{<:Number},<:Function,Nothing}=nothing
- ω_A::Union{Vector{<:Number},<:Function,Nothing}=nothing
- vdot_A::Union{Vector{<:Number},<:Function,Nothing}=nothing
- ωdot_A::Union{Vector{<:Number},<:Function,Nothing}=nothing
"""
function set_motion_basis_A!(; model::Model,u_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,v_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,ω_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,vdot_A::Union{Vector{<:Number},<:Function,Nothing}=nothing,ωdot_A::Union{Vector{<:Number},<:Function,Nothing}=nothing)

    # Reset values that were not input back to nothing on model
    @pack! model = u_A,v_A,ω_A,vdot_A,ωdot_A 

    # Update model
    update_model!(model)
    
end
export set_motion_basis_A!


"""
plot_undeformed_assembly(model::Model)

Plots the nodal coordinates of the assembly of beams

# Arguments
- model::Model
"""
function plot_undeformed_assembly(model::Model,view=(45,45))

    # Initialize backend
    gr()

    # Initialize plot
    plt = plot(;xlabel="\$x_1\$",ylabel="\$x_2\$",zlabel="\$x_3\$",title="Undeformed assembly",camera=view,aspect_ratio=:equal,grid=:true)

    # Initialize plot limits
    x1min,x1max,x2min,x2max,x3min,x3max=0,0,0,0,0,0

    # Loop over beams
    for beam in model.beams
        
        # Nodal coordinates for the current beam
        @unpack r_n = beam

        # Extract x1, x2, and x3 coordinates from the nodal coordinates
        x1 = [point[1] for point in r_n]
        x2 = [point[2] for point in r_n]
        x3 = [point[3] for point in r_n]

        # Update plot limits
        x1ext, x2ext, x3ext = extrema(x1), extrema(x2), extrema(x3)
        x1min,x1max = min(x1min,x1ext[1]),max(x1max,x1ext[2])
        x2min,x2max = min(x2min,x2ext[1]),max(x2max,x2ext[2])
        x3min,x3max = min(x3min,x3ext[1]),max(x3max,x3ext[2])
        # plot!(xlims=(x1min,x1max), ylims=(x2min,x2max), zlims=(x3min,x3max))
        
        # Plot nodes
        scatter!(x1, x2, x3, c=:black, ms=3, label=false)
        
        # Plot lines 
        for i in 1:length(r_n)-1
            plot!([r_n[i][1], r_n[i+1][1]], [r_n[i][2], r_n[i+1][2]], [r_n[i][3], r_n[i+1][3]], c=:black, lw=2, label=false)
        end

    end

    display(plt)

    return plt

end
export plot_undeformed_assembly
