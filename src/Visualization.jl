using Plots, ColorSchemes

"""
    plot_undeformed_assembly(model::Model; kwars...)

Plots the undeformed assembly of beams in the model

# Arguments
- `model::Model`

# Keyword arguments
- `view::Tuple{<:Number,<:Number}` = view angles
- `equalAspectRatio::Bool` = flag to set equal aspect ratio plot
- `plotNodes::Bool` = flag to plot nodes
- `nodesColor=:black` = color of nodes
- `linesColor=:black` = color of lines (beams)
"""
function plot_undeformed_assembly(model::Model; view::Tuple{Int64,Int64}=(45,45),equalAspectRatio::Bool=true,plotNodes::Bool=true,nodesColor=:black,linesColor=:black)

    # Set backend
    gr()

    # Initialize plot
    plt = plot(xlabel="\$x_1\$",ylabel="\$x_2\$",zlabel="\$x_3\$",title="Undeformed assembly",camera=view,aspect_ratio=:equal,grid=:true)

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
        
        # Plot nodes
        if plotNodes
            scatter!(x1, x2, x3, c=nodesColor, ms=3, label=false)
        end
        
        # Plot lines 
        for i in 1:length(r_n)-1
            plot!([r_n[i][1], r_n[i+1][1]], [r_n[i][2], r_n[i+1][2]], [r_n[i][3], r_n[i+1][3]], c=linesColor, lw=2, label=false)
        end

        # Update plot limits, if applicable
        if equalAspectRatio
            x1ext, x2ext, x3ext = extrema(x1), extrema(x2), extrema(x3)
            x1min,x1max = min(x1min,x1ext[1]),max(x1max,x1ext[2])
            x2min,x2max = min(x2min,x2ext[1]),max(x2max,x2ext[2])
            x3min,x3max = min(x3min,x3ext[1]),max(x3max,x3ext[2])
            lowerLim = min(x1min,x2min,x3min)
            upperLim = max(x1max,x2max,x3max)
            plot!(xlims=(lowerLim,upperLim), ylims=(lowerLim,upperLim), zlims=(lowerLim,upperLim))
        end

    end

    return plt
end
export plot_undeformed_assembly


"""
    plot_steady_deformation(problem::Problem; kwargs...)

Plots the initial and final deformed states for the model in the given problem

# Arguments
- `problem::Problem` = problem

# Keyword arguments
- `plotBCs::Bool` = flag to plot BCs
- `plotDistLoads::Bool` = flag to plot distributed loads (includes gravitational and aerodynamic loads)
- `view::Union{Nothing,Tuple{Real,Real}}` = view angles
- `scale::Number` = displacements and rotations scale
- `lw::Number` = linewidth
- `colorUndef` = color for undeformed assembly
- `colorDef` = color for deformed assembly
- `grid::Bool` = flag for grid
- `legendPos` = legend position
- `tolPlane::Number` = displacement tolerance to plot as plane
- `plotAeroSurf` = flag to plot aerodynamic surfaces
- `surfα::Float64` = transparency factor of aerodynamic surfaces 
- `ΔuDef::Vector{<:Number}` = displacement vector for first node of deformed assembly relative to the undeformed one
- `plotLimits::Union{Nothing,Vector{Tuple{T1,T2}}}` = plot axis limits    
- `showScale::Bool` = flag to show scale on plot
- `scalePos::Vector{<:Number}` = position of scale on plot
- `save::Bool` = flag to save the figure
- `savePath::String` = relative path on which to save the figure
"""
function plot_steady_deformation(problem::Problem; plotBCs::Bool=true,plotUndeformed::Bool=true,plotDistLoads::Bool=true,view::Union{Nothing,Tuple{Real,Real}}=nothing,scale::Number=1,lw::Number=1,colorUndef=:black,colorDef=:blue,grid::Bool=true,legendPos=:best,tolPlane::Number=1e-8,plotAeroSurf::Bool=true,surfα::Float64=0.5,ΔuDef::Vector{<:Number}=zeros(3),plotLimits::Union{Nothing,Vector{Tuple{T1,T2}}}=nothing,showScale::Bool=false,scalePos::Vector{<:Real}=[0,0],save::Bool=false,savePath::String="/test/outputs/figures/fig.pdf") where {T1<:Real,T2<:Real}

    # Validate
    @assert typeof(problem) in [SteadyProblem,TrimProblem,EigenProblem]
    @assert scale > 0
    @assert lw > 0
    @assert tolPlane > 0
    @assert 0 < surfα <= 1
    @assert length(ΔuDef) == 3

    # Set backend
    pyplot()

    # Unpack
    @unpack elements,nElementsTotal,units = problem.model

    # Initialize plot
    plt = plot(grid=grid)
    if plotUndeformed
        plot!([NaN],[NaN], c=colorUndef, lw=lw, label="Undeformed")
    end
    plot!([NaN],[NaN], c=colorDef, lw=lw, label="Deformed")
    plot!(legend=legendPos)

    # Initialize arrays
    x1Undef = Array{Vector{Float64}}(undef,nElementsTotal)
    x2Undef = Array{Vector{Float64}}(undef,nElementsTotal)
    x3Undef = Array{Vector{Float64}}(undef,nElementsTotal)
    x1Def = Array{Vector{Float64}}(undef,nElementsTotal)
    x2Def = Array{Vector{Float64}}(undef,nElementsTotal)
    x3Def = Array{Vector{Float64}}(undef,nElementsTotal)
    undefAirfoilCoords_n1 = Array{Union{Nothing,Matrix{Float64}}}(undef,nElementsTotal)
    undefAirfoilCoords_n2 = Array{Union{Nothing,Matrix{Float64}}}(undef,nElementsTotal)
    defAirfoilCoords_n1 = Array{Union{Nothing,Matrix{Float64}}}(undef,nElementsTotal)
    defAirfoilCoords_n2 = Array{Union{Nothing,Matrix{Float64}}}(undef,nElementsTotal)

    # Loop over elements
    for (e,element) in enumerate(elements)     
        @unpack r_n1,r_n2,nodalStates,aero = element
        @unpack u_n1,u_n2,p_n1,p_n2 = nodalStates
        # Set undeformed coordinates
        x1Undef[e] = [r_n1[1]; r_n2[1]]
        x2Undef[e] = [r_n1[2]; r_n2[2]]
        x3Undef[e] = [r_n1[3]; r_n2[3]]
        # Set deformed coordinates
        x1Def[e] = x1Undef[e] .+ scale*[u_n1[1]; u_n2[1]] .+ ΔuDef[1]
        x2Def[e] = x2Undef[e] .+ scale*[u_n1[2]; u_n2[2]] .+ ΔuDef[2]
        x3Def[e] = x3Undef[e] .+ scale*[u_n1[3]; u_n2[3]] .+ ΔuDef[3]
        # Skip elements without aero surfaces, or if those are not to be plotted
        if !plotAeroSurf || isnothing(aero)
            undefAirfoilCoords_n1[e],undefAirfoilCoords_n2[e],defAirfoilCoords_n1[e],defAirfoilCoords_n2[e] = nothing,nothing,nothing,nothing
            continue
        end
        # Rotation tensors
        R_n1,_ = rotation_tensor_WM(scale*p_n1)
        R_n2,_ = rotation_tensor_WM(scale*p_n2)
        # Undeformed nodal airfoil coordinates
        undefAirfoilCoords_n1[e],undefAirfoilCoords_n2[e] = get_undeformed_airfoil_coords(element)
        # Deformed nodal airfoil coordinates ( bring to origin (-r) and resolve in basis A (R_n*), then throw back to initial position (+r) and add scaled displacements (+scale*u))
        defAirfoilCoords_n1[e] = R_n1*(undefAirfoilCoords_n1[e] .- r_n1) .+ r_n1 .+ scale*u_n1 .+ ΔuDef
        defAirfoilCoords_n2[e] = R_n2*(undefAirfoilCoords_n2[e] .- r_n2) .+ r_n2 .+ scale*u_n2 .+ ΔuDef
    end

    # Set TFs for plane views
    x1Plane = maximum(abs.(vcat(x1Def...))) < tolPlane ? true : false
    x2Plane = maximum(abs.(vcat(x2Def...))) < tolPlane ? true : false
    x3Plane = maximum(abs.(vcat(x3Def...))) < tolPlane ? true : false

    # Plot undeformed and deformed beam assemblies according to view
    if isnothing(view)
        if x1Plane
            plot!(xlabel=string("\$x_2\$ [",units.length,"]"),ylabel=string("\$x_3\$ [",units.length,"]"))
            if plotUndeformed
                plot!(x2Undef, x3Undef, color=colorUndef,lw=lw,label=false)
            end
            plot!(x2Def, x3Def, color=colorDef,aspect_ratio=:equal,lw=lw,label=false)
        elseif x2Plane
            plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_3\$ [",units.length,"]"))
            if plotUndeformed
                plot!(x1Undef, x3Undef, color=colorUndef,lw=lw,label=false)
            end
            plot!(x1Def, x3Def, color=colorDef,aspect_ratio=:equal,lw=lw,label=false)
        elseif x3Plane
            plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"))
            if plotUndeformed
                plot!(x1Undef, x2Undef, color=colorUndef,lw=lw,label=false)
            end
            plot!(x1Def, x2Def, color=colorDef,aspect_ratio=:equal,lw=lw,label=false)
        else
            plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"),zlabel=string("\$x_3\$ [",units.length,"]"))
            if plotUndeformed
                plot!(x1Undef, x2Undef, x3Undef, camera=(45,45),color=colorUndef,lw=lw,label=false)
            end
            plot!(x1Def, x2Def, x3Def, camera=(45,45),color=colorDef,lw=lw,label=false)
        end
    else
        plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"),zlabel=string("\$x_3\$ [",units.length,"]"))
        if plotUndeformed
            plot!(x1Undef, x2Undef, x3Undef, camera=view,color=colorUndef,lw=lw,label=false)
        end
        plot!(x1Def, x2Def, x3Def, camera=view,color=colorDef,lw=lw,label=false)
    end

    # Plot aerodynamic surfaces, if applicable
    for e in 1:nElementsTotal
        if isnothing(defAirfoilCoords_n1[e])
            continue
        end
        # Plot undeformed aerodynamic surfaces
        if plotUndeformed
            Xu = [undefAirfoilCoords_n1[e][1,:] undefAirfoilCoords_n2[e][1,:]]
            Yu = [undefAirfoilCoords_n1[e][2,:] undefAirfoilCoords_n2[e][2,:]]
            Zu = [undefAirfoilCoords_n1[e][3,:] undefAirfoilCoords_n2[e][3,:]]
            surface!(Xu,Yu,Zu, color=palette([colorUndef,colorUndef],2), colorbar=false, alpha=surfα)
        end
        # Plot deformed aerodynamic surfaces
        Xd = [defAirfoilCoords_n1[e][1,:] defAirfoilCoords_n2[e][1,:]]
        Yd = [defAirfoilCoords_n1[e][2,:] defAirfoilCoords_n2[e][2,:]]
        Zd = [defAirfoilCoords_n1[e][3,:] defAirfoilCoords_n2[e][3,:]]
        surface!(Xd,Yd,Zd, color=palette([colorDef,colorDef],2), colorbar=false, alpha=surfα)
    end

    # Plot distributed loads (including aerodynamic and gravitational), if applicable
    if plotDistLoads
        plt = plot_distributed_loads!(plt,problem,x1Def,x2Def,x3Def,x1Plane,x2Plane,x3Plane,view,1)
    end

    # Plot generalized displacements and concentrated loads, if applicable
    if plotBCs
        plt = plot_BCs!(plt,problem,x1Def,x2Def,x3Def,x1Plane,x2Plane,x3Plane,view,0,1)
    end

    # Initialize TF for plane plot
    isPlane = isnothing(view) && (x1Plane || x2Plane || x3Plane)

    # Compute plot limits
    plotMax = max(maximum(reduce(vcat, x1Def)),maximum(reduce(vcat, x2Def)),maximum(reduce(vcat, x3Def)),maximum(reduce(vcat, x1Undef)),maximum(reduce(vcat, x2Undef)),maximum(reduce(vcat, x3Undef)))
    plotMin = min(minimum(reduce(vcat, x1Def)),minimum(reduce(vcat, x2Def)),minimum(reduce(vcat, x3Def)),minimum(reduce(vcat, x1Undef)),minimum(reduce(vcat, x2Undef)),minimum(reduce(vcat, x3Undef)))

    # Set the same plot extrema across all axis (equal aspect ratio), if applicable
    if !isPlane
        plot!(xlims=[plotMin,plotMax],ylims=[plotMin,plotMax],zlims=[plotMin,plotMax])
    end

    # Plot scale, if applicable
    if showScale
        # Display scale
        scaleString = "Deformation scale: $(string(scale))×"
        if isPlane
            annotate!(scalePos[1], scalePos[2], text(scaleString, 8))
        else
            plot!([NaN], [NaN], c=:white, lw=0, label=scaleString)
        end
    end

    # Save, if applicable
    if save
        savefig(string(pwd(),savePath))
    end

    return plt
end
export plot_steady_deformation


"""
    plot_steady_outputs(problem::Problem; kwargs...)

Plots outputs of a steady problem

# Arguments
- `problem::Problem` = problem

# Keyword arguments
- `outputs::Vector{String}` = list of outputs
- `beamGroups` = list of beams in each group for which arclengths are concatenated
- `lw::Number` = linewidth
- `colorScheme` = color scheme
- `legendPos` = legend position
- `save::Bool` = flag to save figures
- `saveFolder::String` = relative path of folder where to save figures
- `figureExtension::String` = figure extension
"""
function plot_steady_outputs(problem::Problem; outputs::Vector{String}=["u","p","F","M","V","Ω","α","cn","cm","ct","cl","cd"],beamGroups=1:length(problem.model.beams),lw::Number=1,colorScheme=:rainbow,legendPos=:best,save::Bool=false,saveFolder::String="/test/outputs/figures/",figureExtension::String=".pdf")

    # Validate
    validOutputs = ["u","u1","u2","u3","p","p1","p2","p3","F","F1","F2","F3","M","M1","M2","M3","V","V1","V2","V3","Ω","Ω1","Ω2","Ω3","γ","γ1","γ2","γ3","κ","κ1","κ2","κ3","P","P1","P2","P3","H","H1","H2","H3","α","cn","cm","ct","cl","cd"]
    @assert all(out -> out in validOutputs, outputs) "set outputs as one of $(join(validOutputs, ','))"

    # Set backend
    gr()

    # Initialize elemental and nodal variables arrays for each beam
    Nbeams = length(problem.model.beams)
    x1Nodes = [Vector{Float64}() for _ in 1:Nbeams]
    x1Elems = [Vector{Float64}() for _ in 1:Nbeams]
    uNodes = [Vector{Float64}() for _ in 1:Nbeams]
    pNodes = [Vector{Float64}() for _ in 1:Nbeams]
    FNodes = [Vector{Float64}() for _ in 1:Nbeams]
    MNodes = [Vector{Float64}() for _ in 1:Nbeams]
    VElems = [Vector{Float64}() for _ in 1:Nbeams]
    ΩElems = [Vector{Float64}() for _ in 1:Nbeams]
    γElems = [Vector{Float64}() for _ in 1:Nbeams]
    κElems = [Vector{Float64}() for _ in 1:Nbeams]
    PElems = [Vector{Float64}() for _ in 1:Nbeams]
    HElems = [Vector{Float64}() for _ in 1:Nbeams]
    αₑElems = [Vector{Float64}() for _ in 1:Nbeams]
    cnElems = [Vector{Float64}() for _ in 1:Nbeams]
    cmElems = [Vector{Float64}() for _ in 1:Nbeams]
    ctElems = [Vector{Float64}() for _ in 1:Nbeams]
    clElems = [Vector{Float64}() for _ in 1:Nbeams]
    cdElems = [Vector{Float64}() for _ in 1:Nbeams]

    # Loop over elements
    for (e,element) in enumerate(problem.model.elements)
        # Parent beam ID
        beamID = element.parent.ID
        # Elemental and nodal arclength positions
        x1_n = vcat(element.x1_n1,element.x1_n2)
        x1_e = element.x1
        # Elemental and nodal states
        u = vcat(problem.nodalStatesOverσ[end][e].u_n1,problem.nodalStatesOverσ[end][e].u_n2)
        p = vcat(problem.nodalStatesOverσ[end][e].p_n1,problem.nodalStatesOverσ[end][e].p_n2)
        F = vcat(problem.nodalStatesOverσ[end][e].F_n1,problem.nodalStatesOverσ[end][e].F_n2)
        M = vcat(problem.nodalStatesOverσ[end][e].M_n1,problem.nodalStatesOverσ[end][e].M_n2)
        V = problem.elementalStatesOverσ[end][e].V
        Ω = problem.elementalStatesOverσ[end][e].Ω
        γ = problem.compElementalStatesOverσ[end][e].γ
        κ = problem.compElementalStatesOverσ[end][e].κ
        P = problem.compElementalStatesOverσ[end][e].P
        H = problem.compElementalStatesOverσ[end][e].H
        # Add to arrays
        append!(x1Nodes[beamID],x1_n)
        append!(x1Elems[beamID],x1_e)
        append!(uNodes[beamID],u)
        append!(pNodes[beamID],p)
        append!(FNodes[beamID],F)
        append!(MNodes[beamID],M)
        append!(VElems[beamID],V)
        append!(ΩElems[beamID],Ω)
        append!(γElems[beamID],γ)
        append!(κElems[beamID],κ)
        append!(PElems[beamID],P)
        append!(HElems[beamID],H)
        # Skip elements without aero
        if isnothing(element.aero)
            continue
        end
        # Aerodynamic coefficients
        αₑ = problem.aeroVariablesOverσ[end][e].flowAnglesAndRates.αₑ
        cn = problem.aeroVariablesOverσ[end][e].aeroCoefficients.cn
        cm = problem.aeroVariablesOverσ[end][e].aeroCoefficients.cm
        ct = problem.aeroVariablesOverσ[end][e].aeroCoefficients.ct
        cl = cn .* cos(αₑ) .+ ct .* sin(αₑ)
        cd = cn .* sin(αₑ) .- ct .* cos(αₑ)
        # Add to arrays
        append!(αₑElems[beamID],αₑ)
        append!(cnElems[beamID],cn)
        append!(cmElems[beamID],cm)
        append!(ctElems[beamID],ct)
        append!(clElems[beamID],cl)
        append!(cdElems[beamID],cd)
    end

    # Plot outputs
    @unpack units = problem.model
    if "u" in outputs || "u1" in outputs
        plt_u1 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=uNodes,ind=1,units=units,YLabel="u_1",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"u1",figureExtension))
        end
    end
    if "u" in outputs || "u2" in outputs
        plt_u2 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=uNodes,ind=2,units=units,YLabel="u_2",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"u2",figureExtension))
        end
    end
    if "u" in outputs || "u3" in outputs
        plt_u3 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=uNodes,ind=3,units=units,YLabel="u_3",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"u3",figureExtension))
        end
    end
    if "p" in outputs || "p1" in outputs
        plt_p1 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=pNodes,ind=1,units=units,YLabel="p_1",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"p1",figureExtension))
        end
    end
    if "p" in outputs || "p2" in outputs
        plt_p2 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=pNodes,ind=2,units=units,YLabel="p_2",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"p2",figureExtension))
        end
    end 
    if "p" in outputs || "p3" in outputs
        plt_p3 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=pNodes,ind=3,units=units,YLabel="p_3",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"p3",figureExtension))
        end
    end
    if "F" in outputs || "F1" in outputs
        plt_F1 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=FNodes,ind=1,units=units,YLabel="F_1^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"F1",figureExtension))
        end
    end
    if "F" in outputs || "F2" in outputs
        plt_F2 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=FNodes,ind=2,units=units,YLabel="F_2^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"F2",figureExtension))
        end
    end
    if "F" in outputs || "F3" in outputs
        plt_F3 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=FNodes,ind=3,units=units,YLabel="F_3^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"F3",figureExtension))
        end
    end
    if "M" in outputs || "M1" in outputs
        plt_M1 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=MNodes,ind=1,units=units,YLabel="M_1^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"M1",figureExtension))
        end
    end
    if "M" in outputs || "M2" in outputs
        plt_M2 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=MNodes,ind=2,units=units,YLabel="M_2^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"M2",figureExtension))
        end
    end
    if "M" in outputs || "M3" in outputs
        plt_M3 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=MNodes,ind=3,units=units,YLabel="M_3^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"M3",figureExtension))
        end
    end
    if "V" in outputs || "V1" in outputs
        plt_V1 = plot_output_of_x1(beamGroups,x1=x1Elems,output=VElems,ind=1,units=units,YLabel="V_1^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"V1",figureExtension))
        end
    end
    if "V" in outputs || "V2" in outputs
        plt_V2 = plot_output_of_x1(beamGroups,x1=x1Elems,output=VElems,ind=2,units=units,YLabel="V_2^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"V2",figureExtension))
        end
    end
    if "V" in outputs || "V3" in outputs
        plt_V3 = plot_output_of_x1(beamGroups,x1=x1Elems,output=VElems,ind=3,units=units,YLabel="V_3^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"V3",figureExtension))
        end
    end
    if "Ω" in outputs || "Ω1" in outputs
        plt_Ω1 = plot_output_of_x1(beamGroups,x1=x1Elems,output=ΩElems,ind=1,units=units,YLabel="\\Omega_1^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"Ω1",figureExtension))
        end
    end
    if "Ω" in outputs || "Ω2" in outputs
        plt_Ω2 = plot_output_of_x1(beamGroups,x1=x1Elems,output=ΩElems,ind=2,units=units,YLabel="\\Omega_2^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"Ω2",figureExtension))
        end
    end
    if "Ω" in outputs || "Ω3" in outputs
        plt_Ω3 = plot_output_of_x1(beamGroups,x1=x1Elems,output=ΩElems,ind=3,units=units,YLabel="\\Omega_3^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"Ω3",figureExtension))
        end
    end
    if "γ" in outputs || "γ1" in outputs
        plt_γ1 = plot_output_of_x1(beamGroups,x1=x1Elems,output=γElems,ind=1,units=units,YLabel="\\gamma_1^+",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"γ1",figureExtension))
        end
    end
    if "γ" in outputs || "γ2" in outputs
        plt_γ2 = plot_output_of_x1(beamGroups,x1=x1Elems,output=γElems,ind=2,units=units,YLabel="\\gamma_2^+",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"γ2",figureExtension))
        end
    end
    if "γ" in outputs || "γ3" in outputs
        plt_γ3 = plot_output_of_x1(beamGroups,x1=x1Elems,output=γElems,ind=3,units=units,YLabel="\\gamma_3^+",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"γ3",figureExtension))
        end
    end
    if "κ" in outputs || "κ1" in outputs
        plt_κ1 = plot_output_of_x1(beamGroups,x1=x1Elems,output=κElems,ind=1,units=units,YLabel="\\kappa_1^+",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"κ1",figureExtension))
        end
    end
    if "κ" in outputs || "κ2" in outputs
        plt_κ2 = plot_output_of_x1(beamGroups,x1=x1Elems,output=κElems,ind=2,units=units,YLabel="\\kappa_2^+",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"κ2",figureExtension))
        end
    end
    if "κ" in outputs || "κ3" in outputs
        plt_κ3 = plot_output_of_x1(beamGroups,x1=x1Elems,output=κElems,ind=3,units=units,YLabel="\\kappa_3^+",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"κ3",figureExtension))
        end
    end
    if "P" in outputs || "P1" in outputs
        plt_P1 = plot_output_of_x1(beamGroups,x1=x1Elems,output=PElems,ind=1,units=units,YLabel="P_1^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"P1",figureExtension))
        end
    end
    if "P" in outputs || "P2" in outputs
        plt_P2 = plot_output_of_x1(beamGroups,x1=x1Elems,output=PElems,ind=2,units=units,YLabel="P_2^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"P2",figureExtension))
        end
    end
    if "P" in outputs || "P3" in outputs
        plt_P3 = plot_output_of_x1(beamGroups,x1=x1Elems,output=PElems,ind=3,units=units,YLabel="P_3^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"P3",figureExtension))
        end
    end
    if "H" in outputs || "H1" in outputs
        plt_H1 = plot_output_of_x1(beamGroups,x1=x1Elems,output=HElems,ind=1,units=units,YLabel="H_1^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"H1",figureExtension))
        end
    end
    if "H" in outputs || "H2" in outputs
        plt_H2 = plot_output_of_x1(beamGroups,x1=x1Elems,output=HElems,ind=2,units=units,YLabel="H_2^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"H2",figureExtension))
        end
    end
    if "H" in outputs || "H3" in outputs
        plt_H3 = plot_output_of_x1(beamGroups,x1=x1Elems,output=HElems,ind=3,units=units,YLabel="H_3^*",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"H3",figureExtension))
        end
    end
    if "α" in outputs 
        plt_α = plot_output_of_x1(beamGroups,x1=x1Elems,output=αₑElems,ind=0,units=units,YLabel="\\alpha",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"α",figureExtension))
        end
    end
    if "cn" in outputs 
        plt_cn = plot_output_of_x1(beamGroups,x1=x1Elems,output=cnElems,ind=0,units=units,YLabel="c_n",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"cn",figureExtension))
        end
    end
    if "cm" in outputs 
        plt_cm = plot_output_of_x1(beamGroups,x1=x1Elems,output=cmElems,ind=0,units=units,YLabel="c_m",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"cm",figureExtension))
        end
    end
    if "ct" in outputs 
        plt_ct = plot_output_of_x1(beamGroups,x1=x1Elems,output=ctElems,ind=0,units=units,YLabel="c_t",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"ct",figureExtension))
        end
    end
    if "cl" in outputs 
        plt_cl = plot_output_of_x1(beamGroups,x1=x1Elems,output=clElems,ind=0,units=units,YLabel="c_l",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"cl",figureExtension))
        end
    end
    if "cd" in outputs 
        plt_cd = plot_output_of_x1(beamGroups,x1=x1Elems,output=cdElems,ind=0,units=units,YLabel="c_d",lw=lw,colorScheme=colorScheme,legendPos=legendPos)
        if save
            savefig(string(pwd(),saveFolder,"cd",figureExtension))
        end
    end

    return nothing
end
export plot_steady_outputs


# Plots output along the arclength coordinate (x1) for each beam group
function plot_output_of_x1(beamGroups; x1,output,ind,units,YLabel,colorScheme=:rainbow,lw=1,legendPos=:best)

    # Initialize multiplication factor
    γ = 1

    # Set ylabel unit
    if YLabel == "\\alpha"
        γ = 180/π
        yLabelUnit = "deg"
    elseif YLabel in ["c_n","c_m","c_t","c_l","c_d"]
        yLabelUnit = " "   
    elseif occursin("u",YLabel)
        yLabelUnit = units.length
    elseif occursin("p",YLabel)
        yLabelUnit = " "
    elseif occursin("F",YLabel)
        yLabelUnit = units.force
    elseif occursin("M",YLabel)
        yLabelUnit = string(units.force,".",units.length)
    elseif occursin("V",YLabel)
        yLabelUnit = string(units.length,"/s")
    elseif occursin("\\Omega",YLabel)
        γ = units.angle == "deg" ? 180/π : 1
        yLabelUnit = string(units.angle,"/s")
    elseif occursin("\\gamma",YLabel)
        yLabelUnit = " "
    elseif occursin("\\kappa",YLabel)
        yLabelUnit = " "
    elseif occursin("P",YLabel)
        yLabelUnit = string(units.mass,"/s")
    elseif occursin("H",YLabel)
        yLabelUnit = string(units.mass,".",units.length,"/s")       
    end

    # Define pallete
    p = palette(colorScheme,length(beamGroups))

    # Initialize plot 
    plt = plot(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$",YLabel,"\$ [",yLabelUnit,"]"))

    # Loop beam groups
    for (i,beamGroup) in enumerate(beamGroups)
        # Initialize x1 end coordinate from previous beam in the group
        x1Previous = 0
        # Initialize arrays for plot
        x1BeamGroup = Vector{Float64}()
        outputBeamGroup = Vector{Float64}()
        # Loop over beams
        for (b,beamID) in enumerate(beamGroup)
            # Range of current x1 array
            Nx1 = length(x1[beamID])
            # Adjust x1 coordinate
            append!(x1BeamGroup,x1[beamID])
            if b > 1
                x1Previous += x1[b-1][end]
                x1BeamGroup[end-Nx1+1:end] += x1Previous
            end
            # Set range of outputs array to be plotted
            range = ind > 0 ? (ind:3:length(output[beamID])) : (1:length(output[beamID]))
            # Concatenate outputs for beam group
            append!(outputBeamGroup,γ*output[beamID][range])
        end
        # Plot output
        beamGroupLabel = length(beamGroups) > 1 ? string("Beam group ",i) : false
        plt = plot!(plt, x1BeamGroup,outputBeamGroup,lw=lw,palette=p,label=beamGroupLabel,legend=legendPos)
    end

    # Display
    display(plt)

    return plt
end


"""
    plot_mode_shapes(problem::Problem; kwargs...)

Plots the mode shapes of the model in the given problem

# Arguments
- `problem::Problem` = problem

# Keyword arguments
- `plotBCs::Bool` = flag to plot BCs
- `view::Union{Nothing,Tuple{Real,Real}}` = view angles
- `nModes::Union{Nothing,Int64}` = number of modes to plot
- `scale::Number` = displacements and rotations scale
- `frequencyLabel::String` = option for frequency label (only frequency or frequency and damping)
- `lw::Number` = linewidth
- `colorSteady` = color for steadily deformed assembly
- `modalColorScheme` = color scheme for mode shapes
- `grid::Bool` = flag for grid
- `legendPos` = legend position
- `tolPlane::Number` = displacement tolerance to plot as plane
- `plotAeroSurf` = flag to plot aerodynamic surfaces
- `surfα::Float64` = transparency factor of aerodynamic surfaces 
- `save::Bool` = flag to save the figure
- `savePath::String` = relative path on which to save the figure
"""
function plot_mode_shapes(problem::Problem; plotBCs::Bool=true,view::Union{Nothing,Tuple{Real,Real}}=nothing,nModes::Union{Nothing,Int64}=nothing,scale::Number=1,frequencyLabel::String="frequency&damping",lw::Number=1,colorSteady=:black,modalColorScheme=:jet1,grid::Bool=true,legendPos=:best,tolPlane::Number=1e-8,plotAeroSurf::Bool=true,surfα::Float64=0.5,save::Bool=false,savePath::String="/test/outputs/figures/fig.pdf")

    # Validate
    @assert problem isa EigenProblem
    @assert lw > 0
    @assert tolPlane > 0
    @assert 0 < surfα <= 1
    nModes = isnothing(nModes) ? problem.nModes : nModes
    @assert nModes <= problem.nModes

    # Set backend
    pyplot()

    # Unpack
    @unpack modeShapesAbs = problem
    @unpack elements,nElementsTotal,units = problem.model
    damps = problem.dampingsOscillatory
    freqs = problem.frequenciesOscillatory 

    # Set mode colors
    modeColors = nModes == 1 ? get(colorschemes[modalColorScheme], LinRange(0,1,nModes+1)) : get(colorschemes[modalColorScheme], LinRange(0,1,nModes))

    # Set frequency multiplier
    if units.frequency == "rad/s"
        γ = 1
    elseif units.frequency == "Hz"
        γ = 1/(2π)
    elseif units.frequency == "rpm"
        γ = 60/2π
    end

    # Initialize plot
    plt = plot(grid=grid)
    plot!([NaN],[NaN], c=colorSteady, lw=lw, label="Steady", legend=legendPos)
    for m in 1:nModes
        if frequencyLabel == "frequency"
            label = string("Mode $(m): ω = $(round(γ*freqs[m],sigdigits=3)) ", units.frequency)
        elseif frequencyLabel == "frequency&damping"
            label = string("Mode $(m): σ = $(round(γ*damps[m],sigdigits=3)) ± $(round(γ*freqs[m],sigdigits=3))i ", units.frequency)
        end
        plot!([NaN],[NaN], c=modeColors[m], lw=lw, label=label)
    end

    # Initialize arrays
    x1Undef = Array{Vector{Float64}}(undef,nElementsTotal)
    x2Undef = Array{Vector{Float64}}(undef,nElementsTotal)
    x3Undef = Array{Vector{Float64}}(undef,nElementsTotal)
    x1Steady = Array{Vector{Float64}}(undef,nElementsTotal)
    x2Steady = Array{Vector{Float64}}(undef,nElementsTotal)
    x3Steady = Array{Vector{Float64}}(undef,nElementsTotal)
    undefAirfoilCoords_n1 = Array{Union{Nothing,Matrix{Float64}}}(undef,nElementsTotal)
    undefAirfoilCoords_n2 = Array{Union{Nothing,Matrix{Float64}}}(undef,nElementsTotal)
    steadyAirfoilCoords_n1 = Array{Union{Nothing,Matrix{Float64}}}(undef,nElementsTotal)
    steadyAirfoilCoords_n2 = Array{Union{Nothing,Matrix{Float64}}}(undef,nElementsTotal)

    # Loop over elements
    for (e,element) in enumerate(elements)     
        @unpack r_n1,r_n2,nodalStates,aero = element
        @unpack u_n1,u_n2,p_n1,p_n2 = nodalStates
        # Set undeformed coordinates
        x1Undef[e] = [r_n1[1]; r_n2[1]]
        x2Undef[e] = [r_n1[2]; r_n2[2]]
        x3Undef[e] = [r_n1[3]; r_n2[3]]
        # Set steady deformed coordinates
        x1Steady[e] = x1Undef[e] .+ [u_n1[1]; u_n2[1]]
        x2Steady[e] = x2Undef[e] .+ [u_n1[2]; u_n2[2]]
        x3Steady[e] = x3Undef[e] .+ [u_n1[3]; u_n2[3]]
        # Skip elements without aero surfaces, or if those are not to be plotted
        if !plotAeroSurf || isnothing(aero)
            undefAirfoilCoords_n1[e],undefAirfoilCoords_n2[e],steadyAirfoilCoords_n1[e],steadyAirfoilCoords_n2[e] = nothing,nothing,nothing,nothing
            continue
        end
        # Rotation tensors
        R_n1,_ = rotation_tensor_WM(p_n1)
        R_n2,_ = rotation_tensor_WM(p_n2)
        # Undeformed nodal airfoil coordinates
        undefAirfoilCoords_n1[e],undefAirfoilCoords_n2[e] = get_undeformed_airfoil_coords(element)
        # Steady deformed nodal airfoil coordinates ( bring to origin (-r) and resolve in basis A (R_n*), then throw back to initial position (+r) and add displacements (+u))
        steadyAirfoilCoords_n1[e] = R_n1*(undefAirfoilCoords_n1[e] .- r_n1) .+ r_n1 .+ u_n1    
        steadyAirfoilCoords_n2[e] = R_n2*(undefAirfoilCoords_n2[e] .- r_n2) .+ r_n2 .+ u_n2
    end

    # Initialize arrays
    x1Modal = Array{Vector{Float64}}(undef,nModes,nElementsTotal)
    x2Modal = Array{Vector{Float64}}(undef,nModes,nElementsTotal)
    x3Modal = Array{Vector{Float64}}(undef,nModes,nElementsTotal)
    modalAirfoilCoords_n1 = Array{Union{Nothing,Matrix{Float64}}}(undef,nModes,nElementsTotal)
    modalAirfoilCoords_n2 = Array{Union{Nothing,Matrix{Float64}}}(undef,nModes,nElementsTotal)

    # Initialize plot limits
    plotMin = plotMax = 0

    # Loop over modes
    for m in 1:nModes
        # Loop over elements
        for (e,element) in enumerate(elements)
            # Nodal displacements and rotations    
            u_n1 = modeShapesAbs[m].nodalStates[e].u_n1
            u_n2 = modeShapesAbs[m].nodalStates[e].u_n2
            p_n1 = modeShapesAbs[m].nodalStates[e].p_n1
            p_n2 = modeShapesAbs[m].nodalStates[e].p_n2
            # Set modal coordinates
            x1Modal[m,e] = x1Steady[e] .+ scale*[u_n1[1]; u_n2[1]]
            x2Modal[m,e] = x2Steady[e] .+ scale*[u_n1[2]; u_n2[2]]
            x3Modal[m,e] = x3Steady[e] .+ scale*[u_n1[3]; u_n2[3]]
            # Skip elements without aero surfaces, or if those are not to be plotted
            if !plotAeroSurf || isnothing(element.aero)
                modalAirfoilCoords_n1[m,e],modalAirfoilCoords_n2[m,e] = nothing,nothing
                continue
            end
            @unpack r_n1,r_n2 = element
            # Rotation tensors
            R_n1,_ = rotation_tensor_WM(scale*p_n1)
            R_n2,_ = rotation_tensor_WM(scale*p_n2)
            # Modal airfoil coordinates ( bring to origin (-r) and resolve in basis A (R_n*), then throw back to initial position (+r), add scaled modal displacements (+scale*u) and steady position displacement)
            modalAirfoilCoords_n1[m,e] = R_n1*(undefAirfoilCoords_n1[e] .- r_n1) .+ r_n1 .+ scale*u_n1 .+ (steadyAirfoilCoords_n1[e] - undefAirfoilCoords_n1[e])
            modalAirfoilCoords_n2[m,e] = R_n2*(undefAirfoilCoords_n2[e] .- r_n2) .+ r_n2 .+ scale*u_n2 .+ (steadyAirfoilCoords_n2[e] - undefAirfoilCoords_n2[e])
        end

        # Set TFs for plane views
        x1Plane = maximum(abs.(vcat(x1Modal[m,:]...))) < tolPlane ? true : false
        x2Plane = maximum(abs.(vcat(x2Modal[m,:]...))) < tolPlane ? true : false
        x3Plane = maximum(abs.(vcat(x3Modal[m,:]...))) < tolPlane ? true : false

        # Plot steady and modal beam assemblies according to view
        if isnothing(view)
            if x1Plane
                plot!(xlabel=string("\$x_2\$ [",units.length,"]"),ylabel=string("\$x_3\$ [",units.length,"]"))
                plot!(x2Steady, x3Steady, color=colorSteady,lw=lw,label=false)
                plot!(x2Modal[m,:], x3Modal[m,:], color=modeColors[m],aspect_ratio=:equal,lw=lw,label=false)
            elseif x2Plane
                plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_3\$ [",units.length,"]"))
                plot!(x1Steady, x3Steady, color=colorSteady,lw=lw,label=false)
                plot!(x1Modal[m,:], x3Modal[m,:], color=modeColors[m],aspect_ratio=:equal,lw=lw,label=false)
            elseif x3Plane
                plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"))
                plot!(x1Steady, x2Steady, color=colorSteady,lw=lw,label=false)
                plot!(x1Modal[m,:], x2Modal[m,:], color=modeColors[m],aspect_ratio=:equal,lw=lw,label=false)
            else
                view = (45,45)
                plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"),zlabel=string("\$x_3\$ [",units.length,"]"))
                plot!(x1Steady, x2Steady, x3Steady, camera=view,color=colorSteady,lw=lw,label=false)
                plot!(x1Modal[m,:], x2Modal[m,:], x3Modal[m,:], camera=view,color=modeColors[m],lw=lw,label=false)
            end
        else
            plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"),zlabel=string("\$x_3\$ [",units.length,"]"))
            plot!(x1Steady, x2Steady, x3Steady, camera=view,color=colorSteady,lw=lw,label=false)
            plot!(x1Modal[m,:], x2Modal[m,:], x3Modal[m,:], camera=view,color=modeColors[m],lw=lw,label=false)
        end

        # Plot generalized displacements and concentrated loads, if applicable
        if plotBCs && m==1
            plt = plot_BCs!(plt,problem,x1Steady,x2Steady,x3Steady,x1Plane,x2Plane,x3Plane,view,0,1)
        end

        # Plot aerodynamic surfaces, if applicable
        for e in 1:nElementsTotal
            if isnothing(modalAirfoilCoords_n1[m,e])
                continue
            end
            # Plot steady deformed aerodynamic surfaces
            Xu = [steadyAirfoilCoords_n1[e][1,:] steadyAirfoilCoords_n2[e][1,:]]
            Yu = [steadyAirfoilCoords_n1[e][2,:] steadyAirfoilCoords_n2[e][2,:]]
            Zu = [steadyAirfoilCoords_n1[e][3,:] steadyAirfoilCoords_n2[e][3,:]]
            surface!(Xu,Yu,Zu, color=palette([colorSteady,colorSteady],2), colorbar=false, alpha=surfα)
            # Plot modal aerodynamic surfaces
            Xd = [modalAirfoilCoords_n1[m,e][1,:] modalAirfoilCoords_n2[m,e][1,:]]
            Yd = [modalAirfoilCoords_n1[m,e][2,:] modalAirfoilCoords_n2[m,e][2,:]]
            Zd = [modalAirfoilCoords_n1[m,e][3,:] modalAirfoilCoords_n2[m,e][3,:]]
            surface!(Xd,Yd,Zd, color=palette([modeColors[m],modeColors[m]],2), colorbar=false, alpha=surfα)
        end

        # Set the same plot extrema across all axis (equal aspect ratio), if applicable
        if !(x1Plane || x2Plane || x3Plane) || !isnothing(view)
            plotMax = max(plotMax,maximum(reduce(vcat, x1Modal[m,:])),maximum(reduce(vcat, x2Modal[m,:])),maximum(reduce(vcat, x3Modal[m,:])),maximum(reduce(vcat, x1Steady)),maximum(reduce(vcat, x2Steady)),maximum(reduce(vcat, x3Steady)))
            plotMin = min(plotMin,minimum(reduce(vcat, x1Modal[m,:])),minimum(reduce(vcat, x2Modal[m,:])),minimum(reduce(vcat, x3Modal[m,:])),minimum(reduce(vcat, x1Steady)),minimum(reduce(vcat, x2Steady)),minimum(reduce(vcat, x3Steady)))
            plot!(xlims=[plotMin,plotMax],ylims=[plotMin,plotMax],zlims=[plotMin,plotMax])
        end
    end

    # Save, if applicable
    if save
        savefig(string(pwd(),savePath))
    end

    return plt
end
export plot_mode_shapes


"""
    plot_dynamic_deformation(problem::Problem; kwargs...)

Plots the animated deformation of the model in the given problem

# Arguments
- `problem::Problem` = problem

# Keyword arguments
- `refBasis::String` = reference observer basis for plot
- `plotFrequency::Int64` = frequency of time steps to plot
- `plotUndeformed::Bool` = flag to plot undeformed assembly
- `plotBCs::Bool` = flag to plot BCs
- `plotDistLoads::Bool` = flag to plot distributed loads (includes gravitational and aerodynamic loads)
- `view::Union{Nothing,Tuple{Real,Real}}` = view angles
- `fps::Number` = frame rate for gif
- `scale::Number` = displacements and rotations scale
- `frequencyLabel::String` = option for frequency label (only frequency or frequency and damping)
- `lw::Number` = linewidth
- `colorUndef` = color for undeformed assembly
- `colorDef` = color for deformed assembly
- `grid::Bool` = flag for grid
- `legendPos` = legend position
- `tolPlane::Number` = displacement tolerance to plot as plane
- `plotAeroSurf` = flag to plot aerodynamic surfaces
- `surfα::Float64` = transparency factor of aerodynamic surfaces 
- `plotLimits::Union{Nothing,Vector{Tuple{T1,T2}}}` = plot axis limits
- `save::Bool` = flag to save the figure
- `savePath::String` = relative path on which to save the figure
- `showScale::Bool` = flag to show scale on plot
- `showTimeStamp::Bool` = flag to show time stamp on plot
- `scalePos::Vector{<:Number}` = position of scale on plot
- `timeStampPos::Vector{<:Number}` = position of time stamp on plot
- `displayProgress::Bool` = flag to display progress of gif creation
"""
function plot_dynamic_deformation(problem::Problem; refBasis::String="A", plotFrequency::Int64=1,plotUndeformed::Bool=false,plotBCs::Bool=true,plotDistLoads::Bool=true,view::Union{Nothing,Tuple{Real,Real}}=nothing,fps::Number=30,scale::Number=1,lw::Number=1,colorUndef=:black,colorDef=:blue,grid::Bool=true,legendPos=:best,tolPlane::Number=1e-8,plotAeroSurf::Bool=true,surfα::Float64=0.5,plotLimits::Union{Nothing,Vector{Tuple{T1,T2}}}=nothing,save::Bool=false,savePath::String="/test/outputs/figures/fig.gif",showScale::Bool=true,showTimeStamp::Bool=true,scalePos::Vector{<:Number}=[0.1;0.05;0.05],timeStampPos::Vector{<:Number}=[0.5;0.05;0.05],displayProgress::Bool=false) where {T1<:Number,T2<:Number}

    # Validate
    @assert problem isa DynamicProblem
    @assert refBasis in ["I","A"]
    @assert plotFrequency >= 1
    @assert scale > 0
    @assert lw > 0
    @assert tolPlane > 0
    @assert 0 < surfα <= 1
    @assert fps > 0
    @assert endswith(savePath,".gif")
    @assert length(scalePos) == 3
    @assert length(timeStampPos) == 3
    if !isnothing(plotLimits)
        @assert length(plotLimits) == 3
        @assert all((plotLimits[a][1] < plotLimits[a][2]) for a in 1:3)
    end

    # Set backend
    pyplot()

    # Unpack
    @unpack timeVector,sizeOfTime = problem
    @unpack elements,nElementsTotal,units,u_A,R_A_ofTime = problem.model

    # Initialize plot
    plt = plot(grid=grid)

    # Initialize arrays
    x1Undef = Array{Vector{Float64}}(undef,nElementsTotal)
    x2Undef = Array{Vector{Float64}}(undef,nElementsTotal)
    x3Undef = Array{Vector{Float64}}(undef,nElementsTotal)
    x1Def = Array{Vector{Float64}}(undef,sizeOfTime,nElementsTotal)
    x2Def = Array{Vector{Float64}}(undef,sizeOfTime,nElementsTotal)
    x3Def = Array{Vector{Float64}}(undef,sizeOfTime,nElementsTotal)
    undefAirfoilCoords_n1 = Array{Union{Nothing,Matrix{Float64}}}(undef,nElementsTotal)
    undefAirfoilCoords_n2 = Array{Union{Nothing,Matrix{Float64}}}(undef,nElementsTotal)
    defAirfoilCoords_n1 = Array{Union{Nothing,Matrix{Float64}}}(undef,sizeOfTime,nElementsTotal)
    defAirfoilCoords_n2 = Array{Union{Nothing,Matrix{Float64}}}(undef,sizeOfTime,nElementsTotal)

    # Compute undeformed assembly variables
    for (e,element) in enumerate(elements)     
        @unpack r_n1,r_n2,nodalStates,aero = element
        @unpack u_n1,u_n2,p_n1,p_n2 = nodalStates
        # Set undeformed coordinates
        x1Undef[e] = [r_n1[1]; r_n2[1]]
        x2Undef[e] = [r_n1[2]; r_n2[2]]
        x3Undef[e] = [r_n1[3]; r_n2[3]]
        # Skip elements without aero surfaces, or if those are not to be plotted
        if !plotAeroSurf || isnothing(aero)
            undefAirfoilCoords_n1[e],undefAirfoilCoords_n2[e] = nothing,nothing
            continue
        end
        # Undeformed nodal airfoil coordinates
        undefAirfoilCoords_n1[e],undefAirfoilCoords_n2[e] = get_undeformed_airfoil_coords(element)
    end

    # Initialize plot limits, if applicable
    if isnothing(plotLimits)
        plotMin = plotMax = 0
    end

    # Loop over time
    anim = @animate for t=1:plotFrequency:length(timeVector)

        # Set current time
        timeNow = timeVector[t]

        # Set current displacement vector and rotation tensor from basis A to reference basis (either basis I or A)
        if refBasis == "I"
            uRefNow = u_A(timeNow)
            R = R_A_ofTime[t]'
        elseif refBasis == "A"
            uRefNow = zeros(3)
            R = I3
        end

        # Loop over elements
        for (e,element) in enumerate(elements)     
            @unpack r_n1,r_n2,aero = element
            u_n1 = problem.nodalStatesOverTime[t][e].u_n1
            u_n2 = problem.nodalStatesOverTime[t][e].u_n2
            # Set deformed nodal position vectors
            posNode1 = [x1Undef[e][1]+scale*u_n1[1]; x2Undef[e][1]+scale*u_n1[2]; x3Undef[e][1]+scale*u_n1[3]]
            posNode2 = [x1Undef[e][2]+scale*u_n2[1]; x2Undef[e][2]+scale*u_n2[2]; x3Undef[e][2]+scale*u_n2[3]]
            # Set deformed coordinates (including movement of reference basis) 
            x1Def[t,e] = uRefNow[1] .+ [(R*posNode1)[1]; (R*posNode2)[1]]
            x2Def[t,e] = uRefNow[2] .+ [(R*posNode1)[2]; (R*posNode2)[2]]
            x3Def[t,e] = uRefNow[3] .+ [(R*posNode1)[3]; (R*posNode2)[3]]
        end

        # Set TFs for plane views
        x1Plane = maximum(abs.(vcat(x1Def[t,:]...))) < tolPlane ? true : false
        x2Plane = maximum(abs.(vcat(x2Def[t,:]...))) < tolPlane ? true : false
        x3Plane = maximum(abs.(vcat(x3Def[t,:]...))) < tolPlane ? true : false

        # Plot undeformed and deformed beam assemblies according to view
        if isnothing(view)
            if x1Plane
                if plotUndeformed
                    plot(x2Undef, x3Undef, color=colorUndef,lw=lw,label=false)
                    plot!(x2Def[t,:], x3Def[t,:], color=colorDef,aspect_ratio=:equal,lw=lw,label=false)
                    plot!(xlabel=string("\$x_2\$ [",units.length,"]"),ylabel=string("\$x_3\$ [",units.length,"]"))
                else
                    plot([NaN], [NaN], label=false)
                    plot!(x2Def[t,:], x3Def[t,:], color=colorDef,aspect_ratio=:equal,lw=lw,label=false)
                    plot(xlabel=string("\$x_2\$ [",units.length,"]"),ylabel=string("\$x_3\$ [",units.length,"]"))
                end
            elseif x2Plane
                if plotUndeformed
                    plot(x1Undef, x3Undef, color=colorUndef,lw=lw,label=false)
                    plot!(x1Def[t,:], x3Def[t,:], color=colorDef,aspect_ratio=:equal,lw=lw,label=false)
                    plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_3\$ [",units.length,"]"))
                else
                    plot([NaN], [NaN], label=false)
                    plot!(x1Def[t,:], x3Def[t,:], color=colorDef,aspect_ratio=:equal,lw=lw,label=false)
                    plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_3\$ [",units.length,"]"))
                end
            elseif x3Plane
                if plotUndeformed
                    plot(x1Undef, x2Undef, color=colorUndef,lw=lw,label=false)
                    plot!(x1Def[t,:], x2Def[t,:], color=colorDef,aspect_ratio=:equal,lw=lw,label=false)
                    plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"))
                else
                    plot([NaN], [NaN], label=false)
                    plot!(x1Def[t,:], x2Def[t,:], color=colorDef,aspect_ratio=:equal,lw=lw,label=false)
                    plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"))
                end
            else
                if plotUndeformed
                    plot(x1Undef, x2Undef, x3Undef, camera=(45,45),color=colorUndef,lw=lw, label=false)
                    plot!(x1Def[t,:], x2Def[t,:], x3Def[t,:], camera=(45,45),color=colorDef,lw=lw,label=false)
                    plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"),zlabel=string("\$x_3\$ [",units.length,"]"))
                else
                    plot([NaN], [NaN], label=false)
                    plot!(x1Def[t,:], x2Def[t,:], x3Def[t,:], camera=(45,45),color=colorDef,lw=lw,label=false)
                    plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"),zlabel=string("\$x_3\$ [",units.length,"]"))
                end
            end
        else
            if plotUndeformed
                plot(x1Undef, x2Undef, x3Undef, camera=view,color=colorUndef,lw=lw, label=false)
                plot!(x1Def[t,:], x2Def[t,:], x3Def[t,:], camera=view,color=colorDef,lw=lw,label=false)
                plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"),zlabel=string("\$x_3\$ [",units.length,"]"))
            else
                plot([NaN], [NaN], label=false)
                plot!(x1Def[t,:], x2Def[t,:], x3Def[t,:], camera=view,color=colorDef,lw=lw,label=false)
                plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"),zlabel=string("\$x_3\$ [",units.length,"]"))
            end
        end

        # Initialize TF for plane plot
        isPlane = isnothing(view) && (x1Plane || x2Plane || x3Plane)

        # Plot aerodynamic surfaces, if applicable
        for (e,element) in enumerate(elements)
            if isnothing(undefAirfoilCoords_n1[e])
                continue
            end
            # Update flag
            isPlane = false
            # Unpack
            @unpack r_n1,r_n2,aero = element
            u_n1 = problem.nodalStatesOverTime[t][e].u_n1
            u_n2 = problem.nodalStatesOverTime[t][e].u_n2
            p_n1 = problem.nodalStatesOverTime[t][e].p_n1
            p_n2 = problem.nodalStatesOverTime[t][e].p_n2
            # Rotation tensors
            R_n1,_ = rotation_tensor_WM(scale*p_n1)
            R_n2,_ = rotation_tensor_WM(scale*p_n2)
            # Deformed nodal airfoil coordinates ( bring to origin (-r) and resolve in reference basis (R*R_n*), then throw back to initial position (+r) and add scaled displacements (+scale*u))
            defAirfoilCoords_n1[t,e] = R*R_n1*(undefAirfoilCoords_n1[e] .- r_n1) .+ r_n1 .+ scale*u_n1 .+ uRefNow
            defAirfoilCoords_n2[t,e] = R*R_n2*(undefAirfoilCoords_n2[e] .- r_n2) .+ r_n2 .+ scale*u_n2 .+ uRefNow
            # Current coordinates
            Xd = [defAirfoilCoords_n1[t,e][1,:] defAirfoilCoords_n2[t,e][1,:]]
            Yd = [defAirfoilCoords_n1[t,e][2,:] defAirfoilCoords_n2[t,e][2,:]]
            Zd = [defAirfoilCoords_n1[t,e][3,:] defAirfoilCoords_n2[t,e][3,:]]
            # Plot undeformed aerodynamic surfaces, if applicable
            if plotUndeformed
                Xu = [undefAirfoilCoords_n1[e][1,:] undefAirfoilCoords_n2[e][1,:]]
                Yu = [undefAirfoilCoords_n1[e][2,:] undefAirfoilCoords_n2[e][2,:]]
                Zu = [undefAirfoilCoords_n1[e][3,:] undefAirfoilCoords_n2[e][3,:]]
                surface!(Xu,Yu,Zu, color=palette([colorUndef,colorUndef],2), colorbar=false, alpha=surfα)
            end
            # Plot deformed aerodynamic surfaces
            surface!(Xd,Yd,Zd, color=palette([colorDef,colorDef],2), colorbar=false, alpha=surfα)
        end

        # Plot BCs (generalized displacements and concentrated loads), if applicable
        if plotBCs
            plt = plot_BCs!(plt,problem,x1Def[t,:],x2Def[t,:],x3Def[t,:],x1Plane,x2Plane,x3Plane,view,timeNow,t)
        end

        # Plot distributed loads (including aerodynamic), if applicable
        if plotDistLoads
            # Plot distributed loads (including aerodynamic)
            plt = plot_distributed_loads!(plt,problem,x1Def[t,:],x2Def[t,:],x3Def[t,:],x1Plane,x2Plane,x3Plane,view,t)
        end

        # Set the same plot extrema across all axis (equal aspect ratio), if applicable
        if isnothing(plotLimits)
            plotMax = max(plotMax,maximum(reduce(vcat, x1Def[t,:])),maximum(reduce(vcat, x2Def[t,:])),maximum(reduce(vcat, x3Def[t,:])),maximum(reduce(vcat, x1Undef)),maximum(reduce(vcat, x2Undef)),maximum(reduce(vcat, x3Undef)))
            plotMin = min(plotMin,minimum(reduce(vcat, x1Def[t,:])),minimum(reduce(vcat, x2Def[t,:])),minimum(reduce(vcat, x3Def[t,:])),minimum(reduce(vcat, x1Undef)),minimum(reduce(vcat, x2Undef)),minimum(reduce(vcat, x3Undef)))
            plot!(xlims=[plotMin,plotMax],ylims=[plotMin,plotMax],zlims=[plotMin,plotMax])
        else
            plot!(xlims=[plotLimits[1][1],plotLimits[1][2]],ylims=[plotLimits[2][1],plotLimits[2][2]],zlims=[plotLimits[3][1],plotLimits[3][2]])
        end

        # Plot time stamp, if applicable
        if showTimeStamp
            # Set rounded time value
            roundedTime = round(timeNow,sigdigits=2)
            # Set positions
            if isnothing(plotLimits)
                XPos = plotMin + timeStampPos[1]*(plotMax-plotMin)
                YPos = plotMax + timeStampPos[2]*(plotMax-plotMin)
                ZPos = plotMax + timeStampPos[3]*(plotMax-plotMin)
            else
                XPos = plotLimits[1][1] + timeStampPos[1]*(plotLimits[1][2]-plotLimits[1][1])
                YPos = plotLimits[2][2] + timeStampPos[2]*(plotLimits[2][2]-plotLimits[2][1])
                ZPos = plotLimits[3][2] + timeStampPos[3]*(plotLimits[3][2]-plotLimits[3][1])
            end
            # Display time stamp
            if isPlane
                annotate!(XPos, YPos, text("Time: $(roundedTime) s", 8))
            else
                plot!([NaN], [NaN], c=:white, lw=0, label="Time: $(roundedTime) s")
            end
        end

        # Plot scale, if applicable
        if showScale
            # Set positions
            if isnothing(plotLimits)
                XPos = plotMin + scalePos[1]*(plotMax-plotMin)
                YPos = plotMax + scalePos[2]*(plotMax-plotMin)
                ZPos = plotMax + scalePos[3]*(plotMax-plotMin)
            else
                XPos = plotLimits[1][1] + scalePos[1]*(plotLimits[1][2]-plotLimits[1][1])
                YPos = plotLimits[2][2] + scalePos[2]*(plotLimits[2][2]-plotLimits[2][1])
                ZPos = plotLimits[3][2] + scalePos[3]*(plotLimits[3][2]-plotLimits[3][1])
            end
            # Display scale
            scaleString = "Scale: $(string(scale))×"
            if isPlane
                annotate!(XPos, YPos, text(scaleString, 8))
            else
                plot!([NaN], [NaN], c=:white, lw=0, label=scaleString)
            end
        end

        # Plot legend, if applicable
        if plotUndeformed
            plot!([NaN],[NaN], c=colorUndef, lw=lw, label="Undeformed")
            plot!([NaN],[NaN], c=colorDef, lw=lw, label="Deformed")
            plot!(legend=legendPos)
        end
        
        # Show progress, if applicable
        if displayProgress
            progress = round(timeNow/timeVector[end]*100,digits=1)
            println("Animation progress: $progress %")
        end
    end

    # Save, if applicable
    if save
        gif(anim, pkgdir(AeroBeams)*savePath, fps=fps)
    end

    return plt
end
export plot_dynamic_deformation


"""
    plot_time_outputs(problem::Problem; kwargs...)

Plots outputs of a dynamic problem

# Arguments
- `problem::Problem` = problem

# Keyword arguments
- nodes::Vector{Tuple{Int64,Int64}} = global IDs of nodes for which to plot
- elements::Vector{Int64} = global IDs of elements for which to plot
- nodalOutputs::Vector{String} = nodal outputs to plot
- `elementalOutputs::Vector{String}` = elemental outputs to plot
- `lw::Number` = linewidth
- `colorScheme` = color scheme
- `showLegend::Bool` = flag to show legend
- `legendPos` = legend position
- `save::Bool` = flag to save figures
- `saveFolder::String` = relative path of folder where to save figures
- `figureExtension::String` = figure extension
"""
function plot_time_outputs(problem::Problem; nodes::Vector{Tuple{Int64,Int64}}=Vector{Tuple{Int64,Int64}}(),elements::Vector{Int64}=Vector{Int64}(),nodalOutputs::Vector{String}=["u","p","F","M"],elementalOutputs::Vector{String}=["u","p","F","M","V","Ω","α","cn","cm","ct","cl","cd"],lw::Number=1,colorScheme=:rainbow,showLegend::Bool=true,legendPos=:best,save::Bool=false,saveFolder::String="/test/outputs/figures/",figureExtension::String=".pdf")

    # Validate
    validNodalOutputs = ["u","u1","u2","u3","p","p1","p2","p3","F","F1","F2","F3","M","M1","M2","M3"]
    validElementalOutputs = ["u","u1","u2","u3","p","p1","p2","p3","F","F1","F2","F3","M","M1","M2","M3","V","V1","V2","V3","Ω","Ω1","Ω2","Ω3","γ","γ1","γ2","γ3","κ","κ1","κ2","κ3","P","P1","P2","P3","H","H1","H2","H3","α","cn","cm","ct","cl","cd","udot","udot1","udot2","udot3","pdot","pdot1","pdot2","pdot3","Vdot","Vdot1","Vdot2","Vdot3","Ωdot","Ωdot1","Ωdot2","Ωdot3"]
    @assert all(out -> out in validNodalOutputs, nodalOutputs) "set 'nodalOutputs' as one of $(join(validNodalOutputs, ','))"
    @assert all(out -> out in validElementalOutputs, elementalOutputs) "set 'elementalOutputs' as one of $(join(validElementalOutputs, ','))"
    @assert all(n -> 0 < n[1] <= problem.model.nElementsTotal, nodes) "element of 'nodes' must be between 0 and total number of elements in the model"
    @assert all(n -> n[2] in [1,2], nodes) "set 'nodes' as a tuple in which the first value is the element and the second is its local node (1 or 2)"

    # Set backend
    gr()

    # Unpack 
    @unpack timeVector,sizeOfTime = problem
    @unpack units = problem.model

    # Initialize plots
    plt_u1 = plot()
    plt_u2 = plot()
    plt_u3 = plot()
    plt_p1 = plot()
    plt_p2 = plot()
    plt_p3 = plot()
    plt_F1 = plot()
    plt_F2 = plot()
    plt_F3 = plot()
    plt_M1 = plot()
    plt_M2 = plot()
    plt_M3 = plot()
    plt_V1 = plot()
    plt_V2 = plot()
    plt_V3 = plot()
    plt_Ω1 = plot()
    plt_Ω2 = plot()
    plt_Ω3 = plot()
    plt_γ1 = plot()
    plt_γ2 = plot()
    plt_γ3 = plot()
    plt_κ1 = plot()
    plt_κ2 = plot()
    plt_κ3 = plot()
    plt_P1 = plot()
    plt_P2 = plot()
    plt_P3 = plot()
    plt_H1 = plot()
    plt_H2 = plot()
    plt_H3 = plot()
    plt_udot1 = plot()
    plt_udot2 = plot()
    plt_udot3 = plot()
    plt_pdot1 = plot()
    plt_pdot2 = plot()
    plt_pdot3 = plot()
    plt_Vdot1 = plot()
    plt_Vdot2 = plot()
    plt_Vdot3 = plot()
    plt_Ωdot1 = plot()
    plt_Ωdot2 = plot()
    plt_Ωdot3 = plot()
    plt_α = plot()
    plt_cn = plot()
    plt_cm = plot()
    plt_ct = plot()
    plt_cl = plot()
    plt_cd = plot()

    # Loop over elements
    Ne = length(elements)
    for (i,e) in enumerate(elements)
        # u1
        if "u" in elementalOutputs || "u1" in elementalOutputs
            u1 = [problem.elementalStatesOverTime[t][e].u[1] for t in 1:sizeOfTime]
            plt_u1 = plot_output_of_time!(plt_u1,t=timeVector,output=u1,element=e,units=units,YLabel="u_1",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"u1",figureExtension))
            end
            if i == Ne
                display(plt_u1)
            end
        end
        # u2
        if "u" in elementalOutputs || "u2" in elementalOutputs
            u2 = [problem.elementalStatesOverTime[t][e].u[2] for t in 1:sizeOfTime]
            plt_u2 = plot_output_of_time!(plt_u2,t=timeVector,output=u2,element=e,units=units,YLabel="u_2",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"u2",figureExtension))
            end
            if i == Ne
                display(plt_u2)
            end
        end
        # u3
        if "u" in elementalOutputs || "u3" in elementalOutputs
            u3 = [problem.elementalStatesOverTime[t][e].u[3] for t in 1:sizeOfTime]
            plt_u3 = plot_output_of_time!(plt_u3,t=timeVector,output=u3,element=e,units=units,YLabel="u_3",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"u3",figureExtension))
            end
            if i == Ne
                display(plt_u3)
            end
        end
        # p1
        if "p" in elementalOutputs || "p1" in elementalOutputs
            p1 = [problem.elementalStatesOverTime[t][e].p[1] for t in 1:sizeOfTime]
            plt_p1 = plot_output_of_time!(plt_p1,t=timeVector,output=p1,element=e,units=units,YLabel="p_1",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"p1",figureExtension))
            end
            if i == Ne
                display(plt_p1)
            end
        end
        # p2
        if "p" in elementalOutputs || "p2" in elementalOutputs
            p2 = [problem.elementalStatesOverTime[t][e].p[2] for t in 1:sizeOfTime]
            plt_p2 = plot_output_of_time!(plt_p2,t=timeVector,output=p2,element=e,units=units,YLabel="p_2",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"p2",figureExtension))
            end
            if i == Ne
                display(plt_p2)
            end
        end
        # p3
        if "p" in elementalOutputs || "p3" in elementalOutputs
            p3 = [problem.elementalStatesOverTime[t][e].p[3] for t in 1:sizeOfTime]
            plt_p3 = plot_output_of_time!(plt_p3,t=timeVector,output=p3,element=e,units=units,YLabel="p_3",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"p3",figureExtension))
            end
            if i == Ne
                display(plt_p3)
            end
        end
        # F1
        if "F" in elementalOutputs || "F1" in elementalOutputs
            F1 = [problem.elementalStatesOverTime[t][e].F[1] for t in 1:sizeOfTime]
            plt_F1 = plot_output_of_time!(plt_F1,t=timeVector,output=F1,element=e,units=units,YLabel="F_1^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"F1",figureExtension))
            end
            if i == Ne
                display(plt_F1)
            end
        end
        # F2
        if "F" in elementalOutputs || "F2" in elementalOutputs
            F2 = [problem.elementalStatesOverTime[t][e].F[2] for t in 1:sizeOfTime]
            plt_F2 = plot_output_of_time!(plt_F2,t=timeVector,output=F2,element=e,units=units,YLabel="F_2^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"F2",figureExtension))
            end
            if i == Ne
                display(plt_F2)
            end
        end
        # F3
        if "F" in elementalOutputs || "F3" in elementalOutputs
            F3 = [problem.elementalStatesOverTime[t][e].F[3] for t in 1:sizeOfTime]
            plt_F3 = plot_output_of_time!(plt_F3,t=timeVector,output=F3,element=e,units=units,YLabel="F_3^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"F3",figureExtension))
            end
            if i == Ne
                display(plt_F3)
            end
        end
        # M1
        if "M" in elementalOutputs || "M1" in elementalOutputs
            M1 = [problem.elementalStatesOverTime[t][e].M[1] for t in 1:sizeOfTime]
            plt_M1 = plot_output_of_time!(plt_M1,t=timeVector,output=M1,element=e,units=units,YLabel="M_1^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"M1",figureExtension))
            end
            if i == Ne
                display(plt_M1)
            end
        end
        # M2
        if "M" in elementalOutputs || "M2" in elementalOutputs
            M2 = [problem.elementalStatesOverTime[t][e].M[2] for t in 1:sizeOfTime]
            plt_M2 = plot_output_of_time!(plt_M2,t=timeVector,output=M2,element=e,units=units,YLabel="M_2^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"M2",figureExtension))
            end
            if i == Ne
                display(plt_M1)
            end
        end
        # M3
        if "M" in elementalOutputs || "M3" in elementalOutputs
            M3 = [problem.elementalStatesOverTime[t][e].M[3] for t in 1:sizeOfTime]
            plt_M3 = plot_output_of_time!(plt_M3,t=timeVector,output=M3,element=e,units=units,YLabel="M_3^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"M3",figureExtension))
            end
            if i == Ne
                display(plt_M3)
            end
        end
        # V1
        if "V" in elementalOutputs || "V1" in elementalOutputs
            V1 = [problem.elementalStatesOverTime[t][e].V[1] for t in 1:sizeOfTime]
            plt_V1 = plot_output_of_time!(plt_V1,t=timeVector,output=V1,element=e,units=units,YLabel="V_1^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"V1",figureExtension))
            end
            if i == Ne
                display(plt_V1)
            end
        end
        # V2
        if "V" in elementalOutputs || "V2" in elementalOutputs
            V2 = [problem.elementalStatesOverTime[t][e].V[2] for t in 1:sizeOfTime]
            plt_V2 = plot_output_of_time!(plt_V2,t=timeVector,output=V2,element=e,units=units,YLabel="V_2^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"V2",figureExtension))
            end
            if i == Ne
                display(plt_V2)
            end
        end
        # V3
        if "V" in elementalOutputs || "V3" in elementalOutputs
            V3 = [problem.elementalStatesOverTime[t][e].V[3] for t in 1:sizeOfTime]
            plt_V3 = plot_output_of_time!(plt_V3,t=timeVector,output=V3,element=e,units=units,YLabel="V_3^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"V3",figureExtension))
            end
            if i == Ne
                display(plt_V3)
            end
        end
        # Ω1
        if "Ω" in elementalOutputs || "Ω1" in elementalOutputs
            Ω1 = [problem.elementalStatesOverTime[t][e].Ω[1] for t in 1:sizeOfTime]
            plt_Ω1 = plot_output_of_time!(plt_Ω1,t=timeVector,output=Ω1,element=e,units=units,YLabel="\\Omega_1^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"Ω1",figureExtension))
            end
            if i == Ne
                display(plt_Ω1)
            end
        end
        # Ω2
        if "Ω" in elementalOutputs || "Ω2" in elementalOutputs
            Ω2 = [problem.elementalStatesOverTime[t][e].Ω[2] for t in 1:sizeOfTime]
            plt_Ω2 = plot_output_of_time!(plt_Ω2,t=timeVector,output=Ω2,element=e,units=units,YLabel="\\Omega_2^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"Ω2",figureExtension))
            end
            if i == Ne
                display(plt_Ω2)
            end
        end
        # Ω3
        if "Ω" in elementalOutputs || "Ω3" in elementalOutputs
            Ω3 = [problem.elementalStatesOverTime[t][e].Ω[3] for t in 1:sizeOfTime]
            plt_Ω3 = plot_output_of_time!(plt_Ω3,t=timeVector,output=Ω3,element=e,units=units,YLabel="\\Omega_3^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"Ω3",figureExtension))
            end
            if i == Ne
                display(plt_Ω3)
            end
        end
        # γ1
        if "γ" in elementalOutputs || "γ1" in elementalOutputs
            γ1 = [problem.compElementalStatesOverTime[t][e].γ[1] for t in 1:sizeOfTime]
            plt_γ1 = plot_output_of_time!(plt_γ1,t=timeVector,output=γ1,element=e,units=units,YLabel="\\gamma_1^+",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"γ1",figureExtension))
            end
            if i == Ne
                display(plt_γ1)
            end
        end
        # γ2
        if "γ" in elementalOutputs || "γ2" in elementalOutputs
            γ2 = [problem.compElementalStatesOverTime[t][e].γ[2] for t in 1:sizeOfTime]
            plt_γ2 = plot_output_of_time!(plt_γ2,t=timeVector,output=γ2,element=e,units=units,YLabel="\\gamma_2^+",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"γ2",figureExtension))
            end
            if i == Ne
                display(plt_γ2)
            end
        end
        # γ3
        if "γ" in elementalOutputs || "γ3" in elementalOutputs
            γ3 = [problem.compElementalStatesOverTime[t][e].γ[3] for t in 1:sizeOfTime]
            plt_γ3 = plot_output_of_time!(plt_γ3,t=timeVector,output=γ3,element=e,units=units,YLabel="\\gamma_3^+",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"γ3",figureExtension))
            end
            if i == Ne
                display(plt_γ3)
            end
        end
        # κ1
        if "κ" in elementalOutputs || "κ1" in elementalOutputs
            κ1 = [problem.compElementalStatesOverTime[t][e].κ[1] for t in 1:sizeOfTime]
            plt_κ1 = plot_output_of_time!(plt_κ1,t=timeVector,output=κ1,element=e,units=units,YLabel="\\kappa_1^+",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"κ1",figureExtension))
            end
            if i == Ne
                display(plt_κ1)
            end
        end
        # κ2
        if "κ" in elementalOutputs || "κ2" in elementalOutputs
            κ2 = [problem.compElementalStatesOverTime[t][e].κ[2] for t in 1:sizeOfTime]
            plt_κ2 = plot_output_of_time!(plt_κ2,t=timeVector,output=κ2,element=e,units=units,YLabel="\\kappa_2^+",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"κ2",figureExtension))
            end
            if i == Ne
                display(plt_κ2)
            end
        end
        # κ3
        if "κ" in elementalOutputs || "κ3" in elementalOutputs
            κ3 = [problem.compElementalStatesOverTime[t][e].κ[3] for t in 1:sizeOfTime]
            plt_κ3 = plot_output_of_time!(plt_κ3,t=timeVector,output=κ3,element=e,units=units,YLabel="\\kappa_3^+",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"κ3",figureExtension))
            end
            if i == Ne
                display(plt_κ3)
            end
        end
        # P1
        if "P" in elementalOutputs || "P1" in elementalOutputs
            P1 = [problem.compElementalStatesOverTime[t][e].P[1] for t in 1:sizeOfTime]
            plt_P1 = plot_output_of_time!(plt_P1,t=timeVector,output=P1,element=e,units=units,YLabel="P_1^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"P1",figureExtension))
            end
            if i == Ne
                display(plt_P1)
            end
        end
        # P2
        if "P" in elementalOutputs || "P2" in elementalOutputs
            P2 = [problem.compElementalStatesOverTime[t][e].P[2] for t in 1:sizeOfTime]
            plt_P2 = plot_output_of_time!(plt_P2,t=timeVector,output=P2,element=e,units=units,YLabel="P_2^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"P2",figureExtension))
            end
            if i == Ne
                display(plt_P2)
            end
        end
        # P3
        if "P" in elementalOutputs || "P3" in elementalOutputs
            P3 = [problem.compElementalStatesOverTime[t][e].P[3] for t in 1:sizeOfTime]
            plt_P3 = plot_output_of_time!(plt_P3,t=timeVector,output=P3,element=e,units=units,YLabel="P_3^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"P3",figureExtension))
            end
            if i == Ne
                display(plt_P3)
            end
        end
        # H1
        if "H" in elementalOutputs || "H1" in elementalOutputs
            H1 = [problem.compElementalStatesOverTime[t][e].H[1] for t in 1:sizeOfTime]
            plt_H1 = plot_output_of_time!(plt_H1,t=timeVector,output=H1,element=e,units=units,YLabel="H_1^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"H1",figureExtension))
            end
            if i == Ne
                display(plt_H1)
            end
        end
        # H2
        if "H" in elementalOutputs || "H2" in elementalOutputs
            H2 = [problem.compElementalStatesOverTime[t][e].H[2] for t in 1:sizeOfTime]
            plt_H2 = plot_output_of_time!(plt_H2,t=timeVector,output=H2,element=e,units=units,YLabel="H_2^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"H2",figureExtension))
            end
            if i == Ne
                display(plt_H2)
            end
        end
        # H3
        if "H" in elementalOutputs || "H3" in elementalOutputs
            H3 = [problem.compElementalStatesOverTime[t][e].H[3] for t in 1:sizeOfTime]
            plt_H3 = plot_output_of_time!(plt_H3,t=timeVector,output=H3,element=e,units=units,YLabel="H_3^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"H3",figureExtension))
            end
            if i == Ne
                display(plt_H3)
            end
        end
        # udot1
        if "udot" in elementalOutputs || "udot1" in elementalOutputs
            udot1 = [problem.elementalStatesRatesOverTime[t][e].udot[1] for t in 1:sizeOfTime]
            plt_udot1 = plot_output_of_time!(plt_udot1,t=timeVector,output=udot1,element=e,units=units,YLabel="\\dot{u}_1",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"udot1",figureExtension))
            end
            if i == Ne
                display(plt_udot1)
            end
        end
        # udot2
        if "udot" in elementalOutputs || "udot2" in elementalOutputs
            udot2 = [problem.elementalStatesRatesOverTime[t][e].udot[2] for t in 1:sizeOfTime]
            plt_udot2 = plot_output_of_time!(plt_udot2,t=timeVector,output=udot2,element=e,units=units,YLabel="\\dot{u}_2",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"udot2",figureExtension))
            end
            if i == Ne
                display(plt_udot2)
            end
        end
        # udot3
        if "udot" in elementalOutputs || "udot3" in elementalOutputs
            udot3 = [problem.elementalStatesRatesOverTime[t][e].udot[3] for t in 1:sizeOfTime]
            plt_udot3 = plot_output_of_time!(plt_udot3,t=timeVector,output=udot3,element=e,units=units,YLabel="\\dot{u}_3",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"udot3",figureExtension))
            end
            if i == Ne
                display(plt_udot3)
            end
        end
        # pdot1
        if "pdot" in elementalOutputs || "pdot1" in elementalOutputs
            pdot1 = [problem.elementalStatesRatesOverTime[t][e].pdot[1] for t in 1:sizeOfTime]
            plt_pdot1 = plot_output_of_time!(plt_pdot1,t=timeVector,output=pdot1,element=e,units=units,YLabel="\\dot{p}_1",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"pdot1",figureExtension))
            end
            if i == Ne
                display(plt_pdot1)
            end
        end
        # pdot2
        if "pdot" in elementalOutputs || "pdot2" in elementalOutputs
            pdot2 = [problem.elementalStatesRatesOverTime[t][e].pdot[2] for t in 1:sizeOfTime]
            plt_pdot2 = plot_output_of_time!(plt_pdot2,t=timeVector,output=pdot2,element=e,units=units,YLabel="\\dot{p}_2",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"pdot2",figureExtension))
            end
            if i == Ne
                display(plt_pdot2)
            end
        end
        # pdot3
        if "pdot" in elementalOutputs || "pdot3" in elementalOutputs
            pdot3 = [problem.elementalStatesRatesOverTime[t][e].pdot[3] for t in 1:sizeOfTime]
            plt_pdot3 = plot_output_of_time!(plt_pdot3,t=timeVector,output=pdot3,element=e,units=units,YLabel="\\dot{p}_3",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"pdot3",figureExtension))
            end
            if i == Ne
                display(plt_pdot3)
            end
        end   
        # Vdot1
        if "Vdot" in elementalOutputs || "Vdot1" in elementalOutputs
            Vdot1 = [problem.elementalStatesRatesOverTime[t][e].Vdot[1] for t in 1:sizeOfTime]
            plt_Vdot1 = plot_output_of_time!(plt_Vdot1,t=timeVector,output=Vdot1,element=e,units=units,YLabel="\\dot{V}_1",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"Vdot1",figureExtension))
            end
            if i == Ne
                display(plt_Vdot1)
            end
        end
        # Vdot2
        if "Vdot" in elementalOutputs || "Vdot2" in elementalOutputs
            Vdot2 = [problem.elementalStatesRatesOverTime[t][e].Vdot[2] for t in 1:sizeOfTime]
            plt_Vdot2 = plot_output_of_time!(plt_Vdot2,t=timeVector,output=Vdot2,element=e,units=units,YLabel="\\dot{V}_2",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"Vdot2",figureExtension))
            end
            if i == Ne
                display(plt_Vdot2)
            end
        end
        # Vdot3
        if "Vdot" in elementalOutputs || "Vdot3" in elementalOutputs
            Vdot3 = [problem.elementalStatesRatesOverTime[t][e].Vdot[3] for t in 1:sizeOfTime]
            plt_Vdot3 = plot_output_of_time!(plt_Vdot3,t=timeVector,output=Vdot3,element=e,units=units,YLabel="\\dot{V}_3",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"Vdot3",figureExtension))
            end
            if i == Ne
                display(plt_Vdot3)
            end
        end
        # Ωdot1
        if "Ωdot" in elementalOutputs || "Ωdot1" in elementalOutputs
            Ωdot1 = [problem.elementalStatesRatesOverTime[t][e].Ωdot[1] for t in 1:sizeOfTime]
            plt_Ωdot1 = plot_output_of_time!(plt_Ωdot1,t=timeVector,output=Ωdot1,element=e,units=units,YLabel="\\dot{\\Omega}_1",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"Ωdot1",figureExtension))
            end
            if i == Ne
                display(plt_Ωdot1)
            end
        end
        # Ωdot2
        if "Ωdot" in elementalOutputs || "Ωdot2" in elementalOutputs
            Ωdot2 = [problem.elementalStatesRatesOverTime[t][e].Ωdot[2] for t in 1:sizeOfTime]
            plt_Ωdot2 = plot_output_of_time!(plt_Ωdot2,t=timeVector,output=Ωdot2,element=e,units=units,YLabel="\\dot{\\Omega}_2",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"Ωdot2",figureExtension))
            end
            if i == Ne
                display(plt_Ωdot2)
            end
        end
        # Ωdot3
        if "Ωdot" in elementalOutputs || "Ωdot3" in elementalOutputs
            Ωdot3 = [problem.elementalStatesRatesOverTime[t][e].Ωdot[3] for t in 1:sizeOfTime]
            plt_Ωdot3 = plot_output_of_time!(plt_Ωdot3,t=timeVector,output=Ωdot3,element=e,units=units,YLabel="\\dot{\\Omega}_3",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"Ωdot3",figureExtension))
            end
            if i == Ne
                display(plt_Ωdot3)
            end
        end
        # Skip elements without aero
        if isnothing(problem.model.elements[e].aero)
            continue
        end
        # Aerodynamic variables
        αₑ = [problem.aeroVariablesOverTime[t][e].flowAnglesAndRates.αₑ for t in 1:sizeOfTime]
        cn = [problem.aeroVariablesOverTime[t][e].aeroCoefficients.cn for t in 1:sizeOfTime]
        cm = [problem.aeroVariablesOverTime[t][e].aeroCoefficients.cm for t in 1:sizeOfTime]
        ct = [problem.aeroVariablesOverTime[t][e].aeroCoefficients.ct for t in 1:sizeOfTime]
        cl = cn .* cos.(αₑ) .+ ct .* sin.(αₑ)
        cd = cn .* sin.(αₑ) .- ct .* cos.(αₑ)
        # α
        if "α" in elementalOutputs
            plt_α = plot_output_of_time!(plt_α,t=timeVector,output=αₑ,element=e,units=units,YLabel="\\alpha",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"α",figureExtension))
            end
            if i == Ne
                display(plt_α)
            end
        end
        # cn
        if "cn" in elementalOutputs
            plt_cn = plot_output_of_time!(plt_cn,t=timeVector,output=cn,element=e,units=units,YLabel="c_n",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"cn",figureExtension))
            end
            if i == Ne
                display(plt_cn)
            end
        end
        # cm
        if "cm" in elementalOutputs
            plt_cm = plot_output_of_time!(plt_cm,t=timeVector,output=cm,element=e,units=units,YLabel="c_m",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"cm",figureExtension))
            end
            if i == Ne
                display(plt_cm)
            end
        end
        # ct
        if "ct" in elementalOutputs
            plt_ct = plot_output_of_time!(plt_ct,t=timeVector,output=ct,element=e,units=units,YLabel="c_t",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"ct",figureExtension))
            end
            if i == Ne
                display(plt_ct)
            end
        end
        # cl
        if "cl" in elementalOutputs
            plt_cl = plot_output_of_time!(plt_cl,t=timeVector,output=cl,element=e,units=units,YLabel="c_l",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"cl",figureExtension))
            end
            if i == Ne
                display(plt_cl)
            end
        end
        # cd
        if "cd" in elementalOutputs
            plt_cd = plot_output_of_time!(plt_cd,t=timeVector,output=cd,element=e,units=units,YLabel="c_d",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"cd",figureExtension))
            end
            if i == Ne
                display(plt_cd)
            end
        end
    end

    # Loop over nodal tuples
    Nn = length(nodes)
    for (i,tup) in enumerate(nodes)
        # Extract element and local node
        e = tup[1]
        localNode = tup[2]
        # Extract global node ID
        nodeGlobalID = problem.model.elements[e].nodesGlobalID[localNode]
        # Set fields
        if localNode == 1
            uField = :u_n1
            pField = :p_n1
            FField = :F_n1
            MField = :N_n1
        else
            uField = :u_n2
            pField = :p_n2
            FField = :F_n2
            MField = :M_n2
        end
        # u1
        if "u" in nodalOutputs || "u1" in nodalOutputs
            u1 = [getfield(problem.nodalStatesOverTime[t][e],uField)[1] for t in 1:sizeOfTime]
            plt_u1 = plot_output_of_time!(plt_u1,t=timeVector,output=u1,node=nodeGlobalID,units=units,YLabel="u_1",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"u1",figureExtension))
            end
            if i == Nn
                display(plt_u1)
            end
        end
        # u2
        if "u" in nodalOutputs || "u2" in nodalOutputs
            u2 = [getfield(problem.nodalStatesOverTime[t][e],uField)[2] for t in 1:sizeOfTime]
            plt_u2 = plot_output_of_time!(plt_u2,t=timeVector,output=u2,node=nodeGlobalID,units=units,YLabel="u_2",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"u2",figureExtension))
            end
            if i == Nn
                display(plt_u2)
            end
        end
        # u3
        if "u" in nodalOutputs || "u3" in nodalOutputs
            u3 = [getfield(problem.nodalStatesOverTime[t][e],uField)[3] for t in 1:sizeOfTime]
            plt_u3 = plot_output_of_time!(plt_u3,t=timeVector,output=u3,node=nodeGlobalID,units=units,YLabel="u_3",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"u3",figureExtension))
            end
            if i == Nn
                display(plt_u3)
            end
        end
        # p1
        if "p" in nodalOutputs || "p1" in nodalOutputs
            p1 = [getfield(problem.nodalStatesOverTime[t][e],pField)[1] for t in 1:sizeOfTime]
            plt_p1 = plot_output_of_time!(plt_p1,t=timeVector,output=p1,node=nodeGlobalID,units=units,YLabel="p_1",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"p1",figureExtension))
            end
            if i == Nn
                display(plt_p1)
            end
        end
        # p2
        if "p" in nodalOutputs || "p2" in nodalOutputs
            p2 = [getfield(problem.nodalStatesOverTime[t][e],pField)[2] for t in 1:sizeOfTime]
            plt_p2 = plot_output_of_time!(plt_p2,t=timeVector,output=p2,node=nodeGlobalID,units=units,YLabel="p_2",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"p2",figureExtension))
            end
            if i == Nn
                display(plt_p2)
            end
        end
        # p3
        if "p" in nodalOutputs || "p3" in nodalOutputs
            p3 = [getfield(problem.nodalStatesOverTime[t][e],pField)[3] for t in 1:sizeOfTime]
            plt_p3 = plot_output_of_time!(plt_p3,t=timeVector,output=u3,node=nodeGlobalID,units=units,YLabel="p_3",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"p3",figureExtension))
            end
            if i == Nn
                display(plt_p3)
            end
        end
        # F1
        if "F" in nodalOutputs || "F1" in nodalOutputs
            F1 = [getfield(problem.nodalStatesOverTime[t][e],FField)[1] for t in 1:sizeOfTime]
            plt_F1 = plot_output_of_time!(plt_F1,t=timeVector,output=F1,node=nodeGlobalID,units=units,YLabel="F_1^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"F1",figureExtension))
            end
            if i == Nn
                display(plt_F1)
            end
        end
        # F2
        if "F" in nodalOutputs || "F2" in nodalOutputs
            F2 = [getfield(problem.nodalStatesOverTime[t][e],FField)[2] for t in 1:sizeOfTime]
            plt_F2 = plot_output_of_time!(plt_F2,t=timeVector,output=F2,node=nodeGlobalID,units=units,YLabel="F_2^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"F2",figureExtension))
            end
            if i == Nn
                display(plt_F2)
            end
        end
        # F3
        if "F" in nodalOutputs || "F3" in nodalOutputs
            F3 = [getfield(problem.nodalStatesOverTime[t][e],FField)[3] for t in 1:sizeOfTime]
            plt_F3 = plot_output_of_time!(plt_F3,t=timeVector,output=F3,node=nodeGlobalID,units=units,YLabel="F_3^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"F3",figureExtension))
            end
            if i == Nn
                display(plt_F3)
            end
        end
        # M1
        if "M" in nodalOutputs || "M1" in nodalOutputs
            M1 = [getfield(problem.nodalStatesOverTime[t][e],MField)[1] for t in 1:sizeOfTime]
            plt_M1 = plot_output_of_time!(plt_M1,t=timeVector,output=M1,node=nodeGlobalID,units=units,YLabel="M_1^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"M1",figureExtension))
            end
            if i == Nn
                display(plt_M1)
            end
        end
        # M2
        if "M" in nodalOutputs || "M2" in nodalOutputs
            M2 = [getfield(problem.nodalStatesOverTime[t][e],MField)[2] for t in 1:sizeOfTime]
            plt_M2 = plot_output_of_time!(plt_M2,t=timeVector,output=M2,node=nodeGlobalID,units=units,YLabel="M_2^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"M2",figureExtension))
            end
            if i == Nn
                display(plt_M2)
            end
        end
        # M3
        if "M" in nodalOutputs || "M3" in nodalOutputs
            M3 = [getfield(problem.nodalStatesOverTime[t][e],MField)[3] for t in 1:sizeOfTime]
            plt_M3 = plot_output_of_time!(plt_M3,t=timeVector,output=M3,node=nodeGlobalID,units=units,YLabel="M_3^*",lw=lw,colorScheme=colorScheme,showLegend=showLegend,legendPos=legendPos)
            if save
                savefig(string(pwd(),saveFolder,"M3",figureExtension))
            end
            if i == Nn
                display(plt_M3)
            end
        end
    end

    return nothing
end
export plot_time_outputs


# Plots output over time
function plot_output_of_time!(plt; t,output,element=NaN,node=NaN,units,YLabel,lw=1,colorScheme=:rainbow,showLegend=true,legendPos=:best)

    # Validate: either "element" or "node" is NaN
    @assert isnan(element) || isnan(node)
    @assert !(isnan(element) && isnan(node))

    # Initialize multiplication factor
    γ = 1

    # Set ylabel unit
    if YLabel == "\\alpha"
        γ = 180/π
        yLabelUnit = "deg"
    elseif YLabel in ["c_n","c_m","c_t","c_l","c_d"]
        yLabelUnit = " "   
    elseif occursin("\\dot{u}",YLabel)
        yLabelUnit = string(units.length,"/s")
    elseif occursin("\\dot{p}",YLabel)
        yLabelUnit = "1/s"
    elseif occursin("\\dot{V}",YLabel)
        yLabelUnit = string(units.length,"\$/s^2\$")
    elseif occursin("\\dot{\\Omega}",YLabel)
        γ = units.angle == "deg" ? 180/π : 1
        yLabelUnit = string(units.angle,"\$/s^2\$")    
    elseif occursin("u",YLabel)
        yLabelUnit = units.length
    elseif occursin("p",YLabel)
        yLabelUnit = " "
    elseif occursin("F",YLabel)
        yLabelUnit = units.force
    elseif occursin("M",YLabel)
        yLabelUnit = string(units.force,".",units.length)
    elseif occursin("V",YLabel)
        yLabelUnit = string(units.length,"/s")
    elseif occursin("\\Omega",YLabel)
        γ = units.angle == "deg" ? 180/π : 1
        yLabelUnit = string(units.angle,"/s")
    elseif occursin("\\gamma",YLabel)
        yLabelUnit = " "
    elseif occursin("\\kappa",YLabel)
        yLabelUnit = " "
    elseif occursin("P",YLabel)
        yLabelUnit = string(units.mass,"/s")
    elseif occursin("H",YLabel)
        yLabelUnit = string(units.mass,".",units.length,"/s")       
    end

    # Define pallete
    p = palette(colorScheme,length(plt.series_list)+1)

    # Set axis labels
    plt = plot!(plt, xlabel=string("\$t\$ [s]"),ylabel=string("\$",YLabel,"\$ [",yLabelUnit,"]"))

    # Set legend
    if showLegend
        currentLegend = !isnan(element) ? string("Element ", element) : string("Node ", node)
    else
        currentLegend = false
    end

    # Plot output
    plt = plot!(plt, t,γ*output,lw=lw,palette=p,label=currentLegend,legend=legendPos)

    # Update colors from previous plots
    [(plt.series_list[i][:linecolor] = p[i]) for i in 1:length(plt.series_list)]

    return plt
end


# Computes the undeformed nodal airfoil coordinates
function get_undeformed_airfoil_coords(element::Element)

    @unpack x1_n1,x1_n2,r_n1,r_n2,R0_n1,R0_n2 = element
    @unpack c,normSparPos = element.parent.aeroSurface
    @unpack airfoil,Rw = element.aero

    # Airfoil coordinates in the X-plane
    Y,Z = airfoil.coordinates[:,1],airfoil.coordinates[:,2]

    # Number of points
    N = length(Y)

    # Nodal chords and normalized spar positions
    if c isa Number
        c1 = c2 = c
    else
        c1,c2 = c(x1_n1),c(x1_n2)
    end
    if normSparPos isa Number
        normSparPos1 = normSparPos2 = normSparPos
    else
        normSparPos1,normSparPos2 = normSparPos(x1_n1),normSparPos(x1_n2)
    end

    # Rotate to face forward
    Y = 1 .- Y

    # Offset to spar position
    Y1,Y2 = Y .- (1 .- normSparPos1), Y .- (1 .- normSparPos2)
    Z1 = Z2 = Z

    # Scale by chord
    Y1 *= c1
    Z1 *= c1
    Y2 *= c2
    Z2 *= c2

    # Rotate to match initial nodal orientation and translate to initial position
    XYZRot1,XYZRot2 = Rw*R0_n1*[zeros(N)'; Y1'; Z1'], Rw*R0_n2*[zeros(N)'; Y2'; Z2']           
    X1,X2 = XYZRot1[1,:]' .+ r_n1[1], XYZRot2[1,:]' .+ r_n2[1]
    Y1,Y2 = XYZRot1[2,:]' .+ r_n1[2], XYZRot2[2,:]' .+ r_n2[2]
    Z1,Z2 = XYZRot1[3,:]' .+ r_n1[3], XYZRot2[3,:]' .+ r_n2[3]

    # Set undeformed coordinates of both nodes
    X = vcat(X1,X2)
    Y = vcat(Y1,Y2)
    Z = vcat(Z1,Z2)
    undefAirfoilCoords = [X; Y; Z]

    # Separate by node
    undefAirfoilCoords_n1 = Matrix(hcat(undefAirfoilCoords[1,:],undefAirfoilCoords[3,:],undefAirfoilCoords[5,:])')
    undefAirfoilCoords_n2 = Matrix(hcat(undefAirfoilCoords[2,:],undefAirfoilCoords[4,:],undefAirfoilCoords[6,:])')

    return undefAirfoilCoords_n1,undefAirfoilCoords_n2
end


# Plots all boundary conditions at the current time
function plot_BCs!(plt,problem,x1Def,x2Def,x3Def,x1Plane,x2Plane,x3Plane,view,timeNow,tIndNow)

    @unpack BCs = problem.model
    
    # Loop BCs
    for BC in BCs
        # Update BCs
        update_BC_data!(BC,timeNow)
        @unpack isLoad,deadLoadsOnA,followerLoadsOnA,deadLoadsOnb,followerLoadsOnb,beam,node,R0_n,Fmax,Mmax = BC
        # Reference length is beam length
        L = beam.length
        # Find global ID of first element that contains the BC'ed node
        eLocalID = node == 1 ? 1 : node-1
        localNode = node == 1 ? 1 : 2
        element = beam.elements[eLocalID]
        eGlobalID = element.globalID
        # Nodal rotation tensor from basis b to basis B
        if problem isa DynamicProblem
            R_n,_ = localNode == 1 ? rotation_tensor_WM(problem.nodalStatesOverTime[tIndNow][eGlobalID].p_n1) : rotation_tensor_WM(problem.nodalStatesOverTime[tIndNow][eGlobalID].p_n2)
        else
            R_n,_ = localNode == 1 ? rotation_tensor_WM(element.nodalStates.p_n1) : rotation_tensor_WM(element.nodalStates.p_n2)
        end
        # Position where BC is applied
        P2 = [x1Def[eGlobalID][localNode]; x2Def[eGlobalID][localNode]; x3Def[eGlobalID][localNode]]
        # Draw BC
        plt = draw_BC!(plt,isLoad=isLoad,deadLoadsOnA=deadLoadsOnA,followerLoadsOnA=followerLoadsOnA,deadLoadsOnb=deadLoadsOnb,followerLoadsOnb=followerLoadsOnb,R0_n=R0_n,R_n=R_n,Fmax=Fmax,Mmax=Mmax,L=L,P2=P2,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
    end

    return plt

end


# Draws boundary condition
function draw_BC!(plt; isLoad,deadLoadsOnA,followerLoadsOnA,deadLoadsOnb,followerLoadsOnb,R0_n,R_n,Fmax,Mmax,L,P2,x1Plane,x2Plane,x3Plane,view)

    # Loop DOFs
    for DOF in 1:6
        # Draw generalized displacement BC
        if !isLoad[DOF]
            # Displacement DOFs
            if DOF in [1,2,3]
                plt = draw_generalized_displacement!(plt,axis=DOF,P2=P2,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
            # Rotation DOFs    
            else
                Δ = -[1 0 0; 0 1 0; 0 0 1][:, DOF-3]
                plt = draw_generalized_displacement!(plt,axis=DOF-3,P2=P2,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
                plt = draw_generalized_displacement!(plt,axis=DOF-3,P2=P2,Δ=Δ,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view,color=:orange)
            end
        # Draw generalized forces BC    
        else
            # Dead loads initially resolved on basis A
            if !iszero(deadLoadsOnA[DOF])
                if 1 <= DOF <= 3
                    plt = draw_concentrated_force!(plt,DOF=DOF,F=deadLoadsOnA[DOF],R=I3,P2=P2,Fmax=Fmax,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
                else
                    plt = draw_concentrated_moment!(plt,DOF=DOF,M=deadLoadsOnA[DOF],R=I3,P2=P2,Mmax=Mmax,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
                end
            end
            # Follower loads initially resolved on basis A
            if !iszero(followerLoadsOnA[DOF])
                if 1 <= DOF <= 3
                    plt = draw_concentrated_force!(plt,DOF=DOF,F=followerLoadsOnA[DOF],R=R_n,P2=P2,Fmax=Fmax,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
                else
                    plt = draw_concentrated_moment!(plt,DOF=DOF,M=followerLoadsOnA[DOF],R=R_n,P2=P2,Mmax=Mmax,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
                end
            end
            # Dead loads initially resolved on basis b
            if !iszero(deadLoadsOnb[DOF])
                if 1 <= DOF <= 3
                    plt = draw_concentrated_force!(plt,DOF=DOF,F=deadLoadsOnb[DOF],R=R0_n,P2=P2,Fmax=Fmax,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
                else
                    plt = draw_concentrated_moment!(plt,DOF=DOF,M=deadLoadsOnb[DOF],R=R0_n,P2=P2,Mmax=Mmax,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
                end
            end
            # Follower loads initially resolved on basis b
            if !iszero(followerLoadsOnb[DOF])
                if 1 <= DOF <= 3
                    plt = draw_concentrated_force!(plt,DOF=DOF,F=followerLoadsOnb[DOF],R=R_n*R0_n,P2=P2,Fmax=Fmax,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
                else
                    plt = draw_concentrated_moment!(plt,DOF=DOF,M=followerLoadsOnb[DOF],R=R_n*R0_n,P2=P2,Mmax=Mmax,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
                end
            end
        end
    end

    return plt
end


# Draws concentrated force
function draw_concentrated_force!(plt; DOF,F,R,P2,Fmax,L,x1Plane,x2Plane,x3Plane,view,color=:green,λ=1/10)

    # Get direction vector
    if DOF == 1
        dir = a1
    elseif DOF == 2
        dir = a2
    elseif DOF == 3
        dir = a3
    end   

    # Vectors for quiver
    P1 = P2 .- λ*L*F/Fmax*(R*dir)
    ΔP = P2 .- P1

    # Plot
    if x1Plane && isnothing(view)
        quiver!([P1[2]], [P1[3]], quiver=([ΔP[2]], [ΔP[3]]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x2Plane && isnothing(view)
        quiver!([P1[1]], [P1[3]], quiver=([ΔP[1]], [ΔP[3]]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x3Plane && isnothing(view)
        quiver!([P1[1]], [P1[2]], quiver=([ΔP[1]], [ΔP[2]]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
    else
        quiver!([P1[1]], [P1[2]], [P1[3]], quiver=([ΔP[1]], [ΔP[2]], [ΔP[3]]), color=color, lw=1, quiverhead=0.5)
    end

    return plt
end


# Draws concentrated moment
function draw_concentrated_moment!(plt; DOF,M,R,P2,Mmax,L,x1Plane,x2Plane,x3Plane,view,color=:green,λ=1)

    # Get direction vector
    if DOF == 4
        dir = a1
    elseif DOF == 5
        dir = a2
    elseif DOF == 6
        dir = a3
    end   

    # Plot
    plt = draw_circular_vector!(plt, origin=P2,M=M/Mmax,R=R,L=L*λ,axis=DOF-3,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view,color=color)

    return plt
end


# Draws generalized displacement boundary condition
function draw_generalized_displacement!(plt; axis,P2,L,x1Plane,x2Plane,x3Plane,view,Δ=zeros(3),λ=1/200,Ndiv=21,color=:red)

    # Set x and y triangle vertices for 2D view, if applicable
    if isnothing(view) && (x1Plane || x2Plane || x3Plane)
        if axis == 1
            xpoints = [0; -1; -1]
            ypoints = [0; 1/2; -1/2]
        elseif axis == 2
            xpoints = cos.(LinRange(0, 2π, Ndiv))./2
            ypoints = sin.(LinRange(0, 2π, Ndiv))./2
        elseif axis == 3
            xpoints = [0; 1/2; -1/2]
            ypoints = [0; -1; -1]
        end
    end

    # Set triangle for 2D view or cone for 3D view
    if x1Plane && isnothing(view)
        x = P2[2] + π*λ*L*Δ[2] .+ 2π*L*λ*xpoints
        y = P2[3] + π*λ*L*Δ[3] .+ 2π*L*λ*ypoints
    elseif x2Plane && isnothing(view)
        x = P2[1] + π*λ*L*Δ[1]  .+ 2π*L*λ*xpoints
        y = P2[3] + π*λ*L*Δ[3]  .+ 2π*L*λ*ypoints
    elseif x3Plane && isnothing(view)
        x = P2[1] + π*λ*L*Δ[1]  .+ 2π*L*λ*xpoints
        y = P2[2] + π*λ*L*Δ[2]  .+ 2π*L*λ*ypoints
    else
        # Cone drawing variables
        r = LinRange(0, 2π, Ndiv)
        θ = LinRange(0, 2π, Ndiv)
        H = repeat(r', length(θ), 1)
        T = repeat(θ, 1, length(r))
        # Cone coordinates
        x1Cone = -L*λ/2 * H .* cos.(T)
        x2Cone = -L*λ/2 * H .* sin.(T)
        x3Cone = -L*λ * H
        # Change according to axis
        x1Cone, x2Cone, x3Cone = if axis == 1
            x3Cone, x2Cone, x1Cone
        elseif axis == 2
            x2Cone, x3Cone, x1Cone
        elseif axis == 3
            x1Cone, x2Cone, x3Cone
        else
            error("Unexpected axis: $axis")
        end
    end
    
    # Plot
    if isnothing(view) && (x1Plane || x2Plane || x3Plane)
        plot!(x, y, seriestype=:shape, fillcolor=color, linecolor=color, label=false)
    else
        surface!(P2[1] + π*λ*L*Δ[1] .+ x1Cone, P2[2] + π*λ*L*Δ[2] .+ x2Cone, P2[3] + π*λ*L*Δ[3] .+ x3Cone, color=palette([color,color],2), colorbar=false)
    end

    return plt
end


# Plots all distributed loads at the current time
function plot_distributed_loads!(plt,problem::Problem,x1Def,x2Def,x3Def,x1Plane,x2Plane,x3Plane,view,tIndNow)

    @unpack timeNow,maxAeroForce,maxAeroMoment = problem
    @unpack gravityVector,elements = problem.model

    # Loop elements
    for (e,element) in enumerate(elements)

        @unpack aero,R0,R0_n1,R0_n2,μ,hasDistributedDeadForcesBasisA,hasDistributedDeadMomentsBasisA,hasDistributedDeadForcesBasisb,hasDistributedDeadMomentsBasisb,hasDistributedFollowerForcesBasisA,hasDistributedFollowerMomentsBasisA,hasDistributedFollowerForcesBasisb,hasDistributedFollowerMomentsBasisb,f_A_of_ζt,m_A_of_ζt,f_b_of_ζt,m_b_of_ζt,ff_A_of_ζt,mf_A_of_ζt,ff_b_of_ζt,mf_b_of_ζt,globalID = element

        # Reference length = beam length
        L = element.parent.length

        # Skip if there are no distributed loads on current element
        if iszero([hasDistributedDeadForcesBasisA; hasDistributedDeadMomentsBasisA; hasDistributedDeadForcesBasisb; hasDistributedDeadMomentsBasisb; hasDistributedFollowerForcesBasisA; hasDistributedFollowerMomentsBasisA; hasDistributedFollowerForcesBasisb; hasDistributedFollowerMomentsBasisb]) && isnothing(aero) && iszero(gravityVector)
            continue
        end

        # Interpolated coordinates of deformed structure at element's nodes and midpoint
        x1DefInterp = [x1Def[e][1]; (x1Def[e][1]+x1Def[e][2])/2; x1Def[e][2]]
        x2DefInterp = [x2Def[e][1]; (x2Def[e][1]+x2Def[e][2])/2; x2Def[e][2]]
        x3DefInterp = [x3Def[e][1]; (x3Def[e][1]+x3Def[e][2])/2; x3Def[e][2]]
        P2 = hcat(x1DefInterp, x2DefInterp, x3DefInterp)

        # Interpolated local coordinate at element's nodes and midpoint
        ζInterp = [-1; 0; 1]

        # Number of divisions for vectors (element's nodes and midpoint)
        Ndiv = length(ζInterp)

        # Interpolation arclength coordinate for estimation of maximum loads
        x1Interp = LinRange(0,L,element.parent.nElements*3)

        # Compute rotation tensors at element's nodes and midpoint, if applicable
        if !iszero([hasDistributedFollowerForcesBasisA; hasDistributedFollowerMomentsBasisA; hasDistributedFollowerForcesBasisb; hasDistributedFollowerMomentsBasisb]) || !isnothing(aero)
            if problem isa DynamicProblem
                p = problem.elementalStatesOverTime[tIndNow][e].p
                p_n1 = problem.nodalStatesOverTime[tIndNow][e].p_n1
                p_n2 = problem.nodalStatesOverTime[tIndNow][e].p_n2
                R_e,_ = rotation_tensor_WM(p)
                R_n1,_ = rotation_tensor_WM(p_n1)
                R_n2,_ = rotation_tensor_WM(p_n2)
            else
                @unpack R = element
                @unpack p_n1,p_n2 = element.nodalStates
                R_e = R
                R_n1,_ = rotation_tensor_WM(p_n1)
                R_n2,_ = rotation_tensor_WM(p_n2)
            end
        end

        ## Gravity loads
        # ----------------------------------------------------------------------
        if !iszero(gravityVector) && μ > 0
            for (dir,ai) in zip(1:3,[a1, a2, a3])
                Fi = repeat([gravityVector[dir]],3)
                if maximum(abs.(Fi)) == 0
                    continue
                end
                plt = draw_distributed_forces!(plt, F=Fi,R=[I3,I3,I3],ai=ai,P2=P2,Ndiv=Ndiv,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view,color=:gold)
            end
        end

        ## Dead distributed forces initially resolved in basis A
        # ----------------------------------------------------------------------
        if hasDistributedDeadForcesBasisA
            currentLoad = f_A_of_ζt.(ζInterp,timeNow)
            Fmax = maximum(abs.(element.parent.f_A_of_x1t(x1Interp,timeNow)))
            for (dir,ai) in zip(1:3,[a1, a2, a3])
                Fi = [currentLoad[d][dir] for d in 1:3]
                if maximum(abs.(Fi)) == 0
                    continue
                end
                plt = draw_distributed_forces!(plt, F=Fi,Fmax=Fmax,R=[I3,I3,I3],ai=ai,P2=P2,Ndiv=Ndiv,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
            end
        end

        ## Dead distributed moments initially resolved in basis A
        # ----------------------------------------------------------------------
        if hasDistributedDeadMomentsBasisA
            currentLoad = m_A_of_ζt.(ζInterp,timeNow)
            Mmax = maximum(abs.(element.parent.m_A_of_x1t(x1Interp,timeNow)))
            for (dir,ai) in zip(1:3,[a1, a2, a3])
                Mi = [currentLoad[d][dir] for d in 1:3]
                if maximum(abs.(Mi)) == 0
                    continue
                end
                plt = draw_distributed_moments!(plt, M=Mi,Mmax=Mmax,R=[I3,I3,I3],axis=dir,P2=P2,Ndiv=Ndiv,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
            end
        end

        ## Dead distributed forces initially resolved in basis b
        # ----------------------------------------------------------------------
        if hasDistributedDeadForcesBasisb
            currentLoad = f_b_of_ζt.(ζInterp,timeNow)
            Fmax = maximum(abs.(element.parent.f_b_of_x1t(x1Interp,timeNow)))
            for (dir,ai) in zip(1:3,[a1, a2, a3])
                Fi = [currentLoad[d][dir] for d in 1:3]
                if maximum(abs.(Fi)) == 0
                    continue
                end
                plt = draw_distributed_forces!(plt, F=Fi,Fmax=Fmax,R=[R0,R0,R0],ai=ai,P2=P2,Ndiv=Ndiv,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
            end
        end

        ## Dead distributed moments initially resolved in basis b
        # ----------------------------------------------------------------------
        if hasDistributedDeadMomentsBasisb
            currentLoad = m_b_of_ζt.(ζInterp,timeNow)
            Mmax = maximum(abs.(element.parent.m_b_of_x1t(x1Interp,timeNow)))
            for (dir,ai) in zip(1:3,[a1, a2, a3])
                Mi = [currentLoad[d][dir] for d in 1:3]
                if maximum(abs.(Mi)) == 0
                    continue
                end
                plt = draw_distributed_moments!(plt, M=Mi,Mmax=Mmax,R=[R0,R0,R0],axis=dir,P2=P2,Ndiv=Ndiv,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
            end
        end

        ## Distributed follower forces initially resolved in basis A
        # ----------------------------------------------------------------------
        if hasDistributedFollowerForcesBasisA
            currentLoad = ff_A_of_ζt.(ζInterp,timeNow)
            Fmax = maximum(abs.(element.parent.ff_A_of_x1t(x1Interp,timeNow)))
            for (dir,ai) in zip(1:3,[a1, a2, a3])
                Fi = [currentLoad[d][dir] for d in 1:3]
                if maximum(abs.(Fi)) == 0
                    continue
                end
                plt = draw_distributed_forces!(plt, F=Fi,Fmax=Fmax,R=[R_n1,R_e,R_n2],ai=ai,P2=P2,Ndiv=Ndiv,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
            end
        end

        ## Distributed follower moments initially resolved in basis A
        # ----------------------------------------------------------------------
        if hasDistributedFollowerMomentsBasisA
            currentLoad = mf_A_of_ζt.(ζInterp,timeNow)
            Mmax = maximum(abs.(element.parent.mf_A_of_x1t(x1Interp,timeNow)))
            for (dir,ai) in zip(1:3,[a1, a2, a3])
                Mi = [currentLoad[d][dir] for d in 1:3]
                if maximum(abs.(Mi)) == 0
                    continue
                end
                plt = draw_distributed_moments!(plt, M=Mi,Mmax=Mmax,R=[R_n1,R_e,R_n2],axis=dir,P2=P2,Ndiv=Ndiv,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
            end
        end

        ## Distributed follower forces initially resolved in basis b
        # ----------------------------------------------------------------------
        if hasDistributedFollowerForcesBasisb
            currentLoad = ff_b_of_ζt.(ζInterp,timeNow)
            Fmax = maximum(abs.(element.parent.ff_b_of_x1t(x1Interp,timeNow)))
            for (dir,ai) in zip(1:3,[a1, a2, a3])
                Fi = [currentLoad[d][dir] for d in 1:3]
                if maximum(abs.(Fi)) == 0
                    continue
                end
                plt = draw_distributed_forces!(plt, F=Fi,Fmax=Fmax,R=[R_n1*R0_n1,R_e*R0,R_n2*R0_n2],ai=ai,P2=P2,Ndiv=Ndiv,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
            end
        end

        ## Distributed follower moments initially resolved in basis b
        # ----------------------------------------------------------------------
        if hasDistributedFollowerMomentsBasisb
            currentLoad = mf_b_of_ζt.(ζInterp,timeNow)
            Mmax = maximum(abs.(element.parent.mf_b_of_x1t(x1Interp,timeNow)))
            for (dir,ai) in zip(1:3,[a1, a2, a3])
                Mi = [currentLoad[d][dir] for d in 1:3]
                if maximum(abs.(Mi)) == 0
                    continue
                end
                plt = draw_distributed_moments!(plt, M=Mi,Mmax=Mmax,R=[R_n1*R0_n1,R_e*R0,R_n2*R0_n2],axis=dir,P2=P2,Ndiv=Ndiv,L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)
            end
        end

        ## Aerodynamic loads
        # ----------------------------------------------------------------------
        if !isnothing(aero)
            # Current values
            if problem isa DynamicProblem
                ct = problem.aeroVariablesOverTime[tIndNow][globalID].aeroCoefficients.ct
                cn = problem.aeroVariablesOverTime[tIndNow][globalID].aeroCoefficients.cn
                cm = problem.aeroVariablesOverTime[tIndNow][globalID].aeroCoefficients.cm
            else
                ct = problem.aeroVariablesOverσ[end][globalID].aeroCoefficients.ct
                cn = problem.aeroVariablesOverσ[end][globalID].aeroCoefficients.cn
                cm = problem.aeroVariablesOverσ[end][globalID].aeroCoefficients.cm
            end
            # Mormalized values of aerodynamic loads
            ctNorm = ct / maxAeroForce
            cnNorm = cn / maxAeroForce
            cmNorm = cm / maxAeroMoment
            # Draw aerodynamic loads (at element's midpoint)
            plt = draw_aero_loads!(plt,ctNorm=ctNorm,cnNorm=cnNorm,cmNorm=cmNorm,R=R_e*R0,P2=P2[2,:],L=L,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view)   
        end
    end

    return plt
end


# Draws distributed forces
function draw_distributed_forces!(plt; F,Fmax=maximum(abs.(F)),R,ai,P2,Ndiv,L,x1Plane,x2Plane,x3Plane,view,color=:purple,λ=1/10)

    # Normalized vector of distributed forces in current direction  
    FNorm = F / Fmax
    
    # Vectors for quiver
    P1 = [P2[j,:] .- L * λ * (R[j] * ai * FNorm[j]) for j in 1:Ndiv]
    P1 = hcat(P1...)'
    ΔP = P2 .- P1
    
    # Plot
    if x1Plane && isnothing(view)
        quiver!(P1[:,2], P1[:,3], quiver=(ΔP[:,2], ΔP[:,3]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x2Plane && isnothing(view)
        quiver!(P1[:,1], P1[:,3], quiver=(ΔP[:,1], ΔP[:,3]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x3Plane && isnothing(view)
        quiver!(P1[:,1], P1[:,2], quiver=(ΔP[:,1], ΔP[:,2]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
    else
        quiver!(P1[:,1], P1[:,2], P1[:,3], quiver=(ΔP[:,1], ΔP[:,2], ΔP[:,3]), color=color, lw=1, quiverhead=0.5)
    end
    
    return plt
end


# Draws distributed moments
function draw_distributed_moments!(plt; M,Mmax,R,axis,P2,Ndiv,L,x1Plane,x2Plane,x3Plane,view,color=:purple,λ=1/2)

    # Get normalized vector of distributed forces    
    MNorm = M / Mmax

    # Loop over nodes and midpoint
    for j=1:Ndiv
        # Draw moment 
        plt = draw_circular_vector!(plt,origin=P2[:,j],M=MNorm[j],R=R[j],L=L*λ, axis=axis,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view,color=color)
    end
    
    return plt
end


# Draws aerodynamic loads
function draw_aero_loads!(plt; ctNorm,cnNorm,cmNorm,R,P2,L,x1Plane,x2Plane,x3Plane,view,color=:green,λ=1/4)
        
    # Vectors for quiver
    ctDir = a2
    cnDir = a3
    P1_ct = P2 .- L * λ * ctNorm * (R * ctDir)
    P1_cn = P2 .- L * λ * cnNorm * (R * cnDir)
    ΔP_ct = P2 .- P1_ct
    ΔP_cn = P2 .- P1_cn
    
    # Plot ct and cn loads
    if x1Plane && isnothing(view)
        quiver!([P1_ct[2]], [P1_ct[3]], quiver=([ΔP_ct[2]], [ΔP_ct[3]]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
        quiver!([P1_cn[2]], [P1_cn[3]], quiver=([ΔP_cn[2]], [ΔP_cn[3]]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x2Plane && isnothing(view)
        quiver!([P1_ct[1]], [P1_ct[3]], quiver=([ΔP_ct[1]], [ΔP_ct[3]]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
        quiver!([P1_cn[1]], [P1_cn[3]], quiver=([ΔP_cn[1]], [ΔP_cn[3]]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x3Plane && isnothing(view)
        quiver!([P1_ct[1]], [P1_ct[2]], quiver=([ΔP_ct[1]], [ΔP_ct[2]]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
        quiver!([P1_cn[1]], [P1_cn[2]], quiver=([ΔP_cn[1]], [ΔP_cn[2]]), color=color, lw=1, quiverhead=0.5, aspect_ratio=:equal)
    else
        quiver!([P1_ct[1]], [P1_ct[2]], [P1_ct[3]], quiver=([ΔP_ct[1]], [ΔP_ct[2]], [ΔP_ct[3]]), color=color, lw=1, quiverhead=0.5)
        quiver!([P1_cn[1]], [P1_cn[2]], [P1_cn[3]], quiver=([ΔP_cn[1]], [ΔP_cn[2]], [ΔP_cn[3]]), color=color, lw=1, quiverhead=0.5)
    end

    # Plot cm load
    plt = draw_circular_vector!(plt,origin=P2,M=cmNorm,R=R,L=L, axis=1,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view,color=color)
    
    return plt
end


# Draws aerodynamic loads
function draw_circular_vector!(plt; origin,M,R,L,axis,x1Plane,x2Plane,x3Plane,view,color=:red,r=L*abs(M)/20,angle=9π/5,ah=0.3*r,divisions=30)
    
    # Circle definition
    circFun = if axis == 1
        θ -> [0; r*cos(θ); r*sin(θ)]
    elseif axis == 2
        θ -> [r*cos(θ); 0; r*sin(θ)]
    elseif axis == 3
        θ -> [r*cos(θ); r*sin(θ); 0]
    else
        error("Choose axis as 1, 2 or 3")
    end
    θRange = range(0,stop=angle,length=divisions)
    f = hcat(circFun.(θRange)...)
    circumference = R * f .+ origin

    # Arrow parameters
    dir = sign(M)    
    if axis == 1
        tt = dir == 1 ? angle : 0
        arrow = [0       0  0;
                 -ah/2    0  ah/2;
                 -dir*ah  0  -dir*ah]
        Rtt = [1       0        0;
               0  cos(tt) -sin(tt);
               0 sin(tt)  cos(tt)]
    elseif axis == 2
        tt = dir == 1 ? 0 : angle
        arrow = [ah/2 0 -ah/2;
                 0 0 0;
                 dir*ah  0 dir*ah]
        Rtt = [cos(tt) 0 -sin(tt);
               0       1        0;
               sin(tt) 0 cos(tt)]
    elseif axis == 3
        tt = dir == 1 ? angle : 0
        arrow = [-ah/2 0 ah/2;
                 -dir*ah  0 -dir*ah;
                 0 0 0]
        Rtt = [cos(tt) -sin(tt) 0;
               sin(tt)  cos(tt) 0;
               0        0       1]
    end
    arrowhead = R * (circFun(tt) .+ Rtt * arrow) .+ origin
    
    # Plot circumference and arrowhead
    if x1Plane && isnothing(view)
        plot!(circumference[2,:], circumference[3,:], color=color, lw=1, aspect_ratio=:equal, label=false)
        plot!(arrowhead[2,:], arrowhead[3,:], color=color, lw=1, aspect_ratio=:equal, label=false)
    elseif x2Plane && isnothing(view)
        plot!(circumference[1,:], circumference[3,:], color=color, lw=1, aspect_ratio=:equal, label=false)
        plot!(arrowhead[1,:], arrowhead[3,:], color=color, lw=1, aspect_ratio=:equal, label=false)
    elseif x3Plane && isnothing(view)
        plot!(circumference[1,:], circumference[2,:], color=color, lw=1, aspect_ratio=:equal, label=false)
        plot!(arrowhead[1,:], arrowhead[2,:], color=color, lw=1, aspect_ratio=:equal, label=false)
    else
        plot!(circumference[1,:], circumference[2,:], circumference[3,:], color=color, lw=1, label=false)
        plot!(arrowhead[1,:], arrowhead[2,:], arrowhead[3,:], color=color, lw=1, label=false)
    end

    return plt
end