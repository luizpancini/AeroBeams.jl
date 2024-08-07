"""
plot_undeformed_assembly(model::Model)

Plots the nodal coordinates of the assembly of beams

# Arguments
- `model::Model`
"""
function plot_undeformed_assembly(model::Model,view=(45,45))

    # Initialize backend
    pyplot()

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
            plot!([r_n[i][1], r_n[i+1][1]], [r_n[i][2], r_n[i+1][2]], [r_n[i][3], r_n[i+1][3]], c=:black, linewidth=2, label=false)
        end

    end

    display(plt)

    return plt
end
export plot_undeformed_assembly


"""
plot_steady_deformation(problem::Problem)

Plots the initial and final deformed states for the model in the given problem

# Arguments
- `problem::Problem`
- `view` = view angles
"""
function plot_steady_deformation(problem::Problem; plotBCs::Bool=true,view::Union{Nothing,Tuple{Int64,Int64}}=nothing,scale::Number=1,linewidth::Number=1,colorUndef=:black,colorDef=:blue,grid::Bool=true,tolPlane::Number=1e-8,plotAeroSurf::Bool=true,surfα::Float64=0.5,save::Bool=false,savePath::String="/test/outputs/figures/fig.pdf")

    # Validate
    @assert typeof(problem) in [SteadyProblem,TrimProblem,EigenProblem]
    @assert scale > 0
    @assert linewidth > 0
    @assert tolPlane > 0
    @assert 0 < surfα <= 1

    # Set backend
    pyplot()

    # Unpack
    @unpack elements,nElementsTotal,units = problem.model

    # Initialize plot
    plt = plot(grid=grid)
    plot!([NaN],[NaN], c=colorUndef, linewidth=linewidth, label="Undeformed")
    plot!([NaN],[NaN], c=colorDef, linewidth=linewidth, label="Deformed")

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
        x1Def[e] = x1Undef[e] .+ scale*[u_n1[1]; u_n2[1]]
        x2Def[e] = x2Undef[e] .+ scale*[u_n1[2]; u_n2[2]]
        x3Def[e] = x3Undef[e] .+ scale*[u_n1[3]; u_n2[3]]
        # Skip elements without aero surfaces, or if those are not to be plotted
        if !plotAeroSurf || isnothing(aero)
            undefAirfoilCoords_n1[e],undefAirfoilCoords_n2[e],defAirfoilCoords_n1[e],defAirfoilCoords_n2[e] = nothing,nothing,nothing,nothing
            continue
        end
        @unpack Rw = aero
        # Rotation tensors
        R_n1,_ = rotation_tensor_WM(p_n1)
        R_n2,_ = rotation_tensor_WM(p_n2)
        # Undeformed nodal airfoil coordinates
        undefAirfoilCoords_n1[e],undefAirfoilCoords_n2[e] = get_undeformed_airfoil_coords(element)
        # Deformed nodal airfoil coordinates ( bring to origin (-r) and resolve in basis A (R_n*Rw*), then throw back to initial position (+r) and add scaled displacements (+scale*u))
        defAirfoilCoords_n1[e] = R_n1*Rw*(undefAirfoilCoords_n1[e] .- r_n1) .+ r_n1 .+ scale*u_n1    
        defAirfoilCoords_n2[e] = R_n2*Rw*(undefAirfoilCoords_n2[e] .- r_n2) .+ r_n2 .+ scale*u_n2
    end

    # Set TFs for plane views
    x1Plane = maximum(abs.(vcat(x1Def...))) < tolPlane ? true : false
    x2Plane = maximum(abs.(vcat(x2Def...))) < tolPlane ? true : false
    x3Plane = maximum(abs.(vcat(x3Def...))) < tolPlane ? true : false

    # Plot undeformed and deformed beam assemblies actording to view
    if isnothing(view)
        if x1Plane
            plot!(xlabel=string("\$x_2\$ [",units.length,"]"),ylabel=string("\$x_3\$ [",units.length,"]"))
            plot!(x2Undef, x3Undef, color=colorUndef,linewidth=linewidth,label=false)
            plot!(x2Def, x3Def, color=colorDef,aspect_ratio=:equal,linewidth=linewidth,label=false)
        elseif x2Plane
            plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_3\$ [",units.length,"]"))
            plot!(x1Undef, x3Undef, color=colorUndef,linewidth=linewidth,label=false)
            plot!(x1Def, x3Def, color=colorDef,aspect_ratio=:equal,linewidth=linewidth,label=false)
        elseif x3Plane
            plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"))
            plot!(x1Undef, x2Undef, color=colorUndef,linewidth=linewidth,label=false)
            plot!(x1Def, x2Def, color=colorDef,aspect_ratio=:equal,linewidth=linewidth,label=false)
        else
            view = (45,45)
            plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"),zlabel=string("\$x_3\$ [",units.length,"]"))
            plot!(x1Undef, x2Undef, x3Undef, camera=view,color=colorUndef,linewidth=linewidth,label=false)
            plot!(x1Def, x2Def, x3Def, camera=view,color=colorDef,linewidth=linewidth,label=false)
        end
    else
        plot!(xlabel=string("\$x_1\$ [",units.length,"]"),ylabel=string("\$x_2\$ [",units.length,"]"),zlabel=string("\$x_3\$ [",units.length,"]"))
        plot!(x1Undef, x2Undef, x3Undef, camera=view,color=colorUndef,linewidth=linewidth,label=false)
        plot!(x1Def, x2Def, x3Def, camera=view,color=colorDef,linewidth=linewidth,label=false)
    end

    # Plot aerodynamic surfaces, if applicable
    for e in 1:nElementsTotal
        if isnothing(defAirfoilCoords_n1[e])
            continue
        end
        # Plot undeformed aerodynamic surfaces
        Xu = [undefAirfoilCoords_n1[e][1,:] undefAirfoilCoords_n2[e][1,:]]
        Yu = [undefAirfoilCoords_n1[e][2,:] undefAirfoilCoords_n2[e][2,:]]
        Zu = [undefAirfoilCoords_n1[e][3,:] undefAirfoilCoords_n2[e][3,:]]
        surface!(Xu,Yu,Zu, color=palette([colorUndef,colorUndef],2), colorbar=false, alpha=surfα)
        # Plot deformed aerodynamic surfaces
        Xd = [defAirfoilCoords_n1[e][1,:] defAirfoilCoords_n2[e][1,:]]
        Yd = [defAirfoilCoords_n1[e][2,:] defAirfoilCoords_n2[e][2,:]]
        Zd = [defAirfoilCoords_n1[e][3,:] defAirfoilCoords_n2[e][3,:]]
        surface!(Xd,Yd,Zd, color=palette([colorDef,colorDef],2), colorbar=false, alpha=surfα)
    end

    # Plot BCs and distributed loads, if applicable
    if plotBCs
        # Plot generalized displacements and concentrated loads
        plt = plot_BCs!(plt,problem,x1Def,x2Def,x3Def,x1Plane,x2Plane,x3Plane,view,0,1)
        # Compute maximum aero loads
        maxAeroForce,maxAeroMoment = maximum_aero_loads(problem)
        # Plot distributed loads (including aerodynamic)
        plt = plot_distributed_loads!(plt,problem,x1Def,x2Def,x3Def,x1Plane,x2Plane,x3Plane,view,maxAeroForce,maxAeroMoment,1)
    end

    # Set the same plot extrema across all axis (equal aspect ratio), if applicable
    if !(x1Plane || x2Plane || x3Plane) || !isnothing(view)
        plotMax = max(maximum(reduce(vcat, x1Def)),maximum(reduce(vcat, x2Def)),maximum(reduce(vcat, x3Def)),maximum(reduce(vcat, x1Undef)),maximum(reduce(vcat, x2Undef)),maximum(reduce(vcat, x3Undef)))
        plotMin = min(minimum(reduce(vcat, x1Def)),minimum(reduce(vcat, x2Def)),minimum(reduce(vcat, x3Def)),minimum(reduce(vcat, x1Undef)),minimum(reduce(vcat, x2Undef)),minimum(reduce(vcat, x3Undef)))
        plot!(xlims=[plotMin,plotMax],ylims=[plotMin,plotMax],zlims=[plotMin,plotMax])
    end

    # Save, if applicable
    if save
        savefig(string(pwd(),savePath))
    end

    return plt
end
export plot_steady_deformation


"""
plot_steady_outputs(problem::Problem)

Plots outputs of a steady problem

# Arguments
- `problem::Problem`
"""
function plot_steady_outputs(problem::Problem; outputs::Vector{String}=["u","p","F","M","V","Ω","α","cn","cm","ct","cl","cd"],beamGroups=1:length(problem.model.beams),lw::Number=1,save::Bool=false,saveFolder::String="/test/outputs/figures/",saveExtension::String=".pdf")

    # Validate
    validOutputs = ["u","u1","u2","u3","p","p1","p2","p3","F","F1","F2","F3","M","M1","M2","M3","V","V1","V2","V3","Ω","Ω1","Ω2","Ω3","α","cn","cm","ct","cl","cd"]
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
        # Add to arrays
        append!(x1Nodes[beamID],x1_n)
        append!(x1Elems[beamID],x1_e)
        append!(uNodes[beamID],u)
        append!(pNodes[beamID],p)
        append!(FNodes[beamID],F)
        append!(MNodes[beamID],M)
        append!(VElems[beamID],V)
        append!(ΩElems[beamID],Ω)
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
        plt_u1 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=uNodes,ind=1,units=units,YLabel="u_1",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"u1",saveExtension))
        end
    end
    if "u" in outputs || "u2" in outputs
        plt_u2 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=uNodes,ind=2,units=units,YLabel="u_2",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"u2",saveExtension))
        end
    end
    if "u" in outputs || "u3" in outputs
        plt_u3 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=uNodes,ind=3,units=units,YLabel="u_3",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"u3",saveExtension))
        end
    end
    if "p" in outputs || "p1" in outputs
        plt_p1 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=pNodes,ind=1,units=units,YLabel="p_1",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"p1",saveExtension))
        end
    end
    if "p" in outputs || "p2" in outputs
        plt_p2 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=pNodes,ind=2,units=units,YLabel="p_2",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"p2",saveExtension))
        end
    end 
    if "p" in outputs || "p3" in outputs
        plt_p3 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=pNodes,ind=3,units=units,YLabel="p_3",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"p3",saveExtension))
        end
    end
    if "F" in outputs || "F1" in outputs
        plt_F1 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=FNodes,ind=1,units=units,YLabel="F_1^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"F1",saveExtension))
        end
    end
    if "F" in outputs || "F2" in outputs
        plt_F2 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=FNodes,ind=2,units=units,YLabel="F_2^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"F2",saveExtension))
        end
    end
    if "F" in outputs || "F3" in outputs
        plt_F3 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=FNodes,ind=3,units=units,YLabel="F_3^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"F3",saveExtension))
        end
    end
    if "M" in outputs || "M1" in outputs
        plt_M1 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=MNodes,ind=1,units=units,YLabel="M_1^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"M1",saveExtension))
        end
    end
    if "M" in outputs || "M2" in outputs
        plt_M2 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=MNodes,ind=2,units=units,YLabel="M_2^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"M2",saveExtension))
        end
    end
    if "M" in outputs || "M3" in outputs
        plt_M3 = plot_output_of_x1(beamGroups,x1=x1Nodes,output=MNodes,ind=3,units=units,YLabel="M_3^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"M3",saveExtension))
        end
    end
    if "V" in outputs || "V1" in outputs
        plt_V1 = plot_output_of_x1(beamGroups,x1=x1Elems,output=VElems,ind=1,units=units,YLabel="V_1^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"V1",saveExtension))
        end
    end
    if "V" in outputs || "V2" in outputs
        plt_V2 = plot_output_of_x1(beamGroups,x1=x1Elems,output=VElems,ind=2,units=units,YLabel="V_2^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"V2",saveExtension))
        end
    end
    if "V" in outputs || "V3" in outputs
        plt_V3 = plot_output_of_x1(beamGroups,x1=x1Elems,output=VElems,ind=3,units=units,YLabel="V_3^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"V3",saveExtension))
        end
    end
    if "Ω" in outputs || "Ω1" in outputs
        plt_Ω1 = plot_output_of_x1(beamGroups,x1=x1Elems,output=ΩElems,ind=1,units=units,YLabel="Ω_1^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"Ω1",saveExtension))
        end
    end
    if "Ω" in outputs || "Ω2" in outputs
        plt_Ω2 = plot_output_of_x1(beamGroups,x1=x1Elems,output=ΩElems,ind=2,units=units,YLabel="Ω_2^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"Ω2",saveExtension))
        end
    end
    if "Ω" in outputs || "Ω3" in outputs
        plt_Ω3 = plot_output_of_x1(beamGroups,x1=x1Elems,output=ΩElems,ind=3,units=units,YLabel="Ω_3^*",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"Ω3",saveExtension))
        end
    end
    if "α" in outputs 
        plt_α = plot_output_of_x1(beamGroups,x1=x1Elems,output=αₑElems,ind=0,units=units,YLabel="\\alpha",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"α",saveExtension))
        end
    end
    if "cn" in outputs 
        plt_cn = plot_output_of_x1(beamGroups,x1=x1Elems,output=cnElems,ind=0,units=units,YLabel="c_n",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"cn",saveExtension))
        end
    end
    if "cm" in outputs 
        plt_cm = plot_output_of_x1(beamGroups,x1=x1Elems,output=cmElems,ind=0,units=units,YLabel="c_m",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"cm",saveExtension))
        end
    end
    if "ct" in outputs 
        plt_ct = plot_output_of_x1(beamGroups,x1=x1Elems,output=ctElems,ind=0,units=units,YLabel="c_t",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"ct",saveExtension))
        end
    end
    if "cl" in outputs 
        plt_cl = plot_output_of_x1(beamGroups,x1=x1Elems,output=clElems,ind=0,units=units,YLabel="c_l",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"cl",saveExtension))
        end
    end
    if "cd" in outputs 
        plt_cd = plot_output_of_x1(beamGroups,x1=x1Elems,output=cdElems,ind=0,units=units,YLabel="c_d",lw=lw)
        if save
            savefig(string(pwd(),saveFolder,"cd",saveExtension))
        end
    end

    return nothing
end
export plot_steady_outputs


"""
plot_output_of_x1(beamGroups,x1,output,ind,colors=get(colorschemes[:rainbow], LinRange(0,1,length(beamGroups))))

Plots output along the arclength coordinate for each beam group

# Arguments
- `problem::Problem`
"""
function plot_output_of_x1(beamGroups; x1,output,ind,units,YLabel,colors=get(colorschemes[:rainbow],LinRange(0,1,length(beamGroups)+1)),lw=1)

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
    elseif occursin("Ω",YLabel)
        γ = units.angle == "deg" ? 180/π : 1
        yLabelUnit = string(units.angle,"/s") 
    end

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
        plt = plot!(plt, x1BeamGroup,outputBeamGroup,lw=lw,color=colors[i],label="Beam group $(i)")
    end

    # Display
    display(plt)

    return plt
end


"""
get_undeformed_airfoil_coords(element::Element)

Computes the undeformed nodal airfoil coordinates

# Arguments
- `element::Element`
"""
function get_undeformed_airfoil_coords(element::Element)

    @unpack x1_n1,x1_n2,r_n1,r_n2,R0_n1,R0_n2 = element
    @unpack c,normSparPos = element.parent.aeroSurface
    @unpack airfoil = element.aero

    # Airfoil coordinates in the X-plane
    Y,Z = airfoil.coordinates[:,1],airfoil.coordinates[:,2]

    # Number of points
    N = length(Y)

    # Nodal chords and normalized spar positions
    if c isa Float64
        c1 = c2 = c
    else
        c1,c2 = c(x1_n1),c(x1_n2)
    end
    if normSparPos isa Float64
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
    XYZRot1,XYZRot2 = R0_n1*[zeros(N)'; Y1'; Z1'], R0_n2*[zeros(N)'; Y2'; Z2']           
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


"""
plot_BCs!(plt,problem::Problem,element::Element)

Plots all boundary conditions at the current time

# Arguments
- `element::Element`
"""
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
        plt = draw_BC!(plt,isLoad=isLoad,deadLoadsOnA=deadLoadsOnA,followerLoadsOnA=followerLoadsOnA,deadLoadsOnb=deadLoadsOnb,followerLoadsOnb=followerLoadsOnb,R0_n=R0_n,R_n=R_n,Fmax=Fmax,Mmax=Mmax,L=L,P2=P2,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view,tIndNow=tIndNow)
    end

    return plt

end


"""
draw_BC!(plt,problem::Problem,element::Element)

Draws boundary condition

# Arguments
- 
"""
function draw_BC!(plt; isLoad,deadLoadsOnA,followerLoadsOnA,deadLoadsOnb,followerLoadsOnb,R0_n,R_n,Fmax,Mmax,L,P2,x1Plane,x2Plane,x3Plane,view,tIndNow)

    # Loop DOFs
    for DOF in 1:6
        # Draw generalized displacement BC
        if !isLoad[DOF] && tIndNow == 1
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


"""
draw_concentrated_force!(plt,problem::Problem,element::Element)

Draws concentrated force

# Arguments
- 
"""
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
        plt = quiver!(plt, [P1[2]], [P1[3]], quiver=([ΔP[2]], [ΔP[3]]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x2Plane && isnothing(view)
        plt = quiver!(plt, [P1[1]], [P1[3]], quiver=([ΔP[1]], [ΔP[3]]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x3Plane && isnothing(view)
        plt = quiver!(plt, [P1[1]], [P1[2]], quiver=([ΔP[1]], [ΔP[2]]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    else
        plt = quiver!(plt, [P1[1]], [P1[2]], [P1[3]], quiver=([ΔP[1]], [ΔP[2]], [ΔP[3]]), color=color, linewidth=1, quiverhead=0.5)
    end

    return plt
end


"""
draw_concentrated_moment!(plt,problem::Problem,element::Element)

Draws concentrated moment

# Arguments
- 
"""
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


"""
draw_generalized_displacement!(plt,problem::Problem,element::Element)

Draws generalized displacement boundary condition

# Arguments
- 
"""
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
        plt = plot!(plt, x, y, seriestype=:shape, fillcolor=color, linecolor=color, label=false)
    else
        plt = surface!(plt, P2[1] + π*λ*L*Δ[1] .+ x1Cone, P2[2] + π*λ*L*Δ[2] .+ x2Cone, P2[3] + π*λ*L*Δ[3] .+ x3Cone, color=palette([color,color],2), colorbar=false)
    end

    return plt
end


"""
plot_distributed_loads!(plt,problem::Problem,element::Element)

Plots all distributed loads at the current time

# Arguments
- 
"""
function plot_distributed_loads!(plt,problem::Problem,x1Def,x2Def,x3Def,x1Plane,x2Plane,x3Plane,view,maxAeroForce,maxAeroMoment,tIndNow)

    @unpack timeNow = problem
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


"""
draw_distributed_forces!(element::Element)

Draws distributed forces

# Arguments
- element::Element
"""
function draw_distributed_forces!(plt; F,Fmax=maximum(abs.(F)),R,ai,P2,Ndiv,L,x1Plane,x2Plane,x3Plane,view,color=:purple,λ=1/10)

    # Normalized vector of distributed forces in current direction  
    FNorm = F / Fmax
    
    # Vectors for quiver
    P1 = [P2[j,:] .- L * λ * (R[j] * ai * FNorm[j]) for j in 1:Ndiv]
    P1 = hcat(P1...)'
    ΔP = P2 .- P1
    
    # Plot
    if x1Plane && isnothing(view)
        plt = quiver!(plt, P1[:,2], P1[:,3], quiver=(ΔP[:,2], ΔP[:,3]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x2Plane && isnothing(view)
        plt = quiver!(plt, P1[:,1], P1[:,3], quiver=(ΔP[:,1], ΔP[:,3]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x3Plane && isnothing(view)
        plt = quiver!(plt, P1[:,1], P1[:,2], quiver=(ΔP[:,1], ΔP[:,2]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    else
        plt = quiver!(plt, P1[:,1], P1[:,2], P1[:,3], quiver=(ΔP[:,1], ΔP[:,2], ΔP[:,3]), color=color, linewidth=1, quiverhead=0.5)
    end
    
    return plt
end


"""
draw_distributed_moments!(element::Element)

Draws distributed moments

# Arguments
- element::Element
"""
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


"""
draw_aero_loads!(element::Element)

Draws aerodynamic loads

# Arguments
- element::Element
"""
function draw_aero_loads!(plt; ctNorm,cnNorm,cmNorm,R,P2,L,x1Plane,x2Plane,x3Plane,view,color=:green,λ=1/4)
        
    # Vectors for quiver
    ctDir = [0; 1; 0]
    cnDir = [0; 0; 1]
    P1_ct = P2 .- L * λ * ctNorm * (R * ctDir)
    P1_cn = P2 .- L * λ * cnNorm * (R * cnDir)
    ΔP_ct = P2 .- P1_ct
    ΔP_cn = P2 .- P1_cn
    
    # Plot ct load
    if x1Plane && isnothing(view)
        plt = quiver!(plt, [P1_ct[2]], [P1_ct[3]], quiver=([ΔP_ct[2]], [ΔP_ct[3]]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x2Plane && isnothing(view)
        plt = quiver!(plt, [P1_ct[1]], [P1_ct[3]], quiver=([ΔP_ct[1]], [ΔP_ct[3]]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x3Plane && isnothing(view)
        plt = quiver!(plt, [P1_ct[1]], [P1_ct[2]], quiver=([ΔP_ct[1]], [ΔP_ct[2]]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    else
        plt = quiver!(plt, [P1_ct[1]], [P1_ct[2]], [P1_ct[3]], quiver=([ΔP_ct[1]], [ΔP_ct[2]], [ΔP_ct[3]]), color=color, linewidth=1, quiverhead=0.5)
    end

    # Plot cn load
    if x1Plane && isnothing(view)
        plt = quiver!(plt, [P1_cn[2]], [P1_cn[3]], quiver=([ΔP_cn[2]], [ΔP_cn[3]]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x2Plane && isnothing(view)
        plt = quiver!(plt, [P1_cn[1]], [P1_cn[3]], quiver=([ΔP_cn[1]], [ΔP_cn[3]]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    elseif x3Plane && isnothing(view)
        plt = quiver!(plt, [P1_cn[1]], [P1_cn[2]], quiver=([ΔP_cn[1]], [ΔP_cn[2]]), color=color, linewidth=1, quiverhead=0.5, aspect_ratio=:equal)
    else
        plt = quiver!(plt, [P1_cn[1]], [P1_cn[2]], [P1_cn[3]], quiver=([ΔP_cn[1]], [ΔP_cn[2]], [ΔP_cn[3]]), color=color, linewidth=1, quiverhead=0.5)
    end

    # Plot cm load
    plt = draw_circular_vector!(plt,origin=P2,M=cmNorm,R=R,L=L, axis=1,x1Plane=x1Plane,x2Plane=x2Plane,x3Plane=x3Plane,view=view,color=color)
    
    return plt
end


"""
draw_aero_loads!(element::Element)

Draws aerodynamic loads

# Arguments
- element::Element
"""
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
    
    # Plot circumference
    if x1Plane && isnothing(view)
        plt = plot!(plt, circumference[2,:], circumference[3,:], color=color, linewidth=1, aspect_ratio=:equal, label=false)
    elseif x2Plane && isnothing(view)
        plt = plot!(plt, circumference[1,:], circumference[3,:], color=color, linewidth=1, aspect_ratio=:equal, label=false)
    elseif x3Plane && isnothing(view)
        plt = plot!(plt, circumference[1,:], circumference[2,:], color=color, linewidth=1, aspect_ratio=:equal, label=false)
    else
        plt = plot!(plt, circumference[1,:], circumference[2,:], circumference[3,:], color=color, linewidth=1, label=false)
    end

    # Plot arrowhead
    if x1Plane && isnothing(view)
        plt = plot!(plt, arrowhead[2,:], arrowhead[3,:], color=color, linewidth=1, aspect_ratio=:equal, label=false)
    elseif x2Plane && isnothing(view)
        plt = plot!(plt, arrowhead[1,:], arrowhead[3,:], color=color, linewidth=1, aspect_ratio=:equal, label=false)
    elseif x3Plane && isnothing(view)
        plt = plot!(plt, arrowhead[1,:], arrowhead[2,:], color=color, linewidth=1, aspect_ratio=:equal, label=false)
    else
        plt = plot!(plt, arrowhead[1,:], arrowhead[2,:], arrowhead[3,:], color=color, linewidth=1, label=false)
    end

    return plt
end


"""
maximum_aero_loads(problem::Problem)

Computes the maximum absolute value of aerodynamic loads (over time or load factor)

# Arguments
- problem::Problem
"""
function maximum_aero_loads(problem::Problem)

    # Compute maximum aerodynamic loads over time or load factor
    # --------------------------------------------------------------------------
    # Get aerodynamic variables of all elements
    elementsAeroVariables = problem isa DynamicProblem ? [problem.aeroVariablesOverTime[i] for i in 1:length(problem.t)] : problem.aeroVariablesOverσ[end]
    # Filter out Nothings (from elements without aero)
    elementsAeroVariablesFiltered = filter(x -> x !== nothing, elementsAeroVariables)
    
    if !isempty(elementsAeroVariablesFiltered)
        # Maxima
        ctMax = maximum(vcat(aero.aeroCoefficients.ct for aero in elementsAeroVariablesFiltered)...)
        cnMax = maximum(vcat(aero.aeroCoefficients.cn for aero in elementsAeroVariablesFiltered)...)
        cmMax = maximum(vcat(aero.aeroCoefficients.cm for aero in elementsAeroVariablesFiltered)...)
        # Maximum force and moment coefficients
        maxAeroForce = max(ctMax,cnMax)
        maxAeroMoment = cmMax
    else
        maxAeroForce = maxAeroMoment = 0
    end

    return maxAeroForce,maxAeroMoment
end