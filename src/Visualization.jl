using Plots, ColorSchemes, Measures

## OPTIONS, BACKEND AND SAVING
# ------------------------------------------------------------------------------

# Validates the options that are common to the plotting functions
function _check_options(backend,DPI,plotLimits=nothing,surfα=0.5,plotSize=(600,400); positives...)
    @assert backend in [:gr, :plotlyjs] "backend must be either :gr or :plotlyjs"
    @assert 100 <= DPI <= 600 "set DPI between 100 and 600"
    @assert 0 < surfα <= 1 "set surfα in the interval (0,1]"
    for (name,value) in positives
        @assert value > 0 "set $(name) as a positive number"
    end
    @assert length(plotSize) == 2 && all(>(0),plotSize) "set plotSize as a pair of positive numbers (width,height) in pixels"
    if !isnothing(plotLimits)
        @assert all(plotLimits[a][1] < plotLimits[a][2] for a in 1:3) "set each of the plotLimits as an ascending pair"
    end
    return nothing
end

# Sets the plotting backend
_set_backend(backend) = getfield(Plots, backend)()

# Full path of a file, given its path relative to the package's directory (with or without a leading separator)
_full_path(relPath) = joinpath(pkgdir(AeroBeams), lstrip(relPath,['/','\\']))

# Saves a figure at a path relative to the package's directory
function _save_figure(fig, relPath)
    fullPath = _full_path(relPath)
    mkpath(dirname(fullPath))
    savefig(fig, fullPath)
    return nothing
end


## BASIC TYPES AND CONSTANTS
# ------------------------------------------------------------------------------

const _I3 = SMatrix{3,3,Float64,9}(1,0,0,0,1,0,0,0,1)
const _O3 = SVector{3,Float64}(0,0,0)
const _NAN3 = SVector{3,Float64}(NaN,NaN,NaN)
_unit(i) = SVector{3,Float64}(i == 1, i == 2, i == 3)
_nanmax(a, xs...) = foldl((m,x) -> isnan(x) ? m : max(m,x), xs; init=a)

# Rigid-body transformation x ↦ u + R*x, from basis A to the basis of the observer
struct _Xform
    u::SVector{3,Float64}
    R::SMatrix{3,3,Float64,9}
end
_Xform() = _Xform(_O3,_I3)
(T::_Xform)(x) = T.u + T.R*x

# Displacements and rotation tensors at the nodes and midpoint of an element (all scaled by the plot's scale)
struct _Pose
    u1::SVector{3,Float64}
    u2::SVector{3,Float64}
    R1::SMatrix{3,3,Float64,9}
    R2::SMatrix{3,3,Float64,9}
    Rm::SMatrix{3,3,Float64,9}
end
_pose(u1,u2,p1,p2,pm) = _Pose(SVector{3,Float64}(u1), SVector{3,Float64}(u2), rotation_tensor_WM(SVector{3,Float64}(p1)), rotation_tensor_WM(SVector{3,Float64}(p2)), rotation_tensor_WM(SVector{3,Float64}(pm)))
const _NO_POSE = _Pose(_O3,_O3,_I3,_I3,_I3)

# Undeformed geometry of a model: element end positions and airfoil sections, and the reference length for loads (the largest beam length, limited by the extent of the assembly)
struct _Geometry
    model::Model
    beams::Vector{Vector{Int}}
    r1::Vector{SVector{3,Float64}}
    r2::Vector{SVector{3,Float64}}
    S1::Vector{Union{Nothing,Matrix{Float64}}}
    S2::Vector{Union{Nothing,Matrix{Float64}}}
    Lref::Float64
    anchor::Tuple{Int,Int}
end

function _geometry(model::Model; surfaces::Bool=true)
    elements = model.assemblyData.elements
    r1 = [SVector{3,Float64}(element.assemblyData.r_n1) for element in elements]
    r2 = [SVector{3,Float64}(element.assemblyData.r_n2) for element in elements]
    S1 = Vector{Union{Nothing,Matrix{Float64}}}(nothing,length(elements))
    S2 = Vector{Union{Nothing,Matrix{Float64}}}(nothing,length(elements))
    if surfaces
        for (e,element) in enumerate(elements)
            isnothing(element.aero.surface) && continue
            S1[e],S2[e] = get_undeformed_airfoil_coords(element)
        end
    end
    points = vcat(r1,r2)
    extent = maximum(maximum(p[a] for p in points) - minimum(p[a] for p in points) for a in 1:3)
    ℓmax = maximum(beam.geometry.ℓ for beam in model.beams)
    Lref = extent > 0 ? min(extent,ℓmax) : ℓmax
    # Anchor: the element end (element, node) closest to the origin of basis A, i.e. the point that initially coincides with it
    ends = [(e,n) for e in eachindex(r1) for n in 1:2]
    anchor = ends[argmin([norm(n == 1 ? r1[e] : r2[e]) for (e,n) in ends])]
    return _Geometry(model, [copy(beam.assemblyData.elementsGlobalIDs) for beam in model.beams], r1, r2, S1, S2, Lref, anchor)
end

# A configuration (undeformed or deformed) of the assembly, as seen by the observer
struct _Config
    geom::_Geometry
    poses::Vector{_Pose}
    T::_Xform
    X1::Vector{SVector{3,Float64}}
    X2::Vector{SVector{3,Float64}}
    S1::Vector{Union{Nothing,Matrix{Float64}}}
    S2::Vector{Union{Nothing,Matrix{Float64}}}
end

# Current position of the anchor of the assembly, as seen by the observer
function _anchor_position(cfg::_Config)
    e,n = cfg.geom.anchor
    return n == 1 ? cfg.X1[e] : cfg.X2[e]
end

_undeformed_config(geom::_Geometry) = _Config(geom, fill(_NO_POSE,length(geom.r1)), _Xform(), geom.r1, geom.r2, geom.S1, geom.S2)

function _configuration(geom::_Geometry, poses::Vector{_Pose}, T::_Xform=_Xform())
    n = length(poses)
    X1 = Vector{SVector{3,Float64}}(undef,n)
    X2 = Vector{SVector{3,Float64}}(undef,n)
    S1 = Vector{Union{Nothing,Matrix{Float64}}}(nothing,n)
    S2 = Vector{Union{Nothing,Matrix{Float64}}}(nothing,n)
    for e in 1:n
        pose = poses[e]
        X1[e] = T(geom.r1[e] + pose.u1)
        X2[e] = T(geom.r2[e] + pose.u2)
        isnothing(geom.S1[e]) && continue
        S1[e] = _move_section(geom.S1[e], geom.r1[e], pose.u1, pose.R1, T)
        S2[e] = _move_section(geom.S2[e], geom.r2[e], pose.u2, pose.R2, T)
    end
    return _Config(geom,poses,T,X1,X2,S1,S2)
end

# Moves the undeformed section points S (3×N) of the node at r to its deformed position r+u and rotation R
function _move_section(S, r, u, R, T)
    M = T.R*R
    return M*S .+ (T.u + T.R*(r + u - R*r))
end

# Displacements and rotation parameters of an element state (or modal state) at its nodes and midpoint: u1, u2, p1, p2 and pm
_state_values(s) = (SVector{3,Float64}(s.nodal.u1), SVector{3,Float64}(s.nodal.u2), SVector{3,Float64}(s.nodal.p1), SVector{3,Float64}(s.nodal.p2), SVector{3,Float64}(s.midpoint.p))

# Poses of all elements from their states, with the displacements and rotation parameters scaled (and superposed to those of a mode shape, times `factor`, if given). The assembly is translated so that the midpoint of the element `centralize` stays at its undeformed position
function _poses(states, scale; modal=nothing, factor=0.0, centralize=nothing)
    Δu = isnothing(centralize) ? _O3 : -scale*SVector{3,Float64}(states[centralize].midpoint.u)
    poses = Vector{_Pose}(undef,length(states))
    for e in eachindex(states)
        x = scale .* _state_values(states[e])
        isnothing(modal) || (x = x .+ factor .* _state_values(modal[e]))
        poses[e] = _pose(x[1]+Δu, x[2]+Δu, x[3], x[4], x[5])
    end
    return poses
end


## VIEW, EXTENT AND LIMITS
# ------------------------------------------------------------------------------

# Bounding box of the points that are plotted
struct _Extent
    lo::Vector{Float64}
    hi::Vector{Float64}
end
_Extent() = _Extent(fill(Inf,3),fill(-Inf,3))

# Expands the bounding box with a point, the points of a matrix of sections (3×N) or a configuration, all relative to `shift`
function _expand!(ext::_Extent, x::AbstractVector, shift=_O3)
    any(isnan,x) && return ext
    for a in 1:3
        ext.lo[a] = min(ext.lo[a],x[a]-shift[a])
        ext.hi[a] = max(ext.hi[a],x[a]-shift[a])
    end
    return ext
end

function _expand!(ext::_Extent, S::AbstractMatrix, shift=_O3)
    any(isnan,S) && return ext
    for a in 1:3
        ext.lo[a] = min(ext.lo[a],minimum(@view S[a,:])-shift[a])
        ext.hi[a] = max(ext.hi[a],maximum(@view S[a,:])-shift[a])
    end
    return ext
end

function _expand!(ext::_Extent, cfg::_Config, shift=_O3)
    for e in eachindex(cfg.X1)
        _expand!(ext,cfg.X1[e],shift)
        _expand!(ext,cfg.X2[e],shift)
        isnothing(cfg.S1[e]) && continue
        _expand!(ext,cfg.S1[e],shift)
        _expand!(ext,cfg.S2[e],shift)
    end
    return ext
end

# Plot view: plotted axes (two for planar views, three otherwise), camera angles (3D) and unit vector normal to the plane (planar views)
struct _View
    axes::Vector{Int}
    camera::Tuple{Real,Real}
    out::Union{Nothing,SVector{3,Float64}}
end
_is3D(v::_View) = length(v.axes) == 3

# Selects the view: 3D if requested or if surfaces are drawn, planar if the structure is flat in one of the axes
function _select_view(ext::_Extent, view, tolPlane, force3D)
    camera = something(view,(45,45))
    if isnothing(view) && !force3D
        for a in 1:3
            ext.hi[a]-ext.lo[a] < tolPlane && return _View([i for i in 1:3 if i != a],camera,_unit(a))
        end
    end
    return _View([1,2,3],camera,nothing)
end

# Coordinates of the points along the plotted axes
_coords(points, v::_View) = Tuple([p[a] for p in points] for a in v.axes)

# Axis limits for the plotted axes: user-defined limits (expanded to a cube about their center in 3D views, so that the scales of the axes are equal), a cube centered on the bounding box (3D), or padded limits (planar views, if fixed)
function _limits(ext::_Extent, v::_View, plotLimits; avoidDistortion=false, pad=0.1, fixed=true)
    if !isnothing(plotLimits)
        lims = [collect(Float64,plotLimits[a]) for a in v.axes]
        if _is3D(v) && avoidDistortion
            h = maximum(l[2]-l[1] for l in lims)/2
            lims = [(l[1]+l[2])/2 .+ [-h,h] for l in lims]
        end
        return lims
    end
    size = ext.hi - ext.lo
    if _is3D(v)
        c = (ext.lo + ext.hi)/2
        h = (0.5+pad)*maximum(size)
        h > 0 || (h = 1.0)
        return [[c[a]-h,c[a]+h] for a in 1:3]
    elseif fixed
        m = maximum(size[v.axes])
        m > 0 || (m = 1.0)
        return [[ext.lo[a]-pad*m,ext.hi[a]+pad*m] for a in v.axes]
    end
    return nothing
end

function _set_limits!(plt, lims)
    isnothing(lims) && return plt
    plot!(plt, xlims=lims[1], ylims=lims[2])
    length(lims) == 3 && plot!(plt, zlims=lims[3])
    return plt
end

# Creates an empty plot for the view
function _new_plot(model::Model, v::_View; plotAxes=true, plotGrid=true, legendPos=:topright, DPI=300, plotSize=(600,400), kwargs...)
    plt = plot(grid=plotGrid, axis=plotAxes, legend=legendPos, dpi=DPI, size=plotSize; kwargs...)
    plotAxes || plot!(plt, framestyle=:none)
    # In 3D views the limits are cubic (equal scales), and an equal aspect ratio would shrink the plot area as a function of the limits
    _is3D(v) || plot!(plt, aspect_ratio=:equal)
    _is3D(v) && plot!(plt, camera=v.camera)
    if plotAxes
        labels = [string("\$x_",a,"\$ [",model.units.length,"]") for a in v.axes]
        plot!(plt, xlabel=labels[1], ylabel=labels[2])
        # The label of the vertical axis of planar plots may be clipped by the backend, which underestimates the width of the tick labels
        _is3D(v) ? plot!(plt, zlabel=labels[3]) : plot!(plt, left_margin=3mm)
    end
    return plt
end


## STRUCTURE
# ------------------------------------------------------------------------------

"""
    get_undeformed_airfoil_coords(element::Element)

Returns the coordinates of the undeformed airfoil sections at the nodes of the element, as 3×N matrices (nodes 1 and 2), resolved in basis A

# Arguments
- `element::Element`: element (must have an aerodynamic surface)
"""
function get_undeformed_airfoil_coords(element::Element)

    @assert !isnothing(element.aero.surface) "the element has no aerodynamic surface"

    @unpack x1_n1,x1_n2,R0_n1,R0_n2 = element.geometry
    @unpack r_n1,r_n2 = element.assemblyData
    @unpack airfoil,c,normSparPos = element.aero.surface.geometry
    @unpack Rw = element.aero.geometry

    # Airfoil coordinates, rotated to face forward
    coordinates = airfoil.database.coordinates
    Y = 1 .- coordinates[:,1]
    Z = coordinates[:,2]
    N = length(Y)

    # Nodal chords and normalized spar positions
    c1,c2 = c isa Real ? (c,c) : (c(x1_n1),c(x1_n2))
    s1,s2 = normSparPos isa Real ? (normSparPos,normSparPos) : (normSparPos(x1_n1),normSparPos(x1_n2))

    # Sections offset to the spar, scaled by the chord, rotated to the nodal orientation and translated to the nodal position
    section(R0_n,r_n,c_n,s_n) = (Rw*R0_n)*vcat(zeros(1,N),((Y .- (1-s_n))*c_n)',(Z*c_n)') .+ r_n

    return section(R0_n1,r_n1,c1,s1),section(R0_n2,r_n2,c2,s2)
end
export get_undeformed_airfoil_coords

# Draws the beams of a configuration as one polyline per beam
function _draw_beams!(plt, cfg::_Config, v::_View; color, lw, ls=:solid, label=false)
    for (i,ids) in enumerate(cfg.geom.beams)
        points = Vector{SVector{3,Float64}}(undef,2*length(ids))
        for (j,e) in enumerate(ids)
            points[2j-1] = cfg.X1[e]
            points[2j] = cfg.X2[e]
        end
        plot!(plt, _coords(points,v)..., color=color, lw=lw, ls=ls, label=i == 1 ? label : false)
    end
    return plt
end

# Draws a translucent mesh of polygonal faces (tuples of 1-based vertex indices)
function _mesh!(plt, x, y, z, faces; color, alpha=1)
    plot!(plt, x, y, z, seriestype=:mesh3d, connections=faces, fillcolor=color, linecolor=color, lw=0, alpha=alpha, label=false)
    return plt
end

# Draws the aerodynamic surfaces of a configuration as a single mesh, with a strip of quadrilateral faces between the sections at the nodes of each element
function _draw_surfaces!(plt, cfg::_Config, v::_View; color, surfα)
    _is3D(v) || return plt
    x,y,z = Float64[],Float64[],Float64[]
    faces = NTuple{4,Int}[]
    for e in eachindex(cfg.S1)
        isnothing(cfg.S1[e]) && continue
        N = size(cfg.S1[e],2)
        offset = length(x)
        for S in (cfg.S1[e],cfg.S2[e])
            append!(x,@view S[1,:])
            append!(y,@view S[2,:])
            append!(z,@view S[3,:])
        end
        for i in 1:N-1
            push!(faces,(offset+i,offset+i+1,offset+N+i+1,offset+N+i))
        end
    end
    isempty(faces) || _mesh!(plt, x, y, z, faces, color=color, alpha=surfα)
    return plt
end

# Draws a configuration: beams and, optionally, aerodynamic surfaces
function _draw_configuration!(plt, cfg::_Config, v::_View; color, lw, ls=:solid, label=false, surfaces=false, surfα=0.5)
    _draw_beams!(plt, cfg, v, color=color, lw=lw, ls=ls, label=label)
    surfaces && _draw_surfaces!(plt, cfg, v, color=color, surfα=surfα)
    return plt
end

# Saved states, aerodynamic variables, etc. of a problem and their time indices
_history(problem::DynamicProblem) = problem.solutionHistory
_history(problem) = problem.loadstepHistory
_R_A(problem::DynamicProblem, ti) = SMatrix{3,3,Float64,9}(problem.model.kinematicsBasisA.R_AVec[ti])
_R_A(problem, ti) = SMatrix{3,3,Float64,9}(problem.model.kinematicsBasisA.state.R_A)

# Time of the saved step k of a problem (steady problems are drawn at the initial time)
_step_time(problem::DynamicProblem, k) = problem.solutionHistory.t[k]
_step_time(problem, k) = 0.0


## LOADS AND BOUNDARY CONDITIONS
# ------------------------------------------------------------------------------

# Distributed loads of an element: function field, flag field, whether it is resolved in basis b, follower, and moment
const _DIST_LOADS = (
    (fun=:f_A_of_ζt,  has=:has_f_A,  inb=false, follower=false, moment=false),
    (fun=:f_b_of_ζt,  has=:has_f_b,  inb=true,  follower=false, moment=false),
    (fun=:ff_A_of_ζt, has=:has_ff_A, inb=false, follower=true,  moment=false),
    (fun=:ff_b_of_ζt, has=:has_ff_b, inb=true,  follower=true,  moment=false),
    (fun=:m_A_of_ζt,  has=:has_m_A,  inb=false, follower=false, moment=true),
    (fun=:m_b_of_ζt,  has=:has_m_b,  inb=true,  follower=false, moment=true),
    (fun=:mf_A_of_ζt, has=:has_mf_A, inb=false, follower=true,  moment=true),
    (fun=:mf_b_of_ζt, has=:has_mf_b, inb=true,  follower=true,  moment=true),
)

# Normalization of the sizes of the loads drawn: the largest load of each category is drawn with the nominal size, times the user's scaler λ
struct _LoadScales
    λ::Float64
    Fbc::Float64
    Mbc::Float64
    Fdist::Float64
    Mdist::Float64
    cF::Float64
    cM::Float64
end

# Concentrated loads of a BC at time t (evaluated from the prescribed functions, as the BC histories are overwritten by every problem solved on the model): dead and follower loads, resolved in the bases A and b
function _bc_loads(bc, t)
    funs = bc.loadFuns
    return map((funs.deadLoadsOnAofTime,funs.followerLoadsOnAofTime,funs.deadLoadsOnbofTime,funs.followerLoadsOnbofTime)) do f
        SVector{6,Float64}(ntuple(i -> f[i](t),6))
    end
end

function _load_scales(problem, geom::_Geometry, λ, steps; BCs=true, dist=true)
    model = geom.model
    times = [_step_time(problem,k) for k in steps]
    Fbc = Mbc = Fdist = Mdist = cF = cM = 0.0
    if BCs
        for bc in model.BCs, t in times, loads in _bc_loads(bc,t)
            Fbc = max(Fbc,maximum(abs,@view loads[1:3]))
            Mbc = max(Mbc,maximum(abs,@view loads[4:6]))
        end
    end
    if dist
        for element in model.assemblyData.elements, kind in _DIST_LOADS
            getfield(element.distLoads,kind.has) || continue
            f = getfield(element.distLoads,kind.fun)
            m = maximum(norm(f(ζ,t)) for ζ in (0.0,0.5,1.0) for t in times)
            kind.moment ? (Mdist = max(Mdist,m)) : (Fdist = max(Fdist,m))
        end
        for k in steps, av in _history(problem).aeroVariables[k]
            isnothing(av) && continue
            @unpack ct,cn,cm = av.aeroCoefficients
            cF = _nanmax(cF,abs(ct),abs(cn))
            cM = _nanmax(cM,abs(cm))
        end
    end
    return _LoadScales(λ,Fbc,Mbc,Fdist,Mdist,cF,cM)
end

# Adds a stroke (polyline) of the given color to the strokes, which are drawn at once by color
function _stroke!(strokes, color, points)
    stroke = get!(() -> SVector{3,Float64}[], strokes, color)
    append!(stroke,points)
    push!(stroke,_NAN3)
    return nothing
end

# Unit vector perpendicular to the given unit vector
function _perpendicular(t̂)
    ref = abs(t̂[1]) < 0.9 ? _unit(1) : _unit(2)
    return normalize(cross(t̂,ref))
end

# Adds an arrow of vector d that ends (tip) at the point of application, with a head of wings perpendicular to d (or in the plane of the plot, if given its normal `out`)
function _arrow!(strokes, color, tip, d, out)
    ℓ = norm(d)
    ℓ > 0 || return nothing
    t̂ = d/ℓ
    planar = !isnothing(out) && abs(dot(t̂,out)) < 0.99
    n̂ = planar ? normalize(cross(t̂,out)) : _perpendicular(t̂)
    h = ℓ/4
    base = tip - h*t̂
    _stroke!(strokes,color,(tip-d,tip))
    _stroke!(strokes,color,(base+0.4h*n̂,tip,base-0.4h*n̂))
    if isnothing(out)
        m̂ = cross(t̂,n̂)
        _stroke!(strokes,color,(base+0.4h*m̂,tip,base-0.4h*m̂))
    end
    return nothing
end

# Adds a circular arrow (arc of 324° with a V head) of radius r about the axis of the basis R, centered at C. The sense of the arrow follows the sign of the moment M
function _moment!(strokes, color, C, R, axis, M, r)
    r > 0 || return nothing
    s = sign(M)
    e1 = R[:,mod1(axis+1,3)]
    e2 = R[:,mod1(axis+2,3)]
    ϑ = range(0,9π/5,length=30)
    arc = [C + r*(cos(a)*e1 + s*sin(a)*e2) for a in ϑ]
    tangent = -sin(ϑ[end])*e1 + s*cos(ϑ[end])*e2
    radial = cos(ϑ[end])*e1 + s*sin(ϑ[end])*e2
    ah = 0.3r
    _stroke!(strokes,color,arc)
    _stroke!(strokes,color,(arc[end]-ah*tangent+ah/2*radial,arc[end],arc[end]-ah*tangent-ah/2*radial))
    return nothing
end

# Draws a symbol of support (triangle or circle in planar views, cone in 3D) for a generalized displacement BC at the point P about the axis `axis`
function _support!(plt, P, axis, rotation, T, v, size)
    S = 2π*size/200
    _support_symbol!(plt,P,axis,T,v,S,:red)
    rotation && _support_symbol!(plt,P-S/2*T.R[:,axis],axis,T,v,S,:orange)
    return plt
end

function _support_symbol!(plt, P, axis, T, v, S, color)
    d = T.R[:,axis]
    if _is3D(v)
        # Cone with apex at P and base (radius S/2) at distance S along -d
        p1 = T.R[:,mod1(axis+1,3)]
        p2 = T.R[:,mod1(axis+2,3)]
        n = 20
        points = vcat([P], [P - S*d + S/2*(cos(θ)*p1 + sin(θ)*p2) for θ in range(0,2π,length=n+1)[1:n]], [P - S*d])
        faces = NTuple{3,Int}[]
        for i in 1:n
            j = mod1(i+1,n)
            push!(faces,(1,1+i,1+j))
            push!(faces,(n+2,1+j,1+i))
        end
        _mesh!(plt, (getindex.(points,a) for a in 1:3)..., faces, color=color)
    else
        out = v.out
        if abs(dot(d,out)) > 0.5
            points = [P + S/2*(cos(θ)*_unit(v.axes[1]) + sin(θ)*_unit(v.axes[2])) for θ in range(0,2π,length=21)]
        else
            inPlane = normalize(d - dot(d,out)*out)
            n̂ = cross(out,inPlane)
            points = [P,P-S*inPlane+S/2*n̂,P-S*inPlane-S/2*n̂]
        end
        x,y = _coords(points,v)
        plot!(plt, x, y, seriestype=:shape, fillcolor=color, linecolor=color, label=false)
    end
    return plt
end

# Adds the concentrated loads of the BCs at time t to the strokes, and draws the supports of the generalized displacement BCs
function _bc_loads!(plt, strokes, cfg::_Config, v::_View, scales::_LoadScales, t)
    Lref = cfg.geom.Lref
    λ = scales.λ
    TR = cfg.T.R
    for bc in cfg.geom.model.BCs
        @unpack beam,node = bc.definition
        @unpack R0_n,isLoad = bc.deducedData
        e = beam.elements[node == 1 ? 1 : node-1].assemblyData.globalID
        P,Rn = node == 1 ? (cfg.X1[e],cfg.poses[e].R1) : (cfg.X2[e],cfg.poses[e].R2)
        deadA,followerA,deadb,followerb = _bc_loads(bc,t)
        loads = ((deadA,_I3),(followerA,Rn),(deadb,R0_n),(followerb,Rn*R0_n))
        for DOF in 1:6
            axis = mod1(DOF,3)
            if !isLoad[DOF]
                _support!(plt,P,axis,DOF > 3,cfg.T,v,Lref*λ)
                continue
            end
            for (values,R) in loads
                value = values[DOF]
                value == 0 && continue
                if DOF <= 3
                    _arrow!(strokes,:green,P,Lref*λ/10*value/scales.Fbc*(TR*R)[:,axis],v.out)
                else
                    _moment!(strokes,:green,P,TR*R,axis,value,Lref*λ/20*abs(value)/scales.Mbc)
                end
            end
        end
    end
    return nothing
end

# Adds the distributed (including gravitational and aerodynamic) loads at the saved step k (time index ti, time t) to the strokes
function _distributed_loads!(strokes, cfg::_Config, v::_View, problem, scales::_LoadScales, k, ti, t)
    model = cfg.geom.model
    Lref,λ,TR = cfg.geom.Lref,scales.λ,cfg.T.R
    g = SVector{3,Float64}(model.gravityVector)
    gravity = TR*(_R_A(problem,ti)'*g)
    aeroVariables = _history(problem).aeroVariables[k]
    for (e,element) in enumerate(model.assemblyData.elements)
        pose = cfg.poses[e]
        P = (cfg.X1[e],(cfg.X1[e]+cfg.X2[e])/2,cfg.X2[e])
        R = (pose.R1,pose.Rm,pose.R2)
        R0 = (element.geometry.R0_n1,element.geometry.R0,element.geometry.R0_n2)
        if !iszero(g) && element.sectional.μ > 0
            for j in 1:3
                _arrow!(strokes,:gold,P[j],Lref*λ/10/norm(g)*gravity,v.out)
            end
        end
        for kind in _DIST_LOADS
            getfield(element.distLoads,kind.has) || continue
            f = getfield(element.distLoads,kind.fun)
            for (j,ζ) in enumerate((0.0,0.5,1.0))
                value = SVector{3,Float64}(f(ζ,t))
                M = TR*(kind.inb ? (kind.follower ? R[j]*R0[j] : R0[j]) : (kind.follower ? R[j] : _I3))
                if kind.moment
                    for axis in 1:3
                        _moment!(strokes,:purple,P[j],M,axis,value[axis],Lref*λ/40*abs(value[axis])/scales.Mdist)
                    end
                else
                    _arrow!(strokes,:purple,P[j],Lref*λ/10/scales.Fdist*(M*value),v.out)
                end
            end
        end
        aero = aeroVariables[e]
        if !isnothing(aero) && !isnothing(element.aero.surface)
            @unpack ct,cn,cm = aero.aeroCoefficients
            M = TR*R[2]*element.aero.geometry.RwR0
            _arrow!(strokes,:green,P[2],Lref*λ/4*ct/scales.cF*M[:,2],v.out)
            _arrow!(strokes,:green,P[2],Lref*λ/4*cn/scales.cF*M[:,3],v.out)
            _moment!(strokes,:green,P[2],M,1,cm,Lref*λ/20*abs(cm)/scales.cM)
        end
    end
    return nothing
end

# Loads to be drawn on a configuration: concentrated (BCs) and/or distributed, at the saved step k (time index ti, time t) of the problem
struct _LoadSpec
    problem::Problem
    scales::_LoadScales
    k::Int
    ti::Int
    t::Float64
    BCs::Bool
    dist::Bool
end

# Draws the loads and supports of the specification
function _draw_loads!(plt, cfg::_Config, v::_View, spec::_LoadSpec)
    strokes = Dict{Symbol,Vector{SVector{3,Float64}}}()
    spec.BCs && _bc_loads!(plt,strokes,cfg,v,spec.scales,spec.t)
    spec.dist && _distributed_loads!(strokes,cfg,v,spec.problem,spec.scales,spec.k,spec.ti,spec.t)
    for (color,points) in strokes
        plot!(plt, _coords(points,v)..., color=color, lw=1, label=false)
    end
    return plt
end


## FRAMES
# ------------------------------------------------------------------------------

# A layer of a frame: a configuration, drawn with the given attributes (no shape is drawn if the color is nothing), and the loads to be drawn on it
struct _Layer
    cfg::_Config
    color
    ls::Symbol
    label
    surfaces::Bool
    loads::Union{Nothing,_LoadSpec}
end
_Layer(cfg::_Config; color=nothing, ls=:solid, label=false, surfaces=false, loads=nothing) = _Layer(cfg,color,ls,label,surfaces,loads)

# Expands the bounding box with the shapes drawn in the layers
function _expand!(ext::_Extent, layers::Vector{_Layer}, shift=_O3)
    for layer in layers
        isnothing(layer.color) || _expand!(ext,layer.cfg,shift)
    end
    return ext
end

# Whether aerodynamic surfaces are drawn in the layers
_has_surfaces(layers) = any(layer -> layer.surfaces && any(!isnothing,layer.cfg.S1), layers)

# Adds a legend entry (without data) to the plot
function _legend_entry!(plt, v::_View, color, lw, label)
    plot!(plt, (fill(NaN,1) for _ in v.axes)..., color=color, lw=lw, label=label)
    return plt
end

# Position of an annotation in a planar view, relative to the limits: from the lower limit along x1 and from the upper limit along x2 and x3
function _note_position(v::_View, lims, pos)
    a,b = v.axes
    x = lims[1][a == 1 ? 1 : 2] + pos[a]*(lims[1][2]-lims[1][1])
    y = lims[2][2] + pos[b]*(lims[2][2]-lims[2][1])
    return x,y
end

# Draws the layers in a new plot, with the given axis limits. The keyword arguments for the plot are those of _new_plot
function _draw_frame(model::Model, layers::Vector{_Layer}, v::_View, lims; lw, surfα, plotKwargs...)
    plt = _new_plot(model, v; plotKwargs...)
    for layer in layers
        isnothing(layer.color) || _draw_configuration!(plt, layer.cfg, v, color=layer.color, lw=lw, ls=layer.ls, label=layer.label, surfaces=layer.surfaces, surfα=surfα)
        isnothing(layer.loads) || _draw_loads!(plt, layer.cfg, v, layer.loads)
    end
    _set_limits!(plt, lims)
    return plt
end

# Plots a single frame. The view and the limits are determined by the layers, and the plot is decorated by decorate(plt,v), if given
function _plot_frame(model::Model, layers::Vector{_Layer}; view, tolPlane, plotLimits, avoidDistortion=false, decorate=nothing, drawKwargs...)
    ext = _expand!(_Extent(),layers)
    v = _select_view(ext, view, tolPlane, _has_surfaces(layers))
    lims = _limits(ext, v, plotLimits; avoidDistortion=avoidDistortion, fixed=false)
    plt = _draw_frame(model, layers, v, lims; drawKwargs...)
    isnothing(decorate) || decorate(plt,v)
    return plt
end

# Animates the frames, each with the layers given by layersOf(frame) and decorated by decorate(plt,v,lims,frame), where lims are the axis limits of the plot. The view and limits are fixed over the animation (or follow the movement of basis A), which requires a first pass over the frames: the configurations are recomputed to be drawn, rather than stored
function _animate(model::Model, frames, layersOf; view, tolPlane, plotLimits, avoidDistortion, follow, fps, save, savePath, displayProgress, decorate, drawKwargs...)

    @assert !isempty(frames) "there are no frames to animate"

    # Following basis A, the limits are shifted by the translation of the observer's transformation; following the assembly, by the position of its anchor (both of the last layer, if there are many)
    shiftOf(layers) = follow == :basisA ? layers[end].cfg.T.u : follow == :assembly ? _anchor_position(layers[end].cfg) : _O3

    # View and limits, from the extent over all frames
    ext = _Extent()
    force3D = false
    for f in frames
        layers = layersOf(f)
        _expand!(ext,layers,shiftOf(layers))
        force3D |= _has_surfaces(layers)
    end
    v = _select_view(ext, view, tolPlane, force3D)
    lims = _limits(ext, v, plotLimits; avoidDistortion=avoidDistortion)

    # Frames
    anim = Animation()
    for (i,f) in enumerate(frames)
        layers = layersOf(f)
        shift = shiftOf(layers)
        limsNow = [lims[j] .+ shift[a] for (j,a) in enumerate(v.axes)]
        plt = _draw_frame(model, layers, v, limsNow; drawKwargs...)
        decorate(plt,v,limsNow,f)
        frame(anim,plt)
        displayProgress && println("Animation progress: ",round(i/length(frames)*100,digits=1)," %")
    end

    # Gif
    if save
        path = _full_path(savePath)
        mkpath(dirname(path))
        return gif(anim,path,fps=fps)
    else
        return gif(anim,fps=fps)
    end
end


## STEADY PROBLEMS
# ------------------------------------------------------------------------------

"""
    plot_undeformed_assembly(model::Model; kwargs...)

Plots the undeformed assembly of beams in the model

# Arguments
- `model::Model`

# Keyword arguments
- `backend::Symbol`: backend for plot
- `view::Tuple{Real,Real}`: view angles
- `equalAspectRatio::Bool`: flag to set equal aspect ratio plot
- `plotNodes::Bool`: flag to plot nodes
- `nodesColor=:black`: color of nodes
- `linesColor=:black`: color of lines (beams)
- `DPI::Int`: pixel density
- `plotSize::Tuple{Real,Real}`: size of the plot (width, height) in pixels
- `plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}`: plot axis limits
- `avoidDistortion::Bool`: flag to avoid distortion in 3D plots
"""
function plot_undeformed_assembly(model::Model; backend::Symbol=:gr, view::Tuple{Real,Real}=(45,45), equalAspectRatio::Bool=true, plotNodes::Bool=true, nodesColor=:black, linesColor=:black, DPI::Int=300, plotSize::Tuple{Real,Real}=(600,400), plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}=nothing, avoidDistortion::Bool=false)

    _check_options(backend,DPI,plotLimits,0.5,plotSize)
    _set_backend(backend)

    geom = _geometry(model, surfaces=false)
    cfg = _undeformed_config(geom)
    ext = _expand!(_Extent(),cfg)
    v = _View([1,2,3],view,nothing)

    plt = _new_plot(model, v, legendPos=false, DPI=DPI, plotSize=plotSize, title="Undeformed assembly")
    _draw_beams!(plt, cfg, v, color=linesColor, lw=2)
    plotNodes && scatter!(plt, _coords(model.assemblyData.r_n,v)..., c=nodesColor, ms=3, label=false)
    _set_limits!(plt, equalAspectRatio || !isnothing(plotLimits) ? _limits(ext,v,plotLimits; avoidDistortion=avoidDistortion) : nothing)

    return plt
end
export plot_undeformed_assembly

# A deformed configuration to be plotted: the saved step of a problem (steady, trim or eigen), with its plot attributes
struct _SteadyEntry
    problem::Problem
    step::Int
    centralize::Union{Nothing,Int}
    color
    ls::Symbol
    label
    final::Bool
end

# Writes the deformation scale on the plot (as a legend entry in 3D views, at the given position of the axes otherwise)
function _show_scale!(plt, v::_View, scale, pos)
    label = string("Deformation scale: ",scale,"×")
    _is3D(v) ? _legend_entry!(plt, v, :white, 1, label) : annotate!(plt, pos[1], pos[2], text(label,8))
    return plt
end

# Checks that the problem is a steady-like problem which has been solved
function _check_steady_problem(problem)
    @assert problem isa Union{SteadyProblem,TrimProblem,EigenProblem} "the problem must be a steady, trim or eigen problem"
    @assert !isempty(problem.loadstepHistory.states) "solve the problem before plotting"
    return nothing
end

# Plots the deformed configurations of the entries, over the undeformed assembly of the first one. Loads are plotted for the final entries only
function _plot_steady_entries(entries::Vector{_SteadyEntry}; plotBCs, plotUndeformed, plotDistLoads, view, scale, lw, colorUndef, undefLabel, plotAxes, plotGrid, legendPos, tolPlane, plotAeroSurf, surfα, plotLimits, avoidDistortion, showScale, scalePos, loadsSizeScaler, DPI, plotSize, save, savePath)

    geoms = IdDict{Any,_Geometry}()
    geometry(problem) = get!(() -> _geometry(problem.model, surfaces=plotAeroSurf), geoms, problem.model)

    # Layers: undeformed assembly and deformed configurations
    layers = _Layer[]
    plotUndeformed && push!(layers, _Layer(_undeformed_config(geometry(entries[1].problem)), color=colorUndef, label=undefLabel, surfaces=plotAeroSurf))
    for en in entries
        geom = geometry(en.problem)
        cfg = _configuration(geom, _poses(_history(en.problem).states[en.step],scale,centralize=en.centralize))
        loads = en.final && (plotBCs || plotDistLoads) ?
        _LoadSpec(en.problem, _load_scales(en.problem,geom,loadsSizeScaler,(en.step,),BCs=plotBCs,dist=plotDistLoads), en.step, 1, 0.0, plotBCs, plotDistLoads) :
        nothing
        push!(layers, _Layer(cfg, color=en.color, ls=en.ls, label=en.label, surfaces=en.final && plotAeroSurf, loads=loads))
    end

    plt = _plot_frame(entries[1].problem.model, layers, view=view, tolPlane=tolPlane, plotLimits=plotLimits, avoidDistortion=avoidDistortion, decorate=(plt,v) -> showScale && _show_scale!(plt,v,scale,scalePos), lw=lw, surfα=surfα, plotAxes=plotAxes, plotGrid=plotGrid, legendPos=legendPos, DPI=DPI, plotSize=plotSize)
    save && _save_figure(plt,savePath)

    return plt
end

"""
    plot_steady_deformation(problem::Problem; kwargs...)

Plots the deformed state of the model in the given problem

# Arguments
- `problem::Problem`: problem

# Keyword arguments
- `backend::Symbol`: backend for plot
- `plotBCs::Bool`: flag to plot BCs
- `plotUndeformed::Bool`: flag to plot undeformed assembly
- `plotDistLoads::Bool`: flag to plot distributed loads (includes gravitational and aerodynamic loads)
- `view::Union{Nothing,Tuple{Real,Real}}`: view angles (by default, the plot is planar if the assembly lies on a plane of the basis A, and 3D otherwise)
- `scale::Real`: displacements and rotations scale
- `lw::Real`: linewidth
- `colorUndef`: color for undeformed assembly
- `colorDef`: color for deformed assembly
- `plotLegend::Bool`: flag to plot legend
- `plotAxes::Bool`: flag to plot axes
- `plotGrid::Bool`: flag to plot grid
- `legendPos`: legend position
- `tolPlane::Real`: tolerance of the extent of the assembly along an axis to plot as plane
- `plotAeroSurf::Bool`: flag to plot aerodynamic surfaces (this forces a 3D plot)
- `surfα::Float64`: transparency factor of aerodynamic surfaces
- `element2centralize::Union{Nothing,Int}`: element for displacement offset
- `σ2plot::Vector{<:Real}`: (tentative) load factors to plot
- `partialLoadStepsLineStyle::Symbol`: linestyle for partial load steps
- `plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}`: plot axis limits
- `showScale::Bool`: flag to show scale on plot
- `scalePos::Vector{<:Real}`: position of scale on plot
- `loadsSizeScaler::Real`: scaler for plot size of distributed and concentrated loads
- `DPI::Int`: pixel density
- `plotSize::Tuple{Real,Real}`: size of the plot (width, height) in pixels
- `avoidDistortion::Bool`: flag to avoid distortion in 3D plots
- `save::Bool`: flag to save the figure
- `savePath::String`: relative path (from the package's directory) on which to save the figure
"""
function plot_steady_deformation(problem::Problem; backend::Symbol=:gr, plotBCs::Bool=true, plotUndeformed::Bool=true, plotDistLoads::Bool=true, view::Union{Nothing,Tuple{Real,Real}}=nothing, scale::Real=1, lw::Real=1, colorUndef=:black, colorDef=:blue, plotLegend::Bool=true, plotAxes::Bool=true, plotGrid::Bool=true, legendPos=:topright, tolPlane::Real=1e-8, plotAeroSurf::Bool=true, surfα::Float64=0.5, element2centralize::Union{Nothing,Int}=nothing, σ2plot::Vector{<:Real}=[1], partialLoadStepsLineStyle::Symbol=:dot, plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}=nothing, showScale::Bool=false, scalePos::Vector{<:Real}=[0,0], loadsSizeScaler::Real=1, DPI::Int=300, plotSize::Tuple{Real,Real}=(600,400), avoidDistortion::Bool=false, save::Bool=false, savePath::String="/test/outputs/figures/fig.pdf")

    _check_steady_problem(problem)
    _check_options(backend,DPI,plotLimits,surfα,plotSize, scale=scale, lw=lw, tolPlane=tolPlane, loadsSizeScaler=loadsSizeScaler)
    @assert isnothing(element2centralize) || element2centralize <= problem.model.assemblyData.nElementsTotal "element2centralize must be at most the number of elements"
    @assert !isempty(σ2plot) "set at least one load factor in σ2plot"
    _set_backend(backend)

    # Closest saved load steps, the last of which is the final one
    steps = [argmin(abs.(problem.loadstepHistory.σ .- σ)) for σ in σ2plot]
    entries = [_SteadyEntry(problem, k, element2centralize, colorDef, i == length(steps) ? :solid : partialLoadStepsLineStyle, i == length(steps) && plotLegend ? "Deformed" : false, i == length(steps)) for (i,k) in enumerate(steps)]

    return _plot_steady_entries(entries, plotBCs=plotBCs, plotUndeformed=plotUndeformed, plotDistLoads=plotDistLoads, view=view, scale=scale, lw=lw, colorUndef=colorUndef, undefLabel=plotLegend ? "Undeformed" : false, plotAxes=plotAxes, plotGrid=plotGrid, legendPos=legendPos, tolPlane=tolPlane, plotAeroSurf=plotAeroSurf, surfα=surfα, plotLimits=plotLimits, avoidDistortion=avoidDistortion, showScale=showScale, scalePos=scalePos, loadsSizeScaler=loadsSizeScaler, DPI=DPI, plotSize=plotSize, save=save, savePath=savePath)
end
export plot_steady_deformation

"""
    plot_steady_deformations(problems::Array{<:Problem}; kwargs...)

Plots the deformed states of the models in the given problems, over the undeformed assembly of the first problem

# Arguments
- `problems::Array{<:Problem}`: problems

# Keyword arguments
- `backend::Symbol`: backend for plot
- `plotUndeformed::Bool`: flag to plot undeformed assembly
- `plotBCs::Bool`: flag to plot BCs
- `plotDistLoads::Bool`: flag to plot distributed loads (includes gravitational and aerodynamic loads)
- `view::Union{Nothing,Tuple{Real,Real}}`: view angles (by default, the plot is planar if the assemblies lie on a plane of the basis A, and 3D otherwise)
- `scale::Real`: displacements and rotations scale
- `lw::Real`: linewidth
- `colorUndef`: color for undeformed assembly
- `colorScheme`: color scheme for deformed assemblies
- `plotAxes::Bool`: flag to plot axes
- `plotGrid::Bool`: flag to plot grid
- `legendPos`: legend position
- `legendEntries::Union{Nothing,Array{String}}`: legend entries (one for each problem)
- `tolPlane::Real`: tolerance of the extent of the assemblies along an axis to plot as plane
- `plotAeroSurf::Bool`: flag to plot aerodynamic surfaces (this forces a 3D plot)
- `surfα::Float64`: transparency factor of aerodynamic surfaces
- `elements2centralize::Union{Nothing,AbstractArray{Int}}`: elements for displacement offset (one for each problem)
- `plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}`: plot axis limits
- `loadsSizeScaler::Real`: scaler for plot size of distributed and concentrated loads
- `DPI::Int`: pixel density
- `plotSize::Tuple{Real,Real}`: size of the plot (width, height) in pixels
- `avoidDistortion::Bool`: flag to avoid distortion in 3D plots
- `save::Bool`: flag to save the figure
- `savePath::String`: relative path (from the package's directory) on which to save the figure
"""
function plot_steady_deformations(problems::Array{<:Problem}; backend::Symbol=:gr, plotUndeformed::Bool=true, plotBCs::Bool=true, plotDistLoads::Bool=true, view::Union{Nothing,Tuple{Real,Real}}=nothing, scale::Real=1, lw::Real=2, colorUndef=:black, colorScheme=:rainbow, plotAxes::Bool=true, plotGrid::Bool=true, legendPos=:topright, legendEntries::Union{Nothing,Array{String}}=nothing, tolPlane::Real=1e-8, plotAeroSurf::Bool=true, surfα::Float64=0.5, elements2centralize::Union{Nothing,AbstractArray{Int}}=nothing, plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}=nothing, loadsSizeScaler::Real=1, DPI::Int=300, plotSize::Tuple{Real,Real}=(600,400), avoidDistortion::Bool=false, save::Bool=false, savePath::String="/test/outputs/figures/fig.pdf")

    @assert !isempty(problems) "set at least one problem"
    foreach(_check_steady_problem, problems)
    _check_options(backend,DPI,plotLimits,surfα,plotSize, scale=scale, lw=lw, tolPlane=tolPlane, loadsSizeScaler=loadsSizeScaler)
    @assert isnothing(legendEntries) || length(legendEntries) == length(problems) "specify legendEntries with the same size as problems"
    @assert isnothing(elements2centralize) || length(elements2centralize) == length(problems) "specify elements2centralize with the same size as problems"
    @assert isnothing(elements2centralize) || all(elements2centralize[i] <= problem.model.assemblyData.nElementsTotal for (i,problem) in enumerate(problems)) "elements2centralize must be at most the number of elements"
    _set_backend(backend)

    colors = cgrad(colorScheme, length(problems), categorical=true)
    entries = [_SteadyEntry(problem, length(problem.loadstepHistory.states), isnothing(elements2centralize) ? nothing : elements2centralize[i], colors[i], :solid, isnothing(legendEntries) ? false : legendEntries[i], true) for (i,problem) in enumerate(problems)]

    return _plot_steady_entries(entries, plotBCs=plotBCs, plotUndeformed=plotUndeformed, plotDistLoads=plotDistLoads, view=view, scale=scale, lw=lw, colorUndef=colorUndef, undefLabel=isnothing(legendEntries) ? false : "Undeformed", plotAxes=plotAxes, plotGrid=plotGrid, legendPos=legendPos, tolPlane=tolPlane, plotAeroSurf=plotAeroSurf, surfα=surfα, plotLimits=plotLimits, avoidDistortion=avoidDistortion, showScale=false, scalePos=[0,0], loadsSizeScaler=loadsSizeScaler, DPI=DPI, plotSize=plotSize, save=save, savePath=savePath)
end
export plot_steady_deformations


## MODE SHAPES
# ------------------------------------------------------------------------------

# Conversion factors from rad/s to the frequency units of the model
const _FREQUENCY_FACTORS = Dict("rad/s" => 1.0, "Hz" => 1/(2π), "rpm" => 60/(2π))

# A mode to be plotted, with its color and label
struct _ModeEntry
    shape::ModeShape
    color
    label
end

# Selects the modes to be plotted, with their colors and labels. By default, these are the oscillatory modes (those of positive frequency, numbered in ascending order of frequency) of the problem; alternatively, the mode shapes (e.g., from `filter_tracked_modes`) and their numbers may be given
function _select_modes(problem::EigenProblem, nModes, modes2plot, modeLabels, modalColorScheme, numberOfColors, nullDampingTol, showLegend; modeShapes=nothing, modeNumbers=nothing)

    if isnothing(modeShapes)
        @assert isnothing(modeNumbers) "modeNumbers can only be input together with modeShapes"
        shapes = filter(mode -> mode.frequency > eps(), problem.eigensolution.modeShapes)
        numbers = collect(1:length(shapes))
        @assert !isempty(shapes) "solve the problem before plotting (it must have oscillatory modes)"
        @assert isnothing(nModes) || isnothing(modes2plot) "input either nModes or modes2plot"
        modes = isnothing(modes2plot) ? collect(1:something(nModes,length(shapes))) : sort(collect(modes2plot))
        @assert !isempty(modes) && modes[1] >= 1 && modes[end] <= length(shapes) "the problem has $(length(shapes)) oscillatory modes"
    else
        @assert isnothing(nModes) && isnothing(modes2plot) "nModes and modes2plot cannot be input together with modeShapes"
        @assert isnothing(modeNumbers) || length(modeNumbers) == length(modeShapes) "set one mode number for each of the $(length(modeShapes)) mode shapes"
        # Absent mode shapes (e.g., tracks outside the frequency bounds at this point) are skipped
        present = findall(!isnothing, modeShapes)
        @assert !isempty(present) "modeShapes has no mode shapes to plot"
        shapes = ModeShape[modeShapes[i] for i in present]
        numbers = isnothing(modeNumbers) ? present : collect(modeNumbers[present])
        modes = collect(1:length(shapes))
    end
    @assert isnothing(modeLabels) || length(modeLabels) == length(modes) "set one label for each of the $(length(modes)) modes plotted"
    @assert ismissing(numberOfColors) || numberOfColors >= length(modes) "set numberOfColors as at least the number of modes plotted"

    # Colors by position among the modes plotted or, if the number of colors is given, by mode number
    nColors = ismissing(numberOfColors) ? length(modes) : numberOfColors
    colors = cgrad(modalColorScheme, nColors, categorical=true)

    # Labels with the frequency (and damping, if not null)
    unit = problem.model.units.frequency
    γ = _FREQUENCY_FACTORS[unit]
    function label(m,mode)
        name = isnothing(modeLabels) ? "Mode $(numbers[mode])" : modeLabels[m]
        ω = round(γ*shapes[mode].frequency,digits=2)
        abs(γ*shapes[mode].damping) < nullDampingTol && return string(name,": ω = ",ω," ",unit)
        return string(name,": σ = ",round(γ*shapes[mode].damping,sigdigits=3)," ± ",ω,"i ",unit)
    end

    return [_ModeEntry(shapes[mode], colors[ismissing(numberOfColors) ? m : mod1(numbers[mode],nColors)], showLegend ? label(m,mode) : false) for (m,mode) in enumerate(modes)]
end

# Returns the function that gives the layers of the steady configuration (with the BCs, if applicable) and the modes, given the amplitudes of the modes
function _mode_layers(problem::EigenProblem, entries::Vector{_ModeEntry}; plotSteady, colorSteady, steadyLabel, plotBCs, plotAeroSurf, element2centralize)
    geom = _geometry(problem.model, surfaces=plotAeroSurf)
    k = length(problem.loadstepHistory.states)
    states = problem.loadstepHistory.states[k]
    loads = plotBCs ? _LoadSpec(problem, _load_scales(problem,geom,1.0,(k,),dist=false), k, 1, 0.0, true, false) : nothing
    steady = _Layer(_configuration(geom,_poses(states,1,centralize=element2centralize)), color=plotSteady ? colorSteady : nothing, label=plotSteady ? steadyLabel : false, surfaces=plotAeroSurf, loads=loads)
    return amplitudes -> vcat(steady, [_Layer(_configuration(geom,_poses(states,1,modal=entry.shape.state,factor=a,centralize=element2centralize)), color=entry.color, label=entry.label, surfaces=plotAeroSurf) for (entry,a) in zip(entries,amplitudes)])
end

"""
    plot_mode_shapes(problem::EigenProblem; kwargs...)

Plots the mode shapes of the model in the given eigenproblem, over its steady configuration

# Arguments
- `problem::EigenProblem`: problem

# Keyword arguments
- `backend::Symbol`: backend for plot
- `plotBCs::Bool`: flag to plot BCs
- `view::Union{Nothing,Tuple{Real,Real}}`: view angles (by default, the plot is planar if the assembly and its modes lie on a plane of the basis A, and 3D otherwise)
- `nModes::Union{Nothing,Int}`: number of oscillatory modes to plot (first ones)
- `modes2plot::Union{Nothing,AbstractArray{Int}}`: array of oscillatory modes to plot
- `modeShapes::Union{Nothing,AbstractVector{<:Union{ModeShape,Nothing}}}`: mode shapes to plot, instead of the problem's oscillatory modes (e.g., a column of the mode shapes output by [`filter_tracked_modes`](@ref filter_tracked_modes), whose `nothing` entries are skipped). Cannot be input with `nModes` or `modes2plot`
- `modeNumbers::Union{Nothing,AbstractVector{Int}}`: numbers of the mode shapes in `modeShapes`, used in the labels and colors (e.g., the kept tracks output by [`filter_tracked_modes`](@ref filter_tracked_modes)). By default, the positions in `modeShapes`
- `scale::Real`: displacements and rotations scale
- `element2centralize::Union{Nothing,Int}`: element for displacement offset
- `nullDampingTol::Real`: tolerance for consideration of null damping
- `plotSteady::Bool`: flag to plot steady solution
- `lw::Real`: linewidth
- `colorSteady`: color for steadily deformed assembly
- `modalColorScheme`: color scheme for mode shapes
- `numberOfColors::Union{Missing,Int}`: number of colors in legend color gradient (if given, the modes are colored by their numbers)
- `plotAxes::Bool`: flag to plot axes
- `plotGrid::Bool`: flag to plot grid
- `plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}`: plot axis limits
- `legendPos`: legend position
- `modeLabels::Union{Nothing,Vector{String}}`: labels for each mode
- `tolPlane::Real`: tolerance of the extent of the assembly along an axis to plot as plane
- `plotAeroSurf::Bool`: flag to plot aerodynamic surfaces (this forces a 3D plot)
- `surfα::Float64`: transparency factor of aerodynamic surfaces
- `DPI::Int`: pixel density
- `plotSize::Tuple{Real,Real}`: size of the plot (width, height) in pixels
- `avoidDistortion::Bool`: flag to avoid distortion in 3D plots
- `save::Bool`: flag to save the figure
- `savePath::String`: relative path (from the package's directory) on which to save the figure
"""
function plot_mode_shapes(problem::EigenProblem; backend::Symbol=:gr, plotBCs::Bool=true, view::Union{Nothing,Tuple{Real,Real}}=nothing, nModes::Union{Nothing,Int}=nothing, modes2plot::Union{Nothing,AbstractArray{Int}}=nothing, modeShapes::Union{Nothing,AbstractVector{<:Union{ModeShape,Nothing}}}=nothing, modeNumbers::Union{Nothing,AbstractVector{Int}}=nothing, scale::Real=1, element2centralize::Union{Nothing,Int}=nothing, nullDampingTol::Real=1e-5, plotSteady::Bool=true, lw::Real=2, colorSteady=:black, modalColorScheme=:rainbow, numberOfColors::Union{Missing,Int}=missing, plotAxes::Bool=true, plotGrid::Bool=true, plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}=nothing, legendPos=:topright, modeLabels::Union{Nothing,Vector{String}}=nothing, tolPlane::Real=1e-8, plotAeroSurf::Bool=true, surfα::Float64=0.5, DPI::Int=300, plotSize::Tuple{Real,Real}=(600,400), avoidDistortion::Bool=false, save::Bool=false, savePath::String="/test/outputs/figures/fig.pdf")

    _check_steady_problem(problem)
    _check_options(backend,DPI,plotLimits,surfα,plotSize, lw=lw, tolPlane=tolPlane, nullDampingTol=nullDampingTol)
    @assert isnothing(element2centralize) || element2centralize <= problem.model.assemblyData.nElementsTotal "element2centralize must be at most the number of elements"
    entries = _select_modes(problem, nModes, modes2plot, modeLabels, modalColorScheme, numberOfColors, nullDampingTol, true, modeShapes=modeShapes, modeNumbers=modeNumbers)
    _set_backend(backend)

    layersOf = _mode_layers(problem, entries, plotSteady=plotSteady, colorSteady=colorSteady, steadyLabel="Steady", plotBCs=plotBCs, plotAeroSurf=plotAeroSurf, element2centralize=element2centralize)
    plt = _plot_frame(problem.model, layersOf(fill(scale,length(entries))), view=view, tolPlane=tolPlane, plotLimits=plotLimits, avoidDistortion=avoidDistortion, lw=lw, surfα=surfα, plotAxes=plotAxes, plotGrid=plotGrid, legendPos=legendPos, DPI=DPI, plotSize=plotSize)
    save && _save_figure(plt,savePath)

    return plt
end
export plot_mode_shapes

"""
    plot_mode_shapes_animation(problem::EigenProblem; kwargs...)

Plots the animation of the mode shapes of the model in the given eigenproblem, over its steady configuration

# Arguments
- `problem::EigenProblem`: problem

# Keyword arguments
- `backend::Symbol`: backend for plot
- `plotBCs::Bool`: flag to plot BCs
- `view::Union{Nothing,Tuple{Real,Real}}`: view angles (by default, the plot is planar if the assembly and its modes lie on a plane of the basis A, and 3D otherwise)
- `nModes::Union{Nothing,Int}`: number of oscillatory modes to plot (first ones)
- `modes2plot::Union{Nothing,AbstractArray{Int}}`: array of oscillatory modes to plot
- `modeShapes::Union{Nothing,AbstractVector{<:Union{ModeShape,Nothing}}}`: mode shapes to plot, instead of the problem's oscillatory modes (e.g., a column of the mode shapes output by [`filter_tracked_modes`](@ref filter_tracked_modes), whose `nothing` entries are skipped). Cannot be input with `nModes` or `modes2plot`
- `modeNumbers::Union{Nothing,AbstractVector{Int}}`: numbers of the mode shapes in `modeShapes`, used in the labels and colors (e.g., the kept tracks output by [`filter_tracked_modes`](@ref filter_tracked_modes)). By default, the positions in `modeShapes`
- `matchModeFrequency::Bool`: flag to match frequencies of each mode in animation
- `scale::Real`: displacements and rotations scale
- `nFramesPerCycle::Int`: number of frames over vibration cycle
- `fps::Int`: frame rate for gif
- `element2centralize::Union{Nothing,Int}`: element for displacement offset
- `nullDampingTol::Real`: tolerance for consideration of null damping
- `plotSteady::Bool`: flag to plot steady solution
- `lw::Real`: linewidth
- `colorSteady`: color for steadily deformed assembly
- `modalColorScheme`: color scheme for mode shapes
- `numberOfColors::Union{Missing,Int}`: number of colors in legend color gradient (if given, the modes are colored by their numbers)
- `plotAxes::Bool`: flag to plot axes
- `plotGrid::Bool`: flag to plot grid
- `plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}`: plot axis limits
- `showLegend::Bool`: flag to show legend
- `legendFontSize::Int`: legend font size
- `legendPos`: legend position
- `modeLabels::Union{Nothing,Vector{String}}`: labels for each mode
- `tolPlane::Real`: tolerance of the extent of the assembly along an axis to plot as plane
- `plotAeroSurf::Bool`: flag to plot aerodynamic surfaces (this forces a 3D plot)
- `surfα::Float64`: transparency factor of aerodynamic surfaces
- `DPI::Int`: pixel density
- `plotSize::Tuple{Real,Real}`: size of the plot (width, height) in pixels
- `avoidDistortion::Bool`: flag to avoid distortion in 3D plots
- `save::Bool`: flag to save the figure
- `savePath::String`: relative path (from the package's directory) on which to save the figure (must have the .gif extension)
- `displayProgress::Bool`: flag to display progress of gif creation
"""
function plot_mode_shapes_animation(problem::EigenProblem; backend::Symbol=:gr, plotBCs::Bool=true, view::Union{Nothing,Tuple{Real,Real}}=nothing, nModes::Union{Nothing,Int}=nothing, modes2plot::Union{Nothing,AbstractArray{Int}}=nothing, modeShapes::Union{Nothing,AbstractVector{<:Union{ModeShape,Nothing}}}=nothing, modeNumbers::Union{Nothing,AbstractVector{Int}}=nothing, matchModeFrequency::Bool=true, scale::Real=1, nFramesPerCycle::Int=51, fps::Int=30, element2centralize::Union{Nothing,Int}=nothing, nullDampingTol::Real=1e-5, plotSteady::Bool=true, lw::Real=2, colorSteady=:black, modalColorScheme=:rainbow, numberOfColors::Union{Missing,Int}=missing, plotAxes::Bool=true, plotGrid::Bool=true, plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}=nothing, showLegend::Bool=true, legendFontSize::Int=12, legendPos=:topright, modeLabels::Union{Nothing,Vector{String}}=nothing, tolPlane::Real=1e-8, plotAeroSurf::Bool=true, surfα::Float64=0.5, DPI::Int=300, plotSize::Tuple{Real,Real}=(600,400), avoidDistortion::Bool=false, save::Bool=false, savePath::String="/test/outputs/figures/anim.gif", displayProgress::Bool=false)

    _check_steady_problem(problem)
    _check_options(backend,DPI,plotLimits,surfα,plotSize, lw=lw, tolPlane=tolPlane, nullDampingTol=nullDampingTol, fps=fps)
    @assert nFramesPerCycle >= 3 "nFramesPerCycle must be greater than or equal to 3"
    @assert endswith(savePath,".gif") "set a savePath with the .gif extension"
    @assert isnothing(element2centralize) || element2centralize <= problem.model.assemblyData.nElementsTotal "element2centralize must be at most the number of elements"
    entries = _select_modes(problem, nModes, modes2plot, modeLabels, modalColorScheme, numberOfColors, nullDampingTol, showLegend, modeShapes=modeShapes, modeNumbers=modeNumbers)
    _set_backend(backend)

    layersOf = _mode_layers(problem, entries, plotSteady=plotSteady, colorSteady=colorSteady, steadyLabel=showLegend ? "Steady" : false, plotBCs=plotBCs, plotAeroSurf=plotAeroSurf, element2centralize=element2centralize)

    # Frames over a vibration cycle (without repeating the extremes), with the amplitude of each mode between -1 and 1 and, optionally, with the time scaled so that the modes oscillate at their relative frequencies
    frames = vcat(1:nFramesPerCycle,(nFramesPerCycle-1):-1:2)
    ratios = [matchModeFrequency ? entry.shape.frequency/entries[1].shape.frequency : NaN for entry in entries]
    function amplitudes(frame)
        ζ = 2(frame-1)/(nFramesPerCycle-1) - 1
        return [scale*(matchModeFrequency ? sin(ratio*asin(ζ)) : ζ) for ratio in ratios]
    end

    return _animate(problem.model, frames, frame -> layersOf(amplitudes(frame)), view=view, tolPlane=tolPlane, plotLimits=plotLimits, avoidDistortion=avoidDistortion, follow=:none, fps=fps, save=save, savePath=savePath, displayProgress=displayProgress, decorate=(plt,v,lims,frame) -> nothing, lw=lw, surfα=surfα, plotAxes=plotAxes, plotGrid=plotGrid, legendPos=showLegend ? legendPos : false, DPI=DPI, plotSize=plotSize, legendfontsize=legendFontSize)
end
export plot_mode_shapes_animation


## DYNAMIC PROBLEMS
# ------------------------------------------------------------------------------

# Number of time steps saved in the solution history of a dynamic problem
function _saved_steps(problem::DynamicProblem)
    nSaved = problem.solutionHistory.index - 1
    @assert nSaved > 0 "solve the problem with trackingTimeSteps=true before plotting"
    return nSaved
end

# Transformation from the basis A to the basis of the observer, at the time index ti. The observer is inertial and at rest (:none, :assembly) or attached to basis A, with its origin at the inertial origin (:basisA, so that the axes values read the position relative to it)
function _observer_xform(problem::DynamicProblem, ti, follow::Symbol)
    uA = SVector{3,Float64}(problem.model.kinematicsBasisA.u_AVec[ti])
    RA = _R_A(problem,ti)
    return follow == :basisA ? _Xform(RA'*uA,_I3) : _Xform(uA,RA)
end

# Layer with the configuration of a dynamic problem at its saved step k (of time index ti), seen by the inertial observer at rest or by the observer that follows the assembly
function _dynamic_layer(problem::DynamicProblem, geom::_Geometry, scales::_LoadScales, k; follow=:none, scale, color, ls=:solid, label=false, surfaces, plotBCs, plotDistLoads)
    ti = 1 + (k-1)*problem.trackingFrequency
    T = _observer_xform(problem, ti, follow)
    cfg = _configuration(geom, _poses(problem.solutionHistory.states[k],scale), T)
    loads = plotBCs || plotDistLoads ? _LoadSpec(problem, scales, k, ti, problem.solutionHistory.t[k], plotBCs, plotDistLoads) : nothing
    return _Layer(cfg, color=color, ls=ls, label=label, surfaces=surfaces, loads=loads)
end

"""
    plot_dynamic_deformation(problem::DynamicProblem; kwargs...)

Plots the animated deformation of the model in the given problem

# Arguments
- `problem::DynamicProblem`: problem

# Keyword arguments
- `backend::Symbol`: backend for plot
- `plotFrequency::Int`: frequency of time steps to plot
- `plotUndeformed::Bool`: flag to plot undeformed assembly
- `plotBCs::Bool`: flag to plot BCs
- `plotDistLoads::Bool`: flag to plot distributed loads (includes gravitational and aerodynamic loads)
- `plotAxes::Bool`: flag to plot axes
- `view::Union{Nothing,Tuple{Real,Real}}`: view angles (by default, the plot is planar if the assembly lies on a plane of the observer's basis over the whole animation, and 3D otherwise)
- `fps::Int`: frame rate for gif
- `scale::Real`: displacements and rotations scale
- `lw::Real`: linewidth
- `colorUndef`: color for undeformed assembly
- `colorDef`: color for deformed assembly
- `plotGrid::Bool`: flag to plot grid
- `legendPos`: legend position
- `tolPlane::Real`: tolerance of the extent of the assembly along an axis to plot as plane
- `plotAeroSurf::Bool`: flag to plot aerodynamic surfaces (this forces a 3D plot)
- `surfα::Float64`: transparency factor of aerodynamic surfaces
- `loadsSizeScaler::Real`: scaler for plot size of distributed and concentrated loads
- `plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}`: plot axis limits
- `follow::Symbol`: observer that the animation follows. With `:none`, the plot box is fixed and the motion is seen by an inertial observer at rest. With `:basisA`, the observer is attached to basis A and the plot box (of the size given by `plotLimits`) moves with it, while the axes values read the position relative to the inertial origin (resolved in basis A). With `:assembly`, the plot box is centered on the current position of the point of the assembly that initially coincides with the origin of basis A, while the axes values read the position seen by the inertial observer
- `DPI::Int`: pixel density
- `plotSize::Tuple{Real,Real}`: size of the plot (width, height) in pixels
- `avoidDistortion::Bool`: flag to avoid distortion in 3D plots
- `save::Bool`: flag to save the figure
- `savePath::String`: relative path (from the package's directory) on which to save the figure (must have the .gif extension)
- `showScale::Bool`: flag to show scale on plot
- `showTimeStamp::Bool`: flag to show time stamp on plot
- `scalePos::Vector{<:Real}`: position of scale on plot (relative to the plot limits)
- `timeStampPos::Vector{<:Real}`: position of time stamp on plot (relative to the plot limits)
- `displayProgress::Bool`: flag to display progress of gif creation
"""
function plot_dynamic_deformation(problem::DynamicProblem; backend::Symbol=:gr, plotFrequency::Int=1, plotUndeformed::Bool=false, plotBCs::Bool=true, plotDistLoads::Bool=true, plotAxes::Bool=true, view::Union{Nothing,Tuple{Real,Real}}=nothing, fps::Int=30, scale::Real=1, lw::Real=1, colorUndef=:black, colorDef=:blue, plotGrid::Bool=true, legendPos=:topright, tolPlane::Real=1e-8, plotAeroSurf::Bool=true, surfα::Float64=0.5, loadsSizeScaler::Real=1, plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}=nothing, follow::Symbol=:basisA, DPI::Int=300, plotSize::Tuple{Real,Real}=(600,400), avoidDistortion::Bool=false, save::Bool=false, savePath::String="/test/outputs/figures/fig.gif", showScale::Bool=true, showTimeStamp::Bool=true, scalePos::Vector{<:Real}=[0.1;0.05;0.05], timeStampPos::Vector{<:Real}=[0.5;0.05;0.05], displayProgress::Bool=false)

    _check_options(backend,DPI,plotLimits,surfα,plotSize, scale=scale, lw=lw, tolPlane=tolPlane, loadsSizeScaler=loadsSizeScaler, fps=fps)
    @assert plotFrequency >= 1 "plotFrequency must be at least 1"
    @assert follow in (:none,:basisA,:assembly) "follow must be :none, :basisA or :assembly"
    @assert endswith(savePath,".gif") "set a savePath with the .gif extension"
    @assert length(scalePos) == 3 && length(timeStampPos) == 3 "set scalePos and timeStampPos with 3 components"
    nSaved = _saved_steps(problem)
    _set_backend(backend)

    geom = _geometry(problem.model, surfaces=plotAeroSurf)
    scales = _load_scales(problem, geom, loadsSizeScaler, 1:nSaved, BCs=plotBCs, dist=plotDistLoads)
    undeformed(T) = _Layer(_configuration(geom, fill(_NO_POSE,length(geom.r1)), T), color=colorUndef, label="Undeformed", surfaces=plotAeroSurf)
    function layersOf(k)
        deformed = _dynamic_layer(problem, geom, scales, k, follow=follow, scale=scale, color=colorDef, label=plotUndeformed ? "Deformed" : false, surfaces=plotAeroSurf, plotBCs=plotBCs, plotDistLoads=plotDistLoads)
        # The undeformed assembly stays at the origin of the inertial basis, unless the observer follows the assembly
        return plotUndeformed ? [undeformed(follow == :none ? _Xform() : deformed.cfg.T),deformed] : [deformed]
    end

    # Time stamp and scale: in the title of 3D plots and as annotations in planar plots (which may be positioned above the axes, thus the top margin)
    function decorate(plt, v, lims, k)
        t = round(problem.solutionHistory.t[k], sigdigits=3)
        if _is3D(v)
            lines = String[]
            showTimeStamp && push!(lines, string("Time = ",t," s"))
            showScale && push!(lines, string("Scale = ",scale,"×"))
            isempty(lines) || title!(plt, join(lines,"\n"))
        else
            showTimeStamp && annotate!(plt, _note_position(v,lims,timeStampPos)..., text(string("Time: ",t," s"),8))
            showScale && annotate!(plt, _note_position(v,lims,scalePos)..., text(string("Scale: ",scale,"×"),8))
        end
        return nothing
    end
    margin = (showTimeStamp || showScale) ? (top_margin=10mm,) : (;)

    return _animate(problem.model, 1:plotFrequency:nSaved, layersOf, view=view, tolPlane=tolPlane, plotLimits=plotLimits, avoidDistortion=avoidDistortion, follow=follow, fps=fps, save=save, savePath=savePath, displayProgress=displayProgress, decorate=decorate, lw=lw, surfα=surfα, plotAxes=plotAxes, plotGrid=plotGrid, legendPos=plotUndeformed ? legendPos : false, DPI=DPI, plotSize=plotSize; margin...)
end
export plot_dynamic_deformation

"""
    plot_dynamic_deformations(problems::Array{DynamicProblem}; kwargs...)

Plots the animated deformation of the models in the given problems

# Arguments
- `problems::Array{DynamicProblem}`: problems

# Keyword arguments
- `backend::Symbol`: backend for plot
- `plotFrequency::Int`: frequency of time steps to plot (counted over the time steps saved by the problems)
- `timeTol::Real`: tolerance to identify time step
- `plotBCs::Bool`: flag to plot BCs
- `plotDistLoads::Bool`: flag to plot distributed loads (includes gravitational and aerodynamic loads)
- `plotAxes::Bool`: flag to plot axes
- `view::Union{Nothing,Tuple{Real,Real}}`: view angles (by default, the plot is planar if the assemblies lie on a plane of the observer's basis over the whole animation, and 3D otherwise)
- `fps::Int`: frame rate for gif
- `lw::Real`: linewidth
- `colorScheme::Symbol`: colors for assemblies in problems
- `plotGrid::Bool`: flag to plot grid
- `legendPos`: legend position
- `legendEntries::Union{Nothing,Array{String}}`: legend entries for problems (one for each problem)
- `tolPlane::Real`: tolerance of the extent of the assemblies along an axis to plot as plane
- `plotAeroSurf::Bool`: flag to plot aerodynamic surfaces (this forces a 3D plot)
- `surfα::Float64`: transparency factor of aerodynamic surfaces
- `loadsSizeScaler::Real`: scaler for plot size of distributed and concentrated loads
- `plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}`: plot axis limits
- `follow::Symbol`: observer that the animation follows. With `:none`, the plot box is fixed and the motion is seen by an inertial observer at rest. With `:basisA`, the observer is attached to basis A and the plot box (of the size given by `plotLimits`) moves with it, while the axes values read the position relative to the inertial origin (resolved in basis A). With `:assembly`, the plot box is centered on the current position of the point of the (last) assembly that initially coincides with the origin of basis A, while the axes values read the position seen by the inertial observer
- `DPI::Int`: pixel density
- `plotSize::Tuple{Real,Real}`: size of the plot (width, height) in pixels
- `avoidDistortion::Bool`: flag to avoid distortion in 3D plots
- `save::Bool`: flag to save the figure
- `savePath::String`: relative path (from the package's directory) on which to save the figure (must have the .gif extension)
- `showTimeStamp::Bool`: flag to show time stamp on plot
- `displayProgress::Bool`: flag to display progress of gif creation
"""
function plot_dynamic_deformations(problems::Array{DynamicProblem}; backend::Symbol=:gr, plotFrequency::Int=1, timeTol::Real=1e-8, plotBCs::Bool=true, plotDistLoads::Bool=true, plotAxes::Bool=true, view::Union{Nothing,Tuple{Real,Real}}=nothing, fps::Int=30, lw::Real=1, colorScheme::Symbol=:rainbow, plotGrid::Bool=true, legendPos=:topright, legendEntries::Union{Nothing,Array{String}}=nothing, tolPlane::Real=1e-8, plotAeroSurf::Bool=true, surfα::Float64=0.5, loadsSizeScaler::Real=1, plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}=nothing, follow::Symbol=:basisA, DPI::Int=300, plotSize::Tuple{Real,Real}=(600,400), avoidDistortion::Bool=false, save::Bool=false, savePath::String="/test/outputs/figures/fig.gif", showTimeStamp::Bool=true, displayProgress::Bool=false)

    @assert !isempty(problems) "set at least one problem"
    _check_options(backend,DPI,plotLimits,surfα,plotSize, lw=lw, tolPlane=tolPlane, loadsSizeScaler=loadsSizeScaler, fps=fps, timeTol=timeTol)
    @assert plotFrequency >= 1 "plotFrequency must be at least 1"
    @assert follow in (:none,:basisA,:assembly) "follow must be :none, :basisA or :assembly"
    @assert endswith(savePath,".gif") "set a savePath with the .gif extension"
    @assert isnothing(legendEntries) || length(legendEntries) == length(problems) "specify legendEntries with the same size as problems"
    nSaved = [_saved_steps(problem) for problem in problems]
    @assert all(nSaved .>= 2) "the problems must have at least two saved time steps"
    _set_backend(backend)

    # Time grid, with the step of the (saved) time vectors of the problems, and the saved step of each problem at each instant of the grid (zero if none)
    times = [problem.solutionHistory.t[1:n] for (problem,n) in zip(problems,nSaved)]
    Δt = plotFrequency*gcd_floats_approx([t[2]-t[1] for t in times])
    grid = range(minimum(first.(times)), maximum(last.(times)), step=Δt)
    function savedStep(t, τ)
        k = searchsortedfirst(t, τ-timeTol)
        return k <= length(t) && t[k] <= τ+timeTol ? k : 0
    end
    steps = [[savedStep(t,τ) for t in times] for τ in grid]
    frames = [i for i in eachindex(grid) if any(>(0),steps[i])]
    @assert !isempty(frames) "the problems have no saved time step in common with the time grid; increase timeTol"

    colors = cgrad(colorScheme, length(problems), categorical=true)
    geoms = [_geometry(problem.model, surfaces=plotAeroSurf) for problem in problems]
    scales = [_load_scales(problem, geom, loadsSizeScaler, 1:n, BCs=plotBCs, dist=plotDistLoads) for (problem,geom,n) in zip(problems,geoms,nSaved)]
    layersOf(i) = [_dynamic_layer(problems[p], geoms[p], scales[p], k, follow=follow, scale=1, color=colors[p], surfaces=plotAeroSurf, plotBCs=plotBCs, plotDistLoads=plotDistLoads) for (p,k) in enumerate(steps[i]) if k > 0]

    # Legend entries of the problems, which are always present, and time stamp
    function decorate(plt, v, lims, i)
        isnothing(legendEntries) || foreach(p -> _legend_entry!(plt, v, colors[p], lw, legendEntries[p]), eachindex(problems))
        showTimeStamp && title!(plt, string("\$t=",round(grid[i],sigdigits=6),"\$ s"))
        return nothing
    end

    return _animate(problems[1].model, frames, layersOf, view=view, tolPlane=tolPlane, plotLimits=plotLimits, avoidDistortion=avoidDistortion, follow=follow, fps=fps, save=save, savePath=savePath, displayProgress=displayProgress, decorate=decorate, lw=lw, surfα=surfα, plotAxes=plotAxes, plotGrid=plotGrid, legendPos=isnothing(legendEntries) ? false : legendPos, DPI=DPI, plotSize=plotSize, top_margin=10mm)
end
export plot_dynamic_deformations

"""
    plot_snapshots(problem::DynamicProblem; kwargs...)

Plots the deformed states for the model in the given dynamic problem at specified time instances (snapshots)

# Arguments
- `problem::DynamicProblem`: problem

# Keyword arguments
- `backend::Symbol`: backend for plot
- `plotBCs::Bool`: flag to plot BCs
- `plotDistLoads::Bool`: flag to plot distributed loads (includes gravitational and aerodynamic loads)
- `view::Union{Nothing,Tuple{Real,Real}}`: view angles (by default, the plot is planar if the assembly lies on a plane of the inertial basis over all snapshots, and 3D otherwise)
- `scale::Real`: displacements and rotations scale
- `lw::Real`: linewidth
- `ls::Symbol`: line style
- `color`: color
- `plotAxes::Bool`: flag to plot axes
- `plotGrid::Bool`: flag to plot grid
- `tolPlane::Real`: tolerance of the extent of the assembly along an axis to plot as plane
- `plotAeroSurf::Bool`: flag to plot aerodynamic surfaces (this forces a 3D plot)
- `surfα::Float64`: transparency factor of aerodynamic surfaces
- `snapshots::Vector{<:Real}`: (tentative) time instances to plot, in ascending order
- `plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}`: plot axis limits
- `showScale::Bool`: flag to show scale on plot
- `scalePos::Vector{<:Real}`: position of scale on plot (in the coordinates of the axes of planar plots)
- `DPI::Int`: pixel density
- `plotSize::Tuple{Real,Real}`: size of the plot (width, height) in pixels
- `avoidDistortion::Bool`: flag to avoid distortion in 3D plots
- `save::Bool`: flag to save the figure
- `savePath::String`: relative path (from the package's directory) on which to save the figure
"""
function plot_snapshots(problem::DynamicProblem; backend::Symbol=:gr, plotBCs::Bool=true, plotDistLoads::Bool=true, view::Union{Nothing,Tuple{Real,Real}}=nothing, scale::Real=1, lw::Real=1, ls::Symbol=:dash, color=:black, plotAxes::Bool=true, plotGrid::Bool=true, tolPlane::Real=1e-8, plotAeroSurf::Bool=true, surfα::Float64=0.5, snapshots::Vector{<:Real}, plotLimits::Union{Nothing,Tuple{Vector{<:Real},Vector{<:Real},Vector{<:Real}}}=nothing, showScale::Bool=false, scalePos::Vector{<:Real}=[0,0], DPI::Int=300, plotSize::Tuple{Real,Real}=(600,400), avoidDistortion::Bool=false, save::Bool=false, savePath::String="/test/outputs/figures/fig.pdf")

    _check_options(backend,DPI,plotLimits,surfα,plotSize, scale=scale, lw=lw, tolPlane=tolPlane)
    nSaved = _saved_steps(problem)
    t = problem.solutionHistory.t[1:nSaved]
    @assert !isempty(snapshots) && issorted(snapshots) "set the snapshots in ascending order"
    @assert snapshots[1] >= t[1] && snapshots[end] <= t[end] "the snapshots must be within the time interval of the solution: [$(t[1]),$(t[end])]"
    _set_backend(backend)

    geom = _geometry(problem.model, surfaces=plotAeroSurf)
    scales = _load_scales(problem, geom, 1.0, 1:nSaved, BCs=plotBCs, dist=plotDistLoads)
    layers = [_dynamic_layer(problem, geom, scales, argmin(abs.(t .- snapshot)), scale=scale, color=color, ls=ls, surfaces=plotAeroSurf, plotBCs=plotBCs, plotDistLoads=plotDistLoads) for snapshot in snapshots]

    decorate(plt, v) = showScale && (_is3D(v) ? title!(plt, string("Scale = ",scale,"×")) : annotate!(plt, scalePos[1], scalePos[2], text(string("Deformation scale: ",scale,"×"),8)))
    plt = _plot_frame(problem.model, layers, view=view, tolPlane=tolPlane, plotLimits=plotLimits, avoidDistortion=avoidDistortion, decorate=decorate, lw=lw, surfα=surfα, plotAxes=plotAxes, plotGrid=plotGrid, legendPos=false, DPI=DPI, plotSize=plotSize)
    save && _save_figure(plt,savePath)

    return plt
end
export plot_snapshots

## OUTPUTS OVER THE ARCLENGTH AND OVER TIME
# ------------------------------------------------------------------------------

# Output quantity: label of the vertical axis (# stands for the component), its unit and conversion factor (functions of the units of the model), whether it is aerodynamic (a scalar of the elements with aerodynamic surface; the others have three components), its value at the midpoint of an element (from the element's state and aerodynamic variables) and, for nodal quantities, at a local node of an element
struct _Quantity
    label::String
    unit::Function
    factor::Function
    aero::Bool
    value::Function
    nodal::Union{Nothing,Function}
end

_no_unit(units) = " "
_no_factor(units) = 1.0
_angle_factor(units) = units.angle == "deg" ? 180/π : 1.0
_nodal_value(field1, field2) = (state, node) -> getfield(state.nodal, node == 1 ? field1 : field2)
_lift(a) = a.aeroCoefficients.cn*cos(a.unsteadyKin.αₑ) + a.aeroCoefficients.ct*sin(a.unsteadyKin.αₑ)
_drag(a) = a.aeroCoefficients.cn*sin(a.unsteadyKin.αₑ) - a.aeroCoefficients.ct*cos(a.unsteadyKin.αₑ)

const _QUANTITIES = Dict{String,_Quantity}(
    "u"    => _Quantity("u_#",units -> units.length, _no_factor, false, (s,a) -> s.midpoint.u, _nodal_value(:u1,:u2)),
    "p"    => _Quantity("p_#", _no_unit, _no_factor, false, (s,a) -> s.midpoint.p, _nodal_value(:p1,:p2)),
    "F"    => _Quantity("F_#^*", units -> units.force, _no_factor, false, (s,a) -> s.midpoint.F, _nodal_value(:F1,:F2)),
    "M"    => _Quantity("M_#^*", units -> string(units.force,".",units.length), _no_factor,  false, (s,a) -> s.midpoint.M, _nodal_value(:M1,:M2)),
    "V"    => _Quantity("V_#^*", units -> string(units.length,"/s"), _no_factor, false, (s,a) -> s.midpoint.V, nothing),
    "Ω"    => _Quantity("\\Omega_#^*", units -> string(units.angle,"/s"), _angle_factor, false, (s,a) -> s.midpoint.Ω, nothing),
    "γ"    => _Quantity("\\gamma_#^+", _no_unit, _no_factor, false, (s,a) -> s.midpoint.γ, nothing),
    "κ"    => _Quantity("\\kappa_#^+", _no_unit, _no_factor, false, (s,a) -> s.midpoint.κ, nothing),
    "P"    => _Quantity("P_#^*", units -> string(units.mass,"/s"),  _no_factor,    false, (s,a) -> s.midpoint.P, nothing),
    "H"    => _Quantity("H_#^*", units -> string(units.mass,".",units.length,"/s"), _no_factor, false, (s,a) -> s.midpoint.H, nothing),
    "udot" => _Quantity("\\dot{u}_#", units -> string(units.length,"/s"),_no_factor,    false, (s,a) -> s.rates.udot, nothing),
    "pdot" => _Quantity("\\dot{p}_#", units -> "1/s",_no_factor, false, (s,a) -> s.rates.pdot, nothing),
    "Vdot" => _Quantity("\\dot{V}_#", units -> string(units.length,"\$/s^2\$"), _no_factor, false, (s,a) -> s.rates.Vdot, nothing),
    "Ωdot" => _Quantity("\\dot{\\Omega}_#", units -> string(units.angle,"\$/s^2\$"), _angle_factor, false, (s,a) -> s.rates.Ωdot, nothing),
    "α"    => _Quantity("\\alpha",units -> "deg", units -> 180/π, true, (s,a) -> a.unsteadyKin.αₑ, nothing),
    "cn"   => _Quantity("c_n", _no_unit, _no_factor, true, (s,a) -> a.aeroCoefficients.cn, nothing),
    "cm"   => _Quantity("c_m", _no_unit, _no_factor, true, (s,a) -> a.aeroCoefficients.cm, nothing),
    "ct"   => _Quantity("c_t", _no_unit, _no_factor, true, (s,a) -> a.aeroCoefficients.ct, nothing),
    "cl"   => _Quantity("c_l", _no_unit, _no_factor, true, (s,a) -> _lift(a), nothing),
    "cd"   => _Quantity("c_d", _no_unit, _no_factor, true, (s,a) -> _drag(a), nothing),
)

const _STEADY_OUTPUTS = ["u","p","F","M","V","Ω","γ","κ","P","H","α","cn","cm","ct","cl","cd"]
const _TIME_NODAL_OUTPUTS = ["u","p","F","M"]
const _TIME_ELEMENTAL_OUTPUTS = vcat(_STEADY_OUTPUTS, ["udot","pdot","Vdot","Ωdot"])

# Valid names of outputs: the quantities and, for the ones with three components, each of their components
function _valid_outputs(names)
    valid = String[]
    for name in names
    push!(valid, name)
    _QUANTITIES[name].aero || append!(valid, string.(name,1:3))
    end
    return valid
end

# (name,component) of the requested outputs, in order and without repetitions (the component is 0 for aerodynamic quantities)
function _expand_outputs(outputs)
    expanded = Tuple{String,Int}[]
    for out in outputs
        if haskey(_QUANTITIES,out)
            _QUANTITIES[out].aero ? push!(expanded,(out,0)) : append!(expanded,[(out,c) for c in 1:3])
        else
            push!(expanded,(chop(out),parse(Int,last(out))))
        end
    end
    return unique(expanded)
end

_output_file(name, c) = c == 0 ? name : string(name,c)
_output_ylabel(q::_Quantity, c, units) = string("\$", c == 0 ? q.label : replace(q.label,"#"=>string(c)), "\$ [", q.unit(units), "]")
_component(value, c) = c == 0 ? value : value[c]

# File names of the figures of the outputs, which must differ even if the file system is case insensitive (e.g. p1 and P1)
function _figure_files(outputs)
    files = String[]
    for (name,c) in outputs
        file = _output_file(name,c)
        while lowercase(file) in lowercase.(files)
            file = string(file,"_")
        end
        push!(files, file)
    end
    return files
end

# Shows a figure of an output and, optionally, saves it in the folder (relative to the package's directory)
function _finish_output_figure(plt, file, save, saveFolder, figureExtension)
    display(plt)
    save && _save_figure(plt, joinpath(saveFolder, string(file,figureExtension)))
    return nothing
end

# Arclength positions and values of an output (component c) at the last load step, along a beam: at the nodes of its elements for nodal quantities, at their midpoints otherwise (only for the elements with aerodynamic surface in the case of aerodynamic quantities)
function _beam_profile(problem, beam, q::_Quantity, c)
    elements = problem.model.assemblyData.elements
    states = problem.loadstepHistory.states[end]
    aeroVars = problem.loadstepHistory.aeroVariables[end]
    x = Float64[]
    y = Float64[]
    for e in beam.assemblyData.elementsGlobalIDs
        geometry = elements[e].geometry
        if q.aero
            isnothing(aeroVars[e]) && continue
            push!(x, geometry.x1)
            push!(y, q.value(states[e],aeroVars[e]))
        elseif isnothing(q.nodal)
            push!(x, geometry.x1)
            push!(y, q.value(states[e],aeroVars[e])[c])
        else
            push!(x, geometry.x1_n1, geometry.x1_n2)
            push!(y, q.nodal(states[e],1)[c], q.nodal(states[e],2)[c])
        end
    end
    return x,y
end

"""
    plot_steady_outputs(problem::Problem; kwargs...)

Plots outputs of a steady problem, over the arclength of the beams, at its last load step

# Arguments
- `problem::Problem`: problem

# Keyword arguments
- `backend::Symbol`: backend for plot
- `outputs::Vector{String}`: list of outputs: `"u"`, `"p"`, `"F"`, `"M"`, `"V"`, `"Ω"`, `"γ"`, `"κ"`, `"P"`, `"H"` (all three components, or a single one, e.g. `"F3"`) and, for the elements with aerodynamic surface, `"α"`, `"cn"`, `"cm"`, `"ct"`, `"cl"`, `"cd"`
- `beamGroups`: list of beams in each group for which arclengths are concatenated
- `lw::Real`: linewidth
- `colorScheme`: color scheme
- `legendPos`: legend position
- `save::Bool`: flag to save figures
- `DPI::Int`: pixel density
- `saveFolder::String`: relative path (from the package's directory) of the folder where to save figures
- `figureExtension::String`: figure extension
"""
function plot_steady_outputs(problem::Problem; backend::Symbol=:gr, outputs::Vector{String}=["u","p","F","M","V","Ω","α","cn","cm","ct","cl","cd"], beamGroups=1:length(problem.model.beams), lw::Real=1, colorScheme=:rainbow, legendPos=:best, save::Bool=false, DPI::Int=300, saveFolder::String="/test/outputs/figures/", figureExtension::String=".pdf")

    _check_steady_problem(problem)
    _check_options(backend,DPI,lw=lw)
    valid = _valid_outputs(_STEADY_OUTPUTS)
    @assert all(out -> out in valid, outputs) "set outputs as one of $(join(valid, ','))"
    nBeams = length(problem.model.beams)
    @assert !isempty(beamGroups) && all(b -> 1 <= b <= nBeams, vcat(beamGroups...)) "set beamGroups with beam IDs between 1 and $nBeams"
    _set_backend(backend)

    beams = problem.model.beams
    units = problem.model.units
    multipleGroups = length(beamGroups) > 1
    colors = palette(colorScheme, length(beamGroups))
    requested = _expand_outputs(outputs)

    for ((name,c),file) in zip(requested, _figure_files(requested))
        q = _QUANTITIES[name]
        factor = q.factor(units)
        plt = plot(xlabel=string("\$x_1\$ [",units.length,"]"), ylabel=_output_ylabel(q,c,units), legend=multipleGroups ? legendPos : false, dpi=DPI)
        hasData = false
        for (i,group) in enumerate(beamGroups)
            x = Float64[]
            y = Float64[]
            offset = 0.0
            for b in group
                xb,yb = _beam_profile(problem, beams[b], q, c)
                append!(x, xb .+ offset)
                append!(y, factor .* yb)
                offset += beams[b].geometry.ℓ
            end
            isempty(x) && continue
            plot!(plt, x, y, lw=lw, c=colors[i], label=multipleGroups ? string("Beam group ",i) : false)
            hasData = true
        end
        hasData && _finish_output_figure(plt, file, save, saveFolder, figureExtension)
    end
    return nothing
end
export plot_steady_outputs

"""
    plot_time_outputs(problem::DynamicProblem; kwargs...)

Plots outputs of a dynamic problem, over time

# Arguments
- `problem::DynamicProblem`: problem

# Keyword arguments
- `backend::Symbol`: backend for plot
- `nodes::Vector{Tuple{Int,Int}}`: nodes for which to plot, each given by the global ID of an element and the local node (1 or 2) of that element
- `elements::Vector{Int}`: global IDs of elements for which to plot
- `nodalOutputs::Vector{String}`: nodal outputs to plot: `"u"`, `"p"`, `"F"`, `"M"` (all three components, or a single one, e.g. `"u2"`)
- `elementalOutputs::Vector{String}`: elemental outputs to plot: `"u"`, `"p"`, `"F"`, `"M"`, `"V"`, `"Ω"`, `"γ"`, `"κ"`, `"P"`, `"H"`, `"udot"`, `"pdot"`, `"Vdot"`, `"Ωdot"` (all three components, or a single one, e.g. `"F3"`) and, for the elements with aerodynamic surface, `"α"`, `"cn"`, `"cm"`, `"ct"`, `"cl"`, `"cd"`
- `lw::Real`: linewidth
- `colorScheme`: color scheme
- `showLegend::Bool`: flag to show legend
- `legendPos`: legend position
- `DPI::Int`: pixel density
- `save::Bool`: flag to save figures
- `saveFolder::String`: relative path (from the package's directory) of the folder where to save figures
- `figureExtension::String`: figure extension
"""
function plot_time_outputs(problem::DynamicProblem; backend::Symbol=:gr, nodes::Vector{Tuple{Int,Int}}=Vector{Tuple{Int,Int}}(), elements::Vector{Int}=Vector{Int}(), nodalOutputs::Vector{String}=["u","p","F","M"], elementalOutputs::Vector{String}=["u","p","F","M","V","Ω","α","cn","cm","ct","cl","cd"], lw::Real=1, colorScheme=:rainbow, showLegend::Bool=true, legendPos=:best, DPI::Int=300, save::Bool=false, saveFolder::String="/test/outputs/figures/", figureExtension::String=".pdf")

    _check_options(backend,DPI,lw=lw)
    validNodal = _valid_outputs(_TIME_NODAL_OUTPUTS)
    validElemental = _valid_outputs(_TIME_ELEMENTAL_OUTPUTS)
    @assert all(out -> out in validNodal, nodalOutputs) "set 'nodalOutputs' as one of $(join(validNodal, ','))"
    @assert all(out -> out in validElemental, elementalOutputs) "set 'elementalOutputs' as one of $(join(validElemental, ','))"
    nElements = problem.model.assemblyData.nElementsTotal
    @assert all(e -> 0 < e <= nElements, elements) "set 'elements' between 1 and the total number of elements in the model"
    @assert all(n -> 0 < n[1] <= nElements, nodes) "element of 'nodes' must be between 1 and the total number of elements in the model"
    @assert all(n -> n[2] in (1,2), nodes) "set 'nodes' as a tuple in which the first value is the element and the second is its local node (1 or 2)"
    nSaved = _saved_steps(problem)
    _set_backend(backend)

    units = problem.model.units
    modelElements = problem.model.assemblyData.elements
    time = problem.solutionHistory.t[1:nSaved]
    states = problem.solutionHistory.states
    aeroVars = problem.solutionHistory.aeroVariables

    # Series to plot (elements first, then nodes) and the outputs of each kind
    series = vcat([(string("Element ",e),e,0) for e in elements],
    [(string("Node ",modelElements[e].assemblyData.nodesGlobalIDs[n]),e,n) for (e,n) in nodes])
    elementalPairs = isempty(elements) ? Tuple{String,Int}[] : _expand_outputs(elementalOutputs)
    nodalPairs = isempty(nodes) ? Tuple{String,Int}[] : _expand_outputs(nodalOutputs)
    colors = palette(colorScheme, max(length(series),1))
    requested = unique(vcat(elementalPairs,nodalPairs))

    for ((name,c),file) in zip(requested, _figure_files(requested))
        q = _QUANTITIES[name]
        factor = q.factor(units)
        plt = plot(xlabel="\$t\$ [s]", ylabel=_output_ylabel(q,c,units), legend=showLegend ? legendPos : false, dpi=DPI)
        hasData = false
        for (i,(label,e,n)) in enumerate(series)
            (name,c) in (n == 0 ? elementalPairs : nodalPairs) || continue
            q.aero && isnothing(aeroVars[1][e]) && continue
            y = [factor*_component(n == 0 ? q.value(states[k][e],aeroVars[k][e]) : q.nodal(states[k][e],n), c) for k in 1:nSaved]
            plot!(plt, time, y, lw=lw, c=colors[i], label=showLegend ? label : false)
            hasData = true
        end
        hasData && _finish_output_figure(plt, file, save, saveFolder, figureExtension)
    end
    return nothing
end
export plot_time_outputs
