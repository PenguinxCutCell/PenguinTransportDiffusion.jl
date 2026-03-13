module PenguinTransportDiffusion

using LinearAlgebra
using SparseArrays
using StaticArrays

using CartesianOperators
using PenguinBCs
using PenguinSolverCore

using PenguinDiffusion
using PenguinTransport

export AdvDiffModelMono, AdvDiffModelDiph
export MovingAdvDiffModelMono, MovingAdvDiffModelDiph
export AdvDiffCoupledModelMono
export assemble_steady_mono!, assemble_unsteady_mono!
export assemble_steady_diph!, assemble_unsteady_diph!
export assemble_unsteady_mono_moving!, assemble_unsteady_diph_moving!
export solve_steady!, solve_unsteady!, solve_unsteady_moving!
export rebuild!, update_advection_ops!

mutable struct AdvDiffModelMono{N,T,DM,VTω,VTγ,BCT,SCH}
    diff::DM
    uω::VTω
    uγ::VTγ
    bc::BCT
    scheme::SCH
    periodic::NTuple{N,Bool}
end

mutable struct AdvDiffModelDiph{N,T,DM,V1ω,V1γ,V2ω,V2γ,BCT,SCH}
    diff::DM
    u1ω::V1ω
    u1γ::V1γ
    u2ω::V2ω
    u2γ::V2γ
    bc::BCT
    scheme::SCH
    periodic::NTuple{N,Bool}
end

mutable struct MovingAdvDiffModelMono{N,T,DM,VTω,VTγ,WTγ,BCT,SCH}
    diff::DM
    uω::VTω
    uγ::VTγ
    wγ::WTγ
    bc::BCT
    scheme::SCH
    periodic::NTuple{N,Bool}
end

mutable struct MovingAdvDiffModelDiph{N,T,DM,V1ω,V1γ,V2ω,V2γ,WTγ,BCT,SCH}
    diff::DM
    u1ω::V1ω
    u1γ::V1γ
    u2ω::V2ω
    u2γ::V2γ
    wγ::WTγ
    bc::BCT
    scheme::SCH
    periodic::NTuple{N,Bool}
end

function _as_velocity_components(u, N::Int)
    if u isa Tuple
        length(u) == N || throw(ArgumentError("velocity tuple must contain $N components"))
        return u
    elseif u isa AbstractVector
        length(u) == N || throw(ArgumentError("velocity container must contain $N components"))
        return Tuple(u)
    end
    throw(ArgumentError("velocity field must be a tuple/vector of $N components"))
end

function _eval_fun_or_const(v, x::SVector{N,T}, t::T) where {N,T}
    if v isa Number
        return convert(T, v)
    elseif v isa Ref
        return _eval_fun_or_const(v[], x, t)
    elseif v isa Function
        if applicable(v, x..., t)
            return convert(T, v(x..., t))
        elseif applicable(v, x...)
            return convert(T, v(x...))
        end
    end
    throw(ArgumentError("callback/value must be numeric, Ref, (x...), or (x..., t)"))
end

function _component_values(cap::AssembledCapacity{N,T}, c, xpts::Vector{SVector{N,T}}, t::T) where {N,T}
    nt = cap.ntotal
    if c isa AbstractVector
        length(c) == nt || throw(ArgumentError("velocity vector length must be $nt"))
        out = Vector{T}(undef, nt)
        @inbounds for i in 1:nt
            xi = xpts[i]
            if all(isfinite, xi)
                vi = convert(T, c[i])
                out[i] = isfinite(vi) ? vi : zero(T)
            else
                out[i] = zero(T)
            end
        end
        return out
    elseif c isa Number || c isa Ref || c isa Function
        out = Vector{T}(undef, nt)
        @inbounds for i in 1:nt
            xi = xpts[i]
            out[i] = all(isfinite, xi) ? _eval_fun_or_const(c, xi, t) : zero(T)
        end
        return out
    end
    throw(ArgumentError("unsupported velocity component type $(typeof(c))"))
end

function _velocity_tuple_values(cap::AssembledCapacity{N,T}, u, xpts::Vector{SVector{N,T}}, t::T) where {N,T}
    comps = _as_velocity_components(u, N)
    return ntuple(d -> _component_values(cap, comps[d], xpts, t), N)
end

function _velocity_values(cap::AssembledCapacity{N,T}, uω, uγ, t::T) where {N,T}
    return _velocity_tuple_values(cap, uω, cap.C_ω, t), _velocity_tuple_values(cap, uγ, cap.C_γ, t)
end

function _split_diph_source(source, ::Type{T}) where {T}
    if source isa Tuple && length(source) == 2
        return source[1], source[2]
    elseif source isa Function
        s1 = (args...) -> begin
            s = applicable(source, args...) ? source(args...) : source(args[1:(end - 1)]...)
            s[1]
        end
        s2 = (args...) -> begin
            s = applicable(source, args...) ? source(args...) : source(args[1:(end - 1)]...)
            s[2]
        end
        return s1, s2
    end
    throw(ArgumentError("diph source must be a function returning a tuple or a tuple of two callbacks/constants"))
end

function _value_time_dependent(v, x::SVector{N,T}) where {N,T}
    if v isa Ref
        return _value_time_dependent(v[], x)
    end
    return v isa Function && applicable(v, x..., zero(T))
end

function _vel_time_dependent(u, x::SVector{N,T}) where {N,T}
    comps = _as_velocity_components(u, N)
    for c in comps
        _value_time_dependent(c, x) && return true
    end
    return false
end

function _source_time_dependent(source, x::SVector{N,T}) where {N,T}
    return source isa Function && applicable(source, x..., zero(T))
end

function _adv_border_time_dependent(bc::BorderConditions, x::SVector{N,T}) where {N,T}
    for side_bc in values(bc.borders)
        if side_bc isa Inflow
            _value_time_dependent(side_bc.value, x) && return true
        end
    end
    return false
end

function _side_names(N::Int)
    if N == 1
        return (:left, :right)
    elseif N == 2
        return (:left, :right, :bottom, :top)
    elseif N == 3
        return (:left, :right, :bottom, :top, :backward, :forward)
    end
    throw(ArgumentError("unsupported dimension N=$N; expected 1, 2, or 3"))
end

function _split_borderconditions(bc::BorderConditions, N::Int)
    bdiff = Dict{Symbol,AbstractBoundary}()
    badv = Dict{Symbol,AbstractBoundary}()
    for side in _side_names(N)
        side_bc = get(bc.borders, side, nothing)
        if side_bc === nothing
            continue
        elseif side_bc isa Periodic
            bdiff[side] = side_bc
            badv[side] = side_bc
        elseif side_bc isa Dirichlet || side_bc isa Neumann || side_bc isa Robin
            bdiff[side] = side_bc
            badv[side] = Outflow()
        elseif side_bc isa Inflow || side_bc isa Outflow
            bdiff[side] = Neumann(0.0)
            badv[side] = side_bc
        else
            throw(ArgumentError("unsupported mixed BC type $(typeof(side_bc)) on side $side"))
        end
    end
    return BorderConditions(; bdiff...), BorderConditions(; badv...)
end

function _insert_block!(A::SparseMatrixCSC{T,Int}, rows::UnitRange{Int}, cols::UnitRange{Int}, B::SparseMatrixCSC{T,Int}) where {T}
    size(B, 1) == length(rows) || throw(DimensionMismatch("block rows do not match target range"))
    size(B, 2) == length(cols) || throw(DimensionMismatch("block cols do not match target range"))
    @inbounds for j in 1:size(B, 2)
        for p in nzrange(B, j)
            i = B.rowval[p]
            A[rows[i], cols[j]] = A[rows[i], cols[j]] + B.nzval[p]
        end
    end
    return A
end

function _insert_vec!(b::Vector{T}, rows::UnitRange{Int}, v::Vector{T}) where {T}
    length(v) == length(rows) || throw(DimensionMismatch("vector block length mismatch"))
    @inbounds for i in eachindex(v)
        b[rows[i]] += v[i]
    end
    return b
end

function _scale_rows!(A::SparseMatrixCSC{T,Int}, rows::UnitRange{Int}, α::T) where {T}
    α == one(T) && return A
    r1 = first(rows)
    r2 = last(rows)
    @inbounds for j in 1:size(A, 2)
        for p in nzrange(A, j)
            i = A.rowval[p]
            if r1 <= i <= r2
                A.nzval[p] *= α
            end
        end
    end
    return A
end

function _theta_from_scheme(::Type{T}, scheme) where {T}
    if scheme isa Symbol
        if scheme === :BE
            return one(T)
        elseif scheme === :CN
            return convert(T, 0.5)
        end
        throw(ArgumentError("unknown scheme `$scheme`; expected :BE or :CN"))
    elseif scheme isa Real
        θ = convert(T, scheme)
        (zero(T) <= θ <= one(T)) || throw(ArgumentError("numeric θ must satisfy 0 ≤ θ ≤ 1"))
        return θ
    end
    throw(ArgumentError("scheme must be a Symbol (:BE/:CN) or a numeric theta"))
end

function _validate_mono_layout(cap::AssembledCapacity, lay)
    nt = cap.ntotal
    length(lay.ω) == nt || throw(ArgumentError("layout ω length ($(length(lay.ω))) must match cap.ntotal ($nt)"))
    length(lay.γ) == nt || throw(ArgumentError("layout γ length ($(length(lay.γ))) must match cap.ntotal ($nt)"))
    return nothing
end

function _validate_diph_layout(cap1::AssembledCapacity, cap2::AssembledCapacity, lay)
    nt1 = cap1.ntotal
    nt2 = cap2.ntotal
    nt1 == nt2 || throw(ArgumentError("cap1/cap2 ntotal mismatch"))
    length(lay.ω1) == nt1 || throw(ArgumentError("layout ω1 length mismatch"))
    length(lay.γ1) == nt1 || throw(ArgumentError("layout γ1 length mismatch"))
    length(lay.ω2) == nt2 || throw(ArgumentError("layout ω2 length mismatch"))
    length(lay.γ2) == nt2 || throw(ArgumentError("layout γ2 length mismatch"))
    return nothing
end

function AdvDiffModelMono(
    cap::AssembledCapacity{N,T},
    D,
    uω,
    uγ;
    source=((args...) -> zero(T)),
    bc::BorderConditions=BorderConditions(),
    bc_interface_diff::Union{Nothing,PenguinBCs.Robin}=nothing,
    layout::UnknownLayout=layout_mono(cap.ntotal),
    coeff_mode::Symbol=:harmonic,
    scheme::AdvectionScheme=Centered(),
) where {N,T}
    bc_diff, bc_adv = _split_borderconditions(bc, N)
    pflags_diff = periodic_flags(bc_diff, N)
    pflags_adv = periodic_flags(bc_adv, N)

    opsD = DiffusionOps(cap; periodic=pflags_diff)
    diff = PenguinDiffusion.DiffusionModelMono(
        cap,
        opsD,
        D;
        source=source,
        bc_border=bc_diff,
        bc_interface=bc_interface_diff,
        layout=layout,
        coeff_mode=coeff_mode,
    )

    return AdvDiffModelMono{N,T,typeof(diff),typeof(uω),typeof(uγ),typeof(bc),typeof(scheme)}(
        diff,
        uω,
        uγ,
        bc,
        scheme,
        pflags_adv,
    )
end

function AdvDiffModelDiph(
    cap1::AssembledCapacity{N,T},
    D1,
    u1ω,
    u1γ,
    cap2::AssembledCapacity{N,T},
    D2,
    u2ω,
    u2γ;
    source=((args...) -> (zero(T), zero(T))),
    bc::BorderConditions=BorderConditions(),
    ic::Union{Nothing,InterfaceConditions}=nothing,
    bc_interface::Union{Nothing,InterfaceConditions}=nothing,
    layout::UnknownLayout=layout_diph(cap1.ntotal),
    coeff_mode::Symbol=:harmonic,
    scheme::AdvectionScheme=Centered(),
) where {N,T}
    bc_diff, bc_adv = _split_borderconditions(bc, N)
    pflags_diff = periodic_flags(bc_diff, N)
    pflags_adv = periodic_flags(bc_adv, N)

    ops1 = DiffusionOps(cap1; periodic=pflags_diff)
    ops2 = DiffusionOps(cap2; periodic=pflags_diff)
    s1, s2 = _split_diph_source(source, T)

    diff = PenguinDiffusion.DiffusionModelDiph(
        cap1, ops1, D1, s1,
        cap2, ops2, D2, s2;
        bc_border=bc_diff,
        ic=ic,
        bc_interface=bc_interface,
        layout=layout,
        coeff_mode=coeff_mode,
    )

    return AdvDiffModelDiph{
        N,T,typeof(diff),typeof(u1ω),typeof(u1γ),typeof(u2ω),typeof(u2γ),typeof(bc),typeof(scheme)
    }(
        diff,
        u1ω,
        u1γ,
        u2ω,
        u2γ,
        bc,
        scheme,
        pflags_adv,
    )
end

function AdvDiffModelDiph(
    cap::AssembledCapacity{N,T},
    D1,
    D2,
    u1ω,
    u1γ,
    u2ω,
    u2γ;
    source=((args...) -> (zero(T), zero(T))),
    bc::BorderConditions=BorderConditions(),
    ic::Union{Nothing,InterfaceConditions}=nothing,
    bc_interface::Union{Nothing,InterfaceConditions}=nothing,
    layout::UnknownLayout=layout_diph(cap.ntotal),
    coeff_mode::Symbol=:harmonic,
    scheme::AdvectionScheme=Centered(),
) where {N,T}
    return AdvDiffModelDiph(
        cap,
        D1,
        u1ω,
        u1γ,
        cap,
        D2,
        u2ω,
        u2γ;
        source=source,
        bc=bc,
        ic=ic,
        bc_interface=bc_interface,
        layout=layout,
        coeff_mode=coeff_mode,
        scheme=scheme,
    )
end

function MovingAdvDiffModelMono(
    diff::PenguinDiffusion.MovingDiffusionModelMono{N,T},
    uω,
    uγ;
    wγ=ntuple(_ -> zero(T), N),
    bc::BorderConditions=diff.bc_border,
    scheme::AdvectionScheme=Centered(),
) where {N,T}
    _, bc_adv = _split_borderconditions(bc, N)
    return MovingAdvDiffModelMono{N,T,typeof(diff),typeof(uω),typeof(uγ),typeof(wγ),typeof(bc),typeof(scheme)}(
        diff,
        uω,
        uγ,
        wγ,
        bc,
        scheme,
        periodic_flags(bc_adv, N),
    )
end

function MovingAdvDiffModelMono(
    grid::PenguinDiffusion.CartesianGrid{N,T},
    body,
    D,
    uω,
    uγ;
    wγ=ntuple(_ -> zero(T), N),
    source=((args...) -> zero(T)),
    bc::BorderConditions=BorderConditions(),
    bc_interface_diff::Union{Nothing,PenguinBCs.Robin}=nothing,
    layout::UnknownLayout=layout_mono(prod(grid.n)),
    coeff_mode::Symbol=:harmonic,
    scheme::AdvectionScheme=Centered(),
    geom_method::Symbol=:vofijul,
) where {N,T}
    bc_diff, _ = _split_borderconditions(bc, N)
    diff = PenguinDiffusion.MovingDiffusionModelMono(
        grid,
        body,
        D;
        source=source,
        bc_border=bc_diff,
        bc_interface=bc_interface_diff,
        layout=layout,
        coeff_mode=coeff_mode,
        geom_method=geom_method,
    )
    return MovingAdvDiffModelMono(diff, uω, uγ; wγ=wγ, bc=bc, scheme=scheme)
end

function MovingAdvDiffModelDiph(
    diff::PenguinDiffusion.MovingDiffusionModelDiph{N,T},
    u1ω,
    u1γ,
    u2ω,
    u2γ;
    wγ=ntuple(_ -> zero(T), N),
    bc::BorderConditions=diff.bc_border,
    scheme::AdvectionScheme=Centered(),
) where {N,T}
    _, bc_adv = _split_borderconditions(bc, N)
    return MovingAdvDiffModelDiph{
        N,T,typeof(diff),typeof(u1ω),typeof(u1γ),typeof(u2ω),typeof(u2γ),typeof(wγ),typeof(bc),typeof(scheme)
    }(
        diff,
        u1ω,
        u1γ,
        u2ω,
        u2γ,
        wγ,
        bc,
        scheme,
        periodic_flags(bc_adv, N),
    )
end

function MovingAdvDiffModelDiph(
    grid::PenguinDiffusion.CartesianGrid{N,T},
    body1,
    D1,
    u1ω,
    u1γ,
    D2,
    u2ω,
    u2γ;
    wγ=ntuple(_ -> zero(T), N),
    source=((args...) -> (zero(T), zero(T))),
    body2=nothing,
    bc::BorderConditions=BorderConditions(),
    ic::Union{Nothing,InterfaceConditions}=nothing,
    bc_interface::Union{Nothing,InterfaceConditions}=nothing,
    layout::UnknownLayout=layout_diph(prod(grid.n)),
    coeff_mode::Symbol=:harmonic,
    scheme::AdvectionScheme=Centered(),
    geom_method::Symbol=:vofijul,
) where {N,T}
    bc_diff, _ = _split_borderconditions(bc, N)
    diff = PenguinDiffusion.MovingDiffusionModelDiph(
        grid,
        body1,
        D1,
        D2;
        source=source,
        body2=body2,
        bc_border=bc_diff,
        ic=ic,
        bc_interface=bc_interface,
        layout=layout,
        coeff_mode=coeff_mode,
        geom_method=geom_method,
    )
    return MovingAdvDiffModelDiph(diff, u1ω, u1γ, u2ω, u2γ; wγ=wγ, bc=bc, scheme=scheme)
end

function _velocity_values_diph(model::AdvDiffModelDiph{N,T}, t::T) where {N,T}
    cap1 = model.diff.cap1
    cap2 = model.diff.cap2
    u1ωv = _velocity_tuple_values(cap1, model.u1ω, cap1.C_ω, t)
    u1γv = _velocity_tuple_values(cap1, model.u1γ, cap1.C_γ, t)
    u2ωv = _velocity_tuple_values(cap2, model.u2ω, cap2.C_ω, t)
    u2γv = _velocity_tuple_values(cap2, model.u2γ, cap2.C_γ, t)
    return u1ωv, u1γv, u2ωv, u2γv
end

function _ops_for_time(model::AdvDiffModelMono{N,T}, t::T) where {N,T}
    cap = model.diff.cap
    uωv, uγv = _velocity_values(cap, model.uω, model.uγ, t)
    ops = advection_ops(cap, uωv, uγv; periodic=model.periodic, scheme=model.scheme)
    return ops, uωv, uγv
end

function _ops_for_time(model::AdvDiffModelDiph{N,T}, t::T) where {N,T}
    u1ωv, u1γv, u2ωv, u2γv = _velocity_values_diph(model, t)
    opsA1 = advection_ops(model.diff.cap1, u1ωv, u1γv; periodic=model.periodic, scheme=model.scheme)
    opsA2 = advection_ops(model.diff.cap2, u2ωv, u2γv; periodic=model.periodic, scheme=model.scheme)
    return opsA1, opsA2, u1ωv, u1γv, u2ωv, u2γv
end

function update_advection_ops!(model::AdvDiffModelMono{N,T}; t::T=zero(T)) where {N,T}
    _ops_for_time(model, t)
    return model
end

function update_advection_ops!(model::AdvDiffModelDiph{N,T}; t::T=zero(T)) where {N,T}
    _ops_for_time(model, t)
    return model
end

function rebuild!(model::AdvDiffModelMono{N,T}, moments; bc=zero(T), t::T=zero(T)) where {N,T}
    CartesianOperators.rebuild!(model.diff.cap, moments; bc=bc)

    bc_diff, bc_adv = _split_borderconditions(model.bc, N)
    pflags_diff = periodic_flags(bc_diff, N)
    opsD = DiffusionOps(model.diff.cap; periodic=pflags_diff)
    model.diff = PenguinDiffusion.DiffusionModelMono(
        model.diff.cap,
        opsD,
        model.diff.D;
        source=model.diff.source,
        bc_border=bc_diff,
        bc_interface=model.diff.bc_interface,
        layout=model.diff.layout,
        coeff_mode=model.diff.coeff_mode,
    )

    model.periodic = periodic_flags(bc_adv, N)
    update_advection_ops!(model; t=t)
    return model
end

function rebuild!(model::AdvDiffModelDiph{N,T}, moments; bc=zero(T), t::T=zero(T)) where {N,T}
    cap1 = model.diff.cap1
    cap2 = model.diff.cap2
    if moments isa Tuple && length(moments) == 2
        CartesianOperators.rebuild!(cap1, moments[1]; bc=bc)
        cap1 === cap2 || CartesianOperators.rebuild!(cap2, moments[2]; bc=bc)
    else
        CartesianOperators.rebuild!(cap1, moments; bc=bc)
        cap1 === cap2 || CartesianOperators.rebuild!(cap2, moments; bc=bc)
    end

    bc_diff, bc_adv = _split_borderconditions(model.bc, N)
    pflags_diff = periodic_flags(bc_diff, N)
    ops1 = DiffusionOps(cap1; periodic=pflags_diff)
    ops2 = DiffusionOps(cap2; periodic=pflags_diff)
    model.diff = PenguinDiffusion.DiffusionModelDiph(
        cap1, ops1, model.diff.D1, model.diff.source1,
        cap2, ops2, model.diff.D2, model.diff.source2;
        bc_border=bc_diff,
        ic=model.diff.ic,
        layout=model.diff.layout,
        coeff_mode=model.diff.coeff_mode,
    )

    model.periodic = periodic_flags(bc_adv, N)
    update_advection_ops!(model; t=t)
    return model
end

function assemble_steady_mono!(sys::LinearSystem{T}, model::AdvDiffModelMono{N,T}, t::T) where {N,T}
    PenguinDiffusion.assemble_steady_mono!(sys, model.diff, t)

    cap = model.diff.cap
    lay = model.diff.layout.offsets
    _validate_mono_layout(cap, lay)

    opsA, uωv, _ = _ops_for_time(model, t)
    conv_bulk = reduce(+, opsA.C)
    conv_iface = convert(T, 0.5) * reduce(+, opsA.K)

    _insert_block!(sys.A, lay.ω, lay.ω, conv_bulk + conv_iface)
    _insert_block!(sys.A, lay.ω, lay.γ, conv_iface)

    _, bc_adv = _split_borderconditions(model.bc, N)
    PenguinTransport.apply_box_bc_transport_mono!(
        sys.A,
        sys.b,
        cap,
        uωv,
        bc_adv,
        model.scheme;
        t=t,
        ωrows=lay.ω,
    )

    active_rows = PenguinTransport._mono_row_activity(cap, lay)
    sys.A, sys.b = PenguinTransport._apply_row_identity_constraints!(sys.A, sys.b, active_rows)
    sys.cache = nothing
    return sys
end

function assemble_steady_diph!(sys::LinearSystem{T}, model::AdvDiffModelDiph{N,T}, t::T) where {N,T}
    PenguinDiffusion.assemble_steady_diph!(sys, model.diff, t)

    cap1 = model.diff.cap1
    cap2 = model.diff.cap2
    lay = model.diff.layout.offsets
    _validate_diph_layout(cap1, cap2, lay)

    opsA1, opsA2, u1ωv, _, u2ωv, _ = _ops_for_time(model, t)

    conv_bulk1 = reduce(+, opsA1.C)
    conv_iface1 = convert(T, 0.5) * reduce(+, opsA1.K)
    conv_bulk2 = reduce(+, opsA2.C)
    conv_iface2 = convert(T, 0.5) * reduce(+, opsA2.K)

    _insert_block!(sys.A, lay.ω1, lay.ω1, conv_bulk1 + conv_iface1)
    _insert_block!(sys.A, lay.ω1, lay.γ1, conv_iface1)
    _insert_block!(sys.A, lay.ω2, lay.ω2, conv_bulk2 + conv_iface2)
    _insert_block!(sys.A, lay.ω2, lay.γ2, conv_iface2)

    _, bc_adv = _split_borderconditions(model.bc, N)
    PenguinTransport.apply_box_bc_transport_mono!(
        sys.A,
        sys.b,
        cap1,
        u1ωv,
        bc_adv,
        model.scheme;
        t=t,
        ωrows=lay.ω1,
    )
    PenguinTransport.apply_box_bc_transport_mono!(
        sys.A,
        sys.b,
        cap2,
        u2ωv,
        bc_adv,
        model.scheme;
        t=t,
        ωrows=lay.ω2,
    )

    active_rows = PenguinDiffusion._diph_row_activity(cap1, cap2, lay)
    sys.A, sys.b = PenguinTransport._apply_row_identity_constraints!(sys.A, sys.b, active_rows)
    sys.cache = nothing
    return sys
end

function _init_unsteady_state_mono(model::AdvDiffModelMono{N,T}, u0) where {N,T}
    lay = model.diff.layout.offsets
    nt = model.diff.cap.ntotal
    nsys = maximum((last(lay.ω), last(lay.γ)))
    u = zeros(T, nsys)
    if length(u0) == nsys
        u .= Vector{T}(u0)
    elseif length(u0) == nt
        u[lay.ω] .= Vector{T}(u0)
    else
        throw(DimensionMismatch("u0 length must be $nt (ω block) or $nsys (full system)"))
    end
    return u
end

function _init_unsteady_state_diph(model::AdvDiffModelDiph{N,T}, u0) where {N,T}
    lay = model.diff.layout.offsets
    nt = model.diff.cap1.ntotal
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))
    u = zeros(T, nsys)
    if length(u0) == nsys
        u .= Vector{T}(u0)
    elseif length(u0) == 2 * nt
        u0v = Vector{T}(u0)
        u[lay.ω1] .= u0v[1:nt]
        u[lay.ω2] .= u0v[(nt + 1):(2 * nt)]
    else
        throw(DimensionMismatch("u0 length must be $(2 * nt) (ω1+ω2) or $nsys (full system)"))
    end
    return u
end

function _init_unsteady_state_mono_moving(model::MovingAdvDiffModelMono{N,T}, u0) where {N,T}
    lay = model.diff.layout.offsets
    nt = prod(model.diff.grid.n)
    nsys = maximum((last(lay.ω), last(lay.γ)))
    u = zeros(T, nsys)
    if length(u0) == nsys
        u .= Vector{T}(u0)
    elseif length(u0) == nt
        u[lay.ω] .= Vector{T}(u0)
        u[lay.γ] .= u[lay.ω]
    else
        throw(DimensionMismatch("u0 length must be $nt (ω block) or $nsys (full system)"))
    end
    return u
end

function _init_unsteady_state_diph_moving(model::MovingAdvDiffModelDiph{N,T}, u0) where {N,T}
    lay = model.diff.layout.offsets
    nt = prod(model.diff.grid.n)
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))
    u = zeros(T, nsys)
    if length(u0) == nsys
        u .= Vector{T}(u0)
    elseif length(u0) == 2 * nt
        u0v = Vector{T}(u0)
        u[lay.ω1] .= u0v[1:nt]
        u[lay.ω2] .= u0v[(nt + 1):(2 * nt)]
        u[lay.γ1] .= u[lay.ω1]
        u[lay.γ2] .= u[lay.ω2]
    else
        throw(DimensionMismatch("u0 length must be $(2 * nt) (ω1+ω2) or $nsys (full system)"))
    end
    return u
end

function assemble_unsteady_mono!(
    sys::LinearSystem{T},
    model::AdvDiffModelMono{N,T},
    uⁿ,
    t::T,
    dt::T,
    scheme_or_theta,
) where {N,T}
    dt > zero(T) || throw(ArgumentError("dt must be positive"))
    θ = _theta_from_scheme(T, scheme_or_theta)
    assemble_steady_mono!(sys, model, t + θ * dt)

    cap = model.diff.cap
    lay = model.diff.layout.offsets
    _validate_mono_layout(cap, lay)

    nt = cap.ntotal
    nsys = maximum((last(lay.ω), last(lay.γ)))
    ufull = if length(uⁿ) == nsys
        Vector{T}(uⁿ)
    elseif length(uⁿ) == nt
        v = zeros(T, nsys)
        v[lay.ω] .= Vector{T}(uⁿ)
        v
    else
        v = zeros(T, nsys)
        v[lay.ω] .= Vector{T}(uⁿ[lay.ω])
        v
    end

    if θ != one(T)
        Aω_prev = sys.A[lay.ω, :]
        corr = Aω_prev * ufull
        _scale_rows!(sys.A, lay.ω, θ)
        _insert_vec!(sys.b, lay.ω, (-(one(T) - θ)) .* corr)
    end

    M = cap.buf.V ./ dt
    sys.A = sys.A + sparse(lay.ω, lay.ω, M, nsys, nsys)
    _insert_vec!(sys.b, lay.ω, M .* Vector{T}(ufull[lay.ω]))

    active_rows = PenguinTransport._mono_row_activity(cap, lay)
    sys.A, sys.b = PenguinTransport._apply_row_identity_constraints!(sys.A, sys.b, active_rows)
    sys.cache = nothing
    return sys
end

function assemble_unsteady_diph!(
    sys::LinearSystem{T},
    model::AdvDiffModelDiph{N,T},
    uⁿ,
    t::T,
    dt::T,
    scheme_or_theta,
) where {N,T}
    dt > zero(T) || throw(ArgumentError("dt must be positive"))
    θ = _theta_from_scheme(T, scheme_or_theta)
    assemble_steady_diph!(sys, model, t + θ * dt)

    lay = model.diff.layout.offsets
    nt = model.diff.cap1.ntotal
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))

    ufull = if length(uⁿ) == nsys
        Vector{T}(uⁿ)
    elseif length(uⁿ) == 2 * nt
        v = zeros(T, nsys)
        u0v = Vector{T}(uⁿ)
        v[lay.ω1] .= u0v[1:nt]
        v[lay.ω2] .= u0v[(nt + 1):(2 * nt)]
        v
    else
        v = zeros(T, nsys)
        v[lay.ω1] .= Vector{T}(uⁿ[lay.ω1])
        v[lay.ω2] .= Vector{T}(uⁿ[lay.ω2])
        v
    end

    if θ != one(T)
        Aω1_prev = sys.A[lay.ω1, :]
        Aω2_prev = sys.A[lay.ω2, :]
        corr1 = Aω1_prev * ufull
        corr2 = Aω2_prev * ufull
        _scale_rows!(sys.A, lay.ω1, θ)
        _scale_rows!(sys.A, lay.ω2, θ)
        _insert_vec!(sys.b, lay.ω1, (-(one(T) - θ)) .* corr1)
        _insert_vec!(sys.b, lay.ω2, (-(one(T) - θ)) .* corr2)
    end

    M1 = model.diff.cap1.buf.V ./ dt
    M2 = model.diff.cap2.buf.V ./ dt
    rows = vcat(collect(lay.ω1), collect(lay.ω2))
    vals = vcat(M1, M2)
    sys.A = sys.A + sparse(rows, rows, vals, nsys, nsys)
    _insert_vec!(sys.b, lay.ω1, M1 .* Vector{T}(ufull[lay.ω1]))
    _insert_vec!(sys.b, lay.ω2, M2 .* Vector{T}(ufull[lay.ω2]))

    active_rows = PenguinDiffusion._diph_row_activity(model.diff.cap1, model.diff.cap2, lay)
    sys.A, sys.b = PenguinTransport._apply_row_identity_constraints!(sys.A, sys.b, active_rows)
    sys.cache = nothing
    return sys
end

function assemble_unsteady_mono_moving!(
    sys::LinearSystem{T},
    model::MovingAdvDiffModelMono{N,T},
    uⁿ,
    t::T,
    dt::T,
    scheme_or_theta,
) where {N,T}
    dt > zero(T) || throw(ArgumentError("dt must be positive"))
    θ = _theta_from_scheme(T, scheme_or_theta)

    PenguinDiffusion.assemble_unsteady_mono_moving!(sys, model.diff, uⁿ, t, dt; scheme=θ)

    cap = something(model.diff.cap_slab)
    lay = model.diff.layout.offsets
    _validate_mono_layout(cap, lay)

    τ = t + θ * dt
    uωv, uγv = _velocity_values(cap, model.uω, model.uγ, τ)
    wγv = _velocity_tuple_values(cap, model.wγ, cap.C_γ, τ)
    uγrel = PenguinTransport._relative_interface_velocity(uγv, wγv)

    opsA = PenguinTransport._advection_ops_moving(
        cap,
        uωv,
        uγrel;
        periodic=model.periodic,
        scheme=model.scheme,
    )

    conv_bulk = reduce(+, opsA.C)
    conv_iface = convert(T, 0.5) * reduce(+, opsA.K)

    _insert_block!(sys.A, lay.ω, lay.ω, conv_bulk + conv_iface)
    _insert_block!(sys.A, lay.ω, lay.γ, conv_iface)

    _, bc_adv = _split_borderconditions(model.bc, N)
    PenguinTransport.apply_box_bc_transport_mono!(
        sys.A,
        sys.b,
        cap,
        uωv,
        bc_adv,
        model.scheme;
        t=τ,
        ωrows=lay.ω,
    )

    active_rows = PenguinTransport._mono_row_activity(cap, lay)
    sys.A, sys.b = PenguinTransport._apply_row_identity_constraints!(sys.A, sys.b, active_rows)
    sys.cache = nothing
    return sys
end

function assemble_unsteady_diph_moving!(
    sys::LinearSystem{T},
    model::MovingAdvDiffModelDiph{N,T},
    uⁿ,
    t::T,
    dt::T,
    scheme_or_theta,
) where {N,T}
    dt > zero(T) || throw(ArgumentError("dt must be positive"))
    θ = _theta_from_scheme(T, scheme_or_theta)

    PenguinDiffusion.assemble_unsteady_diph_moving!(sys, model.diff, uⁿ, t, dt; scheme=θ)

    cap1 = something(model.diff.cap1_slab)
    cap2 = something(model.diff.cap2_slab)
    lay = model.diff.layout.offsets
    _validate_diph_layout(cap1, cap2, lay)

    τ = t + θ * dt
    u1ωv = _velocity_tuple_values(cap1, model.u1ω, cap1.C_ω, τ)
    u1γv = _velocity_tuple_values(cap1, model.u1γ, cap1.C_γ, τ)
    u2ωv = _velocity_tuple_values(cap2, model.u2ω, cap2.C_ω, τ)
    u2γv = _velocity_tuple_values(cap2, model.u2γ, cap2.C_γ, τ)
    wγv = _velocity_tuple_values(cap1, model.wγ, cap1.C_γ, τ)

    u1γrel = PenguinTransport._relative_interface_velocity(u1γv, wγv)
    u2γrel = PenguinTransport._relative_interface_velocity(u2γv, wγv)
    opsA1 = PenguinTransport._advection_ops_moving(cap1, u1ωv, u1γrel; periodic=model.periodic, scheme=model.scheme)
    opsA2 = PenguinTransport._advection_ops_moving(cap2, u2ωv, u2γrel; periodic=model.periodic, scheme=model.scheme)

    conv_bulk1 = reduce(+, opsA1.C)
    conv_iface1 = convert(T, 0.5) * reduce(+, opsA1.K)
    conv_bulk2 = reduce(+, opsA2.C)
    conv_iface2 = convert(T, 0.5) * reduce(+, opsA2.K)

    _insert_block!(sys.A, lay.ω1, lay.ω1, conv_bulk1 + conv_iface1)
    _insert_block!(sys.A, lay.ω1, lay.γ1, conv_iface1)
    _insert_block!(sys.A, lay.ω2, lay.ω2, conv_bulk2 + conv_iface2)
    _insert_block!(sys.A, lay.ω2, lay.γ2, conv_iface2)

    _, bc_adv = _split_borderconditions(model.bc, N)
    PenguinTransport.apply_box_bc_transport_mono!(
        sys.A,
        sys.b,
        cap1,
        u1ωv,
        bc_adv,
        model.scheme;
        t=τ,
        ωrows=lay.ω1,
    )
    PenguinTransport.apply_box_bc_transport_mono!(
        sys.A,
        sys.b,
        cap2,
        u2ωv,
        bc_adv,
        model.scheme;
        t=τ,
        ωrows=lay.ω2,
    )

    active_rows = PenguinDiffusion._diph_row_activity(cap1, cap2, lay)
    sys.A, sys.b = PenguinTransport._apply_row_identity_constraints!(sys.A, sys.b, active_rows)
    sys.cache = nothing
    return sys
end

function PenguinSolverCore.assemble!(sys::LinearSystem{T}, model::AdvDiffModelMono{N,T}, t, dt) where {N,T}
    assemble_unsteady_mono!(sys, model, sys.x, convert(T, t), convert(T, dt), one(T))
end

function PenguinSolverCore.assemble!(sys::LinearSystem{T}, model::AdvDiffModelDiph{N,T}, t, dt) where {N,T}
    assemble_unsteady_diph!(sys, model, sys.x, convert(T, t), convert(T, dt), one(T))
end

function PenguinSolverCore.assemble!(sys::LinearSystem{T}, model::MovingAdvDiffModelMono{N,T}, t, dt) where {N,T}
    assemble_unsteady_mono_moving!(sys, model, sys.x, convert(T, t), convert(T, dt), one(T))
end

function PenguinSolverCore.assemble!(sys::LinearSystem{T}, model::MovingAdvDiffModelDiph{N,T}, t, dt) where {N,T}
    assemble_unsteady_diph_moving!(sys, model, sys.x, convert(T, t), convert(T, dt), one(T))
end

function solve_steady!(model::AdvDiffModelMono{N,T}; t::T=zero(T), method::Symbol=:direct, kwargs...) where {N,T}
    lay = model.diff.layout.offsets
    n = maximum((last(lay.ω), last(lay.γ)))
    sys = LinearSystem(spzeros(T, n, n), zeros(T, n))
    assemble_steady_mono!(sys, model, t)
    solve!(sys; method=method, reuse_factorization=false, kwargs...)
    return (system=sys, solution=copy(sys.x))
end

function solve_steady!(model::AdvDiffModelDiph{N,T}; t::T=zero(T), method::Symbol=:direct, kwargs...) where {N,T}
    lay = model.diff.layout.offsets
    n = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))
    sys = LinearSystem(spzeros(T, n, n), zeros(T, n))
    assemble_steady_diph!(sys, model, t)
    solve!(sys; method=method, reuse_factorization=false, kwargs...)
    return (system=sys, solution=copy(sys.x))
end

function solve_unsteady_moving!(
    model::MovingAdvDiffModelMono{N,T},
    u0,
    tspan::Tuple{T,T};
    dt::T,
    scheme=:BE,
    method::Symbol=:direct,
    save_history::Bool=true,
    kwargs...,
) where {N,T}
    t0, tend = tspan
    tend >= t0 || throw(ArgumentError("tspan must satisfy tend >= t0"))
    dt > zero(T) || throw(ArgumentError("dt must be positive"))
    θ = _theta_from_scheme(T, scheme)

    u = _init_unsteady_state_mono_moving(model, u0)
    lay = model.diff.layout.offsets
    nsys = maximum((last(lay.ω), last(lay.γ)))

    times = T[t0]
    states = Vector{Vector{T}}()
    save_history && push!(states, copy(u))

    sys = LinearSystem(spzeros(T, nsys, nsys), zeros(T, nsys); x=copy(u))
    t = t0
    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))

    while t < tend - tol
        dt_step = min(dt, tend - t)
        assemble_unsteady_mono_moving!(sys, model, u, t, dt_step, θ)
        solve!(sys; method=method, reuse_factorization=false, kwargs...)
        u .= sys.x
        t += dt_step
        push!(times, t)
        save_history && push!(states, copy(u))
    end

    if !save_history
        states = [copy(u)]
        times = T[t]
    end
    return (times=times, states=states, system=sys, reused_constant_operator=false)
end

function solve_unsteady_moving!(
    model::MovingAdvDiffModelDiph{N,T},
    u0,
    tspan::Tuple{T,T};
    dt::T,
    scheme=:BE,
    method::Symbol=:direct,
    save_history::Bool=true,
    kwargs...,
) where {N,T}
    t0, tend = tspan
    tend >= t0 || throw(ArgumentError("tspan must satisfy tend >= t0"))
    dt > zero(T) || throw(ArgumentError("dt must be positive"))
    θ = _theta_from_scheme(T, scheme)

    u = _init_unsteady_state_diph_moving(model, u0)
    lay = model.diff.layout.offsets
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))

    times = T[t0]
    states = Vector{Vector{T}}()
    save_history && push!(states, copy(u))

    sys = LinearSystem(spzeros(T, nsys, nsys), zeros(T, nsys); x=copy(u))
    t = t0
    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))

    while t < tend - tol
        dt_step = min(dt, tend - t)
        assemble_unsteady_diph_moving!(sys, model, u, t, dt_step, θ)
        solve!(sys; method=method, reuse_factorization=false, kwargs...)
        u .= sys.x
        t += dt_step
        push!(times, t)
        save_history && push!(states, copy(u))
    end

    if !save_history
        states = [copy(u)]
        times = T[t]
    end
    return (times=times, states=states, system=sys, reused_constant_operator=false)
end

function _matrix_time_dependent(model::AdvDiffModelMono{N,T}) where {N,T}
    xω = model.diff.cap.C_ω[1]
    xγ = model.diff.cap.C_γ[1]
    return PenguinDiffusion._mono_matrix_time_dependent(model.diff) ||
           _vel_time_dependent(model.uω, xω) ||
           _vel_time_dependent(model.uγ, xγ)
end

function _rhs_time_dependent(model::AdvDiffModelMono{N,T}) where {N,T}
    xω = model.diff.cap.C_ω[1]
    xγ = model.diff.cap.C_γ[1]
    return PenguinDiffusion._mono_rhs_time_dependent(model.diff) ||
           _adv_border_time_dependent(model.bc, xω) ||
           _source_time_dependent(model.diff.source, xω) ||
           _vel_time_dependent(model.uω, xω) ||
           _vel_time_dependent(model.uγ, xγ)
end

function _matrix_time_dependent(model::AdvDiffModelDiph{N,T}) where {N,T}
    x1ω = model.diff.cap1.C_ω[1]
    x1γ = model.diff.cap1.C_γ[1]
    x2ω = model.diff.cap2.C_ω[1]
    x2γ = model.diff.cap2.C_γ[1]
    return PenguinDiffusion._diph_matrix_time_dependent(model.diff) ||
           _vel_time_dependent(model.u1ω, x1ω) ||
           _vel_time_dependent(model.u1γ, x1γ) ||
           _vel_time_dependent(model.u2ω, x2ω) ||
           _vel_time_dependent(model.u2γ, x2γ)
end

function _rhs_time_dependent(model::AdvDiffModelDiph{N,T}) where {N,T}
    x1ω = model.diff.cap1.C_ω[1]
    x2ω = model.diff.cap2.C_ω[1]
    x1γ = model.diff.cap1.C_γ[1]
    x2γ = model.diff.cap2.C_γ[1]
    return PenguinDiffusion._diph_rhs_time_dependent(model.diff) ||
           _adv_border_time_dependent(model.bc, x1ω) ||
           _adv_border_time_dependent(model.bc, x2ω) ||
           _source_time_dependent(model.diff.source1, x1ω) ||
           _source_time_dependent(model.diff.source2, x2ω) ||
           _vel_time_dependent(model.u1ω, x1ω) ||
           _vel_time_dependent(model.u1γ, x1γ) ||
           _vel_time_dependent(model.u2ω, x2ω) ||
           _vel_time_dependent(model.u2γ, x2γ)
end

function _prepare_constant_unsteady_mono(model::AdvDiffModelMono{N,T}, t0::T, dt::T, θ::T) where {N,T}
    lay = model.diff.layout.offsets
    nsys = maximum((last(lay.ω), last(lay.γ)))
    sys0 = LinearSystem(spzeros(T, nsys, nsys), zeros(T, nsys))
    assemble_steady_mono!(sys0, model, t0 + θ * dt)
    Asteady = sys0.A
    bsteady = copy(sys0.b)
    Aconst = copy(Asteady)
    θ != one(T) && _scale_rows!(Aconst, lay.ω, θ)
    M = model.diff.cap.buf.V ./ dt
    Aconst = Aconst + sparse(lay.ω, lay.ω, M, nsys, nsys)
    Aω_prev = Asteady[lay.ω, :]
    return Aconst, bsteady, Aω_prev, M
end

function _prepare_constant_unsteady_diph(model::AdvDiffModelDiph{N,T}, t0::T, dt::T, θ::T) where {N,T}
    lay = model.diff.layout.offsets
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))
    sys0 = LinearSystem(spzeros(T, nsys, nsys), zeros(T, nsys))
    assemble_steady_diph!(sys0, model, t0 + θ * dt)
    Asteady = sys0.A
    bsteady = copy(sys0.b)
    Aconst = copy(Asteady)
    θ != one(T) && _scale_rows!(Aconst, lay.ω1, θ)
    θ != one(T) && _scale_rows!(Aconst, lay.ω2, θ)
    M1 = model.diff.cap1.buf.V ./ dt
    M2 = model.diff.cap2.buf.V ./ dt
    rows = vcat(collect(lay.ω1), collect(lay.ω2))
    vals = vcat(M1, M2)
    Aconst = Aconst + sparse(rows, rows, vals, nsys, nsys)
    Aω1_prev = Asteady[lay.ω1, :]
    Aω2_prev = Asteady[lay.ω2, :]
    return Aconst, bsteady, Aω1_prev, Aω2_prev, M1, M2
end

function _set_constant_rhs_mono!(
    b::Vector{T},
    bsteady::Vector{T},
    Aω_prev::SparseMatrixCSC{T,Int},
    M::Vector{T},
    lay,
    u::Vector{T},
    θ::T,
) where {T}
    copyto!(b, bsteady)
    if θ != one(T)
        corr = Aω_prev * u
        _insert_vec!(b, lay.ω, (-(one(T) - θ)) .* corr)
    end
    _insert_vec!(b, lay.ω, M .* Vector{T}(u[lay.ω]))
    return b
end

function _set_constant_rhs_diph!(
    b::Vector{T},
    bsteady::Vector{T},
    Aω1_prev::SparseMatrixCSC{T,Int},
    Aω2_prev::SparseMatrixCSC{T,Int},
    M1::Vector{T},
    M2::Vector{T},
    lay,
    u::Vector{T},
    θ::T,
) where {T}
    copyto!(b, bsteady)
    if θ != one(T)
        corr1 = Aω1_prev * u
        corr2 = Aω2_prev * u
        _insert_vec!(b, lay.ω1, (-(one(T) - θ)) .* corr1)
        _insert_vec!(b, lay.ω2, (-(one(T) - θ)) .* corr2)
    end
    _insert_vec!(b, lay.ω1, M1 .* Vector{T}(u[lay.ω1]))
    _insert_vec!(b, lay.ω2, M2 .* Vector{T}(u[lay.ω2]))
    return b
end

function solve_unsteady!(
    model::AdvDiffModelMono{N,T},
    u0,
    tspan::Tuple{T,T};
    dt::T,
    scheme=:BE,
    method::Symbol=:direct,
    save_history::Bool=true,
    kwargs...,
) where {N,T}
    t0, tend = tspan
    tend >= t0 || throw(ArgumentError("tspan must satisfy tend >= t0"))
    dt > zero(T) || throw(ArgumentError("dt must be positive"))
    θ = _theta_from_scheme(T, scheme)

    u = _init_unsteady_state_mono(model, u0)
    lay = model.diff.layout.offsets
    nsys = maximum((last(lay.ω), last(lay.γ)))

    matrix_dep = _matrix_time_dependent(model)
    rhs_dep = _rhs_time_dependent(model)
    constant_operator = !matrix_dep && !rhs_dep

    times = T[t0]
    states = Vector{Vector{T}}()
    save_history && push!(states, copy(u))

    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))
    t = t0

    if constant_operator
        Aconst, bsteady, Aω_prev, M = _prepare_constant_unsteady_mono(model, t0, dt, θ)
        sys = LinearSystem(Aconst, copy(bsteady); x=copy(u))
        while t < tend - tol
            dt_step = min(dt, tend - t)
            if abs(dt_step - dt) <= tol
                _set_constant_rhs_mono!(sys.b, bsteady, Aω_prev, M, lay, u, θ)
                solve!(sys; method=method, reuse_factorization=true, kwargs...)
            else
                assemble_unsteady_mono!(sys, model, u, t, dt_step, θ)
                solve!(sys; method=method, reuse_factorization=false, kwargs...)
            end
            u .= sys.x
            t += dt_step
            push!(times, t)
            save_history && push!(states, copy(u))
        end
        if !save_history
            states = [copy(u)]
            times = T[t]
        end
        return (times=times, states=states, system=sys, reused_constant_operator=true)
    end

    sys = LinearSystem(spzeros(T, nsys, nsys), zeros(T, nsys); x=copy(u))
    while t < tend - tol
        dt_step = min(dt, tend - t)
        assemble_unsteady_mono!(sys, model, u, t, dt_step, θ)
        solve!(sys; method=method, reuse_factorization=false, kwargs...)
        u .= sys.x
        t += dt_step
        push!(times, t)
        save_history && push!(states, copy(u))
    end

    if !save_history
        states = [copy(u)]
        times = T[t]
    end
    return (times=times, states=states, system=sys, reused_constant_operator=false)
end

function solve_unsteady!(
    model::AdvDiffModelDiph{N,T},
    u0,
    tspan::Tuple{T,T};
    dt::T,
    scheme=:BE,
    method::Symbol=:direct,
    save_history::Bool=true,
    kwargs...,
) where {N,T}
    t0, tend = tspan
    tend >= t0 || throw(ArgumentError("tspan must satisfy tend >= t0"))
    dt > zero(T) || throw(ArgumentError("dt must be positive"))
    θ = _theta_from_scheme(T, scheme)

    u = _init_unsteady_state_diph(model, u0)
    lay = model.diff.layout.offsets
    nsys = maximum((last(lay.ω1), last(lay.γ1), last(lay.ω2), last(lay.γ2)))

    matrix_dep = _matrix_time_dependent(model)
    rhs_dep = _rhs_time_dependent(model)
    constant_operator = !matrix_dep && !rhs_dep

    times = T[t0]
    states = Vector{Vector{T}}()
    save_history && push!(states, copy(u))

    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))
    t = t0

    if constant_operator
        Aconst, bsteady, Aω1_prev, Aω2_prev, M1, M2 = _prepare_constant_unsteady_diph(model, t0, dt, θ)
        sys = LinearSystem(Aconst, copy(bsteady); x=copy(u))
        while t < tend - tol
            dt_step = min(dt, tend - t)
            if abs(dt_step - dt) <= tol
                _set_constant_rhs_diph!(sys.b, bsteady, Aω1_prev, Aω2_prev, M1, M2, lay, u, θ)
                solve!(sys; method=method, reuse_factorization=true, kwargs...)
            else
                assemble_unsteady_diph!(sys, model, u, t, dt_step, θ)
                solve!(sys; method=method, reuse_factorization=false, kwargs...)
            end
            u .= sys.x
            t += dt_step
            push!(times, t)
            save_history && push!(states, copy(u))
        end
        if !save_history
            states = [copy(u)]
            times = T[t]
        end
        return (times=times, states=states, system=sys, reused_constant_operator=true)
    end

    sys = LinearSystem(spzeros(T, nsys, nsys), zeros(T, nsys); x=copy(u))
    while t < tend - tol
        dt_step = min(dt, tend - t)
        assemble_unsteady_diph!(sys, model, u, t, dt_step, θ)
        solve!(sys; method=method, reuse_factorization=false, kwargs...)
        u .= sys.x
        t += dt_step
        push!(times, t)
        save_history && push!(states, copy(u))
    end

    if !save_history
        states = [copy(u)]
        times = T[t]
    end
    return (times=times, states=states, system=sys, reused_constant_operator=false)
end

include("coupling.jl")

end # module
