module PenguinTransportDiffusion

using LinearAlgebra
using SparseArrays
using StaticArrays

using CartesianOperators
using PenguinBCs
using PenguinSolverCore

using PenguinDiffusion
using PenguinTransport

export AdvDiffModelMono
export assemble_steady_mono!, assemble_unsteady_mono!
export solve_steady!, solve_unsteady!
export rebuild!, update_advection_ops!

mutable struct AdvDiffModelMono{N,T,DM,VTω,VTγ,BCT,SCH}
    diff::DM
    uω::VTω
    uγ::VTγ
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
        elseif side_bc isa Dirichlet || side_bc isa Neumann
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
        return convert(T, scheme)
    end
    throw(ArgumentError("scheme must be a Symbol (:BE/:CN) or a numeric theta"))
end

function _validate_mono_layout(cap::AssembledCapacity, lay)
    nt = cap.ntotal
    length(lay.ω) == nt || throw(ArgumentError("layout ω length ($(length(lay.ω))) must match cap.ntotal ($nt)"))
    length(lay.γ) == nt || throw(ArgumentError("layout γ length ($(length(lay.γ))) must match cap.ntotal ($nt)"))
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

function _ops_for_time(model::AdvDiffModelMono{N,T}, t::T) where {N,T}
    cap = model.diff.cap
    uωv, uγv = _velocity_values(cap, model.uω, model.uγ, t)
    ops = advection_ops(cap, uωv, uγv; periodic=model.periodic, scheme=model.scheme)
    return ops, uωv, uγv
end

function update_advection_ops!(model::AdvDiffModelMono{N,T}; t::T=zero(T)) where {N,T}
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
        layout=model.diff.layout,
    )

    active_rows = PenguinTransport._mono_row_activity(cap, lay)
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

function PenguinSolverCore.assemble!(sys::LinearSystem{T}, model::AdvDiffModelMono{N,T}, t, dt) where {N,T}
    assemble_unsteady_mono!(sys, model, sys.x, convert(T, t), convert(T, dt), one(T))
end

function solve_steady!(model::AdvDiffModelMono{N,T}; t::T=zero(T), method::Symbol=:direct, kwargs...) where {N,T}
    lay = model.diff.layout.offsets
    n = maximum((last(lay.ω), last(lay.γ)))
    sys = LinearSystem(spzeros(T, n, n), zeros(T, n))
    assemble_steady_mono!(sys, model, t)
    solve!(sys; method=method, reuse_factorization=false, kwargs...)
    return (system=sys, solution=copy(sys.x))
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

    times = T[t0]
    states = Vector{Vector{T}}()
    save_history && push!(states, copy(u))

    sys = LinearSystem(spzeros(T, nsys, nsys), zeros(T, nsys); x=copy(u))
    tol = sqrt(eps(T)) * max(one(T), abs(t0), abs(tend))
    t = t0

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

end # module
