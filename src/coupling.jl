mutable struct AdvDiffCoupledModelMono
    base_model::Any
    velocity_projection::Any
    metadata::Dict{Symbol,Any}
end

function AdvDiffCoupledModelMono(base_model::AdvDiffModelMono;
                                 velocity_projection = nothing,
                                 metadata::AbstractDict{Symbol,<:Any}=Dict{Symbol,Any}())
    return AdvDiffCoupledModelMono(base_model, velocity_projection, Dict{Symbol,Any}(metadata))
end

mutable struct AdvDiffCoupledStateMono
    concentration
    incoming_velocity
    last_time
    refresh_count::Int
end

_advdiff_scalar_type(model::AdvDiffModelMono{N,T}) where {N,T} = T

function _advdiff_system_size(model::AdvDiffModelMono)
    lay = model.diff.layout.offsets
    return maximum((last(lay.ω), last(lay.γ)))
end

function _init_get(init, key::Symbol, default)
    if init === nothing
        return default
    elseif init isa NamedTuple
        return hasproperty(init, key) ? getproperty(init, key) : default
    elseif init isa AbstractDict
        return get(init, key, default)
    end
    return default
end

function _without_keys(kwargs::Base.Pairs, drop::Tuple{Vararg{Symbol}})
    kept = Symbol[]
    vals = Any[]
    for key in keys(kwargs)
        key in drop && continue
        push!(kept, key)
        push!(vals, kwargs[key])
    end
    return NamedTuple{Tuple(kept)}(Tuple(vals))
end

_clone_velocity_component(v) = v isa AbstractArray ? copy(v) : deepcopy(v)

function _sanitize_scalar_field(data)
    if data isa AbstractArray
        clean = copy(data)
        @inbounds for i in eachindex(clean)
            if !isfinite(clean[i])
                clean[i] = zero(eltype(clean))
            end
        end
        return clean
    end
    return deepcopy(data)
end

function _normalize_velocity_components(data, N::Int)
    comps = if data isa Tuple
        length(data) == N || throw(ArgumentError("velocity tuple must have $N components."))
        data
    elseif data isa AbstractVector
        length(data) == N || throw(ArgumentError("velocity container must have $N components."))
        Tuple(data)
    else
        throw(ArgumentError("velocity data must be tuple/vector with $N components."))
    end

    return ntuple(d -> _clone_velocity_component(comps[d]), N)
end

_component_symbol(d::Int) = d == 1 ? :x : d == 2 ? :y : :z

function _components_from_namedtuple(data::NamedTuple, N::Int)
    return ntuple(d -> begin
        name = _component_symbol(d)
        hasproperty(data, name) || throw(ArgumentError("named velocity must provide component `:$name`."))
        _clone_velocity_component(getproperty(data, name))
    end, N)
end

function _default_velocity_projection(model::AdvDiffModelMono{N}, data) where {N}
    if data isa NamedTuple && hasproperty(data, :ω) && hasproperty(data, :γ)
        uω = _normalize_velocity_components(getproperty(data, :ω), N)
        uγ = _normalize_velocity_components(getproperty(data, :γ), N)
        return uω, uγ
    elseif data isa NamedTuple
        comps = _components_from_namedtuple(data, N)
        return comps, deepcopy(comps)
    end

    comps = _normalize_velocity_components(data, N)
    return comps, deepcopy(comps)
end

function _project_velocity(coupled::AdvDiffCoupledModelMono, data)
    callback = coupled.velocity_projection
    callback === nothing && return _default_velocity_projection(coupled.base_model, data)

    if applicable(callback, coupled.base_model, data)
        return callback(coupled.base_model, data)
    elseif applicable(callback, data)
        return callback(data)
    end

    throw(ArgumentError(
        "velocity_projection must accept (base_model, data) or (data)."))
end

function _rebuild_advdiff_with_velocity(model::AdvDiffModelMono, uω, uγ)
    return AdvDiffModelMono(
        model.diff.cap,
        model.diff.D,
        uω,
        uγ;
        source=model.diff.source,
        bc=model.bc,
        bc_interface_diff=model.diff.bc_interface,
        layout=model.diff.layout,
        coeff_mode=model.diff.coeff_mode,
        scheme=model.scheme,
    )
end

function _set_model_velocity!(coupled::AdvDiffCoupledModelMono, uω, uγ)
    base = coupled.base_model
    if uω isa typeof(base.uω) && uγ isa typeof(base.uγ)
        base.uω = uω
        base.uγ = uγ
        return base
    end

    coupled.base_model = _rebuild_advdiff_with_velocity(base, uω, uγ)
    return coupled.base_model
end

function PenguinSolverCore.initialize_state(model::AdvDiffCoupledModelMono, init)
    base = model.base_model
    T = _advdiff_scalar_type(base)
    nsys = _advdiff_system_size(base)

    concentration = _init_get(init, :concentration, zeros(T, nsys))
    velocity = _init_get(init, :velocity, (ω=deepcopy(base.uω), γ=deepcopy(base.uγ)))
    last_time = _init_get(init, :time, nothing)

    return AdvDiffCoupledStateMono(concentration, velocity, last_time, 0)
end

function PenguinSolverCore.copy_state(state::AdvDiffCoupledStateMono)
    return AdvDiffCoupledStateMono(
        deepcopy(state.concentration),
        deepcopy(state.incoming_velocity),
        deepcopy(state.last_time),
        state.refresh_count,
    )
end

function PenguinSolverCore.advance_steady!(block::CoupledBlock{M,S,C}; kwargs...) where {M<:AdvDiffCoupledModelMono,S<:AdvDiffCoupledStateMono,C}
    coupled = block.model
    state = block.state
    T = _advdiff_scalar_type(coupled.base_model)

    tval = haskey(kwargs, :t) ? convert(T, kwargs[:t]) : zero(T)
    if state.incoming_velocity !== nothing
        uω, uγ = _project_velocity(coupled, state.incoming_velocity)
        _set_model_velocity!(coupled, uω, uγ)
    end

    update_advection_ops!(coupled.base_model; t=tval)
    state.refresh_count += 1

    solve_kwargs = _without_keys(kwargs, (:t,))
    sol = solve_steady!(coupled.base_model; t=tval, solve_kwargs...)

    state.concentration = _sanitize_scalar_field(sol.solution)
    state.last_time = tval

    if block.cache isa AbstractDict
        block.cache[:advdiff_system] = sol.system
        block.cache[:advection_refreshed_at] = tval
    end

    return block
end

function PenguinSolverCore.advance_unsteady!(block::CoupledBlock{M,S,C}, t, dt; kwargs...) where {M<:AdvDiffCoupledModelMono,S<:AdvDiffCoupledStateMono,C}
    coupled = block.model
    state = block.state
    T = _advdiff_scalar_type(coupled.base_model)

    t0 = convert(T, t)
    dt_step = convert(T, dt)

    if state.incoming_velocity !== nothing
        uω, uγ = _project_velocity(coupled, state.incoming_velocity)
        _set_model_velocity!(coupled, uω, uγ)
    end

    update_advection_ops!(coupled.base_model; t=t0)
    state.refresh_count += 1

    step_signature = (t0, dt_step)
    u0 = if block.cache isa AbstractDict
        if get(block.cache, :step_signature, nothing) != step_signature
            base_state = state.concentration === nothing ? zeros(T, _advdiff_system_size(coupled.base_model)) : copy(state.concentration)
            block.cache[:step_signature] = step_signature
            block.cache[:step_base_concentration] = base_state
        end
        block.cache[:step_base_concentration]
    else
        state.concentration === nothing ? zeros(T, _advdiff_system_size(coupled.base_model)) : state.concentration
    end

    solve_kwargs = _without_keys(kwargs, (:t, :dt, :save_history))
    sol = solve_unsteady!(
        coupled.base_model,
        u0,
        (t0, t0 + dt_step);
        dt=dt_step,
        save_history=false,
        solve_kwargs...,
    )

    state.concentration = _sanitize_scalar_field(sol.system.x)
    state.last_time = t0 + dt_step

    if block.cache isa AbstractDict
        block.cache[:advdiff_system] = sol.system
        block.cache[:advection_refreshed_at] = t0
    end

    return block
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:concentration}) where {M<:AdvDiffCoupledModelMono,S<:AdvDiffCoupledStateMono,C}
    return _sanitize_scalar_field(block.state.concentration)
end

function PenguinSolverCore.get_coupling_field(block::CoupledBlock{M,S,C}, ::Val{:velocity}) where {M<:AdvDiffCoupledModelMono,S<:AdvDiffCoupledStateMono,C}
    return block.state.incoming_velocity
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{M,S,C}, ::Val{:velocity}, data) where {M<:AdvDiffCoupledModelMono,S<:AdvDiffCoupledStateMono,C}
    block.state.incoming_velocity = deepcopy(data)
    return block
end

function PenguinSolverCore.set_coupling_field!(block::CoupledBlock{M,S,C}, ::Val{:concentration}, data) where {M<:AdvDiffCoupledModelMono,S<:AdvDiffCoupledStateMono,C}
    block.state.concentration = _sanitize_scalar_field(data)
    return block
end

function PenguinSolverCore.block_summary(block::CoupledBlock{M,S,C}) where {M<:AdvDiffCoupledModelMono,S<:AdvDiffCoupledStateMono,C}
    return "AdvDiffCoupledModelMono(name=$(block.name), refresh_count=$(block.state.refresh_count))"
end
