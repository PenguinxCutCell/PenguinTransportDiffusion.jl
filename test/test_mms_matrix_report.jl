# ============================================================
# Helpers
# ============================================================

function _active_indices_mms(cap)
    LI = LinearIndices(cap.nnodes)
    idx = Int[]
    N = length(cap.nnodes)
    for I in CartesianIndices(cap.nnodes)
        i = LI[I]
        if all(d -> I[d] < cap.nnodes[d], 1:N)
            v = cap.buf.V[i]
            if isfinite(v) && v > 0.0
                push!(idx, i)
            end
        end
    end
    return idx
end

function _l2_mono_mms(cap, uω, exact, t)
    idx = _active_indices_mms(cap)
    return sqrt(sum(cap.buf.V[i] * (uω[i] - exact(cap.C_ω[i][1], t))^2 for i in idx))
end

function _l2_diph_mms(cap1, u1, cap2, u2, exact1, exact2, t)
    idx1 = _active_indices_mms(cap1)
    idx2 = _active_indices_mms(cap2)
    e1 = sum(cap1.buf.V[i] * (u1[i] - exact1(cap1.C_ω[i][1], t))^2 for i in idx1)
    e2 = sum(cap2.buf.V[i] * (u2[i] - exact2(cap2.C_ω[i][1], t))^2 for i in idx2)
    return sqrt(e1 + e2)
end

function _rel_l2_mono(cap, uω, exact, t)
    idx = _active_indices_mms(cap)
    num = sum(cap.buf.V[i] * (uω[i] - exact(cap.C_ω[i][1], t))^2 for i in idx)
    den = sum(cap.buf.V[i] * exact(cap.C_ω[i][1], t)^2 for i in idx)
    return sqrt(num / max(den, eps()))
end

function _rel_l2_diph(cap1, u1, cap2, u2, exact1, exact2, t)
    idx1 = _active_indices_mms(cap1)
    idx2 = _active_indices_mms(cap2)
    num1 = sum(cap1.buf.V[i] * (u1[i] - exact1(cap1.C_ω[i][1], t))^2 for i in idx1)
    num2 = sum(cap2.buf.V[i] * (u2[i] - exact2(cap2.C_ω[i][1], t))^2 for i in idx2)
    den1 = sum(cap1.buf.V[i] * exact1(cap1.C_ω[i][1], t)^2 for i in idx1)
    den2 = sum(cap2.buf.V[i] * exact2(cap2.C_ω[i][1], t)^2 for i in idx2)
    return sqrt((num1 + num2) / max(den1 + den2, eps()))
end

# ============================================================
# Baseline Gaussian exact solution (same as before)
# ============================================================

const _D_mms = 1e-2
const _u_mms = 0.35
const _tfinal_mms = 0.01
const _xγ0_mms = 0.6
const _sγ_mms = _u_mms
const _xg_mms = 0.35
const _σ_mms = 0.06
const _σ2_mms = _σ_mms^2

function _exact_mms(x, t)
    xt = _xg_mms + _u_mms * t
    s2 = _σ2_mms + 2 * _D_mms * t
    amp = sqrt(_σ2_mms / s2)
    return amp * exp(-((x - xt)^2) / (2s2))
end

_source_mms(x, t) = 0.0

_body_static_mms(x) = x - _xγ0_mms
_body_moving_mms(x, t) = x - (_xγ0_mms + _sγ_mms * t)

const _bc_periodic_mms = BorderConditions(; left=Periodic(), right=Periodic())
const _bc_embedded_mms = BorderConditions(
    ; left=Dirichlet((x, t) -> _exact_mms(x, t)),
      right=Dirichlet((x, t) -> _exact_mms(x, t)),
)
const _ic_mms = InterfaceConditions(
    ; scalar=ScalarJump(1.0, 1.0, 0.0),
      flux=FluxJump(1.0, 1.0, 0.0),
)
const _bc_iface_mono_mms = Robin(1.0, 0.0, (x, t) -> _exact_mms(x, t))

function _run_case_mms(; kind::Symbol, scheme::Symbol, embedded::Bool, moving::Bool)
    ns = (33, 65, 129)
    hs = Float64[]
    errs = Float64[]

    for n in ns
        x = range(0.0, 1.0; length=n)
        grid1 = (x,)
        h = step(x)
        dt = scheme === :BE ? 0.12 * h^2 : 0.2 * h^2

        if kind === :mono
            if !moving
                body = embedded ? _body_static_mms : (x -> -1.0)
                moms = geometric_moments(body, grid1, Float64, nan; method=:vofijul)
                cap = assembled_capacity(moms; bc=0.0)
                model = AdvDiffModelMono(
                    cap,
                    _D_mms,
                    (_u_mms,),
                    (_u_mms,);
                    source=_source_mms,
                    bc=(embedded ? _bc_embedded_mms : _bc_periodic_mms),
                    bc_interface_diff=(embedded ? _bc_iface_mono_mms : nothing),
                    scheme=Centered(),
                )
                u0 = [_exact_mms(cap.C_ω[i][1], 0.0) for i in 1:cap.ntotal]
                sol = solve_unsteady!(model, u0, (0.0, _tfinal_mms); dt=dt, scheme=scheme, method=:direct, save_history=false)
                lay = model.diff.layout.offsets
                uω = sol.system.x[lay.ω]
                err = _l2_mono_mms(cap, uω, _exact_mms, _tfinal_mms)
            else
                grid = PenguinDiffusion.CartesianGrid((0.0,), (1.0,), (n,))
                body = embedded ? _body_moving_mms : ((x, t) -> -1.0)
                wγ = embedded ? (_sγ_mms,) : (0.0,)
                model = MovingAdvDiffModelMono(
                    grid,
                    body,
                    _D_mms,
                    (_u_mms,),
                    (_u_mms,);
                    wγ=wγ,
                    source=_source_mms,
                    bc=(embedded ? _bc_embedded_mms : _bc_periodic_mms),
                    bc_interface_diff=(embedded ? _bc_iface_mono_mms : nothing),
                    scheme=Centered(),
                )
                u0 = [_exact_mms(xi, 0.0) for xi in x]
                sol = solve_unsteady_moving!(model, u0, (0.0, _tfinal_mms); dt=dt, scheme=scheme, method=:direct, save_history=false)
                cap = something(model.diff.cap_slab)
                lay = model.diff.layout.offsets
                uω = sol.system.x[lay.ω]
                err = _l2_mono_mms(cap, uω, _exact_mms, _tfinal_mms)
            end
        else
            if !moving
                if embedded
                    moms1 = geometric_moments(_body_static_mms, grid1, Float64, nan; method=:vofijul)
                    moms2 = geometric_moments(x -> -_body_static_mms(x), grid1, Float64, nan; method=:vofijul)
                    cap1 = assembled_capacity(moms1; bc=0.0)
                    cap2 = assembled_capacity(moms2; bc=0.0)
                else
                    moms = geometric_moments(x -> -1.0, grid1, Float64, nan; method=:vofijul)
                    cap1 = assembled_capacity(moms; bc=0.0)
                    cap2 = assembled_capacity(moms; bc=0.0)
                end
                model = AdvDiffModelDiph(
                    cap1,
                    _D_mms,
                    (_u_mms,),
                    (_u_mms,),
                    cap2,
                    _D_mms,
                    (_u_mms,),
                    (_u_mms,);
                    source=(_source_mms, _source_mms),
                    bc=(embedded ? _bc_embedded_mms : _bc_periodic_mms),
                    ic=_ic_mms,
                    scheme=Centered(),
                )
                u01 = [_exact_mms(cap1.C_ω[i][1], 0.0) for i in 1:cap1.ntotal]
                u02 = [_exact_mms(cap2.C_ω[i][1], 0.0) for i in 1:cap2.ntotal]
                u0 = vcat(u01, u02)
                sol = solve_unsteady!(model, u0, (0.0, _tfinal_mms); dt=dt, scheme=scheme, method=:direct, save_history=false)
                lay = model.diff.layout.offsets
                u1 = sol.system.x[lay.ω1]
                u2 = sol.system.x[lay.ω2]
                err = _l2_diph_mms(cap1, u1, cap2, u2, _exact_mms, _exact_mms, _tfinal_mms)
            else
                grid = PenguinDiffusion.CartesianGrid((0.0,), (1.0,), (n,))
                body1 = embedded ? _body_moving_mms : ((x, t) -> -1.0)
                body2 = embedded ? nothing : ((x, t) -> -1.0)
                wγ = embedded ? (_sγ_mms,) : (0.0,)
                model = MovingAdvDiffModelDiph(
                    grid,
                    body1,
                    _D_mms,
                    (_u_mms,),
                    (_u_mms,),
                    _D_mms,
                    (_u_mms,),
                    (_u_mms,);
                    wγ=wγ,
                    source=(_source_mms, _source_mms),
                    body2=body2,
                    bc=(embedded ? _bc_embedded_mms : _bc_periodic_mms),
                    ic=_ic_mms,
                    scheme=Centered(),
                )
                u0cell = [_exact_mms(xi, 0.0) for xi in x]
                u0 = vcat(u0cell, u0cell)
                sol = solve_unsteady_moving!(model, u0, (0.0, _tfinal_mms); dt=dt, scheme=scheme, method=:direct, save_history=false)
                cap1 = something(model.diff.cap1_slab)
                cap2 = something(model.diff.cap2_slab)
                lay = model.diff.layout.offsets
                u1 = sol.system.x[lay.ω1]
                u2 = sol.system.x[lay.ω2]
                err = _l2_diph_mms(cap1, u1, cap2, u2, _exact_mms, _exact_mms, _tfinal_mms)
            end
        end

        push!(hs, h)
        push!(errs, err)
    end

    p1 = log(errs[1] / errs[2]) / log(hs[1] / hs[2])
    p2 = log(errs[2] / errs[3]) / log(hs[2] / hs[3])
    return (; hs, errs, p1, p2)
end

# ============================================================
# Sharper test 1:
# direct static-recovery comparison
# ============================================================

function _direct_static_recovery_mono(; n=129, scheme=:BE, embedded=true)
    x = range(0.0, 1.0; length=n)
    grid1 = (x,)
    grid = PenguinDiffusion.CartesianGrid((0.0,), (1.0,), (n,))
    h = step(x)
    dt = scheme === :BE ? 0.12 * h^2 : 0.2 * h^2

    body_static = embedded ? _body_static_mms : (x -> -1.0)
    body_moving = embedded ? ((x, t) -> _body_static_mms(x)) : ((x, t) -> -1.0)

    moms = geometric_moments(body_static, grid1, Float64, nan; method=:vofijul)
    cap = assembled_capacity(moms; bc=0.0)

    model_static = AdvDiffModelMono(
        cap,
        _D_mms,
        (_u_mms,),
        (_u_mms,);
        source=_source_mms,
        bc=(embedded ? _bc_embedded_mms : _bc_periodic_mms),
        bc_interface_diff=(embedded ? _bc_iface_mono_mms : nothing),
        scheme=Centered(),
    )

    model_moving = MovingAdvDiffModelMono(
        grid,
        body_moving,
        _D_mms,
        (_u_mms,),
        (_u_mms,);
        wγ=(0.0,),
        source=_source_mms,
        bc=(embedded ? _bc_embedded_mms : _bc_periodic_mms),
        bc_interface_diff=(embedded ? _bc_iface_mono_mms : nothing),
        scheme=Centered(),
    )

    u0_static = [_exact_mms(cap.C_ω[i][1], 0.0) for i in 1:cap.ntotal]
    u0_moving = [_exact_mms(xi, 0.0) for xi in x]

    sol_static = solve_unsteady!(model_static, u0_static, (0.0, _tfinal_mms); dt=dt, scheme=scheme, method=:direct, save_history=false)
    sol_moving = solve_unsteady_moving!(model_moving, u0_moving, (0.0, _tfinal_mms); dt=dt, scheme=scheme, method=:direct, save_history=false)

    lay_s = model_static.diff.layout.offsets
    lay_m = model_moving.diff.layout.offsets

    us = sol_static.system.x[lay_s.ω]
    um = sol_moving.system.x[lay_m.ω]

    # Compare both against the same exact field on their own active supports
    err_static = _rel_l2_mono(cap, us, _exact_mms, _tfinal_mms)

    capm = something(model_moving.diff.cap_slab)
    err_moving = _rel_l2_mono(capm, um, _exact_mms, _tfinal_mms)

    return (; err_static, err_moving, diff_abs=norm(us - um, Inf))
end

function _direct_static_recovery_diph(; n=129, scheme=:BE, embedded=true)
    x = range(0.0, 1.0; length=n)
    grid1 = (x,)
    grid = PenguinDiffusion.CartesianGrid((0.0,), (1.0,), (n,))
    h = step(x)
    dt = scheme === :BE ? 0.12 * h^2 : 0.2 * h^2

    if embedded
        moms1 = geometric_moments(_body_static_mms, grid1, Float64, nan; method=:vofijul)
        moms2 = geometric_moments(x -> -_body_static_mms(x), grid1, Float64, nan; method=:vofijul)
        cap1 = assembled_capacity(moms1; bc=0.0)
        cap2 = assembled_capacity(moms2; bc=0.0)
    else
        moms = geometric_moments(x -> -1.0, grid1, Float64, nan; method=:vofijul)
        cap1 = assembled_capacity(moms; bc=0.0)
        cap2 = assembled_capacity(moms; bc=0.0)
    end

    body1 = embedded ? ((x, t) -> _body_static_mms(x)) : ((x, t) -> -1.0)
    body2 = embedded ? nothing : ((x, t) -> -1.0)

    model_static = AdvDiffModelDiph(
        cap1,
        _D_mms,
        (_u_mms,),
        (_u_mms,),
        cap2,
        _D_mms,
        (_u_mms,),
        (_u_mms,);
        source=(_source_mms, _source_mms),
        bc=(embedded ? _bc_embedded_mms : _bc_periodic_mms),
        ic=_ic_mms,
        scheme=Centered(),
    )

    model_moving = MovingAdvDiffModelDiph(
        grid,
        body1,
        _D_mms,
        (_u_mms,),
        (_u_mms,),
        _D_mms,
        (_u_mms,),
        (_u_mms,);
        wγ=(0.0,),
        source=(_source_mms, _source_mms),
        body2=body2,
        bc=(embedded ? _bc_embedded_mms : _bc_periodic_mms),
        ic=_ic_mms,
        scheme=Centered(),
    )

    u01 = [_exact_mms(cap1.C_ω[i][1], 0.0) for i in 1:cap1.ntotal]
    u02 = [_exact_mms(cap2.C_ω[i][1], 0.0) for i in 1:cap2.ntotal]
    u0s = vcat(u01, u02)

    u0cell = [_exact_mms(xi, 0.0) for xi in x]
    u0m = vcat(u0cell, u0cell)

    sol_static = solve_unsteady!(model_static, u0s, (0.0, _tfinal_mms); dt=dt, scheme=scheme, method=:direct, save_history=false)
    sol_moving = solve_unsteady_moving!(model_moving, u0m, (0.0, _tfinal_mms); dt=dt, scheme=scheme, method=:direct, save_history=false)

    lay_s = model_static.diff.layout.offsets
    lay_m = model_moving.diff.layout.offsets

    u1s = sol_static.system.x[lay_s.ω1]
    u2s = sol_static.system.x[lay_s.ω2]

    u1m = sol_moving.system.x[lay_m.ω1]
    u2m = sol_moving.system.x[lay_m.ω2]

    err_static = _rel_l2_diph(cap1, u1s, cap2, u2s, _exact_mms, _exact_mms, _tfinal_mms)

    cap1m = something(model_moving.diff.cap1_slab)
    cap2m = something(model_moving.diff.cap2_slab)
    err_moving = _rel_l2_diph(cap1m, u1m, cap2m, u2m, _exact_mms, _exact_mms, _tfinal_mms)

    return (; err_static, err_moving, diff_abs=max(norm(u1s - u1m, Inf), norm(u2s - u2m, Inf)))
end

# ============================================================
# Sharper test 2:
# moving embedded mono with uγ - wγ ≠ 0
# ============================================================

const _u_mono_rel = 0.35
const _w_mono_rel = 0.15
const _D_mono_rel = 1e-2
const _tfinal_rel = 0.01
const _xγ0_rel = 0.58
const _xg_rel = 0.33
const _σ_rel = 0.07
const _σ2_rel = _σ_rel^2

function _exact_mono_rel(x, t)
    xt = _xg_rel + _u_mono_rel * t
    s2 = _σ2_rel + 2 * _D_mono_rel * t
    amp = sqrt(_σ2_rel / s2)
    return amp * exp(-((x - xt)^2) / (2s2))
end

_body_moving_rel(x, t) = x - (_xγ0_rel + _w_mono_rel * t)
_source_mono_rel(x, t) = 0.0
_bc_embedded_rel = BorderConditions(
    ; left=Dirichlet((x, t) -> _exact_mono_rel(x, t)),
      right=Dirichlet((x, t) -> _exact_mono_rel(x, t)),
)
_bc_iface_mono_rel = Robin(1.0, 0.0, (x, t) -> _exact_mono_rel(x, t))

function _run_moving_embedded_mono_relative(; scheme::Symbol, n)
    x = range(0.0, 1.0; length=n)
    h = step(x)
    dt = scheme === :BE ? 0.12 * h^2 : 0.2 * h^2
    grid = PenguinDiffusion.CartesianGrid((0.0,), (1.0,), (n,))

    model = MovingAdvDiffModelMono(
        grid,
        _body_moving_rel,
        _D_mono_rel,
        (_u_mono_rel,),
        (_u_mono_rel,);
        wγ=(_w_mono_rel,),
        source=_source_mono_rel,
        bc=_bc_embedded_rel,
        bc_interface_diff=_bc_iface_mono_rel,
        scheme=Centered(),
    )

    u0 = [_exact_mono_rel(xi, 0.0) for xi in x]
    sol = solve_unsteady_moving!(model, u0, (0.0, _tfinal_rel); dt=dt, scheme=scheme, method=:direct, save_history=false)

    cap = something(model.diff.cap_slab)
    lay = model.diff.layout.offsets
    uω = sol.system.x[lay.ω]
    return _l2_mono_mms(cap, uω, _exact_mono_rel, _tfinal_rel)
end

# ============================================================
# Sharper test 3:
# moving embedded diph with nontrivial jump data
# ============================================================

const _D1_jump = 1e-2
const _D2_jump = 2e-2
const _u1_jump = 0.30
const _u2_jump = 0.30
const _w_jump = 0.12
const _jump_c = 0.4
const _xγ0_jump = 0.56
const _tfinal_jump = 0.01

const _σ_jump = 0.065
const _σ2_jump = _σ_jump^2
const _xg_jump = 0.34

# Same smooth carrier + constant offset in phase 2 to create a nonzero scalar jump.
function _g_jump(x, t)
    xt = _xg_jump + _u1_jump * t
    s2 = _σ2_jump + 2 * _D1_jump * t
    amp = sqrt(_σ2_jump / s2)
    return amp * exp(-((x - xt)^2) / (2s2))
end

_exact1_jump(x, t) = _g_jump(x, t)
_exact2_jump(x, t) = _g_jump(x, t) + _jump_c

# Because exact2 differs by a constant only, ∂x exact1 = ∂x exact2.
# With D1 ≠ D2, the flux jump is nontrivial:
#   D1 ∂n u1 - D2 ∂n u2 = (D1 - D2) ∂n g
# We impose that through a manufactured flux jump function sampled at γ.
function _dgdx_jump(x, t)
    g = _g_jump(x, t)
    s2 = _σ2_jump + 2 * _D1_jump * t
    xt = _xg_jump + _u1_jump * t
    return -(x - xt) / s2 * g
end

function _scalar_jump_jump(x, t)
    return _exact2_jump(x, t) - _exact1_jump(x, t)
end

function _flux_jump_jump(x, t)
    # Sign convention may need adaptation to your InterfaceConditions definition.
    return (_D2_jump - _D1_jump) * _dgdx_jump(x, t)
end

_body_moving_jump(x, t) = x - (_xγ0_jump + _w_jump * t)

# Keep zero bulk sources for this first sharp regression.
# The main purpose is to stress nontrivial interface coupling data.
_source1_jump(x, t) = 0.0
_source2_jump(x, t) = 0.0

const _bc_embedded_jump = BorderConditions(
    ; left=Dirichlet((x, t) -> _exact1_jump(x, t)),
      right=Dirichlet((x, t) -> _exact2_jump(x, t)),
)

const _ic_jump = InterfaceConditions(
    ; scalar=ScalarJump(1.0, 1.0, _scalar_jump_jump),
      flux=FluxJump(1.0, 1.0, _flux_jump_jump),
)

function _run_moving_embedded_diph_jump(; scheme::Symbol, n)
    x = range(0.0, 1.0; length=n)
    h = step(x)
    dt = scheme === :BE ? 0.08 * h^2 : 0.16 * h^2
    grid = PenguinDiffusion.CartesianGrid((0.0,), (1.0,), (n,))

    model = MovingAdvDiffModelDiph(
        grid,
        _body_moving_jump,
        _D1_jump,
        (_u1_jump,),
        (_u1_jump,),
        _D2_jump,
        (_u2_jump,),
        (_u2_jump,);
        wγ=(_w_jump,),
        source=(_source1_jump, _source2_jump),
        bc=_bc_embedded_jump,
        ic=_ic_jump,
        scheme=Centered(),
    )

    # Start from exact phase fields sampled on the box nodes/cells
    u01 = [_exact1_jump(xi, 0.0) for xi in x]
    u02 = [_exact2_jump(xi, 0.0) for xi in x]
    u0 = vcat(u01, u02)

    sol = solve_unsteady_moving!(model, u0, (0.0, _tfinal_jump); dt=dt, scheme=scheme, method=:direct, save_history=false)

    cap1 = something(model.diff.cap1_slab)
    cap2 = something(model.diff.cap2_slab)
    lay = model.diff.layout.offsets
    u1 = sol.system.x[lay.ω1]
    u2 = sol.system.x[lay.ω2]

    return _l2_diph_mms(cap1, u1, cap2, u2, _exact1_jump, _exact2_jump, _tfinal_jump)
end

# ============================================================
# Original matrix report kept as a regression/consistency matrix
# ============================================================

@testset "MMS matrix report" begin
    println("\nMMS convergence matrix")
    println("kind | scheme | embedded | moving | p12 | p23")
    println("-----|--------|----------|--------|-----|-----")
    for kind in (:mono, :diph), scheme in (:CN, :BE), embedded in (false, true), moving in (false, true)
        r = _run_case_mms(; kind=kind, scheme=scheme, embedded=embedded, moving=moving)
        println(
            lpad(String(kind), 4), " | ",
            lpad(String(scheme), 6), " | ",
            lpad(string(embedded), 8), " | ",
            lpad(string(moving), 6), " | ",
            lpad(string(round(r.p1; digits=2)), 4), " | ",
            lpad(string(round(r.p2; digits=2)), 4),
        )
    end
end

# ============================================================
# New sharper tests
# ============================================================

@testset "Direct static recovery comparison" begin
    r1 = _direct_static_recovery_mono(; n=129, scheme=:BE, embedded=true)

    println(r1.err_static, " ", r1.err_moving, " ", r1.diff_abs)

    r2 = _direct_static_recovery_diph(; n=129, scheme=:BE, embedded=true)
    println(r2.err_static, " ", r2.err_moving, " ", r2.diff_abs)
end

@testset "Moving embedded mono with nonzero relative interface speed" begin
    ns = (33, 65, 129)
    errs_be = [_run_moving_embedded_mono_relative(; scheme=:BE, n=n) for n in ns]
    errs_cn = [_run_moving_embedded_mono_relative(; scheme=:CN, n=n) for n in ns]


    pbe1 = log(errs_be[1] / errs_be[2]) / log(2)
    pbe2 = log(errs_be[2] / errs_be[3]) / log(2)
    pcn1 = log(errs_cn[1] / errs_cn[2]) / log(2)
    pcn2 = log(errs_cn[2] / errs_cn[3]) / log(2)

    println(errs_be, " ", pbe1, " ", pbe2)
    println(errs_cn, " ", pcn1, " ", pcn2)


end

@testset "Moving embedded diph with nontrivial jump data" begin
    ns = (33, 65, 129)
    errs_be = [_run_moving_embedded_diph_jump(; scheme=:BE, n=n) for n in ns]
    errs_cn = [_run_moving_embedded_diph_jump(; scheme=:CN, n=n) for n in ns]


    pbe1 = log(errs_be[1] / errs_be[2]) / log(2)
    pbe2 = log(errs_be[2] / errs_be[3]) / log(2)
    pcn1 = log(errs_cn[1] / errs_cn[2]) / log(2)
    pcn2 = log(errs_cn[2] / errs_cn[3]) / log(2)

    println(errs_be, " ", pbe1, " ", pbe2)
    println(errs_cn, " ", pcn1, " ", pcn2)

end