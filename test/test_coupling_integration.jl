import PenguinDarcy
using PenguinSolverCore

function _build_coupling_case(; n=49, λ0=1.0, α=0.0, D=2e-3)
    grid = (range(0.0, 1.0; length=n),)
    cap = assembled_capacity(full_moments(grid); bc=0.0)

    bc_darcy = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
    ops = DiffusionOps(cap; periodic=periodic_flags(bc_darcy, 1))
    darcy_base = PenguinDarcy.DarcyModelMono(cap, ops, λ0; source=(x, t) -> 0.0, bc_border=bc_darcy)

    bc_transport = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
    uzero = (zeros(cap.ntotal),)
    adv_base = AdvDiffModelMono(cap, D, uzero, uzero; source=0.0, bc=bc_transport, scheme=Centered())

    lay = adv_base.diff.layout.offsets
    nsys = maximum((last(lay.ω), last(lay.γ)))
    c0 = zeros(Float64, nsys)
    @inbounds for i in 1:cap.ntotal
        c0[lay.ω[i]] = cap.C_ω[i][1]
    end

    mobility_cb = function (base_model, concentration, t)
        laym = base_model.layout.offsets
        cω = length(concentration) >= last(laym.ω) ? concentration[laym.ω] : concentration
        csum = 0.0
        nfinite = 0
        @inbounds for value in cω
            isfinite(value) || continue
            csum += value
            nfinite += 1
        end
        cavg = nfinite == 0 ? 0.0 : csum / nfinite
        return max(1e-8, λ0 * (1 + α * cavg))
    end

    darcy_coupled = PenguinDarcy.DarcyCoupledModelMono(darcy_base; mobility_update! = mobility_cb)
    adv_coupled = AdvDiffCoupledModelMono(adv_base)

    darcy_block = CoupledBlock(:darcy, darcy_coupled; init=nothing, cache=Dict{Symbol,Any}())
    transport_block = CoupledBlock(:transport, adv_coupled; init=(concentration=c0,), cache=Dict{Symbol,Any}())

    return (cap=cap,
            darcy_base=darcy_base,
            adv_base=adv_base,
            c0=c0,
            darcy_block=darcy_block,
            transport_block=transport_block,
            λ0=λ0,
            D=D)
end

@testset "Coupling integration: one-way equivalence" begin
    dt = 0.01
    case = _build_coupling_case(; α=0.0)

    sysd = PenguinDarcy.solve_steady!(case.darcy_base; method=:direct)
    vel = PenguinDarcy.recover_velocity(case.darcy_base, sysd.x)
    adv_ref = AdvDiffModelMono(case.adv_base.diff.cap, case.D, (vel.x,), (vel.x,);
                               source=case.adv_base.diff.source,
                               bc=case.adv_base.bc,
                               bc_interface_diff=case.adv_base.diff.bc_interface,
                               layout=case.adv_base.diff.layout,
                               coeff_mode=case.adv_base.diff.coeff_mode,
                               scheme=case.adv_base.scheme)

    ref = solve_unsteady!(adv_ref, case.c0, (0.0, dt); dt=dt, method=:direct, save_history=false)
    cref = ref.system.x

    one_way = CoupledProblem(
        [case.darcy_block, case.transport_block],
        OneWayCoupling([:darcy, :transport]);
        maps=[CouplingMap(:darcy, :transport, :velocity)],
    )

    step_coupled!(one_way, 0.0, dt; method=:direct)
    ccoupled = case.transport_block.state.concentration

    lay = case.adv_base.diff.layout.offsets
    idx = active_physical_indices(case.cap)
    @test isapprox(ccoupled[lay.ω][idx], cref[lay.ω][idx]; atol=1e-10, rtol=1e-10)
end

@testset "Coupling integration: two-way residual decrease" begin
    dt = 0.01
    case = _build_coupling_case(; α=0.6)

    mode = TwoWayCoupling(
        [:darcy, :transport];
        maxiter=30,
        atol=1e-10,
        rtol=1e-8,
        relaxation=Dict(:darcy => 0.7, :transport => 1.0),
        norm_fields=[:velocity, :concentration],
        sweep=:GaussSeidel,
        verbose=false,
    )

    problem = CoupledProblem(
        [case.darcy_block, case.transport_block],
        mode;
        maps=[
            CouplingMap(:darcy, :transport, :velocity),
            CouplingMap(:transport, :darcy, :concentration),
        ],
    )

    _, history = step_coupled!(problem, 0.0, dt; method=:direct, return_history=true)

    @test length(history.iterations) >= 2
    @test history.residuals[end] < history.residuals[1]

    threshold = mode.atol + mode.rtol * max(history.residuals[1], eps())
    @test history.residuals[end] <= threshold
end

@testset "Coupling integration: zero-feedback degeneracy" begin
    dt = 0.01

    case_one = _build_coupling_case(; α=0.0)
    one_way = CoupledProblem(
        [case_one.darcy_block, case_one.transport_block],
        OneWayCoupling([:darcy, :transport]);
        maps=[CouplingMap(:darcy, :transport, :velocity)],
    )
    step_coupled!(one_way, 0.0, dt; method=:direct)
    c_one = copy(case_one.transport_block.state.concentration)

    case_two = _build_coupling_case(; α=0.0)
    two_mode = TwoWayCoupling(
        [:darcy, :transport];
        maxiter=40,
        atol=1e-8,
        rtol=1e-6,
        relaxation=Dict(:darcy => 0.7, :transport => 1.0),
        norm_fields=[:velocity],
        sweep=:GaussSeidel,
        verbose=false,
    )

    two_way = CoupledProblem(
        [case_two.darcy_block, case_two.transport_block],
        two_mode;
        maps=[
            CouplingMap(:darcy, :transport, :velocity),
            CouplingMap(:transport, :darcy, :concentration),
        ],
    )

    _, history = step_coupled!(two_way, 0.0, dt; method=:direct, return_history=true, throw_on_nonconvergence=false)
    c_two = case_two.transport_block.state.concentration

    @test history.residuals[end] <= history.residuals[1]

    lay = case_two.adv_base.diff.layout.offsets
    idx = active_physical_indices(case_two.cap)
    @test isapprox(c_two[lay.ω][idx], c_one[lay.ω][idx]; atol=1e-6, rtol=1e-6)
end
