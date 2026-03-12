using PenguinSolverCore

@testset "AdvDiffCoupledModelMono SolverCore adapter" begin
    grid = (range(0.0, 1.0; length=49),)
    cap = assembled_capacity(full_moments(grid); bc=0.0)
    nt = cap.ntotal

    bc = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
    D = 2e-3
    uzeros = (zeros(nt),)

    base = AdvDiffModelMono(cap, D, uzeros, uzeros; source=0.0, bc=bc, scheme=Centered())
    coupled = AdvDiffCoupledModelMono(base)

    lay = base.diff.layout.offsets
    nsys = maximum((last(lay.ω), last(lay.γ)))
    c0 = zeros(nsys)

    block = CoupledBlock(:transport, coupled; init=(concentration=c0,), cache=Dict{Symbol,Any}())

    vel1 = (x=fill(0.4, nt),)
    set_coupling_field!(block, Val(:velocity), vel1)
    @test get_coupling_field(block, Val(:velocity)) == vel1

    c_before = copy(get_coupling_field(block, Val(:concentration)))
    advance_unsteady!(block, 0.0, 0.02; method=:direct, scheme=:BE)

    c_after_1 = get_coupling_field(block, Val(:concentration))
    @test length(c_after_1) == nsys
    @test block.state.refresh_count == 1
    @test norm(c_after_1 - c_before) > 0.0

    vel2 = (x=fill(1.1, nt),)
    set_coupling_field!(block, Val(:velocity), vel2)
    c_prev = copy(c_after_1)
    advance_unsteady!(block, 0.02, 0.02; method=:direct, scheme=:BE)

    c_after_2 = get_coupling_field(block, Val(:concentration))
    @test block.state.refresh_count == 2
    @test norm(c_after_2 - c_prev) > 0.0
    @test haskey(block.cache, :advection_refreshed_at)

    s = PenguinSolverCore.block_summary(block)
    @test occursin("AdvDiffCoupledModelMono", s)
end
