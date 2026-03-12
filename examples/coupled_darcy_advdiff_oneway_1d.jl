using CartesianGeometry: geometric_moments, nan
using CartesianOperators
using PenguinBCs
using PenguinDarcy
using PenguinSolverCore
using PenguinTransportDiffusion

full_moments(grid) = geometric_moments((args...) -> -1.0, grid, Float64, nan; method=:vofijul)

n = 65
grid = (range(0.0, 1.0; length=n),)
cap = assembled_capacity(full_moments(grid); bc=0.0)
nt = cap.ntotal

bc_darcy = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
ops = DiffusionOps(cap; periodic=periodic_flags(bc_darcy, 1))
darcy_base = DarcyModelMono(cap, ops, 1.0; source=(x, t) -> 0.0, bc_border=bc_darcy)
darcy_block = CoupledBlock(:darcy, DarcyCoupledModelMono(darcy_base); init=nothing, cache=Dict{Symbol,Any}())

bc_transport = BorderConditions(; left=Dirichlet(1.0), right=Dirichlet(0.0))
uzero = (zeros(nt),)
adv_base = AdvDiffModelMono(cap, 2e-3, uzero, uzero; source=0.0, bc=bc_transport, scheme=Centered())

lay = adv_base.diff.layout.offsets
nsys = maximum((last(lay.ω), last(lay.γ)))
c0 = zeros(Float64, nsys)
@inbounds for i in 1:nt
    c0[lay.ω[i]] = cap.C_ω[i][1]
end

transport_block = CoupledBlock(:transport, AdvDiffCoupledModelMono(adv_base);
                               init=(concentration=c0,),
                               cache=Dict{Symbol,Any}())

problem = CoupledProblem(
    [darcy_block, transport_block],
    OneWayCoupling([:darcy, :transport]);
    maps=[CouplingMap(:darcy, :transport, :velocity)],
)

dt = 0.01
step_coupled!(problem, 0.0, dt; method=:direct)
c_coupled = transport_block.state.concentration

# Reference run with manually injected Darcy velocity.
sys_d = solve_steady!(darcy_base; method=:direct)
vel = recover_velocity(darcy_base, sys_d.x)
adv_ref = AdvDiffModelMono(cap, 2e-3, (vel.x,), (vel.x,); source=0.0, bc=bc_transport, scheme=Centered())
sol_ref = solve_unsteady!(adv_ref, c0, (0.0, dt); dt=dt, method=:direct, save_history=false)
c_ref = sol_ref.system.x

println("One-way Darcy -> AdvDiff 1D")
println("  dt = $dt")
println("  max|c_coupled - c_ref| = ", maximum(abs.(c_coupled .- c_ref)))
