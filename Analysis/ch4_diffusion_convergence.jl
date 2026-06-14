using Pkg
Pkg.activate(@__DIR__)

using GeothermalWells
import GeothermalWells: get_thermal_conductivity,
    get_volumetric_heat_capacity,
    eigen_estimator_get_dmax,
    advection!,
    create_advection_index_lists,
    ADI_and_ADV_callback!,
    AbstractInletModel
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2, DiscreteCallback
using KernelAbstractions: CPU, adapt, @kernel, @index, zeros
# using CUDA: CUDABackend
using LaTeXStrings
using Plots

include(joinpath("..", "paths.jl"))
include("theme.jl")
ensure_plots_dir()


# define new structs so one can dispatch to different thermal_conductivity / rho_c functions with constant diffusion
struct MaterialPropertiesConvergence{RealT<:Real} <: AbstractMaterialProperties{RealT}
    D::RealT
end

struct BoreholeConvergence{RealT<:Real} <: AbstractBorehole{RealT}
    h::RealT
    # dummy variables:
    xc::RealT
    yc::RealT
    r_inner::RealT
    t_inner::RealT
    r_outer::RealT
end

# dispatch to functions for constant diffusion coefficient and unit volumetric heat capacity
@inline get_thermal_conductivity(x, y, z, boreholes::Tuple{Vararg{BoreholeConvergence}}, materials::MaterialPropertiesConvergence) = materials.D
@inline eigen_estimator_get_dmax(materials::MaterialPropertiesConvergence) = materials.D
@inline get_volumetric_heat_capacity(x, y, z, boreholes::Tuple{Vararg{BoreholeConvergence}}, materials::MaterialPropertiesConvergence) = one(typeof(materials.D))

# no advection for convergence test
@inline advection!(temp, dt, t, cache, boreholes::Tuple{Vararg{BoreholeConvergence}}) = nothing
# @inline create_advection_index_lists(backend, gridx, gridy, gridz, boreholes::Tuple{Vararg{BoreholeConvergence}}) = (NaN, NaN, NaN, NaN, NaN, NaN, 1, NaN)


function diffusion_simulation(N, Δt, backend=backend)
    D = 5e-3 # diffusion coefficient
    materials_convergence = MaterialPropertiesConvergence(D)
    borehole_convergence = BoreholeConvergence{typeof(D)}(-1, #= rest are dummy variables to not trip up create_advection_index_lists =# 0, 0, 0, 0, 0)
    boreholes_convergence = (borehole_convergence,)


    # analytical solution for pure diffusion from a point source
    # dont put t = 0
    @kernel function gaussian_diffusion_kernel!(ϕ, gridx, gridy, gridz, t, D)
        k, j, i = @index(Global, NTuple)
        ϕ[k, j, i] = (4π * D * t)^(-3 / 2) * exp(-(gridx[i]^2 + gridy[j]^2 + gridz[k]^2) / (4D * t))
    end

    xmin, xmax = -1, 1
    ymin, ymax = -1, 1
    zmin, zmax = -1, 1

    gridx = range(xmin, xmax; length=N)
    gridy = range(ymin, ymax; length=N)
    gridz = range(zmin, zmax; length=N)

    t0 = 1.0 # cant be 0
    ϕ_convergence_initial = zeros(backend, Float64, N, N, N)
    gaussian_diffusion_kernel!(backend)(ϕ_convergence_initial, gridx, gridy, gridz, t0, D, ndrange=(N, N, N))
    inlet = ConstantInlet(20.0) # dummy inlet model
    cache = create_cache(backend=backend, gridx=gridx, gridy=gridy, gridz=gridz, materials=materials_convergence, boreholes=boreholes_convergence, inlet_model=inlet)

    callback = DiscreteCallback((u, t, integrator) -> true, ADI_and_ADV_callback!, save_positions=(false, false))

    tspan = (0.0, 1.0)

    prob = ODEProblem(rhs_diffusion_z!, ϕ_convergence_initial, tspan, cache)
    println("Starting simulation for N = $N and dt = $Δt")
    @time sol_xyz = solve(prob, ROCK2(max_stages=100, eigen_est=eigen_estimator), save_everystep=false, callback=callback, adaptive=false, dt=Δt)

    t_end = tspan[2] + t0
    # overwrite the inital condition to be the analytical solution at t = t_end
    gaussian_diffusion_kernel!(backend)(ϕ_convergence_initial, gridx, gridy, gridz, t_end, D, ndrange=(N, N, N))
    rms_error = sqrt(sum((sol_xyz[end] .- ϕ_convergence_initial) .^ 2) / length(sol_xyz[end]))
    println("rms error:", rms_error)
    flush(stdout)
    return rms_error
end

backend = CPU() # or ROCBackend() or CUDABackend() or CPU()

Δxs = [2^-5, 2^-6, 2^-7, 2^-8, 2^-9]
Ns = [length(range(-1, 1, step=dx)) for dx in Δxs]
Δt = 2^-7
# RMS_space = [diffusion_simulation(N, Δt) for N in Ns] # to actually run the simulation

# result from running the simulation
RMS_space = [0.0047825670272605945, 0.0012045108612489402, 0.00030285130786611864, 7.605028883277006e-5, 1.91992478527979e-5]


ΔTs = [2^-1, 2^-2, 2^-3, 2^-4, 2^-5]
N = 2^10
# RMS_time = [diffusion_simulation(N, Δt) for Δt in ΔTs] # to actually run the simulation

# result from running the simulation
RMS_time = [0.006427694581949278, 0.0011750223859080672, 0.0002646728101217797, 6.83103535361796e-5, 2.6586952607580914e-5]

begin
p1 = plot(Δxs, RMS_space,
    title="",
    label="Numerical Error",
    marker=:o,
    xlabel=L"\Delta x",
    ylabel=L"\Vert T_{\mathrm{ana}} - T_{\mathrm{num}} \; \Vert_{L^2}",
    yscale=:log10, xscale=:log10,
    xticks=(Δxs, [L"2^{-%$i}" for i in 5:5+length(Δxs)-1]),
    legend=:topleft,
    grid=true,
    linewidth=3, gridlinewidth=2,
    markersize=4, markerstrokewidth=0.1,
    xtickfontsize=14, ytickfontsize=14,
    xguidefontsize=16, yguidefontsize=16,
    ztickfontsize=14, zguidefontsize=16,
    legendfontsize=14,
    box=:on,
    size=(600, 400),
    fontfamily="Computer modern",
    dpi=300,
    left_margin=2Plots.mm,
)

ylim = ylims(p1)
xlim = xlims(p1)

Ns_ref = range(0.8 * Ns[1], 1.2 * Ns[end], length=100);
ref_line = Ns[1]^2.0 * 0.8 * RMS_space[1] .* (Ns_ref .^ -2.0);
Δxs_ref = range(1.2 * Δxs[1], 0.8 * Δxs[end], length=100);
ref_line_space = Δxs[end]^-2.0 * 0.8 * RMS_space[end] .* (Δxs_ref .^ 2.0);

plot!(p1, Δxs_ref, ref_line_space, linewidth=3,
    label=L"\mathcal{O}(\Delta x)^2",
    linestyle=:dash, xlims=xlim, ylims=ylim)


savefig(p1, joinpath(plots_dir(), "diffusion_spatial_convergence.pdf"))
end

begin
p2 = plot(ΔTs, RMS_time,
    title="",
    fontfamily="Computer modern",
    label="Numerical error",
    marker=:o,
    xlabel=L"\Delta t",
    ylabel=L"\Vert T_{\mathrm{ana}} - T_{\mathrm{num}} \; \Vert_{L^2}",
    yscale=:log10, xscale=:log10,
    xticks=(ΔTs, [L"2^{-%$i}" for i in 1:1+length(ΔTs)-1]),
    ms=4,
    legend=:topleft,
    grid=true,
    linewidth=3, gridlinewidth=2,
    markersize=4, markerstrokewidth=0.1,
    xtickfontsize=14, ytickfontsize=14,
    xguidefontsize=16, yguidefontsize=16,
    ztickfontsize=14, zguidefontsize=16,
    legendfontsize=14,
    box=:on,
    size=(600, 400),
    dpi=300,
    left_margin=2Plots.mm,
)
ylim = ylims(p2)
xlim = xlims(p2)

Δt_ref = range(1.2 * ΔTs[1], 0.8 * ΔTs[end], length=100);
ref_line_time = ΔTs[end]^-2.0 * 0.8 * RMS_time[end] .* (Δt_ref .^ 2.0);
plot!(p2, Δt_ref, ref_line_time, linewidth=3,
    label=L"\mathcal{O}(\Delta t)^2",
    linestyle=:dash, xlims=xlim, ylims=ylim)

savefig(p2, joinpath(plots_dir(), "diffusion_temporal_convergence.pdf"))
end