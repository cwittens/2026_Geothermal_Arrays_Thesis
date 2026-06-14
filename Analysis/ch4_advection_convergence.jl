using Pkg
Pkg.activate(@__DIR__)

using GeothermalWells
import GeothermalWells: get_thermal_conductivity,
    get_volumetric_heat_capacity,
    eigen_estimator_get_dmax,
    ADI_and_ADV_callback!,
    AbstractInletModel,
    advection!
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2, DiscreteCallback
using KernelAbstractions: CPU, adapt, @kernel, @index, zeros
# using CUDA: CUDABackend
using LaTeXStrings
using Plots

include(joinpath("..", "paths.jl"))
include("theme.jl")
ensure_plots_dir()


# =============================================================================
# zero diffusion via multiply dispatch
# =============================================================================
struct NoDiffusionMaterials{RealT<:Real} <: AbstractMaterialProperties{RealT} end

@inline get_thermal_conductivity(x, y, z, boreholes, materials::NoDiffusionMaterials{RealT}) where {RealT} = zero(RealT)
@inline get_volumetric_heat_capacity(x, y, z, boreholes, materials::NoDiffusionMaterials{RealT}) where {RealT} = one(RealT)
@inline eigen_estimator_get_dmax(materials::NoDiffusionMaterials{RealT}) where {RealT} = eps(RealT)

# =============================================================================
# GPU-safe advection override for countxy_inner==0 cases
# =============================================================================
@inline function advection!(ϕ, dt, t, cache, boreholes::Tuple{Vararg{Borehole{T}}}) where {T}
    (; u_tmp, Idx_list, Idx_list_Inner, Idx_list_Outer, count_outer_per_bh, countxy_inner, countxy_outer, countz, gridx, gridy, gridz, backend, inlet_model, T_outlet, T_outlet_counter, T_turnaround_mean) = cache

    fill!(T_outlet, 0)
    fill!(T_outlet_counter, 0)
    fill!(T_turnaround_mean, 0)

    # only difference to package advection function.
    if countxy_inner > 0
        GeothermalWells.kernel_accumulate_outlet!(backend)(T_outlet, T_outlet_counter, ϕ, Idx_list_Inner, gridz, boreholes, dt, ndrange=(countxy_inner))
        T_outlet ./= T_outlet_counter
    end

    GeothermalWells.kernel_accumulate_turnaround_mean!(backend)(T_turnaround_mean, ϕ, Idx_list_Outer, gridz, count_outer_per_bh, boreholes, ndrange=(countz, countxy_outer))

    GeothermalWells.kernel_advection!(backend)(u_tmp, ϕ, gridx, gridy, gridz, Idx_list, Idx_list_Outer, T_turnaround_mean, countxy_inner, dt, t, boreholes, inlet_model, T_outlet, ndrange=(countz, countxy_inner + countxy_outer))

    GeothermalWells.kernel_copy_advection!(backend)(ϕ, u_tmp, Idx_list, ndrange=(countz, countxy_inner + countxy_outer))

    return nothing
end

# =============================================================================
# Inlet and analytical profiles
# =============================================================================
struct SineInlet{RealT<:Real} <: AbstractInletModel
    T_mean::RealT
    amplitude::RealT
    velocity::RealT
    wavelength::RealT
end

@inline T_sine(z, t, T_mean, amplitude, velocity, wavelength) =
    T_mean + amplitude * sin(2π * (z - velocity * t) / wavelength)

@inline (inlet::SineInlet)(bh_idx, T_outlet_values, t) =
    T_sine(0, t, inlet.T_mean, inlet.amplitude, inlet.velocity, inlet.wavelength)




@kernel function sine_advection_kernel!(ϕ, gridz, t, T_mean, amplitude, velocity, wavelength)
    k, j, i = @index(Global, NTuple)
    ϕ[k, j, i] = T_sine(gridz[k], t, T_mean, amplitude, velocity, wavelength)
end


# =============================================================================
# Setup helper for a 1D-in-z outer-pipe advection problem
# =============================================================================
function outer_pipe_borehole(FloatT, h, velocity)
    r_outer = 1.0
    ρ_water = 998.2
    ṁ = velocity * π * r_outer^2 * ρ_water
    # borehole struct does: v_outer = ṁ / (A_outer * ρ_water)
    return Borehole{FloatT}(
        0.0,        # xc
        0.0,        # yc
        h,          # h
        0.0,        # r_inner
        0.0,        # t_inner
        r_outer,    # r_outer
        0.0,        # t_outer
        r_outer,    # r_backfill
        ṁ,          # mass flow rate chosen so v_outer == velocity
        0.0,        # insulation_depth
    )
end

function outer_pipe_grid(FloatT, backend, zmax, Δz)
    # This is basically a 1D problem in z. The 2x2 cross-section is only there
    # because the package's ADI needs at least two x/y grid points.
    # All four xy nodes are classified as outer-pipe cells.
    gridx = adapt(backend, FloatT[-0.25, 0.25])
    gridy = adapt(backend, FloatT[-0.25, 0.25])
    gridz = adapt(backend, (FloatT(0):FloatT(Δz):FloatT(zmax)))

    return gridx, gridy, gridz
end


function advection_simulation(Δz, Δt, backend=backend)
    FloatT = Float64
    velocity = 1.0
    zmax = 10
    t_end = 5
    T_mean = 20.0
    amplitude = 3.0
    wavelength = 0.5

    materials = NoDiffusionMaterials{FloatT}()
    borehole = outer_pipe_borehole(FloatT, zmax, velocity)
    boreholes = (borehole,)

    gridx, gridy, gridz = outer_pipe_grid(FloatT, backend, zmax, Δz)

    ϕ_convergence_initial = zeros(backend, FloatT, length(gridz), length(gridy), length(gridx))
    sine_advection_kernel!(backend)(ϕ_convergence_initial, gridz, 0.0, T_mean, amplitude, velocity, wavelength, ndrange=size(ϕ_convergence_initial))

    inlet = SineInlet(T_mean, amplitude, velocity, wavelength)
    cache = create_cache(backend=backend, gridx=gridx, gridy=gridy, gridz=gridz, materials=materials, boreholes=boreholes, inlet_model=inlet)

    callback = DiscreteCallback((u, t, integrator) -> true, ADI_and_ADV_callback!, save_positions=(false, false))

    tspan = (0.0, t_end)

    prob = ODEProblem(rhs_diffusion_z!, ϕ_convergence_initial, tspan, cache)
    println("Starting simulation for Δz = $Δz and dt = $Δt")
    @time sol_z = solve(prob, ROCK2(max_stages=100, eigen_est=eigen_estimator), save_everystep=false, callback=callback, adaptive=false, dt=Δt)

    # overwrite the initial condition to be the analytical solution at t = t_end
    ϕ_convergence_analytical = zeros(backend, FloatT, length(gridz), length(gridy), length(gridx))
    sine_advection_kernel!(backend)(ϕ_convergence_analytical, gridz, t_end, T_mean, amplitude, velocity, wavelength, ndrange=size(ϕ_convergence_analytical))


    ϕ_num_cpu = adapt(CPU(), sol_z[end])
    ϕ_ana_cpu = adapt(CPU(), ϕ_convergence_analytical)

    rms_error = sqrt(sum((ϕ_num_cpu .- ϕ_ana_cpu) .^ 2) / length(ϕ_num_cpu))
    println("rms error:", rms_error)
    return rms_error
end

backend = CPU() # or ROCBackend() or CUDABackend() or CPU()

Δzs = [2^-5, 2^-6, 2^-7, 2^-8, 2^-9, 2^-10, 2^-11, 2^-12, 2^-13, 2^-14, 2^-15]
cfl_substep = 0.37
Δts = [16 * cfl_substep * Δz for Δz in Δzs]
# RMS_space = [advection_simulation(Δz, Δt) for (Δz, Δt) in zip(Δzs, Δts)] # to actually run the simulation

# result from running the simulation
RMS_space = [2.3342483182670066, 1.6420741287476908, 0.9764416606561254, 0.5336765221832545, 0.27907582931952246, 0.14276488825204736, 0.07222434876713346, 0.0363400343014762, 0.018225517767242464, 0.009127164250281745, 0.004567282549925141]


orders = [log(RMS_space[i] / RMS_space[i-1]) / log(Δzs[i] / Δzs[i-1]) for i in 2:length(Δzs)]

begin
    p1 = plot(Δzs, RMS_space,
        title="",
        label="Numerical error",
        marker=:o,
        xlabel=L"\Delta z",
        ylabel=L"\Vert T_{\mathrm{ana}} - T_{\mathrm{num}} \; \Vert_{L^2}",
        yscale=:log10, xscale=:log10,
        xticks=(Δzs[1:2:end], [L"2^{-%$i}" for i in 5:2:5+length(Δzs)-1]),
        ms=4,
        legend=:topleft,
        grid=true,
        linewidth=3, gridlinewidth=2,
        markersize=5, markerstrokewidth=0.1,
        xtickfontsize=14, ytickfontsize=14,
        xguidefontsize=16, yguidefontsize=16,
        ztickfontsize=14, zguidefontsize=16,
        legendfontsize=14,
        box=:on,
        size=(600, 400),
        dpi=300,
        fontfamily="Computer modern",
        left_margin=2Plots.mm,
    )

    ylim = ylims(p1)
    xlim = xlims(p1)

    Δzs_ref = range(1.2 * Δzs[1], 0.8 * Δzs[end], length=100)
    ref_line_space = Δzs[end]^-1.0 * 0.8 * RMS_space[end] .* (Δzs_ref .^ 1.0)

    plot!(p1, Δzs_ref, ref_line_space, linewidth=3,
        label=L"\mathcal{O}(\Delta z)",
        linestyle=:dash, xlims=xlim, ylims=ylim)

    savefig(p1, joinpath(plots_dir(), "advection_spatial_convergence.pdf"))
end