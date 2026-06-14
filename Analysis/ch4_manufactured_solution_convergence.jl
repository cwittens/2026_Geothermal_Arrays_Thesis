using Pkg
Pkg.activate(@__DIR__)

using GeothermalWells
import GeothermalWells: get_thermal_conductivity,
    get_volumetric_heat_capacity,
    eigen_estimator_get_dmax,
    ADI_and_ADV_callback!,
    advection!
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2, DiscreteCallback
using DiffEqCallbacks: CallbackSet
using KernelAbstractions: CPU, adapt, @kernel, @index, zeros
# using CUDA: CUDABackend
using LaTeXStrings
using Plots

include(joinpath("..", "paths.jl"))
include("theme.jl")
ensure_plots_dir()

# =============================================================================
# Variable coefficients k(x,y,z) and ρc(x,y,z) via dispatch
# =============================================================================
struct MMSMaterials{RealT<:Real} <: AbstractMaterialProperties{RealT}
    k0::RealT
    ax::RealT
    ay::RealT
    az::RealT
    rx::RealT
    ry::RealT
    rz::RealT
end

@inline get_thermal_conductivity(x, y, z, boreholes, m::MMSMaterials) =
    m.k0 * (1 + m.ax * x + m.ay * y + m.az * z)

@inline get_volumetric_heat_capacity(x, y, z, boreholes, m::MMSMaterials) =
    1 + m.rx * x + m.ry * y + m.rz * z

@inline function eigen_estimator_get_dmax(m::MMSMaterials)
    # max k / min ρc on the unit cube (coefficients may be negative)
    kmax = m.k0 * (1 + max(m.ax, 0) + max(m.ay, 0) + max(m.az, 0))
    ρcmin = 1 + min(m.rx, 0) + min(m.ry, 0) + min(m.rz, 0)
    return kmax / ρcmin
end

# =============================================================================
# GPU-safe advection override for countxy_inner == 0 (same as advection test)
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
# Exact solution
# =============================================================================
@inline q_mms(s) = 64 * s^3 * (1 - s)^3

@inline T_exact(x, y, z, t, Tm, A, λ) = Tm + A * exp(-λ * t) * q_mms(x) * q_mms(y) * q_mms(z)

@kernel function exact_solution_kernel!(ϕ, @Const(gridx), @Const(gridy), @Const(gridz), t, Tm, A, λ)
    k, j, i = @index(Global, NTuple)
    ϕ[k, j, i] = T_exact(gridx[i], gridy[j], gridz[k], t, Tm, A, λ)
end

# =============================================================================
# Manufactured source term
# =============================================================================

# calculated the source terms using:
#=
using Symbolics

@variables x y z t
Dt = Differential(t)
Dx = Differential(x)
Dy = Differential(y)
Dz = Differential(z)

@variables λ A Tm k0 ax ay az rx ry rz vz

q(s) = 64 * s^3 * (1-s)^3
# choose q(s) too have no problems at the boundary,
# because q(0)=q(1)=q'(0)=q'(1)=q''(0)=q''(1)=0.

k = k0 * (1 + ax * x + ay * y + az * z)
ρc = 1 + rx*x + ry*y + rz*z

T = Tm + A*exp(-λ*t) * q(x)*q(y)*q(z)
S = Dt(T) - (1/ρc)*(Dx(k*Dx(T)) + Dy(k*Dy(T)) + Dz(k*Dz(T))) + vz*Dz(T)

S_expr = expand_derivatives(S)
=#

# results from symbolics calculation:
@inline function source_julia(t, x, y, z, λ, A, k0, ax, ay, az, rx, ry, rz, vz)
    return ((-(786432A * (x^2) * ((1 - x)^3) * (y^3) * ((1 - y)^3) * (z^3) * ((1 - z)^3) * exp(-t * λ) - 786432A * ((1 - x)^2) * (x^3) * (y^3) * ((1 - y)^3) * (z^3) * ((1 - z)^3) * exp(-t * λ)) * ax * k0 - (-4718592A * (x^3) * ((1 - x)^3) * (y^3) * ((1 - y)^3) * (z^2) * ((1 - z)^2) * exp(-t * λ) + 1572864A * (x^3) * ((1 - x)^3) * (y^3) * ((1 - y)^3) * z * ((1 - z)^3) * exp(-t * λ) + 1572864A * (x^3) * ((1 - x)^3) * (y^3) * ((1 - y)^3) * (1 - z) * (z^3) * exp(-t * λ)) * (1 + ax * x + ay * y + az * z) * k0 - (-4718592A * (x^3) * ((1 - x)^3) * (y^2) * ((1 - y)^2) * (z^3) * ((1 - z)^3) * exp(-t * λ) + 1572864A * (x^3) * ((1 - x)^3) * y * ((1 - y)^3) * (z^3) * ((1 - z)^3) * exp(-t * λ) + 1572864A * (x^3) * ((1 - x)^3) * (1 - y) * (y^3) * (z^3) * ((1 - z)^3) * exp(-t * λ)) * (1 + ax * x + ay * y + az * z) * k0 - (-4718592A * (x^2) * ((1 - x)^2) * (y^3) * ((1 - y)^3) * (z^3) * ((1 - z)^3) * exp(-t * λ) + 1572864A * x * ((1 - x)^3) * (y^3) * ((1 - y)^3) * (z^3) * ((1 - z)^3) * exp(-t * λ) + 1572864A * (1 - x) * (x^3) * (y^3) * ((1 - y)^3) * (z^3) * ((1 - z)^3) * exp(-t * λ)) * (1 + ax * x + ay * y + az * z) * k0 - (786432A * (x^3) * ((1 - x)^3) * (y^2) * ((1 - y)^3) * (z^3) * ((1 - z)^3) * exp(-t * λ) - 786432A * (x^3) * ((1 - x)^3) * ((1 - y)^2) * (y^3) * (z^3) * ((1 - z)^3) * exp(-t * λ)) * ay * k0 - (786432A * (x^3) * ((1 - x)^3) * (y^3) * ((1 - y)^3) * (z^2) * ((1 - z)^3) * exp(-t * λ) - 786432A * (x^3) * ((1 - x)^3) * (y^3) * ((1 - y)^3) * ((1 - z)^2) * (z^3) * exp(-t * λ)) * az * k0) / (1 + rx * x + ry * y + rz * z) + (786432A * (x^3) * ((1 - x)^3) * (y^3) * ((1 - y)^3) * (z^2) * ((1 - z)^3) * exp(-t * λ) - 786432A * (x^3) * ((1 - x)^3) * (y^3) * ((1 - y)^3) * ((1 - z)^2) * (z^3) * exp(-t * λ)) * vz - 262144A * (x^3) * ((1 - x)^3) * (y^3) * ((1 - y)^3) * (z^3) * ((1 - z)^3) * exp(-t * λ) * λ)
end

@kernel function source_kernel!(ϕ, @Const(gridx), @Const(gridy), @Const(gridz),
                                t_eval, dt, λ, A, k0, ax, ay, az, rx, ry, rz, vz)
    k, j, i = @index(Global, NTuple)
    ϕ[k, j, i] += dt * source_julia(t_eval, gridx[i], gridy[j], gridz[k],
                                    λ, A, k0, ax, ay, az, rx, ry, rz, vz)
end

function make_source_callback(λ, A, k0, ax, ay, az, rx, ry, rz, vz)
    function source_affect!(integrator)
        (; backend, gridx, gridy, gridz, Nx, Ny, Nz) = integrator.p
        dt = integrator.t - integrator.tprev
        t_mid = integrator.tprev + dt / 2 # midpoint quadrature of the source integral
        source_kernel!(backend)(integrator.u, gridx, gridy, gridz, t_mid, dt,
                                λ, A, k0, ax, ay, az, rx, ry, rz, vz,
                                ndrange=(Nz, Ny, Nx))
        return nothing
    end
    return DiscreteCallback((u, t, integrator) -> true, source_affect!, save_positions=(false, false))
end

# =============================================================================
# Borehole covering the whole unit square so all cells are outer-pipe cells
# =============================================================================
function mms_borehole(FloatT, velocity)
    r_outer = 1.0 # corner of [0,1]^2 has r = sqrt(0.5) < 1, so all cells are outer pipe
    ρ_water = 998.2
    ṁ = velocity * π * r_outer^2 * ρ_water # so that v_outer == velocity
    return Borehole{FloatT}(
        0.5,     # xc
        0.5,     # yc
        1.0,     # h == zmax, must be in gridz
        0.0,     # r_inner
        0.0,     # t_inner
        r_outer, # r_outer
        0.0,     # t_outer
        r_outer, # r_backfill
        ṁ,
        0.0,     # insulation_depth
    )
end


# =============================================================================
# MMS simulation
# =============================================================================
function mms_simulation(Δ, Δt, backend=backend; v = 1.0)
    FloatT = Float64

    # parameters (must match the ones baked into validation below)
    Tm = 20.0
    A = 1.0
    λ = 1.0
    k0 = 5e-3
    ax, ay, az = 0.5, 0.3, 0.4
    rx, ry, rz = 0.4, 0.2, 0.3
    t_end = 1.0

    materials = MMSMaterials{FloatT}(k0, ax, ay, az, rx, ry, rz)
    borehole = mms_borehole(FloatT, v)
    @assert borehole.v_outer ≈ v
    boreholes = (borehole,)

    gridx = adapt(backend, collect(FloatT, 0:Δ:1))
    gridy = adapt(backend, collect(FloatT, 0:Δ:1))
    gridz = adapt(backend, collect(FloatT, 0:Δ:1)) # contains h = 1.0 exactly for Δ = 2^-n

    Nx, Ny, Nz = length(gridx), length(gridy), length(gridz)

    ϕ0 = zeros(backend, FloatT, Nz, Ny, Nx)
    exact_solution_kernel!(backend)(ϕ0, gridx, gridy, gridz, 0.0, Tm, A, λ, ndrange=(Nz, Ny, Nx))

    inlet = ConstantInlet(Tm)
    cache = create_cache(backend=backend, gridx=gridx, gridy=gridy, gridz=gridz,
                         materials=materials, boreholes=boreholes, inlet_model=inlet)

    ADI_cb = DiscreteCallback((u, t, integrator) -> true, ADI_and_ADV_callback!, save_positions=(false, false))
    source_cb = make_source_callback(λ, A, k0, ax, ay, az, rx, ry, rz, v)
    callback = CallbackSet(ADI_cb, source_cb) # order matters: source after ADI + advection

    prob = ODEProblem(rhs_diffusion_z!, ϕ0, (0.0, t_end), cache)
    println("Starting MMS simulation for Δ = $Δ and dt = $Δt")
    @time sol = solve(prob, ROCK2(max_stages=100, eigen_est=eigen_estimator),
                      save_everystep=false, callback=callback, adaptive=false, dt=Δt)

    ϕ_exact = zeros(backend, FloatT, Nz, Ny, Nx)
    exact_solution_kernel!(backend)(ϕ_exact, gridx, gridy, gridz, t_end, Tm, A, λ, ndrange=(Nz, Ny, Nx))

    ϕ_num_cpu = adapt(CPU(), sol[end])
    ϕ_ex_cpu = adapt(CPU(), ϕ_exact)

    rms_error = sqrt(sum((ϕ_num_cpu .- ϕ_ex_cpu) .^ 2) / length(ϕ_num_cpu))
    println("rms error: ", rms_error)
    return rms_error
end

backend = CPU() # or ROCBackend() or CUDABackend() or CPU()

Δs = [2^-5, 2^-6, 2^-7, 2^-8, 2^-9, 2^-10]
Δts = [16 * Δ for Δ in Δs] # v = 1: each of the 4 semi-Lagrangian substeps jumps 4 cells
# RMS = [mms_simulation(Δ, Δt) for (Δ, Δt) in zip(Δs, Δts)] # to actually run the simulation

# result from running the simulation
RMS = [0.11645697622450904, 0.03999090569623872, 0.01894645825955841, 0.00930306508241787, 0.0046144429031319215, 0.002298508566624391]
orders = [log(RMS[i] / RMS[i-1]) / log(Δs[i] / Δs[i-1]) for i in 2:length(Δs)]
# orders = [1.542053205019253, 1.0777437704325936, 1.026150159969764, 1.01154963619842, 1.005458440178769]



begin
p1 = plot(Δs, RMS,
    title="",
    label="Numerical error",
    marker=:o,
    xlabel=L"\Delta x",
    ylabel=L"\Vert T_{\mathrm{man}} - T_{\mathrm{num}} \; \Vert_{L^2}",
    yscale=:log10, xscale=:log2,
    xticks=(Δs[1:1:end], [L"2^{-%$i}" for i in 5:1:5+length(Δs)-1]),
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

Δzs_ref = range(1.2 * Δs[1], 0.8 * Δs[end], length=100);
ref_line_space = Δs[end]^-1.0 * 0.86 * RMS[end] .* (Δzs_ref .^ 1.0);

plot!(p1, Δzs_ref, ref_line_space, linewidth=3,
    label=L"\mathcal{O}(\Delta x)",
    linestyle=:dash, xlims=xlim, ylims=ylim)
    
@info savefig(p1, joinpath(plots_dir(),  "manufactured_solution_convergence.pdf"))
end


# for 2nd order convergence (no advection aka v=0.0)

Δs = [2^-3, 2^-4, 2^-5, 2^-6, 2^-7]
# RMS_no_advection = [mms_simulation(Δ, 2^-11, v=0.005) for Δ in Δs] # to actually run the simulation

# result from running the simulation
RMS_no_advection = [0.003458950891894509, 0.0010373470200951845, 0.000266052023352635, 6.666546842876306e-5, 1.72073260429515e-5]
orders_no_advection = [log(RMS_no_advection[i] / RMS_no_advection[i-1]) / log(Δs[i] / Δs[i-1]) for i in 2:length(Δs)]
# orders_no_advection = [1.7374359366850733, 1.9631183141368076, 1.996696806737184, 1.9539167384194078]

