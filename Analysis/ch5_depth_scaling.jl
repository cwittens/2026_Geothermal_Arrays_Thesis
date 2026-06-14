using Pkg
Pkg.activate(@__DIR__)

using JLD2: @load
using GeothermalWells
using Statistics: mean
using LaTeXStrings
using LsqFit: curve_fit
using Plots

include(joinpath("..", "paths.jl"))
include("theme.jl")
ensure_plots_dir()

# =============================================================================
# Load simulation results
# =============================================================================
data_dir = joinpath(simulation_dir(), "scaling_study")

files = [
    joinpath(data_dir, "scaling_Q_study_01.jld2"),
    joinpath(data_dir, "scaling_Q_study_02.jld2"),
    joinpath(data_dir, "scaling_Q_study_03.jld2"),
    joinpath(data_dir, "scaling_Q_study_04.jld2"),
    joinpath(data_dir, "scaling_Q_study_05.jld2"),
    joinpath(data_dir, "scaling_Q_study_06.jld2"),
    joinpath(data_dir, "scaling_Q_study_07.jld2"),
    joinpath(data_dir, "scaling_Q_study_08.jld2"),
    joinpath(data_dir, "scaling_Q_study_09.jld2"),
    joinpath(data_dir, "scaling_Q_study_10.jld2"),
]
 
depths = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]
 
# from simulation setup
T_inlet = 20.0
ṁ = 998.0 * 0.005  # kg/s
c_water = 4179.0   # J/(kg·K)
 
T_outlets = Float64[]
ΔTs = Float64[]
Qs_kW = Float64[]
 
for (i, file) in enumerate(files)
    @load file saved_values cache_cpu
 
    # Final snapshot
    T_final = saved_values.saveval[end]
    t_final = saved_values.t[end]
 
    # Extract outlet temperature (inner pipe at surface)
    T_inner, T_outer, gz = GeothermalWells.get_temperatures_along_z_single_well(T_final, cache_cpu)
    T_out = T_inner[1]
 
    ΔT = T_out - T_inlet
    Q = ṁ * c_water * ΔT / 1000  # kW
 
    push!(T_outlets, T_out)
    push!(ΔTs, ΔT)
    push!(Qs_kW, Q)
 
    println("Depth $(depths[i])m: T_outlet = $(round(T_out, digits=2))°C, " *
            "ΔT = $(round(ΔT, digits=2))K, Q = $(round(Q, digits=2)) kW " *
            "(t_final = $(round(t_final / 3600 / 24 / 365, digits=1)) years)")
end


# =============================================================================
# Fit models 
# =============================================================================
mask = 000 .<= depths .<= 3000
d_fit = depths[mask]
Q_fit = Qs_kW[mask]
 
# Model 1: Q = a*d^2 + b*d 
# Model 1 changed: Q = a*d^2 + 0
model_quad(d, p) = p[1] .* d.^2 .+ p[2] .* d .* 0
fit_quad = curve_fit(model_quad, d_fit, Q_fit, [3e-5, -0.02])
a, b = fit_quad.param

# # Model 2: Q = c * d^α
# model_power(d, p) = p[1] .* d.^p[2]
# fit_power = curve_fit(model_power, d_fit, Q_fit, [1e-4, 2.0])
# c_pow, α = fit_power.param
# println("Power law fit: Q = $(c_pow) d^$(round(α, digits=3))")


# =============================================================================
# Plot
# =============================================================================
d_dense = range(minimum(depths), maximum(depths), 200)

p1 = plot(
    grid=true,
    box=:on,
    size=(600, 400),
    dpi=300,
    titlefont=font(16),
    linewidth=3, gridlinewidth=2,
    markersize=4, markerstrokewidth=0.1,
    xtickfontsize=14, ytickfontsize=14,
    xguidefontsize=16, yguidefontsize=16,
    ztickfontsize=14, zguidefontsize=16,
    legendfontsize=14,
    fontfamily="Computer modern",
)

scatter!(p1, depths, Qs_kW, 
    markersize=6, color=1, label="Simulation data")
plot!(p1, d_dense, model_quad(d_dense, fit_quad.param), z_order=1,
    linewidth=3, linestyle=:dash, color=2,
    label=L"Quadratic fit: $ad^2$ (a=%$(round(a, sigdigits=3)))")#, b=$(round(b, sigdigits=3)))")

# plot!(p1, d_dense, model_power(d_dense, fit_power.param),
#     linewidth=2, linestyle=:dot, color=:orange,
#     label="Fit: c·d^α (α=$(round(α, digits=2)))")

plot!(p1,
    xlabel="Borehole Depth [m]",
    ylabel="Heat Extraction Rate [kW]",
    # title="Heat Extraction vs Borehole Depth (20-year, constant inlet)",
    legend=:topleft,
    grid=true,
)

savefig(p1, joinpath(plots_dir(), "scaling_Q_with_fits.pdf"))
