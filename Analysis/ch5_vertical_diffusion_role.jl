using Pkg
Pkg.activate(@__DIR__)

using JLD2: @load, load
using GeothermalWells
using Statistics: mean
using Plots

include(joinpath("..", "paths.jl"))
include("theme.jl")
ensure_plots_dir()


# =============================================================================
# Load the two runs: identical Brown-derived single-well setup (3000 m, Q = 150 kW,
# 20 years), differing only in whether operator B (vertical diffusion / ROCK2) is applied.
# =============================================================================
with_diffusion_file = joinpath(simulation_dir(), "array_studies_2x2", "array_study_01.jld2")
without_diffusion_file = joinpath(simulation_dir(), "no_z_diffusion_test.jld2")

saved_with = load(with_diffusion_file, "saved_values")
cache_cpu  = load(with_diffusion_file, "cache_cpu")   # grids shared by both runs
saved_without = load(without_diffusion_file, "saved_values")

# Final-time temperature fields (20 years)
T_with    = saved_with.saveval[end]
T_without = saved_without.saveval[end]

# Grid and borehole location (identical for both runs)
gridx, gridy, gridz = cache_cpu.gridx, cache_cpu.gridy, cache_cpu.gridz
xc, yc = cache_cpu.boreholes[1].xc, cache_cpu.boreholes[1].yc

begin
# =============================================================================
# ΔT = T_with − T_without  (positive ⇒ vertical diffusion keeps the rock warmer)
# Radial profiles at several depths; the 3000 m slice is the well bottom.
# =============================================================================
depths = [3000, 2000, 1000]          # [m]; 3000 m = bottom of the well
labels = ["Depth: 3000 m", "Depth: 2000 m", "Depth: 1000 m"]
linestyles = [:solid,  :solid, :dot,]

p = plot(
    xlabel = "Radial distance from borehole [m]",
    ylabel = "\$T_{\\mathrm{with}} - T_{\\mathrm{without}}\$  [°C]",
    xlims  = (-100, 100),
    legend = :topright,
    legendfontsize = 11,
    box = :on,
    gridlinewidth = 3,
    xtickfontsize = 14, ytickfontsize = 14,
    xguidefontsize = 16, yguidefontsize = 16,
    size = (600, 400),
    dpi = 300,
    fontfamily="Computer modern",
    right_margin = 3Plots.mm
)

for (i, depth) in enumerate(depths)
    r, T1 = GeothermalWells.extract_x_profile(T_with,    gridx, gridy, gridz, depth, xc, yc, true)
    _, T0 = GeothermalWells.extract_x_profile(T_without, gridx, gridy, gridz, depth, xc, yc, true)
    plot!(p, r, T1 .- T0,
        label = labels[i],
        linewidth = 4,
        linestyle = linestyles[i],
    )
end
@info savefig(p, joinpath(plots_dir(), "no_z_diffusion_delta_T.pdf"))
end


begin
# =============================================================================
# ΔT along the borehole: inlet (annulus) and outlet (inner pipe).
# Same two runs, differing only in operator B (vertical diffusion / ROCK2).
# Shows that the fluid temperatures we actually report are essentially
# untouched except right at the well bottom.
# =============================================================================
T_in_with,  T_out_with,  gridz_adj = GeothermalWells.get_temperatures_along_z_single_well(T_with,    cache_cpu)
T_in_no,    T_out_no,    _          = GeothermalWells.get_temperatures_along_z_single_well(T_without, cache_cpu)

dT_in  = T_in_with  .- T_in_no     # inlet  (annulus)
dT_out = T_out_with .- T_out_no    # outlet (inner pipe)

p_z = plot(
    xlabel = "\$T_{\\mathrm{with}} - T_{\\mathrm{without}}\$  [°C]",
    ylabel = "Depth [m]",
    legend = :bottomright,
    legendfontsize = 11,
    box = :on,
    gridlinewidth = 2,
    xtickfontsize = 14, ytickfontsize = 14,
    xguidefontsize = 16, yguidefontsize = 16,
    size = (600, 450),
    dpi = 300,
    fontfamily="Computer modern",
    right_margin = 3Plots.mm,
)

plot!(p_z, dT_in,  -gridz_adj, label = "Inlet",  color = 1, linewidth = 4)
plot!(p_z, dT_out, -gridz_adj, label = "Outlet", color = 2, linewidth = 4)

# Match the depth-axis convention of the Brown/Hu profile figures
yticks!(p_z, 0:-1000:-3000, string.(0:1000:3000))
p
# dont save it as it is not used in the thesis, but this is where the 0.02C value comes from.
end