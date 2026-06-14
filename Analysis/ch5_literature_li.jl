using Pkg
Pkg.activate(@__DIR__)

using GeothermalWells
using KernelAbstractions: CPU, adapt
using Statistics: mean
using JLD2: load
using Plots

include(joinpath("..", "paths.jl"))
include("theme.jl")
ensure_plots_dir()

# =============================================================================
# Setup directories
# =============================================================================
simulation_Li_path() = joinpath(simulation_dir(), "Li_et_al_simulation_data.jld2")


# Li et al. <https://doi.org/10.1016/j.renene.2021.01.036>
let
    saved_values_Li = load(simulation_Li_path(), "saved_values")
    cache_cpu_Li = load(simulation_Li_path(), "cache_cpu")

    # Extract grids and final temperature field
    gridx_cpu = cache_cpu_Li.gridx
    gridy_cpu = cache_cpu_Li.gridy
    gridz_cpu = cache_cpu_Li.gridz
    xc = cache_cpu_Li.boreholes[1].xc
    yc = cache_cpu_Li.boreholes[1].yc

    T = saved_values_Li.saveval[end]           # final temperature field [°C]
    t = saved_values_Li.t[end] / 3600 / 24    # final time [days]

    # Depths at which to compare with Li et al. data (matches their Figure 5)
    depths = [500.0, 1000.0, 1500.0, 2000.0]
    colors = SERIES

    # =============================================================================
    # Create comparison plot
    # =============================================================================
    p = plot(
        xlabel="Radial distance from borehole [m]",
        ylabel="Temperature of rock layers [°C]",
        xlims=(-2, 60),
        ylims=(20, 60),
        legend=:bottomright,
        box=:on,
        gridlinewidth=2,
        xtickfontsize=14, ytickfontsize=14,
        xguidefontsize=16, yguidefontsize=16,
        legendfontsize=12,
        size=(600, 450),
        fontfamily="Computer modern",
        dpi=300
    )

    # Add legend entries (plotted first with empty data)
    plot!(p, [], [], label="GeothermalWells.jl", color=:black, linewidth=2)
    scatter!(p, [], [], label="Numerical data (Li et al.)",
        color=:black, markersize=3, markerstrokewidth=0, markershape=:diamond)

    # Plot simulation results for each depth
    for (i, depth) in enumerate(depths)
        # Extract radial temperature profile at this depth
        r, T_profile = GeothermalWells.extract_x_profile(T, gridx_cpu, gridy_cpu, gridz_cpu, depth, xc, yc)
        plot!(p, r, T_profile,
            label="",
            color=colors[i],
            linewidth=4)
    end

    # Add Li et al. numerical data as scatter points
    for i in 1:length(depths)
        r_num, T_num = data_li(i)
        scatter!(p, r_num, T_num,
            label="",
            color=colors[i],
            markersize=4,
            markershape=:diamond)
    end

    # Add depth annotations to identify each curve
    # Positions chosen to avoid overlap with curves
    annotate!(p, 50, 58.1 - 2, text("2000m", color=p.series_list[6][:linecolor], :left, 15))
    annotate!(p, 50, 51.0 - 2, text("1500m", color=p.series_list[5][:linecolor], :left, 15))
    annotate!(p, 50, 42.5 - 2, text("1000m", color=p.series_list[4][:linecolor], :left, 15))
    annotate!(p, 50, 29.6 + 2, text("500m", color=p.series_list[3][:linecolor], :left, 15))

    # Save figure
    @info savefig(p, joinpath(plots_dir(), "Li_et_al_radial_profile.pdf"))
end
