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
simulation_Hu_path() = joinpath(simulation_dir(), "Hu_et_al_simulation_data.jld2")


# Hu et al. <https://doi.org/10.1016/j.renene.2019.09.141>
let
    saved_values_Hu = load(simulation_Hu_path(), "saved_values")
    cache_cpu_Hu = load(simulation_Hu_path(), "cache_cpu")


    # Extract grids from cache
    gridx_cpu = cache_cpu_Hu.gridx
    gridy_cpu = cache_cpu_Hu.gridy
    gridz_cpu = cache_cpu_Hu.gridz

    # =============================================================================
    # Visualization: Figure 7 - Temperature profiles at 1, 5, 10, and 25 years
    # =============================================================================


    # Create individual plots for each time point
    # Indices correspond to: 1 year, 5 years, 10 years, 25 years
    time_indices = [1, 5, 10, 25] .* 2 .+ 1
    plots_array = []

    for idx in time_indices
        # Extract temperatures along borehole depth
        T_inner, T_outer, gridz_adjusted = GeothermalWells.get_temperatures_along_z_single_well(
            saved_values_Hu.saveval[idx], cache_cpu_Hu
        )

        year = Int(saved_values_Hu.t[idx] / (3600 * 24 * 365))

        # Create temperature profile plot (inlet down, outlet up)
        p = scatter(
            vcat(T_inner, reverse(T_outer)),
            vcat(-gridz_adjusted, reverse(-gridz_adjusted)),
            label="GeothermalWells.jl",
            ylabel="Depth [m]",
            xlabel="Temperature [°C]",
            xlims=(19.519284893437614, 36.545016559411344),
            legend=false,  # Only show legend on one plot
            title="after $year year" * (year == 1 ? " " : "s"),
            markersize=4, markerstrokewidth=0.1,
            titlefont=font(16),
            box=:on,
            xtickfontsize=15, ytickfontsize=15,
            xguidefontsize=16, yguidefontsize=16,
            ztickfontsize=14, zguidefontsize=16,
            legendfontsize=10,
            gridlinewidth=2,
            size=(600, 400),
            fontfamily="Computer modern",
            dpi=300,
        )

        # Manually set y-axis ticks (since yflip affects arrow directions)
        yticks!(p, 0:-1000:-3000, string.(0:1000:3000))

        # Add Hu et al. numerical data
        depth, temperature = data_hu(year)
        scatter!(p, temperature, -depth,
            label="Numerical data (Hu et al.)",
            color=2,
            markersize=4, markerstrokewidth=0.1,
            markershape=:diamond)

        # Add inlet/outlet annotations with arrows
        color = :black
        if year == 1
            annotate!(p, [(26.8, -1590, text("Inlet", color, 14, "Computer Modern")),
                (31.2, -1500, text("Outlet", color, 14, "Computer Modern"))])
            plot!(p, [25, 27], [-1550, -1950], arrow=true, color=color, linewidth=2, label="")
            plot!(p, [33, 32.2], [-1750, -1050], arrow=true, color=color, linewidth=2, label="")
        elseif year == 5
            annotate!(p, [(26.3, -1590, text("Inlet", color, 14, "Computer Modern")),
                (29.7, -1500, text("Outlet", color, 14, "Computer Modern"))])
            plot!(p, [24.5, 26.5], [-1550, -1950], arrow=true, color=color, linewidth=2, label="")
            plot!(p, [31.5, 30.7], [-1750, -1050], arrow=true, color=color, linewidth=2, label="")
        elseif year == 10
            annotate!(p, [(25.8, -1590, text("Inlet", color, 14, "Computer Modern")),
                (29.2, -1500, text("Outlet", color, 14, "Computer Modern"))])
            plot!(p, [24, 26], [-1550, -1950], arrow=true, color=color, linewidth=2, label="")
            plot!(p, [31, 30.2], [-1750, -1050], arrow=true, color=color, linewidth=2, label="")
        elseif year == 25
            annotate!(p, [(25.8, -1590, text("Inlet", color, 14, "Computer Modern")),
                (28.5, -1500, text("Outlet", color, 14, "Computer Modern"))])
            plot!(p, [23.7, 25.7], [-1550, -1950], arrow=true, color=color, linewidth=2, label="")
            plot!(p, [30.3, 29.5], [-1750, -1050], arrow=true, color=color, linewidth=2, label="")
        end

        push!(plots_array, p)
    end

    # Adjust labels for combined plot layout
    # Remove redundant y-labels and x-labels from interior plots
    yticks!(plots_array[2], [0, -1000, -2000, -3000], [""])
    yticks!(plots_array[4], [0, -1000, -2000, -3000], [""])
    ylabel!(plots_array[2], "")
    ylabel!(plots_array[4], "")
    # xticks!(plots_array[1], [20, 25, 30, 35], [""])
    xlabel!(plots_array[1], "")
    # xticks!(plots_array[2], [20, 25, 30, 35], [""])
    xlabel!(plots_array[2], "")

    # Add legend entries (plotted first with empty data)
    p_legend = scatter([], [], label="GeothermalWells.jl", color=1,
        markersize=4, markerstrokewidth=0.1,
        legend=:top, framestyle=:none, legendfontsize=13, legend_column=2,
    )


    scatter!(p_legend, [], [], label="Numerical data (Hu et al.)",
        color=2, markersize=4, markerstrokewidth=0.1, markershape=:diamond)


    push!(plots_array, p_legend)


    # Combine into 2x2 layout
    p_combined = plot(
        plots_array...,
        layout=@layout([a b; c d; g{0.04h}]),
        size=(960, 720 + 80),
        dpi=300,
        left_margin=2Plots.mm,
        # top_margin=2Plots.mm,
        # bottom_margin=-1Plots.mm
    )

    # Save figure
    @info savefig(p_combined, joinpath(plots_dir(), "Hu_et_al_in_outlet.pdf"))
end
