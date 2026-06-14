using Pkg
Pkg.activate(@__DIR__)

using JLD2
using GeothermalWells
using Plots

include(joinpath("..", "paths.jl"))
include("theme.jl")
ensure_plots_dir()


function plot_simulation_times(plot_idx, nx1_idx=[])

    configs = [(1, 1), (2, 2), (3, 1), (1, 3), (5, 1), (1, 5), (3, 3), (10, 1), (1, 10), (3, 4), (4, 4), (4, 5), (5, 5), (6, 6), (7, 7), (8, 8), (9, 9), (10, 10)]

    times_10d = [
        24.15,
        55.05,
        50.28,
        50.43,
        77.05,
        79.03,
        104.93,
        161.77,
        160.74,
        136.85,
        182.86,
        237.56,
        316.37,
        536.67,
        864.74,
        1351.7,
        2027.16,
        3058.0
    ]

    # new Fable 5 improved version of GeothermalWells.jl (not included in Master Thesis)
    times_10d_new = [7.7, 17.94, 15.86, 15.69, 23.38, 23.01, 32.45, 43.09, 43.04, 40.29, 50.29, 62.24, 74.21, 101.87, 132.11, 169.76, 210.1, 256.66,]

    # new Fable 5 improved version_2 of GeothermalWells.jl (not included in Master Thesis)
    times_10d_new = [6.87, 14.18, 12.35, NaN, 19.44, NaN, 24.01, 30.86, NaN, NaN, 41.88, 49.28, 55.69, 75.49, 97.19, 130.85, 161.63, 189.97,]

    # times_10d = times_10d_new

    times_1y = 1 * times_10d .* (365 / 10) # extrapolate to 1 year
    times_1y_new = 1 * times_10d_new .* (365 / 10) # extrapolate to 1 year

    n_wells = [r * c for (r, c) in configs]
    labels = ["$(r)×$(c)" for (r, c) in configs]

    # Sort by number of wells for a clean line plot

    n_wells_sorted = n_wells[plot_idx]
    times_sorted = times_1y[plot_idx]
    times_sorted_new = times_1y_new[plot_idx]
    labels_sorted = labels[plot_idx]

    # Plot total simulation time
    p = plot(n_wells_sorted, times_sorted ./ 3600,
        xlabel="Number of Wells",
        ylabel="Time [hours]",
        label="Wall-clock Time for 1 Year of Operation",
        fontfamily="Computer modern",
        legend=:topleft,
        marker=:o, markersize=6,
        grid=true, box=:on,
        size=(600, 400), dpi=300,
        linewidth=3, gridlinewidth=2,
        xtickfontsize=14, ytickfontsize=14,
        xguidefontsize=16, yguidefontsize=16,
        ztickfontsize=14, zguidefontsize=16,
        legendfontsize=13,)

    # plot!(n_wells_sorted, times_sorted_new ./ 3600, label="Fable 5 Improvement", marker=:o, markersize=6, color = 3, lw=3)


    max_time = 4maximum(times_sorted ./ 3600)
    offset = max_time * 0.02

    if length(plot_idx) < 10
        # Annotate each point with its array configuration
        for i in eachindex(n_wells_sorted)
            # Alternate above/below to reduce overlap
            if i != 1
                annotate!(p, n_wells_sorted[i], times_sorted[i] / 3600 - offset,
                    text(labels_sorted[i], 11, :bottom))
            else
                annotate!(p, n_wells_sorted[i], times_sorted[i] / 3600 + offset,
                    text(labels_sorted[i], 11, :top))
            end
        end

    else # only plot 1x1, 3x3, 5x5, 10x10 and annotate them
        for i in plot_idx
            if configs[i] in [(5, 5), (6, 6), (7, 7), (8, 8), (9, 9), (10, 10)]
                annotate!(p, n_wells[i], times_1y[i] / 3600 - offset,
                    text(labels[i], 11, :bottom))
            elseif configs[i] in [(1, 1), (3, 3),]
                annotate!(p, n_wells[i], times_1y[i] / 3600 + offset,
                    text(labels[i], 11, :top))
            end
        end

    end



    n_wells_nx1 = n_wells[nx1_idx]
    times_nx1 = times_1y[nx1_idx]
    times_nx1_new = times_1y_new[nx1_idx]
    labels_nx1 = labels[nx1_idx]

    scatter!(p, n_wells_nx1, times_nx1 ./ 3600,
        marker=:o, markersize=6, color=2,
        label="",
    )

    # scatter!(p, n_wells_nx1, times_nx1_new ./ 3600,
    #     marker=:o, markersize=6, color=3,
    #     label="",
    # )

    for i in eachindex(n_wells_nx1)
        annotate!(p, n_wells_nx1[i], times_nx1[i] / 3600 + offset,
            text(labels_nx1[i], 11, :top))
    end



    return p
end



plot_idx = [1, 2, 7, 10, 11, 12, 13]
nx1_idx = [3, 5, 8]

p_25 = plot_simulation_times(plot_idx, nx1_idx)

p_100 = plot_simulation_times([plot_idx..., 14, 15, 16, 17, 18])

@info savefig(p_25, joinpath(plots_dir(), "simulation_time_vs_number_of_wells_25.pdf"))
@info savefig(p_100, joinpath(plots_dir(), "simulation_time_vs_number_of_wells_100.pdf"))


begin
    # Percentages of total profiler wall time
    N = 8
    n_wells = (1:N) .^ 2
    thomas = [75.9, 69.9, 55.4, 42.9, 32.1, 24.4, 18.5, 14.2, 11.1, 11.9][1:N]
    advection = [3.9, 11.0, 22.4, 38.2, 52.2, 62.5, 71.1, 77.0, 81.7, 82.3][1:N]
    adi_explicit = [3.2, 5.5, 6.1, 5.9, 5.2, 4.4, 3.6, 3.1, 2.5, 2.1][1:N]
    rock2 = [2.9, 5.2, 5.8, 5.7, 4.9, 4.2, 3.5, 2.9, 2.4, 2.0][1:N]
    other = [3.4, 4.5, 4.8, 4.6, 4.0, 3.4, 2.8, 2.3, 1.9, 1.6][1:N]
    gpu_util = [89.3, 95.9, 94.6, 97.2, 98.3, 98.9, 99.5, 99.4, 99.8, 99.8][1:N]

    # new Fable 5 improved version of GeothermalWells.jl (not included in Master Thesis)
    # thomas = [46.7, 44.0, 38.8, 34.9, 30.8, 29.1, 27.8, 26.5, 25.8, 25.4][1:N]
    # advection = [6.4, 4.7, 4.0, 3.7, 3.6, 3.6, 3.6, 3.6, 3.6, 3.6][1:N]
    # adi_explicit = [10.2, 16.8, 20.0, 21.8, 22.3, 23.4, 24.1, 24.5, 24.8, 25.1][1:N]
    # rock2 = [9.3, 16.1, 18.9, 20.8, 21.2, 22.3, 22.9, 23.3, 23.6, 23.9][1:N]
    # other = [9.7, 13.4, 15.5, 16.5, 17.0, 17.6, 18.3, 18.5, 18.9, 19.1][1:N]
    # gpu_util = [82.2, 95.0, 97.1, 97.6, 94.9, 96.0, 96.6, 96.4, 96.6, 96.9][1:N]

    # Cumulative boundaries for stacked fill plot
    y_thomas = thomas
    y_advection = y_thomas .+ advection
    y_adi_explicit = y_advection .+ adi_explicit
    y_rock2 = y_adi_explicit .+ rock2
    y_other = y_rock2 .+ other

    configs_label = [(1, 1), (3, 3), (4, 4), (5, 5), (6, 6), (7, 7), (8, 8), (9, 9), (10, 10)][1:N-1]
    n_wells_label = [r * c for (r, c) in configs_label]
    labels = ["$(r)×$(c)" for (r, c) in configs_label]
    labels[1] = "   1x1"
    labels[end] = "$N×$N       "
    p = plot(
        xlabel="Array Configurations",
        ylabel="Fraction of Total Wall Time [%]",
        fontfamily="Computer modern",
        legend=:right,
        grid=true, box=:on,
        size=(600, 400), dpi=300,
        linewidth=3, gridlinewidth=2,
        xtickfontsize=12, ytickfontsize=13,
        xguidefontsize=14, yguidefontsize=14,
        legendfontsize=11,
        xlims=(1, N^2),
        ylims=(0, 100),
        xticks=(n_wells_label, labels),)

    colors = SERIES
    x_dummy = [-100.0, -99.0]
    y_dummy = [-100.0, -100.0]

    # Dummy legend series
    plot!(p, x_dummy, y_dummy,
        fillrange=y_dummy, fillalpha=0.8,
        linewidth=3,
        color=colors[1],
        label="ADI: Thomas Solves",
    )

    plot!(p, x_dummy, y_dummy,
        fillrange=y_dummy, fillalpha=0.8,
        linewidth=3,
        color=colors[2],
        label="Advection",
    )

    plot!(p, x_dummy, y_dummy,
        fillrange=y_dummy, fillalpha=0.8,
        linewidth=3,
        color=colors[3],
        label="ADI: Explicit Diffusion",
    )

    plot!(p, x_dummy, y_dummy,
        fillrange=y_dummy, fillalpha=0.8,
        linewidth=3,
        color=colors[4],
        label="Vertical Diffusion",
    )

    plot!(p, x_dummy, y_dummy,
        fillrange=y_dummy, fillalpha=0.8,
        linewidth=3,
        color=colors[5],
        label="Other Kernels",
    )

    # Actual areas: reverse drawing order
    plot!(p, n_wells, y_other,
        fillrange=y_rock2, fillalpha=0.8,
        linewidth=3, marker=:o, markersize=5,
        color=colors[5], label=false,
    )

    plot!(p, n_wells, y_rock2,
        fillrange=y_adi_explicit, fillalpha=0.8,
        linewidth=3, marker=:o, markersize=5,
        color=colors[4], label=false,
    )

    plot!(p, n_wells, y_adi_explicit,
        fillrange=y_advection, fillalpha=0.8,
        linewidth=3, marker=:o, markersize=5,
        color=colors[3], label=false,
    )

    plot!(p, n_wells, y_advection,
        fillrange=y_thomas, fillalpha=0.8,
        linewidth=3, marker=:o, markersize=5,
        color=colors[2], label=false,
    )

    plot!(p, n_wells, y_thomas,
        fillrange=0, fillalpha=0.8,
        linewidth=3, marker=:o, markersize=5,
        color=colors[1], label=false,
    )

    @info savefig(p, joinpath(plots_dir(), "kernel_time_breakdown.pdf"))
end
