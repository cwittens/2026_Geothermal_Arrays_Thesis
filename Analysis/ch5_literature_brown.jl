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
simulation_Brown_path() = joinpath(simulation_dir(), "Brown_et_al", "Brown_et_al_simulation_data.jld2")
simulation_Brown_Array_path(spacing) = joinpath(simulation_dir(), "Brown_et_al", "Brown_et_al_line_array_$(spacing)m_eastwest.jld2")


# Brown et al. <https://doi.org/10.1016/j.renene.2022.11.100>
let
    saved_values_Brown = load(simulation_Brown_path(), "saved_values")
    cache_cpu_Brown = load(simulation_Brown_path(), "cache_cpu")


    # =============================================================================
    # Visualization: Figure 4b - Temperature profiles along borehole depth
    # =============================================================================

    # Load reference data from Brown et al. Figure 4b
    (z_beier, T_beier), (z_brown, T_brown) = data_brown_single_well_b()

    # Extract temperatures along borehole depth from simulation
    T_inner, T_outer, gridz_adjusted = GeothermalWells.get_temperatures_along_z_single_well(saved_values_Brown.saveval[end], cache_cpu_Brown)

    move_back = 0.1
    # Create plot comparing inlet/outlet temperatures
    p1 = plot(legend=false)
    plot!(p1, Shape([1.4 - move_back, 4, 4, 1.4 - move_back], [-890, -890, -720, -720]),
        fillcolor=:white, linecolor=:black, label="")

    scatter!(p1, [1.55 - move_back], [-770], color=1, markersize=6, label="")
    scatter!(p1, [1.55 - move_back], [-840], color=2, markershape=:diamond, markersize=6, label="")

    annotate!(p1, 1.7 - move_back, -770, text("GeothermalWells.jl", :left, 11, "Computer Modern"))
    annotate!(p1, 1.7 - move_back, -840, text("Numerical data (Brown et al.)", :left, 11, "Computer Modern"))


    scatter!(p1,
        T_inner, -gridz_adjusted,
        label="GeothermalWells.jl",
        color=1,
        ylabel="Depth [m]",
        xlabel="Temperature [°C]",
        markersize=4,
        xtickfontsize=14, ytickfontsize=14,
        xguidefontsize=16, yguidefontsize=16,
        ztickfontsize=14, zguidefontsize=16,
        box=:on,
        gridlinewidth=2,
        size=(600, 400),
        fontfamily="Computer modern",
        dpi=300,
        bottom_margin=2Plots.mm,
        legend=false,
    )
    scatter!(p1, T_outer, -gridz_adjusted, label="", color=1, markersize=4,)

    # Add Brown et al. numerical data
    scatter!(p1, T_brown, -z_brown,
        label="Numerical data (Brown et al.)",
        color=2,
        markershape=:diamond,
        markersize=4,
    )

    # Add annotations to identify inlet and outlet curves
    color = :black
    annotate!(p1, [(2.85, -300, text("Inlets", color, 14, "Computer Modern")),
        (4.7, -260, text("Outlets", color, 14, "Computer Modern"))])

    # Add arrows pointing to the curves
    plot!(p1, [2.4, 3.2], [-300, -500], arrow=true, color=color, linewidth=2, label="")
    plot!(p1, [4.32, 4.12], [-300, -80], arrow=true, color=color, linewidth=2, label="")

    # Manually set y-axis ticks (since yflip affects arrow directions)
    yticks!(p1, -0:-200:-800, string.(0:200:800))

    # Save figure
    @info savefig(p1, joinpath(plots_dir(), "Brown_et_al_in_outlet.pdf"))


    # =============================================================================
    # Visualization: Figure 4d - Radial temperature profiles at different depths
    # =============================================================================

    # Extract grids and final temperature field
    gridx_cpu = cache_cpu_Brown.gridx
    gridy_cpu = cache_cpu_Brown.gridy
    gridz_cpu = cache_cpu_Brown.gridz
    xc = cache_cpu_Brown.boreholes[1].xc
    yc = cache_cpu_Brown.boreholes[1].yc

    T = saved_values_Brown.saveval[end]             # final temperature field [°C]
    t = saved_values_Brown.t[end] / 3600 / 24 / 365  # final time [years]

    # Depths at which to compare with Brown et al. data (matches their Figure 4d)
    depths = [300.0, 600.0, 920.0]  # [m]
    colors = SERIES

    # fake legend, bottom-left corner — coords now in (radial distance, temperature)
    lx0, lx1 = -94.0, 29      # box left / right
    ly0, ly1 = -1.8, 5.8      # box bottom / top
    row_top, row_bot = 3.5, 0.5  # y of the two entries
    sw_x = lx0 + 7              # x of swatches
    txt_x = lx0 + 14             # x where text starts
    p2 = plot(legend=false)
    plot!(p2, Shape([lx0, lx1, lx1, lx0], [ly0, ly0, ly1, ly1]),
        fillcolor=:white, linecolor=:black, label="")

    # line swatch for the simulation curve
    plot!(p2, [sw_x - 5, sw_x + 5], [row_top, row_top], color=:black, linewidth=4, label="")
    # diamond swatch for Brown et al.
    scatter!(p2, [sw_x + 0.1], [row_bot], color=:black, markershape=:diamond, markersize=8, label="")

    annotate!(p2, txt_x, row_top, text("GeothermalWells.jl", :left, 13, "Computer Modern"))
    annotate!(p2, txt_x, row_bot, text("Numerical data (Brown et al.)", :left, 13, "Computer Modern"))

    # Create comparison plot
    plot!(p2,
        xlabel="Radial distance from borehole [m]",
        ylabel="Temperature of rock layers [°C]",
        xlims=(-100, 100),
        ylims=(-3, 42),
        box=:on,
        gridlinewidth=2,
        legend=:bottomleft,
        legendfontsize=12,
        markersize=4,
        xtickfontsize=14, ytickfontsize=14,
        xguidefontsize=16, yguidefontsize=16,
        ztickfontsize=14, zguidefontsize=16,
        left_margin=1Plots.mm,
        right_margin=3Plots.mm,
        size=(600, 450),
        fontfamily="Computer modern",
        dpi=300,)

    # Add legend entries (plotted first with empty data)
    # plot!(p2, [], [], label="GeothermalWells.jl", color=:black, linewidth=2)
    # scatter!(p2, [], [], label="Numerical data (Brown et al.)",
    #     color=:black, markersize=3, markerstrokewidth=0, markershape=:diamond)

    # Plot simulation results for each depth
    for (i, depth) in enumerate(depths)
        # Extract radial temperature profile at this depth
        # Note: Using full_profile=true to get both sides of borehole
        r, T_profile = GeothermalWells.extract_x_profile(T, gridx_cpu, gridy_cpu, gridz_cpu, depth, xc, yc, true)
        plot!(p2, r, T_profile,
            label="",
            color=colors[i],
            linewidth=4.5)
    end

    # Add Brown et al. numerical data as scatter points
    for i in 1:length(depths)
        r_num, T_num = data_brown_single_well_d(i)
        scatter!(p2, r_num, T_num,
            label="",
            color=colors[i],
            markersize=4.3,
            markershape=:diamond)
    end

    # Add depth annotations to identify each curve
    # Positions chosen to avoid overlap with curves
    annotate!(p2, 65 + 5, 37.2, text("920m", color=p2.series_list[8-2][:linecolor], :left, 15, "Computer Modern"))
    annotate!(p2, 65 + 5, 26.5, text("600m", color=p2.series_list[7-2][:linecolor], :left, 15, "Computer Modern"))
    annotate!(p2, 65 + 5, 16.5, text("300m", color=p2.series_list[6-2][:linecolor], :left, 15, "Computer Modern"))

    # Save figure
    @info savefig(p2, joinpath(plots_dir(), "Brown_et_al_radial_profile.pdf"))


    # =============================================================================
    # Visualization: Figure 6 - Radial temperature profiles for arrays at  different depths
    # =============================================================================



    P = []
    for spacing in [20, 30, 40, 50]

        saved_values_array = load(simulation_Brown_Array_path(spacing), "saved_values")
        cache_cpu_array = load(simulation_Brown_Array_path(spacing), "cache_cpu")

        # Extract grids and final temperature field
        gridx_cpu = cache_cpu_array.gridx
        gridy_cpu = cache_cpu_array.gridy
        gridz_cpu = cache_cpu_array.gridz
        xc = cache_cpu_array.boreholes[3].xc
        yc = cache_cpu_array.boreholes[1].yc

        T = saved_values_array.saveval[end]             # final temperature field [°C]
        t = saved_values_array.t[end] / 3600 / 24 / 365  # final time [years]


        # Depths at which to compare with Brown et al. data (matches their Figure 4d)
        depths = [300, 600, 920]  # [m]
        colors = SERIES

        if spacing == 20
            xlims = (-110, 110)
        elseif spacing == 30
            xlims = (-130, 130)
        elseif spacing == 40
            xlims = (-150, 150)
        else
            xlims = (-200, 200)
        end
        # Create comparison plot
        p = plot(
            xlabel="Distance from central borehole [m]",
            ylabel="Temperature [°C]",
            title="Array Spacing of $spacing m",
            xlims=xlims,
            ylims=(-3, 42),
            box=:on,
            gridlinewidth=2,
            legend=:bottomleft,
            legendfontsize=10,
            markersize=3,
            xtickfontsize=14, ytickfontsize=14,
            xguidefontsize=16, yguidefontsize=16,
            ztickfontsize=14, zguidefontsize=16,
            # size=(600, 450),
            fontfamily="Computer modern",
            dpi=300
        )


        # Plot simulation results for each depth
        for (i, depth) in enumerate(depths)
            # Extract radial temperature profile at this depth
            # Note: Using full_profile=true to get both sides of borehole
            r, T_profile = GeothermalWells.extract_x_profile(T, gridx_cpu, gridy_cpu, gridz_cpu, depth, xc, yc, true)
            plot!(p, r, T_profile,
                label="",
                color=colors[i],
                linewidth=4)
        end

        # Add Brown et al. numerical data as scatter points
        for i in 1:length(depths)
            r_num, T_num = data_brown_array(spacing, i)
            scatter!(p, r_num, T_num,
                label="",
                color=colors[i],
                markersize=4,
                markerstrokewidth=0.1,
                markershape=:diamond)
        end

        # Add depth annotations to identify each curve
        # Positions chosen to avoid overlap with curves
        if spacing == 50
            annotate!(p, 150, 37.2 - 1, text("920m", color=p.series_list[3][:linecolor], :left, 13, "Computer Modern"))
            annotate!(p, 150, 26.5 - 1, text("600m", color=p.series_list[2][:linecolor], :left, 13, "Computer Modern"))
            annotate!(p, 150, 16.5 - 1, text("300m", color=p.series_list[1][:linecolor], :left, 13, "Computer Modern"))
        end
        push!(P, p)
    end


    # Add legend entries (plotted first with empty data)
    p_legend = plot([], [], label="GeothermalWells.jl", color=:black, linewidth=2,
        legend=:top, framestyle=:none, legendfontsize=14, legend_column=2,
    )

    # plot!(p_legend, [], [],label = "920m", color=3, linewidth=2)

    scatter!(p_legend, [], [], label="Numerical data (Brown et al.)",
        color=:black, markersize=3, markerstrokewidth=0, markershape=:diamond)


    push!(P, p_legend)


    xlabel!(P[1], "")
    xlabel!(P[2], "")
    ylabel!(P[2], "")
    ylabel!(P[4], "")
    p_array = plot(P..., layout=@layout([a b; c d; g{0.04h}]), size=(960, 720), left_margin=3Plots.mm,
        right_margin=3Plots.mm,
    )
    @info savefig(p_array, joinpath(plots_dir(), "Brown_et_al_array_radial_profile.pdf"))
end
