using Pkg
Pkg.activate(@__DIR__)

using GeothermalWells
using KernelAbstractions: CPU
using Plots

include(joinpath("..", "paths.jl"))
include("theme.jl")
ensure_plots_dir()

function make_boreholes(XC, YC)
    return tuple(
        (Borehole{Float64}(
            xc,
            yc,
            2000.0,
            0.0511,
            0.0114,
            0.0885,
            0.00833,
            0.115,
            42 * 998.2 / 3600,
            0.0,
        ) for xc in XC, yc in YC)...,
    )
end

function make_grids(boreholes)
    xmin, xmax = -100.0, 100.0
    ymin, ymax = -100.0, 100.0
    dx_fine = 0.0025
    growth_factor = 1.3
    dx_max = 10.0

    gridx = create_adaptive_grid_1d(
        xmin=xmin, xmax=xmax,
        dx_fine=dx_fine, growth_factor=growth_factor, dx_max=dx_max,
        boreholes=boreholes, backend=CPU(), Float_used=Float64, direction=:x,
    )

    gridy = create_adaptive_grid_1d(
        xmin=ymin, xmax=ymax,
        dx_fine=dx_fine, growth_factor=growth_factor, dx_max=dx_max,
        boreholes=boreholes, backend=CPU(), Float_used=Float64, direction=:y,
    )

    return gridx, gridy
end

function generate_grid_overview(filename, XC, YC; zoom_center=(XC[1], YC[1]), zoom_halfwidth=0.17)
    boreholes = make_boreholes(XC, YC)
    gridx, gridy = make_grids(boreholes)

    default(
    grid=true,
    box=:on,
    dpi=300,
    titlefont=font(16),
    linewidth=1.2, gridlinewidth=2,
    markersize=4, markerstrokewidth=0.1,
    xtickfontsize=10, ytickfontsize=10,
    xguidefontsize=13, yguidefontsize=13,
    ztickfontsize=14, zguidefontsize=16,
    legendfontsize=10,
    fontfamily="Computer modern",
)

    p1 = plot_grid(
        gridx, gridy;
        size=(450, 450),
        boreholes=boreholes,
        legend=false,
        annotate=false,
    )

    p2 = plot_grid(
        gridx, gridy;
        size=(450, 450),
        boreholes=boreholes,
        legend=:bottomleft,
        xlims=(zoom_center[1] - zoom_halfwidth, zoom_center[1] + zoom_halfwidth),
        ylims=(zoom_center[2] - zoom_halfwidth, zoom_center[2] + zoom_halfwidth),
        annotate=false,
    )

    p = plot(p1, p2; layout=(1, 2), size=(900, 450), dpi=300, left_margin=4Plots.mm)

    savefig(joinpath(plots_dir(), filename))
end

generate_grid_overview("Grid_overview_single.pdf", [0.0], [0.0]; zoom_center=(0.0, 0.0))
generate_grid_overview("Grid_overview_2x2.pdf", [-30, 30.0], [-30, 30.0]; zoom_center=(30.0, 30.0))
generate_grid_overview("Grid_overview_3x3.pdf", [-60, 0, 60], [-60, 0, 60]; zoom_center=(60.0, 60.0))

