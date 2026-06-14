using Pkg
Pkg.activate(@__DIR__)

# Plot ROCK2 stability regions used in the thesis.
using OrdinaryDiffEqStabilizedRK
using SciMLBase: ODEProblem, solve, remake
using Plots

include(joinpath("..", "paths.jl"))
include("theme.jl")
ensure_plots_dir()


for s in [5, 9] # , 20

    # ROCK2 stability interval is approximately 0.81 s^2 for damping η = 0.95
    # Abdulle & Medovikov report this behavior for their damped ROCK2 polynomials.
    L = 0.81 * s^2

    # Dahlquist test equation: u' = z u, with h = 1, so z = hλ = λ
    function f!(du, u, z, t)
        du[1] = z * u[1]
    end

    prob = ODEProblem(f!, ComplexF64[1.0+0im], (0.0, 1.0), 0.0 + 0im)

    # Force a fixed ROCK2 degree/stage count.
    # eigen_est is supplied because ROCK2 normally uses spectral-radius estimates.
    eigen_est = integrator -> (integrator.eigen_est = abs(integrator.p))
    alg = ROCK2(min_stages=s, max_stages=s, eigen_est=eigen_est)

    function rock2_absR(z)
        sol = solve(
            remake(prob; p=z),
            alg;
            dt=1.0,
            adaptive=false,
            save_everystep=false,
            save_start=false,
        )
        return abs(sol.u[end][1])
    end

    # Plot window. Make the imaginary range much smaller than the real range:
    # ROCK2 is designed mainly for the negative real axis.

    if s == 5
        ymax = 4
    elseif s == 9
        ymax = 4
    elseif s == 20
        ymax = 8
    end

    xs = range(-1.01L, 0.01L, length=5 * 500)
    ys = range(-ymax, ymax, length=5 * 260)

    A = Matrix{Float64}(undef, length(ys), length(xs))

    @info "Computing ROCK2 stability grid for s = $s ..."
    for i in eachindex(ys), j in eachindex(xs)
        z = xs[j] + im * ys[i]
        A[i, j] = rock2_absR(z)
    end

    # Plot |R(z)| <= 1.

    # --- parameters ---
    # s = 20
    L = 0.81 * s^2   # approximate real-axis extent

    # stable mask
    stable = Float64.(A .<= 1.0)


    p = heatmap(
        xs, ys, stable;
        # xlims = (-1.01L, 0.5),
        # ylims = (-ymax, ymax),
        box=:on,
        clim=(0, 1),
        colorbar=false,
        c=cgrad([:white, colorant"#9B001F"]),   # red fill
        xlabel="Re(z)",
        ylabel="Im(z)",
        # title = "$s stages",
        legend=false,
        grid=false,
        framestyle=:box,
        size=(800, 320),
        dpi=300,
        left_margin=4Plots.mm,
        bottom_margin=6Plots.mm,
        linewidth=3, gridlinewidth=2,
        markersize=4, markerstrokewidth=0.1,
        xtickfontsize=14, ytickfontsize=14,
        xguidefontsize=16, yguidefontsize=16,
        ztickfontsize=14, zguidefontsize=16,
        legendfontsize=14,
        fontfamily="Computer modern",
        # guidefontsize = 13,
        # tickfontsize = 11,
        # titlefontsize = 16,
    )

    # Stability boundary
    contour!(
        p, xs, ys, stable;
        levels=[0.99999],
        color=:black,
        linewidth=2.2,
    )

    # Horizontal axis
    # hline!(
    #     p, [0.0];
    #     color = :gray40,
    #     linestyle = :dash,
    #     linewidth = 1.2,
    # )
    savefig(joinpath(plots_dir(), "rock2_stability_region_s_$s.pdf"))
end

