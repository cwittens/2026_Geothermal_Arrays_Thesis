# Reproduction of Brown et al. (2023), Figure 6c
# Reference: Investigating scalability of deep borehole heat exchangers:
#            Numerical modelling of arrays with varied modes of operation
# https://doi.org/10.1016/j.renene.2022.11.100
#
# This script reproduces the 30 m spacing case from Figure 6:
# a 5-DBHE line array, operated for 20 years at a constant 50 kW
# heat load per DBHE. To avoid switching the extraction direction from
# x to y, the line array is rotated so that it runs east-west, i.e. along
# the x-axis. Because the material model, domain size, and boundary
# conditions are symmetric in x and y, this is the same physical setup as
# the north-south line in the paper, only rotated by 90 degrees.

using Pkg
Pkg.activate(@__DIR__)
# Pkg.instantiate()

using GeothermalWells
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2
using KernelAbstractions: CPU
using JLD2: @save, @load
using CUDA: CUDABackend

# =============================================================================
# Setup directories
# =============================================================================
include(joinpath("..", "paths.jl"))
simulation_data_dir() = joinpath(simulation_dir(), "Brown_et_al")
!isdir(simulation_data_dir()) && mkdir(simulation_data_dir())

# =============================================================================
# Backend and precision
# =============================================================================
# Choose backend: CPU() for testing, or CUDABackend()/ROCBackend() for GPU
backend = CUDABackend()
Float_used = Float64

# =============================================================================
# Material properties
# =============================================================================
# Brown et al. Table 1. The rock is represented by the weighted-average
# homogeneous properties used in the paper.
materials = HomogenousMaterialProperties{Float_used}(
    2.55,                    # k_rock [W/(m·K)]
    2.356e6,                 # rho_c_rock [J/(m³·K)] = 2480 * 950
    0.59,                    # k_water [W/(m·K)]
    998 * 4179,              # rho_c_water [J/(m³·K)]
    52.7,                    # k_steel, outer pipe [W/(m·K)]
    7850 * 475,              # rho_c_steel [J/(m³·K)] (not specified in Brown et al.)
    0.45,                    # k_inner_pipe, polyethylene [W/(m·K)]
    941 * 1800,              # rho_c_inner_pipe [J/(m³·K)] (not specified in Brown et al.)
    1.05,                    # k_backfill / grout [W/(m·K)]
    995 * 1200               # rho_c_backfill / grout [J/(m³·K)]
)

# =============================================================================
# Borehole geometry: 5-DBHE east-west line array, 30 m spacing
# =============================================================================
# Brown et al. Table 1:
# - DBHE depth = 922 m
# - borehole diameter = 0.216 m -> borehole radius = 0.108 m
# - outer diameter of inner pipe = 0.1005 m
# - inner pipe wall thickness = 0.00688 m
# - outer pipe wall thickness = 0.0081 m
# - grout thickness = 0.01905 m
#
# Derived radii:
# - r_inner = 0.1005/2 - 0.00688 = 0.04337 m
# - r_outer = 0.108 - 0.01905 - 0.0081 = 0.08085 m

borehole_spacing = 30.0
XC = collect(-2borehole_spacing:borehole_spacing:2borehole_spacing)  # [-60, -30, 0, 30, 60] m

boreholes = tuple(
    (Borehole{Float_used}(
        xc,                      # xc [m]
        0.0,                     # yc [m]
        922.0,                   # h, borehole depth [m]
        0.04337,                 # r_inner, inner radius of central pipe [m]
        0.00688,                 # t_inner, inner pipe wall thickness [m]
        0.08085,                 # r_outer, inner radius of outer pipe [m]
        0.0081,                  # t_outer, outer pipe wall thickness [m]
        0.108,                   # r_backfill, borehole radius [m]
        998.0 * 0.005,           # m_dot [kg/s], 5 L/s water
        0.0                      # insulation_depth [m]
    ) for xc in XC)...
)

# =============================================================================
# Grid setup
# =============================================================================
# Brown et al. state a minimum model domain of 500 x 500 x 1300 m and
# a DBHE depth of about 920 m. The line array is placed along x, so the
# Figure 6 profile is extracted directly with `extract_x_profile` at y = 0.
xmin, xmax = -250.0, 250.0
ymin, ymax = -250.0, 250.0
zmin, zmax = 0.0, 1300.0

# These are numerical-resolution choices for GeothermalWells.jl, not
# directly tabulated in Brown et al.
dx_fine = 0.003       # fine spacing near borehole [m]
growth_factor = 1.3   # geometric growth rate
dx_max = 10.0         # maximum spacing far from borehole [m]
dz = 30.0             # vertical spacing [m]

gridx = create_adaptive_grid_1d(
    xmin=xmin, xmax=xmax,
    dx_fine=dx_fine, growth_factor=growth_factor, dx_max=dx_max,
    boreholes=boreholes, backend=backend, Float_used=Float_used, direction=:x
)

gridy = create_adaptive_grid_1d(
    xmin=ymin, xmax=ymax,
    dx_fine=dx_fine, growth_factor=growth_factor, dx_max=dx_max,
    boreholes=boreholes, backend=backend, Float_used=Float_used, direction=:y
)

gridz = create_uniform_gridz_with_borehole_depths(
    zmin=zmin, zmax=zmax, dz=dz,
    boreholes=boreholes, backend=backend
)

# =============================================================================
# Checkpoint/restart configuration
# =============================================================================

checkpoint_id = splitext(basename(@__FILE__))[1]
checkpoint_dir = joinpath(simulation_data_dir(), "checkpoints")

# =============================================================================
# Initial condition + restart
# =============================================================================
# Brown et al. initial condition:
# surface temperature = 9 °C, geothermal gradient = 33.4 °C/km.
# The basal heat flux reported in the paper is 85.17 mW/m², consistent with
# k_rock * gradient = 2.55 * 0.0334 ≈ 0.08517 W/m².
T0_fresh = initial_condition_thermal_gradient(
    backend, Float_used, gridx, gridy, gridz;
    T_surface=9.0,
    gradient=0.0334
)

tspan_full = (0, 3600 * 24 * 365 * 20)  # 20 years [s]
saveat_full = range(tspan_full..., 41)  # save every half year

T0, tspan, saveat = prepare_restart(
    T0_fresh, tspan_full, saveat_full;
    checkpoint_dir=checkpoint_dir,
    checkpoint_id=checkpoint_id,
    backend=backend
)

# =============================================================================
# Inlet model
# =============================================================================
# Brown et al.: P_DBHE = 50 kW per borehole, flow rate = 5 L/s.
# The package inlet model imposes the corresponding temperature drop:
# ΔT = P_DBHE / (m_dot * c_water).
P_DBHE = 50e3
c_water = 4179.0
ΔT_heat_exchanger = P_DBHE / (boreholes[1].ṁ * c_water)

inlet_model = HeatExchangerInlet{Float_used}(ΔT_heat_exchanger)

# =============================================================================
# Create simulation cache
# =============================================================================

cache = create_cache(
    backend=backend,
    gridx=gridx,
    gridy=gridy,
    gridz=gridz,
    materials=materials,
    boreholes=boreholes,
    inlet_model=inlet_model
)

# =============================================================================
# Time integration
# =============================================================================

prob = ODEProblem(rhs_diffusion_z!, T0, tspan, cache)

Δt = 160  # [s]

callback, saved_values = get_simulation_callback(
    saveat=saveat,
    print_every_n=100_000,
    checkpoint_dir=checkpoint_dir,
    checkpoint_id=checkpoint_id,
    checkpoint_every_n=300_000
)

println("Simulating Brown et al. Figure 6c setup")
println("  Array: 5 DBHE east-west line array")
println("  Spacing: $(borehole_spacing) m")
println("  Heat load: $(P_DBHE / 1e3) kW per DBHE")
println("  ΔT_heat_exchanger = $(round(ΔT_heat_exchanger, digits=4)) K")
println("  tspan = $(tspan)")
println("  Δt = $(Δt) s")
println("  Grid size: $(length(gridx)) x $(length(gridy)) x $(length(gridz))")
flush(stdout)

t_elapsed = @elapsed solve(
    prob,
    ROCK2(max_stages=100, eigen_est=eigen_estimator),
    save_everystep=false,
    callback=callback,
    adaptive=false,
    dt=Δt,
    maxiters=Int(1e10)
)

# Assemble full snapshot history from disk
reload_snapshots!(saved_values, checkpoint_dir, checkpoint_id)

println("Simulation completed in $(round(t_elapsed / 3600, digits=2)) hours")

# =============================================================================
# Save simulation data
# =============================================================================
# Create CPU cache for analysis and plotting.
cache_cpu = create_cache(
    backend=CPU(),
    gridx=gridx,
    gridy=gridy,
    gridz=gridz,
    materials=materials,
    boreholes=boreholes,
    inlet_model=inlet_model
)

data_file = joinpath(simulation_data_dir(), "Brown_et_al_line_array_30m_eastwest.jld2")
@save data_file saved_values Δt cache_cpu t_elapsed
println("Simulation data saved to $(data_file)")
