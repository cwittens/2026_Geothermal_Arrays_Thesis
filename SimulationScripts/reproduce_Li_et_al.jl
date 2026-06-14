# Reproduction of Li et al. (2021) validation case
# Reference: Heat extraction model and characteristics of coaxial deep borehole heat exchanger
# https://doi.org/10.1016/j.renene.2021.01.036
# Compares simulation results with numerical data from Figure 5 at 500m, 1000m, 1500m, and 2000m depths

using Pkg
Pkg.activate(@__DIR__)

using GeothermalWells
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2
using KernelAbstractions: CPU, adapt
using Statistics: mean
using JLD2: @save, @load
using CUDA: CUDABackend

# =============================================================================
# Setup directories
# =============================================================================
include(joinpath("..", "paths.jl"))
simulation_data_dir() = simulation_dir()
!isdir(simulation_data_dir()) && mkdir(simulation_data_dir())

# =============================================================================
# Backend and precision
# =============================================================================
# Choose backend: CPU() for testing, or CUDABackend()/ROCBackend() for GPU
backend = CUDABackend()
Float_used = Float64

# =============================================================================
# Material properties (stratified rock)
# =============================================================================
# Material properties from Li et al. Tables 1 and 2.
#
# Li et al. give the four layer properties and prescribe the initial
# geothermal profile through Eq. (11). In that formula the layer bottoms are
# denoted by H_j, but the actual H_j values are not listed in the paper.
#
# The missing layer depths can be reconstructed from the far-field temperatures
# in their Figure 5. Far away from the borehole, the rock is essentially still
# undisturbed after 120 days, so these values should correspond to the initial
# geothermal profile. Using Eq. (11) with
#
#     T_surface = 10 °C,
#     q_g       = 0.075 W/m²,
#     k         = (1.8, 2.6, 3.5, 5.3) W/(m K),
#
# and the far-field temperatures at 500, 1000, 1500, and 2000 m gives a layer
# spacing of about 400 m if we assume equally spaced interfaces,
#
#     H1 = d, H2 = 2d, H3 = 3d,  with  d ≈ 400 m.
#
# This gives the rounded layer bottoms (400, 800, 1200, 2200) m used below.
# If the equal-spacing assumption is dropped, the reconstructed interfaces are
# approximately (401, 796, 1204) m instead, so this choice has no visible effect
# on the comparison.
layer_depths = (400.0, 800.0, 1200.0, 2200.0)

materials_stratified = StratifiedMaterialProperties{4,Float_used}(
    (1.8, 2.6, 3.5, 5.3),                                # k_rock [W/(m·K)]
    (1780 * 1379, 2030 * 1450, 1510 * 1300, 2600 * 878), # rho_c_rock [J/(m³·K)]
    layer_depths,                                        # layer_depths [m]
    0.618,                   # k_water [W/(m·K)]
    4.166e6,                 # rho_c_water [J/(m³·K)]
    41.0,                    # k_steel [W/(m·K)]
    7850 * 475,              # rho_c_steel [J/(m³·K)] (estimated - not specified in paper)
    0.4,                     # k_insulating [W/(m·K)]
    1.955e6,                 # rho_c_insulating [J/(m³·K)] (estimated - not specified in paper)
    1.5,                     # k_backfill [W/(m·K)]
    1.76e6                   # rho_c_backfill [J/(m³·K)] (estimated - not specified in paper)
)

# =============================================================================
# Borehole geometry
# =============================================================================
# Deep coaxial borehole heat exchanger (2000m depth) based on Li et al. Table 2
# No insulation on outer pipe
borehole = Borehole{Float_used}(
    0.0,                      # xc [m]
    0.0,                      # yc [m]
    2000.0,                   # h - borehole depth [m]
    0.0511,                   # r_inner - inner pipe radius [m] (calculated from (125-2×11.4)/2)
    0.0114,                   # t_inner - inner pipe wall thickness [m] (11.4 mm)
    0.0885,                   # r_outer - outer pipe inner radius [m] (calculated from (193.7-2×8.33)/2)
    0.00833,                  # t_outer - outer pipe wall thickness [m] (8.33 mm)
    0.115,                    # r_backfill - outer radius [m] (estimated - not specified in paper)
    42 * 998.2 / 3600,        # ṁ - mass flow rate [kg/s] (42 m³/h converted)
    0.0                       # insulation_depth [m] (no insulation)
)

boreholes = (borehole,)

# =============================================================================
# Grid setup
# =============================================================================
# Domain boundaries
xmin, xmax = -100.0, 100.0
ymin, ymax = -100.0, 100.0
zmin, zmax = 0.0, 2200.0

# Grid parameters
dx_fine = 0.0025      # fine spacing near borehole [m]
growth_factor = 1.3   # geometric growth rate
dx_max = 10.0         # maximum spacing far from borehole [m]
dz = 100.0            # vertical spacing [m]

# Create adaptive grids (fine near borehole, coarse far away)
# The grids are first created on CPU because the custom initial condition below
# is easiest to evaluate there. They are adapted to the chosen backend afterward.
gridx_cpu = create_adaptive_grid_1d(
    xmin=xmin, xmax=xmax,
    dx_fine=dx_fine, growth_factor=growth_factor, dx_max=dx_max,
    boreholes=boreholes, backend=CPU(), Float_used=Float_used, direction=:x
)

gridy_cpu = create_adaptive_grid_1d(
    xmin=ymin, xmax=ymax,
    dx_fine=dx_fine, growth_factor=growth_factor, dx_max=dx_max,
    boreholes=boreholes, backend=CPU(), Float_used=Float_used, direction=:y
)

gridz_cpu = create_uniform_gridz_with_borehole_depths(
    zmin=zmin, zmax=zmax, dz=dz,
    boreholes=boreholes, backend=CPU()
)

# =============================================================================
# Initial condition + restart
# =============================================================================
# Undisturbed geothermal profile from Li et al. Eq. (11).
# The temperature is uniform in the radial direction and varies only with depth.

function segment_length(z, zlo, zhi)
    return max(0.0, min(z, zhi) - zlo)
end

function initial_temperature_Li(z)
    z = max(z, 0.0)

    T_surface = 10.0   # [°C]
    qg = 0.075         # [W/m²]

    k1 = 1.8           # [W/(m·K)]
    k2 = 2.6
    k3 = 3.5
    k4 = 5.3

    return T_surface +
           qg / k1 * segment_length(z, 0.0, 400.0) +
           qg / k2 * segment_length(z, 400.0, 800.0) +
           qg / k3 * segment_length(z, 800.0, 1200.0) +
           qg / k4 * max(z - 1200.0, 0.0)
end

# Check that the reconstructed initial condition agrees with the far-field
# temperatures digitized from Li et al. Figure 5.
n = 5

_, T500  = data_li(1)
_, T1000 = data_li(2)
_, T1500 = data_li(3)
_, T2000 = data_li(4)

@assert abs(initial_temperature_Li(500)  - mean(T500[end-n:end]))  < 0.05
@assert abs(initial_temperature_Li(1000) - mean(T1000[end-n:end])) < 0.05
@assert abs(initial_temperature_Li(1500) - mean(T1500[end-n:end])) < 0.05
@assert abs(initial_temperature_Li(2000) - mean(T2000[end-n:end])) < 0.05


T0_fresh = adapt(
    backend,
    [Float_used(initial_temperature_Li(z)) for z in gridz_cpu, y in gridy_cpu, x in gridx_cpu]
)

# Adapt grids to backend
gridx = adapt(backend, gridx_cpu)
gridy = adapt(backend, gridy_cpu)
gridz = adapt(backend, gridz_cpu)

# =============================================================================
# Checkpoint/restart configuration
# =============================================================================
checkpoint_id = splitext(basename(@__FILE__))[1]
checkpoint_dir = joinpath(simulation_data_dir(), "checkpoints")

tspan_full = (0, 3600 * 24 * 120)  # 120 days [s]
saveat_full = range(tspan_full..., 2)  # save initial and final state

T0, tspan, saveat = prepare_restart(
    T0_fresh, tspan_full, saveat_full;
    checkpoint_dir=checkpoint_dir,
    checkpoint_id=checkpoint_id,
    backend=backend
)

# =============================================================================
# Inlet model
# =============================================================================
# Heat exchanger inlet: T_inlet = T_outlet - Q / C
# where Q is heat extraction rate and C is heat capacity flow rate
Q = 200e3                      # heat extraction rate [W]
c_water = 4174                 # specific heat of water [J/(kg·K)]
inlet_model = HeatExchangerInlet{Float_used}(Q / (borehole.ṁ * c_water))

# =============================================================================
# Create simulation cache
# =============================================================================
cache = create_cache(
    backend=backend,
    gridx=gridx,
    gridy=gridy,
    gridz=gridz,
    materials=materials_stratified,
    boreholes=boreholes,
    inlet_model=inlet_model
)

# =============================================================================
# Time integration
# =============================================================================
prob = ODEProblem(rhs_diffusion_z!, T0, tspan, cache)

# Time step and solver
Δt = 160 # [s]

callback, saved_values = get_simulation_callback(
    saveat=saveat,
    print_every_n=10_000,
    checkpoint_dir=checkpoint_dir,
    checkpoint_id=checkpoint_id,
    checkpoint_every_n=100_000
)

println("Simulating with Δt = $(Δt)s, tspan = $(tspan)")
println("Grid size: $(length(gridx)) x $(length(gridy)) x $(length(gridz))")
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

println("Simulation completed in $(t_elapsed / 60) minutes")
# Simulation completed in 4.2739420605000005 minutes

# =============================================================================
# Save simulation data
# =============================================================================
# Create CPU cache for analysis (makes it easier to not have to deal with GPU arrays)
cache_cpu = create_cache(
    backend=CPU(),
    gridx=gridx,
    gridy=gridy,
    gridz=gridz,
    materials=materials_stratified,
    boreholes=boreholes,
    inlet_model=inlet_model
)

@save joinpath(simulation_data_dir(), "Li_et_al_simulation_data.jld2") saved_values Δt cache_cpu t_elapsed

println("Simulation data saved to $(simulation_data_dir())")