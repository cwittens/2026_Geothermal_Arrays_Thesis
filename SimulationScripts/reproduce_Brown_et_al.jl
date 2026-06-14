# Reproduction of Brown et al. (2023)
# Reference: Investigating scalability of deep borehole heat exchangers: 
#            Numerical modelling of arrays with varied modes of operation
# https://doi.org/10.1016/j.renene.2022.11.100
# Compares simulation results with numerical data from Figures 4b and 4d

using Pkg
Pkg.activate(@__DIR__)

using GeothermalWells
using OrdinaryDiffEqStabilizedRK: ODEProblem, solve, ROCK2
using KernelAbstractions: CPU, adapt
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
# Material properties (homogeneous rock)
# =============================================================================
# Material properties from Brown et al. Table 1
# Using single homogeneous layer since paper uses weighted average
materials = HomogenousMaterialProperties{Float_used}(
    2.55,                    # k_rock - thermal conductivity of rock [W/(m·K)]
    2.356e6,                 # rho_c_rock - volumetric heat capacity [J/(m³·K)]
    0.59,                    # k_water [W/(m·K)]
    998 * 4179,              # rho_c_water [J/(m³·K)]
    52.7,                    # k_steel (outer pipe) [W/(m·K)]
    7850 * 475,              # rho_c_steel [J/(m³·K)] (estimated - not specified in paper)
    0.45,                    # k_insulating (polyethylene inner pipe) [W/(m·K)]
    941 * 1800,              # rho_c_insulating [J/(m³·K)] (estimated - not specified in paper)
    1.05,                    # k_backfill (grout) [W/(m·K)]
    995 * 1200               # rho_c_backfill (grout) [J/(m³·K)]
)

# =============================================================================
# Borehole geometry
# =============================================================================
# Newcastle borehole geometry from Brown et al. Table 1
# Geometry calculations:
# - Inner pipe outer diameter = 0.1005 m → outer radius = 0.05025 m
# - Inner pipe thickness = 0.00688 m → inner radius = 0.05025 - 0.00688 = 0.04337 m
# - Borehole diameter = 0.216 m → borehole radius = 0.108 m
# - Grout thickness = 0.01905 m → r_outer + t_outer = 0.108 - 0.01905 = 0.08895 m
# - Outer pipe thickness = 0.0081 m → r_outer = 0.08895 - 0.0081 = 0.08085 m

borehole = Borehole{Float_used}(
    0.0,                     # xc [m]
    0.0,                     # yc [m]
    922.0,                   # h - borehole depth [m]
    0.04337,                 # r_inner - inner radius of central pipe [m]
    0.00688,                 # t_inner - thickness of inner pipe wall [m]
    0.08085,                 # r_outer - inner radius of outer pipe [m]
    0.0081,                  # t_outer - thickness of outer pipe wall [m]
    0.108,                   # r_backfill - borehole radius [m]
    998.0 * 0.005,           # ṁ - mass flow rate [kg/s] (5 L/s converted)
    0.0                      # insulation_depth [m] (no insulation mentioned in paper)
)

boreholes = (borehole,)

# =============================================================================
# Grid setup
# =============================================================================
# Domain boundaries
xmin, xmax = -100.0, 100.0
ymin, ymax = -100.0, 100.0
zmin, zmax = 0.0, 1100.0

# Grid parameters
dx_fine = 0.003       # fine spacing near borehole [m]
growth_factor = 1.3   # geometric growth rate
dx_max = 10           # maximum spacing far from borehole [m]
dz = 30               # vertical spacing [m]

# Create adaptive grids (fine near borehole, coarse far away)
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
# Linear thermal gradient from Brown et al.
# T(z) = T_surface + gradient * z
T0_fresh = initial_condition_thermal_gradient(
    backend, Float_used, gridx, gridy, gridz;
    T_surface=9.0,      # surface temperature [°C]
    gradient=0.0334     # thermal gradient [K/m]
)

tspan_full = (0, 3600 * 24 * 365 * 20)  # 20 years [s]
saveat_full = range(tspan_full..., 41)   # save every half year

T0, tspan, saveat = prepare_restart(
    T0_fresh, tspan_full, saveat_full;
    checkpoint_dir=checkpoint_dir,
    checkpoint_id=checkpoint_id,
    backend=backend
)

# =============================================================================
# Inlet model
# =============================================================================
# Heat exchanger inlet from Brown et al.
# P_DBHE = 50 kW, ṁ = 4.99 kg/s, c_water = 4179 J/(kg·K)
# ΔT = P / (ṁ * c) = 50000 / (4.99 * 4179) ≈ 2.398 K
Q = 50e3                             # heat extraction rate [W]
c_water = 4179.0                     # specific heat of water [J/(kg·K)]
inlet_model = HeatExchangerInlet{Float_used}(Q / (borehole.ṁ * c_water))

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

# Time step and solver
Δt = 160  # [s]

callback, saved_values = get_simulation_callback(
    saveat=saveat,
    print_every_n=100_000,
    checkpoint_dir=checkpoint_dir,
    checkpoint_id=checkpoint_id,
    checkpoint_every_n=300_000
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

println("Simulation completed in $(t_elapsed / 3600) hours")
# Simulation completed in 3.6831196272019446 hours

# =============================================================================
# Save simulation data
# =============================================================================
# Create CPU cache for analysis (makes it easier to not have to deal with GPU arrays)
cache_cpu = create_cache(
    backend=CPU(),
    gridx=gridx,
    gridy=gridy,
    gridz=gridz,
    materials=materials,
    boreholes=boreholes,
    inlet_model=inlet_model
)

@save joinpath(simulation_data_dir(), "Brown_et_al_simulation_data.jld2") saved_values Δt cache_cpu t_elapsed

println("Simulation data saved to $(simulation_data_dir())")
