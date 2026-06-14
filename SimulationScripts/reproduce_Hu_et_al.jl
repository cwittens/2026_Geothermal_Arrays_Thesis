# Reproduction of Hu et al. (2020) validation case
# Reference: Numerical modeling of a coaxial borehole heat exchanger to exploit 
#            geothermal energy from abandoned petroleum wells in Hinton, Alberta
# https://doi.org/10.1016/j.renene.2019.09.141
# Compares simulation results with numerical data from Figure 7 at 1, 5, 10, and 25 years

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
simulation_data_dir() = simulation_dir()
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
# Material properties from Hu et al. Table 1 (at T = 30°C)
materials = HomogenousMaterialProperties{Float_used}(
    2.8811795365701656,      # k_rock - thermal conductivity [W/(m·K)]
    2.1663426199999996e6,    # rho_c_rock - volumetric heat capacity [J/(m³·K)]
    0.6,                     # k_water [W/(m·K)]
    4.17686808e6,            # rho_c_water [J/(m³·K)]
    44.5,                    # k_steel (outer pipe) [W/(m·K)]
    3.728750e6,              # rho_c_steel [J/(m³·K)]
    0.26,                    # k_insulating (inner pipe) [W/(m·K)]
    1.955e6,                 # rho_c_insulating [J/(m³·K)]
    1.0,                     # k_backfill [W/(m·K)] (not used - no backfill region)
    1.0                      # rho_c_backfill [J/(m³·K)] (not used)
)

# =============================================================================
# Borehole geometry
# =============================================================================
# Hinton borehole geometry from Hu et al.
# Deep coaxial borehole (3500m depth) with insulation on inner pipe down to 1000m
borehole = Borehole{Float_used}(
    0.0,                     # xc [m]
    0.0,                     # yc [m]
    3500.0,                  # h - borehole depth [m]
    0.0381,                  # r_inner - inner pipe radius [m]
    0.01,                    # t_inner - inner pipe wall thickness [m]
    0.0889,                  # r_outer - outer pipe inner radius [m]
    0.01,                    # t_outer - outer pipe wall thickness [m]
    0.0989,                  # r_backfill - outer radius (r_outer + t_outer) [m]
    10.0,                    # ṁ - mass flow rate [kg/s]
    1000.0                   # insulation_depth [m]
)

boreholes = (borehole,)

# =============================================================================
# Grid setup
# =============================================================================
# Domain boundaries
xmin, xmax = -100.0, 100.0
ymin, ymax = -100.0, 100.0
zmin, zmax = 0.0, 3700.0

# Grid parameters
dx_fine = 0.0025      # fine spacing near borehole [m]
growth_factor = 1.3   # geometric growth rate
dx_max = 10.0         # maximum spacing far from borehole [m]
dz = 100.0            # vertical spacing [m]

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
# Linear thermal gradient from Hu et al.
# T(z) = T_surface + gradient * z
T0_fresh = initial_condition_thermal_gradient(
    backend, Float_used, gridx, gridy, gridz;
    T_surface=2.29,     # surface temperature [°C]
    gradient=0.035      # thermal gradient [K/m]
)

tspan_full = (0, 3600 * 24 * 365 * 25)  # 25 years [s]
saveat_full = range(tspan_full..., 51)   # save every half year
# saveat_full = [3600 * 24 * 365 * year for year in [0, 1, 5, 10, 25]]

T0, tspan, saveat = prepare_restart(
    T0_fresh, tspan_full, saveat_full;
    checkpoint_dir=checkpoint_dir,
    checkpoint_id=checkpoint_id,
    backend=backend
)

# =============================================================================
# Inlet model
# =============================================================================
# Constant inlet temperature of 20°C
inlet_model = ConstantInlet{Float_used}(20.0)

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
# Simulation completed in 5.054995181868334 hours

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

@save joinpath(simulation_data_dir(), "Hu_et_al_simulation_data.jld2") saved_values Δt cache_cpu t_elapsed

println("Simulation data saved to $(simulation_data_dir())")
