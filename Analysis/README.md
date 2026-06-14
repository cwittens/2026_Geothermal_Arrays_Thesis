# Analysis Scripts

These scripts read the simulation data in [`../simulation_data/`](../simulation_data/) and
regenerate the thesis figures into [`../figures/`](../figures/). They are cheap: each runs in
minutes on a laptop. The output location is controlled by `plots_dir()` in [`../paths.jl`](../paths.jl).

## Setup

Run everything from this folder, using the local Julia environment (Julia `1.12` or newer):

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()          # one-time: install the exact dependency versions
include("ch5_literature_brown.jl")
```

Each `chX_*.jl` file is self-contained and can be run independently. To regenerate every figure:

```julia
for f in filter(endswith(".jl"), readdir())
    startswith(f, "ch") && include(f)
end
```

## Scripts and the figures they produce

Files are named `ch<chapter>_<topic>.jl` so they sort in the order of the thesis.

| Script | Thesis section | Figures (written to `../figures/`) |
|---|---|---|
| `ch3_grid_overviews.jl` | §3.2 Spatial Discretization | `Grid_overview_single.pdf`, `Grid_overview_2x2.pdf`, `Grid_overview_3x3.pdf` |
| `ch3_rock2_stability.jl` | §3.4 Vertical Diffusion (ROCK2) | `rock2_stability_region_s_5.pdf`, `rock2_stability_region_s_9.pdf` |
| `ch4_benchmarks.jl` | §4.4–4.5 Simulation & Kernel Benchmarks | `simulation_time_vs_number_of_wells_25.pdf`, `simulation_time_vs_number_of_wells_100.pdf`, `kernel_time_breakdown.pdf` |
| `ch4_advection_convergence.jl` | §4.6 Convergence — Advection | `advection_spatial_convergence.pdf` |
| `ch4_diffusion_convergence.jl` | §4.6 Convergence — Diffusion | `diffusion_spatial_convergence.pdf`, `diffusion_temporal_convergence.pdf` |
| `ch4_manufactured_solution_convergence.jl` | §4.6 Convergence — Manufactured Solution | `manufactured_solution_convergence.pdf` |
| `ch5_literature_li.jl` | §5.1 Li et al. validation | `Li_et_al_radial_profile.pdf` |
| `ch5_literature_hu.jl` | §5.2 Hu et al. validation | `Hu_et_al_in_outlet.pdf` |
| `ch5_literature_brown.jl` | §5.3 Brown et al. validation (single well + line array) | `Brown_et_al_in_outlet.pdf`, `Brown_et_al_radial_profile.pdf`, `Brown_et_al_array_radial_profile.pdf` |
| `ch5_vertical_diffusion_role.jl` | §5.4 The Role of Vertical Diffusion | `no_z_diffusion_delta_T.pdf` |
| `ch5_depth_scaling.jl` | §5.5 The Quadratic Reward of Depth | `scaling_Q_with_fits.pdf` |
| `ch5_array_interference.jl` | §5.6 How Close Is Too Close? | `array_2x2_static.pdf`, `interference_sqrt_t.pdf`, `array_3x3_classes.pdf` |

`theme.jl` defines the shared color palette and plot defaults; it is included by every script.

## Which scripts need simulation data?

Most scripts load `.jld2` files from `../simulation_data/` (see
[`../simulation_data/README.md`](../simulation_data/README.md) for the data → figure mapping).

The exceptions are the **convergence studies** (`ch4_*`) and the **grid / stability** plots
(`ch3_*`), which need no external data: they either compute their (cheap) simulation in-script or
plot purely analytical quantities.

## The convergence studies (Chapter 4)

Unlike the array and validation runs, the convergence studies keep their simulation code inside
the analysis script itself. This is just how they ended up organized — these runs are **not**
cheap (the finest grid refinements in particular are expensive). By default the scripts therefore
plot from RMS error values that are recorded in the file, so they reproduce the figures without
rerunning anything.

To actually rerun a convergence simulation, uncomment the line marked
`# to actually run the simulation`:

- `ch4_advection_convergence.jl` — the `RMS_space = [advection_simulation(...) ...]` line.
- `ch4_diffusion_convergence.jl` — the `RMS_space = ...` and `RMS_time = ...` lines.
- `ch4_manufactured_solution_convergence.jl` — the `RMS = ...` line (and `RMS_no_advection = ...`
  for the second-order, advection-off check).

The backend for these is set by `backend = CPU()` near the bottom of each script; switch it to
`CUDABackend()` or `ROCBackend()` if you have a GPU and the matching package installed.
