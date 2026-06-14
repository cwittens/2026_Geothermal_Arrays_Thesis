# Reproducibility Repository — Three-Dimensional Modeling of Deep Borehole Heat Exchanger Arrays

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20597689.svg)](https://doi.org/10.5281/zenodo.20597689)

This repository contains the code and data needed to reproduce every numerical figure in the
M.Sc. thesis

> **Three-Dimensional Modeling of Deep Borehole Heat Exchanger Arrays** <br>
> Collin Wittenstein, Johannes Gutenberg University Mainz, 2026.

The numerical method itself lives in a separate, registered open-source Julia package,
[**GeothermalWells.jl**](https://github.com/cwittens/GeothermalWells.jl). This repository builds on
**GeothermalWells.jl v0.3.1**, and only contains the *experiment* scripts (the exact simulation
setups used in the thesis) and the *analysis* scripts that turn the resulting data into the thesis
figures. The pinned version is recorded in each environment's `Manifest.toml`, so
`Pkg.instantiate()` installs exactly v0.3.1.

## Abstract

Modeling deep borehole heat exchanger (DBHE) arrays in high-fidelity three dimensions is a genuinely hard computational problem: the borehole geometry is only tens of centimeters across, yet the surrounding rock must be resolved over hundreds of meters laterally and several kilometers vertically, while the simulation must span decades of operation. For this reason, fully discretized three-dimensional simulation at the array scale has long been considered prohibitively expensive. This thesis argues the opposite. With numerical methods tailored to the physics and geometry of the problem, and an implementation built to exploit modern parallel hardware, such simulations become practical on a single GPU.

The physical model reduces the system to a single advection-diffusion equation for the temperature field, with rock, grout, pipe walls, and working fluid entering through spatially varying material properties. A wide separation of time scales between advection, horizontal diffusion, and vertical diffusion motivates an operator splitting into three subproblems, each matched to a suitable method: a stabilized explicit Runge–Kutta scheme for the slow vertical diffusion, an Alternating Direction Implicit method for the stiff horizontal diffusion, and an unconditionally stable Semi-Lagrangian method for the advection. The implementation combines custom GPU kernels with Julia's existing solver ecosystem.

The method is validated against three published studies that disagree with one another on how a DBHE should be modeled, ranging from a reduced-dimensionality cylindrical model to a full finite-element model with turbulent flow, with agreement good to excellent. As a first proof of concept, the method is then used for the kind of design study that array deployment requires: how far apart wells must stand before they thermally interfere. It is only one such study in a large parameter space waiting to be optimized, a task left to future work.

The software developed in this work is released as an open-source Julia package, [GeothermalWells.jl](https://github.com/cwittens/GeothermalWells.jl).

## How this repository is organized

The expensive part (running the simulations) is **separated** from the cheap part (making the
figures). There are two top-level Julia projects, each with its own environment:

```text
.
├── SimulationScripts/    # the exact setups that produced the simulation data (GPU, expensive)
├── simulation_data/      # the .jld2 outputs the analysis reads (download these, or regenerate)
├── Analysis/             # scripts that read simulation_data/ and write the thesis figures
├── figures/              # the generated figures (final thesis figures are committed here)
├── paths.jl              # central configuration: where data is read from / figures written to
├── LICENSE
└── README.md             # you are here
```

- [`SimulationScripts/`](SimulationScripts/) — one script per simulation. These need a GPU to be
  practical and are collectively expensive (**~500 GPU-hours on an NVIDIA H200**). The backend is
  configurable; see [`SimulationScripts/README.md`](SimulationScripts/README.md) for the available
  backends and per-experiment costs.
- [`Analysis/`](Analysis/) — one script per thesis figure (or small group of figures). These run
  in minutes on a laptop CPU and only need the data in `simulation_data/`. See
  [`Analysis/README.md`](Analysis/README.md).
- [`simulation_data/`](simulation_data/) — the `.jld2` data products. See
  [`simulation_data/README.md`](simulation_data/README.md).

## Reproducing the figures

You almost certainly want **Path A**.

### Path A — regenerate the figures from existing data (recommended, minutes)

Rerunning all simulations takes roughly 500 GPU-hours, so for reproducing the *figures* it makes
much more sense to use the already-computed simulation data.

1. **Install [Julia](https://julialang.org/) `1.12` or newer.** The data was produced with Julia
   `1.12`.

2. **Get the simulation data.** The full data set is archived on Zenodo at
   [doi.org/10.5281/zenodo.20689818](https://doi.org/10.5281/zenodo.20689818). Download it and
   either unpack it into the repo-local `simulation_data/` folder or point `simulation_dir()` in
   [`paths.jl`](paths.jl) at wherever you unpacked it. Mind the
   [expected folder layout](simulation_data/README.md#expected-folder-layout) — the array and
   scaling runs must sit in their own subfolders.

3. **Instantiate the analysis environment and run a figure script.** From the `Analysis/` folder:

   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()          # one-time: install the exact dependency versions

   include("ch5_literature_brown.jl")   # e.g. the Brown et al. validation figures
   ```

   Each script is self-contained and writes its figures to `figures/` (configurable in
   `paths.jl`). To regenerate **all** figures at once, from the `Analysis/` folder:

   ```julia
   for f in filter(endswith(".jl"), readdir())
       startswith(f, "ch") && include(f)   # run every chX_*.jl script
   end
   ```

### Path B — rerun the simulations from scratch (expensive, GPU)

Only needed if you want to regenerate the underlying data. See
[`SimulationScripts/README.md`](SimulationScripts/README.md) for details and per-experiment
costs. In short, from the `SimulationScripts/` folder:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("reproduce_Brown_et_al.jl")   # writes Brown_et_al_simulation_data.jld2 to simulation_data/
```

The default backend in these scripts is `CUDABackend()` (NVIDIA) and for this it was tested. The underlying package also supports `ROCBackend()` (AMD) and `CPU()`; the convergence studies in `Analysis/` default to `CPU()`.

## Figure → thesis mapping

Analysis scripts are named `ch<chapter>_<topic>.jl` so they sort in thesis order. Each produces
the figures listed below (all written to `figures/`).

| Analysis script | Thesis section | Figure file(s) |
|---|---|---|
| [`ch3_grid_overviews.jl`](Analysis/ch3_grid_overviews.jl) | §3.2 Spatial Discretization | `Grid_overview_single.pdf`, `Grid_overview_2x2.pdf`, `Grid_overview_3x3.pdf` |
| [`ch3_rock2_stability.jl`](Analysis/ch3_rock2_stability.jl) | §3.4 Vertical Diffusion (ROCK2) | `rock2_stability_region_s_5.pdf`, `rock2_stability_region_s_9.pdf` |
| [`ch4_benchmarks.jl`](Analysis/ch4_benchmarks.jl) | §4.4–4.5 Simulation & Kernel Benchmarks | `simulation_time_vs_number_of_wells_25.pdf`, `simulation_time_vs_number_of_wells_100.pdf`, `kernel_time_breakdown.pdf` |
| [`ch4_advection_convergence.jl`](Analysis/ch4_advection_convergence.jl) | §4.6 Convergence — Advection | `advection_spatial_convergence.pdf` |
| [`ch4_diffusion_convergence.jl`](Analysis/ch4_diffusion_convergence.jl) | §4.6 Convergence — Diffusion | `diffusion_spatial_convergence.pdf`, `diffusion_temporal_convergence.pdf` |
| [`ch4_manufactured_solution_convergence.jl`](Analysis/ch4_manufactured_solution_convergence.jl) | §4.6 Convergence — Manufactured Solution | `manufactured_solution_convergence.pdf` |
| [`ch5_literature_li.jl`](Analysis/ch5_literature_li.jl) | §5.1 Li et al. validation | `Li_et_al_radial_profile.pdf` |
| [`ch5_literature_hu.jl`](Analysis/ch5_literature_hu.jl) | §5.2 Hu et al. validation | `Hu_et_al_in_outlet.pdf` |
| [`ch5_literature_brown.jl`](Analysis/ch5_literature_brown.jl) | §5.3 Brown et al. validation (single well + line array) | `Brown_et_al_in_outlet.pdf`, `Brown_et_al_radial_profile.pdf`, `Brown_et_al_array_radial_profile.pdf` |
| [`ch5_vertical_diffusion_role.jl`](Analysis/ch5_vertical_diffusion_role.jl) | §5.4 The Role of Vertical Diffusion | `no_z_diffusion_delta_T.pdf` |
| [`ch5_depth_scaling.jl`](Analysis/ch5_depth_scaling.jl) | §5.5 The Quadratic Reward of Depth | `scaling_Q_with_fits.pdf` |
| [`ch5_array_interference.jl`](Analysis/ch5_array_interference.jl) | §5.6 How Close Is Too Close? | `array_2x2_static.pdf`, `interference_sqrt_t.pdf`, `array_3x3_classes.pdf` |

## A note on the convergence studies (Chapter 4)

The three convergence studies are organized differently from the rest: the simulation code lives
*inside* the analysis script rather than in `SimulationScripts/`. This is purely a matter of how
the studies ended up structured — these runs are **not** cheap, the finest grid refinements in
particular are expensive. Because of that, by default these scripts plot from the RMS error values
recorded in the file, so they reproduce the figures without rerunning anything.

To actually rerun the convergence simulation, uncomment the line marked
`# to actually run the simulation` in the corresponding script:

- [`ch4_advection_convergence.jl`](Analysis/ch4_advection_convergence.jl) — the
  `RMS_space = [advection_simulation(...) ...]` line.
- [`ch4_diffusion_convergence.jl`](Analysis/ch4_diffusion_convergence.jl) — the `RMS_space = ...`
  and `RMS_time = ...` lines.
- [`ch4_manufactured_solution_convergence.jl`](Analysis/ch4_manufactured_solution_convergence.jl)
  — the `RMS = ...` (and, for the second-order check, `RMS_no_advection = ...`) lines.

## Approximate compute cost (Path B)

All production simulations were run on a single **NVIDIA H200** GPU. Regenerating the full data
set costs roughly **500 GPU-hours** in total; the breakdown is approximately:

| Experiment | Simulated time | Approx. wall-clock (H200) |
|---|---|---|
| Li et al. (single well) | 120 days | < 5 min |
| Hu et al. (single well) | 25 years | ~5 h |
| Brown et al. (single well) | 20 years | ~3.7 h |
| Brown et al. line arrays (4 spacings) | 20 years each | ~50 h total |
| Depth scaling (10 depths) | 20 years each | ~45 h total |
| 2×2 and 3×3 interference sweeps (12 spacings each) | 20 years each | ~400 h total |

## Citation

If you use this repository, please cite the thesis and this repository.

```bibtex
@mastersthesis{wittenstein2026dbhearrays,
  title  = {Three-Dimensional Modeling of Deep Borehole Heat Exchanger Arrays},
  author = {Wittenstein, Collin},
  school = {Johannes Gutenberg University Mainz},
  year   = {2026},
  type   = {M.Sc. Thesis}
}

@misc{wittenstein2026dbhearraysRepro,
  title        = {Reproducibility repository for
                  "Three-Dimensional Modeling of Deep Borehole Heat Exchanger Arrays"},
  author       = {Wittenstein, Collin},
  year         = {2026},
  howpublished = {\url{https://github.com/cwittens/2026_Geothermal_Arrays_Thesis}},
  doi          = {10.5281/zenodo.20597689}
}
```

Please also cite the underlying software package
[GeothermalWells.jl](https://github.com/cwittens/GeothermalWells.jl).

## License

The code in this repository is released under the MIT license; see [`LICENSE`](LICENSE).

## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
