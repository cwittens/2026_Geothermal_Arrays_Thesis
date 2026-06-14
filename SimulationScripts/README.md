# Simulation Scripts

These are the exact setups that produced the simulation data used by the thesis figures. Each
script writes one or more `.jld2` files into the folder that `simulation_dir()` in
[`../paths.jl`](../paths.jl) points to (by default the repo-local
[`../simulation_data/`](../simulation_data/)).

> **You probably do not need to run these.** Regenerating the full data set costs roughly
> **500 GPU-hours on an NVIDIA H200**. To just reproduce the *figures*, use the existing
> simulation data and the scripts in [`../Analysis/`](../Analysis/) instead — see the top-level
> [`README.md`](../README.md).

## Setup

Use the local Julia environment in this folder (Julia `1.12` or newer):

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("reproduce_Brown_et_al.jl")
```

## Backend / hardware

The scripts default to `backend = CUDABackend()` (NVIDIA GPU). The underlying package
[GeothermalWells.jl](https://github.com/cwittens/GeothermalWells.jl) is vendor-agnostic, so you
can switch the `backend = ...` line near the top of each script to:

- `ROCBackend()` — AMD GPUs (requires AMDGPU.jl),
- `CPU()` — CPU execution (correct, but far too slow for the production-sized runs).

The scripts checkpoint to a `checkpoints/` subfolder and can be restarted; checkpoint snapshots
are not part of the archived data.

## What each script produces

All paths below are relative to `simulation_dir()` (default `../simulation_data/`).

### Validation against published studies (Chapter 5)

| Script | Output | Used by figure |
|---|---|---|
| `reproduce_Li_et_al.jl` | `Li_et_al_simulation_data.jld2` | Li radial profile |
| `reproduce_Hu_et_al.jl` | `Hu_et_al_simulation_data.jld2` | Hu in/outlet profiles |
| `reproduce_Brown_et_al.jl` | `Brown_et_al/Brown_et_al_simulation_data.jld2` | Brown single-well profiles |
| `reproduce_Brown_et_al_Fig6_20m_eastwest.jl` | `Brown_et_al/Brown_et_al_line_array_20m_eastwest.jld2` | Brown line array |
| `reproduce_Brown_et_al_Fig6_30m_eastwest.jl` | `Brown_et_al/Brown_et_al_line_array_30m_eastwest.jld2` | Brown line array |
| `reproduce_Brown_et_al_Fig6_40m_eastwest.jl` | `Brown_et_al/Brown_et_al_line_array_40m_eastwest.jld2` | Brown line array |
| `reproduce_Brown_et_al_Fig6_50m_eastwest.jl` | `Brown_et_al/Brown_et_al_line_array_50m_eastwest.jld2` | Brown line array |

### Modeling-choice and parametric studies (Chapter 5)

| Script(s) | Output | Used by figure |
|---|---|---|
| `no_z_diffusion_test.jl` | `no_z_diffusion_test.jld2` | Role of vertical diffusion (compared against `array_studies_2x2/array_study_01.jld2`) |
| `scaling_Q_study_01.jl` … `scaling_Q_study_10.jl` | `scaling_study/scaling_Q_study_0X.jld2` (depths 500 m … 5000 m) | Quadratic reward of depth |
| `array_study_01.jl` … `array_study_13.jl` | `array_studies_2x2/array_study_XX.jld2` (2×2 array; spacings 0, 10, … 120 m) | Array interference (2×2) |
| `3x3_array_study_02.jl` … `3x3_array_study_13.jl` | `array_studies_3x3/3x3_array_study_XX.jld2` (3×3 array; spacings 10, … 120 m) | Array interference (3×3) |

Notes on the array sweeps:

- The index encodes the spacing: index `i` corresponds to a spacing of `(i − 1) × 10 m`, so
  `array_study_01.jl` is spacing 0 (a single isolated well) and `array_study_13.jl` is 120 m.
- `array_study_01.jl` (the single well) doubles as the spacing-0 reference for **both** the 2×2
  and 3×3 sweeps, which is why there is no `3x3_array_study_01.jl`.

### Benchmarks (Chapter 4)

`benchmarking.jl` runs each array configuration (1×1 up to 10×10, plus the row arrays) under
`CUDA.@profile` for a short 10-day simulation and saves the raw profiler output to
`benchmarking_kernel_results.jld2`. The wall-clock and per-kernel timing numbers reported in the
thesis were read off these profiler runs on an NVIDIA H200 and are then **hard-coded** into
[`../Analysis/ch4_benchmarks.jl`](../Analysis/ch4_benchmarks.jl), which produces the
`simulation_time_*` and `kernel_time_breakdown` figures. Because the numbers are hardware
specific, the figures are reproduced from those recorded values rather than by rerunning
`benchmarking.jl`.

### Convergence studies (Chapter 4)

The convergence simulations live inside the corresponding `Analysis/ch4_*_convergence.jl` scripts
rather than here — just a matter of how the studies ended up organized, not because they are cheap
(the finest grid refinements are expensive). The analysis scripts plot from recorded RMS values by
default. See [`../Analysis/README.md`](../Analysis/README.md).
