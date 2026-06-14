# Simulation Data

This folder holds the `.jld2` simulation outputs that the analysis scripts read to produce the
thesis figures. Each file stores the saved temperature-field snapshots (`saved_values`), the time
step (`Δt`), a CPU-side cache with the grids and borehole geometry (`cache_cpu`), and a recorded
runtime (`t_elapsed`, which understates the true total for runs that were restarted from a
checkpoint).

The [`Analysis/`](../Analysis/) scripts read from here; the [`SimulationScripts/`](../SimulationScripts/)
scripts write to here. The location is configurable via `simulation_dir()` in
[`../paths.jl`](../paths.jl).

> **Download:** the full data set is archived on Zenodo at
> [doi.org/10.5281/zenodo.20689818](https://doi.org/10.5281/zenodo.20689818). Download it and
> either place the files here or point `simulation_dir()` in `paths.jl` at the folder where you
> unpacked it (see the layout below). Checkpoint/restart snapshots are intentionally **not** part
> of the archive — only the final data products below are needed to reproduce the figures.

## Expected folder layout

The Zenodo archive contains three loose `.jld2` files and four `.zip` archives. Unpack it so that
**the three loose files sit at the top level and each `.zip` becomes a subfolder of the same name**
(the analysis scripts read the array, scaling, and Brown runs from those subfolders). After
unpacking, the folder `simulation_dir()` points to should look like this:

```text
simulation_data/
├── Li_et_al_simulation_data.jld2          # loose file on Zenodo
├── Hu_et_al_simulation_data.jld2          # loose file on Zenodo
├── no_z_diffusion_test.jld2               # loose file on Zenodo
├── Brown_et_al/                           # from Brown_et_al.zip
│   ├── Brown_et_al_simulation_data.jld2
│   ├── Brown_et_al_line_array_20m_eastwest.jld2
│   ├── Brown_et_al_line_array_30m_eastwest.jld2
│   ├── Brown_et_al_line_array_40m_eastwest.jld2
│   └── Brown_et_al_line_array_50m_eastwest.jld2
├── array_studies_2x2/                     # from array_studies_2x2.zip: array_study_01 … 13
├── array_studies_3x3/                     # from array_studies_3x3.zip: 3x3_array_study_02 … 13
└── scaling_study/                         # from scaling_study.zip: scaling_Q_study_01 … 10
```

## Contents

| File / folder | Produced by | Read by |
|---|---|---|
| `Li_et_al_simulation_data.jld2` | `reproduce_Li_et_al.jl` | `ch5_literature_li.jl` |
| `Hu_et_al_simulation_data.jld2` | `reproduce_Hu_et_al.jl` | `ch5_literature_hu.jl` |
| `Brown_et_al/Brown_et_al_simulation_data.jld2` | `reproduce_Brown_et_al.jl` | `ch5_literature_brown.jl` |
| `Brown_et_al/Brown_et_al_line_array_{20,30,40,50}m_eastwest.jld2` | `reproduce_Brown_et_al_Fig6_*_eastwest.jl` | `ch5_literature_brown.jl` |
| `no_z_diffusion_test.jld2` | `no_z_diffusion_test.jl` | `ch5_vertical_diffusion_role.jl` |
| `scaling_study/scaling_Q_study_01…10.jld2` | `scaling_Q_study_01…10.jl` | `ch5_depth_scaling.jl` |
| `array_studies_2x2/array_study_01…13.jld2` | `array_study_01…13.jl` | `ch5_array_interference.jl`, `ch5_vertical_diffusion_role.jl` |
| `array_studies_3x3/3x3_array_study_02…13.jld2` | `3x3_array_study_02…13.jl` | `ch5_array_interference.jl` |

See [`../SimulationScripts/README.md`](../SimulationScripts/README.md) for what each run contains
(depths, spacings, indexing convention).
