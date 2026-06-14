# paths.jl — central place to configure where data is read from and figures are written to.
#
# This file is `include`d by every script in `Analysis/` and `SimulationScripts/`,
# so changing a path here changes it everywhere.

const REPO_ROOT = @__DIR__

# ---------------------------------------------------------------------------
# Where figures are written.
# By default they go into the repo-local `figures/` folder. To regenerate the
# figures directly into the thesis LaTeX project instead, point this at that
# folder (example commented out below).
# ---------------------------------------------------------------------------
plots_dir() = joinpath(REPO_ROOT, "figures")

# ---------------------------------------------------------------------------
# Where the simulation data lives.
#
#  * The `Analysis/` scripts READ the `.jld2` files from here.
#  * The `SimulationScripts/` scripts WRITE their `.jld2` files to here.
#
# By default this is the repo-local `simulation_data/` folder. The full data set
# is large (~500 GPU-hours to regenerate); if you download the archived data
# instead of rerunning the simulations, either drop it into `simulation_data/`
# or point this function at wherever you unpacked it.
#
# Archived simulation data: https://doi.org/10.5281/zenodo.20689818
# ---------------------------------------------------------------------------
simulation_dir() = joinpath(REPO_ROOT, "simulation_data")

function ensure_plots_dir()
    mkpath(plots_dir())
    return plots_dir()
end
