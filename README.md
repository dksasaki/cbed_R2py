# cbed_r2py

A Python wrapper that parallelizes the [CBED](link) diagenetic model across a 2D ocean grid, driven by MOM6-COBALT model output.

## Table of contents

1. [Installation](#installation)
2. [Configuration](#configuration)
3. [Input files](#input-files)
4. [Workflow overview](#workflow-overview)
5. [Scripts](#scripts)

## 1. Installation

Requires [pixi](https://pixi.prefix.dev/latest/installation/). All dependencies are managed via `pixi.toml`.

```bash
git clone https://github.com/dksasaki/cbed_R2py.git
cd cbed_R2py
pixi install
pixi run setup-r
```

`setup-r` installs the R packages `ReacTran` and `marelac` through CRAN. All packages are installed within the pixi environment. Dependencies are pinned in `pixi.toml`.

The `newton_krylov` method (see [CBED method](#cbed-method)) uses the
[`cbed2py`](https://github.com/dksasaki/cbedr2py_1d) package, installed by pixi
from a clone next to this repository (`../cbedr2py_1d`, see `[pypi-dependencies]`
in `pixi.toml`):

```bash
git clone https://github.com/dksasaki/cbedr2py_1d.git ../cbedr2py_1d
git submodule update --init      # src/cbed_R, used by the steady_state method
```



## 2. Configuration

Edit `config.toml` at the project root before running any script. `ncdump` headers for all expected input files are available in `aux/` for reference. **`root_dir` must be the absolute path to the cloned repository** — all R scripts and internal modules are resolved relative to it. All other paths are resolved relative to `config.toml` and can be relative or absolute.

```toml
[paths]
root_dir  = "/path/to/cbed_r2py"      # absolute path to the cloned repository
fpath     = "/path/to/model/output"   # must contain cobalt_btm.nc, cobalt_tr.nc, mom6.nc
ftopo     = "/path/to/ocean_topog.nc"
fgrd      = "/path/to/porosity.grd"
ftemplate = "/path/to/template.nc"
fout      = "data/cache/porosity_out.nc"
cache_dir = "data/cache"

[chunks]
nx = 8
ny = 1

[stitch]
mom6_template = "data/cache/scratch_test/mom6.nc"
cbed_pattern  = "data/cache/cbed*.nc"
output        = "data/cbed_results.nc"

[cbed]
method = "steady_state"   # or "newton_krylov"
model  = "fast"           # cbed2py model for newton_krylov
```

### CBED method

`[cbed] method` selects how `cbed_wrapper.py` computes the sediment state at each point:

| `method` | Model | Forcing | Result |
|----------|-------|---------|--------|
| `steady_state` (default, original behavior) | CBED v1, `src/cbed_R/cbed_v1_func.R` | time mean of the COBALT fields | steady state |
| `newton_krylov` | cbed2py transient model (`model = "fast"` or `"v2"`) | one-year climatology (by day of year, or by month for monthly output) of the COBALT bottom fields, repeated every year | periodic state at Jan 1, found by Jacobian-free Newton-Krylov (`cbed2py.spinup.newton_krylov_spinup`) |

With `newton_krylov`:

- J.OM, O2.w, NO3.w, DIC.w (from `cobalt_btm`), and NH4.w and TAlk.w (from `cobalt_tr`) vary with the seasonal cycle. Temperature and salinity stay at their mean, because `model_reader.py` caches `mom6.nc` already averaged over time.
- `[cbed.params]` overrides `cbed2py` parameters such as `k.NOx`, `k.ODUox` and the `jpoc_mode_*` options. `[cbed.spinup]` passes keywords to `newton_krylov_spinup`, such as `max_newton`, `gmres_iter`, `tol_imbalance`, `tol_f` and `n_sub`.
- Points that fail stay NaN, as with `steady_state`. Points that don't converge keep the best state found, and the per-point diagnostics flag them.

Cost per grid point, measured on one test point with a seasonal J.OM and temperature:

| | Time | Peak memory per worker |
|---|---|---|
| `steady_state` (v1) | 1.8 s | 180 MB |
| `newton_krylov`, `fast` (~50 model years) | 14 s (~8×) | 240 MB |
| `newton_krylov`, `v2` (~90 model years) | 95 s (~55×) | 215 MB |

Raise `--time` in `execute.sh` to match. Memory per worker barely changes.



## 3. Input files

| Variable | File |
|----------|------|
| Bottom temperature, salinity | `mom6.nc` |
| Bottom O2, NO3, DIC, particle fluxes | `cobalt_btm.nc` |
| Bottom NH4, TAlk | `cobalt_tr.nc` |
| Porosity | computed via `porosity_main()` |
| Topography | `ftopo` (set in `config.toml`) |

Reference `ncdump` headers are in `aux/`.


## 4. Workflow overview

```
  config.toml
      │
      ├───────────────────────────────────────┐
      ▼                                       ▼
cbed_wrapper.py  ×  n_chunks            cbed_stitch.py
─────────────────────────────       ──────────────────────────
    MOM6 + COBALT input                 cbed_mom6_001.nc
            │                           cbed_mom6_002.nc
            ▼                               ...
      CBED model (R/rpy2)               cbed_mom6_NNN.nc
      multiprocessing.Pool                     │
            │                                  ▼
            ▼                           cbed_results.nc
  cbed_mom6_<chunk:03d>.nc
```

The domain is split into `nx × ny` tiles (default `8 × 1`). Run `cbed_wrapper.py` once per chunk, then `cbed_stitch.py` once to assemble the final output. For details to run these scripts read below.


## 5. Scripts

### `cbed_wrapper.py`

Runs the CBED model for one domain chunk, distributing grid points across multiple processes. Chunking parameters and the CBED method are read from `config.toml`. Output is saved to `data/cache/cbed_mom6_<chunk:03d>.nc`.

Per-point diagnostics go to `data/cache/spinup_diag_<chunk:03d>.nc`: `status` (1 ok/converged, 0 not converged, -1 failed, -2 land), `wall_time` [s] and, for `newton_krylov`, `model_years`, `newton_steps`, `imbalance` and `rms_G`. 

```bash
python scripts/cbed_wrapper.py <nproc> <chunk>
```

| Argument | Description |
|----------|-------------|
| `nproc` | Number of parallel worker processes |
| `chunk` | Tile index to process (0 = full domain, >0 selects the chunk) |



Example
```bash
python cbed_wrapper.py 200 2  # uses 200 processors and selects the second chunk
```


**Running on a cluster (SLURM)**

A SLURM batch script example is provided at `execute.sh`. Edit the `#SBATCH` directives and the last line to match your chunk index and available resources, then submit with:

```bash
sbatch execute.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=cbed
#SBATCH --partition=sharing
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=120G
#SBATCH --time=01:00:00
#SBATCH --exclude=d3032,d3232,d3203
pixi run python scripts/cbed_wrapper.py 192 5
```

`--cpus-per-task` should match the `nproc` argument passed to `cbed_wrapper.py`.

### `cbed_stitch.py`

Assembles all chunk files into a single analysis-ready NetCDF. All parameters come from `config.toml` — no arguments required. Raises an error if the number of `cbed_mom6_*.nc` files found is not exactly `cbed_mom6_001.nc` … `cbed_mom6_<nx×ny>.nc` (three-digit chunk index, e.g. `cbed_mom6_040.nc`). Output is saved to `data/cbed_results.nc` with dimensions `(l=20, y=ny, x=nx)`.

```bash
python scripts/cbed_stitch.py
```

