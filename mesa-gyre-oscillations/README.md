# MESA to GYRE stellar oscillations

This repository turns a MESA pulsation model into a reproducible adiabatic
GYRE calculation. It creates the GYRE input file, runs the solver, records the
settings, and exports a plain-text mode table alongside the native HDF5
summary and eigenfunction files. The clean RGB workflows target GYRE 9.1.x so
they can record the `chi` and `n_iter` root-convergence diagnostics.

The default calculation matches the quadrupolar modes used in the TDE
project:

- harmonic degree: `l = 2`;
- radial-order search: `n_pg = -20 ... 20`;
- frequency scan: dimensionless angular frequency `0.1 ... 20`;
- reported frequency column: cycles per day;
- inverse-frequency grid with 2000 scan points;
- adiabatic, nonrotating oscillations with a vacuum outer boundary.

All of these choices are command-line options. A different star may require a
different frequency range or grid.

## 1. Prepare the stellar model

The input must be a **MESA GYRE-format pulse file**, often given a `.mesa`
extension. A normal MESA `profile*.data` file does not contain the same
columns and cannot be substituted directly.

To export the required file when a MESA run terminates, add the settings in
[`examples/mesa_export_snippet.inlist`](examples/mesa_export_snippet.inlist)
to the corresponding `&star_job` and `&controls` groups of the MESA inlist.
The important setting is

```fortran
pulse_data_format = 'GYRE'
```

The collaborator can put the resulting file in `models/` or pass a path to a
file stored elsewhere. Stellar models and generated results are ignored by
Git by default.

## 2. Install GYRE and the Python packages

Install and compile GYRE by following the
[official GYRE quick start](https://gyre.readthedocs.io/en/stable/user-guide/quick-start.html),
then expose the installation directory:

```bash
export GYRE_DIR=/path/to/gyre
```

Create a small Python environment for the wrapper and notebook:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The pipeline searches for the executable in this order:

1. `--gyre-executable /explicit/path/to/gyre`;
2. `$GYRE_DIR/bin/gyre`;
3. `gyre` on `$PATH`.

## 3. Run the calculation

From the repository root:

```bash
python scripts/run_gyre/pipeline/run_gyre.py models/my_star.mesa \
  --output-dir results/run_gyre/my_star
```

To inspect the generated configuration before launching GYRE:

```bash
python scripts/run_gyre/pipeline/run_gyre.py models/my_star.mesa \
  --output-dir results/run_gyre/my_star \
  --dry-run
```

Use `--help` to see every numerical option. For example, this scans several
harmonic degrees in microhertz:

```bash
python scripts/run_gyre/pipeline/run_gyre.py models/my_star.mesa \
  --l 0 1 2 \
  --freq-min 10 --freq-max 5000 \
  --freq-units UHZ --grid-type LINEAR \
  --n-freq 4000 \
  --output-dir results/run_gyre/my_star_multi_l
```

The pipeline refuses to replace an existing output directory unless
`--overwrite` is supplied.

Equivalent shortcuts are available through the `Makefile`:

```bash
make dry-run MODEL=models/my_star.mesa
make run MODEL=models/my_star.mesa
```

## Outputs

Each run directory contains:

```text
input_model.mesa       copied input model
gyre.in                complete generated GYRE namelist
gyre.log               solver output
summary.h5             one row per mode
detail.id*.l*.n*.h5   radial eigenfunctions with collision-free mode IDs
modes.csv              easy-to-read mode catalogue
run_metadata.json      paths, settings, checksum, and run status
```

The HDF5 outputs follow GYRE's
[documented output format](https://gyre.readthedocs.io/en/stable/ref-guide/output-files/file-formats.html).
The model is copied into the run directory so that `gyre.in` remains portable
and the calculation can be reproduced later.

## Notebook walkthrough

Open the companion notebook from the repository root:

```bash
jupyter lab scripts/run_gyre/notebook/run_gyre.ipynb
```

The notebook validates a model, displays the generated namelist, optionally
runs GYRE, loads the mode catalogue, and plots the radial and horizontal
eigenfunctions of a selected mode. Its `RUN_GYRE` switch is `False` by
default, so opening the notebook does not launch a calculation unexpectedly.

### Clean RGB notebooks

The clean notebooks contain no saved results or personal absolute paths. Each
one imports the shared case definitions in `pipeline/rgb_cases.py`, performs a
read-only preflight, and keeps the real GYRE run opt-in:

| Case | Notebook | Expected pulse model in `models/` |
|---|---|---|
| 1 M_sun, log(R/R_sun)=0.5 | `GYRE_RGB1Msun_Z0p0142_3p16R_clean.ipynb` | `RGB1Msun_Z0p0142_logR0p5.mesa` |
| 1 M_sun, 10 R_sun | `GYRE_RGB1Msun_Z0p0142_10R_clean.ipynb` | `RGB1Msun_Z0p0142_logR1p0.mesa` |
| 1 M_sun, about 20 R_sun | `GYRE_RGB1Msun_Z0p0142_20R_clean.ipynb` | `RGB1Msun_Z0p0142_logR1p3.mesa` |
| 1 M_sun, about 40 R_sun | `GYRE_RGB1Msun_Z0p0142_40R_clean.ipynb` | `RGB1Msun_Z0p0142_logR1p6.mesa` |
| 3 M_sun, 10 R_sun | `GYRE_RGB3Msun_Z0p0142_10R_clean.ipynb` | `RGB3Msun_Z0p0142_logR1p0.mesa` |

The legacy file called `3.61R` is mislabeled: `log(R/R_sun)=0.5` and its
saved `R_star` both correspond to about 3.16 R_sun. The legacy 3-M_sun
notebook also contains stale saved arrays from a 1-M_sun model. The clean
notebooks preserve neither set of embedded results and verify the mass and
radius from each new GYRE summary.

The required `.mesa` pulse files are deliberately not committed. Obtain them
from the collaborator or export them from MESA, place them under `models/`,
and then open a clean notebook, for example:

```bash
export GYRE_DIR=/path/to/gyre-9.1.1
jupyter lab scripts/run_gyre/notebook/GYRE_RGB1Msun_Z0p0142_10R_clean.ipynb
```

The starting scan bands are case-specific estimates. Production work still
requires increasing `n_freq`, `w_osc`, and `w_ctr`, widening the frequency
and radial-order bounds, and checking that the selected frequencies remain
stable.

## Repository layout

```text
examples/mesa_export_snippet.inlist
models/                              collaborator inputs (ignored)
results/run_gyre/                    generated runs (ignored)
scripts/run_gyre/
  pipeline/gyre_pipeline.py         reusable GYRE runner
  pipeline/rgb_cases.py             shared RGB cases, diagnostics, and plots
  pipeline/build_rgb_notebooks.py   deterministic notebook generator
  notebook/run_gyre.ipynb            guided interactive workflow
  notebook/*_clean.ipynb             clean RGB workflows
```

## Common issues

- **The model is rejected:** confirm that MESA wrote pulse data with
  `pulse_data_format = 'GYRE'`; renaming a standard profile does not convert
  it.
- **No modes are found:** widen the frequency interval, increase `--n-freq`,
  or change between `LINEAR` and `INVERSE` scanning.
- **A run uses too many cores:** pass `--threads 1` (the default) or another
  explicit OpenMP thread count.
- **The frequency range is unclear:** inspect the star's dynamical frequency
  and use a broad preliminary scan before increasing the grid resolution.
- **Rotation or nonadiabatic physics is required:** the current wrapper is an
  adiabatic, nonrotating baseline. Extend the generated `&rot` or
  `&nad_output` configuration only after choosing the intended physics.

The clean RGB workflow is prepared for GYRE 9.1.x and Python 3.9. Configuration
generation, dry runs, diagnostics, and notebook execution can be tested
without launching GYRE; a physical oscillation run still requires a compiled
GYRE executable and the matching MESA pulse model.
