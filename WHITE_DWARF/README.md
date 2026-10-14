# MESA MAMS Profile Tutorial

Date: 2026-07-16
Source: the MESA work directory in `tutorial/` and the automated profile
pipeline in `scripts/mams_profiles/`

This repository teaches a collaborator how to generate middle-age
main-sequence (MAMS) stellar models with MESA. It produces two useful files
for every requested stellar mass:

- a regular MESA structure profile named `MAMS<mass>Msun.profile`;
- a GYRE-format pulsation model named `MAMS<mass>Msun.mesa`.

The notebook creates an independent work directory for each mass. This layout
keeps configurations, logs, profiles, and scheduler jobs separated and makes a
mass grid safe to reproduce.

## Required background

### What MESA does

[Modules for Experiments in Stellar Astrophysics (MESA)](https://docs.mesastar.org/)
solves the one-dimensional stellar-structure and composition-evolution
equations. A MESA *work directory* contains:

- **inlists**, which specify initial conditions, microphysics, numerical
  controls, output controls, and termination conditions;
- **Fortran extras**, which can add custom behavior;
- a locally compiled **`star` executable**, linked to one MESA installation;
- runtime products such as histories, profiles, photos, and terminal logs.

The compiled `star` executable depends on the operating system, compiler,
MESA SDK, and MESA source tree. Each collaborator therefore compiles `star`
on the machine that will run the models.

### What MAMS means in this project

The initial hydrogen mass fraction is determined by the adopted helium and
metal fractions:

$$
X_{\mathrm{initial}}
= 1-Y-Z
= 1-0.2703-0.0142
= 0.7155.
$$

Here, $X_{\mathrm{initial}}$ is the initial hydrogen mass fraction, $Y$ is the
helium mass fraction, and $Z$ is the metal mass fraction. The model evolves
from the pre-main sequence until the central hydrogen abundance satisfies

$$
X_{\mathrm{c}}(\mathrm{H}^1) < 0.35,
$$

where $X_{\mathrm{c}}(\mathrm{H}^1)$ is the central mass fraction of hydrogen-1.
The threshold is close to one half of the initial hydrogen abundance, so this
project uses it as its operational definition of a middle-age main-sequence
model. MESA checks the condition after a timestep, which means the final value
can lie slightly below 0.35.

## Software requirements

The collaborator needs:

1. A Linux or macOS system supported by MESA.
2. A MESA SDK compatible with the chosen MESA release.
3. A local MESA installation.
4. Python 3.9 or newer.
5. JupyterLab for the guided notebook.
6. Slurm only when submitting the prepared models to a Slurm cluster.
7. GYRE only when computing oscillation modes from the generated `.mesa`
   pulsation files. MESA itself writes the GYRE-format files.

The original profiles were generated with the following recorded build:

| Component | Recorded value |
|---|---|
| MESA version identifier | `6d0b0f2` |
| Compiler | `gfortran 13.1.0` |
| MESA SDK | `x86_64-linux-23.7.3` |
| Math backend | `CRMATH` |

Using the same MESA version gives the closest numerical reproduction. A newer
MESA release may rename controls or change numerical details. In that case,
compare the local defaults and update `tutorial/inlist_project` deliberately.

The official installation instructions are at
<https://docs.mesastar.org/en/latest/installation.html>.

## Repository layout

```text
mesa_tutorial/
├── README.md
├── requirements.txt
├── tutorial/
│   ├── inlist
│   ├── inlist_project
│   ├── inlist_pgstar
│   ├── mk, rn, re, clean
│   ├── make/makefile
│   └── src/
│       ├── run.f90
│       └── run_star_extras.f90
└── scripts/mams_profiles/
    ├── pipeline/
    │   ├── __init__.py
    │   └── mams_pipeline.py
    └── notebook/
        └── mams_profiles.ipynb
```

Generated models appear under `results/mams_profiles/`. Git ignores that
directory because the results are reproducible products rather than source
files.

## First-time installation

### 1. Clone this repository

```bash
git clone git@github.com:ZihanZhou26/mesa_tutorial.git
cd mesa_tutorial
```

The repository is private. The owner must grant the collaborator access in the
GitHub repository settings before cloning.

### 2. Install MESA and the MESA SDK

Follow the official MESA installation guide for the collaborator's operating
system. After installation, start a fresh shell and define `MESA_DIR` as the
absolute path to the MESA source tree:

```bash
export MESA_DIR=/absolute/path/to/mesa
```

Many MESA installations provide an SDK initialization command as well. Run
the command specified by the selected SDK before compiling the work directory.
The exact command depends on the SDK release and shell.

Confirm the environment:

```bash
test -d "$MESA_DIR/star"
echo "$MESA_DIR"
```

### 3. Compile the local `star` executable

```bash
cd tutorial
./mk
cd ..
```

Successful compilation creates `tutorial/star`. The Python pipeline links each
isolated run directory to this executable. Compile again after changing files
in `tutorial/src/` or after switching the MESA installation. Changes to inlists
alone take effect on the next run with the existing executable.

Confirm the build:

```bash
test -x tutorial/star
```

### 4. Create a Python environment

Using a virtual environment keeps the notebook dependencies local:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### 5. Start the notebook

From the repository root:

```bash
jupyter lab scripts/mams_profiles/notebook/mams_profiles.ipynb
```

Run the notebook cells from top to bottom. The notebook locates the repository
root automatically when Jupyter starts from the repository or notebook
directory.

## Physics in the supplied model

The authoritative physics template is `tutorial/inlist_project`.

| Setting | Supplied value | Meaning |
|---|---:|---|
| Initial model | pre-main sequence | MESA constructs and relaxes a pre-main-sequence model |
| Nuclear network | `mesa_49.net` | 49-isotope reaction network |
| Initial mass | set by notebook | Stellar mass in solar masses |
| Initial helium | `Y = 0.2703` | Helium mass fraction |
| Initial metallicity | `Z = 0.0142` | Metal mass fraction |
| Type-1 opacity metallicity | `Zbase = 0.0142` | Kept consistent with initial metallicity |
| Mixing-length model | `Cox` | MLT formulation |
| Mixing-length parameter | `2.0` | Ratio of mixing length to pressure scale height |
| Overshoot parameters | `f = 0.0001`, `f0 = 0.0005` | Supplied exponential overshoot controls |
| Mesh coefficient | `1.0` | Spatial resolution control |
| Brunt-Väisälä smoothing cells | `0` | Preserves the unsmoothed profile |
| Termination condition | central H-1 below `0.35` | Operational MAMS definition |
| Pulse format | `GYRE` | Produces an oscillation-ready stellar model |
| Center point | enabled | Adds the stellar center to GYRE data |
| Composition-gradient double points | enabled | Improves representation of sharp gradients |

Any scientific comparison should record changes to these values together with
the MESA version.

## Generate profiles with the notebook

### 1. Edit the parameter cell

The main editable cell contains:

```python
MASSES_MSUN = [1.0, 1.5, 2.0]
INITIAL_Y = 0.2703
INITIAL_Z = 0.0142
CENTRAL_H1_LIMIT = 0.35

PARTITION = "typhon"
WALLTIME = "7-00:00:00"
CPUS_PER_TASK = 1

OVERWRITE = False
PREPARE_RUNS = True
RUN_LOCAL = False
SUBMIT_SLURM = False
```

`MASSES_MSUN` accepts integers or decimal values in solar masses. Examples:

```python
MASSES_MSUN = [0.8, 1.0, 1.2, 1.4]
```

or

```python
MASSES_MSUN = [3.0, 3.01, 3.5, 4.0, 5.0]
```

The output-label convention replaces a decimal point with `p`:

| Mass | Label |
|---:|---|
| `0.15` | `MAMS0p15Msun` |
| `1.0` | `MAMS1Msun` |
| `1.5` | `MAMS1p5Msun` |
| `3.01` | `MAMS3p01Msun` |

### 2. Preview the rendered inlist

The preview cell builds the first inlist in memory and displays only the
settings changed by the pipeline. Check the mass, composition, opacity
metallicity, filenames, and central-hydrogen limit before creating directories.

### 3. Prepare the run directories

With `PREPARE_RUNS = True`, the notebook creates:

```text
results/mams_profiles/
├── manifest.json
├── MAMS1Msun/
│   ├── inlist
│   ├── inlist_project
│   ├── inlist_pgstar
│   ├── metadata.json
│   ├── run.sh
│   ├── submit.slurm
│   └── star -> ../../../tutorial/star
├── MAMS1p5Msun/
└── MAMS2Msun/
```

Running the setup cell again with identical settings reuses the existing
directories. If a directory with the same mass label contains different
settings, the pipeline raises an error and protects the existing calculation.

Set `OVERWRITE = True` only when replacing that generated directory and every
result inside it is intentional. Return it to `False` after replacement.

### 4. Choose an execution mode

#### Local execution

Set:

```python
RUN_LOCAL = True
SUBMIT_SLURM = False
```

The notebook runs the grid sequentially and waits for each model. Local runs
can take substantial time, especially for low-mass stars.

#### Slurm execution

Set:

```python
RUN_LOCAL = False
SUBMIT_SLURM = True
```

Update `PARTITION`, `WALLTIME`, and `CPUS_PER_TASK` for the collaborator's
cluster. Each model receives an independent Slurm script and working directory,
so simultaneous jobs write to separate locations.

The supplied MESA executable is normally a single-task program. Begin with
`CPUS_PER_TASK = 1`. Increase it only after confirming that the selected MESA
build uses threaded libraries that benefit from additional cores.

Useful Slurm commands are:

```bash
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,State,Elapsed,ExitCode
scancel JOB_ID
```

### 5. Inspect completion status

The final notebook section displays a status record for every mass. A completed
model should have both scientific outputs:

```text
MAMS1p5Msun.profile
MAMS1p5Msun.mesa
```

The corresponding `mesa.out` should contain:

```text
termination code: xa_central_lower_limit
```

The final central hydrogen value may be slightly smaller than the requested
threshold because of the last timestep.

## Generate profiles from the command line

The notebook calls a reusable Python module. The same setup can be scripted:

```bash
python3 scripts/mams_profiles/pipeline/mams_pipeline.py \
  --masses 0.8 1.0 1.2 1.4
```

This command prepares directories and a manifest. Add one execution option:

```bash
# Run sequentially on the current machine.
python3 scripts/mams_profiles/pipeline/mams_pipeline.py \
  --masses 0.8 1.0 1.2 --run-local

# Submit independent Slurm jobs.
python3 scripts/mams_profiles/pipeline/mams_pipeline.py \
  --masses 0.8 1.0 1.2 \
  --partition typhon \
  --walltime 7-00:00:00 \
  --cpus-per-task 1 \
  --submit-slurm
```

Show every option with:

```bash
python3 scripts/mams_profiles/pipeline/mams_pipeline.py --help
```

## Understand the generated files

Each mass directory contains:

| Path | Purpose |
|---|---|
| `inlist_project` | Exact mass-specific MESA controls used by the run |
| `metadata.json` | Machine-readable physical and scheduler parameters |
| `run.sh` | Local launcher that records `mesa.out` and `mesa.err` |
| `submit.slurm` | Independent Slurm job script |
| `LOGS/history.data` | Evolution history sampled according to MESA controls |
| `LOGS/profile*.data` | Profiles written by MESA during the run |
| `photos/` | Restart snapshots created by MESA |
| `MAMS<mass>Msun.profile` | Final regular MESA profile |
| `MAMS<mass>Msun.mesa` | Final GYRE-format model |

The first line of a GYRE-format file contains the number of points, stellar
mass in grams, radius in centimeters, luminosity in erg per second, and the
format version. The regular `.profile` file contains named global quantities
followed by zone-by-zone structure and composition data.

## Create a different stellar grid

### Change only the mass

Edit `MASSES_MSUN`. The pipeline handles filenames and isolated directories.
This is the standard way to add another stellar mass while preserving the
supplied physics.

### Change the composition

Edit `INITIAL_Y` and `INITIAL_Z` together with a scientific reason. The
pipeline updates both `initial_z` and opacity `Zbase`, keeping them consistent.
It also checks

$$
X_{\mathrm{initial}} = 1-Y-Z > 0,
$$

where all symbols were defined in the background section. Choose a central
hydrogen limit below the resulting initial hydrogen fraction.

Composition changes can require opacity tables or network considerations.
Record the adopted abundance pattern, because $Z$ alone leaves the mixture
unspecified.

### Change the main-sequence age definition

`CENTRAL_H1_LIMIT` controls the stopping point. Larger positive values stop
earlier; smaller positive values stop later. Approximate exploratory choices
include:

| Central H-1 limit | Qualitative stage |
|---:|---|
| `0.60` | early main sequence |
| `0.35` | this project's MAMS definition |
| `0.05` | late main sequence |

These labels are operational choices rather than universal definitions. For a
precise ZAMS or TAMS definition, select a physical criterion from the research
question and implement the corresponding MESA termination control.

Post-main-sequence, helium-burning, white-dwarf, or compact-object models use
different initial models, networks, mass-loss prescriptions, and termination
conditions. Build those as separate documented pipelines so the meaning of
each output label stays clear.

### Change convection, overshooting, or the network

Edit `tutorial/inlist_project` directly. Important controls include:

```fortran
new_net_name = 'mesa_49.net'
mixing_length_alpha = 2.0
MLT_option = 'Cox'
overshoot_f(:) = 0.0001
overshoot_f0(:) = 0.0005
```

After changing the shared template, use a new output root or intentionally
replace earlier generated directories. Commit the inlist change so every
collaborator can reproduce the new models.

## Troubleshooting

### `MESA_DIR` is unset

Define the variable in the current shell and compile again:

```bash
export MESA_DIR=/absolute/path/to/mesa
cd tutorial
./mk
```

Adding the export command to the shell startup file makes it available in new
terminals.

### `tutorial/star` is missing or lacks executable permission

Run `./mk` inside `tutorial/` with the MESA SDK environment active. Inspect the
compiler output for the first error.

### The inlist contains an unknown control

The collaborator is likely using a different MESA release. Search that
release's defaults files under `$MESA_DIR/star/defaults/` and compare the
official changelog. Reproducing the original build with version identifier
`6d0b0f2` gives the closest control set.

### A prepared directory already exists with different settings

Read its `metadata.json` and preserve useful outputs. Then choose one of two
safe paths:

1. change `OUTPUT_ROOT` to create a new experiment; or
2. set `OVERWRITE = True` for one intentional replacement.

### A Slurm job exits immediately

Read both `slurm-JOB_ID.err` and `mesa.err`. Confirm that compute nodes see the
same repository path, MESA installation, runtime libraries, and environment
initialization as the login node. Cluster module commands can be added near the
top of each generated `submit.slurm` file when required.

### MESA stops before producing both final files

Read the final part of `mesa.out`, then inspect `mesa.err` and
`LOGS/history.data`. Common causes include a timestep limit, convergence
failure, an incompatible inlist control, or missing runtime data. Preserve the
failed directory because its logs contain the evidence needed to diagnose the
model.

### PGSTAR creates display problems on a batch node

Set the PGSTAR window flags to `.false.` in `tutorial/inlist_pgstar`, regenerate
the run directories, and submit again. The scientific evolution and file
outputs remain available in headless mode.

## Reproducibility checklist

For every shared stellar grid, record:

- the Git commit of this repository;
- the MESA version or release;
- the MESA SDK and compiler versions;
- all masses;
- $Y$, $Z$, and the abundance mixture;
- the central-hydrogen termination limit;
- the network, convection, and overshooting settings;
- any local changes to `run_star_extras.f90`;
- whether models ran locally or through a scheduler;
- the contents of `manifest.json` and each `metadata.json`.

Generated profiles can be distributed as a versioned GitHub Release or a
research-data archive. Keeping source code in Git and large model grids in a
data release makes changes easier to review.

## References

- [MESA documentation](https://docs.mesastar.org/)
- [Installing MESA](https://docs.mesastar.org/en/latest/installation.html)
- [Running MESA](https://docs.mesastar.org/en/latest/using_mesa/running.html)
- [GYRE documentation](https://gyre.readthedocs.io/)
