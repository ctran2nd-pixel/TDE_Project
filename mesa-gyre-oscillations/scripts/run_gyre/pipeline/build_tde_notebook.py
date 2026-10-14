#!/usr/bin/env python3
"""Deterministically generate the GYRE-to-TDECalculator handoff notebook.

Run from the repository root (or any directory) with::

    python scripts/run_gyre/pipeline/build_tde_notebook.py

The generated notebook is intentionally output-free and safe to execute in a
fresh clone with neither GYRE results nor a collaborator TDE project present.
"""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf


def repository_root() -> Path:
    """Resolve the project root relative to this generator."""

    return Path(__file__).resolve().parents[3]


def make_notebook() -> nbf.NotebookNode:
    """Build the concise, executable human interface for the bridge."""

    cells = [
        nbf.v4.new_markdown_cell(
            """# GYRE to TDECalculator handoff

This notebook prepares clean GYRE HDF5 output for the collaborator's
unmodified `TDECalculator.py`.  It is an interface adapter, not a second
tidal-disruption implementation.  The calculator assumes adiabatic,
nonrotating quadrupole (`l=2`) modes, reads dimensionless `summary['omega']`,
and loads a legacy `M-gamma.txt` table relative to its own project directory.

Fresh execution is safe in a bare clone: absent GYRE results or absent TDE
project files are reported and the corresponding work is skipped."""
        ),
        nbf.v4.new_markdown_cell(
            """## Contract and scientific limits

The bridge selects only `l=2` rows, requires one `n_pg=0` f mode, rejects
duplicate radial orders, validates the fields and common radial grid required
by the calculator, then stages legacy filenames `detail.l2.n±N.h5`.  It never
changes the GYRE source directory and will not overwrite a staging directory.
The external calculator expects `M_star` and `M_r` in g, `R_star` in cm,
`rho` in g cm^-3, and dimensionless `x=r/R_star`.

Passing this file contract is **not** a validation of the TDE model.  The
collaborator's `M-gamma.txt` is a mass-only calibration tabulated for MAMS
models; it may not apply to an RGB structure even if its HDF5 files stage
successfully.  Confirm or recompute the gamma/criterion with the collaborator
before interpreting RGB results.

The staging manifest records only input handoff provenance. It is not amended
with external calculation results, which remain the caller's responsibility to
save and document.

Dense RGB mixed-mode spectra can also make the external calculation expensive:
its time/mode work is roughly O(`N² × n_modes`) and it allocates arrays over
time, modes, and radial grid.  Use the optional radial-order limits only as
explicit exploratory subsets, then demonstrate convergence both with the
`l=2` GYRE scan and with the calculator's `N` time samples."""
        ),
        nbf.v4.new_code_cell(
            """from pathlib import Path
import sys
import tempfile

def locate_repo_root(start: Path) -> Path:
    for candidate in (start, *start.parents):
        marker = candidate / "scripts" / "run_gyre" / "pipeline" / "tde_bridge.py"
        if marker.is_file():
            return candidate
    raise FileNotFoundError("Run this notebook from inside the cloned repository.")

REPO_ROOT = locate_repo_root(Path.cwd().resolve())
PIPELINE_DIR = REPO_ROOT / "scripts" / "run_gyre" / "pipeline"
if str(PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(PIPELINE_DIR))

from tde_bridge import (
    inspect_gyre_run,
    plot_quadrupole_spectrum,
    prepare_tde_directory,
    preflight_tde_project,
    run_tde_calculator,
)
"""
        ),
        nbf.v4.new_markdown_cell("## Editable parameters"),
        nbf.v4.new_code_cell(
            """# Paths may be repository-relative or absolute. These defaults intentionally
# name no personal machine or collaborator checkout.
GYRE_RUN_DIR = Path("results/run_gyre/my_star")
TDE_ROOT = Path("TDECalculator")

# Leave None to stage in a TemporaryDirectory that this notebook removes at the end.
# Set a new, non-existing path to keep a manifest-backed compatibility directory.
PERSISTENT_OUTPUT_DIR = None
TRANSFER_MODE = "copy"  # "copy" is portable; "hardlink" and "symlink" retain source dependence.

# Optional inclusive, explicitly recorded exploratory subset. Keep None for all l=2 modes.
N_PG_MIN = None
N_PG_MAX = None

# The external calculator is opt-in because it can be costly for dense spectra.
MBH = 1.0e6  # solar masses
Rp = 17.0    # GM_BH units
N = 1000     # orbit time samples; establish convergence before production use
RUN_TDE = False
"""
        ),
        nbf.v4.new_markdown_cell("## Read-only GYRE preflight"),
        nbf.v4.new_code_cell(
            """def resolve_from_repo(path: Path) -> Path:
    return path.expanduser().resolve() if path.is_absolute() else (REPO_ROOT / path).resolve()

SOURCE_DIR = resolve_from_repo(GYRE_RUN_DIR)
inspection = None
if not SOURCE_DIR.is_dir():
    print(f"No GYRE run directory yet: {SOURCE_DIR}")
    print("Run GYRE first, then set GYRE_RUN_DIR to its output directory.")
else:
    try:
        inspection = inspect_gyre_run(SOURCE_DIR, n_pg_min=N_PG_MIN, n_pg_max=N_PG_MAX)
        print(f"Validated {len(inspection.modes)} selected l=2 modes from {inspection.summary_path}")
        print(f"n_pg handoff: {list(inspection.selected_n_pg)}")
        print("TDE uses dimensionless omega; summary display-frequency units are not used.")
    except (OSError, ValueError) as exc:
        print(f"GYRE preflight did not pass: {exc}")
"""
        ),
        nbf.v4.new_markdown_cell("## Mode handoff diagnostic"),
        nbf.v4.new_code_cell(
            """if inspection is None:
    print("No validated quadrupole spectrum to plot.")
else:
    fig, ax, spectrum_report = plot_quadrupole_spectrum(inspection)
    print(
        f"Plotted {spectrum_report.plotted_modes} modes; "
        f"masked non-finite modes: {spectrum_report.masked_nonfinite_modes}."
    )
"""
        ),
        nbf.v4.new_markdown_cell("## TDE project preflight"),
        nbf.v4.new_code_cell(
            """PROJECT_ROOT = resolve_from_repo(TDE_ROOT)
tde_project = None
if not PROJECT_ROOT.is_dir():
    print(f"No TDE project directory yet: {PROJECT_ROOT}")
    print("Set TDE_ROOT to a collaborator clone containing TDECalculator.py and M-gamma.txt.")
else:
    try:
        tde_project = preflight_tde_project(PROJECT_ROOT)
        print(f"TDECalculator module: {tde_project.calculator_path}")
        print(f"Relative gamma table: {tde_project.gamma_table_path}")
        print(f"M-gamma calibration domain: {tde_project.gamma_mass_domain_solar} M_sun")
        print(f"TDECalculator SHA-256: {tde_project.calculator_sha256}")
        print(f"M-gamma SHA-256: {tde_project.gamma_table_sha256}")
    except (OSError, ValueError) as exc:
        print(f"TDE project preflight did not pass: {exc}")
"""
        ),
        nbf.v4.new_markdown_cell("## Stage a fresh legacy-compatible directory"),
        nbf.v4.new_code_cell(
            """prepared = None
_temporary_stage = None
if inspection is None:
    print("Skipping staging until the GYRE preflight succeeds.")
else:
    if PERSISTENT_OUTPUT_DIR is None:
        _temporary_stage = tempfile.TemporaryDirectory(prefix="gyre-tde-")
        stage_dir = Path(_temporary_stage.name) / "star"
        print(f"Using temporary stage: {stage_dir}")
    else:
        stage_dir = resolve_from_repo(Path(PERSISTENT_OUTPUT_DIR))
        print(f"Creating persistent stage (must not already exist): {stage_dir}")
    try:
        prepared = prepare_tde_directory(
            SOURCE_DIR,
            stage_dir,
            n_pg_min=N_PG_MIN,
            n_pg_max=N_PG_MAX,
            transfer_mode=TRANSFER_MODE,
        )
        print(f"Prepared {prepared.directory}")
        print(f"Manifest: {prepared.manifest_path}")
    except (OSError, ValueError) as exc:
        if _temporary_stage is not None:
            _temporary_stage.cleanup()
            _temporary_stage = None
        print(f"Staging did not complete: {exc}")
"""
        ),
        nbf.v4.new_markdown_cell("## Optional external TDE calculation"),
        nbf.v4.new_code_cell(
            """calculator_result = None
if not RUN_TDE:
    print("TDE execution is disabled. Set RUN_TDE=True only after calibration and convergence review.")
elif prepared is None:
    print("Skipping TDE execution because no staged input directory is available.")
elif tde_project is None:
    print("Skipping TDE execution because the TDE project preflight did not pass.")
else:
    # The bridge changes cwd to TDE_ROOT only while constructing the collaborator
    # calculator, so its legacy relative M-gamma.txt lookup works and cwd is restored.
    calculator_result = run_tde_calculator(
        prepared.directory, tde_project.root, MBH=MBH, Rp=Rp, N=N
    )
    print(f"Constructed external TDECalculator for {prepared.directory}")
"""
        ),
        nbf.v4.new_code_cell(
            """# Default staging is deliberately ephemeral. Persistent output is retained only
# when PERSISTENT_OUTPUT_DIR was explicitly set above.
if _temporary_stage is not None:
    _temporary_stage.cleanup()
    _temporary_stage = None
    print("Removed the temporary compatibility directory.")
"""
        ),
    ]
    notebook = nbf.v4.new_notebook(cells=cells)
    # nbformat otherwise generates random cell IDs, which makes a mechanically
    # regenerated notebook change even when its content has not.  Stable IDs
    # keep this tracked human interface byte-reproducible.
    for index, cell in enumerate(cells, start=1):
        cell["id"] = f"gyre-tde-{index:02d}"
    notebook.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.9"},
    }
    return notebook


def main() -> int:
    """Write the tracked, deterministic notebook with no saved execution state."""

    output_path = repository_root() / "scripts" / "run_gyre" / "notebook" / "GYRE_to_TDECalculator.ipynb"
    nbf.write(make_notebook(), output_path)
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
