#!/usr/bin/env python3
"""Regenerate the five collaborator-facing clean RGB GYRE notebooks.

Run from anywhere inside the repository with Python 3.9 or newer::

    python scripts/run_gyre/pipeline/build_rgb_notebooks.py

The generated notebooks contain no execution counts or saved outputs.  Their
real GYRE switch defaults to false, so a missing collaborator-supplied MESA
pulse model produces an explanatory preflight message rather than an error.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import nbformat as nbf


CASES: Dict[str, str] = {
    "rgb_1m_z0p0142_3p16r": "GYRE_RGB1Msun_Z0p0142_3p16R_clean.ipynb",
    "rgb_1m_z0p0142_10r": "GYRE_RGB1Msun_Z0p0142_10R_clean.ipynb",
    "rgb_1m_z0p0142_20r": "GYRE_RGB1Msun_Z0p0142_20R_clean.ipynb",
    "rgb_1m_z0p0142_40r": "GYRE_RGB1Msun_Z0p0142_40R_clean.ipynb",
    "rgb_3m_z0p0142_10r": "GYRE_RGB3Msun_Z0p0142_10R_clean.ipynb",
}


def repository_root() -> Path:
    """Resolve the repository root relative to this generator."""

    return Path(__file__).resolve().parents[3]


def make_notebook(case_key: str) -> nbf.NotebookNode:
    """Build one guided RGB workflow using the shared case module."""

    cells = [
        nbf.v4.new_markdown_cell(
            f"""# Clean RGB oscillation workflow: `{case_key}`

This notebook prepares an adiabatic, nonrotating GYRE calculation from one
MESA pulse-data snapshot. It uses separate radial and nonradial scans for the
asteroseismology option because red-giant mixed modes are far denser than
radial modes. The legacy notebooks were instead oriented toward quadrupolar
tidal/TDE modes (`l=2`); select that intent explicitly in the parameter cell.

`RUN_GYRE` is `False` by default. Opening or running all cells therefore does
not launch GYRE or overwrite results."""
        ),
        nbf.v4.new_markdown_cell(
            """## Expected inputs and outputs

Run the notebook from anywhere inside this clone. The shared implementation
is in `scripts/run_gyre/pipeline/`; generated products go to
`results/run_gyre/rgb_clean/<case>/<scan>/`.

The repository does **not** contain the stellar pulse file. Copy the matching
MESA export into `models/` (or edit `MODEL_PATH` below). It must be written
with `pulse_data_format = 'GYRE'`; an ordinary `profile*.data` is not a
replacement.

After supplying a model and GYRE, set `RUN_GYRE=True`, restart the kernel, and
run all cells. Keep `OVERWRITE=False` unless replacing a pipeline-created run
is intentional."""
        ),
        nbf.v4.new_code_cell(
            """from pathlib import Path
import sys

def locate_repo_root(start):
    for candidate in (start, *start.parents):
        marker = candidate / "scripts" / "run_gyre" / "pipeline" / "rgb_cases.py"
        if marker.is_file():
            return candidate
    raise FileNotFoundError("Run this notebook from inside the cloned repository.")

REPO_ROOT = locate_repo_root(Path.cwd().resolve())
PIPELINE_DIR = REPO_ROOT / "scripts" / "run_gyre" / "pipeline"
if str(PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(PIPELINE_DIR))

from gyre_pipeline import render_gyre_input
from rgb_cases import (
    build_run_configs,
    case_paths,
    configure_plot_style,
    get_case,
    load_case_summaries,
    mode_quality_report,
    plot_eigenfunctions,
    plot_mode_spectrum,
    preflight_case,
    run_case_plan,
    stellar_parameter_report,
)

configure_plot_style()"""
        ),
        nbf.v4.new_markdown_cell(
            """## Editable parameters

Use `SCIENCE_MODE="asteroseismology"` for `l=0,1,2` starting scans. Use
`"quadrupole_tides"` to reproduce the legacy degree choice (`l=2` only).
Neither preset proves completeness: widen the order bounds and increase scan
and spatial resolution until scientifically important frequencies converge."""
        ),
        nbf.v4.new_code_cell(
            f'''CASE_KEY = "{case_key}"
CASE = get_case(CASE_KEY)

# Repository-relative defaults. Absolute paths are also accepted.
MODEL_PATH = Path("models") / CASE.model_filename
OUTPUT_ROOT = Path("results/run_gyre/rgb_clean") / CASE.key
GYRE_EXECUTABLE = None  # Example: Path("/path/to/gyre/bin/gyre"); otherwise use $GYRE_DIR.

RUN_GYRE = False
OVERWRITE = False
SCIENCE_MODE = "asteroseismology"  # Or "quadrupole_tides" for legacy l=2/TDE intent.

# Starting values from the shared case definition; edit for convergence tests.
FREQ_MIN_UHZ, FREQ_MAX_UHZ = CASE.freq_min_uhz, CASE.freq_max_uhz
N_PG_MIN, N_PG_MAX = CASE.n_pg_min, CASE.n_pg_max
RADIAL_N_FREQ = CASE.radial_n_freq
NONRADIAL_N_FREQ = CASE.nonradial_n_freq
NONRADIAL_DEGREES = (1, 2)
NONRADIAL_GRID_TYPE = "INVERSE"
W_OSC, W_EXP, W_CTR = 10.0, 2.0, 100.0
THREADS = 1'''
        ),
        nbf.v4.new_code_cell(
            """MODEL_PATH, OUTPUT_ROOT = case_paths(
    REPO_ROOT, CASE, model_path=MODEL_PATH, output_root=OUTPUT_ROOT
)
print(CASE.label)
print(f"Expected model: {MODEL_PATH}")
print(f"Target radius: {CASE.radius_rsun:.4g} R_sun")
print(f"Starting frequency band: {FREQ_MIN_UHZ:g}--{FREQ_MAX_UHZ:g} microhertz")
print(f"Case note: {CASE.note}")"""
        ),
        nbf.v4.new_markdown_cell("## Read-only preflight"),
        nbf.v4.new_code_cell(
            """preflight = preflight_case(
    REPO_ROOT,
    CASE,
    model_path=MODEL_PATH,
    output_root=OUTPUT_ROOT,
    gyre_executable=GYRE_EXECUTABLE,
)
for key in (
    "model_path", "model_available", "model_valid", "gyre_executable",
    "gyre_available", "output_root",
):
    print(f"{key}: {preflight[key]}")
for message in preflight["messages"]:
    print(f"NOTE: {message}")"""
        ),
        nbf.v4.new_markdown_cell(
            """## Construct and inspect the scan plan

The radial scan is linear and modest. The nonradial scan is much finer and
inverse-frequency by default to better sample the g-mode-rich spectrum.
`w_ctr=100` deliberately resolves the compact RGB core more strongly than the
legacy `w_ctr=10`."""
        ),
        nbf.v4.new_code_cell(
            """configs = build_run_configs(
    CASE,
    science_mode=SCIENCE_MODE,
    freq_min_uhz=FREQ_MIN_UHZ,
    freq_max_uhz=FREQ_MAX_UHZ,
    n_pg_min=N_PG_MIN,
    n_pg_max=N_PG_MAX,
    radial_n_freq=RADIAL_N_FREQ,
    nonradial_n_freq=NONRADIAL_N_FREQ,
    nonradial_degrees=NONRADIAL_DEGREES,
    nonradial_grid_type=NONRADIAL_GRID_TYPE,
    w_osc=W_OSC,
    w_exp=W_EXP,
    w_ctr=W_CTR,
    threads=THREADS,
)
for scan_name, config in configs.items():
    print(f"\\n--- {scan_name} scan ---")
    print(render_gyre_input(config))"""
        ),
        nbf.v4.new_markdown_cell("## Opt-in execution"),
        nbf.v4.new_code_cell(
            """run_results = {}
if not RUN_GYRE:
    print("GYRE execution is disabled. Set RUN_GYRE=True only after preflight passes.")
elif not preflight["model_valid"]:
    raise FileNotFoundError(
        f"Supply a valid MESA GYRE-format pulse model first: {MODEL_PATH}"
    )
elif not preflight["gyre_available"]:
    raise FileNotFoundError("Set GYRE_EXECUTABLE or export GYRE_DIR before running.")
else:
    run_results = run_case_plan(
        MODEL_PATH,
        OUTPUT_ROOT,
        configs,
        gyre_executable=GYRE_EXECUTABLE,
        overwrite=OVERWRITE,
        dry_run=False,
    )
    for scan_name, result in run_results.items():
        print(f"{scan_name}: {result.output_dir}")"""
        ),
        nbf.v4.new_markdown_cell("## Load and validate available mode catalogues"),
        nbf.v4.new_code_cell(
            """summaries = load_case_summaries(OUTPUT_ROOT, tuple(configs))
if not summaries:
    print("No summary files yet; complete a real GYRE run first.")
else:
    for scan_name, summary in summaries.items():
        report = mode_quality_report(summary, configs[scan_name])
        stellar = stellar_parameter_report(summary, CASE)
        print(f"\\n{scan_name}: {report['mode_count']} modes")
        print(f"by degree: {report['by_degree']}")
        print(f"max chi: {report['max_chi']}; max n_iter: {report['max_n_iter']}")
        print(f"model mass: {stellar['mass_msun']} M_sun; radius: {stellar['radius_rsun']} R_sun")
        for warning in stellar["warnings"]:
            print(f"WARNING: {warning}")
        if report["warnings"]:
            for warning in report["warnings"]:
                print(f"WARNING: {warning}")
        else:
            print("No automated quality warnings.")"""
        ),
        nbf.v4.new_markdown_cell(
            """## Mode spectrum

The figure is intentionally title-free and is saved as PDF and 300-dpi PNG
only when real summaries exist. A clean diagnostic requires monotonically
increasing `n_pg` with frequency within each `l`, no order-boundary hits, and
stable frequencies when `n_freq`, `w_osc`, and `w_ctr` are increased."""
        ),
        nbf.v4.new_code_cell(
            """if summaries:
    fig, ax, saved = plot_mode_spectrum(
        summaries, output_prefix=OUTPUT_ROOT / "mode_spectrum"
    )
    print("Saved:", *(str(path) for path in saved), sep="\\n  ")
else:
    print("Spectrum plot skipped: no summaries are available.")"""
        ),
        nbf.v4.new_markdown_cell("## One mode's eigenfunctions"),
        nbf.v4.new_code_cell(
            """detail_files = []
for scan_name in configs:
    detail_files.extend(sorted((OUTPUT_ROOT / scan_name).glob("detail.id*.l*.n*.h5")))

if detail_files:
    DETAIL_INDEX = 0  # Edit to inspect another mode.
    selected_detail = detail_files[DETAIL_INDEX]
    fig, ax, saved = plot_eigenfunctions(
        selected_detail, output_prefix=OUTPUT_ROOT / "eigenfunctions"
    )
    print(f"Selected: {selected_detail}")
    print("Saved:", *(str(path) for path in saved), sep="\\n  ")
else:
    print("Eigenfunction plot skipped: no detail files are available.")"""
        ),
        nbf.v4.new_markdown_cell(
            """## Interpretation checklist

- Confirm `R_star` and `M_star` in the new summary match this case; never infer
  them from a notebook filename or old saved output.
- Treat an `n_pg` boundary hit, duplicate, or nonmonotonic sequence as a reason
  to widen/refine the search.
- Repeat with at least twice `NONRADIAL_N_FREQ` and increased `W_CTR`; retain
  only frequencies stable at the precision needed by the science.
- `quadrupole_tides` is an `l=2` catalogue, not a complete solar-like RGB
  spectrum. Asteroseismic interpretation needs the radial and dipole modes as
  well as careful completeness checks.
- These calculations are adiabatic and nonrotating; they do not predict mode
  damping, growth rates, or rotational splitting."""
        ),
    ]
    for index, cell in enumerate(cells):
        # Stable cell IDs keep notebook regeneration reproducible and avoid
        # meaningless diffs from nbformat's random defaults.
        cell["id"] = f"{case_key}-{index:02d}"
    notebook = nbf.v4.new_notebook(cells=cells)
    notebook.metadata["kernelspec"] = {
        "display_name": "Python 3",
        "language": "python",
        "name": "python3",
    }
    notebook.metadata["language_info"] = {"name": "python", "version": "3.9"}
    return notebook


def main() -> int:
    """Write every declared clean notebook and print the artifact paths."""

    destination = repository_root() / "scripts" / "run_gyre" / "notebook"
    destination.mkdir(parents=True, exist_ok=True)
    for case_key, filename in CASES.items():
        path = destination / filename
        nbf.write(make_notebook(case_key), path)
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
