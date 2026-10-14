"""Shared definitions and analysis helpers for the clean RGB GYRE notebooks.

Inputs
------
One MESA pulse-data model per RGB snapshot, written with
``pulse_data_format = 'GYRE'``, plus a compiled GYRE executable for real runs.

Outputs
-------
Radial and nonradial run directories below ``results/run_gyre/rgb_clean``,
including GYRE HDF5 files, CSV mode tables, metadata, and optional PDF/PNG
figures.  Import this module from the notebooks; it is not a standalone CLI.

The case defaults are deliberately starting points.  Red-giant mixed-mode
spectra can be extremely dense, so production results require scan- and
spatial-grid convergence tests.  The legacy notebooks targeted quadrupolar
(``l=2``) modes for tidal/TDE work.  ``science_mode='asteroseismology'`` adds
separate ``l=0`` and ``l=(1,2)`` scans; ``'quadrupole_tides'`` keeps only the
legacy physical degree without claiming a complete solar-like spectrum.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

try:  # Package import when used as scripts.run_gyre.pipeline.
    from .gyre_pipeline import (
        RunConfig,
        RunResult,
        load_gyre_hdf5,
        load_summary,
        resolve_gyre_executable,
        run_pipeline,
        validate_model_file,
    )
except ImportError:  # Direct import after adding pipeline/ to sys.path.
    from gyre_pipeline import (
        RunConfig,
        RunResult,
        load_gyre_hdf5,
        load_summary,
        resolve_gyre_executable,
        run_pipeline,
        validate_model_file,
    )


@dataclass(frozen=True)
class RGBCase:
    """Physical label, expected input, and numerical starting point."""

    key: str
    label: str
    mass_msun: float
    radius_rsun: float
    metallicity: float
    model_filename: str
    freq_min_uhz: float
    freq_max_uhz: float
    n_pg_min: int
    n_pg_max: int
    radial_n_freq: int
    nonradial_n_freq: int
    note: str


# Radii encoded as logR0p5, logR1p3, and logR1p6 are 10**0.5,
# 10**1.3, and 10**1.6 R_sun.  In particular, the legacy "3.61R" filename
# was a label error: its saved model radius was 2.20109339e11 cm, or 3.164 R_sun.
RGB_CASES: Dict[str, RGBCase] = {
    "rgb_1m_z0p0142_3p16r": RGBCase(
        key="rgb_1m_z0p0142_3p16r",
        label=r"$1\,M_\odot$, $Z=0.0142$, $R\simeq3.16\,R_\odot$ RGB",
        mass_msun=1.0,
        radius_rsun=10.0**0.5,
        metallicity=0.0142,
        model_filename="RGB1Msun_Z0p0142_logR0p5.mesa",
        freq_min_uhz=80.0,
        freq_max_uhz=500.0,
        n_pg_min=-1000,
        n_pg_max=50,
        radial_n_freq=1000,
        nonradial_n_freq=12000,
        note=(
            "Corrected from the legacy notebook label 3.61R; logR=0.5 and "
            "the saved R_star both indicate about 3.16 R_sun."
        ),
    ),
    "rgb_1m_z0p0142_10r": RGBCase(
        key="rgb_1m_z0p0142_10r",
        label=r"$1\,M_\odot$, $Z=0.0142$, $R\simeq10\,R_\odot$ RGB",
        mass_msun=1.0,
        radius_rsun=10.0,
        metallicity=0.0142,
        model_filename="RGB1Msun_Z0p0142_logR1p0.mesa",
        freq_min_uhz=250.0,
        freq_max_uhz=500,
        n_pg_min=-4000,
        n_pg_max=150,
        radial_n_freq=800,
        nonradial_n_freq=20000,
        note="The broad band brackets the expected solar-like power region near tens of microhertz.",
    ),
    "rgb_1m_z0p0142_20r": RGBCase(
        key="rgb_1m_z0p0142_20r",
        label=r"$1\,M_\odot$, $Z=0.0142$, $R\simeq20\,R_\odot$ RGB",
        mass_msun=1.0,
        radius_rsun=10.0**1.3,
        metallicity=0.0142,
        model_filename="RGB1Msun_Z0p0142_logR1p3.mesa",
        freq_min_uhz=2.0,
        freq_max_uhz=30.0,
        n_pg_min=-12000,
        n_pg_max=50,
        radial_n_freq=600,
        nonradial_n_freq=30000,
        note="The exact encoded radius is 19.95 R_sun; mixed modes are denser than in the 10 R_sun case.",
    ),
    "rgb_1m_z0p0142_40r": RGBCase(
        key="rgb_1m_z0p0142_40r",
        label=r"$1\,M_\odot$, $Z=0.0142$, $R\simeq40\,R_\odot$ RGB",
        mass_msun=1.0,
        radius_rsun=10.0**1.6,
        metallicity=0.0142,
        model_filename="RGB1Msun_Z0p0142_logR1p6.mesa",
        freq_min_uhz=0.5,
        freq_max_uhz=10.0,
        n_pg_min=-40000,
        n_pg_max=50,
        radial_n_freq=500,
        nonradial_n_freq=40000,
        note="The exact encoded radius is 39.81 R_sun; the very dense nonradial spectrum can be expensive.",
    ),
    "rgb_3m_z0p0142_10r": RGBCase(
        key="rgb_3m_z0p0142_10r",
        label=r"$3\,M_\odot$, $Z=0.0142$, $R\simeq10\,R_\odot$ RGB",
        mass_msun=3.0,
        radius_rsun=10.0,
        metallicity=0.0142,
        model_filename="RGB3Msun_Z0p0142_logR1p0.mesa",
        freq_min_uhz=1.0,
        freq_max_uhz=500.0,
        n_pg_min=-5000,
        n_pg_max=100,
        radial_n_freq=1000,
        nonradial_n_freq=50000,
        note="Do not use the legacy saved arrays: they report the 1 M_sun model mass and are stale.",
    ),
    "rgb_3m_z0p0142_40r": RGBCase(
        key="rgb_3m_z0p0142_40r",
        label=r"$3\,M_\odot$, $Z=0.0142$, $R\simeq40\,R_\odot$ RGB",
        mass_msun=3.0,
        radius_rsun=40.0,
        metallicity=0.0142,
        model_filename="RGB3Msun_Z0p0142_logR1p6.mesa",
        freq_min_uhz=175.0,
        freq_max_uhz=250.0,
        n_pg_min=-6000,
        n_pg_max=500,
        radial_n_freq=1000,
        nonradial_n_freq=20000,
        note="Do not use the legacy saved arrays: they report the 1 M_sun model mass and are stale.",
    ),
    "rgb_5m_z0p0142_40r": RGBCase(
        key="rgb_5m_z0p0142_40r",
        label=r"$5\,M_\odot$, $Z=0.0142$, $R\simeq40\,R_\odot$ RGB",
        mass_msun=5.0,
        radius_rsun=40.0,
        metallicity=0.0142,
        model_filename="RGB5Msun_Z0p0142_logR1p6.mesa",
        freq_min_uhz=20.0,
        freq_max_uhz=150.0,
        n_pg_min=-5000,
        n_pg_max=200,
        radial_n_freq=1000,
        nonradial_n_freq=40000,
        note="Do not use the legacy saved arrays: they report the 1 M_sun model mass and are stale.",
    ),

    "rgb_8m_z0p0142_100r": RGBCase(
        key="rgb_8m_z0p0142_100r",
        label=r"$8\,M_\odot$, $Z=0.0142$, $R\simeq100\,R_\odot$ RGB",
        mass_msun=8.0,
        radius_rsun=100.0,
        metallicity=0.0142,
        model_filename="RGB8Msun_Z0p0142_logR2p0.mesa",
        freq_min_uhz=15.0,
        freq_max_uhz=80.0,
        n_pg_min=-5000,
        n_pg_max=200,
        radial_n_freq=1000,
        nonradial_n_freq=40000,
        note="Do not use the legacy saved arrays: they report the 1 M_sun model mass and are stale.",
    ),
}


def get_case(key: str) -> RGBCase:
    """Return one named RGB case with an informative error for bad keys."""

    try:
        return RGB_CASES[key]
    except KeyError as exc:
        choices = ", ".join(sorted(RGB_CASES))
        raise KeyError(f"Unknown RGB case {key!r}. Choose one of: {choices}") from exc


def find_repo_root(start: Path | str) -> Path:
    """Find the repository root from the root, notebook directory, or a child."""

    resolved = Path(start).expanduser().resolve()
    for candidate in (resolved, *resolved.parents):
        marker = candidate / "scripts" / "run_gyre" / "pipeline" / "rgb_cases.py"
        if marker.is_file():
            return candidate
    raise FileNotFoundError("Run this notebook from inside the cloned repository.")


def build_run_configs(
    case: RGBCase,
    *,
    science_mode: str = "asteroseismology",
    freq_min_uhz: Optional[float] = None,
    freq_max_uhz: Optional[float] = None,
    n_pg_min: Optional[int] = None,
    n_pg_max: Optional[int] = None,
    radial_n_freq: Optional[int] = None,
    nonradial_n_freq: Optional[int] = None,
    nonradial_degrees: Sequence[int] = (1, 2),
    nonradial_grid_type: str = "INVERSE",
    w_osc: float = 10.0,
    w_exp: float = 2.0,
    w_ctr: float = 100.0,
    threads: int = 1,
) -> Dict[str, RunConfig]:
    """Build separate scan configurations for one RGB snapshot.

    ``asteroseismology`` returns a modest linear radial scan and a much finer
    nonradial scan.  ``quadrupole_tides`` returns only an ``l=2`` nonradial
    scan, matching the legacy notebooks' science intent.  The broad negative
    ``n_pg`` bounds prevent silent truncation but do not prove completeness;
    users must check monotonicity and repeat with higher resolution.
    """

    mode = science_mode.lower()
    if mode not in {"asteroseismology", "quadrupole_tides"}:
        raise ValueError(
            "science_mode must be 'asteroseismology' or 'quadrupole_tides'."
        )
    low = case.freq_min_uhz if freq_min_uhz is None else freq_min_uhz
    high = case.freq_max_uhz if freq_max_uhz is None else freq_max_uhz
    order_min = case.n_pg_min if n_pg_min is None else n_pg_min
    order_max = case.n_pg_max if n_pg_max is None else n_pg_max
    radial_points = case.radial_n_freq if radial_n_freq is None else radial_n_freq
    nonradial_points = (
        case.nonradial_n_freq if nonradial_n_freq is None else nonradial_n_freq
    )

    common = dict(
        n_pg_min=order_min,
        n_pg_max=order_max,
        freq_min=low,
        freq_max=high,
        freq_units="UHZ",
        output_freq_units="UHZ",
        w_osc=w_osc,
        w_exp=w_exp,
        w_ctr=w_ctr,
        threads=threads,
    )
    configs: Dict[str, RunConfig] = {}
    if mode == "asteroseismology":
        configs["radial"] = RunConfig(
            harmonic_degrees=(0,),
            n_freq=radial_points,
            grid_type="LINEAR",
            **common,
        ).validated()
        degrees = tuple(int(degree) for degree in nonradial_degrees)
        if not degrees or any(degree == 0 for degree in degrees):
            raise ValueError("nonradial_degrees must contain positive degrees.")
    else:
        degrees = (2,)
    configs["nonradial"] = RunConfig(
        harmonic_degrees=degrees,
        n_freq=nonradial_points,
        grid_type=nonradial_grid_type,
        **common,
    ).validated()
    return configs


def case_paths(
    repo_root: Path | str,
    case: RGBCase,
    *,
    model_path: Optional[Path | str] = None,
    output_root: Optional[Path | str] = None,
) -> Tuple[Path, Path]:
    """Resolve the expected pulse model and canonical case result directory."""

    root = Path(repo_root).expanduser().resolve()
    model = (
        root / "models" / case.model_filename
        if model_path is None
        else Path(model_path).expanduser()
    )
    output = (
        root / "results" / "run_gyre" / "rgb_clean" / case.key
        if output_root is None
        else Path(output_root).expanduser()
    )
    if not model.is_absolute():
        model = root / model
    if not output.is_absolute():
        output = root / output
    return model.resolve(), output.resolve()


def preflight_case(
    repo_root: Path | str,
    case: RGBCase,
    *,
    model_path: Optional[Path | str] = None,
    output_root: Optional[Path | str] = None,
    gyre_executable: Optional[Path | str] = None,
) -> Dict[str, Any]:
    """Inspect inputs without writing files or launching GYRE."""

    model, output = case_paths(
        repo_root, case, model_path=model_path, output_root=output_root
    )
    report: Dict[str, Any] = {
        "case": case.key,
        "model_path": str(model),
        "model_available": model.is_file(),
        "model_valid": False,
        "output_root": str(output),
        "gyre_available": False,
        "gyre_executable": None,
        "messages": [],
    }
    if model.is_file():
        try:
            report["model_info"] = validate_model_file(model)
        except (OSError, ValueError) as exc:
            report["messages"].append(f"Model validation failed: {exc}")
        else:
            report["model_valid"] = True
    else:
        report["messages"].append(
            "Missing MESA pulse model. Copy a GYRE-format pulse file to "
            f"models/{case.model_filename} or edit MODEL_PATH."
        )
    try:
        executable = resolve_gyre_executable(gyre_executable)
    except FileNotFoundError as exc:
        report["messages"].append(str(exc))
    else:
        report["gyre_available"] = True
        report["gyre_executable"] = str(executable)
    return report


def run_case_plan(
    model_path: Path | str,
    output_root: Path | str,
    configs: Mapping[str, RunConfig],
    *,
    gyre_executable: Optional[Path | str] = None,
    overwrite: bool = True,
    dry_run: bool = False,
) -> Dict[str, RunResult]:
    """Run or prepare each independent radial/nonradial scan."""

    root = Path(output_root)
    results: Dict[str, RunResult] = {}
    for name, config in configs.items():
        results[name] = run_pipeline(
            model_path,
            root / name,
            config,
            gyre_executable=gyre_executable,
            overwrite=overwrite,
            dry_run=dry_run,
        )
    return results


def load_case_summaries(
    output_root: Path | str, scan_names: Sequence[str] = ("radial", "nonradial")
) -> Dict[str, Dict[str, Any]]:
    """Load whichever planned summaries currently exist."""

    root = Path(output_root)
    summaries: Dict[str, Dict[str, Any]] = {}
    for name in scan_names:
        path = root / name / "summary.h5"
        if path.is_file():
            summaries[name] = load_summary(path)
    return summaries


def _real(values: Any) -> np.ndarray:
    """Return a real NumPy view for real or complex GYRE arrays."""

    return np.real(np.asarray(values))


def mode_quality_report(
    summary: Mapping[str, Any],
    config: RunConfig,
    *,
    chi_limit: float = 1.0e-8,
    iteration_limit: int = 20,
) -> Dict[str, Any]:
    """Check scan boundaries, ordering, duplication, and root convergence.

    For adiabatic calculations, ``n_pg`` should increase with increasing
    frequency within each harmonic degree.  A violation is a useful warning
    of missed roots or inadequate resolution, not an automatic repair.
    """

    required = {"l", "n_pg", "freq"}
    missing = sorted(required.difference(summary))
    if missing:
        raise KeyError(f"Summary is missing required datasets: {missing}")

    degree_values = np.asarray(summary["l"], dtype=int)
    orders = np.asarray(summary["n_pg"], dtype=int)
    frequencies = _real(summary["freq"])
    report: Dict[str, Any] = {
        "mode_count": int(len(degree_values)),
        "by_degree": {},
        "boundary_hits": [],
        "max_chi": None,
        "max_n_iter": None,
        "warnings": [],
    }
    if len(degree_values) == 0:
        report["warnings"].append("No modes were found in this scan.")
        return report
    for degree in sorted(set(degree_values.tolist())):
        mask = degree_values == degree
        local_orders = orders[mask]
        local_frequencies = frequencies[mask]
        sequence = local_orders[np.argsort(local_frequencies)]
        differences = np.diff(sequence)
        duplicate_count = int(np.count_nonzero(differences == 0))
        nonmonotonic_count = int(np.count_nonzero(differences <= 0))
        degree_report = {
            "count": int(mask.sum()),
            "n_pg_min_found": int(local_orders.min()),
            "n_pg_max_found": int(local_orders.max()),
            "duplicate_steps": duplicate_count,
            "nonmonotonic_steps": nonmonotonic_count,
        }
        report["by_degree"][int(degree)] = degree_report
        if nonmonotonic_count:
            report["warnings"].append(
                f"l={degree}: n_pg is not strictly increasing with frequency "
                f"at {nonmonotonic_count} step(s)."
            )
        for boundary in (config.n_pg_min, config.n_pg_max):
            if np.any(local_orders == boundary):
                report["boundary_hits"].append({"l": int(degree), "n_pg": boundary})

    if report["boundary_hits"]:
        report["warnings"].append(
            "At least one mode reached an n_pg filter boundary; widen the bounds "
            "before claiming completeness."
        )
    if "chi" in summary and np.asarray(summary["chi"]).size:
        report["max_chi"] = float(np.max(np.abs(np.asarray(summary["chi"]))))
        if report["max_chi"] > chi_limit:
            report["warnings"].append(
                f"Maximum chi={report['max_chi']:.3e} exceeds {chi_limit:.1e}."
            )
    else:
        report["warnings"].append("Summary has no chi dataset; root residuals are unchecked.")
    if "n_iter" in summary and np.asarray(summary["n_iter"]).size:
        report["max_n_iter"] = int(np.max(np.asarray(summary["n_iter"])))
        if report["max_n_iter"] > iteration_limit:
            report["warnings"].append(
                f"Maximum n_iter={report['max_n_iter']} exceeds {iteration_limit}."
            )
    else:
        report["warnings"].append(
            "Summary has no n_iter dataset; iteration convergence is unchecked."
        )
    return report


def stellar_parameter_report(
    summary: Mapping[str, Any],
    case: RGBCase,
    *,
    mass_tolerance_fraction: float = 0.15,
    radius_tolerance_fraction: float = 0.03,
) -> Dict[str, Any]:
    """Compare GYRE's model mass/radius with the case label.

    ``M_star`` and ``R_star`` are cgs values in GYRE summary files.  The mass
    tolerance is intentionally loose because an evolved model's current mass
    need not equal its initial mass after mass loss; a large discrepancy is
    still useful for catching a stale or mismatched pulse file.
    """

    solar_mass_g = 1.98847e33
    solar_radius_cm = 6.957e10
    report: Dict[str, Any] = {
        "mass_msun": None,
        "radius_rsun": None,
        "warnings": [],
    }
    if "M_star" in summary and np.asarray(summary["M_star"]).size:
        mass = float(np.median(_real(summary["M_star"]))) / solar_mass_g
        report["mass_msun"] = mass
        mass_difference = abs(mass - case.mass_msun) / case.mass_msun
        if mass_difference > mass_tolerance_fraction:
            report["warnings"].append(
                f"Summary mass {mass:.4g} M_sun differs from the nominal "
                f"{case.mass_msun:.4g} M_sun case by {mass_difference:.1%}."
            )
    else:
        report["warnings"].append("Summary has no M_star dataset.")
    if "R_star" in summary and np.asarray(summary["R_star"]).size:
        radius = float(np.median(_real(summary["R_star"]))) / solar_radius_cm
        report["radius_rsun"] = radius
        radius_difference = abs(radius - case.radius_rsun) / case.radius_rsun
        if radius_difference > radius_tolerance_fraction:
            report["warnings"].append(
                f"Summary radius {radius:.4g} R_sun differs from the nominal "
                f"{case.radius_rsun:.4g} R_sun case by {radius_difference:.1%}."
            )
    else:
        report["warnings"].append("Summary has no R_star dataset.")
    return report


def configure_plot_style() -> None:
    """Apply a title-free, publication-readable Matplotlib style."""

    plt.rcParams.update(
        {
            "font.size": 16,
            "axes.labelsize": 20,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 15,
            "axes.titlesize": 0,
            "lines.linewidth": 2.2,
            "lines.markersize": 6,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "mathtext.fontset": "stix",
            "font.family": "STIXGeneral",
        }
    )


def plot_mode_spectrum(
    summaries: Mapping[str, Mapping[str, Any]],
    *,
    output_prefix: Optional[Path | str] = None,
) -> Tuple[Any, Any, Tuple[Path, ...]]:
    """Plot mode order against cyclic frequency and optionally save PDF/PNG."""

    configure_plot_style()
    fig, ax = plt.subplots(figsize=(7.5, 4.8), constrained_layout=True)
    plotted = False
    colors = {0: "#222222", 1: "#0072B2", 2: "#D55E00", 3: "#009E73"}
    for summary in summaries.values():
        if not {"l", "n_pg", "freq"}.issubset(summary):
            continue
        degrees = np.asarray(summary["l"], dtype=int)
        frequencies = _real(summary["freq"])
        orders = np.asarray(summary["n_pg"], dtype=int)
        for degree in sorted(set(degrees.tolist())):
            mask = degrees == degree
            ax.scatter(
                frequencies[mask],
                orders[mask],
                s=38,
                alpha=0.8,
                color=colors.get(int(degree)),
                label=rf"$\ell={int(degree)}$",
            )
            plotted = True
    if not plotted:
        plt.close(fig)
        raise ValueError("No l, n_pg, and freq arrays were available to plot.")
    ax.set_xlabel(r"Cyclic frequency $\nu\,(\mu\mathrm{Hz})$")
    ax.set_ylabel(r"Mode order $n_{\mathrm{pg}}$")
    ax.tick_params(which="both", direction="in")
    ax.grid(alpha=0.25)
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), frameon=False)

    saved: Tuple[Path, ...] = ()
    if output_prefix is not None:
        prefix = Path(output_prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        pdf_path = prefix.with_suffix(".pdf")
        png_path = prefix.with_suffix(".png")
        fig.savefig(pdf_path, bbox_inches="tight")
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        saved = (pdf_path, png_path)
    return fig, ax, saved


def plot_eigenfunctions(
    detail_path: Path | str,
    *,
    output_prefix: Optional[Path | str] = None,
) -> Tuple[Any, Any, Tuple[Path, ...]]:
    """Plot real radial and horizontal displacements from one detail file."""

    detail = load_gyre_hdf5(detail_path)
    required = {"x", "xi_r", "xi_h"}
    missing = sorted(required.difference(detail))
    if missing:
        raise KeyError(f"Detail file is missing required datasets: {missing}")
    configure_plot_style()
    fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
    radius = np.asarray(detail["x"])
    ax.plot(radius, _real(detail["xi_r"]), label=r"$\Re(\xi_r)$")
    ax.plot(radius, _real(detail["xi_h"]), label=r"$\Re(\xi_h)$")
    ax.set_xlabel(r"Fractional radius $r/R_\star$")
    ax.set_ylabel(r"Normalized displacement")
    ax.tick_params(which="both", direction="in")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    saved: Tuple[Path, ...] = ()
    if output_prefix is not None:
        prefix = Path(output_prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        pdf_path = prefix.with_suffix(".pdf")
        png_path = prefix.with_suffix(".png")
        fig.savefig(pdf_path, bbox_inches="tight")
        fig.savefig(png_path, dpi=300, bbox_inches="tight")
        saved = (pdf_path, png_path)
    return fig, ax, saved
