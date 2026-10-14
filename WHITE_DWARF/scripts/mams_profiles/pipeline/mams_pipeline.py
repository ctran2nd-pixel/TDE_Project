"""Prepare and launch isolated MESA middle-age main-sequence runs.

The pipeline reads ``tutorial/inlist_project`` as the authoritative physics
template, changes only the requested mass/composition/termination/output
fields, and creates one independent run directory per mass.  Each run keeps
its MESA profile, GYRE model, logs, and metadata separate from every other run.

Inputs
------
* A MESA work directory containing ``inlist``, ``inlist_project``,
  ``inlist_pgstar``, and a compiled ``star`` executable.
* One or more stellar masses and optional composition/evolution controls.

Outputs
-------
Run directories and a JSON manifest under ``results/mams_profiles`` by
default.  Setup is the default action; MESA is launched only with an explicit
``--run-local`` or ``--submit-slurm`` option.

Example
-------
Run from the repository root with Python 3.9 or newer::

    python3 \
      scripts/mams_profiles/pipeline/mams_pipeline.py \
      --masses 1.0 1.5 2.0
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import stat
import subprocess
from dataclasses import asdict, dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable, Sequence


@dataclass(frozen=True)
class MAMSConfig:
    """Physical and scheduler settings for one stellar model."""

    mass_msun: float
    initial_y: float = 0.2703
    initial_z: float = 0.0142
    central_h1_limit: float = 0.35
    partition: str = "typhon"
    walltime: str = "7-00:00:00"
    cpus_per_task: int = 1

    def validate(self) -> None:
        """Reject values that cannot describe the intended MESA run."""

        if self.mass_msun <= 0:
            raise ValueError("mass_msun must be positive")
        if not (0 < self.initial_y < 1):
            raise ValueError("initial_y must lie strictly between 0 and 1")
        if not (0 < self.initial_z < 1):
            raise ValueError("initial_z must lie strictly between 0 and 1")
        if self.initial_y + self.initial_z >= 1:
            raise ValueError("initial_y + initial_z must be less than 1")
        initial_h1 = 1 - self.initial_y - self.initial_z
        if not (0 < self.central_h1_limit < initial_h1):
            raise ValueError(
                "central_h1_limit must lie between zero and the initial "
                f"hydrogen fraction ({initial_h1:.6g})"
            )
        if self.cpus_per_task < 1:
            raise ValueError("cpus_per_task must be at least 1")


@dataclass(frozen=True)
class PreparedRun:
    """Paths and names associated with one prepared MESA run."""

    label: str
    run_dir: Path
    profile_path: Path
    gyre_path: Path
    metadata_path: Path


def find_project_root(start: Path | str = ".") -> Path:
    """Find the nearest parent containing the tutorial MESA work directory."""

    start_path = Path(start).expanduser().resolve()
    candidates = [start_path, *start_path.parents]
    for candidate in candidates:
        if (candidate / "tutorial" / "inlist_project").is_file():
            return candidate
    raise FileNotFoundError(
        f"Could not find tutorial/inlist_project above {start_path}"
    )


def mass_to_token(mass_msun: float | str | Decimal) -> str:
    """Convert a mass to the filename token used by the existing models.

    Examples are ``1 -> '1'``, ``1.5 -> '1p5'``, and
    ``0.15 -> '0p15'``.  Decimal conversion avoids binary floating-point
    artifacts in filenames.
    """

    try:
        value = Decimal(str(mass_msun))
    except InvalidOperation as exc:
        raise ValueError(f"Invalid stellar mass: {mass_msun!r}") from exc
    if not value.is_finite() or value <= 0:
        raise ValueError("stellar mass must be finite and positive")
    normalized = format(value.normalize(), "f")
    return normalized.replace("-", "m").replace(".", "p")


def model_label(mass_msun: float | str | Decimal) -> str:
    """Return the mass-specific basename, such as ``MAMS1p5Msun``."""

    return f"MAMS{mass_to_token(mass_msun)}Msun"


def _replace_assignment(
    text: str,
    key: str,
    value: str,
    *,
    expected_count: int = 1,
) -> str:
    """Replace one Fortran-namelist assignment while retaining comments."""

    pattern = re.compile(
        rf"(?m)^(?P<indent>[ \t]*){re.escape(key)}[ \t]*=[ \t]*"
        rf"(?P<old>[^!\n]*?)(?P<space>[ \t]*)(?P<comment>!.*)?$"
    )

    def replacement(match: re.Match[str]) -> str:
        comment = match.group("comment") or ""
        separator = "  " if comment else ""
        return f"{match.group('indent')}{key} = {value}{separator}{comment}"

    updated, count = pattern.subn(replacement, text)
    if count != expected_count:
        raise ValueError(
            f"Expected {expected_count} assignment(s) for {key!r}, found {count}"
        )
    return updated


def render_inlist_project(template_text: str, config: MAMSConfig) -> str:
    """Render the authoritative project inlist for a single mass."""

    config.validate()
    label = model_label(config.mass_msun)
    replacements = (
        ("filename_for_profile_when_terminate", f"'{label}.profile'"),
        ("save_pulse_data_filename", f"'{label}.mesa'"),
        ("initial_mass", f"{config.mass_msun:.12g}"),
        ("initial_y", f"{config.initial_y:.12g}"),
        ("initial_z", f"{config.initial_z:.12g}"),
        ("xa_central_lower_limit(1)", f"{config.central_h1_limit:.12g}"),
        ("Zbase", f"{config.initial_z:.12g}"),
    )
    rendered = template_text
    for key, value in replacements:
        rendered = _replace_assignment(rendered, key, value)
    return rendered


def _validate_template_dir(template_dir: Path) -> None:
    """Ensure the source MESA work directory has all required inputs."""

    required = ("inlist", "inlist_project", "inlist_pgstar", "star")
    missing = [name for name in required if not (template_dir / name).exists()]
    if missing:
        raise FileNotFoundError(
            f"Template directory {template_dir} is missing: {', '.join(missing)}"
        )
    if not os.access(template_dir / "star", os.X_OK):
        raise PermissionError(f"{template_dir / 'star'} is not executable")


def _write_executable(path: Path, text: str) -> None:
    """Write a small launch script and set its executable bits."""

    path.write_text(text, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _metadata(config: MAMSConfig, template_dir: Path) -> dict[str, object]:
    """Build the deterministic metadata used for idempotency checks."""

    return {
        "config": asdict(config),
        "template_dir": str(template_dir.resolve()),
        "label": model_label(config.mass_msun),
        "profile_filename": f"{model_label(config.mass_msun)}.profile",
        "gyre_filename": f"{model_label(config.mass_msun)}.mesa",
    }


def prepare_run(
    config: MAMSConfig,
    template_dir: Path | str,
    output_root: Path | str,
    *,
    overwrite: bool = False,
    copy_executable: bool = False,
) -> PreparedRun:
    """Create an isolated, directly runnable work directory for one mass.

    Repeating the call with identical settings is safe and returns the
    existing directory.  A directory with different settings is never
    removed unless ``overwrite=True`` is explicitly supplied.
    """

    config.validate()
    template = Path(template_dir).expanduser().resolve()
    output = Path(output_root).expanduser().resolve()
    _validate_template_dir(template)

    label = model_label(config.mass_msun)
    run_dir = output / label
    metadata = _metadata(config, template)
    metadata_path = run_dir / "metadata.json"

    if run_dir.exists() and not overwrite:
        if metadata_path.is_file():
            existing = json.loads(metadata_path.read_text(encoding="utf-8"))
            if existing == metadata:
                return PreparedRun(
                    label=label,
                    run_dir=run_dir,
                    profile_path=run_dir / f"{label}.profile",
                    gyre_path=run_dir / f"{label}.mesa",
                    metadata_path=metadata_path,
                )
        raise FileExistsError(
            f"{run_dir} already exists with different or incomplete settings; "
            "inspect it or rerun with overwrite=True"
        )
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)

    shutil.copy2(template / "inlist", run_dir / "inlist")
    shutil.copy2(template / "inlist_pgstar", run_dir / "inlist_pgstar")
    template_text = (template / "inlist_project").read_text(encoding="utf-8")
    rendered = render_inlist_project(template_text, config)
    (run_dir / "inlist_project").write_text(rendered, encoding="utf-8")

    star_target = run_dir / "star"
    if copy_executable:
        shutil.copy2(template / "star", star_target)
    else:
        relative_star = os.path.relpath(template / "star", start=run_dir)
        star_target.symlink_to(relative_star)

    _write_executable(
        run_dir / "run.sh",
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        'cd "$(dirname "$0")"\n'
        "./star > mesa.out 2> mesa.err\n",
    )
    _write_executable(
        run_dir / "submit.slurm",
        "#!/usr/bin/env bash\n"
        f"#SBATCH --job-name={label}\n"
        "#SBATCH --output=slurm-%j.out\n"
        "#SBATCH --error=slurm-%j.err\n"
        "#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task={config.cpus_per_task}\n"
        f"#SBATCH --time={config.walltime}\n"
        f"#SBATCH --partition={config.partition}\n\n"
        "set -euo pipefail\n"
        'cd "$(dirname "$0")"\n'
        f"export OMP_NUM_THREADS={config.cpus_per_task}\n"
        "./run.sh\n",
    )
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    return PreparedRun(
        label=label,
        run_dir=run_dir,
        profile_path=run_dir / f"{label}.profile",
        gyre_path=run_dir / f"{label}.mesa",
        metadata_path=metadata_path,
    )


def prepare_grid(
    configs: Iterable[MAMSConfig],
    template_dir: Path | str,
    output_root: Path | str,
    *,
    overwrite: bool = False,
    copy_executable: bool = False,
) -> list[PreparedRun]:
    """Prepare a mass grid and save a machine-readable grid manifest."""

    config_list = list(configs)
    labels = [model_label(config.mass_msun) for config in config_list]
    if len(set(labels)) != len(labels):
        raise ValueError("The mass grid contains duplicate output labels")

    runs = [
        prepare_run(
            config,
            template_dir,
            output_root,
            overwrite=overwrite,
            copy_executable=copy_executable,
        )
        for config in config_list
    ]
    output = Path(output_root).expanduser().resolve()
    manifest = {
        "template_dir": str(Path(template_dir).expanduser().resolve()),
        "runs": [
            {
                "label": run.label,
                "run_dir": str(run.run_dir),
                "profile_path": str(run.profile_path),
                "gyre_path": str(run.gyre_path),
            }
            for run in runs
        ],
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return runs


def run_local(run: PreparedRun, *, check: bool = True) -> subprocess.CompletedProcess[str]:
    """Run one prepared model locally and wait for it to finish."""

    return subprocess.run(
        ["./run.sh"], cwd=run.run_dir, check=check, text=True, capture_output=True
    )


def submit_slurm(run: PreparedRun) -> str:
    """Submit one prepared model with Slurm and return its job ID."""

    result = subprocess.run(
        ["sbatch", "--parsable", "submit.slurm"],
        cwd=run.run_dir,
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout.strip().split(";", maxsplit=1)[0]


def inspect_runs(runs: Iterable[PreparedRun]) -> list[dict[str, object]]:
    """Return a compact status table without parsing MESA's numerical data."""

    status = []
    for run in runs:
        status.append(
            {
                "label": run.label,
                "prepared": run.metadata_path.is_file(),
                "mesa_profile": run.profile_path.is_file(),
                "gyre_model": run.gyre_path.is_file(),
                "run_dir": str(run.run_dir),
            }
        )
    return status


def build_configs(
    masses: Sequence[float],
    *,
    initial_y: float = 0.2703,
    initial_z: float = 0.0142,
    central_h1_limit: float = 0.35,
    partition: str = "typhon",
    walltime: str = "7-00:00:00",
    cpus_per_task: int = 1,
) -> list[MAMSConfig]:
    """Construct and validate a consistently configured mass grid."""

    configs = [
        MAMSConfig(
            mass_msun=mass,
            initial_y=initial_y,
            initial_z=initial_z,
            central_h1_limit=central_h1_limit,
            partition=partition,
            walltime=walltime,
            cpus_per_task=cpus_per_task,
        )
        for mass in masses
    ]
    for config in configs:
        config.validate()
    return configs


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--masses", nargs="+", required=True, type=float)
    parser.add_argument("--template-dir", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--initial-y", type=float, default=0.2703)
    parser.add_argument("--initial-z", type=float, default=0.0142)
    parser.add_argument("--central-h1-limit", type=float, default=0.35)
    parser.add_argument("--partition", default="typhon")
    parser.add_argument("--walltime", default="7-00:00:00")
    parser.add_argument("--cpus-per-task", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--copy-executable", action="store_true")
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--run-local", action="store_true")
    action.add_argument("--submit-slurm", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point for setup, optional execution, and status reporting."""

    args = _parser().parse_args(argv)
    root = find_project_root(Path.cwd())
    template_dir = args.template_dir or root / "tutorial"
    output_root = args.output_root or root / "results" / "mams_profiles"
    configs = build_configs(
        args.masses,
        initial_y=args.initial_y,
        initial_z=args.initial_z,
        central_h1_limit=args.central_h1_limit,
        partition=args.partition,
        walltime=args.walltime,
        cpus_per_task=args.cpus_per_task,
    )
    runs = prepare_grid(
        configs,
        template_dir,
        output_root,
        overwrite=args.overwrite,
        copy_executable=args.copy_executable,
    )

    if args.run_local:
        for run in runs:
            print(f"Running {run.label} locally ...", flush=True)
            run_local(run)
    elif args.submit_slurm:
        for run in runs:
            print(f"Submitted {run.label}: job {submit_slurm(run)}")

    print(json.dumps(inspect_runs(runs), indent=2))
    print(f"Manifest: {Path(output_root).expanduser().resolve() / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
