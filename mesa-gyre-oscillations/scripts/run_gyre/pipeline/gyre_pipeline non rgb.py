"""Reusable stages for adiabatic oscillation calculations with GYRE.

Inputs
------
A MESA pulse-data file written with ``pulse_data_format = 'GYRE'`` and a
compiled GYRE executable.

Outputs
-------
The run directory contains a copied stellar model, generated ``gyre.in``,
solver log, native HDF5 summary/detail files, ``modes.csv``, and JSON metadata.

The command-line entry point is ``run_gyre.py`` in this directory. The module
is kept separate so the companion notebook can call each stage directly.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

import h5py
import numpy as np


FREQUENCY_UNITS = (
    "NONE",
    "HZ",
    "UHZ",
    "RAD_PER_SEC",
    "CYC_PER_DAY",
    "ACOUSTIC_DELTA",
    "GRAVITY_DELTA",
    "UPPER_DELTA",
    "LOWER_DELTA",
    "ACOUSTIC_CUTOFF",
    "GRAVITY_CUTOFF",
)


@dataclass(frozen=True)
class RunConfig:
    """Numerical and output choices used to render a GYRE input file.

    ``freq_units`` applies to the scan interval. ``output_freq_units`` only
    controls how discovered frequencies are reported. The distinction is
    important: the TDE workflow scans dimensionless angular frequency but
    reports frequencies in cycles per day.
    """

    harmonic_degrees: Tuple[int, ...] = (2,)
    n_pg_min: int = -20
    n_pg_max: int = 20
    freq_min: float = 0.1
    freq_max: float = 20.0
    n_freq: int = 2000
    freq_units: str = "NONE"
    output_freq_units: str = "CYC_PER_DAY"
    grid_type: str = "INVERSE"
    threads: int = 1
    w_osc: float = 10.0
    w_exp: float = 2.0
    w_ctr: float = 10.0
    n_iter_max: int = 150
    diff_scheme: str = "COLLOC_GL4"
    repair_As: bool = True
    use_n_pg_filter: bool = True

    def validated(self) -> "RunConfig":
        """Return this configuration after checking all user-facing values."""

        if not self.harmonic_degrees:
            raise ValueError("At least one harmonic degree is required.")
        if any(degree < 0 for degree in self.harmonic_degrees):
            raise ValueError("Harmonic degrees must be non-negative integers.")
        if len(set(self.harmonic_degrees)) != len(self.harmonic_degrees):
            raise ValueError("Harmonic degrees must not contain duplicates.")
        if self.n_pg_min > self.n_pg_max:
            raise ValueError("n_pg_min must be less than or equal to n_pg_max.")
        if self.freq_min <= 0 or self.freq_max <= self.freq_min:
            raise ValueError("Require 0 < freq_min < freq_max.")
        if self.n_freq < 2:
            raise ValueError("n_freq must be at least 2.")
        if self.freq_units.upper() not in FREQUENCY_UNITS:
            raise ValueError(f"Unsupported scan frequency units: {self.freq_units}")
        if self.output_freq_units.upper() not in FREQUENCY_UNITS:
            raise ValueError(
                f"Unsupported output frequency units: {self.output_freq_units}"
            )
        if self.grid_type.upper() not in {"LINEAR", "INVERSE"}:
            raise ValueError("grid_type must be LINEAR or INVERSE.")
        if self.threads < 1:
            raise ValueError("threads must be at least 1.")
        return self


@dataclass(frozen=True)
class RunResult:
    """Paths and status returned by :func:`run_pipeline`."""

    output_dir: Path
    input_model: Path
    input_file: Path
    metadata_file: Path
    log_file: Optional[Path]
    summary_file: Optional[Path]
    modes_csv: Optional[Path]
    dry_run: bool


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""

    return datetime.now(timezone.utc).isoformat()


def file_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Calculate a SHA-256 checksum without loading a model fully into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def validate_model_file(model_path: Path | str) -> Dict[str, Any]:
    """Perform inexpensive checks for a MESA GYRE-format pulse file.

    MESA has several versions of this text format, so the check deliberately
    avoids enforcing one column count. It verifies that the initial non-empty
    line is numeric, begins with a plausible number of mesh points, and that
    the file contains more rows than its header alone.
    """

    path = Path(model_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(
            f"Stellar model not found: {path}. Export pulse data from MESA "
            "with pulse_data_format = 'GYRE'."
        )
    if not path.is_file():
        raise ValueError(f"Stellar model path is not a regular file: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Stellar model is empty: {path}")

    first_line = ""
    nonempty_lines = 0
    with path.open("r", encoding="utf-8", errors="strict") as stream:
        for line in stream:
            if line.strip():
                nonempty_lines += 1
                if not first_line:
                    first_line = line.strip()

    fields = first_line.split()
    if len(fields) < 2:
        raise ValueError(
            "The first non-empty line does not resemble a MESA pulse-data header."
        )
    try:
        mesh_points = int(fields[0])
        numeric_header = [
            float(value.replace("D", "E").replace("d", "e"))
            for value in fields[1:]
        ]
    except ValueError as exc:
        raise ValueError(
            "The first non-empty line is not a numeric MESA pulse-data header. "
            "A standard MESA profile*.data file is not a GYRE input model."
        ) from exc
    if mesh_points < 2:
        raise ValueError("The pulse-data header reports fewer than two mesh points.")
    if nonempty_lines <= 1:
        raise ValueError("The pulse-data file contains a header but no model rows.")

    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "mesh_points_header": mesh_points,
        "header_numeric_fields": numeric_header,
        "nonempty_lines": nonempty_lines,
        "sha256": file_sha256(path),
    }


def resolve_gyre_executable(explicit_path: Optional[Path | str] = None) -> Path:
    """Resolve an executable from a CLI path, ``GYRE_DIR``, or ``PATH``."""

    candidates = []
    if explicit_path is not None:
        candidates.append(Path(explicit_path).expanduser())

    gyre_dir = os.environ.get("GYRE_DIR")
    if gyre_dir:
        candidates.append(Path(gyre_dir).expanduser() / "bin" / "gyre")

    on_path = shutil.which("gyre")
    if on_path:
        candidates.append(Path(on_path))

    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.is_file() and os.access(str(resolved), os.X_OK):
            return resolved

    searched = ", ".join(str(path) for path in candidates) or "no candidate paths"
    raise FileNotFoundError(
        "Could not find an executable GYRE frontend. Pass --gyre-executable, "
        f"set GYRE_DIR, or add gyre to PATH. Searched: {searched}"
    )


def render_gyre_input(config: RunConfig, model_filename: str = "input_model.mesa") -> str:
    """Render a complete GYRE v7-compatible adiabatic namelist."""

    cfg = config.validated()
    if "'" in model_filename or "\n" in model_filename:
        raise ValueError("model_filename cannot contain a quote or newline.")

    mode_groups = []

    for degree in cfg.harmonic_degrees:
        lines = [
            "&mode",
            f"  l = {degree}",
        ]
    
        if cfg.use_n_pg_filter:
            lines.extend([
                f"  n_pg_min = {cfg.n_pg_min}",
                f"  n_pg_max = {cfg.n_pg_max}",
            ])
    
        lines.append("/")
    
        mode_groups.append("\n".join(lines))

    return f"""&constants
/

&model
  model_type = 'EVOL'
  file = '{model_filename}'
  file_format = 'MESA'
  repair_As = {str(cfg.repair_As).lower()}
/

{chr(10).join(mode_groups)}

&osc
  outer_bound = 'VACUUM'
/

&rot
/

&num
  diff_scheme = '{cfg.diff_scheme.upper()}'
  n_iter_max = {cfg.n_iter_max}
/

&scan
  grid_type = '{cfg.grid_type.upper()}'
  freq_min = {cfg.freq_min:.16g}
  freq_max = {cfg.freq_max:.16g}
  n_freq = {cfg.n_freq}
  freq_units = '{cfg.freq_units.upper()}'
/

&grid
  w_osc = {cfg.w_osc:.16g}
  w_exp = {cfg.w_exp:.16g}
  w_ctr = {cfg.w_ctr:.16g}
/

&ad_output
  summary_file = 'summary.h5'
  summary_item_list = 'l,n_pg,n_p,n_g,omega,freq,E_norm,R_star,M_star'
  detail_template = 'detail.l%l.n%n.h5'
  detail_item_list = 'l,n_pg,omega,x,xi_r,xi_h,c_1,As,V_2,Gamma_1,rho,M_r'
  freq_units = '{cfg.output_freq_units.upper()}'
/

&nad_output
/
"""


def _json_value(value: Any) -> Any:
    """Convert paths, NumPy scalars, bytes, and complex values for JSON."""

    if isinstance(value, Path):
        return str(value)
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _write_metadata(path: Path, metadata: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(_json_value(metadata), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _complex_array(values: np.ndarray) -> np.ndarray:
    """Convert GYRE's HDF5 ``re``/``im`` compound values to complex arrays."""

    names = values.dtype.names
    if names and "re" in names and "im" in names:
        return values["re"] + 1j * values["im"]
    return values


def load_gyre_hdf5(path: Path | str) -> Dict[str, Any]:
    """Load all root datasets and attributes from a GYRE HDF5 output file."""

    output_path = Path(path)
    if not output_path.is_file():
        raise FileNotFoundError(f"GYRE output not found: {output_path}")

    data: Dict[str, Any] = {}
    with h5py.File(output_path, "r") as handle:
        for name, dataset in handle.items():
            data[name] = _complex_array(np.asarray(dataset))
        data["_attrs"] = {
            name: _json_value(value) for name, value in handle.attrs.items()
        }
    return data


def load_summary(path: Path | str) -> Dict[str, Any]:
    """Load a GYRE summary file and verify that its datasets align by mode."""

    data = load_gyre_hdf5(path)
    lengths = {
        len(value)
        for key, value in data.items()
        if key != "_attrs" and isinstance(value, np.ndarray) and value.ndim > 0
    }
    if len(lengths) > 1:
        raise ValueError(f"Summary datasets have inconsistent lengths: {sorted(lengths)}")
    return data


def export_modes_csv(summary_path: Path | str, csv_path: Path | str) -> Path:
    """Export one CSV row per mode, splitting complex values into two columns."""

    summary = load_summary(summary_path)
    arrays = {
        name: values
        for name, values in summary.items()
        if name != "_attrs" and isinstance(values, np.ndarray) and values.ndim > 0
    }
    if not arrays:
        raise ValueError("GYRE summary contains no mode datasets.")

    preferred = [
        "l",
        "n_pg",
        "n_p",
        "n_g",
        "omega",
        "freq",
        "E_norm",
        "R_star",
        "M_star",
    ]
    names = [name for name in preferred if name in arrays]
    names.extend(sorted(set(arrays) - set(names)))

    columns = []
    for name in names:
        if np.iscomplexobj(arrays[name]):
            columns.extend((f"{name}_real", f"{name}_imag"))
        else:
            columns.append(name)

    destination = Path(csv_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    row_count = len(next(iter(arrays.values())))
    with destination.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for index in range(row_count):
            row: Dict[str, Any] = {}
            for name in names:
                value = arrays[name][index]
                if np.iscomplexobj(arrays[name]):
                    row[f"{name}_real"] = float(np.real(value))
                    row[f"{name}_imag"] = float(np.imag(value))
                else:
                    row[name] = _json_value(value)
            writer.writerow(row)
    return destination


def _prepare_output_directory(
    model_path: Path, output_dir: Path, overwrite: bool
) -> None:
    """Create an empty run directory without risking deletion of the input."""

    resolved_output = output_dir.expanduser().resolve()
    protected_paths = {
        Path(resolved_output.anchor),
        Path.home().resolve(),
        Path.cwd().resolve(),
    }
    if resolved_output in protected_paths:
        raise ValueError(
            f"Refusing to use a protected directory as output: {resolved_output}"
        )
    try:
        model_path.relative_to(resolved_output)
    except ValueError:
        pass
    else:
        raise ValueError("The input model cannot be stored inside the output directory.")

    if resolved_output.exists():
        if not overwrite:
            raise FileExistsError(
                f"Output directory already exists: {resolved_output}. "
                "Choose another path or pass --overwrite."
            )
        if not resolved_output.is_dir():
            raise ValueError(f"Output path exists and is not a directory: {resolved_output}")
        entries = list(resolved_output.iterdir())
        marker = resolved_output / "run_metadata.json"
        if entries and not marker.is_file():
            raise ValueError(
                "Refusing to overwrite a non-empty directory that was not created "
                f"by this pipeline: {resolved_output}"
            )
        shutil.rmtree(resolved_output)
    resolved_output.mkdir(parents=True)


def _extract_gyre_version(log_path: Path) -> Optional[str]:
    """Read the version token printed at the beginning of a GYRE log."""

    pattern = re.compile(r"^gyre \[([^]]+)\]")
    with log_path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = pattern.match(line.strip())
            if match:
                return match.group(1)
    return None


def run_pipeline(
    model_path: Path | str,
    output_dir: Path | str,
    config: RunConfig = RunConfig(),
    *,
    gyre_executable: Optional[Path | str] = None,
    overwrite: bool = False,
    dry_run: bool = False,
) -> RunResult:
    """Prepare, optionally execute, validate, and summarize one GYRE run."""

    cfg = config.validated()
    model_info = validate_model_file(model_path)
    source_model = Path(model_info["path"])
    run_dir = Path(output_dir).expanduser().resolve()
    _prepare_output_directory(source_model, run_dir, overwrite)

    copied_model = run_dir / "input_model.mesa"
    shutil.copy2(source_model, copied_model)
    input_file = run_dir / "gyre.in"
    input_file.write_text(render_gyre_input(cfg, copied_model.name), encoding="utf-8")
    metadata_file = run_dir / "run_metadata.json"

    metadata: Dict[str, Any] = {
        "created_at": utc_now(),
        "status": "prepared" if dry_run else "running",
        "dry_run": dry_run,
        "source_model": model_info,
        "copied_model": str(copied_model),
        "gyre_input": str(input_file),
        "config": asdict(cfg),
        "gyre_executable": str(gyre_executable) if gyre_executable else None,
    }
    _write_metadata(metadata_file, metadata)

    if dry_run:
        return RunResult(
            output_dir=run_dir,
            input_model=copied_model,
            input_file=input_file,
            metadata_file=metadata_file,
            log_file=None,
            summary_file=None,
            modes_csv=None,
            dry_run=True,
        )

    executable = resolve_gyre_executable(gyre_executable)
    metadata["gyre_executable"] = str(executable)
    log_file = run_dir / "gyre.log"
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(cfg.threads)

    try:
        with log_file.open("w", encoding="utf-8") as log_stream:
            completed = subprocess.run(
                [str(executable), input_file.name],
                cwd=run_dir,
                env=environment,
                stdout=log_stream,
                stderr=subprocess.STDOUT,
                check=False,
            )
        metadata["return_code"] = completed.returncode
        metadata["gyre_version"] = _extract_gyre_version(log_file)
        if completed.returncode != 0:
            metadata["status"] = "failed"
            raise RuntimeError(
                f"GYRE exited with code {completed.returncode}. See {log_file}."
            )

        summary_file = run_dir / "summary.h5"
        if not summary_file.is_file() or summary_file.stat().st_size == 0:
            metadata["status"] = "failed"
            raise RuntimeError(f"GYRE completed without a non-empty {summary_file.name}.")

        summary = load_summary(summary_file)
        mode_count = len(summary.get("l", []))
        modes_csv = export_modes_csv(summary_file, run_dir / "modes.csv")
        details = sorted(path.name for path in run_dir.glob("detail.l*.n*.h5"))
        metadata.update(
            {
                "status": "complete",
                "completed_at": utc_now(),
                "mode_count": mode_count,
                "detail_file_count": len(details),
                "summary_file": str(summary_file),
                "modes_csv": str(modes_csv),
                "detail_files": details,
            }
        )
    except Exception as exc:
        metadata.setdefault("completed_at", utc_now())
        metadata["status"] = "failed"
        metadata["error"] = str(exc)
        _write_metadata(metadata_file, metadata)
        raise

    _write_metadata(metadata_file, metadata)
    return RunResult(
        output_dir=run_dir,
        input_model=copied_model,
        input_file=input_file,
        metadata_file=metadata_file,
        log_file=log_file,
        summary_file=summary_file,
        modes_csv=modes_csv,
        dry_run=False,
    )