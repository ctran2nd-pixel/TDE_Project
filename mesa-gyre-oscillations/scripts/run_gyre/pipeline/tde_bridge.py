"""Validate and stage GYRE HDF5 output for the collaborator's TDECalculator.

This module is deliberately an I/O bridge, not a second implementation of
the tidal-disruption calculation.  The collaborator's ``TDECalculator``
expects a directory whose summary contains *only* ``l=2`` modes and whose
detail files are named ``detail.l2.n+0.h5``, ``detail.l2.n-1.h5``, and so on.
The clean GYRE workflow writes collision-safe, ID-bearing filenames instead
(``detail.id%id.l%l.n%n.h5``).  :func:`prepare_tde_directory` validates the
scientific handoff and creates a fresh compatibility directory without
changing the GYRE source run.

Inputs
------
``source_dir`` is a GYRE run directory containing ``summary.h5`` and detail
files.  The selected summary rows must provide ``l``, ``n_pg``, ``omega``,
``M_star``, and ``R_star``.  Every selected detail file must provide ``x``,
``rho``, ``M_r``, ``xi_r``, and ``xi_h``.

Outputs
-------
The staged directory contains a filtered ``summary.h5``, legacy-named detail
files, and ``tde_bridge_manifest.json`` with source checksums and mappings.
The default transfer mode is a byte-for-byte copy, making the result portable.

Example
-------
From the repository root::

    python scripts/run_gyre/pipeline/prepare_tde.py prepare \
        results/run_gyre/my_star --output-dir results/tde_stage/my_star

The bridge validates the file contract only.  In particular it does not
validate whether the collaborator's mass-only M-gamma calibration is suitable
for a non-main-sequence stellar structure.
"""

from __future__ import annotations

import contextlib
import hashlib
import importlib.util
import json
import os
import platform
import re
import shutil
import sys
import uuid
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterator, List, Mapping, Optional, Sequence, Tuple, Type

import h5py
import numpy as np


BRIDGE_VERSION = "1.0"
SUMMARY_FILENAME = "summary.h5"
MANIFEST_FILENAME = "tde_bridge_manifest.json"
REQUIRED_SUMMARY_FIELDS = ("l", "n_pg", "omega", "M_star", "R_star")
REQUIRED_DETAIL_FIELDS = ("x", "rho", "M_r", "xi_r", "xi_h")
TRANSFER_MODES = ("copy", "hardlink", "symlink")
SOLAR_MASS_GRAMS = 1.9884e33  # Exact conversion used by collaborator TDECalculator.py.

_ID_DETAIL_RE = re.compile(
    r"^detail\.id(?P<id>[^.]+)\.l(?P<l>[+-]?\d+)\.n(?P<n>[+-]\d+)\.h5$"
)
_LEGACY_DETAIL_RE = re.compile(
    r"^detail\.l(?P<l>[+-]?\d+)\.n(?P<n>[+-]\d+)\.h5$"
)


@dataclass(frozen=True)
class ModeRecord:
    """One validated quadrupole mode selected for the TDE handoff.

    ``summary_index`` indexes the original, unfiltered GYRE summary.  The
    staged summary preserves the order of these records because the external
    calculator pairs every summary row with its legacy-named detail file.
    """

    summary_index: int
    n_pg: int
    omega_real: float
    omega_imag: float
    source_detail: Path
    source_style: str
    mode_id: Optional[str] = None


@dataclass(frozen=True)
class GyreRunInspection:
    """Validated source-run metadata returned by :func:`inspect_gyre_run`."""

    source_dir: Path
    summary_path: Path
    modes: Tuple[ModeRecord, ...]
    summary_mode_count: int
    l2_mode_count: int
    n_pg_min: Optional[int]
    n_pg_max: Optional[int]
    f_mode_index: int
    detail_grid_size: int
    detail_grid_reference: Path

    @property
    def selected_n_pg(self) -> Tuple[int, ...]:
        """Radial orders in the exact order passed to TDECalculator."""

        return tuple(mode.n_pg for mode in self.modes)


@dataclass(frozen=True)
class PreparedTDERun:
    """Paths and provenance of a completed compatibility staging operation."""

    directory: Path
    summary_path: Path
    manifest_path: Path
    inspection: GyreRunInspection
    transfer_mode: str


@dataclass(frozen=True)
class TDEProject:
    """Validated paths required by the unmodified collaborator calculator."""

    root: Path
    calculator_path: Path
    gamma_table_path: Path
    gamma_mass_domain_solar: Tuple[float, float]
    calculator_sha256: str
    gamma_table_sha256: str


@dataclass(frozen=True)
class SpectrumPlotReport:
    """Accounting for points shown and explicitly masked by the diagnostic."""

    plotted_modes: int
    masked_nonfinite_modes: int


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Return a SHA-256 digest without reading a whole HDF5 file into memory."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _json_value(value: Any) -> Any:
    """Convert paths and NumPy/HDF5 scalar values to JSON-safe values."""

    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _as_real_imag(values: np.ndarray, field_name: str) -> Tuple[np.ndarray, np.ndarray]:
    """Return real and imaginary arrays from native or compound GYRE values."""

    array = np.asarray(values)
    if array.dtype.names:
        names = set(array.dtype.names)
        if {"re", "im"}.issubset(names):
            return np.asarray(array["re"], dtype=float), np.asarray(array["im"], dtype=float)
        if {"real", "imag"}.issubset(names):
            return (
                np.asarray(array["real"], dtype=float),
                np.asarray(array["imag"], dtype=float),
            )
        raise ValueError(
            f"{field_name} has a compound dtype without re/im fields: {array.dtype}"
        )
    if not np.issubdtype(array.dtype, np.number):
        raise ValueError(f"{field_name} must be numeric, got {array.dtype}.")
    complex_array = np.asarray(array, dtype=np.complex128)
    return complex_array.real, complex_array.imag


def _require_vector(
    h5: h5py.File,
    name: str,
    expected_length: Optional[int] = None,
) -> np.ndarray:
    """Load a required one-dimensional root dataset with a clear error."""

    if name not in h5:
        raise ValueError(f"{h5.filename} is missing required dataset {name!r}.")
    item = h5[name]
    if not isinstance(item, h5py.Dataset):
        raise ValueError(f"{h5.filename}:{name} must be a dataset.")
    values = item[...]
    if values.ndim != 1:
        raise ValueError(
            f"{h5.filename}:{name} must be one-dimensional; got shape {values.shape}."
        )
    if expected_length is not None and values.shape[0] != expected_length:
        raise ValueError(
            f"{h5.filename}:{name} has length {values.shape[0]}, expected {expected_length}."
        )
    return values


def _integer_vector(values: np.ndarray, field_name: str) -> np.ndarray:
    """Validate and return an integer-valued vector without silently rounding."""

    real, imag = _as_real_imag(values, field_name)
    if not np.all(np.isfinite(real)) or not np.all(np.isfinite(imag)):
        raise ValueError(f"{field_name} contains non-finite values.")
    if not np.all(imag == 0) or not np.all(real == np.rint(real)):
        raise ValueError(f"{field_name} must contain integer-valued real entries.")
    return np.asarray(np.rint(real), dtype=np.int64)


def _positive_real_vector(values: np.ndarray, field_name: str) -> np.ndarray:
    """Require finite, positive, real values used directly by TDECalculator."""

    real, imag = _as_real_imag(values, field_name)
    if not np.all(np.isfinite(real)) or not np.all(np.isfinite(imag)):
        raise ValueError(f"{field_name} contains non-finite values.")
    if not np.all(imag == 0):
        raise ValueError(f"{field_name} must be real for the collaborator TDECalculator.")
    if not np.all(real > 0):
        raise ValueError(f"{field_name} must be strictly positive for the TDE calculation.")
    return real


def _validate_consistent_stellar_property(values: np.ndarray, name: str) -> None:
    """Reject summary rows that would make TDECalculator depend on row order."""

    if not np.allclose(values, values[0], rtol=1.0e-12, atol=0.0):
        raise ValueError(
            f"Selected l=2 summary rows have inconsistent {name}; TDECalculator uses "
            "the first row's stellar property. Re-run or repair the GYRE output."
        )


def _parse_detail_filename(path: Path) -> Optional[Dict[str, Any]]:
    """Parse either modern ID-bearing or legacy detail filename conventions."""

    match = _ID_DETAIL_RE.match(path.name)
    if match:
        return {
            "style": "id-bearing",
            "l": int(match.group("l")),
            "n_pg": int(match.group("n")),
            "id": match.group("id"),
        }
    match = _LEGACY_DETAIL_RE.match(path.name)
    if match:
        return {
            "style": "legacy",
            "l": int(match.group("l")),
            "n_pg": int(match.group("n")),
            "id": None,
        }
    return None


def _summary_mode_id(values: Optional[np.ndarray], index: int) -> Optional[str]:
    """Turn an optional summary ID value into a filename-comparable string."""

    if values is None:
        return None
    value = values[index]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _metadata_mode_id(value: Any) -> str:
    """Normalize a scalar detail ID using the same convention as summary IDs."""

    scalar = np.asarray(value)
    if scalar.size != 1:
        raise ValueError(f"Detail id metadata must be scalar; got shape {scalar.shape}.")
    return _summary_mode_id(np.asarray([scalar.reshape(-1)[0]]), 0) or ""


def _detail_candidates(source_dir: Path) -> List[Tuple[Path, Dict[str, Any]]]:
    """List only filenames that match a documented TDE/GYRE detail pattern."""

    candidates: List[Tuple[Path, Dict[str, Any]]] = []
    for path in sorted(source_dir.glob("detail*.h5")):
        if not path.is_file():
            continue
        parsed = _parse_detail_filename(path)
        if parsed is not None:
            candidates.append((path.resolve(), parsed))
    return candidates


def _locate_detail(
    candidates: Sequence[Tuple[Path, Dict[str, Any]]],
    n_pg: int,
    mode_id: Optional[str],
    source_dir: Path,
) -> Tuple[Path, str]:
    """Find exactly one source detail, preferring an exact ID-bearing match."""

    matching = [
        (path, metadata)
        for path, metadata in candidates
        if metadata["l"] == 2 and metadata["n_pg"] == n_pg
    ]
    id_matches = [
        (path, metadata)
        for path, metadata in matching
        if metadata["style"] == "id-bearing"
        and mode_id is not None
        and metadata["id"] == mode_id
    ]
    if len(id_matches) == 1:
        return id_matches[0][0], "id-bearing"
    if len(id_matches) > 1:
        raise ValueError(
            f"Ambiguous ID-bearing detail files for l=2, n_pg={n_pg}, id={mode_id} "
            f"in {source_dir}: {[str(path) for path, _ in id_matches]}"
        )

    legacy_matches = [
        (path, metadata) for path, metadata in matching if metadata["style"] == "legacy"
    ]
    if len(legacy_matches) == 1:
        return legacy_matches[0][0], "legacy"
    if len(legacy_matches) > 1:
        raise ValueError(
            f"Ambiguous legacy detail files for l=2, n_pg={n_pg} in {source_dir}: "
            f"{[str(path) for path, _ in legacy_matches]}"
        )
    if len(matching) == 1:
        path, metadata = matching[0]
        if metadata["style"] == "id-bearing" and mode_id is not None:
            raise ValueError(
                f"ID-bearing detail file {path.name} has id={metadata['id']!r}, but "
                f"the summary row for l=2, n_pg={n_pg} has id={mode_id!r}. Refusing "
                "to pair a possibly different eigenfunction."
            )
        return path, str(metadata["style"])
    if len(matching) > 1:
        raise ValueError(
            f"Cannot identify one detail file for l=2, n_pg={n_pg} in {source_dir}. "
            "The summary has no matching id, or several ID-bearing files share this order."
        )
    raise FileNotFoundError(
        f"Missing detail file for selected l=2, n_pg={n_pg} in {source_dir}. "
        "Expected detail.id*.l2.n±N.h5 or detail.l2.n±N.h5."
    )


def _detail_scalar_metadata(detail: h5py.File, name: str, *, required: bool) -> Optional[np.ndarray]:
    """Read one metadata value from an attribute or scalar root dataset."""

    if name in detail.attrs:
        value = np.asarray(detail.attrs[name])
    elif name in detail:
        item = detail[name]
        if not isinstance(item, h5py.Dataset) or item.shape != ():
            raise ValueError(
                f"{detail.filename}:{name} metadata must be an attribute or scalar dataset."
            )
        value = np.asarray(item[()])
    elif required:
        raise ValueError(
            f"{detail.filename} is missing required {name!r} metadata (attribute or scalar dataset)."
        )
    else:
        return None
    if value.size != 1:
        raise ValueError(
            f"{detail.filename}:{name} metadata must be scalar; got shape {value.shape}."
        )
    return value


def _validate_detail_file(
    path: Path,
    expected_n_pg: int,
    expected_omega_real: float,
    expected_omega_imag: float,
    expected_mode_id: Optional[str],
    source_style: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Validate fields and physical grid properties needed by TDECalculator.

    The external calculator uses the f-mode's radial grid and density profile
    for every mode, then adds each mode's eigenfunctions into arrays allocated
    on that grid.  Thus all staged detail files must have equally sized,
    mutually consistent ``x`` grids.  The reference grid is returned to the
    caller for that cross-mode check.
    """

    try:
        with h5py.File(path, "r") as detail:
            x_values = _require_vector(detail, "x")
            size = x_values.shape[0]
            if size < 2:
                raise ValueError(f"{path}: detail grid must contain at least two points.")
            arrays = {"x": x_values}
            arrays.update(
                {
                    name: _require_vector(detail, name, expected_length=size)
                    for name in REQUIRED_DETAIL_FIELDS
                    if name != "x"
                }
            )

            for name, values in arrays.items():
                real, imag = _as_real_imag(values, f"{path}:{name}")
                if not np.all(np.isfinite(real)) or not np.all(np.isfinite(imag)):
                    raise ValueError(f"{path}:{name} contains non-finite values.")

            x, x_imag = _as_real_imag(arrays["x"], f"{path}:x")
            rho, rho_imag = _as_real_imag(arrays["rho"], f"{path}:rho")
            m_r, m_r_imag = _as_real_imag(arrays["M_r"], f"{path}:M_r")
            if not np.all(x_imag == 0) or not np.all(rho_imag == 0) or not np.all(m_r_imag == 0):
                raise ValueError(f"{path}: x, rho, and M_r must be real-valued.")
            if np.any(x < 0) or np.any(x > 1.0 + 1.0e-8):
                raise ValueError(f"{path}: x must be a normalized radius in [0, 1].")
            if not np.all(np.diff(x) > 0):
                raise ValueError(f"{path}: x must be strictly increasing.")
            if np.any(rho < 0):
                raise ValueError(f"{path}: rho must be non-negative.")
            if np.any(m_r < 0) or np.any(np.diff(m_r) < -1.0e-12 * max(1.0, float(np.max(m_r)))):
                raise ValueError(f"{path}: M_r must be non-negative and non-decreasing.")
            for eigenfunction in ("xi_r", "xi_h"):
                _, imaginary = _as_real_imag(arrays[eigenfunction], f"{path}:{eigenfunction}")
                if not np.all(imaginary == 0):
                    raise ValueError(
                        f"{path}:{eigenfunction} must be real for this adiabatic handoff; "
                        "TDECalculator would otherwise silently discard its imaginary part."
                    )

            for name, expected in (("l", 2), ("n_pg", expected_n_pg)):
                metadata = _detail_scalar_metadata(detail, name, required=True)
                actual = _integer_vector(metadata.reshape(1), f"{path}:{name} metadata")
                if int(actual[0]) != expected:
                    raise ValueError(
                        f"{path}: {name} metadata is {int(actual[0])}, expected {expected}."
                    )
            omega_metadata = _detail_scalar_metadata(detail, "omega", required=True)
            omega_real, omega_imag = _as_real_imag(omega_metadata, f"{path}:omega metadata")
            actual_omega_real = float(np.asarray(omega_real).reshape(-1)[0])
            actual_omega_imag = float(np.asarray(omega_imag).reshape(-1)[0])
            if not np.isfinite(actual_omega_real) or not np.isfinite(actual_omega_imag):
                raise ValueError(f"{path}:omega metadata contains non-finite values.")
            if not np.isclose(actual_omega_real, expected_omega_real, rtol=1.0e-10, atol=1.0e-12) or not np.isclose(
                actual_omega_imag, expected_omega_imag, rtol=1.0e-10, atol=1.0e-12
            ):
                raise ValueError(
                    f"{path}:omega metadata ({actual_omega_real}+{actual_omega_imag}j) does not "
                    f"match summary omega ({expected_omega_real}+{expected_omega_imag}j)."
                )
            if source_style == "id-bearing" and expected_mode_id is not None:
                detail_id = _detail_scalar_metadata(detail, "id", required=False)
                if detail_id is not None and _metadata_mode_id(detail_id) != expected_mode_id:
                    raise ValueError(
                        f"{path}:id metadata {_metadata_mode_id(detail_id)!r} does not match "
                        f"summary id {expected_mode_id!r}."
                    )
            return (
                np.asarray(x, dtype=float),
                np.asarray(rho, dtype=float),
                np.asarray(m_r, dtype=float),
                size,
            )
    except OSError as exc:
        raise ValueError(f"Could not read GYRE detail HDF5 file {path}: {exc}") from exc


def _validate_radial_order_bounds(
    n_pg_min: Optional[int], n_pg_max: Optional[int]
) -> Tuple[Optional[int], Optional[int]]:
    """Validate optional inclusive n_pg filters without coercing user input."""

    for name, value in (("n_pg_min", n_pg_min), ("n_pg_max", n_pg_max)):
        if value is not None and (not isinstance(value, (int, np.integer))):
            raise ValueError(f"{name} must be an integer or None.")
    if n_pg_min is not None and n_pg_max is not None and n_pg_min > n_pg_max:
        raise ValueError("n_pg_min must be less than or equal to n_pg_max.")
    return n_pg_min, n_pg_max


def inspect_gyre_run(
    source_dir: Path | str,
    *,
    n_pg_min: Optional[int] = None,
    n_pg_max: Optional[int] = None,
) -> GyreRunInspection:
    """Preflight a GYRE run for the unmodified quadrupole TDECalculator.

    Only ``l=2`` rows are selected.  ``n_pg_min`` and ``n_pg_max`` are
    explicit, inclusive exploratory filters, recorded in the final manifest;
    they never silently truncate a mode set.  The resulting selection must
    contain exactly one f mode (``n_pg=0``), unique radial orders, finite
    nonzero real angular frequencies, and matching radial grids.

    ``summary['freq']`` is intentionally not read: the external calculator
    uses the dimensionless ``summary['omega']`` field, not display frequency
    units.
    """

    n_pg_min, n_pg_max = _validate_radial_order_bounds(n_pg_min, n_pg_max)
    directory = Path(source_dir).expanduser().resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"GYRE run directory does not exist: {directory}")
    summary_path = directory / SUMMARY_FILENAME
    if not summary_path.is_file():
        raise FileNotFoundError(f"GYRE summary not found: {summary_path}")

    try:
        with h5py.File(summary_path, "r") as summary:
            missing = [name for name in REQUIRED_SUMMARY_FIELDS if name not in summary]
            if missing:
                raise ValueError(
                    f"{summary_path} is missing required summary dataset(s): {', '.join(missing)}."
                )
            l_values = _require_vector(summary, "l")
            mode_count = l_values.shape[0]
            if mode_count == 0:
                raise ValueError(f"{summary_path} has no mode rows.")
            n_pg_values = _require_vector(summary, "n_pg", mode_count)
            omega_values = _require_vector(summary, "omega", mode_count)
            m_star_values = _require_vector(summary, "M_star", mode_count)
            r_star_values = _require_vector(summary, "R_star", mode_count)
            id_values = _require_vector(summary, "id", mode_count) if "id" in summary else None

            l_int = _integer_vector(l_values, f"{summary_path}:l")
            n_pg_int = _integer_vector(n_pg_values, f"{summary_path}:n_pg")
            omega_real, omega_imag = _as_real_imag(omega_values, f"{summary_path}:omega")
    except OSError as exc:
        raise ValueError(f"Could not read GYRE summary HDF5 file {summary_path}: {exc}") from exc

    l2_indices = np.flatnonzero(l_int == 2)
    if l2_indices.size == 0:
        raise ValueError(f"{summary_path} has no l=2 quadrupole modes.")
    selected_indices = l2_indices
    if n_pg_min is not None:
        selected_indices = selected_indices[n_pg_int[selected_indices] >= n_pg_min]
    if n_pg_max is not None:
        selected_indices = selected_indices[n_pg_int[selected_indices] <= n_pg_max]
    if selected_indices.size == 0:
        raise ValueError("The requested n_pg filter selects no l=2 modes.")

    selected_n_pg = n_pg_int[selected_indices]
    duplicates = sorted({int(n) for n in selected_n_pg if np.count_nonzero(selected_n_pg == n) > 1})
    if duplicates:
        raise ValueError(
            "Selected l=2 summary contains duplicate n_pg values, so the legacy "
            f"TDE filenames would be ambiguous: {duplicates}."
        )
    f_matches = np.flatnonzero(selected_n_pg == 0)
    if f_matches.size != 1:
        raise ValueError(
            "The selected l=2 handoff must contain exactly one f mode with n_pg=0. "
            "Adjust the explicit n_pg bounds or rerun GYRE."
        )
    selected_omega_real = omega_real[selected_indices]
    selected_omega_imag = omega_imag[selected_indices]
    if not np.all(np.isfinite(selected_omega_real)) or not np.all(
        np.isfinite(selected_omega_imag)
    ):
        raise ValueError(f"{summary_path}: selected l=2 omega contains non-finite values.")
    if not np.all(selected_omega_imag == 0):
        raise ValueError(
            f"{summary_path}: selected l=2 omega must be real for this adiabatic TDE "
            "handoff; TDECalculator would otherwise silently discard its imaginary part."
        )
    if np.any(selected_omega_real <= 0):
        raise ValueError(
            f"{summary_path}: selected l=2 omega must be positive; TDECalculator divides "
            "by its real part."
        )
    m_star = _positive_real_vector(
        m_star_values[selected_indices], f"{summary_path}: selected l=2 M_star"
    )
    r_star = _positive_real_vector(
        r_star_values[selected_indices], f"{summary_path}: selected l=2 R_star"
    )
    _validate_consistent_stellar_property(m_star, "M_star")
    _validate_consistent_stellar_property(r_star, "R_star")

    candidates = _detail_candidates(directory)
    modes: List[ModeRecord] = []
    structures: Dict[int, Tuple[np.ndarray, np.ndarray, np.ndarray, int, Path]] = {}
    for index in selected_indices:
        index_int = int(index)
        n_value = int(n_pg_int[index_int])
        mode_id = _summary_mode_id(id_values, index_int)
        detail_path, source_style = _locate_detail(candidates, n_value, mode_id, directory)
        grid, density, enclosed_mass, size = _validate_detail_file(
            detail_path,
            n_value,
            float(omega_real[index_int]),
            float(omega_imag[index_int]),
            mode_id,
            source_style,
        )
        structures[index_int] = (grid, density, enclosed_mass, size, detail_path)
        modes.append(
            ModeRecord(
                summary_index=index_int,
                n_pg=n_value,
                omega_real=float(omega_real[index_int]),
                omega_imag=float(omega_imag[index_int]),
                source_detail=detail_path,
                source_style=source_style,
                mode_id=mode_id,
            )
        )

    f_summary_index = int(selected_indices[int(f_matches[0])])
    reference_grid, reference_rho, reference_m_r, grid_size, reference_path = structures[
        f_summary_index
    ]
    if not np.isclose(reference_m_r[-1], m_star[0], rtol=1.0e-6, atol=0.0):
        raise ValueError(
            f"{reference_path}: surface M_r={reference_m_r[-1]:.8g} g is inconsistent "
            f"with summary M_star={m_star[0]:.8g} g. The bridge requires agreement "
            "within 1e-6 relative tolerance because TDECalculator combines them as one star."
        )
    for index_int, (grid, density, enclosed_mass, size, detail_path) in structures.items():
        if size != grid_size or not np.allclose(grid, reference_grid, rtol=1.0e-12, atol=1.0e-14):
            raise ValueError(
                f"{detail_path}: x grid differs from f-mode detail {reference_path}. The "
                "external TDECalculator combines every eigenfunction on that grid."
            )
        if not np.allclose(density, reference_rho, rtol=1.0e-12, atol=0.0):
            raise ValueError(
                f"{detail_path}: rho differs from f-mode detail {reference_path}; detail "
                "files may not describe one stellar structure."
            )
        if not np.allclose(enclosed_mass, reference_m_r, rtol=1.0e-12, atol=0.0):
            raise ValueError(
                f"{detail_path}: M_r differs from f-mode detail {reference_path}; detail "
                "files may not describe one stellar structure."
            )
    return GyreRunInspection(
        source_dir=directory,
        summary_path=summary_path,
        modes=tuple(modes),
        summary_mode_count=int(mode_count),
        l2_mode_count=int(l2_indices.size),
        n_pg_min=n_pg_min,
        n_pg_max=n_pg_max,
        f_mode_index=int(f_matches[0]),
        detail_grid_size=grid_size,
        detail_grid_reference=reference_path,
    )


def _copy_attributes(source: h5py.AttributeManager, target: h5py.AttributeManager) -> None:
    """Copy HDF5 attributes without converting their dtypes or values."""

    for key, value in source.items():
        target[key] = value


def _write_filtered_summary(
    source_path: Path, destination_path: Path, selected_indices: Sequence[int]
) -> None:
    """Create a root-dataset-preserving HDF5 summary filtered by row index.

    Dataset values and attributes retain their native HDF5 dtype.  Root
    datasets with leading dimension equal to the source mode count are sliced;
    scalar or unrelated datasets are copied unchanged.  This supports GYRE
    summaries with additional columns while retaining the calculator's exact
    ``pygyre.read_output`` layout.
    """

    indices = np.asarray(selected_indices, dtype=np.intp)
    with h5py.File(source_path, "r") as source, h5py.File(destination_path, "w") as target:
        _copy_attributes(source.attrs, target.attrs)
        mode_count = source["l"].shape[0]
        for name, item in source.items():
            if isinstance(item, h5py.Dataset):
                data = item[...]
                if item.ndim >= 1 and item.shape[0] == mode_count:
                    data = data[indices, ...]
                created = target.create_dataset(name, data=data, dtype=item.dtype)
                _copy_attributes(item.attrs, created.attrs)
            else:
                # Groups are not interpreted by pygyre's summary reader, but retaining
                # them avoids silently discarding unrelated provenance.
                source.copy(name, target)


def _transfer_file(source: Path, destination: Path, transfer_mode: str) -> None:
    """Transfer a validated detail file using the requested explicit mode."""

    if transfer_mode == "copy":
        shutil.copy2(source, destination)
    elif transfer_mode == "hardlink":
        os.link(source, destination)
    elif transfer_mode == "symlink":
        os.symlink(source, destination)
    else:  # Defensive check for callers outside the CLI.
        raise ValueError(f"Unknown transfer mode {transfer_mode!r}; choose from {TRANSFER_MODES}.")


def _mode_manifest_record(mode: ModeRecord, staged_path: Path) -> Dict[str, Any]:
    """Create one JSON-safe source-to-staged mapping record."""

    return {
        "summary_index": mode.summary_index,
        "n_pg": mode.n_pg,
        "mode_id": mode.mode_id,
        "omega": {"real": mode.omega_real, "imag": mode.omega_imag},
        "source_detail": str(mode.source_detail),
        "source_style": mode.source_style,
        "source_sha256": _sha256(mode.source_detail),
        "staged_detail": staged_path.name,
        "staged_sha256": _sha256(staged_path),
    }


def _source_run_provenance(source_dir: Path) -> Dict[str, Any]:
    """Capture optional clean-pipeline metadata without making it mandatory."""

    provenance: Dict[str, Any] = {}
    metadata_path = source_dir / "run_metadata.json"
    gyre_input_path = source_dir / "gyre.in"
    if metadata_path.is_file():
        provenance["run_metadata"] = str(metadata_path)
        provenance["run_metadata_sha256"] = _sha256(metadata_path)
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            provenance["run_metadata_parseable"] = False
        else:
            provenance["run_metadata_parseable"] = True
            if metadata.get("gyre_version") is not None:
                provenance["gyre_version"] = metadata["gyre_version"]
    if gyre_input_path.is_file():
        provenance["gyre_input"] = str(gyre_input_path)
        provenance["gyre_input_sha256"] = _sha256(gyre_input_path)
    return provenance


def _build_manifest(
    inspection: GyreRunInspection,
    temporary_directory: Path,
    transfer_mode: str,
) -> Dict[str, Any]:
    """Assemble reproducibility metadata after every staged file is present."""

    summary_path = temporary_directory / SUMMARY_FILENAME
    mappings = []
    for mode in inspection.modes:
        staged = temporary_directory / legacy_detail_filename(mode.n_pg)
        mappings.append(_mode_manifest_record(mode, staged))
    return {
        "schema_version": 1,
        "bridge": {
            "name": "tde_bridge",
            "version": BRIDGE_VERSION,
            "created_utc": _utc_now(),
        },
        "source": {
            "directory": str(inspection.source_dir),
            "summary": str(inspection.summary_path),
            "summary_sha256": _sha256(inspection.summary_path),
            "run_provenance": _source_run_provenance(inspection.source_dir),
        },
        "selection": {
            "harmonic_degree": 2,
            "n_pg_min": inspection.n_pg_min,
            "n_pg_max": inspection.n_pg_max,
            "requires_f_mode_n_pg": 0,
            "summary_mode_count": inspection.summary_mode_count,
            "available_l2_mode_count": inspection.l2_mode_count,
            "selected_mode_count": len(inspection.modes),
            "selected_n_pg": list(inspection.selected_n_pg),
            "detail_grid_size": inspection.detail_grid_size,
        },
        "transfer": {"mode": transfer_mode, "portable": transfer_mode == "copy"},
        "staged": {
            "summary": SUMMARY_FILENAME,
            "summary_sha256": _sha256(summary_path),
            "detail_mappings": mappings,
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "h5py": h5py.version.version,
        },
        "scientific_scope": {
            "uses_summary_field": "omega",
            "does_not_use_display_frequency": True,
            "warning": (
                "This file-contract validation does not validate applying the "
                "collaborator's mass-only M-gamma calibration to this stellar structure."
            ),
        },
    }


def legacy_detail_filename(n_pg: int) -> str:
    """Return the exact filename built by the collaborator's ``f'n{n:+d}'`` loop."""

    if not isinstance(n_pg, (int, np.integer)):
        raise ValueError("n_pg must be an integer.")
    return f"detail.l2.n{int(n_pg):+d}.h5"


def prepare_tde_directory(
    source_dir: Path | str,
    output_dir: Path | str,
    *,
    n_pg_min: Optional[int] = None,
    n_pg_max: Optional[int] = None,
    transfer_mode: str = "copy",
) -> PreparedTDERun:
    """Stage a fresh, validated, legacy-compatible TDE input directory.

    The destination must not exist.  Source HDF5 files are never modified.
    ``copy`` is the default and only portable transfer choice; ``hardlink``
    and ``symlink`` make source immutability/deletion a continuing concern.
    The result is assembled in a bridge-owned temporary sibling and atomically
    published only after the summary, every detail file, and manifest succeed.
    """

    if transfer_mode not in TRANSFER_MODES:
        raise ValueError(f"transfer_mode must be one of {TRANSFER_MODES}, got {transfer_mode!r}.")
    inspection = inspect_gyre_run(
        source_dir, n_pg_min=n_pg_min, n_pg_max=n_pg_max
    )
    destination = Path(output_dir).expanduser().resolve()
    if destination.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing TDE staging directory: {destination}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.parent / f".{destination.name}.tde-bridge-{uuid.uuid4().hex}"
    try:
        temporary.mkdir()
        _write_filtered_summary(
            inspection.summary_path,
            temporary / SUMMARY_FILENAME,
            [mode.summary_index for mode in inspection.modes],
        )
        for mode in inspection.modes:
            _transfer_file(
                mode.source_detail,
                temporary / legacy_detail_filename(mode.n_pg),
                transfer_mode,
            )
        manifest = _build_manifest(inspection, temporary, transfer_mode)
        manifest_path = temporary / MANIFEST_FILENAME
        manifest_path.write_text(
            json.dumps(_json_value(manifest), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(destination)
    except Exception:
        if temporary.exists():
            shutil.rmtree(temporary)
        raise
    return PreparedTDERun(
        directory=destination,
        summary_path=destination / SUMMARY_FILENAME,
        manifest_path=destination / MANIFEST_FILENAME,
        inspection=inspection,
        transfer_mode=transfer_mode,
    )


def load_manifest(path: Path | str) -> Dict[str, Any]:
    """Load and minimally validate a manifest written by this bridge."""

    manifest_path = Path(path).expanduser().resolve()
    if manifest_path.is_dir():
        manifest_path = manifest_path / MANIFEST_FILENAME
    with manifest_path.open("r", encoding="utf-8") as stream:
        manifest = json.load(stream)
    if manifest.get("schema_version") != 1 or manifest.get("bridge", {}).get("name") != "tde_bridge":
        raise ValueError(f"{manifest_path} is not a recognized tde_bridge manifest.")
    return manifest


def preflight_tde_project(tde_root: Path | str) -> TDEProject:
    """Check an official collaborator clone before any calculator import/run.

    A clone containing only a notebook cannot be run through this bridge: the
    reusable ``TDECalculator.py`` module and its cwd-relative ``M-gamma.txt``
    table are both mandatory for the unmodified external interface.
    """

    root = Path(tde_root).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"TDE project directory does not exist: {root}")
    calculator = root / "TDECalculator.py"
    if not calculator.is_file():
        raise FileNotFoundError(
            f"{root} lacks TDECalculator.py. A notebook-only TDE clone is not "
            "runnable through this bridge; obtain the collaborator's calculator module."
        )
    gamma_table = root / "M-gamma.txt"
    if not gamma_table.is_file():
        raise FileNotFoundError(
            f"{root} lacks M-gamma.txt, required because TDECalculator loads it relative to cwd."
        )
    try:
        table = np.loadtxt(gamma_table, ndmin=2)
    except (OSError, ValueError) as exc:
        raise ValueError(
            f"Could not parse {gamma_table} as the required two-column M-gamma table: {exc}"
        ) from exc
    if table.ndim != 2 or table.shape[0] < 2 or table.shape[1] != 2:
        raise ValueError(
            f"{gamma_table} must contain at least two rows and exactly two columns: "
            "stellar mass [M_sun] and gamma."
        )
    masses, gammas = table[:, 0], table[:, 1]
    if not np.all(np.isfinite(masses)) or not np.all(np.isfinite(gammas)):
        raise ValueError(f"{gamma_table} contains non-finite mass or gamma values.")
    if np.any(masses <= 0) or not np.all(np.diff(masses) > 0):
        raise ValueError(
            f"{gamma_table} masses must be positive and strictly increasing for interpolation."
        )
    if np.any(gammas <= 0):
        raise ValueError(
            f"{gamma_table} gamma values must be finite and strictly positive for the "
            "external disruption criterion."
        )
    return TDEProject(
        root=root,
        calculator_path=calculator,
        gamma_table_path=gamma_table,
        gamma_mass_domain_solar=(float(masses[0]), float(masses[-1])),
        calculator_sha256=_sha256(calculator),
        gamma_table_sha256=_sha256(gamma_table),
    )


@contextlib.contextmanager
def temporary_working_directory(path: Path | str) -> Iterator[Path]:
    """Temporarily change cwd and restore it even if the calculator raises.

    The process working directory is global state, so callers must not run
    this helper concurrently from multiple threads.
    """

    target = Path(path).expanduser().resolve()
    original = Path.cwd()
    os.chdir(target)
    try:
        yield target
    finally:
        os.chdir(original)


def load_tde_calculator(tde_root: Path | str) -> Type[Any]:
    """Import and return the collaborator's unmodified ``TDECalculator`` class.

    Loading occurs with cwd set to the TDE root as a defensive counterpart to
    :func:`run_tde_calculator`; it is restored before this function returns.
    """

    project = preflight_tde_project(tde_root)
    module_name = f"_tde_calculator_bridge_{uuid.uuid4().hex}"
    with temporary_working_directory(project.root):
        specification = importlib.util.spec_from_file_location(module_name, project.calculator_path)
        if specification is None or specification.loader is None:
            raise ImportError(f"Could not create an import specification for {project.calculator_path}")
        module = importlib.util.module_from_spec(specification)
        sys.modules[module_name] = module
        try:
            specification.loader.exec_module(module)
        except Exception:
            sys.modules.pop(module_name, None)
            raise
    calculator = getattr(module, "TDECalculator", None)
    if calculator is None or not callable(calculator):
        raise AttributeError(f"{project.calculator_path} does not define a callable TDECalculator.")
    return calculator


def run_tde_calculator(
    staged_star_dir: Path | str,
    tde_root: Path | str,
    *,
    MBH: float = 1.0e6,
    Rp: float = 17.0,
    N: int = 1000,
    calculator_class: Optional[Type[Any]] = None,
) -> Any:
    """Instantiate the external calculator with cwd-safe legacy table lookup.

    ``staged_star_dir`` must be an existing staged directory. Before import,
    the staged mass is checked against the finite M-gamma interpolation domain
    to prevent the external ``numpy.interp`` call from silently endpoint-
    clamping. The external class is supplied only as an injectable seam for
    tests; normal callers use the collaborator module imported by
    :func:`load_tde_calculator`. No TDE physics is duplicated here.
    """

    staged = Path(staged_star_dir).expanduser().resolve()
    if not staged.is_dir():
        raise FileNotFoundError(f"Staged TDE star directory does not exist: {staged}")
    project = preflight_tde_project(tde_root)
    if not np.isfinite(MBH) or MBH <= 0:
        raise ValueError("MBH must be finite and positive in solar masses.")
    if not np.isfinite(Rp) or Rp <= 0:
        raise ValueError("Rp must be finite and positive in GM_BH units.")
    if not isinstance(N, (int, np.integer)) or N < 2:
        raise ValueError("N must be an integer of at least 2 time samples.")
    staged_summary = staged / SUMMARY_FILENAME
    if not staged_summary.is_file():
        raise FileNotFoundError(f"Staged TDE summary not found: {staged_summary}")
    try:
        with h5py.File(staged_summary, "r") as summary:
            staged_m_star = _positive_real_vector(
                _require_vector(summary, "M_star"), f"{staged_summary}:M_star"
            )
    except OSError as exc:
        raise ValueError(f"Could not read staged summary {staged_summary}: {exc}") from exc
    _validate_consistent_stellar_property(staged_m_star, "M_star")
    staged_mass_solar = float(staged_m_star[0] / SOLAR_MASS_GRAMS)
    
    lower_mass, upper_mass = project.gamma_mass_domain_solar

    # Allow a very small numerical tolerance at the calibration boundaries.
    # The tolerance is in solar masses, not grams.
    mass_tol = 1e-4

    if not (lower_mass - mass_tol <= staged_mass_solar <= upper_mass + mass_tol):
        raise ValueError(
            f"Staged M_star={staged_mass_solar:.8g} M_sun is outside the M-gamma "
            f"table domain [{lower_mass:.8g}, {upper_mass:.8g}] M_sun, "
            f"even allowing a tolerance of {mass_tol:g} M_sun."
        )

    # Explicitly accept endpoint clamping only for values within the tolerance.
    if staged_mass_solar < lower_mass or staged_mass_solar > upper_mass:
        warnings.warn(
            f"Staged M_star={staged_mass_solar:.8g} M_sun is slightly outside "
            f"the M-gamma table domain [{lower_mass:.8g}, {upper_mass:.8g}] M_sun. "
            "The collaborator's interpolation will use the nearest endpoint.",
            RuntimeWarning,
            stacklevel=2,
        )
    
        calculator = calculator_class or load_tde_calculator(project.root)
        with temporary_working_directory(project.root):
            return calculator(str(staged), MBH=float(MBH), Rp=float(Rp), N=int(N))


def plot_quadrupole_spectrum(
    inspection: GyreRunInspection,
    ax: Any = None,
) -> Tuple[Any, Any, SpectrumPlotReport]:
    """Plot selected l=2 g/f/p modes against radial order and ``Re(omega)``.

    The vertical coordinate is GYRE's dimensionless angular frequency, the
    quantity used by the external calculator.  Any non-finite values supplied
    through a manually constructed inspection are explicitly reported in the
    returned :class:`SpectrumPlotReport` and warning, never silently repaired.
    The scoped Matplotlib style does not mutate a caller's global rcParams.
    """

    import matplotlib.pyplot as plt

    style = {
        "font.size": 16,
        "axes.labelsize": 20,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 15,
        "lines.linewidth": 2.2,
        "lines.markersize": 6,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "stix",
        "font.family": "STIXGeneral",
        "text.usetex": False,
    }
    orders = np.asarray([mode.n_pg for mode in inspection.modes], dtype=float)
    omegas = np.asarray([mode.omega_real for mode in inspection.modes], dtype=float)
    finite = np.isfinite(orders) & np.isfinite(omegas)
    report = SpectrumPlotReport(
        plotted_modes=int(np.count_nonzero(finite)),
        masked_nonfinite_modes=int(np.count_nonzero(~finite)),
    )
    if report.masked_nonfinite_modes:
        warnings.warn(
            f"Quadrupole spectrum masked {report.masked_nonfinite_modes} non-finite mode(s).",
            RuntimeWarning,
            stacklevel=2,
        )
    created = ax is None
    with plt.rc_context(style):
        if created:
            figure, axis = plt.subplots(figsize=(7.5, 4.8), constrained_layout=True)
        else:
            axis = ax
            figure = axis.figure
        classes = (
            (orders < 0, "g modes", "o", "#0072B2"),
            (orders == 0, "f mode", "D", "#D55E00"),
            (orders > 0, "p modes", "s", "#009E73"),
        )
        for mask, label, marker, color in classes:
            selected = finite & mask
            if np.any(selected):
                axis.scatter(
                    orders[selected],
                    omegas[selected],
                    s=42,
                    marker=marker,
                    color=color,
                    edgecolor="black",
                    linewidth=0.6,
                    label=label,
                )
        axis.set_xlabel(r"Radial order $n_{\rm pg}$")
        axis.set_ylabel(r"$\mathrm{Re}(\omega)$ [dimensionless]")
        axis.grid(True, alpha=0.22)
        handles, labels = axis.get_legend_handles_labels()
        if handles:
            axis.legend(handles, labels, frameon=False)
    return figure, axis, report
