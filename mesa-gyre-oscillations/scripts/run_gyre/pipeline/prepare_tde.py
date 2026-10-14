#!/usr/bin/env python3
"""Inspect or stage clean GYRE output for the collaborator's TDECalculator.

Examples
--------
Read-only contract check::

    python scripts/run_gyre/pipeline/prepare_tde.py inspect results/run_gyre/my_star

Create a portable compatibility directory::

    python scripts/run_gyre/pipeline/prepare_tde.py prepare results/run_gyre/my_star \
        --output-dir results/tde_stage/my_star

The command never changes the source run and refuses to overwrite a staging
directory.  ``--n-pg-min``/``--n-pg-max`` are explicit exploratory filters;
the selected set still has to contain the f mode (``n_pg=0``).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from tde_bridge import (
    TRANSFER_MODES,
    inspect_gyre_run,
    prepare_tde_directory,
    preflight_tde_project,
)


def _inspection_report(inspection: Any) -> Dict[str, Any]:
    """Make a compact, JSON-friendly read-only preflight report."""

    return {
        "source_dir": str(inspection.source_dir),
        "summary_path": str(inspection.summary_path),
        "summary_mode_count": inspection.summary_mode_count,
        "available_l2_mode_count": inspection.l2_mode_count,
        "selected_mode_count": len(inspection.modes),
        "selected_n_pg": list(inspection.selected_n_pg),
        "n_pg_min": inspection.n_pg_min,
        "n_pg_max": inspection.n_pg_max,
        "detail_grid_size": inspection.detail_grid_size,
        "detail_mappings": [
            {
                "n_pg": mode.n_pg,
                "mode_id": mode.mode_id,
                "source_detail": str(mode.source_detail),
                "source_style": mode.source_style,
            }
            for mode in inspection.modes
        ],
        "warning": (
            "File validation does not validate applying the collaborator's mass-only "
            "M-gamma calibration to this stellar structure."
        ),
    }


def _add_selection_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--n-pg-min",
        type=int,
        help="Inclusive lower radial-order filter; recorded in the manifest.",
    )
    parser.add_argument(
        "--n-pg-max",
        type=int,
        help="Inclusive upper radial-order filter; recorded in the manifest.",
    )
    parser.add_argument(
        "--tde-root",
        type=Path,
        help="Optionally preflight a collaborator TDE project (TDECalculator.py and M-gamma.txt).",
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the two-operation command-line interface."""

    parser = argparse.ArgumentParser(
        description="Validate and stage l=2 GYRE HDF5 files for TDECalculator."
    )
    subparsers = parser.add_subparsers(dest="operation", required=True)

    inspect_parser = subparsers.add_parser("inspect", help="Read-only GYRE/TDE contract preflight.")
    inspect_parser.add_argument("source_dir", type=Path, help="GYRE run directory containing summary.h5.")
    _add_selection_arguments(inspect_parser)

    prepare_parser = subparsers.add_parser(
        "prepare", help="Create a fresh compatibility directory after preflight."
    )
    prepare_parser.add_argument("source_dir", type=Path, help="GYRE run directory containing summary.h5.")
    prepare_parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="New TDE-compatible directory; existing paths are never overwritten.",
    )
    prepare_parser.add_argument(
        "--transfer-mode",
        choices=TRANSFER_MODES,
        default="copy",
        help="Detail-file transfer (default: copy, the portable choice).",
    )
    prepare_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and print the intended handoff without creating output-dir.",
    )
    _add_selection_arguments(prepare_parser)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run an operation with concise errors suitable for collaborator shells."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        inspection = inspect_gyre_run(
            args.source_dir, n_pg_min=args.n_pg_min, n_pg_max=args.n_pg_max
        )
        report = _inspection_report(inspection)
        if args.tde_root is not None:
            project = preflight_tde_project(args.tde_root)
            report["tde_project"] = {
                "root": str(project.root),
                "calculator": str(project.calculator_path),
                "calculator_sha256": project.calculator_sha256,
                "gamma_table": str(project.gamma_table_path),
                "gamma_table_sha256": project.gamma_table_sha256,
                "gamma_mass_domain_solar": list(project.gamma_mass_domain_solar),
            }
        if args.operation == "prepare" and not args.dry_run:
            prepared = prepare_tde_directory(
                args.source_dir,
                args.output_dir,
                n_pg_min=args.n_pg_min,
                n_pg_max=args.n_pg_max,
                transfer_mode=args.transfer_mode,
            )
            report["prepared"] = {
                "directory": str(prepared.directory),
                "summary": str(prepared.summary_path),
                "manifest": str(prepared.manifest_path),
                "transfer_mode": prepared.transfer_mode,
            }
        elif args.operation == "prepare":
            report["dry_run"] = True
            report["would_create"] = str(args.output_dir.expanduser().resolve())
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0
    except (OSError, ValueError, ImportError, AttributeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
