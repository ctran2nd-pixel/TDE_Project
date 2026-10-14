#!/usr/bin/env python3
"""Run a reproducible adiabatic GYRE calculation from a MESA pulse file.

Example
-------
python scripts/run_gyre/pipeline/run_gyre.py models/my_star.mesa \\
  --output-dir results/run_gyre/my_star
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

from gyre_pipeline import FREQUENCY_UNITS, RunConfig, run_pipeline


def build_parser() -> argparse.ArgumentParser:
    """Build the documented command-line interface."""

    parser = argparse.ArgumentParser(
        description=(
            "Generate a GYRE input file, run adiabatic stellar oscillations, "
            "and export HDF5 and CSV mode products."
        )
    )
    parser.add_argument(
        "model",
        type=Path,
        help="MESA pulse-data model written with pulse_data_format='GYRE'.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Run directory (default: results/run_gyre/<model stem>).",
    )
    parser.add_argument(
        "--gyre-executable",
        type=Path,
        help="Explicit path to the compiled gyre frontend.",
    )
    parser.add_argument(
        "--l",
        dest="harmonic_degrees",
        type=int,
        nargs="+",
        default=[2],
        help="One or more harmonic degrees (default: 2).",
    )
    parser.add_argument("--n-pg-min", type=int, default=-20)
    parser.add_argument("--n-pg-max", type=int, default=20)
    parser.add_argument(
        "--freq-min",
        type=float,
        default=0.1,
        help="Minimum scan frequency in --freq-units (default: 0.1).",
    )
    parser.add_argument(
        "--freq-max",
        type=float,
        default=20.0,
        help="Maximum scan frequency in --freq-units (default: 20).",
    )
    parser.add_argument("--n-freq", type=int, default=2000)
    parser.add_argument(
        "--freq-units",
        type=str.upper,
        choices=FREQUENCY_UNITS,
        default="NONE",
        help=(
            "Units of the scan bounds; NONE means dimensionless angular "
            "frequency (default: NONE)."
        ),
    )
    parser.add_argument(
        "--output-freq-units",
        type=str.upper,
        choices=FREQUENCY_UNITS,
        default="CYC_PER_DAY",
        help="Units used for the output freq column (default: CYC_PER_DAY).",
    )
    parser.add_argument(
        "--grid-type",
        type=str.upper,
        choices=("LINEAR", "INVERSE"),
        default="INVERSE",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="OpenMP threads supplied to GYRE (default: 1).",
    )
    parser.add_argument(
        "--w-osc",
        type=float,
        default=10.0,
        help="Spatial-grid oscillatory-region weight (default: 10).",
    )
    parser.add_argument(
        "--w-exp",
        type=float,
        default=2.0,
        help="Spatial-grid evanescent-region weight (default: 2).",
    )
    parser.add_argument(
        "--w-ctr",
        type=float,
        default=10.0,
        help="Spatial-grid central-region weight (default: 10).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace an existing run directory.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and generate inputs without requiring or launching GYRE.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Parse arguments and run the pipeline with concise error reporting."""

    parser = build_parser()
    args = parser.parse_args(argv)
    output_dir = args.output_dir or Path("results/run_gyre") / args.model.stem
    config = RunConfig(
        harmonic_degrees=tuple(args.harmonic_degrees),
        n_pg_min=args.n_pg_min,
        n_pg_max=args.n_pg_max,
        freq_min=args.freq_min,
        freq_max=args.freq_max,
        n_freq=args.n_freq,
        freq_units=args.freq_units,
        output_freq_units=args.output_freq_units,
        grid_type=args.grid_type,
        w_osc=args.w_osc,
        w_exp=args.w_exp,
        w_ctr=args.w_ctr,
        threads=args.threads,
        n_iter_max=N_ITER_MAX,
)

    try:
        result = run_pipeline(
            args.model,
            output_dir,
            config,
            gyre_executable=args.gyre_executable,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
        )
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"Run directory: {result.output_dir}")
    print(f"GYRE input:    {result.input_file}")
    print(f"Metadata:      {result.metadata_file}")
    if result.dry_run:
        print("Status:        dry run; GYRE was not launched")
    else:
        print(f"GYRE log:      {result.log_file}")
        print(f"Summary:       {result.summary_file}")
        print(f"Mode table:    {result.modes_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
