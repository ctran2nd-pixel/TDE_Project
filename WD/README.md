# Fast carbon-oxygen white-dwarf case

This directory contains one runnable MESA `r24.08.1` work directory:

- `make_co_wd_fast/`

The case starts from the bundled `co_core.mod` checkpoint and therefore does
not simulate the pre-main-sequence stage. It also skips the main-sequence,
core-helium-burning, and C/O-core construction stages. The two stages that are
run are:

1. `inlist_remove_env`: remove the AGB envelope and write `co_wd.mod`.
2. `inlist_settle`: enable element diffusion and write `final.mod`.

The checkpoint was produced with MESA `r24.08.1` from the standard
approximately `3.1 Msun`, `Z = 0.02` progenitor. Its stored mass is
approximately `3.0784 Msun` before envelope removal. Do not change
`initial_mass` or `initial_z` and expect this checkpoint to represent a
different progenitor. A different mass or metallicity requires a matching
checkpoint generated in a separate, longer run.

## Run

Use MESA `r24.08.1` and make sure `MESA_DIR` points to that installation:

```bash
export MESA_DIR=/path/to/mesa-r24.08.1
cd /path/to/mesa_project/WD/make_co_wd_fast
./mk
./rn
```

`./mk` is required. It builds the local `star` executable with this case's
`src/run_star_extras.f90`; copying a generic MESA `star` executable will omit
the case-specific envelope-removal logic.

The run writes live terminal output and also saves:

- `remove_env.log` and `LOGS_remove_env/`
- `settle.log` and `LOGS_settle/`
- `co_wd.mod`, the young white dwarf
- `final.mod`, the settled white dwarf

Monitor a stage from another terminal with:

```bash
tail -f remove_env.log
```

If output from an earlier run exists, `./rn` moves it into a timestamped
`archive/` directory before starting.
