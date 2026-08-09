# TADF Stage 4 MOMAP Adapter

This directory contains only the Stage 4 adapter that connects prepared
Gaussian outputs to the repository's MOMAP tooling. It does not include the
library-generation, xTB prescreening, or Gaussian execution stages of a full
TADF screening pipeline.

For the separately maintained end-to-end project, see
[silico-quantum/tadf-screening](https://github.com/silico-quantum/tadf-screening).

## Scope

`stage4_momap.py` accepts S0/S1 and optional T1 Gaussian log files, delegates
the per-molecule calculation to `../momap/tools/tadf.py`, ranks successful
results against an emission window, and writes a Markdown report plus JSON.
MOMAP and Gaussian are separately licensed external programs; this adapter
does not install or validate them.

Each Gaussian log must have a matching `.fchk`, or a matching `.chk` that can
be converted with `formchk`. The adapter stages the log and formatted
checkpoint together in the molecule work directory before invoking MOMAP.

The adjacent `../momap/tools` directory is used by default. If the adapter is
deployed separately, set `MOMAP_TOOLS_DIR` to the directory containing
`tadf.py`; the adapter does not assume a Codex, Claude Code, OpenClaw, or
Copilot installation path.

## Usage

```bash
# Single-molecule mode
python stage4_momap.py --mol-id mol_07566 \
    --s0 logs/s0.log --s1 logs/s1.log --t1 logs/t1.log

# Run ISC only with an explicit, provenance-backed S1-T1 SOC value
python stage4_momap.py --mol-id mol_07566 \
    --s0 logs/s0.log --s1 logs/s1.log --t1 logs/t1.log \
    --hso-cm1 12.5

# Batch mode
python stage4_momap.py candidates.csv --output stage4_output/ --target blue

# CSV columns: mol_id,s0_log,s1_log,t1_log,hso_cm1
```

Use `--window MIN MAX` for a custom range, or `--config PATH` to read an
emission target from an upstream workflow configuration. Run
`python stage4_momap.py --help` for the complete local interface.

ISC is skipped unless `hso_cm1` is supplied either on the command line or in
the CSV row. The JSON result then records `isc.status` as `not_computed`; the
workflow never substitutes a placeholder SOC value. Internally, Stage 4 asks
`tadf.py` to write a dedicated `--json-output` file, so progress messages on
stdout are not interpreted as machine-readable results.

For S1, the reported total energy is the last SCF reference energy plus the
last state-1 excitation energy. The absorption and emission transition
dipoles are state 1 from the first and last TD dipole tables, respectively;
their XYZ vector norm is converted from atomic units to Debye.

Expected outputs under the selected output directory are:

- per-molecule files produced by the MOMAP workflow;
- `stage4_report.md`, unless changed with `--report`;
- `stage4_results.json` with machine-readable results.

The bundled images are presentation assets only; they are not evidence that a
calculation can be reproduced in the current environment. In particular,
`examples/4czipn_clean.png` and `examples/4czipn_bonds.png` are byte-identical
in this snapshot, so they must not be treated as a validated style comparison.
