# Molecular Orbital Analysis

This directory contains an Agent Skill and a PySCF-based command-line helper
for generating explicitly indexed molecular-orbital cube files.

- Read [SKILL.md](SKILL.md) for the input contract, closed-shell and open-shell
  workflow, validation criteria, and reporting requirements.
- Run `scripts/generate_orbital_cubes.py --help` from a Python environment with
  PySCF installed. The script validates one XYZ frame, requires explicit units,
  charge, spin, method, and basis, rejects unconverged or non-finite electronic
  results, validates every cube atom record and scalar payload value, and
  records accepted or rejected artifact provenance in `run.json`. It reads the
  source XYZ once, atomically publishes a run-local `input.xyz`, and accepts a
  run only when fresh non-empty `scf.chk` and `pyscf.log` hashes are present.

No global executable or scientific dependency is bundled. Install PySCF and a
cube-capable viewer separately and record their versions for reproducibility.
