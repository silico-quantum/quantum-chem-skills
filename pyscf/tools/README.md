# Bundled PySCF prototype status

The Python files in this directory are historical prototypes, not supported
entry points. Direct execution is disabled by `_legacy_guard.py` before NumPy,
PySCF, or another optional integration is imported. Their definitions remain
available for static review in an appropriate environment, but no function in
this directory is part of the supported runtime contract.

These prototypes have not been run against a pinned PySCF environment by the
current repository validation. Syntax compilation alone does not establish API
compatibility or scientific validity. Use
`../scripts/run_safe_dft_tda.py` for the one supported bundled calculation: a
closed-shell RKS/TDA run with strict result and provenance checks.

Do not import or reuse their functions in production until each script has:

1. a pinned PySCF/Python environment;
2. a strict input schema with explicit geometry units, charge, and PySCF spin;
3. convergence and finite-value checks for every stage;
4. a no-overwrite output contract and provenance manifest;
5. a small reference calculation compared with official PySCF examples; and
6. an automated regression test for parsing, units, and failure behavior.

Known reasons for quarantine include:

| File | Blocking issue examples |
|---|---|
| `analysis.py` | checkpoint reconstruction and population/dipole APIs are unvalidated |
| `cascf.py` | advertises perturbation interfaces that are not established by official core examples |
| `ccsd.py` | incomplete open-shell/reporting and downstream convergence handling |
| `dft.py` | unvalidated functional recommendations and grid API; `compare_functionals` still contains the quarantined unsupported `output_file` call |
| `geometry.py` | optimizer, transition-state, and frequency interfaces require optional version-specific integrations |
| `mp2.py` | input and convergence/reporting contract is incomplete |
| `pes.py` | geometry mutation, constraints, state tracking, and failed-point policy are unvalidated |
| `scf.py` | legacy CLI lacks strict XYZ and accepted-result output contracts |
| `spectrum.py` | contains unestablished imports and quarantined `td.f` access plus unsafe broadening paths |
| `tddft.py` | contains quarantined `td.f` access and incomplete spin/root handling; use `td.oscillator_strength()` only in a validated runner |

The same direct-execution quarantine applies to
`../scripts/dft_calculation.py` and `../references/benzene-dft-tddft.py`; both
stop before optional scientific packages are imported. Other reference
documents remain review-only. Use the compact workflow in `../SKILL.md`, the
supported runner, and official documentation as the starting point.

## Historical reference documents

The long documents below are retained for review but are not authoritative API
documentation and must not be executed as recipes without line-by-line
verification against the installed version:

- `../references/theory/pyscf-api-reference.md`
- `../references/theory/pyscf-advanced.md`
- `../references/theory/pyscf-jax-integration.md`
- `../references/practice/2d-potential-energy-surface.md`
- `../references/practice/emission-spectrum-guide.md`
- `../references/practice/emission-spectrum-workflow.md`
- `../references/benzene-dft-tddft.py`

In particular, converting PySCF NumPy arrays to JAX arrays does not make an SCF
calculation end-to-end differentiable. Any differentiable or accelerator
integration requires a dedicated, documented backend and its own validation.
