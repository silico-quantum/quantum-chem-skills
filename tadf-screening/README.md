# 🧪 TADF Emitter Screening Pipeline

> **High-throughput computational screening of Thermally Activated Delayed Fluorescence (TADF) emitters.**

[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/)
[![MOMAP](https://img.shields.io/badge/MOMAP-2024A-orange.svg)](http://www.momap.cn)

## 🚀 Four-Stage Screening Pipeline

```
 Stage 1           Stage 2           Stage 3              Stage 4 ✨
┌──────────┐     ┌──────────┐     ┌──────────────┐     ┌──────────────┐
│ Library  │ →   │  xTB     │ →   │  Gaussian    │ →   │   MOMAP      │
│ RDKit    │     │ GFN2-xTB │     │ B3LYP/6-31G* │     │  TVCF/ISC    │
│ D-A enum │     │ sTDA pre │     │ S0/T1/S1 opt │     │  k_RISC, Φ   │
└──────────┘     └──────────┘     └──────────────┘     └──────────────┘
  20k→5k           5k→74             74→TBD               TBD→TOP 🏆
```

### Stage 1 — Library Generation
Stochastic Donor-Acceptor fragment assembly via RDKit/DeepChem.
**Input:** D/A fragment pools  **→  Output:** 3D structures (XYZ/SDF)

### Stage 2 — xTB Pre-screening
GFN2-xTB geometry optimization + sTDA excitation.
**Filter:** S1 in 350–500 nm, f ≥ 0.05  **→  Output:** Shortlisted candidates

### Stage 3 — Gaussian TDDFT
B3LYP/6-31G(d) S0/T1 optimization + S1 TDDFT (marcus cluster).
**Filter:** Normal termination × 2 per molecule  **→  Output:** `.log` + `.fchk` files

### Stage 4 — MOMAP Photophysics ✨
Full photophysics via MOMAP 2024A TVCF method.
**Computes:** EVC → fluorescence spectrum → ISC rate → quantum yield
**Filter:** Peak emission in 450–490 nm + k_RISC > k_r

```bash
# Single molecule
python stage4_momap.py --mol-id mol_07566 \
    --s0 logs/s0.log --s1 logs/s1.log --t1 logs/t1.log

# Batch from CSV
python stage4_momap.py candidates.csv --output stage4_output/

# Columns: mol_id,s0_log,s1_log,t1_log
```

**Output:**
- `stage4_output/<mol_id>/` — MOMAP EVC, spectrum, ISC results
- `stage4_output/stage4_report.md` — Ranked report with 🔵 blue window flags
- `stage4_output/stage4_results.json` — Machine-readable results

---

## 📈 Example: Blue TADF Discovery

| Candidate ID | λ_emi (nm) | ΔE_ST (eV) | f | Blue Window | Stage 3 → 4 |
|:---:|:---:|:---:|:---:|:---:|:---:|
| mol_07566 | 450.4 | 0.XX | 0.0938 | 🔵 YES | Stage 3 |
| mol_06765 | 448.3 | 0.XX | 0.1104 | 🟡 near | Stage 3 |
| mol_07825 | 362.6 | 0.XX | 0.0295 | ❌ UV | → filtered out |

*Full Stage 4 report in `stage4_report.md`*

---

## 🛠️ Pipeline Dependencies

| Stage | Tools | Server |
|-------|-------|--------|
| 1 | RDKit, DeepChem | Local |
| 2 | xTB 6.7.1 | marcus2 |
| 3 | Gaussian 16, formchk | marcus (Slurm) |
| 4 | MOMAP 2024A | marcus2 / Slurm |

---
**Silico (硅灵)** 🔮 — AI Research Partner
