# TADF Molecular Generator (tadf-gen)

Generative deep learning skill for TADF (Thermally Activated Delayed Fluorescence) molecular design. Uses graph-based RNN generation models trained on high-throughput xTB/DFT screening data to produce novel blue narrowband TADF candidates.

## Overview

This skill integrates generative models with the blue_1m TADF screening pipeline (see `../tadf-screening/SKILL.md`). Instead of exhaustive combinatorial enumeration of donor-acceptor fragment pairs, the generative model learns the distribution of successful TADF molecules and samples novel candidates from the learned manifold.

### Pipeline Integration

```
tadf-screening (combinatorial)          tadf-gen (generative)
─────────────────────────────          ──────────────────────
Stage 1-2: Fragment enumeration
Stage 3: 3D gen (800k XYZ)            
Stage 4: xTB SP (1.26M)               ← Training data
Stage 5: xTB OPT (434k)               ← Property labels
                                       ↓
                              Train Gen-DL model
                              Generate novel candidates
                              Filter × property predictor
                                       ↓
Stage 6: sTDA (top-N)                 ← Validate generated
Stage 7-8: DFT + MOMAP
```

## Generative Model: DeepMoleculeGen

### Source

- **Repository**: [spark8ku/DeepMoleculeGen](https://github.com/spark8ku/DeepMoleculeGen)
- **Model**: Gen-DL — RNN-based molecular graph generator with graph convolutions
- **Framework**: Apache MXNet 1.8.0
- **Reference**: [Scientific Data 7, 289 (2020)](https://www.nature.com/articles/s41597-020-00634-8)
- **Web**: <http://deep4chem.korea.ac.kr/DeepMoleculeGen>

### Architecture

```
Input (atom init) → Atom Embedding → Graph Conv (F_skip) → Dense → Policy
                        ↓                                          ↓
                   Bond Embedding → RNN (N_rnn=3) → Atom + Bond prediction
```

Key components:
- **Embedding**: Atom types (N_A=7), bond types (N_B), mask tokens
- **Graph Convolution**: Message-passing layers with feature dimensions F_h=[32,64,128,128,256,256]
- **Policy Network**: Fh_policy=128, outputs next atom + bond predictions
- **Autoregressive generation**: Builds molecular graph atom-by-atom, bond-by-bond

### Requirements

```
mxnet==1.8.0
numpy==1.20.3
scipy==1.10.1
networkx==3.1
rdkit (for molecular graph operations)
```

### Pre-trained Model

- Location: `tools/DeepMoleculeGen/ckpt/DeepMoleculeGen/`
- `config.json`: Model hyperparameters (N_C=7, F_e=16, N_rnn=3, F_skip=256)
- `ckpt.7z`: Pre-trained weights (Git LFS, ~46MB)
- Note: Checkpoint trained on general organic molecules with optical properties — to be fine-tuned on TADF-specific data

## Installation

### On GPU Node (c20, 8×RTX 5090)

```bash
# Python 3.11 is available at /opt/Python-3.11.4
export PATH=/opt/Python-3.11.4/bin:$PATH

# Install MXNet with GPU support
pip install mxnet-cu112  # CUDA 12.8 needs compatible build

# Install dependencies
cd skills/tadf-gen/tools/DeepMoleculeGen
pip install -r requirement.txt
pip install rdkit-pypi
```

### Local (Mac)

```bash
pip install mxnet  # CPU only
pip install -r requirement.txt
pip install rdkit
```

## Usage

### 1. Prepare Training Data from Screening Pipeline

```bash
# Extract SMILES + properties from Stage 4 CSVs
python3 scripts/prepare_training_data.py \
  --csv-dir tadf-screening/data/blue_1m/logs/ \
  --output data/tadf_train.json \
  --property hl_gap_ev \
  --min-gap 1.5 --max-gap 3.5
```

### 2. Fine-tune Pre-trained Model

```bash
python3 scripts/train_tadf.py \
  --config tools/DeepMoleculeGen/ckpt/DeepMoleculeGen/config.json \
  --checkpoint tools/DeepMoleculeGen/ckpt/DeepMoleculeGen/ckpt.7z \
  --train-data data/tadf_train.json \
  --epochs 100 \
  --device gpu
```

### 3. Generate Novel Candidates

```bash
python3 scripts/generate.py \
  --model models/tadf-gen-params.npz \
  --num-molecules 10000 \
  --output data/generated_tadf.smi \
  --device gpu
```

### 4. Validate Generated Molecules

```bash
python3 scripts/validate.py \
  --smi-file data/generated_tadf.smi \
  --xtb-path ~/software/xtb-6.7.1/bin/xtb \
  --hl-gap-range 1.5-3.5 \
  --output data/validated_candidates.json
```

## Property Prediction

The model can be extended with a property prediction head to directly estimate:
- HOMO-LUMO gap (xTB level)
- S₁ excitation energy (sTDA level)
- Oscillator strength
- ΔE_ST (singlet-triplet gap)

Using the 1.26M xTB SP results from Stage 4 as training labels.

## Scripts

| Script | Purpose |
|--------|---------|
| `scripts/prepare_training_data.py` | Extract SMILES + properties from Stage 4 CSVs |
| `scripts/train_tadf.py` | Fine-tune Gen-DL on TADF data |
| `scripts/generate.py` | Sample novel molecular graphs |
| `scripts/validate.py` | Screen generated candidates with xTB |

## Data Flow

```
Stage 4 CSVs (1.26M rows)
    ↓ scripts/prepare_training_data.py
tadf_train.json (SMILES + HL gap + energy)
    ↓ scripts/train_tadf.py (GPU: c20)
tadf-gen-params.npz (fine-tuned weights)
    ↓ scripts/generate.py
generated_tadf.smi (novel candidates)
    ↓ scripts/validate.py (xTB SP)
validated_candidates.json (filtered hits)
    ↓ Stage 6 sTDA pipeline
```

## References

1. DeepMoleculeGen: Gen-DL model for molecular properties, *Sci Data* 7, 289 (2020).  
   <https://doi.org/10.1038/s41597-020-00634-8>
2. `spark8ku/DeepMoleculeGen` — <https://github.com/spark8ku/DeepMoleculeGen>
3. blue_1m TADF screening pipeline — see `../tadf-screening/SKILL.md`

## License

DeepMoleculeGen is distributed under its original license (see `tools/DeepMoleculeGen/LICENSE`).  
This skill wrapper and training scripts are part of the quantum-chem-skills monorepo.
