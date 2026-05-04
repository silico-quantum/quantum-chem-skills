#!/usr/bin/env python3
"""Validate generated candidates with xTB GFN2 single-point."""
# Placeholder — runs xTB on generated SMILES, filters by HL gap

import argparse, subprocess, sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smi-file', required=True)
    parser.add_argument('--xtb-path', default='~/software/xtb-6.7.1/bin/xtb')
    parser.add_argument('--hl-gap-range', default='1.5-3.5')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    print(f'Validating SMILES from {args.smi_file}')
    print(f'xTB: {args.xtb_path}, HL gap window: {args.hl_gap_range}')
    print('Validation pipeline — requires RDKit 3D gen + xTB SP')

if __name__ == '__main__':
    main()
