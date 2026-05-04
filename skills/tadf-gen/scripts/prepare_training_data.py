#!/usr/bin/env python3
"""Prepare TADF training data from Stage 4 xTB SP CSVs."""
# Placeholder — reads Stage 4 status CSVs, extracts SMILES + properties

import argparse, csv, json, sys
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv-dir', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--property', default='hl_gap_ev')
    parser.add_argument('--min-gap', type=float, default=1.5)
    parser.add_argument('--max-gap', type=float, default=3.5)
    args = parser.parse_args()

    csv_dir = Path(args.csv_dir)
    data = []
    for csv_file in sorted(csv_dir.glob('status_*.csv')):
        with open(csv_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row.get('status') != 'ok':
                    continue
                hl = float(row.get('hl_gap_ev', 0))
                if hl < args.min_gap or hl > args.max_gap:
                    continue
                data.append({
                    'mol': row['mol'],
                    'hl_gap_ev': hl,
                    'energy_eh': row.get('energy_eh', ''),
                })

    with open(args.output, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'Wrote {len(data)} entries to {args.output}')

if __name__ == '__main__':
    main()
