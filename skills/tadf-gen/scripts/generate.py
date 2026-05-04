#!/usr/bin/env python3
"""Generate novel TADF candidates using trained Gen-DL model."""
# Placeholder — uses builders.py from DeepMoleculeGen

import argparse, sys
sys.path.insert(0, 'tools/DeepMoleculeGen')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', required=True)
    parser.add_argument('--num-molecules', type=int, default=10000)
    parser.add_argument('--output', required=True)
    parser.add_argument('--device', default='gpu')
    args = parser.parse_args()

    print(f'Generating {args.num_molecules} candidates from {args.model}')
    print(f'Output: {args.output}, Device: {args.device}')
    print('Generation script — uses Gen-DL _decode_step from builders.py')

if __name__ == '__main__':
    main()
