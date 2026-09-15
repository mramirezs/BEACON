#!/usr/bin/env python3
"""
BEACON - Parsear salida de DiscoTope-3.0 al formato unificado por residuo.

DiscoTope-3.0 genera CSVs en: <out_dir>/<pdb_name>/<pdb_name>_<chain>_discotope3.csv
Columnas: pdb, chain, res_id, residue, DiscoTope-3.0_score, calibrated_score,
          epitope, rsa, pLDDTs, length, alphafold_struc_flag

Este script busca el CSV recursivamente y genera:
  res_id, residue, discotope3_score, rsa, prediction

Uso:
    python3 parse_discotope3.py <dt3_output_dir> <output_csv>
    python3 parse_discotope3.py results/ag85b/discotope3/ results/ag85b/discotope3/discotope3_parsed.csv
"""

import sys
import os
import csv
import glob


THRESHOLD = -3.7  # DiscoTope-3.0 default threshold


def find_dt3_csv(output_dir):
    """Find DiscoTope-3.0 output CSV, searching subdirectories."""
    # DiscoTope creates: <out_dir>/<pdb_name>/<pdb_name>_<chain>_discotope3.csv
    patterns = [
        os.path.join(output_dir, "**", "*_discotope3.csv"),
        os.path.join(output_dir, "*_discotope3.csv"),
        os.path.join(output_dir, "**", "*.csv"),
    ]
    for pattern in patterns:
        files = glob.glob(pattern, recursive=True)
        # Exclude our own parsed output
        files = [f for f in files if 'parsed' not in os.path.basename(f)]
        if files:
            return files[0]
    return None


def parse_dt3_csv(csv_file):
    """Parse DiscoTope-3.0 output CSV."""
    residues = []
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames

        # Find score column
        score_col = None
        for c in cols:
            if 'discotope' in c.lower() and 'score' in c.lower():
                score_col = c
                break
        if score_col is None:
            # Try calibrated_score as fallback
            for c in cols:
                if 'calibrated' in c.lower():
                    score_col = c
                    break
        if score_col is None:
            print(f"ERROR: No se encontro columna de score")
            print(f"  Columnas: {cols}")
            sys.exit(1)

        print(f"  Columna de score: '{score_col}'")

        for row in reader:
            score = float(row[score_col])
            res_id = int(row['res_id'])
            residue = row['residue']
            rsa = float(row.get('rsa', 0))

            # Prediction: use 'epitope' column if available, otherwise threshold
            if 'epitope' in row:
                pred = 1 if row['epitope'].strip().lower() == 'true' else 0
            else:
                pred = 1 if score >= THRESHOLD else 0

            residues.append({
                'res_id': res_id,
                'residue': residue,
                'discotope3_score': round(score, 6),
                'rsa': round(rsa, 5),
                'prediction': pred
            })

    return residues


def main():
    if len(sys.argv) < 3:
        print("Uso: python3 parse_discotope3.py <dt3_output_dir> <output_csv>")
        sys.exit(1)

    dt3_dir = sys.argv[1]
    output_csv = sys.argv[2]

    csv_file = find_dt3_csv(dt3_dir)
    if csv_file is None:
        print(f"ERROR: No se encontro CSV de DiscoTope-3.0 en {dt3_dir}")
        print("  Esperado: <dir>/<pdb_name>/<pdb_name>_<chain>_discotope3.csv")
        sys.exit(1)

    print(f"Parseando: {csv_file}")
    residues = parse_dt3_csv(csv_file)

    # Write output
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['res_id', 'residue', 'discotope3_score', 'rsa', 'prediction'])
        writer.writeheader()
        writer.writerows(residues)

    n_epitope = sum(1 for r in residues if r['prediction'] == 1)
    print(f"\nEscrito: {output_csv}")
    print(f"  Total residuos: {len(residues)}")
    print(f"  Rango PDB: {residues[0]['res_id']}–{residues[-1]['res_id']}")
    print(f"  Predichos epitope: {n_epitope} ({100*n_epitope/len(residues):.1f}%)")


if __name__ == "__main__":
    main()
