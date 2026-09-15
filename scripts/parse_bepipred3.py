#!/usr/bin/env python3
"""
BEACON - Parsear salida de BepiPred-3.0 al formato unificado por residuo.

BepiPred-3.0 genera:
  - raw_output.csv: scores por residuo (columnas: Accession, Residue,
    BepiPred-3.0 score, BepiPred-3.0 linear epitope score)
  - Bcell_epitope_preds.fasta: MAYUSCULA=epitope, minuscula=no

Este script lee raw_output.csv (preferido, tiene scores) o el FASTA,
y genera un CSV estandarizado con numeracion PDB.

Uso:
    python3 parse_bepipred3.py <bp3_output_dir> <output_csv> [--start_resid N]
    python3 parse_bepipred3.py results/ag85b/bepipred3/ results/ag85b/bepipred3/bepipred3_parsed.csv --start_resid 2

    --start_resid: posicion PDB del primer residuo (default: 1).
                   Ag85B empieza en 2, Spike RBD en 333, L1 HPV-16 en 20.
"""

import sys
import os
import csv
import glob


THRESHOLD = 0.1496  # BepiPred-3.0 default (percentil 80)


def parse_bp3_csv(csv_file, start_resid=1):
    """Parse raw_output.csv de BepiPred-3.0."""
    residues = []
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        cols = reader.fieldnames
        
        # Find score column — look for partial match
        score_col = None
        for c in cols:
            if 'bepipred' in c.lower() and 'score' in c.lower() and 'linear' not in c.lower():
                score_col = c
                break
        if score_col is None:
            # Fallback: first column with "score" in name
            for c in cols:
                if 'score' in c.lower():
                    score_col = c
                    break
        if score_col is None:
            print(f"ERROR: No se encontro columna de score en {csv_file}")
            print(f"  Columnas disponibles: {cols}")
            sys.exit(1)
        
        # Find residue column
        res_col = None
        for c in cols:
            if c.lower() in ('residue', 'aa', 'amino_acid'):
                res_col = c
                break
        
        print(f"  Columna de score: '{score_col}'")
        print(f"  Columna de residuo: '{res_col}'")
        
        for i, row in enumerate(reader):
            score = float(row[score_col].strip())
            aa = row[res_col].strip() if res_col else '?'
            residues.append({
                'res_id': start_resid + i,
                'residue': aa,
                'bepipred3_score': round(score, 6),
                'prediction': 1 if score >= THRESHOLD else 0
            })
    
    return residues


def parse_bp3_fasta(fasta_file, start_resid=1):
    """Parse FASTA de BepiPred-3.0 (MAYUSCULA=epitope)."""
    residues = []
    seq = ""
    with open(fasta_file) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    
    for i, aa in enumerate(seq):
        residues.append({
            'res_id': start_resid + i,
            'residue': aa.upper(),
            'bepipred3_score': None,
            'prediction': 1 if aa.isupper() else 0
        })
    return residues


def main():
    if len(sys.argv) < 3:
        print("Uso: python3 parse_bepipred3.py <bp3_output_dir> <output_csv> [--start_resid N]")
        sys.exit(1)
    
    bp3_dir = sys.argv[1]
    output_csv = sys.argv[2]
    
    # Parse --start_resid
    start_resid = 1
    if '--start_resid' in sys.argv:
        idx = sys.argv.index('--start_resid')
        start_resid = int(sys.argv[idx + 1])
    
    # Try CSV first (has scores)
    csv_files = glob.glob(os.path.join(bp3_dir, "raw_output*.csv"))
    if not csv_files:
        csv_files = glob.glob(os.path.join(bp3_dir, "*.csv"))
        csv_files = [f for f in csv_files if 'parsed' not in f]  # skip our own output
    
    if csv_files:
        print(f"Parseando CSV: {csv_files[0]}")
        residues = parse_bp3_csv(csv_files[0], start_resid)
    else:
        # Fall back to FASTA
        fasta_files = glob.glob(os.path.join(bp3_dir, "Bcell_epitope_preds.fasta"))
        if not fasta_files:
            fasta_files = glob.glob(os.path.join(bp3_dir, "*.fasta"))
        if not fasta_files:
            print(f"ERROR: No se encontraron archivos de salida en {bp3_dir}")
            sys.exit(1)
        print(f"Parseando FASTA: {fasta_files[0]}")
        residues = parse_bp3_fasta(fasta_files[0], start_resid)
    
    # Write output
    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['res_id', 'residue', 'bepipred3_score', 'prediction'])
        writer.writeheader()
        writer.writerows(residues)
    
    n_epitope = sum(1 for r in residues if r['prediction'] == 1)
    print(f"\nEscrito: {output_csv}")
    print(f"  Total residuos: {len(residues)}")
    print(f"  Rango PDB: {residues[0]['res_id']}–{residues[-1]['res_id']}")
    print(f"  Predichos epitope (score >= {THRESHOLD}): {n_epitope} ({100*n_epitope/len(residues):.1f}%)")


if __name__ == "__main__":
    main()
