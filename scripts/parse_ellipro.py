#!/usr/bin/env python3
"""
BEACON - Parsear salida de ElliPro al formato unificado por residuo.

ElliPro standalone genera:
  - table.txt: protrusion index (PI) por residuo (CSV: No.,Chain ID,Residue Number,Residue Name,Score)
  - output.txt: lista de epitopes lineales y discontinuos

Este script lee table.txt (tiene PI por residuo) y genera:
  res_id, residue, ellipro_score, prediction

Uso:
    python3 parse_ellipro.py <ellipro_dir> <output_csv>
    python3 parse_ellipro.py results/ag85b/ellipro/ results/ag85b/ellipro/ellipro_parsed.csv
"""

import sys
import os
import csv
import glob

AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def parse_table_csv(filepath):
    """
    Parse ElliPro standalone table.txt — CSV with protrusion index per residue.
    Format: No.,Chain ID,Residue Number,Residue Name,Score
    """
    residues = {}
    with open(filepath) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                res_id = int(row.get('Residue Number', row.get('res_id', 0)))
                res_name = row.get('Residue Name', row.get('residue', '?'))
                pi = float(row.get('Score', row.get('score', 0)))
                aa = AA3TO1.get(res_name.upper(), res_name[0] if len(res_name) == 1 else '?')
                if res_id not in residues or pi > residues[res_id][1]:
                    residues[res_id] = (aa, pi)
            except (ValueError, KeyError):
                continue
    return residues


def parse_output_epitopes(filepath):
    """
    Parse ElliPro standalone output.txt — epitope list.
    CSV format with headers like: No.,Structure,Chain,Start Position,End Position,...,Score,Type
    and discontinuous: No.,Structure,Residues,Number of Residues,Score,Type
    """
    residues = {}
    with open(filepath) as f:
        content = f.read()

    # Split into sections (linear and discontinuous)
    sections = content.strip().split('\n\n')
    for section in sections:
        lines = section.strip().split('\n')
        if not lines:
            continue
        header = lines[0]

        if 'Start Position' in header:
            # Linear epitopes
            reader = csv.DictReader(lines)
            for row in reader:
                try:
                    start = int(row['Start Position'])
                    end = int(row['End Position'])
                    score = float(row['Score'])
                    peptide = row.get('Peptide', '')
                    for i, pos in enumerate(range(start, end + 1)):
                        aa = peptide[i] if i < len(peptide) else '?'
                        if pos not in residues or score > residues[pos][1]:
                            residues[pos] = (aa, score)
                except (ValueError, KeyError):
                    continue

        elif 'Residues' in header:
            # Discontinuous epitopes
            reader = csv.DictReader(lines)
            for row in reader:
                try:
                    score = float(row['Score'])
                    res_text = row['Residues']
                    import re
                    for match in re.finditer(r'[A-Z]:([A-Z])(\d+)', res_text):
                        aa = match.group(1)
                        pos = int(match.group(2))
                        if pos not in residues or score > residues[pos][1]:
                            residues[pos] = (aa, score)
                except (ValueError, KeyError):
                    continue

    return residues


def main():
    if len(sys.argv) < 3:
        print("Uso: python3 parse_ellipro.py <ellipro_dir> <output_csv>")
        sys.exit(1)

    ellipro_dir = sys.argv[1]
    output_csv = sys.argv[2]
    threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5

    all_scores = {}

    # Priority 1: table.txt (PI per residue — most complete)
    table_file = os.path.join(ellipro_dir, "ellipro_table.txt")
    if not os.path.exists(table_file):
        table_file = os.path.join(ellipro_dir, "table.txt")
    if os.path.exists(table_file):
        print(f"  Parseando tabla standalone: {table_file}")
        all_scores = parse_table_csv(table_file)
        print(f"  Residuos con PI: {len(all_scores)}")

    # Priority 2: output.txt (epitope list — use if no table)
    if not all_scores:
        output_file = os.path.join(ellipro_dir, "ellipro_output.txt")
        if not os.path.exists(output_file):
            output_file = os.path.join(ellipro_dir, "output.txt")
        if os.path.exists(output_file):
            print(f"  Parseando output standalone: {output_file}")
            all_scores = parse_output_epitopes(output_file)
            print(f"  Residuos mapeados: {len(all_scores)}")

    if not all_scores:
        print(f"ERROR: No se encontraron datos de ElliPro en {ellipro_dir}")
        print("Archivos esperados: ellipro_table.txt o ellipro_output.txt")
        sys.exit(1)

    # Write unified output
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['res_id', 'residue', 'ellipro_score', 'prediction'])
        for pos in sorted(all_scores.keys()):
            aa, score = all_scores[pos]
            pred = 1 if score >= threshold else 0
            writer.writerow([pos, aa, f'{score:.4f}', pred])

    n_total = len(all_scores)
    n_epitope = sum(1 for pos in all_scores if all_scores[pos][1] >= threshold)
    print(f"\nEscrito: {output_csv}")
    print(f"  Residuos con score: {n_total}")
    print(f"  Rango PDB: {min(all_scores.keys())}–{max(all_scores.keys())}")
    print(f"  Predichos epitope (PI >= {threshold}): {n_epitope} ({100*n_epitope/n_total:.1f}%)")


if __name__ == "__main__":
    main()
