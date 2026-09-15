#!/usr/bin/env python3
"""
BEACON - Consenso por voto mayoritario de los tres predictores.

Integra las predicciones parseadas de BepiPred-3.0, ElliPro y DiscoTope-3.0
en un score de consenso por residuo (0-3), donde cada punto representa el
numero de predictores que clasifican ese residuo como epitope.

Genera dos archivos de salida:
  - consensus_per_residue.csv : score de consenso residuo por residuo
  - consensus_regions.csv     : regiones contiguas con consenso >= 2

Uso:
    python3 consensus.py <results_dir> [--min_consensus 2]
    python3 consensus.py results/ag85b/

    Espera encontrar en <results_dir>:
      bepipred3/bepipred3_parsed.csv
      discotope3/discotope3_parsed.csv
      ellipro/ellipro_parsed.csv

    Escribe la salida en:
      <results_dir>/consensus/consensus_per_residue.csv
      <results_dir>/consensus/consensus_regions.csv
"""

import sys
import os
import csv
import glob


MIN_CONSENSUS = 2   # umbral por defecto para formar region (>= 2/3)
MIN_REGION_LEN = 3  # longitud minima de una region (residuos consecutivos)

# Tabla de conversion de 3 letras a 1 letra (para secuencias de region)
AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def find_parsed(results_dir, predictor):
    """Localiza el CSV parseado de un predictor."""
    patterns = [
        os.path.join(results_dir, predictor, f"{predictor}_parsed.csv"),
        os.path.join(results_dir, predictor, "*_parsed.csv"),
    ]
    for pat in patterns:
        files = glob.glob(pat)
        if files:
            return files[0]
    return None


def load_predictor(csv_file, score_col):
    """Carga un CSV parseado: devuelve dict res_id -> (residue, score, prediction)."""
    data = {}
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rid = int(row['res_id'])
            residue = row.get('residue', '?')
            score = row.get(score_col, '')
            pred = int(row['prediction'])
            data[rid] = (residue, score, pred)
    return data


def one_letter(res):
    """Normaliza un residuo a codigo de 1 letra."""
    res = res.strip()
    if len(res) == 1:
        return res.upper()
    return AA3TO1.get(res.upper(), 'X')


def build_consensus(bp3, dt3, elli):
    """
    Combina los tres predictores en un score de consenso por residuo.
    Devuelve una lista de dicts ordenada por res_id.
    """
    all_ids = sorted(set(bp3) | set(dt3) | set(elli))
    rows = []
    for rid in all_ids:
        # Residuo: preferir el de cualquier predictor que lo tenga
        residue = '?'
        for src in (bp3, dt3, elli):
            if rid in src:
                residue = one_letter(src[rid][0])
                break

        bp3_pred = bp3[rid][2] if rid in bp3 else 0
        dt3_pred = dt3[rid][2] if rid in dt3 else 0
        elli_pred = elli[rid][2] if rid in elli else 0

        bp3_score = bp3[rid][1] if rid in bp3 else ''
        dt3_score = dt3[rid][1] if rid in dt3 else ''
        elli_score = elli[rid][1] if rid in elli else ''

        consensus = bp3_pred + dt3_pred + elli_pred

        rows.append({
            'res_id': rid,
            'residue': residue,
            'bepipred3_score': bp3_score,
            'discotope3_score': dt3_score,
            'ellipro_score': elli_score,
            'bp3_pred': bp3_pred,
            'dt3_pred': dt3_pred,
            'elli_pred': elli_pred,
            'consensus_score': consensus,
        })
    return rows


def find_regions(consensus_rows, min_consensus=MIN_CONSENSUS, min_length=MIN_REGION_LEN):
    """
    Identifica regiones contiguas donde consensus_score >= min_consensus.
    Solo se conservan regiones con al menos min_length residuos consecutivos.
    Devuelve lista de regiones con estadisticas.
    """
    regions = []
    current = []

    for row in consensus_rows:
        if row['consensus_score'] >= min_consensus:
            current.append(row)
        else:
            if len(current) >= min_length:
                regions.append(_summarize_region(current))
            current = []
    if len(current) >= min_length:
        regions.append(_summarize_region(current))

    return regions


def _summarize_region(residues):
    """Calcula estadisticas de una region contigua."""
    start = residues[0]['res_id']
    end = residues[-1]['res_id']
    length = end - start + 1
    sequence = ''.join(r['residue'] for r in residues)
    scores = [r['consensus_score'] for r in residues]
    avg_c = round(sum(scores) / len(scores), 2)
    max_c = max(scores)
    bp3_pos = sum(r['bp3_pred'] for r in residues)
    dt3_pos = sum(r['dt3_pred'] for r in residues)
    elli_pos = sum(r['elli_pred'] for r in residues)
    return {
        'start': start, 'end': end, 'length': length,
        'sequence': sequence, 'avg_consensus': avg_c, 'max_consensus': max_c,
        'bp3_positive': bp3_pos, 'dt3_positive': dt3_pos, 'elli_positive': elli_pos,
    }


def main():
    if len(sys.argv) < 2:
        print("Uso: python3 consensus.py <results_dir> [--min_consensus N]")
        sys.exit(1)

    results_dir = sys.argv[1]
    min_consensus = MIN_CONSENSUS
    if '--min_consensus' in sys.argv:
        idx = sys.argv.index('--min_consensus')
        min_consensus = int(sys.argv[idx + 1])

    # Localizar los tres CSV parseados
    bp3_file = find_parsed(results_dir, 'bepipred3')
    dt3_file = find_parsed(results_dir, 'discotope3')
    elli_file = find_parsed(results_dir, 'ellipro')

    for name, f in [('BepiPred-3.0', bp3_file), ('DiscoTope-3.0', dt3_file), ('ElliPro', elli_file)]:
        if f is None:
            print(f"ERROR: No se encontro el CSV parseado de {name} en {results_dir}")
            sys.exit(1)
        print(f"  {name}: {f}")

    # Cargar predicciones
    bp3 = load_predictor(bp3_file, 'bepipred3_score')
    dt3 = load_predictor(dt3_file, 'discotope3_score')
    elli = load_predictor(elli_file, 'ellipro_score')

    # Construir consenso
    consensus_rows = build_consensus(bp3, dt3, elli)

    # Crear carpeta de salida
    out_dir = os.path.join(results_dir, 'consensus')
    os.makedirs(out_dir, exist_ok=True)

    # Escribir consensus_per_residue.csv
    per_res_file = os.path.join(out_dir, 'consensus_per_residue.csv')
    fields = ['res_id', 'residue', 'bepipred3_score', 'discotope3_score',
              'ellipro_score', 'bp3_pred', 'dt3_pred', 'elli_pred', 'consensus_score']
    with open(per_res_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(consensus_rows)

    # Identificar y escribir regiones
    regions = find_regions(consensus_rows, min_consensus)
    reg_file = os.path.join(out_dir, 'consensus_regions.csv')
    reg_fields = ['start', 'end', 'length', 'sequence', 'avg_consensus',
                  'max_consensus', 'bp3_positive', 'dt3_positive', 'elli_positive']
    with open(reg_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=reg_fields)
        writer.writeheader()
        writer.writerows(regions)

    # Resumen
    n = len(consensus_rows)
    c3 = sum(1 for r in consensus_rows if r['consensus_score'] == 3)
    c2plus = sum(1 for r in consensus_rows if r['consensus_score'] >= 2)
    c1plus = sum(1 for r in consensus_rows if r['consensus_score'] >= 1)

    print(f"\nEscrito: {per_res_file}")
    print(f"Escrito: {reg_file}")
    print(f"\n  Total residuos: {n}")
    print(f"  Consenso >= 1/3: {c1plus} ({100*c1plus/n:.1f}%)")
    print(f"  Consenso >= 2/3: {c2plus} ({100*c2plus/n:.1f}%)")
    print(f"  Consenso  = 3/3: {c3} ({100*c3/n:.1f}%)")
    print(f"  Regiones de consenso (>= {min_consensus}/3): {len(regions)}")


if __name__ == "__main__":
    main()
