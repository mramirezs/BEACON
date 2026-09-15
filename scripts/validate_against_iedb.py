#!/usr/bin/env python3
"""
BEACON - Validacion de las predicciones contra epitopes experimentales del IEDB.

Cruza el consenso por residuo (consensus_per_residue.csv) con las etiquetas
experimentales del IEDB (iedb_labels.csv) y calcula metricas de rendimiento
para cada predictor individual y para cada umbral de consenso.

Metricas por metodo: TP, FP, FN, TN, sensibilidad, especificidad, precision,
F1-score y coeficiente de correlacion de Matthews (MCC).

Uso:
    python3 validate_against_iedb.py <consensus_per_residue.csv> <iedb_labels.csv> <output.csv>
    python3 validate_against_iedb.py \\
        results/l1_hpv16/consensus/consensus_per_residue.csv \\
        results/l1_hpv16/validation/iedb_labels.csv \\
        results/l1_hpv16/validation/validation_metrics.csv

Formato esperado de iedb_labels.csv:
    res_id, residue, iedb_positive   (iedb_positive: 1=positivo, 0=negativo)

Nota metodologica: solo se incluyen en el calculo los residuos con etiqueta
experimental (positiva o negativa explicitamente ensayada). Los residuos sin
evidencia experimental deben marcarse aparte y excluirse del computo, ya que
la ausencia de evidencia no equivale a una etiqueta negativa.
"""

import sys
import csv
import math


def load_consensus(path):
    """Carga consensus_per_residue.csv: res_id -> dict de predicciones."""
    data = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            rid = int(row['res_id'])
            data[rid] = {
                'bp3': int(row['bp3_pred']),
                'dt3': int(row['dt3_pred']),
                'elli': int(row['elli_pred']),
                'consensus': int(row['consensus_score']),
            }
    return data


def load_iedb_labels(path):
    """Carga iedb_labels.csv: res_id -> iedb_positive (0/1)."""
    labels = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            rid = int(row['res_id'])
            labels[rid] = int(row['iedb_positive'])
    return labels


def compute_metrics(predictions, labels):
    """
    Calcula TP/FP/FN/TN y metricas derivadas.
    predictions: dict res_id -> 0/1
    labels: dict res_id -> 0/1
    Solo se evaluan los residuos presentes en labels.
    """
    tp = fp = fn = tn = 0
    for rid, true_label in labels.items():
        pred = predictions.get(rid, 0)
        if pred == 1 and true_label == 1:
            tp += 1
        elif pred == 1 and true_label == 0:
            fp += 1
        elif pred == 0 and true_label == 1:
            fn += 1
        else:
            tn += 1

    sens = tp / (tp + fn) if (tp + fn) > 0 else None
    spec = tn / (tn + fp) if (tn + fp) > 0 else None
    prec = tp / (tp + fp) if (tp + fp) > 0 else None
    f1 = (2 * prec * sens / (prec + sens)
          if prec and sens and (prec + sens) > 0 else None)

    # MCC
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / denom) if denom > 0 else None

    return {
        'TP': tp, 'FP': fp, 'FN': fn, 'TN': tn,
        'Sensitivity': sens, 'Specificity': spec,
        'Precision': prec, 'F1': f1, 'MCC': mcc,
    }


def main():
    if len(sys.argv) < 4:
        print("Uso: python3 validate_against_iedb.py <consensus_per_residue.csv> "
              "<iedb_labels.csv> <output.csv>")
        sys.exit(1)

    consensus_file = sys.argv[1]
    iedb_file = sys.argv[2]
    output_file = sys.argv[3]

    consensus = load_consensus(consensus_file)
    labels = load_iedb_labels(iedb_file)

    print(f"  Residuos con consenso: {len(consensus)}")
    print(f"  Residuos con etiqueta IEDB: {len(labels)}")
    pos = sum(1 for v in labels.values() if v == 1)
    neg = sum(1 for v in labels.values() if v == 0)
    print(f"  Positivos IEDB: {pos}  |  Negativos IEDB: {neg}")

    # Definir los metodos a evaluar
    metodos = {
        'BepiPred-3.0': {rid: d['bp3'] for rid, d in consensus.items()},
        'DiscoTope-3.0': {rid: d['dt3'] for rid, d in consensus.items()},
        'ElliPro': {rid: d['elli'] for rid, d in consensus.items()},
        'Score ≥ 1/3': {rid: (1 if d['consensus'] >= 1 else 0) for rid, d in consensus.items()},
        'Score ≥ 2/3': {rid: (1 if d['consensus'] >= 2 else 0) for rid, d in consensus.items()},
        'Score = 3/3': {rid: (1 if d['consensus'] == 3 else 0) for rid, d in consensus.items()},
    }

    # Calcular metricas para cada metodo
    fields = ['label', 'TP', 'FP', 'FN', 'TN', 'Sensitivity',
              'Specificity', 'Precision', 'F1', 'MCC']
    results = []
    for label, preds in metodos.items():
        m = compute_metrics(preds, labels)
        row = {'label': label}
        for k in ['TP', 'FP', 'FN', 'TN']:
            row[k] = m[k]
        for k in ['Sensitivity', 'Specificity', 'Precision', 'F1', 'MCC']:
            row[k] = '' if m[k] is None else m[k]
        results.append(row)

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)

    print(f"\nEscrito: {output_file}")
    print(f"\n{'Metodo':<16}{'TP':>5}{'FP':>5}{'FN':>5}{'TN':>5}{'MCC':>9}")
    print("-" * 45)
    for r in results:
        mcc = f"{float(r['MCC']):.3f}" if r['MCC'] != '' else '—'
        print(f"{r['label']:<16}{r['TP']:>5}{r['FP']:>5}{r['FN']:>5}{r['TN']:>5}{mcc:>9}")


if __name__ == "__main__":
    main()
