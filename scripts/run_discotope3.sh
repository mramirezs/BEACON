#!/bin/bash
# BEACON - Ejecutar BepiPred-3.0 para una proteina
# Uso: bash run_bepipred3.sh <fasta_file> <output_dir>
# Ejemplo: bash run_bepipred3.sh data/1F0N_A.fasta results/ag85b/bepipred3/

set -euo pipefail

FASTA="${1:?Uso: bash run_bepipred3.sh <fasta> <output_dir>}"
OUTDIR="${2:?Especifica directorio de salida}"
BP3_DIR="${BP3_DIR:-$HOME/BEACON/tools/BepiPred-3.0}"

if [ ! -f "$FASTA" ]; then
    echo "ERROR: No se encuentra $FASTA"
    exit 1
fi

mkdir -p "$OUTDIR"

echo "==========================================="
echo "BEACON - BepiPred-3.0"
echo "==========================================="
echo "Entrada:  $FASTA"
echo "Salida:   $OUTDIR"
echo ""

# Activar entorno conda
eval "$(conda shell.bash hook)"
conda activate bepipred3

# Ejecutar BepiPred-3.0
cd "$BP3_DIR"
python3 bepipred3_CLI.py -i "$FASTA" -o "$OUTDIR" -pred mjv_pred

echo ""
echo "BepiPred-3.0 completado."
echo "Resultados en: $OUTDIR"
ls -la "$OUTDIR"
