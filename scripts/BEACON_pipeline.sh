#!/bin/bash
# ============================================================================
# BEACON_pipeline.sh — Pipeline automatizado de prediccion de epitopes B
# ============================================================================
# Ejecuta los tres predictores (BepiPred-3.0, DiscoTope-3.0, ElliPro),
# calcula el consenso por voto mayoritario, y valida contra epitopes IEDB.
#
# Uso:
#   bash BEACON_pipeline.sh --protein ag85b --pdb 1F0N --chain A
#   bash BEACON_pipeline.sh --protein spike_rbd --pdb 6W41 --chain C
#   bash BEACON_pipeline.sh --protein l1_hpv16 --pdb 1DZL --chain A
#   bash BEACON_pipeline.sh --all    # Ejecuta las 3 proteinas del estudio
#
# Requisitos:
#   - Miniconda con entornos: bepipred3, discotope3
#   - BepiPred-3.0 en $HOME/BEACON/tools/BepiPred-3.0/
#   - DiscoTope-3.0 en $HOME/BEACON/tools/DiscoTope-3.0/
#   - ElliPro.jar en $HOME/BEACON/tools/ElliPro.jar
#   - Java (openjdk-17-jre-headless)
#
# Manuel Alain Ramirez Saenz — Master en Bioinformatica, VIU 2025
# ============================================================================

set -euo pipefail

BEACON_DIR="${BEACON_DIR:-$HOME/BEACON}"
BP3_DIR="${BP3_DIR:-$BEACON_DIR/tools/BepiPred-3.0}"
DT3_DIR="${DT3_DIR:-$BEACON_DIR/tools/DiscoTope-3.0}"
ELLIPRO_JAR="${ELLIPRO_JAR:-$BEACON_DIR/tools/ElliPro.jar}"
SCRIPTS_DIR="$BEACON_DIR/scripts"
DATA_DIR="$BEACON_DIR/data"

# Colores para output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_step() { echo -e "${BLUE}[BEACON]${NC} $1"; }
log_ok()   { echo -e "${GREEN}[OK]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_err()  { echo -e "${RED}[ERROR]${NC} $1"; }

# ───────────────────────────────────────────────────────────────
# Funciones del pipeline
# ───────────────────────────────────────────────────────────────

run_bepipred3() {
    local fasta=$1
    local outdir=$2

    log_step "Etapa 1: BepiPred-3.0"
    mkdir -p "$outdir"

    eval "$(conda shell.bash hook)"
    conda activate bepipred3

    cd "$BP3_DIR"
    python3 bepipred3_CLI.py -i "$fasta" -o "$outdir" -pred mjv_pred

    conda deactivate
    cd "$BEACON_DIR"

    python3 "$SCRIPTS_DIR/parse_bepipred3.py" "$outdir" "$outdir/bepipred3_parsed.csv"
    log_ok "BepiPred-3.0 completado: $outdir/bepipred3_parsed.csv"
}

run_discotope3() {
    local pdb=$1
    local outdir=$2

    log_step "Etapa 2: DiscoTope-3.0"
    mkdir -p "$outdir"

    eval "$(conda shell.bash hook)"
    conda activate discotope3

    cd "$DT3_DIR"
    python3 discotope3/main.py \
        --pdb_or_zip_file "$pdb" \
        --out_dir "$outdir" \
        --cpu_only

    conda deactivate
    cd "$BEACON_DIR"

    python3 "$SCRIPTS_DIR/parse_discotope3.py" "$outdir" "$outdir/discotope3_parsed.csv"
    log_ok "DiscoTope-3.0 completado: $outdir/discotope3_parsed.csv"
}

run_ellipro() {
    local pdb=$1
    local chain=$2
    local outdir=$3

    log_step "Etapa 3: ElliPro (standalone)"
    mkdir -p "$outdir"

    if [ ! -f "$ELLIPRO_JAR" ]; then
        log_err "No se encuentra $ELLIPRO_JAR"
        log_err "Descargar de https://tools.iedb.org/ellipro/download/"
        return 1
    fi

    java -jar "$ELLIPRO_JAR" \
        --input-file "$pdb" \
        --chains "$chain" \
        --min-score 0.5 \
        --max-dist 6 \
        --output "$outdir/ellipro_output.txt" \
        --table "$outdir/ellipro_table.txt"

    python3 "$SCRIPTS_DIR/parse_ellipro.py" "$outdir" "$outdir/ellipro_parsed.csv"
    log_ok "ElliPro completado: $outdir/ellipro_parsed.csv"
}

run_consensus() {
    local results_dir=$1

    log_step "Etapa 4: Consenso por voto mayoritario"
    python3 "$SCRIPTS_DIR/consensus.py" "$results_dir"
    log_ok "Consenso completado: $results_dir/consensus/"
}

run_validation() {
    local results_dir=$1

    log_step "Etapa 5: Validacion contra IEDB"

    local consensus_csv="$results_dir/consensus/consensus_per_residue.csv"
    local epitope_csv="$results_dir/validation/iedb_epitopes.csv"

    if [ ! -f "$epitope_csv" ]; then
        log_warn "No se encontro $epitope_csv — omitiendo validacion"
        return 1
    fi

    python3 "$SCRIPTS_DIR/validate_iedb.py" "$consensus_csv" "$epitope_csv"
    log_ok "Validacion completada: $results_dir/validation/"
}

run_visualization() {
    local results_dir=$1
    local protein_name=$2

    log_step "Etapa 6: Visualizacion"
    python3 "$SCRIPTS_DIR/visualize_results.py" "$results_dir" "$protein_name"
    log_ok "Figuras generadas: $results_dir/figures/"
}

# ───────────────────────────────────────────────────────────────
# Pipeline completo para una proteina
# ───────────────────────────────────────────────────────────────

run_pipeline() {
    local protein=$1
    local pdb_id=$2
    local chain=$3
    local protein_label=$4
    local results_dir="$BEACON_DIR/results/$protein"

    echo ""
    echo "============================================================"
    echo "  BEACON Pipeline — $protein_label (PDB: $pdb_id, chain $chain)"
    echo "============================================================"
    echo ""

    # Buscar archivos preparados
    local pdb_clean fasta
    pdb_clean=$(find "$DATA_DIR" -name "${pdb_id}_${chain}*clean.pdb" -o -name "${pdb_id}_${chain}*RBD.pdb" 2>/dev/null | head -1)
    fasta=$(find "$DATA_DIR" -name "${pdb_id}_${chain}*.fasta" 2>/dev/null | head -1)

    if [ -z "$pdb_clean" ] || [ -z "$fasta" ]; then
        log_err "No se encontraron archivos preparados para ${pdb_id}_${chain} en $DATA_DIR"
        log_err "Ejecuta primero: python scripts/prepare_structures.py --pdb $pdb_id --chain $chain"
        return 1
    fi

    echo "  PDB: $pdb_clean"
    echo "  FASTA: $fasta"
    echo ""

    run_bepipred3 "$fasta" "$results_dir/bepipred3"
    run_discotope3 "$pdb_clean" "$results_dir/discotope3"
    run_ellipro "$pdb_clean" "$chain" "$results_dir/ellipro"
    run_consensus "$results_dir"
    run_validation "$results_dir"
    run_visualization "$results_dir" "$protein_label"

    echo ""
    log_ok "Pipeline completo para $protein_label"
    echo ""
}

# ───────────────────────────────────────────────────────────────
# Main
# ───────────────────────────────────────────────────────────────

usage() {
    echo "BEACON Pipeline — Prediccion de epitopes B por consenso"
    echo ""
    echo "Uso:"
    echo "  bash BEACON_pipeline.sh --protein <nombre> --pdb <id> --chain <cadena>"
    echo "  bash BEACON_pipeline.sh --all"
    echo ""
    echo "Proteinas predefinidas (--all):"
    echo "  ag85b      PDB: 1F0N, cadena A  (M. tuberculosis)"
    echo "  spike_rbd  PDB: 6W41, cadena C  (SARS-CoV-2)"
    echo "  l1_hpv16   PDB: 1DZL, cadena A  (HPV-16)"
}

PROTEIN=""
PDB_ID=""
CHAIN=""
RUN_ALL=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --protein) PROTEIN="$2"; shift 2 ;;
        --pdb) PDB_ID="$2"; shift 2 ;;
        --chain) CHAIN="$2"; shift 2 ;;
        --all) RUN_ALL=true; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "Opcion desconocida: $1"; usage; exit 1 ;;
    esac
done

echo "============================================================"
echo "  BEACON v1.0 — B-cell Epitope Analysis CONsensus"
echo "  Pipeline de prediccion por voto mayoritario"
echo "============================================================"
echo "  Directorio: $BEACON_DIR"
echo "  Fecha: $(date)"
echo ""

if $RUN_ALL; then
    log_step "Ejecutando pipeline para las 3 proteinas del estudio..."
    run_pipeline "ag85b" "1F0N" "A" "Ag85B (M. tuberculosis)"
    run_pipeline "spike_rbd" "6W41" "C" "Spike RBD (SARS-CoV-2)"
    run_pipeline "l1_hpv16" "1DZL" "A" "L1 HPV-16"

    echo ""
    echo "============================================================"
    log_ok "BEACON Pipeline completado para las 3 proteinas"
    echo "============================================================"
    echo "Resultados en: $BEACON_DIR/results/"

elif [ -n "$PROTEIN" ] && [ -n "$PDB_ID" ] && [ -n "$CHAIN" ]; then
    run_pipeline "$PROTEIN" "$PDB_ID" "$CHAIN" "$PROTEIN"
else
    log_err "Especifica --all o --protein/--pdb/--chain"
    usage
    exit 1
fi
