# Guía de Instalación y Ejecución de BEACON en Kubuntu
## Pipeline de Predicción de Epítopes B por Consenso (Voto Mayoritario)

**Sistema operativo:** Kubuntu (basado en Ubuntu)
**Requisitos mínimos:** 16 GB RAM, 2+ cores CPU, ~10 GB disco libre
**Autor:** Manuel Alain Ramírez Sáenz — Máster en Bioinformática, VIU 2025

---

## 1. Prerrequisitos del sistema

```bash
sudo apt update && sudo apt install -y \
    build-essential gcc g++ make cmake \
    git wget curl unzip \
    python3-dev libffi-dev libssl-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    python-is-python3
```

## 2. Instalar Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# Seguir instrucciones, aceptar licencia, activar conda init
source ~/.bashrc
```

## 3. Crear estructura de directorios

```bash
mkdir -p ~/BEACON/{data,results/{ag85b,spike_rbd,l1_hpv16}/{bepipred3,discotope3,ellipro,consensus,validation,figures},scripts,tools}
cd ~/BEACON
```

## 4. Instalar BepiPred-3.0

```bash
# Crear entorno conda
conda create -n bepipred3 python=3.10 -y
conda activate bepipred3
conda install pytorch cpuonly -c pytorch -y
pip install fair-esm biopython

# Clonar repositorio
cd ~/BEACON/tools
git clone https://github.com/UberClifford/BepiPred-3.0.git
cd BepiPred-3.0

# Verificar instalación
python bepipred3_CLI.py --help
# Opciones de -pred: mjv_pred (majority vote, recomendado), vt_pred (variable threshold)
conda deactivate
```

> **Nota:** La primera ejecución descargará el modelo ESM-2 (~2.5 GB).
> Esto es normal y solo ocurre una vez.

## 5. Instalar DiscoTope-3.0

```bash
# Crear entorno conda
conda create -n discotope3 python=3.11 -y
conda activate discotope3

# Clonar repositorio
cd ~/BEACON/tools
git clone https://github.com/Magnushhoie/DiscoTope-3.0.git
cd DiscoTope-3.0

# Instalar dependencias y paquete
pip install -r requirements.txt
pip install .

# Descomprimir modelos (IMPORTANTE)
unzip models.zip

# Verificar instalación
python discotope3/main.py --help
conda deactivate
```

> **Nota:** DiscoTope-3.0 usa representaciones latentes ESM-IF1 (inverse folding)
> para predecir epítopes conformacionales. La primera ejecución descargará el modelo.

## 6. Instalar ElliPro (Standalone JAR)

```bash
# Instalar Java (si no está instalado)
sudo apt install openjdk-17-jre-headless -y

# Descargar ElliPro standalone
cd ~/BEACON/tools
wget https://tools.iedb.org/ellipro/download/ElliPro.jar
# (o copiar manualmente desde la pagina de descarga)

# Verificar instalación
java -jar ElliPro.jar --help
```

> **Nota:** ElliPro usa el protrusion index (PI) de Thornton et al.
> para identificar regiones expuestas en la superficie proteica.

## 7. Copiar scripts del pipeline

Copiar los scripts generados al directorio `~/BEACON/scripts/`:

```bash
cd ~/BEACON/scripts

# Scripts de parseo (copiar desde los artifacts del proyecto)
# parse_bepipred3.py    — parsea salida de BepiPred-3.0
# parse_discotope3.py   — parsea salida de DiscoTope-3.0
# parse_ellipro.py      — parsea resultados web de ElliPro
# consensus.py          — calcula consenso por voto mayoritario
# validate_iedb.py      — valida contra epítopes experimentales IEDB
# visualize_results.py  — genera heatmap y barplot de métricas
# BEACON_pipeline.sh    — script maestro del pipeline
# prepare_structures.py — prepara estructuras PDB (limpieza + FASTA)
```

## 8. Ejecución paso a paso — Ag85B (M. tuberculosis)

### 8.1 Preparar estructura

```bash
cd ~/BEACON

# Descargar PDB 1F0N desde RCSB
wget -O data/1F0N.pdb https://files.rcsb.org/download/1F0N.pdb

# Preparar (limpiar cadena A, extraer FASTA)
python3 scripts/prepare_structures.py --pdb 1F0N --chain A --label ag85b
```

Archivos generados:
- `data/1F0N_A_clean.pdb` — estructura limpia (284 residuos, pos 2–285)
- `data/1F0N_A.fasta` — secuencia en formato FASTA

### 8.2 BepiPred-3.0

> **Antes de ejecutar:** verificar que el header del FASTA no contenga
> caracteres especiales (`|`, `:`). Si los tiene, simplificar:
> `sed -i '1s/.*/>(ID)/' archivo.fasta`

```bash
conda activate bepipred3
cd ~/BEACON/tools/BepiPred-3.0

# IMPORTANTE: usar rutas ABSOLUTAS y flag -pred mjv_pred
python bepipred3_CLI.py \
    -i ~/BEACON/data/1F0N_A.fasta \
    -o ~/BEACON/results/ag85b/bepipred3/ \
    -pred mjv_pred

conda deactivate
cd ~/BEACON

# Parsear al formato BEACON
python3 scripts/parse_bepipred3.py \
    results/ag85b/bepipred3/ \
    results/ag85b/bepipred3/bepipred3_parsed.csv
```

### 8.3 DiscoTope-3.0

```bash
conda activate discotope3
cd ~/BEACON/tools/DiscoTope-3.0

python discotope3/main.py \
    --pdb_or_zip_file ~/BEACON/data/1F0N_A_clean.pdb \
    --out_dir ~/BEACON/results/ag85b/discotope3/ \
    --cpu_only

conda deactivate
cd ~/BEACON

# Parsear al formato BEACON
python3 scripts/parse_discotope3.py \
    results/ag85b/discotope3/ \
    results/ag85b/discotope3/discotope3_parsed.csv
```

### 8.4 ElliPro

```bash
mkdir -p results/ag85b/ellipro

java -jar ~/BEACON/tools/ElliPro.jar \
    --input-file data/1F0N_A_clean.pdb \
    --chains A \
    --min-score 0.5 \
    --max-dist 6 \
    --output results/ag85b/ellipro/ellipro_output.txt \
    --table results/ag85b/ellipro/ellipro_table.txt

# Parsear al formato BEACON
python3 scripts/parse_ellipro.py \
    results/ag85b/ellipro/ \
    results/ag85b/ellipro/ellipro_parsed.csv
```

### 8.5 Consenso

```bash
python3 scripts/consensus.py results/ag85b/
```

### 8.6 Validación contra IEDB

```bash
python3 scripts/validate_iedb.py \
    results/ag85b/consensus/consensus_per_residue.csv \
    results/ag85b/validation/iedb_epitopes.csv
```

### 8.7 Visualización

```bash
python3 scripts/visualize_results.py \
    results/ag85b/ \
    "Ag85B (M. tuberculosis)"
```

## 9. Ejecución — Spike RBD (SARS-CoV-2, PDB: 6W41)

> **IMPORTANTE:** 6W41 es un complejo Fab-RBD. La cadena del Spike RBD
> es la **cadena C** (residuos 333–527), no la cadena A.

```bash
# Preparar
wget -O data/6W41.pdb https://files.rcsb.org/download/6W41.pdb
python3 scripts/prepare_structures.py --pdb 6W41 --chain C --label spike_rbd

# BepiPred-3.0
conda activate bepipred3
cd ~/BEACON/tools/BepiPred-3.0
python bepipred3_CLI.py -i ~/BEACON/data/6W41_C_RBD.fasta -o ~/BEACON/results/spike_rbd/bepipred3/ -pred mjv_pred
conda deactivate
cd ~/BEACON
python3 scripts/parse_bepipred3.py results/spike_rbd/bepipred3/ results/spike_rbd/bepipred3/bepipred3_parsed.csv

# DiscoTope-3.0
conda activate discotope3
cd ~/BEACON/tools/DiscoTope-3.0
python discotope3/main.py --pdb_or_zip_file ~/BEACON/data/6W41_C_RBD.pdb --out_dir ~/BEACON/results/spike_rbd/discotope3/ --cpu_only
conda deactivate
cd ~/BEACON
python3 scripts/parse_discotope3.py results/spike_rbd/discotope3/ results/spike_rbd/discotope3/discotope3_parsed.csv

# ElliPro
mkdir -p results/spike_rbd/ellipro
java -jar ~/BEACON/tools/ElliPro.jar --input-file data/6W41_C_RBD.pdb --chains C --min-score 0.5 --max-dist 6 --output results/spike_rbd/ellipro/ellipro_output.txt --table results/spike_rbd/ellipro/ellipro_table.txt
python3 scripts/parse_ellipro.py results/spike_rbd/ellipro/ results/spike_rbd/ellipro/ellipro_parsed.csv

# Consenso + Validación + Figuras
python3 scripts/consensus.py results/spike_rbd/
python3 scripts/validate_iedb.py results/spike_rbd/consensus/consensus_per_residue.csv results/spike_rbd/validation/iedb_epitopes.csv
python3 scripts/visualize_results.py results/spike_rbd/ "Spike RBD (SARS-CoV-2)"
```

## 10. Ejecución — L1 HPV-16 (PDB: 1DZL)

```bash
# Preparar
wget -O data/1DZL.pdb https://files.rcsb.org/download/1DZL.pdb
python3 scripts/prepare_structures.py --pdb 1DZL --chain A --label l1_hpv16

# BepiPred-3.0
conda activate bepipred3
cd ~/BEACON/tools/BepiPred-3.0
python bepipred3_CLI.py -i ~/BEACON/data/1DZL_A.fasta -o ~/BEACON/results/l1_hpv16/bepipred3/ -pred mjv_pred
conda deactivate
cd ~/BEACON
python3 scripts/parse_bepipred3.py results/l1_hpv16/bepipred3/ results/l1_hpv16/bepipred3/bepipred3_parsed.csv

# DiscoTope-3.0
conda activate discotope3
cd ~/BEACON/tools/DiscoTope-3.0
python discotope3/main.py --pdb_or_zip_file ~/BEACON/data/1DZL_A_clean.pdb --out_dir ~/BEACON/results/l1_hpv16/discotope3/ --cpu_only
conda deactivate
cd ~/BEACON
python3 scripts/parse_discotope3.py results/l1_hpv16/discotope3/ results/l1_hpv16/discotope3/discotope3_parsed.csv

# ElliPro
mkdir -p results/l1_hpv16/ellipro
java -jar ~/BEACON/tools/ElliPro.jar --input-file data/1DZL_A_clean.pdb --chains A --min-score 0.5 --max-dist 6 --output results/l1_hpv16/ellipro/ellipro_output.txt --table results/l1_hpv16/ellipro/ellipro_table.txt
python3 scripts/parse_ellipro.py results/l1_hpv16/ellipro/ results/l1_hpv16/ellipro/ellipro_parsed.csv

# Consenso + Validación + Figuras
python3 scripts/consensus.py results/l1_hpv16/
python3 scripts/validate_iedb.py results/l1_hpv16/consensus/consensus_per_residue.csv results/l1_hpv16/validation/iedb_epitopes.csv
python3 scripts/visualize_results.py results/l1_hpv16/ "L1 HPV-16"
```

## 11. Ejecución automatizada (las 3 proteínas)

```bash
bash scripts/BEACON_pipeline.sh --all
```

O para una sola proteína:
```bash
bash scripts/BEACON_pipeline.sh --protein ag85b --pdb 1F0N --chain A
```

## 12. Estructura final del proyecto

```
~/BEACON/
├── scripts/
│   ├── BEACON_pipeline.sh          # Pipeline maestro
│   ├── prepare_structures.py       # Preparación de PDBs
│   ├── run_bepipred3.sh            # Wrapper BepiPred-3.0
│   ├── parse_bepipred3.py          # Parser BepiPred-3.0
│   ├── run_discotope3.sh           # Wrapper DiscoTope-3.0
│   ├── parse_discotope3.py         # Parser DiscoTope-3.0
│   ├── run_ellipro.sh              # Wrapper ElliPro standalone
│   ├── parse_ellipro.py            # Parser ElliPro
│   ├── consensus.py                # Consenso por voto mayoritario
│   ├── validate_iedb.py            # Validación contra IEDB
│   └── visualize_results.py        # Heatmap + barplot
├── data/
│   ├── 1F0N_A_clean.pdb            # Ag85B estructura limpia
│   ├── 1F0N_A.fasta                # Ag85B secuencia
│   ├── 6W41_C_RBD.pdb              # Spike RBD estructura limpia
│   ├── 6W41_C_RBD.fasta            # Spike RBD secuencia
│   ├── 1DZL_A_clean.pdb            # L1 HPV-16 estructura limpia
│   └── 1DZL_A.fasta                # L1 HPV-16 secuencia
├── results/
│   ├── ag85b/
│   │   ├── bepipred3/bepipred3_parsed.csv
│   │   ├── discotope3/discotope3_parsed.csv
│   │   ├── ellipro/ellipro_parsed.csv
│   │   ├── consensus/
│   │   │   ├── consensus_per_residue.csv
│   │   │   └── consensus_regions.csv
│   │   ├── validation/
│   │   │   ├── iedb_epitopes.csv
│   │   │   ├── iedb_residue_labels.csv
│   │   │   └── validation_metrics.csv
│   │   └── figures/
│   │       ├── heatmap_predictions.png
│   │       └── metrics_barplot.png
│   ├── spike_rbd/  (misma estructura)
│   └── l1_hpv16/   (misma estructura)
└── tools/
    ├── BepiPred-3.0/
    └── DiscoTope-3.0/
```

## Troubleshooting

| Problema | Solución |
|----------|----------|
| `ModuleNotFoundError: torch` | Verificar que el entorno conda correcto está activado |
| ESM-2 descarga lenta | Normal la primera vez (~2.5 GB). Usar wifi estable |
| DiscoTope `pybind11` error | `pip install pybind11` antes de instalar discotope3 |
| DiscoTope `biotite`/`torch_geometric` error | Requiere Python ≥ 3.11; recrear entorno con `python=3.11` |
| ElliPro timeout en web | Reducir tamaño del PDB, o reintentar más tarde |
| parse_bepipred3 no encuentra CSV | Verificar ruta de salida; BP3 puede generar solo FASTA |
| `CUDA not available` | Normal con `--cpu_only`; no es error, solo aviso |

## Tiempos estimados de ejecución (CPU)

| Proteína | BepiPred-3.0 | DiscoTope-3.0 | ElliPro (web) |
|----------|-------------|---------------|---------------|
| Ag85B (284 aa) | ~5 min | ~3 min | ~2 min |
| Spike RBD (195 aa) | ~3 min | ~2 min | ~1 min |
| L1 HPV-16 (455 aa) | ~8 min | ~5 min | ~3 min |
