# BEACON

**B-cell Epitope Analysis with CONsensus**

Pipeline bioinformático para la predicción de epítopes de células B mediante consenso por voto mayoritario de tres predictores complementarios: **BepiPred-3.0**, **ElliPro** y **DiscoTope-3.0**.

BEACON integra las predicciones de los tres métodos, calcula un consenso por residuo (escala 0–3) y valida los resultados contra datos experimentales del IEDB. Está diseñado para **priorizar regiones candidatas a epítopo de alta precisión**, orientando la validación experimental posterior.

---

## Tabla de contenidos

- [¿Qué hace BEACON?](#qué-hace-beacon)
- [Requisitos del sistema](#requisitos-del-sistema)
- [Instalación desde cero](#instalación-desde-cero)
  - [1. Preparar el sistema](#1-preparar-el-sistema)
  - [2. Instalar Miniconda](#2-instalar-miniconda)
  - [3. Descargar BEACON](#3-descargar-beacon)
  - [4. Instalar BepiPred-3.0](#4-instalar-bepipred-30)
  - [5. Instalar DiscoTope-3.0](#5-instalar-discotope-30)
  - [6. Instalar ElliPro](#6-instalar-ellipro)
- [Uso](#uso)
- [Estructura del proyecto](#estructura-del-proyecto)
- [Proteínas analizadas](#proteínas-analizadas)
- [Solución de problemas](#solución-de-problemas)
- [Cita](#cita)
- [Licencia](#licencia)

---

## ¿Qué hace BEACON?

El pipeline ejecuta seis etapas para cada proteína:

1. **Preparación estructural** — descarga la estructura del PDB, elimina agua y ligandos, selecciona la cadena de interés y extrae la secuencia FASTA.
2. **BepiPred-3.0** — predicción basada en secuencia (modelo de lenguaje ESM-2).
3. **DiscoTope-3.0** — predicción estructural (representaciones de inverse folding ESM-IF1 + XGBoost).
4. **ElliPro** — predicción estructural por índice de protrusión.
5. **Consenso** — voto mayoritario por residuo: score de 0 a 3 según cuántos predictores coinciden.
6. **Validación** — cruce con epítopes experimentales del IEDB y cálculo de métricas.

Cada predictor binariza su score con un umbral estándar (BepiPred-3.0 ≥ 0.1496, ElliPro ≥ 0.5, DiscoTope-3.0 ≥ −3.7) y el consenso suma las tres predicciones, otorgando igual peso a cada una.

---

## Requisitos del sistema

- **Sistema operativo:** Linux (probado en Kubuntu/Ubuntu). En Windows, usar WSL2.
- **Memoria RAM:** mínimo 16 GB.
- **Procesador:** 2 o más núcleos.
- **Disco libre:** ~10 GB (los modelos de lenguaje ocupan varios GB).
- **Conexión a internet:** necesaria para descargar estructuras, modelos y datos del IEDB.

No se requiere GPU; el pipeline funciona con CPU.

---

## Instalación desde cero

> Estas instrucciones asumen que partes de un sistema Linux limpio, sin nada instalado. Copia y pega cada bloque de comandos en una terminal.

### 1. Preparar el sistema

Instala las herramientas de compilación y utilidades básicas:

```bash
sudo apt update && sudo apt install -y \
    build-essential gcc g++ make cmake \
    git wget curl unzip \
    python3-dev libffi-dev libssl-dev \
    zlib1g-dev libbz2-dev liblzma-dev \
    python-is-python3 \
    openjdk-17-jre-headless
```

> `openjdk-17-jre-headless` es Java, necesario para ElliPro.

### 2. Instalar Miniconda

Miniconda gestiona los entornos de Python aislados que necesita cada predictor:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Sigue las instrucciones en pantalla (acepta la licencia y responde `yes` a la inicialización). Luego recarga la terminal:

```bash
source ~/.bashrc
```

### 3. Descargar BEACON

Clona este repositorio y entra en la carpeta:

```bash
cd ~
git clone https://github.com/TU_USUARIO/BEACON.git
cd BEACON
```

> Reemplaza `TU_USUARIO` por tu nombre de usuario de GitHub.

Crea las carpetas de trabajo:

```bash
mkdir -p data results tools logs
```

### 4. Instalar BepiPred-3.0

```bash
# Crear el entorno
conda create -n bepipred3 python=3.10 -y
conda activate bepipred3
conda install pytorch cpuonly -c pytorch -y
pip install fair-esm biopython

# Descargar la herramienta
cd ~/BEACON/tools
git clone https://github.com/UberClifford/BepiPred-3.0.git

# Verificar
cd BepiPred-3.0
python bepipred3_CLI.py --help
conda deactivate
cd ~/BEACON
```

> **Nota:** la primera ejecución descargará el modelo ESM-2 (~2.5 GB). Es normal y ocurre una sola vez.

### 5. Instalar DiscoTope-3.0

```bash
# Crear el entorno (requiere Python 3.11)
conda create -n discotope3 python=3.11 -y
conda activate discotope3

# Descargar la herramienta
cd ~/BEACON/tools
git clone https://github.com/Magnushhoie/DiscoTope-3.0.git
cd DiscoTope-3.0

# Instalar dependencias y descomprimir los modelos
pip install -r requirements.txt
pip install .
unzip models.zip

# Verificar
python discotope3/main.py --help
conda deactivate
cd ~/BEACON
```

> **Nota:** DiscoTope-3.0 usa representaciones de inverse folding ESM-IF1. La primera ejecución descargará el modelo correspondiente.

### 6. Instalar ElliPro

ElliPro es un archivo Java (`.jar`) que se descarga desde el IEDB:

```bash
cd ~/BEACON/tools
wget https://tools.iedb.org/ellipro/download/ElliPro.jar

# Verificar
java -jar ElliPro.jar --help
```

> Si la descarga directa falla, descarga el `.jar` manualmente desde
> https://tools.iedb.org/ellipro/download/ y colócalo en `~/BEACON/tools/`.

---

## Uso

### Ejecutar las tres proteínas del estudio de una vez

```bash
bash scripts/BEACON_pipeline.sh --all
```

### Ejecutar una sola proteína

```bash
bash scripts/BEACON_pipeline.sh --protein ag85b --pdb 1F0N --chain A
bash scripts/BEACON_pipeline.sh --protein spike_rbd --pdb 6W41 --chain C
bash scripts/BEACON_pipeline.sh --protein l1_hpv16 --pdb 1DZL --chain A
```

### Ejecutar paso a paso (ejemplo con Ag85B)

```bash
# 1. Preparar estructura
python3 scripts/prepare_structures.py --pdb 1F0N --chain A --label ag85b

# 2. BepiPred-3.0
conda activate bepipred3
cd tools/BepiPred-3.0
python bepipred3_CLI.py -i ~/BEACON/data/1F0N_A.fasta \
    -o ~/BEACON/results/ag85b/bepipred3/ -pred mjv_pred
conda deactivate
cd ~/BEACON
python3 scripts/parse_bepipred3.py results/ag85b/bepipred3/ \
    results/ag85b/bepipred3/bepipred3_parsed.csv --start_resid 2

# 3. DiscoTope-3.0
conda activate discotope3
cd tools/DiscoTope-3.0
python discotope3/main.py --pdb_or_zip_file ~/BEACON/data/1F0N_A_clean.pdb \
    --out_dir ~/BEACON/results/ag85b/discotope3/ --cpu_only
conda deactivate
cd ~/BEACON
python3 scripts/parse_discotope3.py results/ag85b/discotope3/ \
    results/ag85b/discotope3/discotope3_parsed.csv

# 4. ElliPro
java -jar tools/ElliPro.jar --input-file data/1F0N_A_clean.pdb \
    --chains A --min-score 0.5 --max-dist 6 \
    --output results/ag85b/ellipro/ellipro_output.txt \
    --table results/ag85b/ellipro/ellipro_table.txt
python3 scripts/parse_ellipro.py results/ag85b/ellipro/ \
    results/ag85b/ellipro/ellipro_parsed.csv

# 5. Consenso
python3 scripts/consensus.py results/ag85b/
```

> **Importante:** el parámetro `--start_resid` en `parse_bepipred3.py` indica la posición del primer residuo según la numeración PDB. Es 2 para Ag85B, 333 para Spike RBD y 20 para L1 HPV-16.

---

## Estructura del proyecto

```
BEACON/
├── scripts/
│   ├── BEACON_pipeline.sh       # Pipeline maestro
│   ├── prepare_structures.py    # Descarga y limpieza de estructuras PDB
│   ├── parse_bepipred3.py       # Parser de BepiPred-3.0
│   ├── parse_discotope3.py      # Parser de DiscoTope-3.0
│   ├── parse_ellipro.py         # Parser de ElliPro
│   ├── consensus.py             # Consenso por voto mayoritario
│   └── run_discotope3.sh        # Wrapper de DiscoTope-3.0
├── data/                        # Estructuras PDB y secuencias FASTA
├── results/                     # Resultados por proteína
│   └── <proteína>/
│       ├── bepipred3/
│       ├── discotope3/
│       ├── ellipro/
│       └── consensus/
├── tools/                       # Predictores externos
│   ├── BepiPred-3.0/
│   ├── DiscoTope-3.0/
│   └── ElliPro.jar
└── README.md
```

---

## Proteínas analizadas

El estudio de referencia aplicó BEACON a tres antígenos de patógenos de relevancia clínica global y un control negativo:

| Proteína | Organismo | PDB | Cadena | Residuos | Resolución |
|----------|-----------|-----|--------|----------|------------|
| Antígeno 85B (Ag85B) | *Mycobacterium tuberculosis* | 1F0N | A | 2–285 | 1.80 Å |
| Spike RBD | SARS-CoV-2 | 6W41 | C | 333–527 | 3.08 Å |
| Cápside mayor L1 | VPH-16 | 1DZL | A | 20–474 | 3.50 Å |
| GAPDH (control negativo) | *Homo sapiens* | 4WNC | A | 3–335 | 1.99 Å |

> En 6W41 (RBD en complejo con el anticuerpo CR3022), la cadena C corresponde al RBD; las cadenas del anticuerpo (H y L) se eliminan en la preparación.

---

## Solución de problemas

| Problema | Solución |
|----------|----------|
| `ModuleNotFoundError: torch` | Verifica que el entorno conda correcto está activado (`conda activate bepipred3` o `discotope3`). |
| Descarga de ESM-2 muy lenta | Es normal la primera vez (~2.5 GB). Usa una conexión estable. |
| DiscoTope error de `pybind11` | Ejecuta `pip install pybind11` antes de instalar DiscoTope. |
| DiscoTope error de `biotite`/`torch_geometric` | Requiere Python ≥ 3.11. Recrea el entorno con `python=3.11`. |
| `parse_bepipred3` no encuentra el CSV | Verifica la ruta de salida; BepiPred puede generar solo el FASTA. |
| Aviso `CUDA not available` | Normal con `--cpu_only`; no es un error. |
| ElliPro no ejecuta | Verifica que Java esté instalado (`java -version`). |

---

## Cita

Si utilizas BEACON en tu investigación, por favor cita este trabajo:

```
Ramírez Sáenz, M. A. (2025). BEACON: Desarrollo y aplicación de un pipeline
de consenso para la predicción computacional de epítopes B en proteínas
antigénicas de patógenos de relevancia clínica global [Trabajo de Fin de
Máster]. Universidad Internacional de Valencia (VIU).
```

### Herramientas integradas

BEACON no sería posible sin las herramientas que integra. Cita también:

- **BepiPred-3.0:** Clifford, J. N., et al. (2022). *Protein Science, 31*(12), e4497.
- **DiscoTope-3.0:** Høie, M. H., et al. (2024). *Frontiers in Immunology, 15*, 1322712.
- **ElliPro:** Ponomarenko, J., et al. (2008). *BMC Bioinformatics, 9*, 514.
- **IEDB:** Vita, R., et al. (2025). *Nucleic Acids Research, 53*(D1), D436–D443.

---

## Licencia

Este proyecto se distribuye bajo la licencia MIT. Consulta el archivo [LICENSE](LICENSE) para más detalles.

---

*Desarrollado por Manuel Alain Ramírez Sáenz — Máster en Bioinformática, Universidad Internacional de Valencia (VIU), 2025.*
