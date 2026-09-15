#!/usr/bin/env python3
"""
BEACON - Preparacion de estructuras PDB para el pipeline de prediccion de epitopes B.

Descarga las estructuras desde RCSB, limpia heteroatomos y agua,
selecciona la cadena de interes, y extrae la secuencia FASTA.

Uso:
    python prepare_structures.py                  # Prepara las 3 proteinas del estudio
    python prepare_structures.py --pdb 1F0N --chain A  # Prepara una proteina especifica
"""

import urllib.request
import argparse
import os
import sys

AA_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M',
}


def download_pdb(pdb_id, outdir="data"):
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{pdb_id}.pdb")
    if not os.path.exists(outfile):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"  Descargando {pdb_id} desde RCSB...")
        urllib.request.urlretrieve(url, outfile)
    return outfile


def clean_pdb(input_pdb, output_pdb, chain='A', res_range=None):
    kept_lines = []
    residues = {}
    with open(input_pdb) as f:
        for line in f:
            if line[:4] == 'ATOM':
                chain_id = line[21]
                if chain_id != chain:
                    continue
                resname = line[17:20].strip()
                if resname not in AA_MAP:
                    continue
                res_seq = int(line[22:26].strip())
                if res_range and (res_seq < res_range[0] or res_seq > res_range[1]):
                    continue
                kept_lines.append(line)
                residues[res_seq] = resname
    kept_lines.append("END\n")
    with open(output_pdb, 'w') as f:
        f.writelines(kept_lines)
    sequence = "".join(AA_MAP.get(residues[r], 'X') for r in sorted(residues))
    return {
        'n_atoms': len(kept_lines) - 1,
        'n_residues': len(residues),
        'first_res': min(residues) if residues else None,
        'last_res': max(residues) if residues else None,
        'sequence': sequence,
    }


def write_fasta(filename, header, sequence, line_width=80):
    with open(filename, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i + line_width] + "\n")


def prepare_protein(pdb_id, chain, label, outdir="data", res_range=None):
    print(f"\nPreparando {label} (PDB: {pdb_id}, cadena {chain})")
    raw_pdb = download_pdb(pdb_id, outdir)
    suffix = f"_{chain}"
    if res_range:
        suffix += f"_{res_range[0]}-{res_range[1]}"
    clean_path = os.path.join(outdir, f"{pdb_id}{suffix}_clean.pdb")
    fasta_path = os.path.join(outdir, f"{pdb_id}{suffix}.fasta")
    info = clean_pdb(raw_pdb, clean_path, chain=chain, res_range=res_range)
    write_fasta(fasta_path, f"{pdb_id}_{chain}_{label}", info['sequence'])
    print(f"  PDB limpio:  {clean_path} ({info['n_atoms']:,} atomos)")
    print(f"  FASTA:       {fasta_path} ({info['n_residues']} residuos, {info['first_res']}-{info['last_res']})")
    return info


def main():
    parser = argparse.ArgumentParser(description="BEACON: Preparacion de estructuras PDB")
    parser.add_argument("--pdb", help="PDB ID especifico")
    parser.add_argument("--chain", default="A", help="Cadena a extraer (default: A)")
    parser.add_argument("--res-start", type=int, help="Residuo inicial")
    parser.add_argument("--res-end", type=int, help="Residuo final")
    parser.add_argument("--label", default="protein", help="Etiqueta FASTA")
    parser.add_argument("--outdir", default="data", help="Directorio de salida")
    args = parser.parse_args()

    if args.pdb:
        res_range = None
        if args.res_start and args.res_end:
            res_range = (args.res_start, args.res_end)
        prepare_protein(args.pdb, args.chain, args.label, args.outdir, res_range)
    else:
        proteins = [
            {"pdb": "1F0N", "chain": "A", "label": "Ag85B_M.tuberculosis", "range": None},
            {"pdb": "6W41", "chain": "C", "label": "Spike_RBD_SARS-CoV-2", "range": None},
            {"pdb": "1DZL", "chain": "A", "label": "L1_HPV-16", "range": None},
        ]
        for p in proteins:
            prepare_protein(p["pdb"], p["chain"], p["label"], args.outdir, p["range"])
        print("\nPreparacion completa.")


if __name__ == "__main__":
    main()
