#!/usr/bin/env python3
import pandas as pd
import os
import re

# Lista de muestras
samples = ["GBMEC5785", "GBMEC5441", "GBMEC5433", "GBMEC5377", "EC-1", "EC-2", "EC-3"]

# Ruta base relativa al script
base_dir = "../results/salmon"

# Diccionario para almacenar TPMs de cada muestra
tpm_dict = {}

for sample in samples:
    quant_path = os.path.join(base_dir, sample, "quant.sf")
    
    # Leer el archivo quant.sf
    df = pd.read_csv(quant_path, sep="\t")
    
    # Limpiar transcript ID: eliminar .1, .2... y quedarnos con transcriptID|GENE
    def clean_name(x):
        parts = x.split("|")
        transcript = re.sub(r"\.\d+$", "", parts[0])  # quita el .1, .2, etc.
        gene = parts[5] if len(parts) > 5 else "NA"
        return f"{transcript}|{gene}"
    
    df["Transcript_Gene"] = df["Name"].apply(clean_name)
    
    # Guardar solo Transcript_Gene y TPM
    tpm_dict[sample] = df.set_index("Transcript_Gene")["TPM"]

# Combinar todos los TPM en una tabla
tpm_matrix = pd.DataFrame(tpm_dict)

# Resetear índice para que Transcript_Gene quede como primera columna
tpm_matrix.reset_index(inplace=True)

# Guardar como TSV
output_path = "../results/TPM_matrix.tsv"
tpm_matrix.to_csv(output_path, sep="\t", index=False)

print(f"✅ Archivo generado: {output_path}")
