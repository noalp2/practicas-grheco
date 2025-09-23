import pandas as pd
import os
import argparse
import random
from urllib.request import urlretrieve

def descargar_pares(tsv_file, output_dir, num_pares=3):
    # Crear carpeta si no existe
    os.makedirs(output_dir, exist_ok=True)

    # Leer archivo TSV
    df = pd.read_csv(tsv_file, sep='\t')

    # Filtrar tumoral y normal
    tumoral = df[df['sample_title'].str.endswith('_C')]
    normal = df[df['sample_title'].str.endswith('_N')]

    # Extraer pares válidos: buscar mismos prefijos antes de _C/_N
    pares = []
    for c_title in tumoral['sample_title']:
        base = c_title[:-2]  # quitar _C
        n_match = normal[normal['sample_title'].str.startswith(base)]
        if not n_match.empty:
            n_row = n_match.iloc[0]
            c_row = df[df['sample_title'] == c_title].iloc[0]
            pares.append((c_row, n_row))

    # Elegir aleatoriamente num_pares
    pares_seleccionados = random.sample(pares, k=min(num_pares, len(pares)))

    # Descargar archivos
    for c_row, n_row in pares_seleccionados:
        for row in [c_row, n_row]:
            urls = row['fastq_ftp'].split(';')
            sample_name = row['sample_title']
            for i, url in enumerate(urls, start=1):
                filename = os.path.join(output_dir, f"{sample_name}_{i}.fastq.gz")
                print(f"Descargando {filename}...")
                urlretrieve(f"ftp://{url}", filename)
    print("Descarga completa.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Descargar muestras pareadas de un TSV y renombrarlas.")
    parser.add_argument("-i", "--input", required=True, help="Archivo TSV con información de muestras")
    parser.add_argument("-o", "--output", required=True, help="Carpeta de destino para las descargas")
    args = parser.parse_args()

    descargar_pares(args.input, args.output)