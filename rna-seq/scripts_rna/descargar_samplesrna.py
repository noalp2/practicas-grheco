import pandas as pd
import os
import argparse
import random
from urllib.request import urlretrieve

def descargar_muestras(tsv_file, output_dir, n=3):
    # Crear carpeta de salida si no existe
    os.makedirs(output_dir, exist_ok=True)

    # Leer archivo TSV
    df = pd.read_csv(tsv_file, sep='\t', engine="python", comment="#")

    # Crear listas de cáncer y normales
    cancer = df[df['sample_title'].str.contains("GBM", case=False, na=False)].copy()
    cancer["grupo"] = "_C"

    normal = df[~df['sample_title'].str.contains("GBM", case=False, na=False)].copy()
    normal["grupo"] = "_N"

    if cancer.empty or normal.empty:
        raise ValueError("No se encontraron suficientes muestras en el archivo TSV")

    # Selección aleatoria
    cancer_sel = cancer.sample(min(n, len(cancer)), random_state=None)
    normal_sel = normal.sample(min(n, len(normal)), random_state=None)

    seleccionadas = pd.concat([cancer_sel, normal_sel])

    # Descargar y renombrar
    for _, row in seleccionadas.iterrows():
        urls = row['fastq_ftp'].split(';')
        sample_name = row['sample_title'].replace(" ", "_") + row['grupo']
        for i, url in enumerate(urls, start=1):
            filename = os.path.join(output_dir, f"{sample_name}_{i}.fastq.gz")
            if os.path.exists(filename):
                print(f"⚠️  Ya existe {filename}, se omite la descarga.")
            else:
                print(f"⬇️ Descargando {filename} ...")
                urlretrieve(f"ftp://{url}", filename)

    # Guardar log de las muestras seleccionadas
    log_file = os.path.join(output_dir, "muestras_descargadas.csv")
    seleccionadas.to_csv(log_file, sep="\t", index=False)
    print(f"✅ Descarga completa. Lista de muestras guardada en: {log_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Descargar muestras de RNAseq y renombrarlas según sample_title")
    parser.add_argument("-i", "--input", required=True, help="Archivo TSV con información de muestras")
    parser.add_argument("-o", "--output", required=True, help="Carpeta de destino para las descargas")
    parser.add_argument("-n", "--num", type=int, default=3, help="Número de muestras de cada grupo (default=3)")
    args = parser.parse_args()

    descargar_muestras(args.input, args.output, args.num)
