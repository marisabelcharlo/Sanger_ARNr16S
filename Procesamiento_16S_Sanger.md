# Generar secuencias limpias de calidad a partir de datos de secuenciación del gen ARNr 16S, obtenidas con la plataforma de Sanger

## Instalar Biophyton

`pip install biopython` o `conda install -c conda-forge biopython`

- Utilizé biophyton v.1.85

## Crear un archivo en el editor de texto nano

`nano trim_sanger_ab1.py`

Tipos de archivos en el directorio de trabajo para cada muestra:
- **.ab1** : Formato binario del equipo. Contiene: señal fluorescente (los “picos” A/C/G/T), llamado de bases (base calling), calidad por base (Phred), metadatos de la corrida.
- **.phd.1** : Texto plano generado por Phred. Lista base por base: nucleótido, calidad (Phred), posición del pico.
- **.txt** : Versión simplificada: solo la secuencia con o sin calidad.
- **.pdf** : Visualización del cromatograma

El script trabaja sobre todos los archivos **.ab1** presentes en el directorio, utilizando un enfoque de ventana móvil sobre los valores de calidad (Phred) de cada base. En este método, se define un tamaño de ventana (por ejemplo, 10 bases) y un umbral mínimo de calidad promedio (por ejemplo, Q20). El algoritmo recorre la secuencia desde el inicio buscando la primera región donde el promedio de calidad dentro de esa ventana alcanza o supera el umbral, estableciendo allí el punto de corte inicial; luego realiza el mismo procedimiento en sentido inverso desde el final para definir el corte terminal. De esta forma, se eliminan automáticamente los extremos de baja calidad característicos de la secuenciación Sanger, manteniendo únicamente el tramo central más confiable. Adicionalmente, se puede definir una longitud mínima de secuencia para descartar lecturas excesivamente cortas tras el recorte, asegurando consistencia y comparabilidad entre muestras (por ejemplo, 100 bp).

```bash 
cat > trim_sanger_ab1.py <<'EOF'
from Bio import SeqIO
from pathlib import Path

INPUT_DIR = Path(".")
OUTPUT_DIR = Path("trimmed")

MIN_QUALITY = 20
WINDOW_SIZE = 10
MIN_LENGTH = 100

OUTPUT_DIR.mkdir(exist_ok=True)

def trim_by_quality(record, min_q=20, window=10):
    qualities = record.letter_annotations["phred_quality"]
    seq_len = len(record.seq)

    start = 0
    end = seq_len

    for i in range(seq_len - window + 1):
        if sum(qualities[i:i+window]) / window >= min_q:
            start = i
            break

    for i in range(seq_len - window, -1, -1):
        if sum(qualities[i:i+window]) / window >= min_q:
            end = i + window
            break

    return record[start:end], start, end

summary = []

for ab1_file in sorted(INPUT_DIR.glob("*.ab1")):
    record = SeqIO.read(ab1_file, "abi")

    trimmed, start, end = trim_by_quality(record, MIN_QUALITY, WINDOW_SIZE)

    sample_name = ab1_file.stem
    trimmed.id = sample_name
    trimmed.name = sample_name
    trimmed.description = ""

    if len(trimmed.seq) >= MIN_LENGTH:
        SeqIO.write(trimmed, OUTPUT_DIR / f"{sample_name}_trimmed.fasta", "fasta")
        SeqIO.write(trimmed, OUTPUT_DIR / f"{sample_name}_trimmed.fastq", "fastq")
        status = "kept"
    else:
        status = "discarded_short"

    summary.append(
        f"{sample_name}\t{len(record.seq)}\t{len(trimmed.seq)}\t{start}\t{end}\t{status}"
    )

with open(OUTPUT_DIR / "trimming_summary.tsv", "w") as f:
    f.write("sample\toriginal_length\ttrimmed_length\tstart\tend\tstatus\n")
    f.write("\n".join(summary))
    f.write("\n")
EOF
```
Para cada archivo .ab1, el script genera dos archivos de salida si la secuencia cumple con la longitud mínima:

- un archivo **.fasta**, con la secuencia trimmeada
- un archivo **.fastq**, con la secuencia trimmeada y sus calidades asociadas.

## Crear un archivo multi-fasta

`cat trimmed/*_trimmed.fasta > todas_27F_trimmed.fasta`

## Verificar cuantas secuencias quedaron

`grep -c "^>" todas_27F_trimmed.fasta`

## Revisar encabezados

`grep "^>" todas_27F_trimmed.fasta`

## Revisar las longitudes

`awk '/^>/ {if (seq) print name, length(seq); name=$0; seq=""; next} {seq=seq$0} END {print name, length(seq)}' todas_27F_trimmed.fasta`

## Elegir una secuencia para copiar y pega en BLASTn

`grep -A 1 ">AMR-3_27F" todas_27F_trimmed.fasta`


