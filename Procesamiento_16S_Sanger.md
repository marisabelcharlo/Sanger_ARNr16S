# Procesamiento de datos Sanger del gen ARNr 16S: trimming, ensamblado y obtención de secuencias consenso

## Instalar Biopython

`pip install biopython` o `conda install -c conda-forge biopython`

- Utilicé Biopython v.1.85

## Crear el script en el editor de texto nano

```bash
nano trim_and_merge_16S.py
```

Tipos de archivos en el directorio de trabajo para cada muestra:
- **.ab1** : Formato binario del equipo. Contiene: señal fluorescente (los “picos” A/C/G/T), llamado de bases (base calling), calidad por base (Phred), metadatos de la corrida.
- **.phd.1** : Texto plano generado por Phred. Lista base por base: nucleótido, calidad (Phred), posición del pico.
- **.txt** : Versión simplificada: solo la secuencia con o sin calidad.
- **.pdf** : Visualización del cromatograma

El script trabaja sobre todos los archivos .ab1 presentes en el directorio y automatiza el procesamiento de las lecturas Sanger correspondientes a los primers 27F y 1492R. En primer lugar, realiza un recorte fijo de los extremos de cada lectura, eliminando 20 nucleótidos del inicio y 100 nucleótidos del final, independientemente de la calidad. Luego aplica un trimming basado en calidad utilizando una ventana móvil sobre los valores Phred de cada base. Para ello, se define un tamaño de ventana de 10 bases y un umbral mínimo de calidad promedio de Q20. El algoritmo recorre la secuencia desde el inicio buscando la primera región donde el promedio de calidad alcanza o supera ese umbral, estableciendo allí el punto de corte inicial; luego realiza el mismo procedimiento desde el final para definir el corte terminal. De esta forma, se eliminan automáticamente los extremos de baja calidad. Posteriormente, las lecturas obtenidas con el primer 1492R son convertidas a su reverso complementario, de modo que queden en la misma orientación que las lecturas 27F. Finalmente, el script identifica el solapamiento entre ambas lecturas mediante alineamiento local y genera una secuencia consenso por muestra, seleccionando en las posiciones discordantes la base con mayor calidad. Como salida, se generan archivos FASTA y FASTQ para las lecturas trimmeadas, archivos FASTA y FASTQ para las secuencias consenso, y tablas resumen del trimming y del ensamblado.

## Parámetros utilizados

- `MIN_QUALITY = 20`: calidad Phred mínima promedio aceptada.
- `WINDOW_SIZE = 10`: número de bases usado para calcular el promedio móvil de calidad.
- `MIN_LENGTH = 100`: longitud mínima para conservar una lectura luego del trimming.
- `FIXED_TRIM_START = 20`: número de bases eliminadas al inicio, independientemente de la calidad.
- `FIXED_TRIM_END = 100`: número de bases eliminadas al final, independientemente de la calidad.
- `MIN_OVERLAP = 50`: longitud mínima del solapamiento entre 27F y 1492R.
- `MIN_IDENTITY = 0.75`: identidad mínima requerida en el solapamiento para aceptar el ensamblado.

```bash 
cat > trim_and_merge_16S.py <<'EOF'
from Bio import SeqIO, pairwise2
from Bio.SeqRecord import SeqRecord
from pathlib import Path

INPUT_DIR = Path(".")
OUT_DIR = Path("trimmed_merged")
TRIM_DIR = OUT_DIR / "trimmed_reads"
CONS_DIR = OUT_DIR / "consensus"

MIN_QUALITY = 20
WINDOW_SIZE = 10
MIN_LENGTH = 100

FIXED_TRIM_START = 20
FIXED_TRIM_END = 100

MIN_OVERLAP = 50
MIN_IDENTITY = 0.75

OUT_DIR.mkdir(exist_ok=True)
TRIM_DIR.mkdir(exist_ok=True)
CONS_DIR.mkdir(exist_ok=True)

def fixed_trim(record, start_trim=20, end_trim=100):
    if len(record.seq) <= start_trim + end_trim:
        return record[0:0]
    return record[start_trim:len(record.seq)-end_trim]

def quality_trim(record, min_q=20, window=10):
    qualities = record.letter_annotations["phred_quality"]
    seq_len = len(record.seq)

    if seq_len < window:
        return record[0:0], 0, 0

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

def sample_id_from_name(name):
    return (
        name.replace("_27F", "")
            .replace("_1492R", "")
            .replace("_trimmed", "")
    )

def find_best_overlap_aln(fwd, rev, min_overlap=50, min_identity=0.75):
    fseq = str(fwd.seq).upper()
    rseq = str(rev.seq).upper()

    alignments = pairwise2.align.localms(
        fseq, rseq,
        2,      # match
        -1,     # mismatch
        -2,     # gap open
        -0.5,   # gap extend
        one_alignment_only=True
    )

    if not alignments:
        return None

    aln_f, aln_r, score, start, end = alignments[0]

    matches = 0
    aligned_length = 0

    for a, b in zip(aln_f, aln_r):
        if a != "-" and b != "-":
            aligned_length += 1
            if a == b and a != "N" and b != "N":
                matches += 1

    if aligned_length < min_overlap:
        return None

    identity = matches / aligned_length

    if identity < min_identity:
        return None

    return aligned_length, identity, aln_f, aln_r

def merge_from_alignment(fwd, rev, aln_f, aln_r, sample):
    fq = fwd.letter_annotations["phred_quality"]
    rq = rev.letter_annotations["phred_quality"]

    consensus_seq = []
    consensus_qual = []

    i_f = 0
    i_r = 0

    for a, b in zip(aln_f, aln_r):
        if a != "-" and b != "-":
            qf = fq[i_f]
            qr = rq[i_r]

            if a == b:
                consensus_seq.append(a)
                consensus_qual.append(max(qf, qr))
            else:
                if qf >= qr:
                    consensus_seq.append(a)
                    consensus_qual.append(qf)
                else:
                    consensus_seq.append(b)
                    consensus_qual.append(qr)

            i_f += 1
            i_r += 1

        elif a != "-" and b == "-":
            consensus_seq.append(a)
            consensus_qual.append(fq[i_f])
            i_f += 1

        elif a == "-" and b != "-":
            consensus_seq.append(b)
            consensus_qual.append(rq[i_r])
            i_r += 1

    consensus = SeqRecord(
        fwd.seq.__class__("".join(consensus_seq)),
        id=sample,
        name=sample,
        description="consensus_SW_alignment"
    )

    consensus.letter_annotations["phred_quality"] = consensus_qual
    return consensus

reads_27F = {}
reads_1492R = {}

summary_trim = []
summary_merge = []

for ab1_file in sorted(INPUT_DIR.glob("*.ab1")):
    record = SeqIO.read(ab1_file, "abi")
    record.id = ab1_file.stem
    record.name = ab1_file.stem
    record.description = ""

    fixed = fixed_trim(record, FIXED_TRIM_START, FIXED_TRIM_END)
    trimmed, q_start, q_end = quality_trim(fixed, MIN_QUALITY, WINDOW_SIZE)

    sample = sample_id_from_name(ab1_file.stem)

    if "1492R" in ab1_file.stem:
        trimmed = trimmed.reverse_complement(
            id=ab1_file.stem,
            name=ab1_file.stem,
            description="trimmed_reverse_complement"
        )
        reads_1492R[sample] = trimmed
        orientation = "1492R_reverse_complement"

    elif "27F" in ab1_file.stem:
        reads_27F[sample] = trimmed
        orientation = "27F"

    else:
        orientation = "unknown"

    status = "kept" if len(trimmed.seq) >= MIN_LENGTH else "discarded_short"

    if status == "kept":
        SeqIO.write(trimmed, TRIM_DIR / f"{ab1_file.stem}_trimmed.fasta", "fasta")
        SeqIO.write(trimmed, TRIM_DIR / f"{ab1_file.stem}_trimmed.fastq", "fastq")

    summary_trim.append(
        f"{ab1_file.stem}\t{orientation}\t{len(record.seq)}\t{len(fixed.seq)}\t{len(trimmed.seq)}\t{status}"
    )

all_samples = sorted(set(reads_27F.keys()) | set(reads_1492R.keys()))

for sample in all_samples:
    if sample not in reads_27F:
        summary_merge.append(f"{sample}\tmissing_27F\tNA\tNA\tNA\tNA\tNA")
        continue

    if sample not in reads_1492R:
        summary_merge.append(f"{sample}\tmissing_1492R\tNA\tNA\tNA\tNA\tNA")
        continue

    fwd = reads_27F[sample]
    rev = reads_1492R[sample]

    best = find_best_overlap_aln(fwd, rev, MIN_OVERLAP, MIN_IDENTITY)

    if best is None:
        summary_merge.append(
            f"{sample}\tno_overlap\t{len(fwd.seq)}\t{len(rev.seq)}\tNA\tNA\tNA"
        )
        continue

    overlap_len, identity, aln_f, aln_r = best
    consensus = merge_from_alignment(fwd, rev, aln_f, aln_r, sample)

    SeqIO.write(consensus, CONS_DIR / f"{sample}_consensus.fasta", "fasta")
    SeqIO.write(consensus, CONS_DIR / f"{sample}_consensus.fastq", "fastq")

    summary_merge.append(
        f"{sample}\tmerged\t{len(fwd.seq)}\t{len(rev.seq)}\t{overlap_len}\t{identity:.4f}\t{len(consensus.seq)}"
    )

with open(OUT_DIR / "trimming_summary.tsv", "w") as f:
    f.write("read\torientation\toriginal_length\tafter_fixed_trim_length\tfinal_trimmed_length\tstatus\n")
    f.write("\n".join(summary_trim))
    f.write("\n")

with open(OUT_DIR / "merge_summary.tsv", "w") as f:
    f.write("sample\tstatus\tlength_27F\tlength_1492R_RC\toverlap\tidentity\tconsensus_length\n")
    f.write("\n".join(summary_merge))
    f.write("\n")
EOF
```

```bash
python trim_and_merge_16S.py
```

```bash
grep -c "merged" trimmed_merged/merge_summary.tsv
grep -c "no_overlap" trimmed_merged/merge_summary.tsv
grep -c "missing_1492R" trimmed_merged/merge_summary.tsv
```

```bash
column -t trimmed_merged/merge_summary.tsv | less -S
```

```bash
awk '/^>/ {if (seq) print name, length(seq); name=$0; seq=""; next} {seq=seq$0} END {print name, length(seq)}' trimmed_merged/consensus/*_consensus.fasta
```

## Salidas generadas

El script genera la carpeta `trimmed_merged/`.

En `trimmed_reads/` se guardan las lecturas individuales luego del recorte. Para cada muestra se genera un archivo `.fasta` y un archivo `.fastq`.

En `consensus/` se guardan las secuencias consenso obtenidas tras empalmar las lecturas 27F y 1492R.

El archivo `trimming_summary.tsv` resume el trimming de cada lectura.

El archivo `merge_summary.tsv` resume el empalme por muestra, incluyendo longitud de cada lectura, tamaño del solapamiento, identidad del solapamiento y longitud final del consenso.

## Criterios de control de calidad

- Consensos de aproximadamente 1300–1500 bp: ensamblado compatible con 16S casi completo.
- Consensos <1000 bp: ensamblado incompleto o lectura parcial; revisar.
- Consensos >1600 bp: posible error de ensamblado; revisar.
- Identidad baja en el solapamiento: posible baja calidad, mezcla o empalme forzado.

