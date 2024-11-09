# Linux obuka

## Vežba SARS-COV-2

Ispod se nalazi kod koji ćemo koristiti na obuci Uvod u Linux u toku sledeće nedelje (utorak, sreda, petak).
Između znakova <> se nalaze delovi koda koji treba da se zamene odgovarajućim fajlom.

### 1. Fastqc

```bash
fastqc -o reports -t 2 fastq/V300114179_L04_80_80_?.fq.gz
```
### 2. Fastp
```bash
zgrep -A2 -B1 --color AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA fastq/V300114179_L04_80_80_1.fq.gz
```
```bash
fastp --average_qual 28 -w 2 -i fastq/V300114179_L04_80_80_1.fq.gz -I fastq/V300114179_L04_80_80_2.fq.gz \
-o clean/uzorak1_1.fq.gz -O clean/uzorak1_2.fq.gz \
-j reports/uzorak1.json -h reports/uzorak1.html \
--adapter_sequence=AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
--adapter_sequence_r2=AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG
```
### 3. Indeksiranje referentnog genoma

```bash
bwa index reference/sequence.fasta
```

### 4. Poravnanje sa referentnim genomom
```bash
bwa mem -t 4 -R
'@RG\tID:uzorak1\tPL:DNBSEQ\tPU:uzorak1\tLB:mutPCR\tSM
:uzorak1' reference/sequence.fasta
clean/uzorak1_?.fq.gz | samtools view -b |samtools
sort > bam/uzorak1_sortiran.bam
```

### 5. Statistika BAM fajla
```bash
samtools flagstat bam/uzorak1_sortiran.bam
```

### 6. Filtriranje nemapiranih očitavanja
```bash
samtools view -f 2 -b bam/uzorak1_sortiran.bam >
bam/uzorak1_filter.bam
```

### 7. Qualimap
```bash
qualimap bamqc -bam bam/uzorak1_filter.bam -outdir
reports/uzorak1
```

### 8. Indeksiranje BAM fajla
```bash
samtools index bam/uzorak1_filter.bam
```

### 9. Provera kvaliteta baza
```bash
fastqc -o reports -t 2 bam/uzorak1_filter.bam
```
