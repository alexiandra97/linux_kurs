# Linux obuka

## Vežba SARS-COV-2

Ispod se nalazi <i> bash </i> kod koji ćemo koristiti na obuci Uvod u Linux u toku sledeće nedelje (utorak, sreda, petak).

Između znakova <> se nalaze delovi koda koji treba da se zamene odgovarajućim fajlom.

### 1. Fastqc

```bash
fastqc -o reports -t 2 fastq/<fastq fajl>
```
### 2. Fastp
```bash
zgrep -A2 -B1 --color AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA fastq/<fastq fajl>
```
```bash
fastp --average_qual 28 -w 2 -i fastq/<fastq fajl> -I fastq/<fastq fajl> \
-o clean/uzorak1_1.fq.gz -O clean/uzorak1_2.fq.gz \
-j reports/uzorak1.json -h reports/uzorak1.html \
--adapter_sequence=AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
--adapter_sequence_r2=AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG
```
### 3. Indeksiranje referentnog genoma

```bash
bwa index reference/<fasta fajl>
```

### 4. Poravnanje sa referentnim genomom
```bash
bwa mem -t 4 -R \
'@RG\tID:uzorak1\tPL:DNBSEQ\tPU:uzorak1\tLB:mutPCR\tSM:uzorak1' reference/<fasta fajl> \
clean/uzorak1_?.fq.gz | samtools view -b |samtools sort \
> bam/uzorak1_sortiran.bam
```

### 5. Statistika BAM fajla
```bash
samtools flagstat bam/uzorak1_sortiran.bam
```

### 6. Filtriranje nemapiranih očitavanja
```bash
samtools view -f 2 -b bam/uzorak1_sortiran.bam > bam/uzorak1_filter.bam
```

### 7. Qualimap
```bash
qualimap bamqc -bam bam/uzorak1_filter.bam -outdir reports/uzorak1
```

### 8. Indeksiranje BAM fajla
```bash
samtools index bam/uzorak1_filter.bam
```

### 9. Provera kvaliteta baza
```bash
fastqc -o reports -t 2 bam/uzorak1_filter.bam
```

### 10. Variant calling
```bash
freebayes -p 1 -q 25 -m 60 --min-coverage 30 -f reference/<fasta fajl> \
bam/uzorak1_filter.bam > vcf/uzorak1_raw.vcf
```
### 11. Zipovanje i indeksiranje vcf fajla
```bash
bgzip vcf/uzorak1_raw.vcf
bcftools index vcf/uzorak1_raw.vcf.gz
```
### 12. Filtriranje vcf fajla
```bash
bcftools filter -i 'QUAL>20 && INFO/DP>10' vcf/uzorak1_raw.vcf.gz > \
vcf/uzorak1_filter.vcf
```

### 13. Brojanje varijanti
```bash
zgrep -c -v "^#" vcf/uzorak1_filter.vcf.gz
```

### 14. Spike protein (filtriranje po poziciji)
```bash
bcftools view -r NC_045512.2:21563-25384 vcf/uzorak1_filter.vcf.gz
```
### 15. Vcf bez header-a
```bash
bcftools view -H vcf/uzorak1_filter.vcf.gz |wc -l
```
### 16. Uključuje varijante sa kvalitetom >=20
```bash
bcftools view -H -i 'QUAL>=20' vcf/uzorak1_filter.vcf.gz | wc -l
```

### 17. Dobijanje konsenzusne sekvence
```bash
bcftools consensus -f reference/<fasta fajl> \
vcf/uzorak1_filter.vcf.gz | sed "s/^>.*$/>uzorak1/g" > fasta/uzorak1.fasta
```
### 18. Filtriranje vcf fajla
```bash
bcftools view -r NC_045512.2:21563-25384 vcf/uzorak1_filter.vcf.gz
```
```bash
bcftools view -H vcf/uzorak1_filter.vcf.gz |wc -l
```
```bash
bcftools view -H -i 'QUAL>=20' vcf/uzorak1_filter.vcf.gz | wc -l
```

### 19. Konverzija vcf u tsv format
```bash
gzip -d -k vcf/uzorak1_filter.vcf.gz
vcf2tsv vcf/uzorak1_filter.vcf vcf/uzorak1_filter.tsv
```
### 20. Preuzimanje gff3
```bash
wget -O reference/sequence.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?save=file&db=nuccore&report=gff3&id=1798174254"
```
### 21. Pravljenje bed fajla
```bash
awk -F "\t" 'BEGIN { OFS="\t" } $3=="gene" { split($9,f,";"); gene=f[3]; gsub(/Name=/,"",gene); print $1,$4,$5,gene }' reference/sequence.gff3 > reference/sequence_genes.bed
```
### 22. Pravljenje header-a za fajl
```bash
echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">' > gene.header
```

### 22. Anotacija
```bash
$ bcftools annotate -a reference/sequence_genes.bed -c CHROM,FROM,TO,GENE -h gene.header vcf/uzorak1_filter.vcf.gz > vcf/uzorak1_anotiran.vcf
```
