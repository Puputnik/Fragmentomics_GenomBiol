# Fragmentomics_GenomBiol
Pipeline to replicate Katsman et. al fragmentomics results

## Reads Preprocessing


### Illumina

trimma le N al 3’
i file li trovi /home/guest/FASTQ_NO_Ns3P/47503_ID1514_2-19-744_S1_L001_R1_001.trimN3p.fastq.gz>

Align with bwa mem and filter out unmapped reads, supplementary and secondary alignments. Sort by coordinates

```
bwa mem reference.fa sample_R1.fastq sample_R2.fastq  | samtools view -F 0x4 -F 0x100 -F 0x800 -b | samtools sort -O BAM -o sample.srt.bam 
```
Mark duplicates using picard

```
java -jar /bin/picard.jar MarkDuplicates I=sample.srt.bam  O=sample.marked_duplicates.srt.bam M=sample.metrics 
```
Filter out duplicates, reads not in proper pair and alignments with mapping quality < 20. Filter out aligments with TLEN >= 700bp

```
samtools view -F 0x400 -q 20 -f 0x2 -h sample.marked_duplicates.srt.bam | awk '( $9 < 700 || $1 ~ /^@/ )' | samtools view -bS -  -o sample.filtered.bam
```

### Nanopore

Basecall HAC

guppy demultiplex trim barcodes
both barcodes version
		keep the ids.

Align with minimap2 and filter out unmapped reads, supplementary and secondary alignments, and alignments with mapping quality < 20. Annotate TLEN field with read length obtained from CIGAR string (custom script provided) and filter out alignments with TLEN >= 700bp 

```
minimap2 -ax map-ont --MD -L reference.mmi sample.fastq.gz | samtools view -h -q 20 -F 0x4 -F 0x100 -F 0x800 | ~/Fragmentomics/samCigarToTlen.pl | awk '( $9 < 700 || $1 ~ /^@/ )' | samtools view -bS -  -o sample.filtered.bam
```




Fragmentomics

motifs
mkdir -p STATS/MOTIF/MOTIF_COUNTS


        samtools view  sample.bam | perl ~/Fragmentomics/stats_maker_0basedstart.pl > STATS/sample.stats

bedtools getfasta  -s  -name -tab -fi reference.fa -bed  STATS/sample.stats | awk '{{ print $1, substr($2,1,4)}}' | awk -v FS='::' -v OFS=' ' '{{print $1,$2,$3}}'   >   STATS/MOTIF/sample.motif

Rscript ~/Fragmentomics/count_motif.R ~/Fragmentomics/chr_list.txt  STATS STATS/MOTIF STATS/MOTIF/MOTIF_COUNTS/  sample  
