# Fragmentomics_GenomBiol
Pipeline to replicate Katsman et. al fragmentomics results

Requires:
guppy
minimap2
bwa
samtools
picard
bedtools
R 
perl

Rpackages: "vroom", "data.table", "ggplot2", "ComplexHeatmap", "circlize", "pals", "ggplot2", "ggbeeswarm", "ggpubr"

## Reads Preprocessing


### Illumina

Trim Ns at 3’ end
TO BE COMPLETED
/home/guest/FASTQ_NO_Ns3P/47503_ID1514_2-19-744_S1_L001_R1_001.trimN3p.fastq.gz>

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
minimap2 -ax map-ont --MD -L reference.mmi sample.fastq.gz | samtools view -h -q 20 -F 0x4 -F 0x100 -F 0x800 | ~/Fragmentomics/Scripts/General/samCigarToTlen.pl | awk '( $9 < 700 || $1 ~ /^@/ )' | samtools view -bS -  -o sample.filtered.bam
```


Then store your filtered bams in a folder of your choice for the subsequent steps (in this tutorial ~/Fragmentomics/Data)
Perform all the subsequent commands from that folder.
The Genome_Biol_Data folder contains .stats files from Katsman et al. to replicate and test the pipeline (original .bam files are not provided here for ethical reasons, but you don't need them as you can perform the whole analysis directly from the files provided here)
The Genome_Biol_Data_analyzed contains all the analysis from Katsman et al.

## Fragmentomics

### 5' 4-mer motifs

Enter the folder in which the filtered bams (from now on "sample.bam") are stored (~/Fragmentomics/Data/) and create the folders for the subsequent analysis

```
mkdir -p STATS/MOTIF/MOTIF_COUNTS
```
Create .stats files (custom format)

```
samtools view  sample.bam | perl ~/Fragmentomics/Scripts/General/stats_maker_0basedstart.pl > STATS/sample.stats
```
.stats files are 15 column, tab separated files which are used for both 4-mer and read-length analysis. 
In order, the columns are:
1) Contig name
2) Alignment start position (0-based)
3) Alignment end position (1-based)
4) Read-id
5) 0 (for formatting purposes)
6) strand
7) sam flag
8) mapping quality
9) Left Hard clip length (from CIGAR)
10) Left Soft clip length (from CIGAR)
11) Aligment length (from CIGAR)
12) Right Soft clip length (from CIGAR)
13) Right Hard clip length (from CIGAR)
14) Mate orientation (F or R, only for illumina)
15) bam TLEN field

```
NC_000014.9	42181158	42181344	00000a94-0391-41f5-9a31-8fa0ba907910	0	+	0	60	0	0	186	1	0	NA	186
NC_000002.12	209229299	209229465	00000edf-2af4-4f11-bfb6-5c9800f27702	0	-	16	60	0	1	166	0	0	NA	166
NC_000007.14	41236402	41236569	00000f90-3b76-4e4a-925b-030dffc521dd	0	+	0	60	0	164	167	54	0	NA	167
...

...
NC_000007.14	30069730	30069868	0000a877-57f5-4eba-859c-81ef1197f442	0	+	0	60	0	0	138	0	0	NA	138

```
Use bedtools to extract template sequences using .stats files as bed. Then use awk to extract the first 4 nucleotides and print read names

```
bedtools getfasta  -s  -name -tab -fi reference.fa -bed  STATS/sample.stats | awk '{{ print $1, substr($2,1,4)}}' | awk -v FS='::' -v OFS=' ' '{{print $1,$2,$3}}'   >   STATS/MOTIF/sample.motif
```
.motif files are 3 column, space separated files which report:
1) Read-id
2) contig:start-end(strand)
3) 5' 4-mer

```
00000a1a-f400-4203-9aaf-5008d975bf73 NC_000003.12:49124212-49124348(+) CCCC 
00000b41-de8c-4fe4-89fd-952d601ccef5 NC_000014.9:82929351-82929599(-) CTGG 
00000f0b-ef56-401b-88ae-54085e19c8b5 NC_000020.11:4622611-4622769(-) cttt 
...

...
0000a096-3c26-46ad-bb20-4ad452db2ea0 NC_000001.11:239968433-239968602(-) CCTT 
```
Obtain the counts of each 4-mer motif using the custom script provided and the .motif files previously produced.
You have to provide a list of the chromosome names you are including in the analysis in a chr_list.txt file (one chromosome per row). The file used for Katsman et. al is provided in the ~/Fragmentomics/Utility folder.

```
Rscript ~/Fragmentomics/Scripts/Motifs/count_motif.R ~/Fragmentomics/Utility/chr_list.txt  STATS STATS/MOTIF STATS/MOTIF/MOTIF_COUNTS/  sample  
```
This will produce .motif.R files, which are R objects (tables) with the raw count of each 4-mer motif.

Use ~/Fragmentomics/Scripts/Motifs/motif_heatmap.R file to make the heatmap from Figure 2
Use ~/Fragmentomics/Scripts/Motifs/plot_motif_CCCA.R file to make the CCCA jitterplot from Figure 2

For both the scripts you have to provide a samples_info_heatmap.tsv file (the one used for Katsman et al. is included in the ~/Fragmentomics/Utility folder) with the following columns:

1) MOTIF_COUNTS folder full path
2) sample basename (file name without the extension)
3) sample printed name
4) Disease annotation: Cancer, Healthy (Can be specified if HU or ISPRO samples)
5) Sample sex (not mandatory)
6) Library type (Illumina, Nanopore, Nanopore_HU)
7) Sample type (Healthy, Lung)
8) Tumor Fraction 

```
~/Fragmentomics/Data/STATS/MOTIF/MOTIF_COUNTS/	BC01.HAC	BC01	Cancer	M	Nanopore	Lung	0.252
~/Fragmentomics/Data/STATS/MOTIF/MOTIF_COUNTS/	BC04.HAC	BC04	ISPRO	M	Nanopore	Healthy	0
~/Fragmentomics/Data/STATS/MOTIF/MOTIF_COUNTS/	BC08_ILL	BC08	Cancer	M	Illumina	Lung	0.105
```

### Read length
create folder for read length counts files
```
mkdir -p STATS/READLENGTH_COUNTS
```
#### Illumina

Obtain the counts of each readlength bin (1bp)
You have to provide a list of the chromosome names you are including in the analysis in a chr_list.txt file (one chromosome per row). The file used for Katsman et. al is provided in the ~/Fragmentomics/Utility folder.

```
/usr/bin/Rscript ~/Fragmentomics/Scripts/Read_length/count_readlength_illumina.R ~/Fragmentomics/Utility/chr_list.txt STATS  STATS/READLENGTH_COUNTS sample
```
#### Nanopore 

Obtain the counts of each readlength bin (1bp)
You have to provide a list of the chromosome names you are including in the analysis in a chr_list.txt file (one chromosome per row). The file used for Katsman et. al is provided in the ~/Fragmentomics/Utility folder.
You can provide a list of read ids which will be used for the analysis. For the paper, for multiplex runs, we used the ids of reads which had barcodes at both ends.
```
/usr/bin/Rscript ~/Fragmentomics/Scripts/Read_length/count_readlength_nanopore_bothbarcodes.R ~/Fragmentomics/Utility/chr_list.txt  STATS STATS/READLENGTH_COUNTS  sample ~/Fragmentomics/Data/BOTH_BARCODES/sample.ids 
```
Use "NO" if you don't want to provide an .ids file (for example for singleplex runs)
```
/usr/bin/Rscript ~/Fragmentomics/Scripts/Read_length/count_readlength_nanopore_bothbarcodes.R ~/Fragmentomics/Utility/chr_list.txt  STATS STATS/READLENGTH_COUNTS  sample  NO 
```
you can create a list of read ids from a fastq (for example obtained via demultiplexing with guppy using the "--both-barcodes" flag) in this way
```
gzip -dc sample.fastq.gz | grep -o  "@........-....-....-....-............"  | sed 's/@//' > ~/Fragmentomics/Data/BOTH_BARCODES/sample.ids
```

Use ~/Fragmentomics/Scripts/Read_length/readlengh_plots.R file to make the read length jitterplots and density plots from figure 2.

You have to provide a samples_info_readlength.tsv file (the one used for Katsman et al. is included in the ~/Fragmentomics/Utility folder) with the following columns:

1) READLENGTH_COUNTS folder full path
2) sample basename (file name without the extension)
3) sample printed name
4) Disease annotation: Cancer, Healthy (Can be specified if HU or ISPRO samples)
5) Sample sex (not mandatory)
6) Library type (Illumina, Nanopore, Nanopore_HU)
7) Sample type (Healthy, Lung)
8) Tumor Fraction 

```
~/Fragmentomics/STATS/READLENGTH_COUNTS/	BC01.HAC	BC01	Cancer	M	Nanopore	Lung	0.252
~/Fragmentomics/STATS/READLENGTH_COUNTS/	BC04.HAC	BC04	ISPRO	M	Nanopore	Healthy	0
~/Fragmentomics/STATS/READLENGTH_COUNTS/	BC08_ILL	BC08	Cancer	M	Illumina	Lung	0.105
```

