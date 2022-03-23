library(vroom)

chr_list <- commandArgs(trailingOnly=TRUE)[1]
stats_path <- commandArgs(trailingOnly=TRUE)[2]
out_path <- commandArgs(trailingOnly=TRUE)[3]
sample_name <- commandArgs(trailingOnly=TRUE)[4]

stats_path <- file.path(stats_path,paste(sample_name,".stats", sep=""))

#### loading chromosome list
chr <- read.table(chr_list, stringsAsFactors = F)

#### loading stat files
stats <- as.data.frame(vroom(stats_path, show_col_types = FALSE, col_names = F))
colnames(stats) <- paste("V", seq(1,length(colnames(stats))), sep="")

#### keeping reads mapped on selected chromosomes, with MAPQ > 20 and without 5' hard/soft clipping
stats <- subset(stats, stats$V8 > 20 & stats$V1 %in% chr[,1] & ((stats$V6 == "+" & stats$V9 == 0 & stats$V10 == 0) | (stats$V6 == "-" & stats$V12 == 0 & stats$V13 == 0)))

#### selecting only one read among the ones that are reported twice (i.e. the DNA fragments for which both the mates have passed the filtering criterias)
stats <- subset(stats, duplicated(stats$V4))

if (!(file.exists(out_path))){
  dir.create(out_path)
}

#### creating count table for each readlength value (calculated from TLEN)
length_counts <- table(abs(stats$V15))
saveRDS(length_counts, file.path(out_path,paste(sample_name,".readlength_counts.R", sep="")))
