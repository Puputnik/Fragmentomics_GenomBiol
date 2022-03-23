library(vroom)

#chr_list <-   "/tank/USB/LIQUID_BIOPSY_MAIN_REVISION/SCRIPTS/chr_list.txt"
#stats_path <- "/tank/USB3/FAST5_AMIR/HAC/700_CUT/STATS/"
#out_path <-   "/tank/USB3/FAST5_AMIR/HAC/700_CUT/STATS/READLENGTH_COUNTS/"
#sample_name <- "BC02.HAC"
###id_path <-   "/tank/USB4/LIQUID_BIOPSY_FINAL/BOTH_BARCODES_IDs"
#id_path <-   "/tank/USB4/LIQUID_BIOPSY_FINAL/BOTH_BARCODES_IDs/BC02.id.txt"

chr_list <-      commandArgs(trailingOnly=TRUE)[1]
stats_path <-    commandArgs(trailingOnly=TRUE)[2]
out_path <-      commandArgs(trailingOnly=TRUE)[3]
sample_name <-   commandArgs(trailingOnly=TRUE)[4]
id_path <-   commandArgs(trailingOnly=TRUE)[5]


stats_path <- file.path(stats_path,paste(sample_name,".stats", sep=""))

#### loading chromosome list
chr <- read.table(chr_list, stringsAsFactors = F)

#### loading stat files
stats <- as.data.frame(vroom(stats_path, show_col_types = FALSE, col_names = F))
colnames(stats) <- paste("V", seq(1,length(colnames(stats))), sep="")

#### keeping reads mapped on selected chromosomes, with MAPQ > 20, readlength < 700 and without 5'/3' hard/soft clipping
stats <- subset(stats, stats$V8 > 20 & stats$V1 %in% chr[,1] & stats$V11 < 700 & stats$V9 == 0 & stats$V10 == 0  & stats$V12 == 0 & stats$V13 == 0) 


#both <- as.data.frame(vroom(file.path(id_path, paste(sample_name, ".id.txt", sep="")), delim="\t", show_col_types = FALSE, col_names = F))

if (id_path != "NO"){
both <- as.data.frame(vroom(id_path, delim="\t", show_col_types = FALSE, col_names = F))
colnames(both) <- "V4"
stats <- merge(stats, both, by="V4", all=F)
}

if (!(file.exists(out_path))){
  dir.create(out_path)
}

#### creating count table for each readlength value
length_counts <- table(stats$V11)
saveRDS(length_counts, file.path(out_path,paste(sample_name,".readlength_counts.R", sep="")))
