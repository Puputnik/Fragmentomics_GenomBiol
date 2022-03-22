library(vroom)

chr_list <- commandArgs(trailingOnly=TRUE)[1]
stats_path <- commandArgs(trailingOnly=TRUE)[2]
motif_path <- commandArgs(trailingOnly=TRUE)[3]
out_path <- commandArgs(trailingOnly=TRUE)[4]
sample_name <- commandArgs(trailingOnly=TRUE)[5]
debug <- c()

stats_path <- file.path(stats_path,paste(sample_name,".stats", sep=""))
motif_path <- file.path(motif_path,paste(sample_name,".motif", sep=""))

#### loading chromosome list
chr <- read.table(chr_list, stringsAsFactors = F)

#### loading stat files
stats <- as.data.frame(vroom(stats_path, show_col_types = FALSE, col_names = F))
colnames(stats) <- paste("V", seq(1,length(colnames(stats))), sep="")
debug <- c(debug, paste("Stats total reads:", length(rownames(stats))))

#### loading motif files
motif <- as.data.frame(vroom(motif_path, show_col_types = FALSE, col_names = F, delim=" "))
motif$X3 <- toupper(motif$X3)
stats <- cbind(stats, motif)
debug <- c(debug, paste("Read names corresponding:", all(stats$V4 == stats$X1)))

#### keeping reads mapped on selected chromosomes, with MAPQ > 20, readlength < 700 and without 5' hard/soft clipping
stats <- subset(stats, stats$V8 > 20 & stats$V1 %in% chr[,1] & stats$V11 < 700 & ((stats$V6 == "+" & stats$V9 == 0 & stats$V10 == 0) | (stats$V6 == "-" & stats$V12 == 0 & stats$V13 == 0)))
debug <- c(debug, paste("Stats filtered reads:", length(rownames(stats))))

####obtaining 4bp motif counts
counts <- table(stats$X3)

if (!(file.exists(out_path))){
  dir.create(out_path)
}

saveRDS(counts, file.path(out_path,paste(sample_name,".motif.R", sep="")))
write.table(debug, file.path(out_path,paste(sample_name,".motif.log", sep="")), row.names = F, col.names = F, quote = F)

