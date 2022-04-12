library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

get_matrix <- function(samples_info){
  c=0
  while (c < length(rownames(samples_info))){
    c=c+1
    path <- samples_info[c,1]
    sname <- samples_info[c,2]
    file_data <- as.data.frame(readRDS(file.path(path,paste(sname,".readlength_counts.R",sep=""))))
    colnames(file_data)[2] <- sname
    file_data[,2] <- as.numeric(file_data[,2]) 
    file_data <- cbind(sname, file_data)
    file_data[,4] <- file_data[,3]/sum(file_data[,3])*100
    file_data[,5] <- samples_info[c,4]
    file_data[,6] <- samples_info[c,7]
    file_data[,7] <- samples_info[c,3]
    file_data[,8] <- samples_info[c,6]
    file_data[,9] <- samples_info[c,8]
    colnames(file_data) <- c("samplename", "readlength","count", "normalizedToCounts", "disease", "type","basename", "library","tumorFrac")
    
    if (c == 1){
      total <- file_data
    } else {
      total <- rbind(total, file_data)
    } 
  }
  return(total)
}
#### load sample annotation file
samples_infob <- read.table("~/Fragmentomics/Utility/samples_info_readlength.tsv", sep="\t", header=F, stringsAsFactors = F)
samples_info <- samples_infob
##### to exclude Illumina and 19_326
#samples_info <- droplevels(samples_info[which(samples_info$V6 != "Illumina"),])
#samples_info <- droplevels(samples_info[which(samples_info$V3 != "19-326"),])

#### create matrix from samples
df <- get_matrix(samples_info)
df$readlength <- as.numeric(as.character(df$readlength))

#### annotate disease field with TF and library info
df$disease[which(df$tumorFrac >= 0.15 & df$disease == "Cancer")] <- "TF>0.15"
df$disease[which(df$tumorFrac < 0.15 & df$disease == "Cancer")] <- "TF<0.15"
df$disease[which(df$library == "Illumina")] <- "Illumina"
df$disease[which(df$basename == "19-326")] <- "19-326"


hashcol <- c("HU" ="#74add1" ,
             "Healthy"    ="#2c7bb6" ,
             "ISPRO"    ="#2c7bb6" ,
             "TF<0.15"    =rgb(211/(255),119/(255),102/(255),1), 
             "TF>0.15"    =rgb(105/(255),58/(255),29/(255),1), 
             "19-326"     = "black"     ,
             "Illumina"=  "#fdc086"         )


#### set output name
#outname <- "~/Fragmentomics/readlength_density.pdf"
outname <- "~/Fragmentomics/readlength_density_all.pdf"

#### print density plot
pdf(outname, 10,5)
p<- ggplot(df,aes(x=readlength,weight=normalizedToCounts, group=samplename, col=disease)) + 
  geom_density(size=0.3, bw=2) + xlim(75,425) + xlab("fragment length")  + 
  scale_color_manual(values = colors <- hashcol) +
  theme(text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size = 3) ) )
print(p)
dev.off()


#### annotate disease field with TF and library info
samples_info <- samples_infob
samples_info$V4[which(samples_info$V8 >= 0.15 & samples_info$V4 == "Cancer")] <- "TF>0.15"
samples_info$V4[which(samples_info$V8 < 0.15  & samples_info$V4 == "Cancer")] <- "TF<0.15"
samples_info$V4[which(samples_info$V6 == "Illumina")] <- "Illumina"
samples_info$V4[which(samples_info$V3 == "19-326")] <- "19-326"

colnames(samples_info) <- c("path","samplename","basename","disease", "run", "library","type","tumorFrac")

get_lengths <- function(samples_info, min_short,max_short,min_long,max_long){
  for (i in rownames(samples_info)){
    dfs <- subset(df, df$samplename == samples_info[i,"samplename"])
    short_count <- sum(dfs$count[which(dfs$readlength >= min_short & dfs$readlength <= max_short)])
    long_count <- sum(dfs$count[which(dfs$readlength >= min_long & dfs$readlength <= max_long)])
    samples_info[i,"ratio"] <- short_count/long_count
  }
  return(samples_info)
}
plotting <- function(data, YL, outname){
  
  #### define printing order
  order_list <- c("Healthy","HU","ISPRO", "TF<0.15", "TF>0.15", "19-326", "Illumina")
  order_list <- order_list[which(order_list %in% data$disease)]
  print(order_list)
  #### plot
  font=20
  
  p<-ggplot(data, aes(x=disease, y=ratio, color=disease))  + geom_quasirandom(dodge.width=0.3, alpha=0.5, bandwidth = 2, varwidth = TRUE, size=4) +
    stat_summary(fun=mean, geom="point", color="black") +
    stat_summary(fun.data=mean_se, geom="errorbar", color="black",width=0.25) +
    scale_color_brewer("",palette="Dark2") +
    theme_classic() +
    ggtitle("Fragment Ratio") + 
    theme(plot.title = element_text(size = 20, face = "bold")) +
    theme(legend.position="none",
          axis.text.x = element_text(size=font*0.9),
          axis.text.y = element_text(size=font),
          axis.title = element_text(size=font),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          
    ) + scale_x_discrete(limits=order_list) + scale_color_manual(values = colors <- hashcol) + ylim(0,YL)
  
  if ("Healthy" %in% data$disease){
    p<- p + stat_compare_means(method = "t.test", size = 8, comparisons = list(c("Healthy", "TF>0.15"))) 
  } else {
    p<- p + stat_compare_means(method = "t.test", size = 8, comparisons = list(c("ISPRO", "TF>0.15"))) 
  }
  
  #### save plot
  pdf(outname, width = 5.5+((length(unique(df$disease))-3)*0.8), height = 5.5)
  print(p  + xlab("  "))
  dev.off()
}

#### get ratio of short fragments (100-150bp) and long fragments (100-220bp) (mononucleosomal peak)
min_short = 100
max_short = 150
min_long = 100
max_long = 220

samples_info_plot <- get_lengths(samples_info, min_short,max_short,min_long,max_long)
#### save mean and sd summary
can <- samples_info_plot$ratio[which(samples_info_plot$disease == "TF>0.15")]
hea <- samples_info_plot$ratio[which(samples_info_plot$disease %in% c("ISPRO", "HU"))]
summary <- rbind(c("Cancers",mean(can, na.rm = T),sd(can, na.rm = T)), c("Healthy",mean(hea, na.rm = T), sd(hea, na.rm = T)))
colnames(summary) <- c("group", "mean", "standard_deviation")
write.table(summary, file="~/Fragmentomics/summary_mononucleosomes.tsv",quote = F, row.names = F, sep="\t")
#### plot
plotting(samples_info_plot, 0.3, "~/Fragmentomics/readlength_jitterplot.pdf")

#### keep data for correlation plot
corr_final <- samples_info_plot[,c("samplename","disease","library","ratio")]
colnames(corr_final)[length(colnames(corr_final))] <- "mononucleosomal"

samples_info_plot$disease[which(samples_info_plot$type == "Healthy")] <- "Healthy"
plotting(samples_info_plot, 0.3, "~/Fragmentomics/readlength_jitterplot_combined.pdf")

#### dinucleosomal peak
min_short = 275
max_short = 325
min_long = 275
max_long = 400

samples_info_plot <- get_lengths(samples_info, min_short,max_short,min_long,max_long)
#samples_info_plot <- droplevels(samples_info_plot[-which(samples_info_plot$disease == "Illumina"),])
#### save mean and sd summary
can <- samples_info_plot$ratio[which(samples_info_plot$disease == "TF>0.15")]
hea <- samples_info_plot$ratio[which(samples_info_plot$disease %in% c("ISPRO", "HU"))]
summary <- rbind(c("Cancers",mean(can, na.rm = T),sd(can, na.rm = T)), c("Healthy",mean(hea, na.rm = T), sd(hea, na.rm = T)))
colnames(summary) <- c("group", "mean", "standard_deviation")
write.table(summary, file="~/Fragmentomics/summary_dinucleosomes.tsv",quote = F, row.names = F, sep="\t")
#### plot
plotting(samples_info_plot, 0.7, "~/Fragmentomics/readlength_jitterplot_long.pdf")

#### keep data for correlation plot
#corr_final <- corr_final[which(corr_final$library != "Illumina"),]
corr_final <- cbind(corr_final, samples_info_plot[,"ratio"])
colnames(corr_final)[length(colnames(corr_final))] <- "dinucleosomal"

samples_info_plot$disease[which(samples_info_plot$type == "Healthy")] <- "Healthy"
plotting(samples_info_plot, 0.7, "~/Fragmentomics/readlength_jitterplot_long_combined.pdf")

#### correlation plot
hashcol <- c("HU" ="#74add1" ,
             "ISPRO"    ="#2c7bb6" ,
             "TF<0.15"    =rgb(211/(255),119/(255),102/(255),1), 
             "TF>0.15"    =rgb(105/(255),58/(255),29/(255),1), 
             "19-326"     = "black"     ,
             "Illumina"=  "#fdc086"         )

p <- ggplot(corr_final, aes(x=mononucleosomal, y=dinucleosomal, col=disease)) + xlim(0,0.67) + ylim(0,0.67) + scale_color_manual(values = colors <- hashcol) + theme(text = element_text(size = 20)) + geom_point(size=2)
pdf("~/Fragmentomics/mononucleosomal_vs_dinucleosomal.pdf",5.5, 4)
print(p)
dev.off()







