library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

get_matrix <- function(samples_info, motifs){
  c=0
  while (c < length(rownames(samples_info))){
    c=c+1
    path <- samples_info[c,1]
    sname <- samples_info[c,2]
    file_data <- as.data.frame(readRDS(file.path(path,paste(sname,".motif.R",sep=""))))
    file_data <- subset(file_data, file_data$Var1 %in% motifs$End.motif)
    colnames(file_data)[2] <- sname
    file_data[,2] <- as.numeric(file_data[,2]) 
    file_data <- cbind(sname, file_data)
    file_data[,4] <- file_data[,3]/sum(file_data[,3])*100
    file_data[,5] <- samples_info[c,4]
    file_data[,6] <- samples_info[c,7]
    #file_data[,7] <- sname
    file_data[,7] <- samples_info[c,3]
    file_data[,8] <- samples_info[c,6]
    file_data[,9] <- samples_info[c,8]
    colnames(file_data) <- c("samplename", "motif","count", "normalizedToCounts", "disease", "type","basename", "library","tumorFrac")
    
    if (c == 1){
      total <- file_data
    } else {
      total <- rbind(total, file_data)
    } 
  }
  return(total)
}
#### load sample annotation file
samples_info <- read.table("~/Fragmentomics/Utility/samples_info_heatmap.tsv", sep="\t", header=F, stringsAsFactors = F)

#### select motif order
motifNames <- read.delim("~/Fragmentomics/Utility/endmotifs.txt", header = T, sep = "\t")

#### set output name
outname <- "~/Fragmentomics/CCCA.pdf"

#### create matrix from samples
df <- get_matrix(samples_info, motifNames)

#### annotate disease field with TF and library info
df$disease[which(df$tumorFrac >= 0.15 & df$disease == "Cancer")] <- "TF>0.15"
df$disease[which(df$tumorFrac < 0.15 & df$disease == "Cancer")] <- "TF<0.15"
#df$disease[which(df$type == "Healthy")] <- "Healthy"
df$disease[which(df$library == "Illumina")] <- "Illumina"
df$disease[which(df$basename == "19-326")] <- "19-326"

data = df[which(df$motif == "CCCA"),]

#### define printing order
order_list <- c("Healthy","HU","ISPRO", "TF<0.15", "TF>0.15", "19-326", "Illumina")
order_list <- order_list[which(order_list %in% data$disease)]

hashcol <- c("HU" ="#74add1" ,
             "Healthy"    ="#2c7bb6" ,
             "ISPRO"    ="#2c7bb6" ,
             "TF<0.15"    =rgb(211/(255),119/(255),102/(255),1), 
             "TF>0.15"    =rgb(105/(255),58/(255),29/(255),1), 
             "19-326"     = "black"     ,
             "Illumina"=  "#fdc086"         )

#### plot
x="disease"
y="normalizedToCounts"
font=20

p<-ggplot(data, aes_string(x=x, y=y, color=x)) + xlab(" ") + ylab("Motif %") +
  geom_quasirandom(dodge.width=0.3, alpha=0.5, bandwidth = 2, varwidth = TRUE, size=4) +#, width = 0.1) +
  stat_summary(fun=mean, geom="point", color="black") +
  stat_summary(fun.data=mean_se, geom="errorbar", color="black",width=0.25) +
  scale_color_brewer("",palette="Dark2") +
  theme_classic() +
  ggtitle("CCCA Motif") + 
  theme(plot.title = element_text(size = 20, face = "bold"))+
  theme(legend.position="none",
        axis.text.x = element_text(size=font*0.9),
        axis.text.y = element_text(size=font),
        axis.title = element_text(size=font),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
  ) + scale_x_discrete(limits=order_list) + ylim(0,2.3) + scale_color_manual(values = colors <- hashcol) 

#### add P-value 
if ("Healthy" %in% data$disease){
  p<- p + stat_compare_means(method = "t.test", size = 8, comparisons = list(c("Healthy", "TF>0.15"))) 
} else {
  p<- p + stat_compare_means(method = "t.test", size = 8, comparisons = list(c("ISPRO", "TF>0.15"))) 
}

#### save plot
pdf(outname, width = 5.5+((length(unique(df$disease))-3)*0.8), height = 5.5)
print(p)
dev.off()


can <- data$normalizedToCounts[which(data$disease == "TF>0.15")]
hea <- data$normalizedToCounts[which(data$disease %in% c("ISPRO", "HU"))]
summary <- rbind(c("Cancers",mean(can, na.rm = T),sd(can, na.rm = T)), c("Healthy",mean(hea, na.rm = T), sd(hea, na.rm = T)))
colnames(summary) <- c("group", "mean", "standard_deviation")
write.table(summary, file="~/Fragmentomics/summary_CCCA.tsv",quote = F, row.names = F, sep="\t")
