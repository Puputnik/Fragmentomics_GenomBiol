library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(pals)

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
samples_info <- read.table("~/Fragmentomics/samples_info_heatmap.tsv", sep="\t", header=F, stringsAsFactors = F)
samples_info$V4[which(samples_info$V4 %in% c("HU", "ISPRO"))] <- "Healthy"

#### select motif order
#### endmotifs.txt is a file containing each 4-mer motif in the same order of Chan et al. 2020
motifNames <- read.delim("~/Fragmentomics/endmotifs.txt", header = T, sep = "\t")

#### set output name
outname <- "~/Fragmentomics/heatmap_motif.pdf"

#### create matrix from samples
df <- get_matrix(samples_info, motifNames)

#### create matrix of "normalizedToCounts" values
groupby = "normalizedToCounts"
dt<-setDT(df, key = "samplename")  
perSample = data.table(motif=unique(dt$motif), key="motif")
for (s in unique(df$samplename)) {
  subTab = dt[as.character(s)]
  groupby=paste0("^",groupby,"$")
  ind=grep(groupby,colnames(dt),value = T)
  perSample[as.character(subTab$motif),as.character(s)] = subTab[, ..ind]
}
perSample<-as.data.frame(perSample)

#### set motif order
order<-motifNames$End.motif
perSample<-perSample[match(order, perSample$motif),]

#### select first 25 motifs
mat<-as.matrix(perSample[1:25,2:(ncol(perSample))])
rownames(mat)<-perSample[1:25,]$motif

#### gather annotation informations from "samples_info"
samples_info_2 <- samples_info
rownames(samples_info_2) <- samples_info_2$V2
name<-colnames(mat)

disease<-as.character(samples_info_2[name,"V4"])
library<-as.character(samples_info_2[name,"V6"])
type<- as.character(samples_info_2[name,"V7"])

tumorFrac<-as.data.frame(dt[dt$motif=="CCCA",c("samplename","tumorFrac")])
tumorFrac<-tumorFrac[tumorFrac$samplename %in% colnames(perSample[,2:ncol(perSample)]), ]
tumorFrac<-tumorFrac$tumorFrac
  
#### define figure graphical parameters
title = "Normalized Motif Counts- both" 

xo= 60  #### figure size multiplier for scaling
font = xo*(100/70)
witdth = xo
height = xo*(53/120)
annosi = xo*(9/120)

#### make annotations
bottomann=rowAnnotation(
  show_legend = T,
                        disease = disease,
                        library = library,
                        #type = type,
                        tumorFrac= tumorFrac,
                        simple_anno_size = unit(annosi, "cm"),
                        annotation_name_gp= gpar(fontsize = font),
                        col = list(disease = c(Healthy =  cols25()[14], Cancer = cols25()[15]),
                                   library = c(Illumina = rgb(198/(255),240/(255),198/(255),1), Nanopore = rgb(246/(255),236/(255),200/(255),1), Nanopore_HU="paleturquoise"),
                                   tumorFrac = colorRamp2(c(min(tumorFrac, na.rm = T), max(tumorFrac,na.rm = T)), c("white", "black"))
                                  )
                              
                        )
#### make heatmap
mat <- t(mat)
row_order <- samples_info$V2[with(samples_info, order(factor(samples_info$V6, levels=c("Nanopore_HU","Nanopore", "Illumina")), samples_info$V8))]
row_labels <- samples_info$V3[match(rownames(mat),samples_info$V2)] 

cols <- c("#E7EBFA", "#B997C6", "#824D99" ,"#4E79C4", "#57A2AC" ,"#7EB875", "#D0B440", "#E67F33" ,"#CE2220" ,"#521913")

h1<-Heatmap(mat, 
            heatmap_legend_param = list(
              legend_direction = "horizontal" 
              ), 
            row_labels = row_labels,
            row_order = row_order,
            name = "Motif",
            column_title = paste0("Motif"), 
            column_title_gp = gpar(fontsize = (font+15),fontface = "bold"),
            row_title_gp = gpar(fontsize=(font+15),fontface = "bold"),
            row_names_gp = gpar(fontsize = font),
            column_names_gp = gpar(fontsize=font),
            col = cols,
            row_names_side = "left",
            cluster_rows = FALSE, 
            cluster_columns = FALSE,
            right_annotation = bottomann,
           )  

#### save pdf
pdf(outname, width = witdth, height = height)
draw(h1,
     padding = unit(xo*(c(30, 100, 10, 60)/60), "mm"),
     heatmap_legend_side = "right"
     )  
dev.off()  
