## Expression analysis
##

#################
## LOAD PACKAGES
#################

library(DESeq2)
library(ggplot2)
library(limma)
library(DSS) 
library(clusterProfiler)
library(ComplexHeatmap)
library(ggbiplot) 
library(caret) 
library(circlize)
library(RColorBrewer)
library(biomaRt) 
library(org.Hs.eg.db) 
library(ReactomePA) 
library(parallel)
library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(venn) 

#############
## FUNCTIONS
#############

## Volcano plot
volcano.plot = function(res, cut.padj, cut.log2FC){
  mut <- as.data.frame(res)
  mutateddf <- mutate(mut, sig=ifelse(mut$padj<cut.padj & mut$log2FoldChange>= cut.log2FC,"UP_regulated", ifelse(mut$padj<cut.padj & mut$log2FoldChange <= -cut.log2FC , "Down_regulated", "Not_different")))
  mutateddf <- na.omit(mutateddf)
  input <- cbind(gene=rownames(mutateddf), mutateddf)
  colnames(input)[8] <- "Significance"
  ggplot(input, aes(log2FoldChange, -log10(pvalue))) + 
    geom_point(aes(col=Significance)) +
    scale_colour_manual(values = c("UP_regulated"= "red", "Down_regulated"="blue",  "Not_different"= "black")) + 
    ggtitle("Volcanoplot DESeq2") +
    theme_bw()
  #volc+geom_text_repel(data=head(input, 20), aes(label=gene))
}

## Intersect
intersect2 <- function(...) {
  args <- list(...)
  nargs <- length(args)
  if(nargs <= 1) {
    if(nargs == 1 && is.list(args[[1]])) {
      do.call("intersect2", args[[1]])
    } else {
      stop("cannot evaluate intersection fewer than 2 arguments")
    }
  } else if(nargs == 2) {
    intersect(args[[1]], args[[2]])
  } else {
    intersect(args[[1]], intersect2(args[-1]))
  }
}


## Scale row for heatmap
scale_rows <- function (x) {
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}
###############################################################################
##############################################################################

tab <- read.table("counts.txt", sep = "\t", header = T, row.names = 1)


# Convert to matrix
countdata <- as.matrix(tab[,1:17])
head(countdata)



# Assign condition (first four are controls, second four contain the expansion)
(condition <- factor(c(rep("SHH", 5), rep("WNT", 5), rep("G3", 5), rep("G4", 5))))

# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)

(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds, parallel = T)


vsd <- vst(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

###########PCA##################

## Plot PCA
p <- p <- plotPCA(vsd)
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
pdf(file = "PCA.pdf", height = 5, width = 10)
print(p)
dev.off()



###########SHH X TODOS

tab <- read.table("counts.txt", sep = "\t", header = T, row.names = 1)

# Assign condition (first four are controls, second four contain the expansion)
(condition <- factor(c(rep("SHH", 5), rep("todos", 15))))

# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)

(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

dds

# Run the DESeq pipeline
dds <- DESeq(dds, parallel = T)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

res <- results(dds, contrast=c("condition","SHH","todos"))
res <- as.data.frame(res)
selected_DEG <- subset(res, abs(res$log2FoldChange) > 2 & res$padj < 0.01)
selected_DEG2 <- selected_DEG[order(selected_DEG$padj), ]
All_DEG_SHH <- rownames(selected_DEG2)

write.table(selected_DEG2, file="selected_DEG2.txt", sep="\t")

All_DEG_SHH_value <- assay(vsd)[All_DEG_SHH, ]
write.table(All_DEG_SHH_value, file="ALL_DEG_SHH_2.txt", sep="\t")



# Read GTF file for gene name conversion
load("GTF_v28.RData")  
convert <- merge(as.data.frame(All_DEG_SHH_value), gtf_v28, by.x="row.names",by.y="gene_id")

########################################
## Heatmap SHHxtodos
########################################

tab_heatmap <- All_DEG_SHH_value
tab_heatmap <- as.matrix(tab_heatmap)
scale_tab <- scale_rows(tab_heatmap)


sampleCondition <- c(rep("SHH", sum(grepl("Shh", colnames(tab_heatmap)))), 
                     rep("WNT", sum(grepl("Wnt", colnames(tab_heatmap)))),
                     rep("Group 3", sum(grepl("G3", colnames(tab_heatmap)))),
                     rep("Group 4", sum(grepl("G4", colnames(tab_heatmap)))))


df <- data.frame(subgroups = sampleCondition, row.names = colnames(tab_heatmap))

ha1 = HeatmapAnnotation(df = df,
                        col = list(subgroups = c("SHH" = "light blue", "WNT" = "purple", "Group 3" = "pink", "Group 4" = "green")))


## Heatmap draw
breaks <- seq(-2,2, by= 0.1)

ht1 <- Heatmap(scale_tab, 
               name = "zscore", column_title = "", width = 1,
               top_annotation = ha1,
               show_row_names = F, show_column_names = F,
               cluster_rows = T, cluster_columns = F,
               col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
               show_column_dend = T, show_row_dend = T, 
               clustering_distance_columns = "pearson", clustering_method_columns = "ward.D2", 
               clustering_distance_rows = "pearson", clustering_method_rows = "ward.D2" 
)
print(ht1)
dev.copy(tiff, "SHH.png", width=8, height=6, res = 500, units = "in")
dev.off()

########################################
## volcano plot SHHxtodos
########################################
png(filename = "ALL_DEG_SHH_Volcano.png", width=10, height=9, res = 500, units = "in")

##pdf(file = "ALL_DEG_SHH_Volcano.pdf")
volcano.plot(res, 0.01, 2)
dev.off()


#################CEMiTool - KEGG########################

## Read GTF file for gene name conversion
load("GTF_v28.RData")  


tab <- read.table("ALL_DEG_SHH_2.txt", sep = "\t", header = T, row.names = 1)
convert <- merge(as.data.frame(tab), gtf_v28, by.x="row.names",by.y="gene_id")

# Adjust a DF with the required format by CEMItool
expressions <- tab[,1:17]
convert <- merge(as.data.frame(tab), gtf_v28, by.x="row.names",by.y="gene_id")
.rowNamesDF(expressions, make.names=TRUE) <- convert$gene_name


## Load packages
library(CEMiTool)


##


# Prepare sample anontation dataframe
annot = data.frame(sampleName = row.names(design), Class = design$subType)

# Read GMT file
gmt_fname <- system.file("extdata", "c2.cp.kegg.v6.2.symbols.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# Read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)

# Creates the cemitool object and perform analysis
cem <- cemitool(as.data.frame(expressions), 
                force_beta = T,
                filter = F,
                apply_vst = F, 
                gmt = gmt_in, 
                interactions = int_df,
                verbose = T)

#
print(cem)

directory <- "<path>"
setwd(directory)

# Create report as pdf and html documents
generate_report(cem, directory = "./Report_SHH", 
                output_format = "html_document", 
                force = T)

# Write analysis results into files
write_files(cem, directory="Report_SHH/Tables", force=TRUE)

# Save all plots
save_plots(cem, "all", directory="Report_SHH/Plots", force = T)


####fazer o mesmo para DEGs em WNTxtodos, G3xtodos e G4xtodos obedecendo a ordem das amostras no arquivo counts######

