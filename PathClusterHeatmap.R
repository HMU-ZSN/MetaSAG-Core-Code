################Draw a clustering heatmap for PathCovCluster.txt####################
####################################################################

library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(dplyr)
library(getopt)

spec <- matrix(c("PathCovClusterFile","i",2,"character","This is PathCovCluster File",
		 "OutPut","o",2,"character","This is the OutPutFolder"),byrow=TRUE,ncol=5)

opt = getopt(spec)
PathCovClusterFile = opt$PathCovClusterFile
OutPutDirt = opt$OutPut


PathCovCluster = read.delim(PathCovClusterFile)
PathCovCluster = acast(PathCovCluster,Path~Cluster,value.var = "Coverage",fill=0)



row_sums=rowSums(PathCovCluster==1)
row_anno <- rowAnnotation(
  bar = anno_barplot(row_sums, gp = gpar(fill = "#5c7ada"), width = unit(1.5, "cm"))
)
col_sums=colSums(PathCovCluster==1)
col_anno <- columnAnnotation(
  bar = anno_barplot(col_sums, gp = gpar(fill = "#5c7ada"), height = unit(1, "cm"))
)



heatmap <- Heatmap(PathCovCluster,
                   name = "Heatmap",
                   right_annotation = row_anno,
                   bottom_annotation = col_anno,
                   col = colorRampPalette(colors = c("white","#b8606c"))(100),
                   rect_gp = gpar(col='grey',lwd=0.5),
                   row_names_gp = gpar(fontsize = 2),
                   column_names_gp = gpar(fontsize = 3),
                   width = unit(2, "mm") * ncol(PathCovCluster),
                   height = unit(0.75, "mm") * nrow(PathCovCluster)

                   )

pdf(paste0(OutPutDirt,'PathCovClusterHeatmap.pdf'))
print(heatmap)
dev.off()
