
library(Seurat)
library(Signac)
library(ggplot2)

normalize_TPM <- function(x){
	x <- x/sum(x)*1e6
	return(x)
}

get_correlation_with_bulk <- function(singleCellSampleName, bulkSampleName, ComparisonName){
	print(ComparisonName)
	single_cell_data <- seurat_obj[[singleCellSampleName]]
	single_cell_data <- as.data.frame(single_cell_data@assays$RNA@counts)
	bulk_data <- bulk_anno[, grep(bulkSampleName, colnames(bulk_anno), value=T), drop=FALSE]

	common_names <- intersect(rownames(single_cell_data), rownames(bulk_anno))

	single_cell_data <- single_cell_data[common_names, , drop=FALSE]
	bulk_data <- bulk_data[common_names, , drop=FALSE]

	single_cell_data_pseudo <- as.data.frame(rowSums(single_cell_data))
	# created_dir
	dir.create('/home/sevastopol/data/gserranos/H1_kaust/Plots/Correlation_with_Bulk', showWarnings = FALSE)

	plot_list <- list()
	for (i in 1:ncol(bulk_data)){
		print(names(bulk_data)[i])
		x <- normalize_TPM(bulk_data[,i, drop=FALSE])
		y <- normalize_TPM(single_cell_data_pseudo[, 1, drop=FALSE])
		corr <- round(cor(x, y, method = "spearman"), 2)

		plot_list[[i]] <- ggplot(data = data.frame(x = -log10(x[,1]), y = -log10(y[,1])), aes(x = x, y = y), alpha=0.8) +
		geom_point() + theme_classic() +
		geom_abline(intercept = 0, slope = 1, color = "red") +
		labs(x = "log10 Bulk TPM", y = "log10 Single cell TPM",
			title = paste0(colnames(bulk_data)[i], "  ", corr)) 
	}
	results <- paste0('/home/sevastopol/data/gserranos/H1_kaust/Plots/Correlation_with_Bulk/', ComparisonName, '.png')
	print(results)
	if(length(plot_list) > 9){plot_list <- plot_list[1:9]}
	png(results)
	cowplot::plot_grid(plotlist = plot_list, ncol = 3)
	dev.off()

}


# get the single cell data
seurat_obj <- list()
seurat_obj['H1_WT']    <- readRDS('/home/sevastopol/data/gserranos/H1_kaust/Data/Single_Results/H1_WT/H1_WT_norm_Peaks_RNA.rds')
seurat_obj['H1X']      <- readRDS('/home/sevastopol/data/gserranos/H1_kaust/Data/Single_Results/H1X/H1X_norm_Peaks_RNA.rds')
seurat_obj['H1_3']     <- readRDS('/home/sevastopol/data/gserranos/H1_kaust/Data/Single_Results/H1_3/H1_3_norm_Peaks_RNA.rds')
seurat_obj['H1X_H1_3'] <- readRDS('/home/sevastopol/data/gserranos/H1_kaust/Data/Single_Results/H1X_H1_3/H1X_H1_3_norm_Peaks_RNA.rds')


# Get the bulk data
bulk <- read.csv('/home/sevastopol/data/gserranos/H1_kaust/Data/Bulk_Azari/all_KOs.csv', row.names=1, header=TRUE)
anno <- read.csv('/home/sevastopol/data/gserranos/H1_kaust/Data/Bulk_Azari/anno.csv', row.names=1, header=TRUE)
bulk_anno <- merge(bulk, anno, by=0)

bulk_anno <- bulk_anno[!bulk_anno$Gene_name %in% c( bulk_anno$Gene_name[duplicated(bulk_anno$Gene_name)]), ]

rownames(bulk_anno) <- bulk_anno$Gene_name
bulk_anno <- bulk_anno[, -c(1, ncol(bulk_anno))]

get_correlation_with_bulk('H1_WT', 'WT', 'WT_correlation_D0')
get_correlation_with_bulk('H1X', 'H1_X', 'H1_X_correlation_D0')
get_correlation_with_bulk('H1_3', 'H1_3', 'H1_3_correlation_D0')
get_correlation_with_bulk('H1X_H1_3', 'H1_X', 'H1X_H1_3_correlation_D0')

singleCellSampleName, bulkSampleName, ComparisonName

for (i in seq(1, length(seurat_obj))){
	print(i)
	print(dim(seurat_obj[[i]]))
}

seurat_obj[['H1_WT']]
wt <- seurat_obj[['H1_WT']]
DefaultAssay(wt) <- 'RNA'
wt[["percent.mt"]] <- PercentageFeatureSet(wt, pattern = "^MT-")
percent.mt < 20

wt <- subset(wt, subset = percent.mt < 20)

plot1 <- FeatureScatter(wt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(paste0('/home/sevastopol/data/gserranos/H1_kaust/Plots/Single_Results/Test.pdf'))
plot1 + plot2 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
FeaturePlot(wt, reduction = 'wnn.umap', features = 'nCount_RNA')
FeaturePlot(wt, reduction = 'wnn.umap', features = 'percent.mt')

dev.off()



DefaultAssay(wt) <- 'peaks'

wt <- LinkPeaks(
  object = wt,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("NANOG", "SOX2")
)

pdf(paste0('/home/sevastopol/data/gserranos/H1_kaust/Plots/Single_Results/Test.pdf'))
CoveragePlot(
  object = seurat_obj[['H1_WT']],
  region = "SOX2",
  features = "SOX2",
  expression.assay = "RNA",
  extend.upstream = 500,
  extend.downstream = 10000
)
dev.off()