
library(Seurat)
library(ggplot2)

get_upset <- function(data_list){
	tmp <- UpSetR::fromList(data_list)
	plot <- UpSetR::upset(tmp, sets = c(names(data_list)),
	order.by = "freq", empty.intersections = "on")
	# plot <- cowplot::plot_grid(NULL, plot$Main_bar, plot$Sizes, plot$Matrix,
	# 		nrow=2, align='hv', rel_heights = c(3,1),
	# 		rel_widths = c(2,3))
	return(plot)
}




my_h5_file <- Read10X_h5("/home/sevastopol/data/gserranos/H1_kaust/Data/Filtered_data/CellBender/WT_GEX_CellBender_out_filtered.h5")

# Create a Seurat object
CellBender_out <- CreateSeuratObject(counts = my_h5_file)
RawSeurat <- readRDS('/home/sevastopol/data/gserranos/H1_kaust/Data/Filtered_data/WT_GEX_filtered_seurat.rds')



is.mito <- base::grep('^MT', rownames(CellBender_out))
qc <- scater::perCellQCMetrics(CellBender_out@assays$RNA@counts, subsets=list(Mito=is.mito))
# this metric with give us an idea of the complexity of our dataset
CellBender_out$log10GenesPerUMI <- log10(CellBender_out$nFeature_RNA) / log10(CellBender_out$nCount_RNA)
CellBender_out$mitoRatio <- PercentageFeatureSet(object = CellBender_out, pattern = "^MT-")
CellBender_out$RPSRatio  <- PercentageFeatureSet(object = CellBender_out, pattern = "^RP[SL]")
CellBender_out$mitoRatio <- CellBender_out@meta.data$mitoRatio / 100
CellBender_out$RPSRatio  <- CellBender_out@meta.data$RPSRatio / 100



metadata <- setNames(CellBender_out@meta.data, c('orig.ident', 'nUMI', 'nGene' ,'log10GenesPerUMI', 
											'mitoRatio', 'RPSRatio'))
CellBender_out@meta.data <- metadata

CellBender_out_flt <- subset(x = CellBender_out, 
					subset= 
					(nUMI >= 500) & #discard cells with less than 500 reads
					(nGene >= 250) & # discard cells with less than 250 genes
					(nGene <= 5000) & # discard cells with more than 5000 genes (duplets or triplets)
					(log10GenesPerUMI > 0.80) & # discard cells with high number of genes detected per UMI
					(mitoRatio < 0.15)) 

CellBender_out_flt <- NormalizeData(CellBender_out_flt)
CellBender_out_flt <- FindVariableFeatures(CellBender_out_flt, selection.method = "vst", nfeatures = 2000)
CellBender_out_flt <- ScaleData(CellBender_out_flt, features = rownames(CellBender_out_flt))
CellBender_out_flt <- RunPCA(CellBender_out_flt, features = VariableFeatures(object = CellBender_out_flt))
CellBender_out_flt <- FindNeighbors(CellBender_out_flt, dims = 1:20)
CellBender_out_flt <- RunUMAP(CellBender_out_flt, dims = 1:50)
CellBender_out_flt <- FindNeighbors(object = CellBender_out_flt, dims = 1:40)
CellBender_out_flt <- FindClusters(object = CellBender_out_flt, resolution = 1.2)

pdf('/home/sevastopol/data/gserranos/H1_kaust/Plots/CellBender/WT_GEX_CellBender_GeneMarkers.pdf')
FeaturePlot(CellBender_out_flt, features = c("NANOG", "SOX2", "CD24", "CD59", "CD9", "POU5F1")) 
VlnPlot(CellBender_out_flt, features = c("NANOG", "SOX2", "CD24", "CD59", "CD9", "POU5F1"))

FeaturePlot(RawSeurat, features = c("NANOG", "SOX2", "CD24", "CD59", "CD9", "POU5F1")) 
VlnPlot(RawSeurat, features = c("NANOG", "SOX2", "CD24", "CD59", "CD9", "POU5F1"))

get_upset(list(CellBender=colnames(CellBender_out_flt), Raw=colnames(RawSeurat)))
dev.off()


pdf('/home/sevastopol/data/gserranos/H1_kaust/Plots/CellBender/RawNotInCB.pdf')
tmp <- subset(x = RawSeurat, cells= c(colnames(RawSeurat)[!colnames(RawSeurat) %in% colnames(CellBender_out_flt)]))
FeaturePlot(tmp, features = c("NANOG", "SOX2", "CD24", "CD59", "CD9", "POU5F1")) 
VlnPlot(tmp, features = c("NANOG", "SOX2", "CD24", "CD59", "CD9", "POU5F1"), pt.size=0)
dev.off()

pdf('/home/sevastopol/data/gserranos/H1_kaust/Plots/CellBender/RawInCB.pdf')
tmp <- subset(x = RawSeurat, cells= c(colnames(RawSeurat)[colnames(RawSeurat) %in% colnames(CellBender_out_flt)]))
FeaturePlot(tmp, features = c("NANOG", "SOX2", "CD24", "CD59", "CD9", "POU5F1")) 
VlnPlot(tmp, features = c("NANOG", "SOX2", "CD24", "CD59", "CD9", "POU5F1"),pt.size=0)
dev.off()


pdf('/home/sevastopol/data/gserranos/H1_kaust/Plots/CellBender/CellBender_markers.pdf')
FeaturePlot(CellBender_out_flt, features = c("NANOG", "SOX2", "CD24", "CD59", "CD9", "POU5F1")) 
VlnPlot(CellBender_out_flt, features = c("NANOG", "SOX2", "CD24", "CD59", "CD9", "POU5F1"), pt.size=0)
dev.off()



# Correlation with BULK

bulk <- read.csv('/home/sevastopol/data/gserranos/H1_kaust/Data/Bulk_Azari/all_KOs.csv', row.names=1, header=TRUE)
anno <- read.csv('/home/sevastopol/data/gserranos/H1_kaust/Data/Bulk_Azari/anno.csv', row.names=1, header=TRUE)
bulk_anno <- merge(bulk, anno, by=0)

bulk_anno <- bulk_anno[!bulk_anno$Gene_name %in% c( bulk_anno$Gene_name[duplicated(bulk_anno$Gene_name)]), ]

rownames(bulk_anno) <- bulk_anno$Gene_name
bulk_anno <- bulk_anno[, -c(1, ncol(bulk_anno))]

common_genes <- intersect(rownames(bulk_anno), rownames(CellBender_out_flt))
bulk_anno_norm_WT <- bulk_anno[, grepl('WT', colnames(bulk_anno))]

bulk_anno_norm_WT <- bulk_anno_norm_WT[common_genes,]

CellBender_out_flt_pseudo <- as.data.frame(CellBender_out_flt@assays$RNA@counts)
CellBender_out_flt_pseudo <- as.data.frame(rowSums(CellBender_out_flt_pseudo))
CellBender_out_flt_pseudo <- CellBender_out_flt_pseudo[common_genes,, drop=FALSE]


normalize_TPM <- function(x){
	x <- x/sum(x)*1e6
	return(x)
}

plot_list <- list()
for (i in 1:ncol(bulk_anno_norm_WT)){

	x <- normalize_TPM(bulk_anno_norm_WT[, i])
	y <- normalize_TPM(CellBender_out_flt_pseudo[, 1])
	corr <- round(cor(x, y, method = "spearman"), 2)

	plot_list[[i]] <- ggplot(data = data.frame(x = log10(x), y = log10(y)), aes(x = x, y = y), alpha=0.8) +
	geom_point() + theme_classic() +
	geom_abline(intercept = 0, slope = 1, color = "red") +
	labs(x = "log10 Bulk counts", y = "log10 Single cell counts",
		title = paste0(colnames(bulk_anno_norm_WT)[i], "  ", corr)) 
}
png('/home/sevastopol/data/gserranos/H1_kaust/Plots/CellBender/WT_correlation.png')
cowplot::plot_grid(plotlist = plot_list, ncol = 3)
dev.off()

