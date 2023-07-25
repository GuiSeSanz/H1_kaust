
library(Seurat)
library(ggplot2)
library(edgeR)


h1_data <- readRDS('/home/sevastopol/data/gserranos/H1_kaust/Data/Filtered_data/H1X_GEX_filtered_seurat.rds')
wt_data <- readRDS('/home/sevastopol/data/gserranos/H1_kaust/Data/Filtered_data/WT_GEX_filtered_seurat.rds')

Idents(h1_data) <- 'Sample'
Idents(wt_data) <- 'Sample'

h1_data_counts <-  as.data.frame(h1_data@assays$RNA@counts)
wt_data_counts <-  as.data.frame(wt_data@assays$RNA@counts)

h1_data_pseudo <- as.data.frame(rowSums(h1_data_counts))
wt_data_pseudo <- as.data.frame(rowSums(wt_data_counts))


bulk <- read.csv('/home/sevastopol/data/gserranos/H1_kaust/Data/Bulk_Azari/all_KOs.csv', row.names=1, header=TRUE)
anno <- read.csv('/home/sevastopol/data/gserranos/H1_kaust/Data/Bulk_Azari/anno.csv', row.names=1, header=TRUE)
bulk_anno <- merge(bulk, anno, by=0)


bulk_anno <- bulk_anno[!bulk_anno$Gene_name %in% c( bulk_anno$Gene_name[duplicated(bulk_anno$Gene_name)]), ]


rownames(bulk_anno) <- bulk_anno$Gene_name
bulk_anno <- bulk_anno[, -c(1, ncol(bulk_anno))]

common_names <- intersect(intersect(rownames(h1_data_pseudo), rownames(bulk_anno)), rownames(wt_data_pseudo))



h1_data_pseudo <- h1_data_pseudo[common_names, , drop=FALSE]
wt_data_pseudo <- wt_data_pseudo[common_names, , drop=FALSE]
bulk_anno <- bulk_anno[common_names, ]


bulk_anno_WT <- bulk_anno[, grepl('WT', colnames(bulk_anno))]
bulk_anno_H1 <- bulk_anno[, grepl('H1_X', colnames(bulk_anno))]


# created_dir
dir.create('/home/sevastopol/data/gserranos/H1_kaust/Plots/Correlation_with_Bulk', showWarnings = FALSE)

normalize_TPM <- function(x){
	x <- x/sum(x)*1e6
	return(x)
}

plot_list <- list()
for (i in 1:ncol(bulk_anno_WT)){

	x <- normalize_TPM(bulk_anno_WT[, i])
	y <- normalize_TPM(wt_data_pseudo[, 1])
	corr <- round(cor(x, y, method = "spearman"), 2)

	plot_list[[i]] <- ggplot(data = data.frame(x = -log10(x), y = -log10(y)), aes(x = x, y = y), alpha=0.8) +
	geom_point() + theme_classic() +
	geom_abline(intercept = 0, slope = 1, color = "red") +
	labs(x = "log10 Bulk TPM", y = "log10 Single cell TPM",
		title = paste0(colnames(bulk_anno_WT)[i], "  ", corr)) 
}
png('/home/sevastopol/data/gserranos/H1_kaust/Plots/Correlation_with_Bulk/WT_correlation.png')
cowplot::plot_grid(plotlist = plot_list, ncol = 3)
dev.off()



plot_list <- list()
for (i in 1:9){

	x <- normalize_TPM(bulk_anno_H1[, i])
	y <- normalize_TPM(wt_data_pseudo[,1])
	corr <- round(cor(x, y, method = "spearman"), 2)

	plot_list[[i]] <- ggplot(data = data.frame(x = -log10(x), y = -log10(y)), aes(x = x, y = y), alpha=0.8) +
	geom_point() + theme_classic() +
	geom_abline(intercept = 0, slope = 1, color = "red") +
	labs(x = "log10 Bulk TPM", y = "log10 Single cell TPM",
		title = paste0(colnames(bulk_anno_H1)[i], "  ", corr))
}

png('/home/sevastopol/data/gserranos/H1_kaust/Plots/Correlation_with_Bulk/H1_correlation.png')
cowplot::plot_grid(plotlist = plot_list, ncol = 3)
dev.off()

