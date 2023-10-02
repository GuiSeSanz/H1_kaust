library(Seurat)
library(Signac)

# Integration attempt




# Easy part, Integrate the RNA together
all_obj <- list()
for (name in c('H1_3', 'H1_WT', 'H1X', 'H1X_H1_3')){
	all_obj[[name]] <- readRDS(paste0('/home/sevastopol/data/gserranos/H1_kaust/Data/Single_Results/',name,'/',name,'_norm_Peaks_RNA_Oct.rds'))
}


for (name in c('H1_3', 'H1_WT', 'H1X', 'H1X_H1_3')){
	DefaultAssay(all_obj[[name]]) <- 'RNA'
	
}


all_obj <- lapply(X = all_obj, FUN = function(x) {
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


##### Only mergint the RNA to test Batch effect
merged_obj <- merge(x = all_obj[['H1_WT']], y = c(all_obj[['H1_3']], all_obj[['H1X']], all_obj[['H1X_H1_3']]), 
					add.cell.id = c('H1_WT', 'H1_3', 'H1X', 'H1X_H1_3'),  project='merged')

VariableFeatures(merged_obj[["RNA"]]) <- rownames(merged_obj[["RNA"]]@scale.data)
merged_obj <-  RunPCA(merged_obj, verbose = FALSE)
PCA_dims <- how_many_PCs(merged_obj, 'merged_PCA_dim')
merged_obj <- RunUMAP(merged_obj, reduction = "pca", dims = 1:PCA_dims)


pdf(paste0(PLOT_PATH, 'Umap_merged_samples.pdf'), width=12, height=10)
	DimPlot(merged_obj, reduction = "umap", group.by = "Sample", cols=colors)
dev.off()




####################################################################################
# RUNNING RPCA <- Best when the query and reference datasets do not exibit substantual batch differences  Stuart*, Butler*, et al., Cell 2019 [Seurat V3]
####################################################################################
# select features that are repeatedly variable across datasets for integration 
features <- SelectIntegrationFeatures(object.list = all_obj)
# run PCA on each dataset using these features
all_obj <- lapply(X = all_obj, FUN = function(x) {
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# Identify anchors 
immune.anchors <- FindIntegrationAnchors(object.list = all_obj, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
int_obj <- IntegrateData(anchorset = immune.anchors)
print("RNA integration DONE!")

DefaultAssay(int_obj) <- "integrated"
int_obj <- ScaleData(int_obj, verbose = TRUE)
int_obj <- RunPCA(int_obj, npcs = 20, verbose = TRUE)
int_obj <- RunUMAP(int_obj, reduction = "pca", dims = 1:20)
int_obj <- FindNeighbors(int_obj, reduction = "pca", dims = 1:20)
int_obj <- FindClusters(int_obj, resolution = 0.1)
DimPlot(int_obj, group.by = "orig.ident") | DimPlot(int_obj, group.by = "seurat_clusters", label = TRUE, label.box = TRUE)
ggplot(int_obj@meta.data, aes(x=seurat_clusters, fill=orig.ident)) + geom_bar(position = "fill")
f <- c("FGFR3", "CCND1", "CCND2", "CCND3","MAF", "MAFB", "MYC", "BRAF", "TP53", "TRAF3", "NSD2")
FeaturePlot(int_obj, f, order = TRUE)
Idents(int_obj) <- "orig.ident"
VlnPlot(int_obj, f, stack = TRUE, flip = TRUE)
#saveRDS(int_obj , "./Riney_organized/seurat_objects/int_5samples_RNAonly.rds")

### SCRNA INTEGRATION: cca
#Set default assay to RNA for all objects in list 
for (i in seq_len(length(all_obj))) {
  DefaultAssay(all_obj[[i]]) <- "RNA"
}

# normalize and identify variable features for each dataset independently
all_obj <- lapply(X = all_obj, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration 
features <- SelectIntegrationFeatures(object.list = all_obj)
# run PCA on each dataset using these features
all_obj <- lapply(X = all_obj, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
# Identify anchors 
immune.anchors <- FindIntegrationAnchors(object.list = all_obj, anchor.features = features, reduction = "cca")
# this command creates an 'integrated' data assay
int_obj <- IntegrateData(anchorset = immune.anchors)
saveRDS(int_obj , "./Riney_organized/seurat_objects/int_5samples_RNAonly_CCA.rds")
print("RNA integration DONE!")
int_obj <- readRDS("./Riney_organized/seurat_objects/int_5samples_RNAonly_CCA.rds")

DefaultAssay(int_obj) <- "integrated"
int_obj <- ScaleData(int_obj, verbose = TRUE)
int_obj <- RunPCA(int_obj, npcs = 20, verbose = TRUE)
int_obj <- RunUMAP(int_obj, reduction = "pca", dims = 1:20)
int_obj <- FindNeighbors(int_obj, reduction = "pca", dims = 1:20)
int_obj <- FindClusters(int_obj, resolution = 0.1)
saveRDS(int_obj , "./Riney_organized/seurat_objects/int_5samples_RNAonly_CCA.rds")

png("./Riney_organized/seurat_objects/int_5samples_RNAonly_CCA.png", w = 1500)
DimPlot(int_obj, group.by = "orig.ident") | DimPlot(int_obj, group.by = "seurat_clusters", label = TRUE, label.box = TRUE) | ggplot(int_obj@meta.data, aes(x=seurat_clusters, fill=orig.ident)) + geom_bar(position = "fill")
dev.off()
table(int_obj$integrated_snn_res.0.1)



# Integrate the ATAC together by retaining the same peaks across all datasets

# integrate RNA and ATAC


# We need to create a common set of peaks across all the datasets to be merged.
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(20, 70, 300), end = c(120, 200, 400)))