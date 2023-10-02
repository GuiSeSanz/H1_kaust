library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)


create_folder <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path)
  }
}

# counts <- Read10X_h5("/home/sevastopol/data/gserranos/H1_kaust/Data/WT_multiome/outs/filtered_feature_bc_matrix.h5")
# fragpath <- "/home/sevastopol/data/gserranos/H1_kaust/Data/WT_multiome/outs/atac_fragments.tsv.gz"
# Sample_Name <- "Test"



Sample_Name <- "H1_3"
counts <- Read10X_h5(paste0("/home/sevastopol/data/gserranos/H1_kaust/Data/Day_0/", Sample_Name, "_multiome/outs/filtered_feature_bc_matrix.h5"))
fragpath <- paste0("/home/sevastopol/data/gserranos/H1_kaust/Data/Day_0/", Sample_Name, "_multiome/outs/atac_fragments.tsv.gz")

print(Sample_Name)

create_folder(paste0('/home/sevastopol/data/gserranos/H1_kaust/Plots/Single_Results/', Sample_Name ))
create_folder(paste0('/home/sevastopol/data/gserranos/H1_kaust/Data/Single_Results/', Sample_Name ))


annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  project = Sample_Name,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(pbmc) <- "RNA"
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)
# pbmc$fragCounts <- CountFragments(fragments = fragpath)
# FRiP(pbmc, "ATAC", total.fragments="fragCounts", col.name = "FRiP", verbose = TRUE)


pdf(paste0('/home/sevastopol/data/gserranos/H1_kaust/Plots/Single_Results/', Sample_Name ,'/QC_prefiltering.pdf'))
VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", 'percent.mt'),
  ncol = 2,
  pt.size = 0
)
dev.off()

pdf(paste0('/home/sevastopol/data/gserranos/H1_kaust/Plots/Single_Results/', Sample_Name ,'/discarded_cells.pdf'))

list_2_plot <- list(
	'nCount_ATAC<0.9' = names(pbmc$nCount_ATAC[pbmc$nCount_ATAC < quantile(pbmc$nCount_ATAC, 0.9)]),
	'nCount_ATAC>0.1' = names(pbmc$nCount_ATAC[pbmc$nCount_ATAC > quantile(pbmc$nCount_ATAC, 0.1)]),
	'nCount_RNA<0.9'  = names(pbmc$nCount_RNA[pbmc$nCount_RNA < quantile(pbmc$nCount_RNA, 0.9)]),
	'nCount_RNA>0.1'  = names(pbmc$nCount_RNA[pbmc$nCount_RNA > quantile(pbmc$nCount_RNA, 0.1)]),
	'nucleosome_signal<2' = names(pbmc$nucleosome_signal[pbmc$nucleosome_signal < 2]),
	'TSS.enrichment>1'    = names(pbmc$TSS.enrichment[pbmc$TSS.enrichment > 1]),
	'mtRatio<0.2'     = names(pbmc$percent.mt[pbmc$percent.mt < 20])
)
# ggvenn::ggvenn(list_2_plot)
# venn::venn(list_2_plot, zcolor = "style")

UpSetR::upset(UpSetR::fromList(list_2_plot), order.by = "freq")
dev.off()



pbmc <- subset(
  x = pbmc,
  subset = 
    nCount_ATAC < quantile(pbmc$nCount_ATAC, 0.9) &
    nCount_ATAC > quantile(pbmc$nCount_ATAC, 0.1) &
    nCount_RNA <  quantile(pbmc$nCount_RNA, 0.9) &
    nCount_RNA >  quantile(pbmc$nCount_RNA, 0.1) &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    percent.mt < 20
)

pdf(paste0('/home/sevastopol/data/gserranos/H1_kaust/Plots/Single_Results/', Sample_Name ,'/QC_postfiltering.pdf'))
VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 2,
  pt.size = 0
)
dev.off()

saveRDS(pbmc, paste0('/home/sevastopol/data/gserranos/H1_kaust/Data/Single_Results/', Sample_Name ,'/', Sample_Name,'_norm_Peaks_RNA.rds'))

# Peak calling
peaks <- CallPeaks(pbmc)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- RunPCA(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


DefaultAssay(pbmc) <- "RNA"
expressed <- c('KLF4', 'SOX2', 'NANOG')
non_expressed <- c('NES', 'PAX6',  'NCAM1', 'FOXA2')
pdf(paste0('/home/sevastopol/data/gserranos/H1_kaust/Plots/Single_Results/', Sample_Name ,'/Gene_Markers.pdf'))
FeaturePlot(pbmc, features = expressed, order=TRUE)
FeaturePlot(pbmc, features = non_expressed, order=TRUE)
dev.off()

saveRDS(pbmc, paste0('/home/sevastopol/data/gserranos/H1_kaust/Data/Single_Results/', Sample_Name ,'/', Sample_Name,'_norm_Peaks_RNA.rds'))


pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

Idents(pbmc) <- "sub.cluster"

p1 <- DimPlot(pbmc, reduction = "umap.rna",  group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap",  group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf(paste0('/home/sevastopol/data/gserranos/H1_kaust/Plots/Single_Results/', Sample_Name ,'/UMAP_integration.pdf'))
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

# saveRDS(pbmc, paste0('/home/sevastopol/data/gserranos/H1_kaust/Data/Single_Results/', Sample_Name ,'/', Sample_Name,'_norm_Peaks_RNA.rds'))

DefaultAssay(pbmc) <- "peaks"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "peaks",
  expression.assay = "RNA"
  # genes.use = c("SOX2", "NANOG")
)

saveRDS(pbmc, paste0('/home/sevastopol/data/gserranos/H1_kaust/Data/Single_Results/', Sample_Name ,'/', Sample_Name,'_norm_Peaks_RNA_Oct.rds'))


pdf(paste0('/home/sevastopol/data/gserranos/H1_kaust/Plots/Single_Results/', Sample_Name ,'/Test.pdf'))
CoveragePlot(
  object = pbmc,
  region = "SOX2",
  features = "SOX2",
  expression.assay = "RNA",
  extend.upstream = 500,
  extend.downstream = 10000
)
dev.off()


# DefaultAssay(pbmc) <- 'RNA'
# rna_barcodes <- colnames(pbmc)
# DefaultAssay(pbmc) <- 'ATAC'
# atac_barcodes <- colnames(pbmc)

