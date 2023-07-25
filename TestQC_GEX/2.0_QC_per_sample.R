library(Seurat)
library(ggplot2)


plot_stats <- function(obj, pdf_name){
  pdf(paste0(PLOT_PATH, pdf_name,'.pdf'))
    # shows the distribution of the transcripts per cells
    print(ggplot(obj@meta.data,aes(color=Phenotype, x=nUMI, fill= Phenotype)) + 
        geom_density(alpha = 0.2) + 
        scale_x_log10() + 
        theme_classic() +
        ylab("Cell density") +
        geom_vline(xintercept = 500) + ggtitle('UMISpercell'))
    print(ggplot(obj@meta.data, aes(color=Phenotype, x=nGene, fill= Phenotype)) + 
        geom_density(alpha = 0.2) + 
        theme_classic() +
        scale_x_log10() + 
        geom_vline(xintercept = 300) + ggtitle('GENESpercell'))
    # Visualize the distribution of genes detected per cell via boxplot
    print(ggplot(obj@meta.data, aes(x=Phenotype, y=log10(nGene), fill=Phenotype)) + 
        geom_boxplot() + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold")) +
        ggtitle("NCells vs NGenes"))
    # correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
    p <- ggplot(obj@meta.data, aes(x=nUMI, y=nGene, color=mitoRatio)) + 
        geom_point() + 
        scale_colour_gradient(low = "gray90", high = "black", limits=c(0,1)) +
        stat_smooth(method=lm) +
        scale_x_log10() + 
        scale_y_log10() + 
        theme_classic() +
        geom_vline(xintercept = 500) +
		geom_hline(yintercept = 5000) +
        geom_hline(yintercept = 250) + ggtitle(paste0('UMIvsGENESpercell  Ncell:: ', ncol(obj)))
    print(p)
    print(p + facet_wrap(~Phenotype))
    print(VlnPlot(obj, features = c("nGene", "nUMI", "mitoRatio", "RPSRatio"), ncol = 4, pt.size=0))
  dev.off()
}

how_many_PCs <- function(obj, pdf_name){
  # determine the correct PC number
  # Determine percent of variation associated with each PC
  pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))
  # Elbow plot to visualize 
  pdf(paste0(PLOT_PATH ,pdf_name,'.pdf'))
  print(ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw())
  dev.off()
  return(pcs)
}


get_CellCycleGenes_per_PCs <- function(PCAs){
  plotter <- data.frame(PC=NULL, CellCycleGenes=NULL)
  for(pc in colnames(PCAs)[1:20]){
    b <- names(sort(abs(PCAs[,pc]), decreasing=TRUE)[1:50])
    plotter <- rbind(plotter,data.frame(PC=pc, CellCycleGenes=sum(b %in% c(g2m_genes, s_genes))))
  }
  p <- ggplot(plotter, aes(x=PC, y=CellCycleGenes, fill=PC)) + geom_bar(alpha=0.8, stat="identity") + 
  ggprism::theme_prism() + theme(legend.position='none', axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values=c(ggthemes::tableau_color_pal('Classic 20')(20) )) + ylim(0,25)
  return(p)
}
get_cell_distribution <- function(dataset){
	n_clusters <- nrow(dataset)
	if( n_clusters < 20){
		colors <- ggthemes::tableau_color_pal('Classic 20')(n_clusters)
	}else{
		colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(n_clusters)
	}
	p <- ggplot(dataset, aes(x=Cluster, y=Ncells, fill=Cluster)) + geom_bar(alpha=0.8, stat="identity") + 
	ggprism::theme_prism()  + theme(legend.position='none') + scale_fill_manual(values=c(colors))
	return(p)
}

get_cells_per_cluster <- function(obj){
  n_cells <- FetchData(obj, vars = c(paste0('SCT_snn_res.', seq(0.2,1.2, by=0.2)))) 
  n_cells <- apply(n_cells, 2, table)
  results <- data.frame(Res=NULL, Cluster=NULL, Ncells=NULL)
  for (name in names(n_cells)){
    tmp <- setNames(as.data.frame(t(n_cells[[name]])), c('Res', 'Cluster', 'Ncells'))
    tmp$Res <- name
    results <- rbind(results, tmp)
  }
  return(results)
}



s_genes   <- c("UBR7","RFC2","RAD51","MCM2","TIPIN","MCM6","UNG","POLD3","WDR76",
"CLSPN","CDC45","CDC6","MSH2","MCM5","POLA1","MCM4","RAD51AP1","GMNN","RPA2",
"CASP8AP2","HELLS","E2F8","GINS2","PCNA","NASP","BRIP1","DSCC1","DTL","CDCA7",
"CENPU","ATAD2","CHAF1B","USP1","SLBP","RRM1","FEN1","RRM2","EXO1","CCNE2",
"TYMS","BLM","PRIM1","UHRF1")
g2m_genes <- c("NCAPD2","ANLN","TACC3","HMMR","GTSE1","NDC80","AURKA","TPX2",
"BIRC5","G2E3","CBX5","RANGAP1","CTCF","CDCA3","TTK","SMC4","ECT2","CENPA",
"CDC20","NEK2","CENPF","TMPO","HJURP","CKS2","DLGAP5","PIMREG","TOP2A","PSRC1",
"CDCA8","CKAP2","NUSAP1","KIF23","KIF11","KIF20B","CENPE","GAS2L3","KIF2C",
"NUF2","ANP32E","LBR","MKI67","CCNB2","CDC25C","HMGB2","CKAP2L","BUB1","CDK1",
"CKS1B","UBE2C","CKAP5","AURKB","CDCA2","TUBB4B","JPT1")


SAMPLE = 'WT_GEX'
DATA_PATH_FILTERED = paste0('/home/sevastopol/data/gserranos/H1_kaust/Data/', SAMPLE, '_experiment/outs/filtered_feature_bc_matrix')

Run_SoupX = FALSE
if (Run_SoupX){
	sc  = SoupX::load10X(paste0('/home/sevastopol/data/gserranos/H1_kaust/Data/', SAMPLE, '_experiment/outs/'))
	sc  = SoupX::autoEstCont(sc, tfidfMin=0.9, soupQuantile=0.85)
	out = SoupX::adjustCounts(sc)
}


PLOT_PATH = '/home/sevastopol/data/gserranos/H1_kaust/Plots/QC/'
PLOT_PATH = paste0(PLOT_PATH, SAMPLE, '/')
# create directory
dir.create(PLOT_PATH, showWarnings = FALSE)

seurat_obj <- CreateSeuratObject(counts = Read10X(data.dir = DATA_PATH_FILTERED), project = SAMPLE, min.cells = 3, min.features = 200)

# Calculate QC metrics
is.mito <- base::grep('^MT', rownames(seurat_obj))
qc <- scater::perCellQCMetrics(seurat_obj@assays$RNA@counts, subsets=list(Mito=is.mito))
# this metric with give us an idea of the complexity of our dataset
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
seurat_obj$RPSRatio  <- PercentageFeatureSet(object = seurat_obj, pattern = "^RP[SL]")
seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
seurat_obj$RPSRatio  <- seurat_obj@meta.data$RPSRatio / 100
seurat_obj$Phenotype <- seurat_obj$orig.ident
seurat_obj$Sample <- SAMPLE
metadata <- setNames(seurat_obj@meta.data, c('orig.ident', 'nUMI', 'nGene' ,'log10GenesPerUMI', 
											'mitoRatio', 'RPSRatio', 'Phenotype', 'Sample'))
seurat_obj@meta.data <- metadata
message('Prefilter plot')
plot_stats(seurat_obj, paste0(SAMPLE, 'Prefilter'))
# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = seurat_obj, 
					subset= 
					(nUMI >= 500) & #discard cells with less than 500 reads
					(nGene >= 250) & # discard cells with less than 250 genes
					(nGene <= 5000) & # discard cells with more than 5000 genes (duplets or triplets)
					(log10GenesPerUMI > 0.80) & # discard cells with high number of genes detected per UMI
					(mitoRatio < 0.15)) 
					#Marina has set the 10% threshold to instead 20
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
message('Post plot')
plot_stats(filtered_seurat, paste0(SAMPLE, 'Postfilter'))
# this normalization is only used to check the influence of the cell cycle phase
seurat_phase <- NormalizeData(filtered_seurat)

seurat_phase <- CellCycleScoring(seurat_phase, 
					g2m.features = g2m_genes, 
					s.features = s_genes)

seurat_phase <- FindVariableFeatures(seurat_phase, 
					selection.method = "vst",
					nfeatures = 2000, 
					verbose = FALSE)
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

PCA_dims <- how_many_PCs(seurat_phase, paste0(SAMPLE, 'NumPCs'))
filtered_seurat <- RunUMAP(seurat_phase, reduction = "pca", dims = 1:50)
message('cellCycle plot')
pdf(paste0(PLOT_PATH ,SAMPLE, '_OriginCellCycle.pdf'))
# Loadings(seurat_phase[['pca']])
# Plot the PCA colored by cell cycle phase
print(DimPlot(seurat_phase,
		reduction = "pca",
		group.by= "Phase"))
print(DimPlot(seurat_phase,
		reduction = "pca",
		group.by= "Phase", split.by='Phenotype'))
print(DimPlot(seurat_phase,
		reduction = "pca",
		group.by= "Phase",
		split.by = "Phase"))
print(get_CellCycleGenes_per_PCs(Loadings(seurat_phase[['pca']])))
print(FeaturePlot(filtered_seurat, features = c( 'S.Score', 'G2M.Score'), order = TRUE) & NoLegend())
dev.off()


filtered_seurat <- SCTransform(filtered_seurat)
filtered_seurat <- RunPCA(object = filtered_seurat)
pdf(paste0(PLOT_PATH ,SAMPLE, '_dims2umap.pdf'), width=12, height=10)
	print(ElbowPlot(filtered_seurat))
	print(DimHeatmap(filtered_seurat, dims = 1:10, cells = 500, balanced = TRUE))
	print(DimHeatmap(filtered_seurat, dims = 11:20, cells = 500, balanced = TRUE))
	print(VizDimLoadings(filtered_seurat, dims = 1:10, reduction = "pca"))
	print(VizDimLoadings(filtered_seurat, dims = 11:20, reduction = "pca"))
dev.off()

message('get PCA dims')
PCA_dims <- how_many_PCs(filtered_seurat, paste0(SAMPLE, 'NumPCs'))
message('get Neighbors')
filtered_seurat <- FindNeighbors(object = filtered_seurat, dims = 1:40)
# Determine the clusters for various resolutions                                
filtered_seurat <- FindClusters(object = filtered_seurat, resolution = c(0.2,0.4, 0.6, 0.8, 1.0, 1.2))
message('RUN UMAP')
filtered_seurat <- RunUMAP(filtered_seurat, reduction = "pca", dims = 1:PCA_dims)

pdf(paste0(PLOT_PATH ,SAMPLE, '_Umap_Stats.pdf'), height = 8)
	# print(clustree::clustree(filtered_seurat@meta.data, prefix = "SCT_snn_res."))
	print(DimPlot(filtered_seurat,  reduction = "umap", group.by = "orig.ident",  cols = c(ggthemes::tableau_color_pal('Classic 10 Medium')(3))) ) 
	features <-  c('nGene', 'nUMI', 'mitoRatio', 'RPSRatio')
	print(FeaturePlot(filtered_seurat, features = features, order = TRUE) & NoLegend())
dev.off()

# message('Res plot')
# pdf(paste0(PLOT_PATH ,SAMPLE, '_allResUmaps.pdf'), height = 8)
# populations <- get_cells_per_cluster(filtered_seurat)
# for (res in (paste0('SCT_snn_res.', seq(0.2,1.2, by=0.2)))){
# 	message(res)
# 	Idents(object = filtered_seurat) <- res
# 	n_clusters <- length(levels(populations[populations$Res == res,'Cluster']))
# 	if( n_clusters < 20){
# 		colors <- ggthemes::tableau_color_pal('Classic 20')(n_clusters)
# 	}else{
# 		colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(n_clusters)
# 	}
# 	print(
# 		cowplot::plot_grid(
# 		DimPlot(filtered_seurat,
# 			reduction = "umap",
# 			label = TRUE,
# 			label.size = 6,
# 			cols= colors) + ggtitle(paste0('Resolution::',res)),
# 		get_cell_distribution(populations[populations$Res == res,]),
# 		nrow=2, rel_heights = c(3,1)
# 		)
# 	)
# }
# dev.off()

saveRDS(filtered_seurat, paste0('/home/sevastopol/data/gserranos/H1_kaust/Data/Filtered_data/' , SAMPLE, '_filtered_seurat.rds'))