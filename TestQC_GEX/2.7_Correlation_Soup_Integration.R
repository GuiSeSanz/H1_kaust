library(Seurat)


soup_counts <- list()
for (SAMPLE in c('WT_GEX', 'H1X_GEX')){
	message(SAMPLE)
	sc  = SoupX::load10X(paste0('/home/sevastopol/data/gserranos/H1_kaust/Data/', SAMPLE, '_experiment/outs/'))
	sc  = SoupX::autoEstCont(sc, tfidfMin=0.9, soupQuantile=0.85)
	out = SoupX::adjustCounts(sc)
	soup_counts[[SAMPLE]] <- out
}




seurat_obj_WT_GEX  <- CreateSeuratObject(counts = soup_counts[['WT_GEX']],  project = 'WT_GEX',  min.cells = 3, min.features = 200)
seurat_obj_WT_GEX <-  PercentageFeatureSet(seurat_obj_WT_GEX, pattern = "^MT-", col.name = "percent.mt")
seurat_obj_WT_GEX <- SCTransform(seurat_obj_WT_GEX)
seurat_obj_H1X_GEX <- CreateSeuratObject(counts = soup_counts[['H1X_GEX']], project = 'H1X_GEX', min.cells = 3, min.features = 200)
seurat_obj_H1X_GEX <- PercentageFeatureSet(seurat_obj_H1X_GEX, pattern = "^MT-", col.name = "percent.mt")
seurat_obj_H1X_GEX <- SCTransform(seurat_obj_H1X_GEX)



seurat_merged <- merge(seurat_obj_WT_GEX, y=seurat_obj_H1X_GEX, project = 'merged_WT_H1', add.cell.ids = c('WT_GEX', 'H1X_GEX'))
seurat_merged$Sample <- seurat_merged$orig.ident


seurat_WT_flt <- readRDS('/home/sevastopol/data/gserranos/H1_kaust/Data/Filtered_data/WT_GEX_filtered_seurat.rds')
seurat_H1_flt <- readRDS('/home/sevastopol/data/gserranos/H1_kaust/Data/Filtered_data/H1X_GEX_filtered_seurat.rds')
