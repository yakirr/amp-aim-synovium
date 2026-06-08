suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    source('/data/srlab/AMP_collab/lakshay-yakir/_common/typing_utils.r')
})

args <- commandArgs(trailingOnly = TRUE)
lineage <- args[1]
cohort <- args[2]

# Step 1: read in AMPp2 and all chunk files 

ampp2_ref <- '/data/srlab/AMP_collab/AMP2_2023_seuratObj_galoz/AMPp2_seuratObj_allcells.202511.rds'
sc <- readRDS(ampp2_ref)
names(sc@meta.data)[1] <- 'sid'
lineage_to_keep = sc$lineage == lineage
sc <- subset(sc, cells = colnames(sc)[lineage_to_keep])

cell_types <- unique(sc$celltypes.med)
basedir <- '/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/'
chunk_files <- list.files(file.path(basedir, 'out_rds', cohort, lineage, cell_types), full.names = TRUE, recursive = FALSE)
chunks <- lapply(chunk_files, function(f) {readRDS(f)})
lineage_obj <- merge(chunks[[1]], y = chunks[-1])
lineage_obj <- JoinLayers(lineage_obj)
rm(chunks)
gc()

# Step 2: Compute marker correlations and save all results 
lineage_obj <- NormalizeData(lineage_obj, normalization.method = "LogNormalize", scale.factor = mean(c(lineage_obj$nCount_RNA, sc$nCount_RNA)))
res <- CompareMarkerCorrelations(lineage_obj,
                                 sc, "celltypes.med", bicluster = F, 
                                 save_path = file.path(basedir, 'diagnostic_plots', cohort, lineage,  
                                                       paste0(lineage, '_allcells_marker_correlations_labeltransfer.png')))

for (ct in cell_types) {
    ct_obj <- subset(lineage_obj, subset = celltypes.med == ct)
    logmsg(paste0(length(colnames(ct_obj)), ' cells were typed as ', ct))
    saveRDS(ct_obj, file.path(basedir, 'out_rds', cohort, lineage, ct, paste0(ct, '_allcells_coarsetypes.rds')))
    
    percell_cts <- setNames(ct_obj$celltypes.med, rownames(ct_obj@meta.data))
    saveRDS(percell_cts, file.path(basedir, 'out_rds', cohort, lineage, ct, paste0(ct, '_allcells_coarsetypesvector.rds')))
} 
                 

