print(paste0(Sys.time(), " | Starting"))

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(ggthemes)
    library(uwot)
    library(pheatmap)
    library(harmony)
    library(tibble)
    library(tidyr)
    library(rlang)
    library(purrr)
    library(glue)
    library(rhdf5)
    library(stringr)
    library(singlecellmethods)
    source('/data/srlab/AMP_collab/lakshay-yakir/_common/typing_utils.r')
    source('/data/srlab/AMP_collab/lakshay/R/utils.R') # TODO: what functions in this file are necessary? Can I remove this line? 
})

start_upR()

args = commandArgs(trailingOnly = TRUE)
lineage <- args[1]
xen_path <- args[2]
cohort <- args[3]
chunk <- args[4]
batch_vars <- if (nchar(args[5]) == 0) NULL else stringr::str_split_1(args[5], ",")

# Step 1: Read in single-cell + Xenium datasets 
   
ampp2_ref <- '/data/srlab/AMP_collab/AMP2_2023_seuratObj_galoz/AMPp2_seuratObj_allcells.202511.rds'
sc <- readRDS(ampp2_ref)
names(sc@meta.data)[1] <- 'sid'
lineage_to_keep = sc$lineage == lineage
sc <- subset(sc, cells = colnames(sc)[lineage_to_keep])
print(paste('Number of sc cells:', dim(sc)[2]))

xen <- readRDS(xen_path)
xen <- subset(xen, features = rownames(xen)[rownames(xen) != 'CCL5']) # Bad probe (may have already been filtered out during segmentation, but just for completeness' sake) 
print(paste('Number of xenium cells:', dim(xen)[2]))

# Step 2: Find DEGs/broad cell type in AMPp2 ct, and build integrated reference 

sc <- subset(sc, features = intersect(rownames(sc), rownames(xen)))
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = median(sc$nCount_RNA), verbose = F)
res <- FindHVGsFromGroups(sc, 'celltypes.med')
hvgs <- res$hvgs

merged <- BuildIntegratedReference(
    xen,
    sc,
    normalization_target = "mean_of_medians",
    hvgs,
    'celltypes.med',
    batch_vars = batch_vars, 
    cca_weights = NULL 
    )

print('Reference built') 

# Step 3: While there exist Xenium cells without labels, subsample (# of AMPp2 cells in the lineage)/iteration and transfer labels from AMPp2 to Xenium, while checking marker gene correlations 

diag_img_dir <- file.path('/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/diagnostic_plots/', cohort, lineage)
dir.create(diag_img_dir, showWarnings = FALSE, recursive = TRUE)
out_filename_stem <- paste0(lineage, "_", chunk)
save_all_xen_umap_path = file.path('/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/out_rds', 
                                   cohort, 
                                   lineage, 
                                   paste0(out_filename_stem, "_umap_after_sc_harmonization.rds"))


merged <- TransferLabelsGraph_new(
    merged = merged,
    annotation_var = "celltypes.med",
    n_neighbors = 30,
    min_ratio = 1.1, 
    max_mult = 1.3, 
    compute_all_xen_umap = TRUE, 
    save_all_xen_umap_path = save_all_xen_umap_path, 
    save_path = file.path(diag_img_dir, paste0(out_filename_stem, '_label_transfer_diagnostics.png'))
    )

res <- CompareMarkerCorrelations(merged[,merged$modality == "xen"],
                                 merged[,merged$modality == "sc"], "celltypes.med", bicluster = F, 
                                 save_path = file.path(diag_img_dir, 
                                                       paste0(out_filename_stem, '_marker_correlations_labeltransfer.png')))

# Step 4: Save results 

colnames(merged) <- sub("^xen_", "", colnames(merged))

out_dir = file.path('/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/out_rds/', cohort, lineage) 
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(merged, file.path(out_dir, paste0(out_filename_stem, '_scandxenium_labeltransfer.rds')))

xen$celltypes.med <- merged$celltypes.med[
  match(colnames(xen), colnames(merged))
]
xen$celltypes.med[!(xen$celltypes.med %in% sc$celltypes.med)] <- 'Untyped'


for (ct in unique(xen$celltypes.med)) {
    dir.create(file.path(out_dir, ct), showWarnings = FALSE, recursive = TRUE)
    xen_sub <- subset(xen, subset = celltypes.med == ct)
    logmsg(paste(length(colnames(xen_sub)), 'cells were typed as', ct))
    saveRDS(xen_sub, 
            file.path(out_dir, ct, paste0(ct, "_", chunk, '_labeltransfer.rds')))
    
    percell_cts <- setNames(xen_sub$celltypes.med, rownames(xen_sub@meta.data))
    saveRDS(percell_cts, file.path(out_dir, ct, paste0(ct, '_', chunk, '_coarsetypesvector.rds')))
}

print(paste(sum(xen$celltypes.med == 'Untyped'), lineage, 'cells did not receive major cell type labels'))
xen_untyped <- subset(xen, cells = colnames(xen)[xen$celltypes.med == 'Untyped'])
saveRDS(xen_untyped, file.path(out_dir, paste0(out_filename_stem, '_untyped_labeltransfer.rds')))

marker_corrs_dir <- file.path('/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/out_analysis', cohort, lineage)
dir.create(marker_corrs_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(res, 
        file.path(marker_corrs_dir, paste0(out_filename_stem, "_tomajorcelltypes_marker_correlations_labeltransfer.RDS")))

print(paste0(Sys.time(), " | Done."))
      
                                 















