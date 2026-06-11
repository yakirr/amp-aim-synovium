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
    source('/data/srlab/AMP_collab/lakshay/R/utils.R')
})

start_upR()

args = commandArgs(trailingOnly = TRUE)
ct <- args[1]
xen_path <- args[2]
cohort <- args[3]
chunk <- args[4]
cca_path <- if (nchar(args[5]) == 0) NULL else args[5]
batch_vars <- if (nchar(args[6]) == 0) NULL else stringr::str_split_1(args[6], ",")

# Step 1: Read in single-cell + Xenium datasets 
   
ampp2_ref <- '/data/srlab/AMP_collab/AMP2_2023_seuratObj_galoz/AMPp2_seuratObj_allcells.202511.rds'
sc <- readRDS(ampp2_ref)
names(sc@meta.data)[1] <- 'sid'
ct_for_seu_obj <- gsub("_", " ", ct)
ct_to_keep = sc$celltypes.med == ct_for_seu_obj
sc <- subset(sc, cells = colnames(sc)[ct_to_keep])
print(paste('Number of sc cells:', dim(sc)[2]))

if (!is.null(cca_path)) {
    cca_res <- readRDS(cca_path)
    cca_weights <- cca_res$xcoef
} else {
    cca_weights <- NULL 
}

xen <- readRDS(xen_path)
xen <- subset(xen, features = rownames(xen)[rownames(xen) != 'CCL5'])
print(paste('Number of xenium cells:', dim(xen)[2]))

# Step 2: Find DEGs/broad cell type in AMPp2 ct, and build integrated reference 

sc <- subset(sc, features = intersect(rownames(sc), rownames(xen)))
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = median(sc$nCount_RNA), verbose = F)
res <- FindHVGsFromGroups(sc, 'cluster_name')
hvgs <- res$hvgs

merged <- BuildIntegratedReference(
    xen,
    sc,
    normalization_target = "mean_of_medians",
    hvgs,
    'cluster_name',
    batch_vars = batch_vars, 
    cca_weights = cca_weights 
    )

print('Reference built') 

# Step 3: While there exist Xenium cells without labels, subsample (# of AMPp2 cells in the ct)/iteration and transfer labels from AMPp2 to Xenium, while checking marker gene correlations 

diag_img_dir <- file.path('/data/srlab/AMP_collab/lakshay-yakir/6.fine_types/diagnostic_plots/', cohort, ct_for_seu_obj)
dir.create(diag_img_dir, showWarnings = FALSE, recursive = TRUE)

out_dir = file.path('/data/srlab/AMP_collab/lakshay-yakir/6.fine_types/out_rds/', cohort, ct_for_seu_obj) 
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_filename_stem <- paste0(ct_for_seu_obj, "_", chunk)
save_all_xen_umap_path = file.path(out_dir, paste0(out_filename_stem, "_umap_after_sc_harmonization.rds"))


merged <- TransferLabelsGraph_new(
    merged = merged,
    annotation_var = "cluster_name",
    n_neighbors = 30,
    min_ratio = 1.1, 
    max_mult = 1.3, 
    compute_all_xen_umap = TRUE, 
    save_all_xen_umap_path = save_all_xen_umap_path, 
    save_path = file.path(diag_img_dir, paste0(out_filename_stem, '_label_transfer_diagnostics.png'))
    )

res <- CompareMarkerCorrelations(merged[,merged$modality == "xen"],
                                 merged[,merged$modality == "sc"], "cluster_name", bicluster = F, 
                                 save_path = file.path(diag_img_dir, 
                                                       paste0(out_filename_stem, '_marker_correlations_labeltransfer.png')))

# Step 4: Save results 

colnames(merged) <- sub("^xen_", "", colnames(merged))
saveRDS(merged, file.path(out_dir, paste0(out_filename_stem, '_scandxenium_labeltransfer.rds')))

xen$cluster_name <- merged$cluster_name[
  match(colnames(xen), colnames(merged))
]
xen$cluster_name[!(xen$cluster_name %in% sc$cluster_name)] <- 'Untyped'
saveRDS(xen, file.path(out_dir, paste0(ct_for_seu_obj, "_", chunk, '_labeltransfer.rds')))

percell_cts <- setNames(xen$cluster_name, rownames(xen@meta.data))
saveRDS(percell_cts, file.path(out_dir, paste0(ct_for_seu_obj, '_', chunk, '_allcells_finetypesvector.rds')))

print(paste(sum(xen$cluster_name == 'Untyped'), ct_for_seu_obj, 'cells did not receive fine cell type labels'))
xen_untyped <- subset(xen, cells = colnames(xen)[xen$cluster_name == 'Untyped'])
saveRDS(xen_untyped, file.path(out_dir, paste0(out_filename_stem, '_untyped_labeltransfer.rds')))

marker_corrs_dir <- file.path('/data/srlab/AMP_collab/lakshay-yakir/6.fine_types/out_analysis', cohort, ct_for_seu_obj)
dir.create(marker_corrs_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(res, 
        file.path(marker_corrs_dir, paste0(out_filename_stem, "_tofinetypes_marker_correlations_labeltransfer.RDS")))

print(paste0(Sys.time(), " | Done."))
      
                                 















