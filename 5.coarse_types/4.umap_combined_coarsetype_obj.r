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
lineage <- args[2]
xen_path <- args[3]
cohort <- args[4]
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
res <- FindHVGsFromGroups(sc, 'sid')
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

# Step 3: Subset Xenium cells, build UMAP, and save results 

xen_harm <- subset(merged, cells = colnames(merged)[merged$modality == 'xen'])
rm(merged)
rm(xen)
gc()
xen_harm <- Run_uwot_umap(
        xen_harm,
        reduction = "harmony",
        spread = 0.8,
        min_dist = 0.3
    )

out_dir = file.path('/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/out_rds/', cohort, lineage, ct_for_seu_obj) 
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_filename_stem <- paste0(ct_for_seu_obj, "_")
save_all_xen_umap_path = file.path(out_dir, paste0(out_filename_stem, "xenonly_umap_after_sc_harmonization.rds"))
saveRDS(xen_harm, save_all_xen_umap_path)

print(paste0(Sys.time(), " | Done."))
      
                                 















