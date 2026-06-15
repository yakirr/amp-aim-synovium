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

args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]
lineage <- args[2]
cohort <- args[3]
cca_path <- if (nchar(args[4]) == 0) NULL else args[4]
batch_vars <- if (nchar(args[5]) == 0) NULL else stringr::str_split_1(args[5], ",")

# Step 1: read in AMPp2 and all chunk files 

ampp2_ref <- '/data/srlab/AMP_collab/AMP2_2023_seuratObj_galoz/AMPp2_seuratObj_allcells.202511.rds'
sc <- readRDS(ampp2_ref)
names(sc@meta.data)[1] <- 'sid'
ct_for_seu_obj <- gsub("_", " ", ct)
ct_to_keep = sc$celltypes.med == ct_for_seu_obj
sc <- subset(sc, cells = colnames(sc)[ct_to_keep])
print(paste('Number of sc cells:', dim(sc)[2]))

basedir <- '/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/'
chunk_files <- list.files(file.path(basedir, 
                                    'out_rds', 
                                    cohort, 
                                    lineage, 
                                    ct_for_seu_obj), 
                          full.names = TRUE, 
                          pattern = '_labeltransfer.rds', 
                          recursive = FALSE)
chunks <- lapply(chunk_files, function(f) {readRDS(f)})
ct_obj <- merge(chunks[[1]], y = chunks[-1])
ct_obj <- JoinLayers(ct_obj)
ct_obj <- subset(ct_obj, features = rownames(ct_obj)[rownames(ct_obj) != 'CCL5']) # Bad probe
rm(chunks)
gc()

if (!is.null(cca_path)) {
    cca_res <- readRDS(cca_path)
    cca_weights <- cca_res$xcoef
} else {
    cca_weights <- NULL 
}

# Step 2: Find DEGs/sid in AMPp2 ct, and build integrated reference 

sc <- subset(sc, features = intersect(rownames(sc), rownames(ct_obj)))
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = median(sc$nCount_RNA), verbose = F)
res <- FindHVGsFromGroups(sc, 'sid')
hvgs <- res$hvgs

merged <- BuildIntegratedReference(
    ct_obj,
    sc,
    normalization_target = "mean_of_medians",
    hvgs,
    'celltypes.med',
    batch_vars = batch_vars, 
    cca_weights = cca_weights, 
    annotate = FALSE
    )
rm(ct_obj)
gc()

print('Reference built')

# Step 3: Compute marker correlations, and generate Xenium-only UMAP 
#res <- CompareMarkerCorrelations(subset(merged, subset = modality == 'xen'),
#                                 subset(merged, subset = modality == 'sc'), 
#                                 "celltypes.med", bicluster = F, 
#                                 save_path = file.path(basedir, 'diagnostic_plots', cohort, lineage, 
#                                                       paste0(lineage, #'_allcells_marker_correlations_labeltransfer.png')))

xen_harm <- subset(merged, cells = colnames(merged)[merged$modality == 'xen'])
rm(merged)
gc()
xen_harm <- Run_uwot_umap(
        xen_harm,
        reduction = "harmony",
        spread = 0.8,
        min_dist = 0.3
    )

# Step 4: save results 

out_dir = file.path('/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/out_rds/', cohort, lineage, ct_for_seu_obj) 
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_filename_stem <- paste0(ct_for_seu_obj, "_")

logmsg(paste0(length(colnames(xen_harm)), ' cells were typed as ', ct))
saveRDS(xen_harm, file.path(out_dir, paste0(out_filename_stem, 'allcells_coarsetypes_umap.rds')))

percell_cts <- setNames(xen_harm$celltypes.med, rownames(xen_harm@meta.data))
saveRDS(percell_cts, file.path(out_dir, paste0(out_filename_stem, 'allcells_coarsetypesvector.rds')))

print(paste0(Sys.time(), " | Done."))


                 

