suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    source('/data/srlab/AMP_collab/lakshay-yakir/_common/typing_utils.r')
})

args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]
ct <- gsub("_", " ", ct) 
cohort  <- args[2]

# ── Step 0: Load data ──────────────────────────────────────────────────────────
ampp2_ref <- '/data/srlab/AMP_collab/AMP2_2023_seuratObj_galoz/AMPp2_seuratObj_allcells.202511.rds'
sc <- readRDS(ampp2_ref)
names(sc@meta.data)[1] <- 'sid'
sc <- subset(sc, cells = colnames(sc)[sc$celltypes.med == ct])

ct_obj <- readRDS('/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/'

basedir <- '/data/srlab/AMP_collab/lakshay-yakir/6.fine_types/'
chunk_files <- list.files(file.path(basedir, 'out_rds', cohort, ct),
                          full.names  = TRUE,
                          recursive   = FALSE,
                          pattern     = paste0(ct, "_\\d+_labeltransfer.rds"))
logmsg("Chunk files found:", length(chunk_files))
chunks  <- lapply(chunk_files, readRDS)
if (length(chunks) == 1) {
    ct_obj <- chunks[[1]]
} else {
    ct_obj <- merge(chunks[[1]], y = chunks[-1])
}
ct_obj  <- JoinLayers(ct_obj)
ct_obj <- NormalizeData(ct_obj, normalization.method = "LogNormalize", scale.factor = median(ct_obj$nCount_RNA))

# ── Step 1: Compute HVGs from sc ──────────────────────────────────────────────
logmsg("Computing HVGs from sc reference...")
hvg_res <- FindHVGsFromGroups(sc, group_var = 'cluster_name')
hvgs    <- hvg_res$hvgs


# ── Step 2: Restrict both objects to HVGs ∩ Xenium genes ──────────────────────
xen_genes    <- rownames(ct_obj)
shared_genes <- intersect(hvgs, xen_genes)
logmsg("Shared genes (HVGs ∩ Xenium panel):", length(shared_genes))

sc_sub  <- sc[shared_genes, ]
xen_sub <- ct_obj[shared_genes, ]

# ── Step 3: wilcoxauc on Xenium object ────────────────────────────────────────
logmsg("Running presto::wilcoxauc on Xenium object...")
xen_expr   <- GetAssayData(xen_sub, assay = DefaultAssay(xen_sub), layer = "data")
xen_groups <- xen_sub[["cluster_name"]][, 1]
xen_auc    <- presto::wilcoxauc(xen_expr, xen_groups)

out_dir <- file.path(basedir, 'out_analysis', cohort, ct)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(xen_auc, file.path(out_dir, paste0(ct, "_xen_wilcoxauc.rds")))
logmsg("Saved Xenium presto output.")

# ── Step 4: wilcoxauc on sc reference object ───────────────────────────────────
logmsg("Running presto::wilcoxauc on sc reference...")
sc_expr   <- GetAssayData(sc_sub, assay = DefaultAssay(sc_sub), layer = "data")
sc_groups <- sc_sub[["cluster_name"]][, 1]
sc_auc    <- presto::wilcoxauc(sc_expr, sc_groups)

saveRDS(sc_auc, file.path(out_dir, paste0(ct, "_sc_wilcoxauc.rds")))
logmsg("Saved sc reference presto output.")

saveRDS(ct_obj, file.path(basedir, 'out_rds', cohort, ct, paste0(ct, '_allcells_finetypes.rds')))

# ── Step 5: wilcoxauc on sc reference object ───────────────────────────────────
logmsg("Computing marker correlations...")                 
corrs <- CompareMarkerCorrelations(xen_sub, 
                                   sc_sub, 
                                   'cluster_name', 
                                   bicluster = FALSE, 
                                   markers.xen = xen_auc, 
                                   markers.sc = sc_auc,   
                                   hvgs_only = FALSE, 
                                   pheatmap_breaks = seq(-1, 1, length.out = 101), 
                                   fontsize = 12, 
                                   show = TRUE, 
                                   save_path = NULL)

saveRDS(corrs, file.path(basedir, 'out_analysis', cohort, ct, paste0(ct, '_tofinetypes_marker_correlations_labeltransfer.RDS')))
