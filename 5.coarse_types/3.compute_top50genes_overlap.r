suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    source('/data/srlab/AMP_collab/lakshay-yakir/_common/typing_utils.r')
})

args <- commandArgs(trailingOnly = TRUE)
lineage <- args[1]
cohort  <- args[2]

# ── Step 0: Load data ──────────────────────────────────────────────────────────
ampp2_ref <- '/data/srlab/AMP_collab/AMP2_2023_seuratObj_galoz/AMPp2_seuratObj_allcells.202511.rds'
sc <- readRDS(ampp2_ref)
names(sc@meta.data)[1] <- 'sid'
sc <- subset(sc, cells = colnames(sc)[sc$lineage == lineage])
cts <- unique(sc$celltypes.med)

basedir <- '/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/'
chunk_files <- list.files(file.path(basedir, 'out_rds', cohort, lineage),
                          full.names  = TRUE,
                          recursive   = TRUE,
                          pattern     = paste0(paste0(cts, "_\\d+_labeltransfer.rds"), collapse = "|"))
logmsg("Chunk files found:", length(chunk_files))
chunks  <- lapply(chunk_files, readRDS)
if (length(chunks) == 1) {
    lineage_obj <- chunks[[1]]
} else {
    lineage_obj <- merge(chunks[[1]], y = chunks[-1])
}
lineage_obj  <- JoinLayers(lineage_obj)
lineage_obj <- NormalizeData(lineage_obj, normalization.method = "LogNormalize", scale.factor = median(lineage_obj$nCount_RNA))

# ── Step 1: Compute HVGs from sc ──────────────────────────────────────────────
logmsg("Computing HVGs from sc reference...")
hvg_res <- FindHVGsFromGroups(sc, group_var = 'celltypes.med')
hvgs    <- hvg_res$hvgs


# ── Step 2: Restrict both objects to HVGs ∩ Xenium genes ──────────────────────
xen_genes    <- rownames(lineage_obj)
shared_genes <- intersect(hvgs, xen_genes)
logmsg("Shared genes (HVGs ∩ Xenium panel):", length(shared_genes))

sc_sub  <- sc[shared_genes, ]
xen_sub <- lineage_obj[shared_genes, ]

# ── Step 3: wilcoxauc on Xenium object ────────────────────────────────────────
logmsg("Running presto::wilcoxauc on Xenium object...")
xen_expr   <- GetAssayData(xen_sub, assay = DefaultAssay(xen_sub), layer = "data")
xen_groups <- xen_sub[["celltypes.med"]][, 1]
xen_auc    <- presto::wilcoxauc(xen_expr, xen_groups)

out_dir <- file.path(basedir, 'out_analysis', cohort, lineage)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(xen_auc, file.path(out_dir, paste0(lineage, "_xen_wilcoxauc.rds")))
logmsg("Saved Xenium presto output.")

# ── Step 4: wilcoxauc on sc reference object ───────────────────────────────────
logmsg("Running presto::wilcoxauc on sc reference...")
sc_expr   <- GetAssayData(sc_sub, assay = DefaultAssay(sc_sub), layer = "data")
sc_groups <- sc_sub[["celltypes.med"]][, 1]
sc_auc    <- presto::wilcoxauc(sc_expr, sc_groups)

saveRDS(sc_auc, file.path(out_dir, paste0(lineage, "_sc_wilcoxauc.rds")))
logmsg("Saved sc reference presto output.")

#for (ct in cell_types) {
#    ct_obj <- subset(lineage_obj, subset = celltypes.med == ct)
#    logmsg(paste0(length(colnames(ct_obj)), ' cells were typed as ', ct))
#    saveRDS(ct_obj, file.path(basedir, 'out_rds', cohort, , ct, paste0(ct, '_allcells_coarsetypes.rds')))
    
#    percell_cts <- setNames(ct_obj$celltypes.med, rownames(ct_obj@meta.data))
#    saveRDS(percell_cts, file.path(basedir, 'out_rds', cohort, , ct, paste0(ct, '_allcells_coarsetypesvector.rds')))
#} 
                 

