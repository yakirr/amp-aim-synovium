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

ct_obj_path <- file.path('/data/srlab/AMP_collab/data/early_disease_synovium/xenium/combined/', 
                          cohort, 
                          'coarsetypes', 
                          paste0(ct, '.rds'))
ct_obj <- readRDS(ct_obj_path)

basedir <- '/data/srlab/AMP_collab/lakshay-yakir/6.fine_types/'
finetype_vector <- file.path(basedir, 'out_rds', cohort, ct, paste0(ct, '_allcells_finetypesvector.rds')) 
finetype_vector <- readRDS(finetype_vector)

finetypes <- setNames(rep('Untyped', ncol(ct_obj)), colnames(ct_obj))
finetypes[names(finetype_vector)] = finetype_vector 
ct_obj$cluster_name <- finetypes

# ── Step 1: Restrict both objects to Xenium genes ─────────────────────────────

logmsg("Restricting to Xenium genes...")

xen_genes    <- rownames(ct_obj)
shared_genes <- intersect(rownames(sc), xen_genes)
logmsg("Shared genes (SC ∩ Xenium panel):", length(shared_genes))

sc_sub  <- sc[shared_genes, ]
xen_sub <- ct_obj[shared_genes, ]

# ── Step 2: Normalize to the mean of median total gene expression ─────────────

norm_constant <- mean(c(sc_sub$nCount_RNA, xen_sub$nCount_RNA))
sc_sub <- NormalizeData(sc_sub, normalization.method = 'LogNormalize', scale.factor = norm_constant)
xen_sub <- NormalizeData(xen_sub, normalization.method = 'LogNormalize', scale.factor = norm_constant)
logmsg("Normalized to mean of median gene expression.")

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

# ── Step 5: save # of overlapping marker genes ─────────────────────────────────

top_ct_markers <- function(df, ct, max_rank = 50) {
    df %>% 
        filter(group == ct) %>% 
        arrange(desc(logFC)) %>% 
        mutate(rank = row_number()) %>% 
        filter(rank <= max_rank)
}

for (cs in unique(sc_sub$cluster_name)) {
    xen_markers <- top_ct_markers(xen_auc, cs)
    sc_markers <- top_ct_markers(sc_auc, cs)
    common_genes <- intersect(xen_markers$feature, sc_markers$feature)
    
    xen_markers %>% 
        filter(feature %in% common_genes) %>% 
        write.csv(file.path(out_dir, paste0(ct, "_overlappingmarkers_xenlogFCs_new.csv")), row.names = FALSE)
    
    sc_markers %>% 
        filter(feature %in% common_genes) %>% 
        write.csv(file.path(out_dir, paste0(ct, "_overlappingmarkers_sclogFCs_new.csv")), row.names = FALSE)
    
    logmsg(paste(length(common_genes), "genes overlap between Xenium and AMPp2", cs, "cells"))
}

# ── Step 5: Compute marker correlations ────────────────────────────────────────

res <- CompareMarkerCorrelations(xen_sub, 
                          sc_sub, 
                          'cluster_name', 
                          bicluster = F, 
                          markers.xen = xen_auc, 
                          markers.sc = sc_auc, 
                          pheatmap_breaks = seq(-1, 1, length.out = 101), 
                          fontsize = 25, 
                          show = FALSE, 
                          save_path = file.path(basedir, 'diagnostic_plots', cohort, ct, 
                                                       paste0(ct, '_allcells_marker_correlations_labeltransfer.png')))

saveRDS(res, file.path(out_dir, paste0(ct, "_tofinetypes_marker_correlations_labeltransfer.RDS")))

xen_mat <- xen_auc %>% 
    group_by(group) %>% 
    arrange(desc(logFC), .by_group = TRUE) %>% 
    mutate(rank = row_number()) %>% 
    filter(rank < 11) 
saveRDS(xen_mat, file.path(out_dir, paste0(ct, "_top10_markergenes_percellstate.csv")))

logmsg("Saved marker correlations.")
