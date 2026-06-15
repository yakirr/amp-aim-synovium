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

lineage_obj_path <- file.path('/data/srlab/AMP_collab/data/early_disease_synovium/xenium/combined/', 
                              cohort, 
                              'lineages', 
                              paste0(lineage, '.rds'))
lineage_obj <- readRDS(lineage_obj_path)

basedir <- '/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/'
coarsetype_vectors <- file.path(basedir, 'out_rds', cohort, lineage, cts, paste0(cts, '_allcells_coarsetypesvector.rds')) 
coarsetype_vectors <- unlist(lapply(coarsetype_vectors, readRDS))
names(coarsetype_vectors) <- gsub("^xen_", "", names(coarsetype_vectors))

coarsetypes <- setNames(rep('Untyped', ncol(lineage_obj)), colnames(lineage_obj))
coarsetypes[names(coarsetype_vectors)] = coarsetype_vectors 
lineage_obj$celltypes.med <- coarsetypes 

# ── Step 1: Restrict both objects to Xenium genes ─────────────────────────────

logmsg("Restricting to Xenium genes...")

xen_genes    <- rownames(lineage_obj)
shared_genes <- intersect(rownames(sc), xen_genes)
logmsg("Shared genes (SC ∩ Xenium panel):", length(shared_genes))

sc_sub  <- sc[shared_genes, ]
xen_sub <- lineage_obj[shared_genes, ]

# ── Step 2: Normalize to the mean of median total gene expression ─────────────

norm_constant <- mean(c(sc_sub$nCount_RNA, xen_sub$nCount_RNA))
sc_sub <- NormalizeData(sc_sub, normalization.method = 'LogNormalize', scale.factor = norm_constant)
xen_sub <- NormalizeData(xen_sub, normalization.method = 'LogNormalize', scale.factor = norm_constant)
logmsg("Normalized to mean of median gene expression.")

# ── Step 3: wilcoxauc on Xenium object ────────────────────────────────────────
logmsg("Running presto::wilcoxauc on Xenium object...")
xen_expr   <- GetAssayData(xen_sub, assay = DefaultAssay(xen_sub), layer = "data")
xen_groups <- xen_sub[["celltypes.med"]][, 1]
xen_auc    <- presto::wilcoxauc(xen_expr, xen_groups)

out_dir <- file.path(basedir, 'out_analysis', cohort, lineage)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(xen_auc, file.path(out_dir, paste0(lineage, "_xen_wilcoxauc_new.rds")))
logmsg("Saved Xenium presto output.")

# ── Step 4: wilcoxauc on sc reference object ───────────────────────────────────
logmsg("Running presto::wilcoxauc on sc reference...")
sc_expr   <- GetAssayData(sc_sub, assay = DefaultAssay(sc_sub), layer = "data")
sc_groups <- sc_sub[["celltypes.med"]][, 1]
sc_auc    <- presto::wilcoxauc(sc_expr, sc_groups)

saveRDS(sc_auc, file.path(out_dir, paste0(lineage, "_sc_wilcoxauc_new.rds")))
logmsg("Saved sc reference presto output.")

# ── Step 5: save # of overlapping marker genes ─────────────────────────────────

top_ct_markers <- function(df, ct, max_rank = 50) {
    df %>% 
        filter(group == ct) %>% 
        arrange(desc(logFC)) %>% 
        mutate(rank = row_number()) %>% 
        filter(rank <= max_rank)
}

for (ct in unique(sc_sub$celltypes.med)) {
    xen_markers <- top_ct_markers(xen_auc, ct)
    sc_markers <- top_ct_markers(sc_auc, ct)
    common_genes <- intersect(xen_markers$feature, sc_markers$feature)
    
    xen_markers %>% 
        filter(feature %in% common_genes) %>% 
        write.csv(file.path(out_dir, paste0(ct, "_overlappingmarkers_xenlogFCs_new.csv")), row.names = FALSE)
    
    sc_markers %>% 
        filter(feature %in% common_genes) %>% 
        write.csv(file.path(out_dir, paste0(ct, "_overlappingmarkers_sclogFCs_new.csv")), row.names = FALSE)
    
    logmsg(paste(length(common_genes), "genes overlap between Xenium and AMPp2", ct, "cells"))
}

# ── Step 6: Compute marker correlations ────────────────────────────────────────

res <- CompareMarkerCorrelations(xen_sub, 
                          sc_sub, 
                          'celltypes.med', 
                          bicluster = F, 
                          markers.xen = xen_auc, 
                          markers.sc = sc_auc, 
                          pheatmap_breaks = seq(-1, 1, length.out = 101), 
                          fontsize = 25, 
                          show = FALSE, 
                          save_path = file.path(basedir, 'diagnostic_plots', cohort, lineage, 
                                                       paste0(lineage, '_allcells_marker_correlations_labeltransfer_new.png')))

saveRDS(res, file.path(out_dir, paste0(lineage, "_tomajorcelltypes_marker_correlations_labeltransfer_new.RDS")))

xen_mat <- xen_auc %>% 
    group_by(group) %>% 
    arrange(desc(logFC), .by_group = TRUE) %>% 
    mutate(rank = row_number()) %>% 
    filter(rank < 11) 
write.csv(xen_mat, file.path(out_dir, paste0(lineage, "_top10_markergenes_percelltype_new.csv")))

logmsg("Saved marker correlations.")