# =============================================================================
# 4.tessera.r
# Run Tessera tiling on Xenium data, cluster tiles, save results
# =============================================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(tessera)
    library(dplyr)
    library(ggthemes)
    library(ggplot2)
    library(tibble)
    library(pheatmap)
    source('../_common/typing_utils.r')
})

samples  <- 'EDP1-EDP2-ARB'
rds_dir  <- paste0("out_rds/", samples)
plot_dir <- paste0("diagnostic_plots/", samples, "/4.tessera")
dir.create(plot_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir,    showWarnings = FALSE, recursive = TRUE)

fast_merge_seurat <- function(seurat_list) {
    message(Sys.time(), " | Merging counts...")
    counts_merged <- do.call(cbind, lapply(seurat_list, function(s) s[["RNA"]]$counts))   
    
    message(Sys.time(), " | Merging metadata...")
    meta_list <- lapply(seurat_list, function(s) s@meta.data)
    meta_merged <- do.call(rbind, unname(meta_list))  # rbind preserves rownames; bind_rows does not
    
    message(Sys.time(), " | Merging PCA...")
    pca_merged <- do.call(rbind, lapply(seurat_list, function(s) s[["pca"]]@cell.embeddings))
    
    message(Sys.time(), " | Creating Seurat object...")
    meta_clean <- meta_merged[, !names(meta_merged) %in% "shape", drop = FALSE]
    merged <- CreateSeuratObject(counts = counts_merged, meta.data = meta_clean)                                        
    merged[["pca"]] <- CreateDimReducObject(
        embeddings = pca_merged,
        key        = "PC_",
        assay      = "RNA"
    )
    message(Sys.time(), " | Done. ", ncol(merged), " cells, ", nrow(merged), " genes.")
    return(list(seurat = merged, shapes = meta_merged[, "shape", drop = FALSE]))
}

# =============================================================================
# Read Xenium object and create tessera_xy reduction
# =============================================================================

message(Sys.time(), " | Reading Xenium RDS...")
xen <- readRDS(file = paste0(rds_dir, "/allcells_qc_harmumapclust_lineage.rds"))

xen[["tessera_xy"]] <- CreateDimReducObject(
    embeddings = as.matrix(xen@meta.data[c('x', 'y')] %>% rename('1' = 'x', '2' = 'y')),
    key   = "TESSERA_",
    assay = DefaultAssay(xen)
)

# =============================================================================
# Run Tessera per cohort
# =============================================================================

cohorts        <- unique(xen@meta.data$cohort)
all_cell_meta  <- list()
all_tiles      <- list()

options(future.globals.maxSize = Inf)

for (coh in cohorts) {
    message(Sys.time(), " | Running Tessera for cohort: ", coh)

    cells_keep <- colnames(xen)[xen@meta.data$cohort == coh]
    xen_coh    <- subset(xen, cells = cells_keep)

    future::plan(future::multicore)
    res <- GetTiles(
        xen_coh,
        'tessera_xy',
        embeddings            = 'harmony',
        group.by              = 'sid',
        dims.use              = 1:25,
        prune_thresh_quantile = 0.99,
        prune_min_cells       = 1,
        max_npts              = 50,
        min_npts              = 5
    )
    future::plan(future::sequential)

    all_cell_meta[[coh]] <- res$obj@meta.data
    all_tiles[[coh]]     <- res$tile_obj

    rm(xen_coh, res)
    gc()
}

# =============================================================================
# Merge cell-to-tile metadata and save
# =============================================================================

message(Sys.time(), " | Merging cell metadata...")
merged_cell_meta <- bind_rows(all_cell_meta)
saveRDS(merged_cell_meta, paste0(rds_dir, "/allcells_qc_harmumapclust_lineage[cellstotiles].rds"))

# =============================================================================
# Merge tile objects, UMAP, cluster, save
# =============================================================================

message(Sys.time(), " | Merging tile Seurat objects...")
result     <- fast_merge_seurat(all_tiles)
tiles      <- result$seurat
tileshapes <- result$shapes
saveRDS(tileshapes, paste0(rds_dir, "/allcells_qc_harmumapclust_lineage[tileshapes].rds"))

set.seed(0)
tiles <- Run_uwot_umap(tiles, reduction = 'pca', spread = 0.8, min_dist = 0.3)
tiles <- FindClusters(tiles, graph.name = 'humap_fgraph', resolution = 0.8, verbose = TRUE)

message(Sys.time(), " | Saving tiles RDS...")
saveRDS(tiles, paste0(rds_dir, "/allcells_qc_harmumapclust_lineage[tiles].rds"))

# =============================================================================
# Diagnostic plots
# =============================================================================

message(Sys.time(), " | Saving diagnostic plots...")

ggsave(paste0(plot_dir, "/umap_tiles_by_cluster.png"),
       DimPlot(tiles, group.by = "seurat_clusters", raster = FALSE) +
           scale_color_tableau("Tableau 20"),
       width = 10, height = 10)

message(Sys.time(), " | Done.")