# =============================================================================
# 2.harmumapclust_markergenes.r
# PCA -> Harmony -> UMAP -> clustering -> marker genes
# =============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(ggthemes)
    library(tidyr)
    library(tibble)
    source('../_common/typing_utils.r')
})

samples <- 'EDP1-EDP2-ARB'

dir.create(paste0("out_analysis/", samples),     showWarnings = FALSE, recursive = TRUE)
dir.create(paste0("diagnostic_plots/", samples), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0("out_rds/", samples), showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# PCA, Harmony, UMAP, cluster
# =============================================================================

message(Sys.time(), " | Reading RDS file...")
seu_qc <- readRDS(file = paste0("out_rds/", samples, "/allcells_qc.rds"))

message(Sys.time(), " | Normalizing...")
seu_qc <- NormalizeData(seu_qc,
                        normalization.method = "LogNormalize",
                        scale.factor = median(seu_qc$nCount_RNA),
                        verbose = FALSE)

message(Sys.time(), " | Finding variable features...")
seu_qc <- FindVariableFeatures(seu_qc, selection.method = "vst", nfeatures = 3000)

message(Sys.time(), " | Scaling data...")
seu_qc <- ScaleData(seu_qc)

message(Sys.time(), " | Running PCA...")
# todo: in future, could add code to subset to ~1M cells prior to PCA and then project loadings back
# onto full dataset. Would run faster and with much less required memory.
set.seed(0)
seu_qc <- RunPCA(seu_qc, verbose = TRUE)

message(Sys.time(), " | Running Harmony...")
set.seed(0)
seu_qc <- harmony::RunHarmony(seu_qc,
                               c("sid", "cohort"),
                               plot_convergence = FALSE,
                               max_iter = 10,
                               early_stop = FALSE)

message(Sys.time(), " | Running UMAP and clustering...")
set.seed(0)
seu_qc <- Run_uwot_umap(seu_qc, spread = 0.8, min_dist = 0.3)
seu_qc <- FindClusters(seu_qc, graph.name = 'humap_fgraph', resolution = 1.5, verbose = TRUE)

message(Sys.time(), " | Dropping non-dense scaled layer...")
seu_qc[["RNA"]]$scale.data <- NULL

message(Sys.time(), " | Saving RDS...")
saveRDS(seu_qc, file = paste0("out_rds/", samples, "/allcells_qc_harmumapclust.rds"))

# =============================================================================
# Diagnostic plots
# =============================================================================

message(Sys.time(), " | Saving diagnostic plots...")
dir.create(paste0("diagnostic_plots/", samples, "/2.harmpcaclust_markergenes"), showWarnings = FALSE, recursive = TRUE)

ggsave(paste0("diagnostic_plots/", samples, "/2.harmpcaclust_markergenes/umap_by_sid.png"),
       DimPlot(seu_qc, group.by = "sid", raster = TRUE),
       width = 10, height = 10)

ggsave(paste0("diagnostic_plots/", samples, "/2.harmpcaclust_markergenes/umap_by_cluster.png"),
       DimPlot(seu_qc, group.by = "seurat_clusters", raster = TRUE) +
           scale_color_tableau("Tableau 20"),
       width = 10, height = 10)

ggsave(paste0("diagnostic_plots/", samples, "/2.harmpcaclust_markergenes/umap_by_cohort.png"),
       DimPlot(seu_qc, group.by = "cohort", raster = FALSE),
       width = 10, height = 10)

m <- list(
    t            = c("CD3E", "CD8A", "CD4", "IL7R", "NKG7", "NCAM1"),
    b.plasma     = c("MS4A1", "CD19", "CD79A", "PRDM1", "MZB1", "XBP1"),
    myeloid      = c("CD14", "CD68", "FCGR3A", "LYVE1", "C1QA", "MARCO"),
    mast         = c("HDC", "KIT"),
    endothelial  = c("PECAM1", "VWF", "CDH5", "FLT1"),
    stromal      = c("PRG4", "CD55", "COL6A1", "CXCL12", "PDPN", "THY1", "PDGFRB", "ACTA2"),
    adipocytes   = c("ADIPOQ", "PLIN1", "PLIN4"),
    proliferating = c("MKI67", "TOP2A"),
    muscles      = c("MYOD1", "CD82")
)

ggsave(paste0("diagnostic_plots/", samples, "/2.harmpcaclust_markergenes/featureplot_markers.png"),
       suppressMessages(
           FeaturePlot(seu_qc, features = unlist(m), ncol = 7, raster = TRUE, order = TRUE) &
               coord_fixed() &
               scico::scale_color_scico(palette = 'batlow', direction = -1)
       ),
       width = 25, height = 30)

# Dot plot of positive control marker genes across clusters
exp_mat <- FetchData(seu_qc, vars = unlist(m), layer = "data")
meta <- seu_qc@meta.data %>%
    select(seurat_clusters) %>%
    bind_cols(exp_mat) %>%
    pivot_longer(-seurat_clusters, names_to = "Gene", values_to = "Expression") %>%
    left_join(
        data.frame(unlist(m)) %>%
            tibble::rownames_to_column("lineage") %>%
            mutate(lineage = gsub("[[:digit:]]", "", lineage)),
        by = join_by(Gene == unlist.m.)
    )

meta_summary <- meta %>%
    group_by(seurat_clusters, Gene, lineage) %>%
    summarise(Avg = mean(Expression),
              Pct = sum(Expression > 0) / length(Expression) * 100,
              .groups = "drop")

dot_plot <- ggplot(meta_summary, aes(x = Gene, y = seurat_clusters)) +
    geom_point(aes(size = Pct, fill = Avg), shape = 21, stroke = NA) +
    scale_size("% detected", range = c(0, 6)) +
    scale_fill_gradientn(colours = viridisLite::mako(100, direction = -1),
                         guide = guide_colorbar(ticks.colour = "black",
                                                frame.colour = "black"),
                         name = "Average\nexpression") +
    ylab("Cluster") + xlab("") +
    facet_grid(. ~ lineage, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 18) +
    theme(axis.text.x  = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
          axis.text.y  = element_text(size = 12, color = "black"),
          axis.title   = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

ggsave(paste0("diagnostic_plots/", samples, "/2.harmpcaclust_markergenes/dotplot_markers_by_cluster.png"),
       dot_plot,
       width = 22, height = 5)

# =============================================================================
# 1. Read RDS
# =============================================================================
message(Sys.time(), " | Saving sub-objects...")
saveRDS(seu_qc[["humap"]]@cell.embeddings,
        paste0("out_rds/", samples, "/allcells_qc_harmumapclust[humap].rds"))
saveRDS(seu_qc[["harmony"]]@cell.embeddings,
        paste0("out_rds/", samples, "/allcells_qc_harmumapclust[hPCs].rds"))
saveRDS(seu_qc@meta.data,
        paste0("out_rds/", samples, "/allcells_qc_harmumapclust[metadata].rds"))

# =============================================================================
# 4. Downsample to 300k cells and save
# =============================================================================

n_cells  <- ncol(seu_qc)
n_target <- 300000
message(Sys.time(), " | Downsampling ", n_cells, " -> ", n_target, " cells...")
set.seed(0)
keep <- sample(colnames(seu_qc), n_target)
seu_ds <- seu_qc[, keep]

message(Sys.time(), " | Saving downsampled RDS...")
saveRDS(seu_ds, paste0("out_rds/", samples, "/allcells_qc_harmumapclust_downsampled.rds"))
message("  Saved: allcells_qc_harmumapclust_downsampled.rds  (", ncol(seu_ds), " cells)")

message(Sys.time(), " | Done.")

# =============================================================================
# Marker genes per cluster in Xenium
# =============================================================================

message(Sys.time(), " | Running wilcoxauc for cluster markers...")
mat <- Seurat::GetAssayData(seu_ds, assay = "RNA", layer = "data")
markers.xen <- presto::wilcoxauc(as.matrix(mat),
                                      seu_ds$seurat_clusters)

write.csv(markers.xen, paste0("out_analysis/", samples, "/xenium_cluster_markers.csv"), row.names = FALSE)

# =============================================================================
# Marker genes per cluster in AMPp2
# =============================================================================

sc <- readRDS('/data/srlab/AMP_collab/AMP2_2023_seuratObj_galoz/AMPp2_seuratObj_allcells.202511.rds')
names(sc@meta.data)[1] <- 'sid'

# filter AMPp2 to the same set of genes in the same order and convert to Seurat
common <- intersect(rownames(sc), VariableFeatures(seu_ds))
sc <- subset(sc, features=common)

# downsample
n_cells  <- ncol(sc)
n_target <- 300000
message(Sys.time(), " | Downsampling ", n_cells, " -> ", n_target, " cells...")
set.seed(0)
keep <- sample(colnames(sc), n_target)
sc <- sc[, keep]

message(Sys.time(), " | Running wilcoxauc on sc cell types...")
mat_sc <- Seurat::GetAssayData(sc, assay = "RNA", layer = "data")
markers.sc <- presto::wilcoxauc(as.matrix(mat_sc), sc@meta.data$cluster_name)

write.csv(markers.sc, paste0("out_analysis/", samples, "/sc_cluster_markers.csv"), row.names = FALSE)

message(Sys.time(), " | Done.")