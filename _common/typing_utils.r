library(dplyr)
library(tidyr)
library(presto)
library(Seurat)
library(pheatmap)
library(gridExtra)

logmsg <- function(...) {
  cat(..., "\n"); flush.console()
}

fig.size <- function(h, w) {
    options(repr.plot.height = h, repr.plot.width = w)
}

dimlist<-function(list_name){
    map(list_name, ~ dim(.x))
}

small_theme <- theme(
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 7),
    plot.title   = element_text(size = 9)
    )

Run_uwot_umap <- function(SeuratObj, reduction = 'harmony', min_dist = 0.01, spread = 1){
    HU <- uwot::umap(SeuratObj@reductions[[reduction]]@cell.embeddings, min_dist = min_dist, 
                 spread = spread, ret_extra = 'fgraph', fast_sgd = FALSE)
    colnames(HU$embedding) = c('HUMAP1', 'HUMAP2')
    rownames(HU$fgraph) = colnames(HU$fgraph) = Cells(SeuratObj)
    SeuratObj[['humap']] <- Seurat::CreateDimReducObject(
        embeddings = HU$embedding,
        assay = 'RNA',
        key = 'HUMAP_',
        global = TRUE
    )
    HU_graph <- Seurat::as.Graph(HU$fgraph)
    DefaultAssay(HU_graph) <- DefaultAssay(SeuratObj)
    SeuratObj[['humap_fgraph']] <- HU_graph
    return(SeuratObj)
}

transfer_labels_graph <- function(graph, query_mask, ref_mask, labels, 
                                  unweighted = FALSE, min_ratio = 1.1) {
  # graph: adjacency matrix (cells x cells)
  # query_mask: boolean mask specifying query cells
  # ref_mask: boolean mask specifying reference cells
  # labels: vector of cell type labels
  # unweighted: if TRUE, convert graph to unweighted (majority vote)

  library(Matrix)

  # Subset adjacency: query -> reference
  A <- graph[query_mask, ref_mask, drop = FALSE]

  if (unweighted) {
    A@x[] <- 1
  }

  ref_labels <- labels[ref_mask]
  label_levels <- sort(unique(ref_labels[!is.na(ref_labels)]))

  # Build sparse one-hot matrix
  label_mat <- sparseMatrix(
    i = seq_along(ref_labels),
    j = match(ref_labels, label_levels),
    x = 1,
    dims = c(length(ref_labels), length(label_levels)),
    dimnames = list(seq_along(ref_labels), label_levels)
  )

  # Neighbor voting
  vote_counts <- A %*% label_mat

  # Assign label with most votes, requiring separation from runner-up
  transferred <- apply(vote_counts, 1, function(v) {
    if (all(v == 0)) return(NA_character_)
  
    ord <- order(v, decreasing = TRUE)
    top1 <- v[ord[1]]
    top2 <- if (length(v) > 1) v[ord[2]] else 0
  
    if (top2 == 0 || (top1 / top2) >= min_ratio) {
      label_levels[ord[1]]
    } else {
      NA_character_
    }
  })

  return(transferred)
}

FindHVGsFromGroups <- function(sc, group_var, logFC_thresh = 0.1) {
    expr <- GetAssayData(sc, assay = DefaultAssay(sc), layer = "data")
    groups <- sc[[group_var]][,1]
    clusters.markers <- presto::wilcoxauc(expr, groups)

    info <- clusters.markers %>%
        dplyr::group_by(group) %>%
        dplyr::filter(padj < 0.05 & auc > 0.5) %>%
        dplyr::arrange(group, dplyr::desc(logFC))
    hvgs <- info %>%
        dplyr::filter(logFC > logFC_thresh)
    hvgs <- unique(hvgs$feature)

    logmsg("Number of HVGs found:", length(hvgs))

    return(list(
        hvgs = hvgs,
        info = info
    ))
}

BuildReference <- function(
    sc, 
    hvgs, 
    batch_vars = NULL, 
    annotation_var, 
    cca_weights = NULL
    ) {

    logmsg("\n[1/4] Finding HVGs, normalizing, and standardizing...")
    logmsg("  Counts dimensions:", nrow(sc), "genes x", ncol(sc), "cells")
    
    sc <- sc %>% 
        NormalizeData(normalization.method = 'LogNormalize', scale.factor = median(sc$nCount_RNA)) %>% 
        ScaleData(features = hvgs) 
    
    logmsg("  Number of highly variable genes: ", length(hvgs))
    
    logmsg("\n[2/4] Running dimensional reduction...")
    
    set.seed(0)
    sc <- RunPCA(sc, features = hvgs, verbose = TRUE)
    reduction_for_harmony <- "pca"
    
    logmsg("  Printing PCA diagnostic plots...")
    #TODO fix plotting issue with last plot
    print(DimPlot(sc, reduction = reduction_for_harmony, group.by = "sid") + NoLegend())
    #print(DimPlot(merged, reduction = reduction_for_harmony, group.by = annotation_var) + NoLegend())
    Sys.sleep(0.5)

    logmsg("\n[3/4] Running Harmony...")
    set.seed(0)
    sc <- harmony::RunHarmony(
        object = sc,
        reduction.use = reduction_for_harmony,
        group.by.vars = c("sid", batch_vars),
        plot_convergence = TRUE,
        max_iter = 10,
        early_stop = FALSE
    )

    logmsg("  Printing Harmony diagnostic plots...")
    print(DimPlot(sc, reduction = "harmony", group.by = "sid") + NoLegend())
    #print(DimPlot(merged, reduction = "harmony", group.by = annotation_var) + NoLegend())
    Sys.sleep(0.5)

    logmsg("\n[4/4] Running UMAP on the Harmony reduction...")
    set.seed(0)
    sc <- Run_uwot_umap(
        sc,
        reduction = "harmony",
        spread = 0.8,
        min_dist = 0.3
    )
    return(sc)
                      }

BuildIntegratedReference <- function(
    xen,
    sc,
    normalization_target = c("mean_of_medians", "sc_median"),
    hvgs,
    annotation_var,
    batch_vars = NULL, 
    cca_weights = NULL
) {
    normalization_target <- match.arg(normalization_target)

    logmsg("\n[1/7] Retrieving count matrices and checking gene overlap...")
    xen_counts <- xen@assays$RNA$counts
    sc_counts  <- sc@assays$RNA$counts

    logmsg("  Xenium counts dimensions:", nrow(xen_counts), "genes x", ncol(xen_counts), "cells")
    logmsg("  Single-cell counts dimensions:", nrow(sc_counts), "genes x", ncol(sc_counts), "cells")

    common_genes <- intersect(rownames(xen_counts), rownames(sc_counts))
    logmsg("  Number of genes shared between objects:", length(common_genes))

    if (!annotation_var %in% colnames(sc@meta.data)) {
        stop("annotation_var not found in sc@meta.data: ", annotation_var)
    }
    if (!"sid" %in% colnames(xen@meta.data) || !"sid" %in% colnames(sc@meta.data)) {
        stop("Both xen and sc must contain a 'sid' column in meta.data.")
    }

    logmsg("\n[2/7] Preparing Seurat objects for merging...")
    xen$modality <- "xen"
    sc$modality  <- "sc"

    xen[[annotation_var]] <- NA_character_
    
    if (!is.null(batch_vars)) {
        keep_cols <- c("sid", "modality", batch_vars, annotation_var)
    } else {
        keep_cols <- c("sid", "modality", annotation_var)
    }    
    xen_meta <- xen@meta.data[, keep_cols, drop = FALSE]
    sc_meta  <- sc@meta.data[, keep_cols, drop = FALSE]

    sc_obj <- CreateSeuratObject(
        counts = sc_counts[common_genes, , drop = FALSE],
        meta.data = sc_meta,
        project = "sc"
    )
    xen_obj <- CreateSeuratObject(
        counts = xen_counts[common_genes, , drop = FALSE],
        meta.data = xen_meta,
        project = "xen"
    )

    merged <- merge(xen_obj, sc_obj, add.cell.ids = c("xen", "sc"))
    merged <- JoinLayers(merged)

    logmsg("  Merged object created.")
    logmsg("  Merged object summary:")
    print(merged)
    logmsg("  Merged object dimensions:", nrow(merged), "genes x", ncol(merged), "cells")
    logmsg("  Cells by modality:")
    print(table(merged$modality))

    logmsg("\n[3/7] Choosing normalization scale factor...")
    xen_median_ncount <- median(merged$nCount_RNA[merged$modality == "xen"])
    sc_median_ncount  <- median(merged$nCount_RNA[merged$modality == "sc"])
    if (normalization_target == "mean_of_medians") {
        scale_factor <- mean(c(xen_median_ncount, sc_median_ncount))
        logmsg("  Using mean of modality-specific medians.")
    } else if (normalization_target == "sc_median") {
        scale_factor <- sc_median_ncount
        logmsg("  Using median nCount_RNA from sc object only.")
    }

    logmsg("  Xen median nCount_RNA:", xen_median_ncount)
    logmsg("  SC median nCount_RNA:", sc_median_ncount)
    logmsg("  Selected normalization scale factor:", scale_factor)

    logmsg("\n[4/7] Normalizing data, setting variable features, and scaling...")
    merged <- NormalizeData(
        merged,
        normalization.method = "LogNormalize",
        scale.factor = scale_factor,
        verbose = FALSE
    )

    hvgs_use <- intersect(hvgs, rownames(merged))
    logmsg("  HVGs supplied:", length(hvgs))
    logmsg("  HVGs present in merged object:", length(hvgs_use))
    VariableFeatures(merged) <- hvgs_use
    merged <- ScaleData(merged, features = hvgs_use, split.by = "modality")
    
    logmsg("\n[5/7] Running dimensional reduction...")
    set.seed(0)

    if (is.null(cca_weights) || (length(cca_weights) == 1 && is.na(cca_weights))) {
        logmsg("  cca_weights is NULL/NA; running PCA.")
        merged <- RunPCA(merged, features = hvgs_use, verbose = TRUE)
        reduction_for_harmony <- "pca"

    } else {
        logmsg("  cca_weights provided; constructing CCA-based reduction.")

        if (is.null(rownames(cca_weights))) {
            stop("cca_weights must have gene names as rownames.")
        }
        hvgs_in_scale <- intersect(hvgs_use, rownames(merged@assays$RNA$scale.data))
        hvgs_in_weights <- intersect(hvgs_use, rownames(cca_weights))
        hvgs_final <- intersect(hvgs_in_scale, hvgs_in_weights)
        logmsg("  HVGs in scale.data:", length(hvgs_in_scale))
        logmsg("  HVGs in cca_weights:", length(hvgs_in_weights))
        logmsg("  HVGs used for CCA projection:", length(hvgs_final))

        cca_embedding <- t(merged@assays$RNA$scale.data[hvgs_final, , drop = FALSE]) %*%
            cca_weights[hvgs_final, , drop = FALSE]

        colnames(cca_embedding) <- paste0("CV_", seq_len(ncol(cca_embedding)))
        merged[["cca"]] <- Seurat::CreateDimReducObject(
            embeddings = cca_embedding,
            assay = "RNA",
            key = "CV_"
        )
        reduction_for_harmony <- "cca"
    }
    
    logmsg("  Printing PCA diagnostic plots...")
    #TODO fix plotting issue with last plot
    print(DimPlot(merged, reduction = reduction_for_harmony, group.by = "modality") + NoLegend())
    print(DimPlot(merged, reduction = reduction_for_harmony, group.by = "sid") + NoLegend())
    #print(DimPlot(merged, reduction = reduction_for_harmony, group.by = annotation_var) + NoLegend())
    Sys.sleep(0.5)

    logmsg("\n[6/7] Running Harmony...")
    set.seed(0)
    merged <- harmony::RunHarmony(
        object = merged,
        reduction.use = reduction_for_harmony,
        group.by.vars = c("modality", "sid", batch_vars),
        plot_convergence = TRUE,
        max_iter = 10,
        early_stop = FALSE
    )

    logmsg("  Printing Harmony diagnostic plots...")
    print(DimPlot(merged, reduction = "harmony", group.by = "modality") + NoLegend())
    print(DimPlot(merged, reduction = "harmony", group.by = "sid") + NoLegend())
    #print(DimPlot(merged, reduction = "harmony", group.by = annotation_var) + NoLegend())
    Sys.sleep(0.5)

    logmsg("\n[7/7] Running UMAP on the Harmony reduction...")
    set.seed(0)
    merged <- Run_uwot_umap(
        merged,
        reduction = "harmony",
        spread = 0.8,
        min_dist = 0.3
    )
    return(merged)
}

ShowTransferDiagnostics <- function(
    merged,
    annotation_var,
    genes_to_visualize = c("KLRD1", "GNLY", "KLRF1"), 
    save_plot = NULL # this should be a file path if not null  
) {
    p1 <- DimPlot(merged, reduction = "humap", group.by = "modality", raster = TRUE, shuffle = TRUE) + NoAxes() + small_theme
    p2 <- DimPlot(
        merged[,merged$modality == "sc"],
        reduction = "humap",
        group.by = annotation_var,
        raster = TRUE,
        cols = Seurat::DiscretePalette(length(unique(merged@meta.data[[annotation_var]]))),
        shuffle = TRUE
    ) + NoAxes() + small_theme
    p12 <- p1 + p2
    f12 <- tempfile(fileext = ".png")
    ggplot2::ggsave(f12, p12, width = 12, height = 5, dpi = 150)
    IRdisplay::display_png(file = f12)

    genes_present <- intersect(genes_to_visualize, rownames(merged))
    genes_missing <- setdiff(genes_to_visualize, rownames(merged))
    if (length(genes_missing) > 0) {
        logmsg("  Genes not found and skipped:", paste(genes_missing, collapse = ", "))
    }

   if (length(genes_present) > 0) {
        p3 <- FeaturePlot(
            merged[,merged$modality == "sc"],
            features = genes_present,
            cols = scico::scico(100, palette = "batlow", direction = -1),
            ncol = min(3, length(genes_present)),
            raster = TRUE,
            order = TRUE
        ) &
        coord_fixed()
        f3 <- tempfile(fileext = ".png")
        ggplot2::ggsave(f3, p3, width = 20, height = 10, dpi = 150)
        IRdisplay::display_png(file = f3)
    } else {
        logmsg("  No requested genes were found, so skipping FeaturePlot.")
    }
    
    if (!(is.null(save_plot))) {
        if (length(genes_present) > 0) {
            p123 <- (p1 + p2) / p3
            ggplot2::ggsave(save_plot, p123, width = 20, height = 20, dpi = 150)
        } else {
            ggplot2::ggsave(save_plot, p12, width = 20, height = 10, dpi = 150)
        }
    
    }

    invisible(NULL)
}




TransferLabelsGraph <- function(
    merged,
    annotation_var,
    graph = NULL,
    n_neighbors = 30,
    min_ratio = 1.1,
    raster = FALSE 
) {
    logmsg("\n[1/4] Checking inputs...")
    if (!"harmony" %in% names(merged@reductions)) {
        stop("The merged object does not contain a 'harmony' reduction.")
    }
    if (!annotation_var %in% colnames(merged@meta.data)) {
        stop("annotation_var not found in merged@meta.data: ", annotation_var)
    }
    if (!all(c("xen", "sc") %in% unique(merged$modality))) {
        stop("merged$modality must contain both 'xen' and 'sc'.")
    }
    logmsg("  Using annotation variable:", annotation_var)
    
    if(is.null(graph)) {
        logmsg("  Using n_neighbors =", n_neighbors, "and min_ratio =", min_ratio)
        logmsg("\n[2/4] Computing nearest-neighbor graph from Harmony embedding...")
        set.seed(0)
        umap_with_graph <- uwot::umap(
            merged@reductions$harmony@cell.embeddings,
            n_neighbors = n_neighbors,
            spread = 0.8,
            min_dist = 0.3,
            ret_extra = "fgraph",
            fast_sgd = FALSE
        )
        logmsg("  Graph computation complete.")
        graph <- umap_with_graph$fgraph
    } else {
        logmsg("\n[2/4] Skipping NNG construction since graph is provided...")
    }
    
    logmsg("\n[3/4] Transferring labels from sc cells to xen cells...")
    query_mask <- merged$modality == "xen"
    ref_mask   <- merged$modality == "sc"
    transferred_labels <- transfer_labels_graph(
        graph = graph,
        query_mask = query_mask,
        ref_mask = ref_mask,
        labels = merged[[annotation_var]][, 1],
        min_ratio = min_ratio
    )
    merged@meta.data[query_mask, annotation_var] <- transferred_labels
    logmsg("  Labels written to:", annotation_var)
    logmsg("  Xen cells without transferred labels:", sum(is.na(transferred_labels)))

    logmsg("\n[4/4] Plotting unlabeled cells and final xen label map...")
    cells_na <- colnames(merged)[
        merged$modality == "xen" & is.na(merged[[annotation_var]][, 1])
    ]
    merged@meta.data[colnames(merged) %in% cells_na, annotation_var] <- 'Unknown'
    p1 <-
        DimPlot(
            merged,
            reduction = "humap",
            cells.highlight = list("Unlabelled" = cells_na),
            cols.highlight = scales::alpha("red", 0.01),
            cols = "grey85",
            shuffle = TRUE, 
            raster = raster
        ) + NoAxes() + small_theme
    p2 <-
        DimPlot(
            subset(merged, subset = modality == "xen"),
            reduction = "humap",
            group.by = annotation_var,
            na.value = "grey80", 
            raster = FALSE
        ) + NoAxes() + small_theme
    p12 <- p1 + p2
    f12 <- tempfile(fileext = ".png")
    ggplot2::ggsave(f12, p12, width = 12, height = 5, dpi = 150)
    IRdisplay::display_png(file = f12)

    return(merged)
}

ShowClusterMarkers <- function(xen, celltype_var, n_show = 30) {
    expr <- GetAssayData(xen, layer = "data")
    clusters.markers <- presto::wilcoxauc(expr, xen[[celltype_var]][,1])

    options(repr.matrix.max.cols = 150, repr.matrix.max.rows = 200)
    show <- clusters.markers %>%
        group_by(group) %>%
        filter(padj < 0.05 & auc > 0.5) %>% 
        arrange(group, desc(logFC)) %>% 
        dplyr::select(feature, group) %>% 
        mutate(row = row_number()) %>% 
        pivot_wider(names_from = group, values_from = feature)

    return(show[1:n_show, ])
}

CompareMarkerCorrelations <- function(xen, sc, celltype_var, bicluster = FALSE, markers.xen = NULL, markers.sc = NULL,
                                      hvgs_only = FALSE, pheatmap_breaks = NULL, fontsize = 12, show = TRUE) {
  # make markers
  hvgs <- rownames(xen[['RNA']]$scale.data)
  if(is.null(markers.xen)) { 
    expr <- GetAssayData(xen, layer = "data")
    if (hvgs_only) {
        expr <- expr[hvgs, ]
    }
    markers.xen <- presto::wilcoxauc(expr, xen[[celltype_var]][, 1])
  }
  if(is.null(markers.sc)) { 
    expr <- GetAssayData(sc, layer = "data")
    if (hvgs_only) {
        expr <- expr[hvgs, ]
    }
    markers.sc <- presto::wilcoxauc(expr, sc[[celltype_var]][, 1])
  }

  # helper: convert marker table to gene x cluster matrix
  marker_table_to_mat <- function(markers, group_col = "group", feature_col = "feature", value_col = "logFC") {
    markers %>%
      dplyr::select(
        group = all_of(group_col),
        feature = all_of(feature_col),
        value = all_of(value_col)
      ) %>%
      distinct(group, feature, .keep_all = TRUE) %>%
      tidyr::pivot_wider(
        names_from = group,
        values_from = value,
        values_fill = 0
      ) %>%
      as.data.frame() -> df

    rownames(df) <- df$feature
    df$feature <- NULL
    as.matrix(df)
  }

  # compare sc to xenium
  xen_mat <- marker_table_to_mat(markers.xen, value_col = "logFC")
  sc_mat  <- marker_table_to_mat(markers.sc,  value_col = "logFC")
  common_genes <- intersect(rownames(xen_mat), rownames(sc_mat))
  xen_mat <- xen_mat[common_genes, , drop = FALSE]
  sc_mat  <- sc_mat[common_genes, , drop = FALSE]

  p1 <- pheatmap::pheatmap(
    cor(sc_mat, sc_mat, method = "pearson", use = "pairwise.complete.obs"),
    cluster_rows = bicluster,
    cluster_cols = bicluster,
    treeheight_row = 0,
    treeheight_col = 0,
    fontsize_row = fontsize,
    fontsize_col = fontsize,
    border_color = NA,
    main = "sc vs sc",
    breaks = pheatmap_breaks, 
    silent = TRUE
  )
  # Extract ordering from p1
  mat2 <- cor(xen_mat, sc_mat, method = "pearson", use = "pairwise.complete.obs")
  if (bicluster) {
    row_order <- p1$tree_row$order
    col_order <- p1$tree_col$order
  } else {
    row_order <- seq_len(nrow(mat2))
    col_order <- seq_len(ncol(mat2))
  }
  mat2 <- mat2[row_order, col_order, drop = FALSE]
  p2 <- pheatmap::pheatmap(
    mat2,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    fontsize_row = fontsize,
    fontsize_col = fontsize,
    border_color = NA,
    main = "xen vs sc",
    breaks = pheatmap_breaks, 
    display_numbers = TRUE, 
    fontsize_number = fontsize, 
    silent = TRUE
  )

  if (show) {
    fig.size(10, 20)
    gridExtra::grid.arrange(p1[[4]], p2[[4]], ncol = 2)
  }
  
  invisible(list(
    markers.xen = markers.xen,
    markers.sc = markers.sc,
    xen_mat = xen_mat,
    sc_mat = sc_mat,
    common_genes = common_genes,
    p1 = p1,
    p2 = p2
  ))
}

# Usage: after clustering Xenium-only embedding at high resolution, collapse Leiden clusters into AMPp2 cell types/states. If marker genes for a cell type are provided, exclude clusters with poor expression of those genes. 

AssignClusterLabels <- function(
    xen,
    sc, 
    leiden_res,
    min_corr, 
    ref_col, 
    hvgs_only = FALSE, 
    ct_marker_genes = NULL, 
    min_ct_marker_score = 0, 
    bicluster = FALSE, 
    markers.xen = NULL, 
    markers.sc = NULL,                                  
    pheatmap_breaks = NULL, 
    fontsize = 12
    ) {
    
    leiden_colname = paste0('humap_fgraph_res.', leiden_res)
    xen@meta.data[[ref_col]] = xen@meta.data[[leiden_colname]]
    pre_collapse <- CompareMarkerCorrelations(xen, sc, ref_col, bicluster, markers.xen, markers.sc, hvgs_only, 
                              pheatmap_breaks, fontsize, show = TRUE)
    cor_mat <- cor(pre_collapse$xen_mat, pre_collapse$sc_mat, method = "pearson", use = "pairwise.complete.obs")
    
    if (!is.null(ct_marker_genes)) {
        ct_score = colSums(xen[['RNA']]$data[ct_marker_genes, ]) < min_ct_marker_score
        other_ct_clusters = as.data.frame(table(ct_score, xen@meta.data[[ref_col]])) %>% 
            rename('cluster' = 'Var2') %>% 
            group_by(cluster) %>% 
            pivot_wider(names_from = ct_score, values_from = Freq, values_fill = 0) %>% 
            filter(`TRUE` > 0.2 * `FALSE`) %>% 
            pull(cluster)
    } else {
        other_ct_clusters = c()
    }
    
    cluster_labels = cor_mat %>% 
        as.data.frame() %>% 
        rownames_to_column('cluster_id') %>% 
        pivot_longer(cols = -cluster_id, names_to = 'cell_state', values_to = 'corr') %>% 
        group_by(cluster_id) %>% 
        filter(corr == max(corr)) %>% 
        ungroup() %>% 
        mutate(cell_state = case_when(
            !(is.null(ct_marker_genes)) & (cluster_id %in% other_ct_clusters) ~ 'Other cell type', 
            corr < min_corr ~ 'Correlation too low', 
            .default = cell_state)) %>% 
        dplyr::select(-corr) 
    
    clusters_failing_correlation = which(cluster_labels$cell_state == 'Correlation too low')
    print(paste(paste(clusters_failing_correlation, collapse = ","), "all had low correlations"))

    xen$seurat_clusters <- xen@meta.data[[leiden_colname]]
    
    if (!is.null(ct_marker_genes)) {
        print(paste(length(other_ct_clusters), 'clusters had cells of other cell types'))
        print(paste('These clusters contain', sum(xen$seurat_clusters %in% other_ct_clusters), 'cells'))
    }
    print(paste(length(clusters_failing_correlation), 'clusters had low correlation'))
    print(paste('These clusters contain', sum(xen$seurat_clusters %in% clusters_failing_correlation), 'cells'))
                                                             
    new_cluster_name = xen[[]] %>% 
        left_join(cluster_labels, by = join_by('seurat_clusters' == 'cluster_id')) %>% 
        mutate(cluster_name_plot = str_split_i(cell_state, ":", 1)) %>% 
        pull(cluster_name_plot) 

    xen$cluster_name_plot <- new_cluster_name 
    
    sc$cluster_name_plot <- sc[[]] %>% 
        mutate(cluster_name_plot = str_split_i(.data[[ref_col]], ":", 1)) %>% 
        pull(cluster_name_plot)
    
    res <- CompareMarkerCorrelations(xen, sc, 'cluster_name_plot', bicluster, markers.xen = NULL, markers.sc = NULL, hvgs_only, 
                              pheatmap_breaks, fontsize, show = TRUE)
    
    invisible(list(
        xen_cluster_name_plot = xen$cluster_name_plot, 
        markers.xen = res$markers.xen,
        markers.sc = res$markers.sc,
        xen_mat = res$xen_mat,
        sc_mat = res$sc_mat,
        common_genes = res$common_genes,
        p1 = res$p1,
        p2 = res$p2
        ))
}
    

PlotStateLogFCComparison <- function(
  markers_sc,
  markers_xen,
  cell_state,
  base_size = 20,
  label_size = 4,
  n_label = 2, 
  fig_size_x = NULL, 
  fig_size_y = NULL 
) {
  library(dplyr)
  library(ggplot2)
  library(ggrepel)

  df <- markers_sc %>%
    dplyr::filter(group == cell_state) %>%
    dplyr::select(feature, logFC_sc = logFC) %>%
    dplyr::inner_join(
      markers_xen %>%
        dplyr::filter(group == cell_state) %>%
        dplyr::select(feature, logFC_xen = logFC),
      by = "feature"
    )

  label_df <- df %>%
    dplyr::filter(
      feature %in% c(
        dplyr::slice_max(df, logFC_sc, n = n_label)$feature,
        dplyr::slice_min(df, logFC_sc, n = n_label)$feature,
        dplyr::slice_max(df, logFC_xen, n = n_label)$feature,
        dplyr::slice_min(df, logFC_xen, n = n_label)$feature
      )
    ) %>%
    dplyr::distinct(feature, .keep_all = TRUE)
  if (!(is.null(fig_size_x) | is.null(fig_size_y))) {
      fig.size(fig_size_x, fig_size_y)
  }
  p <- ggplot(df, aes(x = logFC_xen, y = logFC_sc)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = feature),
      size = label_size,
      max.overlaps = Inf
    ) +
    labs(
      x = "logFC (Xenium)",
      y = "logFC (scRNA-seq)",
      title = paste("logFC comparison:", cell_state)
    ) +
    theme_classic(base_size = base_size)

  print(p)

  r <- cor(df$logFC_xen, df$logFC_sc, use = "pairwise.complete.obs")

  invisible(list(
    plot = p,
    data = df,
    labeled_genes = label_df,
    correlation = r
  ))
}