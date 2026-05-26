#!/usr/bin/env Rscript
# combine_cells.R
# Combines Baysor segmentation outputs into a Seurat object with spatial metadata.
# Equivalent to notebook: 3a_combine_cells.ipynb

suppressPackageStartupMessages({
    library(data.table)
    library(Matrix)
    library(Seurat)
    library(ggplot2)
    library(dplyr)
    library(ggthemes)
    library(tidyr)
    library(tibble)
    library(sf)
    library(sfarrow)
})

# ==============================================================================
# 1. Read and concatenate metadata
# ==============================================================================

meta_files <- list(
  "/data/srlab/AMP_collab/data/early_disease_synovium/Xenium_CTRL-SYN-EDP1_assay manifest.csv",
  "/data/srlab/AMP_collab/data/early_disease_synovium/Xenium_RA-SYN-EDP1_assay manifest.csv",
  "/data/srlab/AMP_collab/data/early_disease_synovium/Xenium_RA-SYN-EDP2_assay manifest.csv",
  "/data/srlab/AMP_collab/data/early_disease_synovium/Xenium_PsD-SYN-EDP1_assay manifest.csv",
  "/data/srlab/AMP_collab/data/early_disease_synovium/Xenium_RA-SYN-ARBITRATE_assay manifest.csv"
)

get_subset_and_edp <- function(path) {
  fname <- basename(path)
  core <- sub("^Xenium_(.*)_assay manifest\\.csv$", "\\1", fname)
  if (core == "RA-SYN-ARBITRATE") {
      list(subset = "RA-SYN", cohort = "ARBITRATE")
  } else {
      edp <- sub(".*-(EDP[12])$", "\\1", core)
      subset <- sub(paste0("-", edp, "$"), "", core)
      list(subset = subset, cohort = edp)
  }
}

meta_list <- lapply(meta_files, function(f) {
  dt <- fread(f)
  info <- get_subset_and_edp(f)
  dt[, subset := info$subset]
  dt[, cohort := info$cohort]
  dt
})

metadata <- rbindlist(meta_list, use.names = TRUE, fill = TRUE)

# ==============================================================================
# 2. Helper: column-bind sparse matrices with mismatched row sets
# ==============================================================================

cbind2_fill <- function(mat_list) {
    rownames_all <- unique(unlist(lapply(mat_list, rownames)))

    add_list <- list()
    for (i in seq_along(mat_list)) {
        cat("."); flush.console()
        mat <- mat_list[[i]]
        add_list[[i]] <- setdiff(rownames_all, rownames(mat))
        toadd <- Matrix::rsparsematrix(length(add_list[[i]]), ncol(mat), 0)
        mat_list[[i]] <- Matrix::rbind2(mat, toadd)
        rownames(mat_list[[i]]) <- c(rownames(mat), add_list[[i]])
        mat_list[[i]] <- mat_list[[i]][rownames_all, , drop = FALSE]
    }

    return(Reduce(Matrix::cbind2, mat_list))
}

# ==============================================================================
# 3. Build sparse count matrix and create Seurat object
# ==============================================================================

t0 <- Sys.time()
sids <- unique(metadata$Sample_ID)
fov_list <- list()

for (sid in sids) {
  cohort <- metadata[metadata$Sample_ID == sid]$cohort[[1]]
  cat("\nProcessing sample:", sid, cohort, "\n")
  flush.console()

  pattern <- file.path("out", sid, "baysor_out", "*", "segmentation_counts.tsv")
  count_files <- Sys.glob(pattern)

  if (length(count_files) == 0) {
    warning(sprintf("No Baysor output found for sid: %s", sid))
    next
  }

  for (i in seq_along(count_files)) {
    cf <- count_files[[i]]
    cat(sprintf("\tFOV %d/%d: %s ", i, length(count_files), cf))
    flush.console()

    df <- read.csv(cf, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    gene_names <- df[[1]]
    cell_ids <- colnames(df)[-1]
    mat_sparse <- Matrix(as.matrix(df[, -1]), sparse = TRUE)

    rownames(mat_sparse) <- gene_names
    colnames(mat_sparse) <- paste0(cell_ids, "_", sid)
    fov_list[[cf]] <- mat_sparse

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("(%.1f sec)\n", elapsed))
    flush.console()
  }
}

cat("Concatenating FOV matrices...\n")
flush.console()
big <- cbind2_fill(fov_list)
big <- Seurat::CreateSeuratObject(
  counts = big,
  project = cohort
)
big[["sid"]] <- sub(".*_", "", rownames(big[[]]))

dir.create("out_rds", showWarnings = FALSE)
saveRDS(big, file = "out_rds/allcells.rds")
cat("Seurat object saved to out_rds/allcells.rds\n")

# ==============================================================================
# 4. Add spatial cell metadata (x, y, confidence, polygon WKT) to Seurat object
# ==============================================================================

t0 <- Sys.time()

for (sid in sids) {
  cohort <- metadata[metadata$Sample_ID == sid, ]$cohort[[1]]
  cat("\nProcessing sample:", sid, cohort, "\n")
  flush.console()

  pattern <- file.path("out", sid, "baysor_out", "*", "segmentation_cell_stats.csv")
  count_files <- Sys.glob(pattern)

  if (length(count_files) == 0) {
    warning(sprintf("No Baysor output found for sid: %s", sid))
    next
  }

  for (i in seq_along(count_files)) {
    cf <- count_files[[i]]
    cat(sprintf("\tFOV %d/%d: %s ", i, length(count_files), cf))
    flush.console()

    df <- read.csv(cf, stringsAsFactors = FALSE, check.names = FALSE)

    poly_file <- file.path(dirname(cf), "segmentation_polygons_2d.json")
    if (!file.exists(poly_file)) {
      warning(sprintf("Missing polygon file for %s", cf))
      next
    }

    geom <- st_read(poly_file, quiet = TRUE)
    geom_df <- data.frame(
      cell = geom$id,
      polygon_wkt = st_as_text(st_geometry(geom)),
      stringsAsFactors = FALSE
    )

    cellmeta <- merge(
      df[, c("cell", "x", "y", "avg_assignment_confidence")],
      geom_df,
      by = "cell",
      all.x = TRUE
    )
    cellmeta$cell_full <- paste0(cellmeta$cell, "_", sid)

    big@meta.data[cellmeta$cell_full,
                  c("x", "y", "avg_assignment_confidence", "polygon_wkt")] <-
      cellmeta[, c("x", "y", "avg_assignment_confidence", "polygon_wkt")]

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("(%.1f sec)\n", elapsed))
    flush.console()
  }
}

cat("Writing final output...\n")
saveRDS(big, file = "out_rds/allcells.rds")
cat("Done. Final object saved to out_rds/allcells.rds\n")