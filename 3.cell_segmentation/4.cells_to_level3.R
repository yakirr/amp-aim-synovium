#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(rhdf5)
    library(arrow)
})

# --- argument parsing ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript 4.cells_to_level3.R <i>  (1-based index into sids)")
}
i <- as.integer(args[1])

# --- load data ---
message("Reading RDS file")
allcells <- readRDS('out_rds/allcells.rds')

counts_mat <- Seurat::GetAssayData(allcells, layer = "counts")
meta       <- allcells@meta.data
sids       <- unique(meta$sid)

if (i < 1 || i > length(sids)) {
    stop(sprintf("Index i=%d is out of range (1-%d)", i, length(sids)))
}

sid <- sids[i]
message("Processing sid #", i, ": ", sid)

# --- output directory ---
dir.create("level3/", showWarnings = FALSE, recursive = TRUE)

# --- subset to this sample ---
cells_use  <- rownames(meta)[meta$sid == sid]
sub_counts <- counts_mat[, cells_use, drop = FALSE]
sub_meta   <- meta[cells_use, , drop = FALSE]

# --- output file paths ---
h5_file <- file.path("level3", paste0(sid, ".h5"))
pq_file <- file.path("level3", paste0(sid, ".parquet"))

if (file.exists(h5_file)) file.remove(h5_file)
if (file.exists(pq_file)) file.remove(pq_file)

# --- write HDF5 ---
rhdf5::h5createFile(h5_file)
rhdf5::h5write(as.matrix(sub_counts), h5_file, "counts")
rhdf5::h5write(colnames(sub_counts),  h5_file, "cells")
rhdf5::h5write(rownames(sub_counts),  h5_file, "genes")

# --- write parquet ---
sub_meta_out <- cbind(cell_id = rownames(sub_meta), sub_meta)
arrow::write_parquet(sub_meta_out, pq_file)

message("Finished sid #", i, ": ", sid)