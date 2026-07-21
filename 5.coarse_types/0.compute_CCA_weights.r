suppressPackageStartupMessages({
    library(CCA)
    library(tidyverse)
    library(Seurat)
    library(Matrix)
    source('/data/srlab/AMP_collab/lakshay-yakir/_common/typing_utils.r')
})

# Step 1: read in data and subset to cell type of interest 
args = commandArgs(trailingOnly = TRUE)
ct <- args

ref_path <- "/data/srlab/AMP_collab/AMP2_2023_seuratObj_galoz/AMPp2_seuratObj_allcells.202511.rds"
ref <- readRDS(ref_path)
ref <- subset(ref, cells = colnames(ref)[ref$celltypes.med %in% ct])

# Step 2: Restrict to Xenium genes and re-QC

xenium_genes <- readLines('/data/srlab/AMP_collab/lakshay/20260603_5100_xenium_genes.txt')
ref <- subset(ref, features = intersect(rownames(ref), xenium_genes))
ngenes_qc <- 100
ncounts_qc <- 150
qc_mask <- (colSums(ref[['RNA']]$counts > 0) > ngenes_qc) & (colSums(ref[['RNA']]$counts) > ncounts_qc)
ref <- subset(ref, cells = colnames(ref)[qc_mask])

adt_exprs_norm <- readRDS("/data/srlab/fzhang/amp/results/2020_01_22_AMP_cite_seq_QC/adt_exprs_norm_filter_2020-02-18.rds") # adt_exprs_norm are CLR normalized and QCed
adt_exprs_norm <- adt_exprs_norm[, colnames(ref)]
ref[['ADT']] <- CreateAssayObject(counts = adt_exprs_norm)
print(paste("Dimensions of gexp are", dim(ref)))
print(paste("Dimensions of ADTs are", dim(adt_exprs_norm)))

# Step 3: Compute HVGs and HVPs 
              
DefaultAssay(ref) <- 'RNA'              
ref$nCount_RNA  <- colSums(ref[['RNA']]$counts)
ref$nFeature_RNA <- colSums(ref[['RNA']]$counts > 0)
ref <- ref %>% 
    NormalizeData(normalization.method = 'LogNormalize', scale.factor = median(ref$nCount_RNA))
ref$sid <- ref$orig.ident 

res <- FindHVGsFromGroups(ref, group_var = "sid", logFC_thresh = 0.1)
hvgs <- res$hvgs
info <- res$info
              
ref <- ref %>% ScaleData(features = hvgs) 
              
compute_kl <- function(prot) {
    high_exp_mask <- cume_dist(adt_exprs_norm[match(prot, rownames(adt_exprs_norm)), ]) > 0.85
    high_exp_dist <- table(ref@meta.data[high_exp_mask, 'cluster_name'])
    high_exp_dist <- high_exp_dist/sum(high_exp_dist)
    full_dist <- table(ref$cluster_name)
    full_dist <- full_dist/sum(full_dist)    
    return(sum(high_exp_dist * log2(high_exp_dist/full_dist + 1e-9)))
}

adt_final <- c()
kls <- c()
for (prot in rownames(adt_exprs_norm)) {
    kl <- compute_kl(prot) 
    kls <- c(kls, kl)
    if (kl > 0.025) {
        adt_final <- c(adt_final, prot)
    }
}
adt_final <- gsub("_", "-", adt_final)
DefaultAssay(ref) <- 'ADT'
ref <- ScaleData(ref, features = adt_final)
              
hvgs <- ref[['RNA']]$scale.data 
hvps <- ref[['ADT']]$scale.data 
print(paste("Dimensions of gexp are", dim(hvgs)))
print(paste("Dimensions of ADTs are", dim(hvps)))
              
# Step 4: Compute CCA 
              
set.seed(0)
system.time({
    res_cca = cc(t(hvgs), t(hvps))
    })
        
outdir <- '/data/srlab/AMP_collab/lakshay-yakir/5.coarse_types/out_ccaweights/'
dir.create(outdir, showWarnings = FALSE)
outpath <- paste0(outdir, paste0(ct, collapse = ""), '_ccaweights.RDS')
saveRDS(res_cca, outpath) 
              
              


