suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    source('/data/srlab/AMP_collab/lakshay-yakir/_common/typing_utils.r')
})

args <- commandArgs(trailingOnly=TRUE)
ct <- args[1]
ct <- gsub("_", " ", ct)
cohort <- args[2]
max_mult <- as.numeric(args[3])

# Step 1: read in Xenium and AMPp2 objects 

ampp2_ref <- '/data/srlab/AMP_collab/AMP2_2023_seuratObj_galoz/AMPp2_seuratObj_allcells.202511.rds'
sc <- readRDS(ampp2_ref)
names(sc@meta.data)[1] <- 'sid'
ct_to_keep = sc$celltypes.med == ct 
sc <- subset(sc, cells = colnames(sc)[ct_to_keep])    
print(paste('Number of sc cells:', dim(sc)[2]))

xen_path <- file.path('/data/srlab/AMP_collab/data/early_disease_synovium/xenium/combined/', 
                      cohort, 
                      'coarsetypes', 
                      paste0(ct, '.rds'))
if (!file.exists(xen_path)) {
  logmsg(paste0("ERROR: path not found: ", xen_path))
  quit(save = "no", status = 1)
}

xen <- readRDS(xen_path)
# Step 2: break Xenium lineage down into smaller chunks, and save each chunk 

i <- 1
sampled_tracker = setNames(character(dim(xen)[2]), colnames(xen))
while (sum(sampled_tracker == "", na.rm = TRUE) > 0) {
    unsampled_cells = which(sampled_tracker == "")
    num_iters = ceiling(length(unsampled_cells) / (max_mult * dim(sc)[2]))

    if (num_iters > 1) {
        sampled_cells = sample(unsampled_cells, size = max_mult * dim(sc)[2])
    } else {
        sampled_cells = unsampled_cells
    }

    logmsg(paste('There are', length(unsampled_cells), 'remaining Xenium cells in this lineage'))
    logmsg(paste('This will require', num_iters, 'more chunks to be saved'))

    xen_sub = subset(xen, cells = colnames(xen)[sampled_cells]) 
    out_dir = file.path('/data/srlab/AMP_collab/lakshay-yakir/6.fine_types/out_rds/', 
                         cohort, 
                         ct)
    dir.create(out_dir, showWarnings = FALSE) 
    saveRDS(xen_sub, file.path(out_dir,  
                               paste0(ct, "_chunk", i, "_untyped.rds"))) 

    sampled_tracker[sampled_cells] = 'Sampled'
    logmsg(paste('Chunk', i, 'saved to',  file.path(out_dir,  
                               paste0(ct, "_chunk", i, "_untyped.rds"))))
    i <- i + 1 
}




