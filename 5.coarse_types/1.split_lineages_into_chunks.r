suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    source('/data/srlab/AMP_collab/lakshay-yakir/_common/typing_utils.r')
})

args <- commandArgs(trailingOnly=TRUE)
lineage <- args[1]
cohort <- args[2]
max_mult <- as.numeric(args[3])

# Step 1: read in Xenium and AMPp2 objects 

ampp2_ref <- '/data/srlab/AMP_collab/AMP2_2023_seuratObj_galoz/AMPp2_seuratObj_allcells.202511.rds'
sc <- readRDS(ampp2_ref)
names(sc@meta.data)[1] <- 'sid'
lineage_to_keep = sc$lineage == lineage 
sc <- subset(sc, cells = colnames(sc)[lineage_to_keep])    
print(paste('Number of sc cells:', dim(sc)[2]))

xen_path <- file.path('/data/srlab/AMP_collab/data/early_disease_synovium/xenium/combined/', 
                      cohort, 
                      'lineages', 
                      paste0(lineage, '.rds'))
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
    dir.create(file.path(dirname(xen_path), paste0(lineage, '_chunks')), showWarnings = FALSE)
    saveRDS(xen_sub, file.path(dirname(xen_path), 
                               paste0(lineage, '_chunks'), 
                               paste0(lineage, "_chunk", i, ".rds"))) 

    sampled_tracker[sampled_cells] = 'Sampled'
    logmsg(paste('Chunk', i, 'saved to',  file.path(dirname(xen_path), 
                               paste0(lineage, '_chunks'), 
                               paste0(lineage, "_chunk", i, ".rds"))))
    i <- i + 1 
}




