# GO ENRICHMENT PER CRYO CONDITION
# SCRIPT A 
# R.v3.5.1 (SEURAT V3)
# DIFFERENTIAL GENE EXPRESSION IN SEURAT V3
# 01 PARSE ARGUMENTS  ------------------------------------------------------------

temp_start_time <- Sys.time()
print(temp_start_time)
temp_args <-
  commandArgs(trailingOnly = T)

# 01 CELL TYPE
temp_sampleIDs  <- 
  temp_args[1]

# 02 SETUP  ------------------------------------------------------------------

library(Seurat)
library(dplyr)

# 03 LOAD SAMPLES -----------------------------------------------------

for(sampleID in temp_sampleIDs) {
  
  temp_objectstring <- "/Output/Rdata/05_seurat_CCA_aligned_processed_downsampled.Rdata"
  
  temp_seurat_10X <- readRDS(paste0("/share/ScratchGeneral/sunwu/projects/cryopreservation_paper/seurat_v3/02_alignment_processing_run02_with_seuratmerge/run05_default_clus_res/output/CCA_",
                                    sampleID, 
                                    temp_objectstring))
  
  temp_seurat_10X@meta.data$condition <- NULL
  
  if(!sampleID == "SCC180161") {
    temp_seurat_10X@meta.data$condition[temp_seurat_10X@meta.data$collection == "Surgery"] <- "Fresh Tissue"
    temp_seurat_10X@meta.data$condition[temp_seurat_10X@meta.data$collection == "Cryopreserved Cell Sus"] <- "Cryopreserved Cell Suspension"
    temp_seurat_10X@meta.data$condition[temp_seurat_10X@meta.data$collection == "Cryopreserved Tissue"] <- "Cryopreserved Tissue"
    
    temp_seurat_10X@meta.data$condition <- factor(temp_seurat_10X@meta.data$condition,
                                                  levels = (c("Fresh Tissue", "Cryopreserved Cell Suspension", "Cryopreserved Tissue")))
  }
  
  if(sampleID == "SCC180161") {
    temp_seurat_10X@meta.data$condition[temp_seurat_10X@meta.data$orig.ident == "SCC180161"] <- "Fresh Tissue"
    temp_seurat_10X@meta.data$condition[temp_seurat_10X@meta.data$orig.ident == "SCC180161SD"] <- "Cryopreserved Tissue"
    temp_seurat_10X@meta.data$condition[temp_seurat_10X@meta.data$orig.ident == "SCC180161ND"] <- "Cryopreserved Overnight"
    
    temp_seurat_10X@meta.data$condition <- factor(temp_seurat_10X@meta.data$condition,
                                                  levels = (c("Fresh Tissue", "Cryopreserved Tissue", "Cryopreserved Overnight")))
    
  }
  
  Idents(object = temp_seurat_10X) <- "condition"
  
  n <- paste0("seurat_10X_", sampleID)
  assign(n, temp_seurat_10X)
  rm(temp_seurat_10X)
}


# 04 SUBSET CONDITIONS ----------------------------------------------------

for(sampleID in temp_sampleIDs) {
  
  temp_seurat_10X <- 
    get(paste0("seurat_10X_", sampleID))
  
  temp_conditions <- 
    unique(temp_seurat_10X@meta.data$condition)
  
  Idents(object = temp_seurat_10X) <- "condition"

  for(condition in temp_conditions){
    
    temp_seurat_10X_subset <- subset(temp_seurat_10X,
                                     idents=condition)
    
    n <- paste0("temp_seurat_10X_subset_",condition)
    assign(n, temp_seurat_10X_subset)
    
  }
}


# 05 DIFFERENTIAL GENE EXPRESSION --------------------------------------------

for(sampleID in temp_sampleIDs) {
  
  for(condition in temp_conditions){
    
    temp_seurat_10X_subset <- 
      get(paste0("temp_seurat_10X_subset_",condition))
    
    Idents(object = temp_seurat_10X_subset) <- "int_PCA_B_res.0.8"
    
    temp_cluster_allmarkers <- FindAllMarkers(
      only.pos = T,
      object = temp_seurat_10X_subset,
      min.pct = 0.25, 
      logfc.threshold = 0.25,
      min.diff.pct = 0.1, 
      test.use = 'MAST', 
      print.bar = T
    ) 
    
    temp_cluster_allmarkers <- arrange(temp_cluster_allmarkers,
                                       (cluster),
                                       desc(avg_logFC))      
    
    # n <- paste0("temp_cluster_allmarkers_", condition)  
    # assign(n, temp_cluster_allmarkers)
    
    write.csv(temp_cluster_allmarkers,
              paste0("01_DGE_",sampleID,"_", condition, ".csv"))
  }
}
  

# FINISH ------------------------------------------------------------------

print("start time")
print(temp_start_time)
print("finish time")
print(Sys.time())
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
