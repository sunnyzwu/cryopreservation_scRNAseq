# kBET 
# v9 compare on a bulk downsampled level > downsampled to the lowest sample number rather than combined downsampled
# run and compute own kNN
#
#
#
#
# 01 PARSE ARGUMENTS  ------------------------------------------------------------

# arguments from command line
temp_start_time <- Sys.time()
print(temp_start_time)
temp_args <-
  commandArgs(trailingOnly = T)

{
# 01 SAMPLE ID
temp_sampleIDs <- 
  temp_args[1]

print(paste0("Sample ID ", temp_sampleIDs))

# 02 NUMBER OF CELLS FOR DOWNSAMPLING
temp_downsampled_num <- 
  temp_args[2]
print(paste0("Downsampled cells (ignored if run by cluster) ", temp_downsampled_num))

if(temp_downsampled_num == "lowest"){
  temp_downsampled_to_lowest <- T
} else temp_downsampled_to_lowest <- F

if(!temp_downsampled_num == "lowest"){
  temp_downsampled_num <- 
    as.numeric(temp_downsampled_num) # this has to be set if parsing in numbers stated in " " from shell script. Will treat as character class and break functions and conditional statements
} 


# 04 RUN KBET BY CLUSTER SUBSTRUCTURE ?
temp_run_by_cluster <- 
  temp_args[3]
temp_run_by_cluster <-
  as.logical(temp_run_by_cluster)

print(paste0("Run by cluster (recommened) ", temp_run_by_cluster))


# 05 NUMBER OF CELLS TO DOWNSAMPLE PER CONDITION PER CLUSTER
temp_downsampled_num_per_clust <- 
  temp_args[4]
temp_downsampled_num_per_clust <- 
  as.numeric(temp_downsampled_num_per_clust) # this has to be set if parsing in numbers stated in " " from shell script. Will treat as character class and break functions and conditional statements

print(paste0("Cells to downsample by cluster ", temp_downsampled_num_per_clust))

# 06 RUN BY INTEGRATED MATRIX
temp_run_by_integrated <- 
  temp_args[5]
temp_run_by_integrated <-
  as.logical(temp_run_by_integrated)
print(paste0("Run by integrated matrix ", temp_run_by_integrated))

# 06 RUN BY SEURAT METRICS ?
temp_run_silloutte_onseurat_object <- 
  temp_args[6]
temp_run_silloutte_onseurat_object <-
  as.logical(temp_run_silloutte_onseurat_object)
print(paste0("Run by seurat metrics ?", temp_run_silloutte_onseurat_object))

  }

# 02 SETUP ------------------------------------------------------------

library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library('FNN')
library(ggfortify)
library(cluster)
library(stringr)

temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 8, 
      height = 8, 
      res = 300, 
      units = 'in'
    )
  }

temp_png_function_2 <-
  function(x) {
    png(
      file = (x), 
      width = 12, 
      height = 6, 
      res = 300, 
      units = 'in'
    )
  }



# 03 LOAD SAMPLES -----------------------------------------------------

for(sampleID in temp_sampleIDs) {
  
  temp_objectstring <- "/Output/Rdata/04_seurat_CCA_aligned_processed.Rdata"
  
  temp_seurat_10X <- readRDS(paste0("/share/ScratchGeneral/sunwu/projects/cryopreservation_paper/seurat_v3/02_alignment_processing_run02_with_seuratmerge/run02_newPCAs_for_merge/output/CCA_",
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
  
  # print cells per sample
  print(table(temp_seurat_10X@meta.data$orig.ident))
  
  temp_df_cellnumber <-  as.data.frame(table(temp_seurat_10X@meta.data$condition))
  temp_df_cellnumber <- arrange(temp_df_cellnumber,
                                Freq)
  temp_lowest_sample_ID <- 
    temp_df_cellnumber[1,1]
  temp_lowest_sample_ID_cellnum <- 
    temp_df_cellnumber[1,2]
  
  # UMAP
  # print UMAP of raw cells
  temp_dimplot <- DimPlot(
    object = temp_seurat_10X,
    group.by = "int_PCAMERGE_C_res.0.4",
    pt.size = 0.5,
    reduction = "UMAPMERGEC",
    label=T
  )
  temp_png_function(paste0("01_UMAP_raw_cells_res_C_0.4.png"))
  print(temp_dimplot)
  dev.off()
  
  temp_dimplot <- DimPlot(
    object = temp_seurat_10X,
    group.by = "condition",
    pt.size = 0.5,
    reduction = "UMAPMERGEC"
  )
  temp_png_function(paste0("01_UMAP_raw_cells_condition_.png"))
  print(temp_dimplot)
  dev.off()
  
  n <- paste0("seurat_10X_", sampleID)
  assign(n, temp_seurat_10X)
  rm(temp_seurat_10X)
}



# 04 GENERATE DATAFRAME FOR DOWNSAMPLING ------------------------------------------

#data: a matrix (rows: samples, columns: features (genes))
#batch: vector or factor with batch label of each cell 
# batch.estimate <- kBET(data, batch)

for(sampleID in temp_sampleIDs) {
  
  temp_seurat_object <- 
    get(paste0("seurat_10X_", sampleID))
  
  if(temp_sampleIDs == "PID20033") {
    Idents(temp_seurat_object) <- "condition"
    temp_seurat_object <- subset(temp_seurat_object,
                                 idents = c("Fresh Tissue", "Cryopreserved Tissue"))
    temp_seurat_object@meta.data$condition <- factor(temp_seurat_object@meta.data$condition,
                                                     levels= unique(temp_seurat_object@meta.data$condition))
    } 
  
  
  temp_df <- data.table(cellbarcode = row.names(temp_seurat_object@meta.data),
                        batch = temp_seurat_object@meta.data$condition,
                        cluster = temp_seurat_object@meta.data$int_PCAMERGE_C_res.0.4)
  
  if(temp_downsampled_to_lowest){
    temp_downsampled_num <- min(table(temp_seurat_object@meta.data$condition))
    temp_downsampled_num <- temp_downsampled_num/2
    print(paste0("downsampling by bulk to ",temp_downsampled_num))
  }
  
  
  n <- paste0("temp_df_",sampleID)
  assign(n, temp_df)
  
}


# 05 STRATIFIED DOWNSUBSAMPLING TO LOWEST BATCH SIZE --------------------------------------------------

# example at https://stackoverflow.com/questions/23479512/stratified-random-sampling-from-data-frame
# set.seed(1)
# n <- 1e4
# d <- data.table(age  = sample(1:5, n, T),
#                 lc   = rbinom(n,   1, .5),
#                 ants = rbinom(n,   1, .7))
# 
# out <- d[, .SD[sample(1:.N, 30)], by=.(age, lc)]

for(sampleID in temp_sampleIDs) {
  
  temp_df <- get(paste0("temp_df_",sampleID))
  
  temp_df_subset <- temp_df
  # print("downsample to lowest condition number")
  # # first downsample to match lowest condition number of cells
  # temp_df_keep_barcodes <- 
  #   temp_df[temp_df$batch %in% temp_lowest_sample_ID,]
  # temp_df <- 
  #   temp_df[!temp_df$batch %in% temp_lowest_sample_ID,]
  # 
  # set.seed(1)
  # temp_df_subset <- 
  #   temp_df[, .SD[sample(1:.N,temp_lowest_sample_ID_cellnum)], by=.(batch)] 
  # 
  # temp_df_subset <- rbind(temp_df_subset,
  #                         temp_df_keep_barcodes)
  print(table(temp_df_subset$batch))
  
  # second downsample to x num cells if not running by cluster
  if(temp_run_by_cluster == F) {
  set.seed(999)
  temp_df_subset <- 
    temp_df_subset[, .SD[sample(1:.N,temp_downsampled_num)], by=.(batch)] 
  
  m <- paste0("temp_df_subset_",sampleID)
  assign(m, temp_df_subset)
  }
  
  if(temp_run_by_cluster) {
    
    temp_clusters_to_use <- 
      as.vector(unique(temp_df_subset$cluster))
    
    # extract each cluster and balance cells from all conditions
    temp_df_subset_combined <- NULL
    temp_clusternum_per_condition <- data.frame(row.names = unique(temp_df_subset$batch))
    for(cluster in temp_clusters_to_use) {
      print(cluster)
      temp_df_subset_cluster <- 
        as.data.frame(temp_df_subset)[temp_df_subset$cluster %in% cluster,]
      # append table of cell numbers per condition
      {
        temp_clusternum_per_condition_subset <- 
          data.frame(table(temp_df_subset_cluster$batch))
        
        temp_clusternum_per_condition_subset <- 
          data.frame(cluster = temp_clusternum_per_condition_subset$Freq)
        
        colnames(temp_clusternum_per_condition_subset) <- cluster
        temp_clusternum_per_condition <- cbind(temp_clusternum_per_condition,
                                               temp_clusternum_per_condition_subset)
        
      }
      temp_df_subset_cluster <- as.data.table(temp_df_subset_cluster)
      
      temp_downsampled_num <- as.numeric(min(table(temp_df_subset_cluster$batch)))
      print(paste0("  smallest cluster per batch = ",temp_downsampled_num, " cells"))
      
      if(temp_downsampled_num > temp_downsampled_num_per_clust) {
        print(paste0("   capping to ",temp_downsampled_num, " cells"))
        temp_downsampled_num <- temp_downsampled_num_per_clust
      } 
      print(paste0("  downsampling to ",temp_downsampled_num, " cells"))
      
      set.seed(423)
      temp_df_subset_cluster <- 
        temp_df_subset_cluster[, .SD[sample(1:.N,temp_downsampled_num)], by=.(batch)] 
    
    print(table(temp_df_subset_cluster$batch))
    
    temp_df_subset_combined <- rbind(temp_df_subset_combined,
                                     temp_df_subset_cluster)
    }
    

    
    # remove all cluster < x cells
    temp_clusternum_per_condition <- 
      t(temp_clusternum_per_condition)
    temp_clusternum_per_condition_filtered <-
      temp_clusternum_per_condition[temp_clusternum_per_condition[,1] > (temp_downsampled_num_per_clust*2),]
    temp_clusternum_per_condition_filtered <-
      temp_clusternum_per_condition_filtered[temp_clusternum_per_condition_filtered[,2] > (temp_downsampled_num_per_clust),]
    if(!temp_sampleIDs == "PID20033") {
      
    temp_clusternum_per_condition_filtered <-
      temp_clusternum_per_condition_filtered[temp_clusternum_per_condition_filtered[,3] > (temp_downsampled_num_per_clust),]
    }
    
    # clusters to use (all > x)
    temp_clusters_to_use <- 
      as.vector(row.names(temp_clusternum_per_condition_filtered))
    temp_df_subset_combined$cluster <- 
      factor(temp_df_subset_combined$cluster,
             levels=temp_clusters_to_use)
    temp_df_subset_combined <- 
      temp_df_subset_combined[temp_df_subset_combined$cluster %in% temp_clusters_to_use,]
    print(paste0("clusters kept based on cell number cutoff ",temp_downsampled_num_per_clust))
    print(temp_clusters_to_use)
    
    
    print(table(temp_df_subset_combined$cluster))
    print(table(temp_df_subset_combined$batch))
    
  j <- paste0("temp_clusters_to_use_", sampleID)
  assign(j, temp_clusters_to_use)
  
  m <- paste0("temp_df_subset_",sampleID)
  assign(m, temp_df_subset_combined)
  }
    
}


# downsample an additional cells from the Fresh tissue dataset
for(sampleID in temp_sampleIDs) {
  
  temp_df_subset <- 
    get(paste0("temp_df_subset_",sampleID))
  
  temp_fresh_barcodes <- 
    temp_df_subset[temp_df_subset$batch %in% "Fresh Tissue",]
  
  temp_fresh_barcodes <- 
    temp_fresh_barcodes$cellbarcode
  
  temp_seurat_object <- 
    get(paste0("seurat_10X_", sampleID))
  
  temp_barcodes <- 
    row.names(temp_seurat_object@meta.data)
  
  temp_barcodes <- 
    temp_barcodes[!temp_barcodes %in% temp_fresh_barcodes]
  
  temp_seurat_object_removed_fresh_barcodes1 <- subset(temp_seurat_object,
                                                       cells = temp_barcodes)
  
  Idents(temp_seurat_object_removed_fresh_barcodes1) <- "condition"
  temp_seurat_object_removed_fresh_barcodes1 <- subset(temp_seurat_object_removed_fresh_barcodes1,
                                                       idents = "Fresh Tissue")
  
  temp_df <- data.table(cellbarcode = row.names(temp_seurat_object_removed_fresh_barcodes1@meta.data),
                        batch = temp_seurat_object_removed_fresh_barcodes1@meta.data$condition,
                        cluster = temp_seurat_object_removed_fresh_barcodes1@meta.data$int_PCAMERGE_C_res.0.4)
  
  if(temp_run_by_cluster == F) {
    set.seed(1)
    temp_df_subset <- 
      temp_df[, .SD[sample(1:.N,temp_downsampled_num)]]
  }
  
  if(temp_run_by_cluster) {
  
    temp_clusters_to_use <- 
      get(paste0("temp_clusters_to_use_", sampleID))
    temp_df$cluster <- 
      factor(temp_df$cluster,
             levels=temp_clusters_to_use)
    # filter for clusters to use
    temp_df <- 
      temp_df[temp_df$cluster %in% temp_clusters_to_use,]
    
    # get downsampled dataframe to identify how many cells to downsample per cluster
    temp_df_subset <- 
      get(paste0("temp_df_subset_",sampleID))
    temp_cellnum_table <- 
      as.data.frame(table(temp_df_subset$cluster))
    temp_cellnum_table$downsampleFT <- 
      temp_cellnum_table$Freq/3
    
    # subset each cluster
    temp_df_subset_combined <- NULL
    for(cluster in temp_clusters_to_use) {
      print(cluster)
      temp_df_subset_cluster <- 
        as.data.frame(temp_df)[temp_df$cluster %in% cluster,]
      
      temp_df_subset_cluster <- as.data.table(temp_df_subset_cluster)
      
      temp_downsampled_num <- 
        temp_cellnum_table$downsampleFT[temp_cellnum_table$Var1 == cluster]
      
      set.seed(572)
      temp_df_subset_cluster <- 
        temp_df_subset_cluster[, .SD[sample(1:.N,temp_downsampled_num)], by=.(batch)] 
      
      print(table(temp_df_subset_cluster$batch))
      
      temp_df_subset_combined <- rbind(temp_df_subset_combined,
                                       temp_df_subset_cluster)
    }
    
    print(table(temp_df_subset_combined$cluster))
    print(table(temp_df_subset_combined$batch))
    temp_df_subset <- temp_df_subset_combined   
    }
    
  
  n <- paste0("temp_df_subset_FT_",sampleID)
  assign(n, temp_df_subset)
  
  rm(temp_seurat_object_removed_fresh_barcodes1)
  
}



# 06 GENERATION OF ALL DATAFRAMES --------------------------------------------

temp_condition_ids <- c("Fresh Tissue", "Cryopreserved Cell Suspension", "Cryopreserved Tissue")
if(temp_sampleIDs == "SCC180161") {
  temp_condition_ids <- c("Fresh Tissue", "Cryopreserved Tissue", "Cryopreserved Overnight") } 
if(temp_sampleIDs == "PID20033") {
  temp_condition_ids <- c("Fresh Tissue", "Cryopreserved Tissue") } 
if(temp_sampleIDs == "CID4523") {
  temp_condition_ids <- c("Fresh Tissue", "Cryopreserved Cell Suspension") } 


for(sampleID in temp_sampleIDs) {
  
  for(type in temp_condition_ids) {
    print(type)
    
    temp_df_subset_FT <- 
      get(paste0("temp_df_subset_FT_",sampleID))
    temp_df_subset <- 
      get(paste0("temp_df_subset_",sampleID))
    
    temp_df <- temp_df_subset[temp_df_subset$batch %in% type,]
    
    # set fresh tissue downsampled data 2 as different batch
    if(type == "Fresh Tissue") {
      temp_df$batch <- "Fresh Tissue 2"
    }
    
    temp_df_combined <- rbind(temp_df_subset_FT,
                              temp_df)
    
    print(table(temp_df_combined$batch))
    print(unique(temp_df_combined$batch))
    print("unique barcodes")
    print(length(unique(temp_df_combined$cellbarcode)))
    print(dim(temp_df_combined))
    
    if(temp_run_by_cluster)  {
      print(table((temp_df_combined$cluster)))
    }
    
    if(type == "Fresh Tissue") {
      temp_name <- "FT"
    }
    if(type == "Cryopreserved Cell Suspension") {
      temp_name <- "CCS"
    }
    if(type == "Cryopreserved Tissue") {
      temp_name <- "CT"
    }
    if(type == "Cryopreserved Overnight") {
      temp_name <- "CO"
    }
    
    n <- paste0("input_kBET_",temp_name)
    assign(n, temp_df_combined)
    
  }
  
}

# 07 EXPORT MATRIX WITH DOWNSAMPLED DATA -------------------------------------

#data: a matrix (rows: samples, columns: features (genes))
#batch: vector or factor with batch label of each cell 
# batch.estimate <- kBET(data, batch)

temp_condition_ids <- c("FT", "CCS", "CT") 

if(temp_sampleIDs == "SCC180161") {
  temp_condition_ids <- c("FT", "CT", "CO") } 

if(temp_sampleIDs == "PID20033") {
  temp_condition_ids <- c("FT", "CT") } 

if(temp_sampleIDs == "CID4523") {
  temp_condition_ids <- c("FT", "CCS") }

for(sampleID in temp_sampleIDs) {
  
  temp_seurat_object <- 
    get(paste0("seurat_10X_", sampleID))
  
  for(type in temp_condition_ids) {
    
    print(type)
    
    temp_df_subset <- get(paste0("input_kBET_",type))
    print(dim(temp_df_subset))
    
    if(temp_run_by_cluster == F)  {
    temp_seurat_object_subset <- subset(temp_seurat_object, 
                                        cells = temp_df_subset$cellbarcode)
    }
    
    if(temp_run_by_cluster)  {
    temp_seurat_object_subset <- subset(temp_seurat_object, 
                                          cells = temp_df_subset$cellbarcode)
      
    temp_clusters_to_use <- 
      get(paste0("temp_clusters_to_use_", sampleID))

    Idents(temp_seurat_object_subset) <- "int_PCAMERGE_C_res.0.4"
    temp_seurat_object_subset <- subset(temp_seurat_object_subset,
                                        ident = temp_clusters_to_use)
    
    temp_seurat_object_subset@meta.data$int_PCAMERGE_C_res.0.4 <-
      factor(temp_seurat_object_subset@meta.data$int_PCAMERGE_C_res.0.4,
             levels = temp_clusters_to_use)
    }
    
    # check FT barcodes are all unique
    if(type == "FT"){
      temp_df_subset <- data.frame(row.names = temp_df_subset$cellbarcode,
                                   condition = temp_df_subset$batch)
      temp_df_subset_sorted <- temp_df_subset[row.names(temp_seurat_object_subset@meta.data),,drop=F]
      print("For FT, matching cellbarcodes (should be TRUE)")
      print(identical(row.names(temp_df_subset_sorted),
                row.names(temp_seurat_object_subset@meta.data)))
      temp_seurat_object_subset <- AddMetaData(temp_seurat_object_subset, 
                                               metadata = temp_df_subset_sorted)
    }
    
    # print UMAP of raw cells
    temp_dimplot <- DimPlot(
      object = temp_seurat_object_subset,
      group.by = "int_PCAMERGE_C_res.0.4",
      pt.size = 0.5,
      reduction = "UMAPMERGEC",
      label=T
    )
    temp_png_function(paste0("02_UMAP_downsampled_cells_",type,"_cluster.png"))
    print(temp_dimplot)
    dev.off()
    
    temp_dimplot <- DimPlot(
      object = temp_seurat_object_subset,
      group.by = "condition",
      pt.size = 0.5,
      reduction = "UMAPMERGEC"
    )
    temp_png_function(paste0("02_UMAP_downsampled_cells_",type,"_condition.png"))
    print(temp_dimplot)
    dev.off()
    
    print(table(temp_seurat_object_subset@meta.data$condition))
    if(temp_run_by_cluster)  {
      print(table((temp_seurat_object_subset@meta.data$int_PCAMERGE_C_res.0.4)))
    }
    
    # export matrix
    if(temp_run_by_integrated == F) {
      print("RNA assay")
    temp_matrix <- GetAssayData(object = temp_seurat_object_subset,
                                assay = "RNA",
                                slot = "data")
    print(paste0("number of features = ", length(rownames(temp_matrix))))
    }
    if(temp_run_by_integrated) {
      print("integrated assay")
      temp_matrix <- GetAssayData(object = temp_seurat_object_subset,
                                  assay = "integrated",
                                  slot = "data")
      print(paste0("number of features = ", length(rownames(temp_matrix))))
    }
    temp_matrix <- t(as.data.frame(temp_matrix))
    
    temp_batch <- temp_seurat_object_subset@meta.data$condition
    
    temp_cluster <- 
      temp_seurat_object_subset@meta.data$int_PCAMERGE_C_res.0.4
    
    h <- paste0("temp_cluster_",sampleID, "_", type)
    assign(h, temp_cluster)

    n <- paste0("temp_matrix_",sampleID, "_", type)
    assign(n, temp_matrix)
    
    m <- paste0("temp_batch_",sampleID, "_", type)
    assign(m, temp_batch)
    
    j <- paste0("temp_seurat_object_subset_",type)
    assign(j, temp_seurat_object_subset)
    
  }
}



# 08 RUN METRICS ON SEURAT PCA ----------------------------------------------------------------------

if(temp_run_silloutte_onseurat_object){
  library(cluster, quietly = TRUE)
  
  reduction <- "pca"
  dims <- 1:100
  
  for(sampleID in temp_sampleIDs) {
    
    temp_df <- NULL
    for(type in temp_condition_ids) {
      
      print(type)
      
      temp_seurat_object_subset <- 
        get(paste0("temp_seurat_object_subset_", type))
      
      # remove clusters < 50 cells
      temp_colname <- "int_PCAMERGE_C_res.0.4"
      print(table(temp_seurat_object_subset@meta.data[,temp_colname]))
      temp_df_cellnum <- as.data.frame(table(temp_seurat_object_subset@meta.data[,temp_colname]))
      temp_df_filtered <- temp_df_cellnum[temp_df_cellnum$Freq > 50,]
      temp_df_filtered <- as.vector(temp_df_filtered$Var1)
      temp_removed <- as.vector(temp_df_cellnum$Var1[!temp_df_cellnum$Var1 %in% temp_df_filtered])
      
      print("trimming samples")
      Idents(temp_seurat_object_subset) <- temp_colname
      temp_seurat_object_subset <- SubsetData(temp_seurat_object_subset,
                                                   ident.use=temp_df_filtered)
      print(temp_removed)
      
      # compute
      dist.matrix <- dist(x = Embeddings(object = temp_seurat_object_subset[[reduction]])[, dims])
      clusters <- temp_seurat_object_subset$condition
      sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
      temp_seurat_object_subset$sil <- sil[, 3]
      
      # mixing metric
      max.k <- 300
      mm <- max.k - MixingMetric(object = temp_seurat_object_subset, 
                                 grouping.var = "condition", 
                                 reduction = reduction, 
                                 dims = dims,
                                 max.k = max.k)
      
      DefaultAssay(object = temp_seurat_object_subset) <- "RNA"
      # Local structure preservation
      ls <- LocalStruct(object = temp_seurat_object_subset,
                        grouping.var = "condition", 
                        reduction = reduction, 
                        reduced.dims = dims, 
                        orig.dims = 1:30)
      ls <- unname(obj = unlist(x = ls))
      
      all.metrics <- list(
        silhouette = temp_seurat_object_subset$sil, 
        mixing.metric = mm,
        local.struct = ls
      )
      
      temp_df_subset <- as.data.frame(temp_seurat_object_subset$sil)
      colnames(temp_df_subset) <- "silhouette"
      temp_df_subset$mixing.metric <-  mm
      temp_df_subset$local.struct <- ls
      temp_df_subset$condition <- type
      temp_df_subset$cluster <- temp_seurat_object_subset@meta.data$int_PCAMERGE_C_res.0.4
      temp_df <- rbind(temp_df,temp_df_subset)
      
    }
  
    # plot
    if(sampleID == "SCC180161") {
      temp_colours <- c("#1b9e77", "#7570b3", "#abd9e9")
    }  else temp_colours <- c("#1b9e77", "#d95f02", "#7570b3")
    
    temp_df$condition <- factor(temp_df$condition, 
                                levels=unique(temp_df$condition))
    temp_df$cluster <- factor(temp_df$cluster, 
                                levels=str_sort(unique(temp_df$cluster), numeric = T))

    for(metrics in c("silhouette","mixing.metric","local.struct")){
      print(metrics)
      temp_df_subset <- temp_df[,c(metrics,"condition", "cluster")]
      colnames(temp_df_subset) <- c("value","condition", "cluster")
      
      temp_ggplot_clusters <- ggplot(temp_df_subset, 
                                     aes(condition, 
                                         value,
                                         label = condition,
                                         fill = condition)) + 
        geom_boxplot(outlier.shape = NA, 
                     size=0.3, 
                     alpha = 0.5) + 
        labs(x='Cluster', 
             y=paste0(metrics, 'Score'),
             title=metrics) +
        theme_bw() +
        # scale_y_continuous(limits=c(-0.5,0.5)) +
        theme(axis.text.x = element_text(angle=25, vjust = .25, size=8)) +
        scale_fill_manual(values = temp_colours) +
        facet_wrap(. ~ cluster,
                   scales="free_x",
                   nrow =1) # free x removes white space
        
      
      temp_png_function_2(paste0("07_",metrics,"_scores_from_seurat_object_by_cluster.png"))
      print(temp_ggplot_clusters)
      dev.off()
      
      temp_ggplot_clusters <- ggplot(temp_df_subset, 
                                     aes(condition, 
                                         value,
                                         label = condition,
                                         fill = condition)) + 
        geom_boxplot(outlier.shape = NA, 
                     size=0.3, 
                     alpha = 0.5) + 
        labs(x='Cluster', 
             y=paste0(metrics, 'Score'),
             title=metrics) +
        theme_bw() +
        # scale_y_continuous(limits=c(-0.5,0.5)) +
        theme(axis.text.x = element_text(angle=25, vjust = .25, size=8)) +
        scale_fill_manual(values = temp_colours) +
        temp_png_function(paste0("07_",metrics,"_scores_from_seurat_object_averaged.png"))
        print(temp_ggplot_clusters)
        dev.off()
      
      
      write.csv(temp_df, "08_silhouette_mixing_LSC_metrics.csv")
    }
    
    
    }
  
}









