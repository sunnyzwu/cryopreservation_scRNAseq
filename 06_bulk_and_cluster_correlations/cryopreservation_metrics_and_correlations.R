# CRYOPRESERVATION PAPER METRICS AND CORRELATION 
# SEURAT V5
#
# MERGE REDUCTIONS AND CLUSTER ID = UMAPMERGEB & int_PCAMERGE_B_res.0.8
# INTEGRATED REDUCTIONS AND CLUSTER ID = UMAPB & int_PCA_B_res.0.8
#
#
# 01 SET UP ------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(cowplot)
library(reshape2)
library(stringr)
library(data.table)
library(dplyr)

# 02 DIRECTORIES -------------------------------------------------------------

setwd("/share/ScratchGeneral/sunwu/projects/cryopreservation_paper/seurat_v3/05_figures")
dir.create("20200326_paper_figures/")
setwd("20200326_paper_figures")

# 03 LOAD DOWNSAMPLED DATA ---------------------------------------------------

# load downsampled data
for(sampleID in temp_sampleIDs[!temp_sampleIDs %in% "CID4523"]) {
  print(sampleID)
  
  temp_seurat_object <- 
    get(paste0("seurat_10X_", sampleID))
  
  temp_orig_IDs <- unique(temp_seurat_object@meta.data$orig.ident)
  
  if(sampleID == "SCC180161"){
    temp_orig_IDs <- c(temp_orig_IDs[1], temp_orig_IDs[3], temp_orig_IDs[2])
  }
  print(temp_orig_IDs)
  
  # load cellranger aggr data
  temp_seurat_aggr <- Read10X(paste0("/share/ScratchGeneral/sunwu/projects/cryopreservation_paper/seurat_v3/06_cellranger_aggr_downsample_reads/output/cellranger_aggr_",sampleID,"/",sampleID,"/outs/raw_gene_bc_matrices_mex/GRCh38/"))
  
  # append new cell barcodes
  temp_grep_function <-
    function(x) {
      (grep(x,
            colnames(temp_seurat_aggr),
            value = T))
    }
  
  if(!sampleID == "PID20033"){
    temp_data_frame <- data.frame(IDs = (c(
      rep(temp_orig_IDs[1],
          times = length(temp_grep_function("-1"))),
      rep(temp_orig_IDs[2],
          times = length(temp_grep_function("-2"))),
      rep(temp_orig_IDs[3],
          times = length(temp_grep_function("-3")))
    )),
    row.names = (
      c(
        temp_grep_function("-1"),
        temp_grep_function("-2"),
        temp_grep_function("-3")
      )
    ))
    
    temp_data_frame$newbarcode <- 
      paste0(temp_data_frame$IDs, 
             "_",
             rownames(temp_data_frame))
    
    print(all.equal(colnames(temp_seurat_aggr),
                    rownames(temp_data_frame)))
    
    temp_data_frame$newbarcode <- gsub("-1",
                                       "",
                                       temp_data_frame$newbarcode )
    temp_data_frame$newbarcode <- gsub("-2",
                                       "",
                                       temp_data_frame$newbarcode )
    temp_data_frame$newbarcode <- gsub("-3",
                                       "",
                                       temp_data_frame$newbarcode )
    
    colnames(x = temp_seurat_aggr) <- 
      temp_data_frame$newbarcode 
  }
  
  if(sampleID == "PID20033"){
    temp_data_frame <- data.frame(IDs = (c(
      rep(temp_orig_IDs[1],
          times = length(temp_grep_function("-1"))),
      rep(temp_orig_IDs[2],
          times = length(temp_grep_function("-2")))
    )),
    row.names = (
      c(
        temp_grep_function("-1"),
        temp_grep_function("-2")
      )
    ))
    
    temp_data_frame$newbarcode <- 
      paste0(temp_data_frame$IDs, 
             "_",
             rownames(temp_data_frame))
    
    print(all.equal(colnames(temp_seurat_aggr),
                    rownames(temp_data_frame)))
    
    temp_data_frame$newbarcode <- gsub("-1",
                                       "",
                                       temp_data_frame$newbarcode )
    temp_data_frame$newbarcode <- gsub("-2",
                                       "",
                                       temp_data_frame$newbarcode )
    
    colnames(x = temp_seurat_aggr) <- 
      temp_data_frame$newbarcode 
  }
  
  temp_seurat_object_aggr <-  CreateSeuratObject(
    counts = temp_seurat_aggr,
    min.cells = 1, 
    min.features = 0
  )
  
  # mito content
  temp_seurat_object_aggr[["percent.mito"]] <- PercentageFeatureSet(temp_seurat_object_aggr, 
                                                                    pattern = "^MT-")
  
  # temp_overlap_barcodes <- colnames(temp_seurat_object_aggr)[colnames(temp_seurat_object_aggr) %in% colnames(temp_seurat_object)]
  
  
  # append downsampled nGenes and UMIs to original object
  temp_seurat_object_aggr_filtered <- subset(temp_seurat_object_aggr,
                                             cells=colnames(temp_seurat_object))
  
  temp_df <- temp_seurat_object_aggr_filtered@meta.data
  temp_df <- temp_df[rownames(temp_seurat_object@meta.data),,drop=F]
  
  print(all.equal(rownames(temp_seurat_object@meta.data),
                  rownames(temp_df)))
  
  temp_df <- temp_df[,colnames(temp_df) %in% c("nCount_RNA", "nFeature_RNA", "percent.mito"),drop=T]
  colnames(temp_df) <- paste0("downsampled_", colnames(temp_df))
  
  temp_seurat_object <- AddMetaData(temp_seurat_object,
                                    temp_df)
  
  
  # append metadata to downsampled object
  temp_df <- temp_seurat_object@meta.data
  temp_df <- temp_df[rownames(temp_seurat_object_aggr_filtered@meta.data),,drop=F]
  
  print(all.equal(rownames(temp_seurat_object_aggr_filtered@meta.data),
                  rownames(temp_df)))
  
  temp_seurat_object_aggr_filtered <- AddMetaData(temp_seurat_object_aggr_filtered,
                                                  temp_df)
  
  temp_seurat_object_aggr_filtered <- NormalizeData(temp_seurat_object_aggr_filtered)
  
  n <- paste0("seurat_10X_downsampledobject_",sampleID)
  assign(n, temp_seurat_object_aggr_filtered)
  
  n <- paste0("seurat_10X_downsampled_",sampleID)
  assign(n, temp_seurat_object)
}


# 04 METRICS GGPLOT DOWNSAMPLED STATS AND PLOT ------------------------------------------------

dir.create("03_nGenes_nUMIs")
# p-val cutoffs
symnum.args <- list(cutpoints = c(0, 0.00001, 0.0001, 0.001, 0.01, 1),
                    symbols = c("****", "***", "**", "*", "ns"))

#
# In other words, we use the following convention for symbols
# indicating statistical significance:
#   
#   * 'ns': p > 0.01
# 
# * '*': p <= 0.01
# 
# * '**': p <= 0.001
# 
# * '***': p <= 0.0001
# 
# * '****': p <= 0.00001

# plot
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x), 
      width = 12, 
      height = 8, 
      useDingbats = F
    )
  }

temp_png_function <-
  function(x) {
    png(
      file = (x),
      width = 12,
      height = 8,
      res = 300,
      units = 'in'
    )
  }
for(plottype in c("Genes", "UMIs", "Mitocondrial")) {
  if(plottype == "Genes") {colname <- "downsampled_nFeature_RNA"}
  if(plottype == "UMIs") {colname <- "downsampled_nCount_RNA"}
  if(plottype == "Mitocondrial") {colname <- "downsampled_percent.mito"}
  
  temp_gene_data_df <- NULL
  for(sampleID in temp_sampleIDs[!temp_sampleIDs %in% "CID4523"]) {
    
    temp_seurat_object <- 
      get(paste0("seurat_10X_downsampled_", sampleID))
    
    temp_df <- data.frame(temp_seurat_object@meta.data[,colname],
                          temp_seurat_object@meta.data[,"orig.ident"],
                          temp_seurat_object@meta.data[,"condition"])
    
    colnames(temp_df) <- c(colname, "orig.ident", "condition")
    
    temp_gene_data_df <- rbind(temp_gene_data_df,
                               temp_df)
  }
  
  # patient
  temp_gene_data_df$patient <- NA
  for(i in unique(temp_gene_data_df$orig.ident)) {
    if(i %in% c("CID4471", "CID4471CELLSUS", "CID4471CHUNKS")){temp_pid = "BC P1"}
    if(i %in% c("CID44971", "CID44971CELLSUS", "CID44971CHUNKS")){temp_pid = "BC P2"}
    if(i %in% c("CID4513", "CID4513CELLSUS", "CID4513CHUNKS")){temp_pid = "BC P3"}
    # if(i %in% c("CID4523", "CID4523CELLSUS", "CCID4523CHUNKS")){temp_pid = "BC P4"}
    if(i %in% c("PID17267", "PID17267CELLSUS", "PID17267CHUNKS")){temp_pid = "PC P1"}
    if(i %in% c("PID20033", "PID20033CELLSUS", "PID20033CHUNKS")){temp_pid = "PC P2"}
    if(i %in% c("SCC180161", "SCC180161SD", "SCC180161ND")){temp_pid = "M P1"}
    temp_gene_data_df$patient[temp_gene_data_df$orig.ident == paste0(i)] <- temp_pid
  }
  
  temp_gene_data_df$condition <- factor(temp_gene_data_df$condition,
                                        levels = c("Fresh Tissue", 
                                                   "Cryopreserved Cell Suspension", 
                                                   "Cryopreserved Tissue",
                                                   "Cryopreserved Overnight"))
  
  temp_gene_data_df$patient <- factor(temp_gene_data_df$patient,
                                      levels = c("BC P1",
                                                 "BC P2",
                                                 "BC P3",
                                                 # "BC P4",
                                                 "PC P1",
                                                 "PC P2",
                                                 "M P1"
                                      ))
  
  # filter out crap samples
  # temp_gene_data_df <- temp_gene_data_df[! temp_gene_data_df$condition == "CID4523CHUNKS",]
  temp_gene_data_df <- temp_gene_data_df[! temp_gene_data_df$orig.ident == "PID20033CELLSUS",]
  
  
  # unique ID
  temp_gene_data_df$ID <- paste0(temp_gene_data_df$patient,
                                 temp_gene_data_df$condition)
  # temp_comparisons <-  unique(temp_gene_data_df$ID)
  temp_gene_data_df$ID <- factor(temp_gene_data_df$ID,
                                 levels = unique(temp_gene_data_df$ID))
  
  # colours 
  # temp_colours <- c("#1b9e77", "#d95f02", "#abd9e9", "#7570b3")
  
  temp_colours <- c("#1b9e77", "#d95f02", "#7570b3", "#abd9e9")
  
  # Comparisons
  my_comparisons <- list(c("BC P1Fresh Tissue", "BC P1Cryopreserved Cell Suspension"),
                         c("BC P1Fresh Tissue", "BC P1Cryopreserved Tissue"), #
                         c("BC P2Fresh Tissue", "BC P2Cryopreserved Cell Suspension"),
                         c("BC P2Fresh Tissue", "BC P2Cryopreserved Tissue"), #
                         c("BC P3Fresh Tissue", "BC P3Cryopreserved Cell Suspension"),
                         c("BC P3Fresh Tissue", "BC P3Cryopreserved Tissue"), #
                         c("PC P1Fresh Tissue", "PC P1Cryopreserved Cell Suspension"),
                         c("PC P1Fresh Tissue", "PC P1Cryopreserved Tissue"), #
                         c("PC P2Fresh Tissue", "PC P2Cryopreserved Tissue"), #
                         c("M P1Fresh Tissue", "M P1Cryopreserved Cell Suspension"),
                         c("M P1Fresh Tissue", "M P1Cryopreserved Overnight")#
  )
  if(plottype == "Genes"){
    temp_num <- 12000
    temp_factor_to_align_stat_bars <- 0.8
    temp_scale <-  c(0,500,1000,2500,5000,12500)
  } 
  if(plottype == "UMIs"){
    temp_num <- 120000
    temp_factor_to_align_stat_bars <- 0.8
    temp_scale <- c(0,1000,5000,10000,25000,125000)
  }
  if(plottype == "Mitocondrial"){
    temp_num <- 25
    temp_factor_to_align_stat_bars <- 0.9
    temp_scale <- c(0,5,10,20,25)
  }
  
  temp_label_y <- c(temp_num, temp_num*temp_factor_to_align_stat_bars,
                    temp_num, temp_num*temp_factor_to_align_stat_bars,
                    temp_num, temp_num*temp_factor_to_align_stat_bars,
                    temp_num, temp_num*temp_factor_to_align_stat_bars,
                    temp_num,
                    temp_num, temp_num*temp_factor_to_align_stat_bars)
  
  # my_comparisons <- list(c("Fresh Tissue", "Cryopreserved Cell Suspension"),
  #                        c("Fresh Tissue", "Cryopreserved Tissue"), #
  #                        c("Fresh Tissue", "Cryopreserved Overnight")
  # )
  
  # 
  # plot graph with p vals
  if(!plottype == "Mitocondrial"){
    temp_gene_boxplot <- ggboxplot(temp_gene_data_df, 
                                   x = "ID",
                                   y= colname,
                                   fill = "condition",
                                   # group = "patient",
                                   # label="patient",
                                   legend = "none",
                                   # facet.by = "patient", 
                                   repel = F,
                                   nrow=1,
                                   size=1, 
                                   # scales ="free_x",
                                   outlier.size = 0.01
                                   # ylim=c(0,100000)
    )+
      stat_compare_means(method = "t.test",
                         paired = F,
                         label = "p.signif", 
                         # ref.group = "Fresh Tissue",
                         hide.ns = F,
                         comparisons = my_comparisons,
                         bracket.size = 0.1,
                         label.y = temp_label_y,
                         symnum.args = symnum.args) +  # Pairwise comparison against all
      # stat_compare_means(method = "anova", 
      #                    label.y = temp_max_label_y) +
      xlab(" ") + ylab(paste0("Number of ",plottype)) + theme(axis.text.x = element_blank(),
                                                              axis.title.y = element_text(size=20)) +
      scale_fill_manual(values = temp_colours) +
      coord_trans(y = "log") +
      scale_y_continuous(breaks=temp_scale,
                         labels=temp_scale)
  }
  if(plottype == "Mitocondrial"){
    temp_gene_boxplot <- ggboxplot(temp_gene_data_df, 
                                   x = "ID",
                                   y= colname,
                                   fill = "condition",
                                   # group = "patient",
                                   # label="patient",
                                   legend = "none",
                                   # facet.by = "patient", 
                                   repel = F,
                                   nrow=1,
                                   size=1, 
                                   # scales ="free_x",
                                   outlier.size = 0.01
                                   # ylim=c(0,100000)
    )+
      stat_compare_means(method = "t.test",
                         paired = F,
                         label = "p.signif", 
                         # ref.group = "Fresh Tissue",
                         hide.ns = F,
                         comparisons = my_comparisons,
                         bracket.size = 0.1,
                         label.y = temp_label_y,
                         symnum.args = symnum.args) +  # Pairwise comparison against all
      # stat_compare_means(method = "anova", 
      #                    label.y = temp_max_label_y) +
      xlab(" ") + ylab("Percentage Mitocondrial") + theme(axis.text.x = element_blank(),
                                                          axis.title.y = element_text(size=20)) +
      scale_fill_manual(values = temp_colours) +
      scale_y_continuous(breaks=temp_scale,
                         labels=temp_scale)
  }
  
  temp_png_function(paste0("03_nGenes_nUMIs/ggboxplot_",plottype,".png"))
  print(temp_gene_boxplot)
  dev.off()
  
  temp_pdf_function(paste0("03_nGenes_nUMIs/ggboxplot_",plottype,".pdf"))
  print(temp_gene_boxplot)
  dev.off()
  
  
}


# print metrics
for(plottype in c("Genes", "UMIs")) {
  print(plottype)
  if(plottype == "Genes") {colname <- "downsampled_nFeature_RNA"}
  if(plottype == "UMIs") {colname <- "downsampled_nCount_RNA"}
  
  temp_gene_data_df <- NULL
  temp_gene_data_df_summary <- NULL
  for(sampleID in temp_sampleIDs) {
    print(sampleID)
    temp_seurat_object <- 
      get(paste0("seurat_10X_downsampled_", sampleID))
    
    temp_df <- data.frame(temp_seurat_object@meta.data[,colname],
                          temp_seurat_object@meta.data[,"orig.ident"],
                          temp_seurat_object@meta.data[,"condition"])
    
    colnames(temp_df) <- c(colname, "orig.ident", "condition")
    
    print("mean")
    temp_conditions <- unique(temp_df$condition)
    for(condition in temp_conditions) {
      temp_df_summary_subset <- data.frame(sample = sampleID,
                                           condition = condition,
                                           mean = mean(temp_df[temp_df$condition == condition, colname]))
      print(condition)
      print(mean(temp_df[temp_df$condition == condition, colname]))
      
      temp_gene_data_df_summary <- rbind(temp_gene_data_df_summary, 
                                         temp_df_summary_subset)
    }
    
    temp_gene_data_df <- rbind(temp_gene_data_df,
                               temp_df)
    
    j <- paste0("temp_df_",plottype, "_", sampleID)
    assign(j, temp_df)
    
  }
  
  # patient
  temp_gene_data_df$patient <- NA
  for(i in unique(temp_gene_data_df$orig.ident)) {
    if(i %in% c("CID4471", "CID4471CELLSUS", "CID4471CHUNKS")){temp_pid = "BC P1"}
    if(i %in% c("CID44971", "CID44971CELLSUS", "CID44971CHUNKS")){temp_pid = "BC P2"}
    if(i %in% c("CID4513", "CID4513CELLSUS", "CID4513CHUNKS")){temp_pid = "BC P3"}
    if(i %in% c("CID4523", "CID4523CELLSUS", "CCID4523CHUNKS")){temp_pid = "BC P4"}
    if(i %in% c("PID17267", "PID17267CELLSUS", "PID17267CHUNKS")){temp_pid = "PC P1"}
    if(i %in% c("PID20033", "PID20033CELLSUS", "PID20033CHUNKS")){temp_pid = "PC P2"}
    if(i %in% c("SCC180161", "SCC180161SD", "SCC180161ND")){temp_pid = "M P1"}
    temp_gene_data_df$patient[temp_gene_data_df$orig.ident == paste0(i)] <- temp_pid
  }
  
  temp_gene_data_df$condition <- factor(temp_gene_data_df$condition,
                                        levels = c("Fresh Tissue", 
                                                   "Cryopreserved Cell Suspension", 
                                                   "Cryopreserved Tissue",
                                                   "Cryopreserved Overnight"))
  
  temp_gene_data_df$patient <- factor(temp_gene_data_df$patient,
                                      levels = c("BC P1",
                                                 "BC P2",
                                                 "BC P3",
                                                 "BC P4",
                                                 "PC P1",
                                                 "PC P2",
                                                 "M P1"
                                      ))
  
  m <- paste0("temp_gene_data_df_",plottype)
  assign(m, temp_gene_data_df)
  
  n <- paste0("temp_gene_data_df_summary_",plottype)
  assign(n, temp_gene_data_df_summary)
  
  write.csv(temp_gene_data_df_summary,
            paste0("03_nGenes_nUMIs/downsampled_mean_",plottype, ".csv"))
}


# 05 SUPP METRICS PER CLUSTER ------------------------------------------------

dir.create("SUPP_03_nGenes_nUMIs")
# p-val cutoffs
symnum.args <- list(cutpoints = c(0, 0.00001, 0.0001, 0.001, 0.01, 1),
                    symbols = c("****", "***", "**", "*", "ns"))

#
# In other words, we use the following convention for symbols
# indicating statistical significance:
#   
#   * 'ns': p > 0.01
# 
# * '*': p <= 0.01
# 
# * '**': p <= 0.001
# 
# * '***': p <= 0.0001
# 
# * '****': p <= 0.00001

# plot
# temp_pdf_function <-
#   function(x) {
#     pdf(
#       file = (x), 
#       width = 14, 
#       height = 8, 
#       useDingbats = F
#     )
#   }

temp_png_function <-
  function(x) {
    png(
      file = (x),
      width = 12,
      height = 6,
      res = 300,
      units = 'in'
    )
  }

for(plottype in c("Genes", "UMIs")) {
  if(plottype == "Genes") {colname <- "downsampled_nFeature_RNA"}
  if(plottype == "UMIs") {colname <- "downsampled_nCount_RNA"}
  if(plottype == "Mitocondrial") {colname <- "downsampled_percent.mito"}
  
  temp_gene_data_df <- NULL
  for(sampleID in temp_sampleIDs[!temp_sampleIDs %in% "CID4523"]) {
    
    temp_seurat_object <- 
      get(paste0("seurat_10X_downsampled_", sampleID))
    
    temp_df <- data.frame(temp_seurat_object@meta.data[,colname],
                          temp_seurat_object@meta.data[,"int_PCA_B_res.0.8"],
                          temp_seurat_object@meta.data[,"orig.ident"],
                          temp_seurat_object@meta.data[,"condition"])
    
    colnames(temp_df) <- c(colname, "clusterID", "orig.ident", "condition")
    
    temp_gene_data_df <- rbind(temp_gene_data_df,
                               temp_df)
  }
  
  # patient
  temp_gene_data_df$patient <- NA
  for(i in unique(temp_gene_data_df$orig.ident)) {
    if(i %in% c("CID4471", "CID4471CELLSUS", "CID4471CHUNKS")){temp_pid = "BC P1"}
    if(i %in% c("CID44971", "CID44971CELLSUS", "CID44971CHUNKS")){temp_pid = "BC P2"}
    if(i %in% c("CID4513", "CID4513CELLSUS", "CID4513CHUNKS")){temp_pid = "BC P3"}
    # if(i %in% c("CID4523", "CID4523CELLSUS", "CCID4523CHUNKS")){temp_pid = "BC P4"}
    if(i %in% c("PID17267", "PID17267CELLSUS", "PID17267CHUNKS")){temp_pid = "PC P1"}
    if(i %in% c("PID20033", "PID20033CELLSUS", "PID20033CHUNKS")){temp_pid = "PC P2"}
    if(i %in% c("SCC180161", "SCC180161SD", "SCC180161ND")){temp_pid = "M P1"}
    temp_gene_data_df$patient[temp_gene_data_df$orig.ident == paste0(i)] <- temp_pid
  }
  
  temp_gene_data_df$condition <- factor(temp_gene_data_df$condition,
                                        levels = c("Fresh Tissue", 
                                                   "Cryopreserved Cell Suspension", 
                                                   "Cryopreserved Tissue",
                                                   "Cryopreserved Overnight"))
  
  temp_gene_data_df$patient <- factor(temp_gene_data_df$patient,
                                      levels = c("BC P1",
                                                 "BC P2",
                                                 "BC P3",
                                                 # "BC P4",
                                                 "PC P1",
                                                 "PC P2",
                                                 "M P1"
                                      ))
  
  # filter out crap samples
  # temp_gene_data_df <- temp_gene_data_df[! temp_gene_data_df$condition == "CID4523CHUNKS",]
  temp_gene_data_df <- temp_gene_data_df[! temp_gene_data_df$orig.ident == "PID20033CELLSUS",]
  
  
  # unique ID
  temp_gene_data_df$ID <- paste0(temp_gene_data_df$patient,
                                 temp_gene_data_df$condition)
  temp_gene_data_df$ID <- factor(temp_gene_data_df$ID,
                                 levels = unique(temp_gene_data_df$ID))
  
  # colours 
  # temp_colours <- c("#1b9e77", "#d95f02", "#abd9e9", "#7570b3")
  
  temp_colours <- c("#1b9e77", "#d95f02", "#7570b3", "#abd9e9")
  
  if(plottype == "Genes"){
    temp_num <- 12000
    temp_factor_to_align_stat_bars <- 0.8
    temp_scale <-  c(0,500,1000,2500,5000,12500)
    temp_label <- "Number of Genes"
  } 
  if(plottype == "UMIs"){
    temp_num <- 120000
    temp_factor_to_align_stat_bars <- 0.8
    temp_scale <- c(0,1000,5000,10000,25000,125000)
    temp_label <- "Number of UMIs"
  }
  if(plottype == "Mitocondrial"){
    temp_num <- 25
    temp_factor_to_align_stat_bars <- 0.9
    temp_scale <- c(0,5,10,20,25)
  }
  
  temp_label_y <- c(temp_num, temp_num*temp_factor_to_align_stat_bars,
                    temp_num, temp_num*temp_factor_to_align_stat_bars,
                    temp_num, temp_num*temp_factor_to_align_stat_bars,
                    temp_num, temp_num*temp_factor_to_align_stat_bars,
                    temp_num,
                    temp_num, temp_num*temp_factor_to_align_stat_bars)
  
  
  
  for(i in unique(temp_gene_data_df$patient)) {
    
    if(i == "PC P2"){
      my_comparisons <- list(c("Fresh Tissue", "Cryopreserved Tissue")
      )
    } 
    if(i == "M P1"){
      my_comparisons <- list(c("Fresh Tissue", "Cryopreserved Cell Suspension"),
                             c("Fresh Tissue", "Cryopreserved Overnight")
      )
    } 
    if(!i == "M P1" && !i == "PC P2"){
      my_comparisons <- list(c("Fresh Tissue", "Cryopreserved Cell Suspension"),
                             c("Fresh Tissue", "Cryopreserved Tissue")
      )
    }
    
    
    
    temp_gene_data_df_subset <- temp_gene_data_df[temp_gene_data_df$patient == i,,drop=F]
    temp_gene_data_df_subset[,colname] <- as.numeric(temp_gene_data_df_subset[,colname])
    temp_gene_data_df_subset$clusterID <- paste0("c",temp_gene_data_df_subset$clusterID)
    
    
    temp_gene_data_df_subset$clusterID <- factor(temp_gene_data_df_subset$clusterID,
                                                 levels = str_sort(unique(temp_gene_data_df_subset$clusterID),
                                                                   numeric=T)
    )
    temp_gene_data_df_subset <- temp_gene_data_df_subset[order(temp_gene_data_df_subset$clusterID),,drop=F]
    
    # unique ID
    temp_gene_data_df_subset$ID <- paste0(temp_gene_data_df_subset$clusterID,
                                          temp_gene_data_df_subset$condition)
    temp_gene_data_df_subset$ID <- factor(temp_gene_data_df_subset$ID,
                                          levels = unique(temp_gene_data_df_subset$ID))
    
    # my_comparisons <- list()
    # for(clusterID in unique(temp_gene_data_df_subset$clusterID)){
    #   temp_gene_data_df_subset_cluster <- temp_gene_data_df_subset[temp_gene_data_df_subset$clusterID == clusterID,,drop=F]
    #   for(number in c(2:length(unique(temp_gene_data_df_subset_cluster$ID)))){
    #     temp_vector <- list(c(as.character(unique(temp_gene_data_df_subset_cluster$ID)[1]),
    #                      as.character(unique(temp_gene_data_df_subset_cluster$ID)[number])
    #                      ))
    # 
    #     my_comparisons <- c(my_comparisons, temp_vector)
    # 
    #   }
    # 
    # }
    # plot graph with p vals
    if(!plottype == "Mitocondrial"){
      temp_gene_boxplot <- ggboxplot(temp_gene_data_df_subset, 
                                     x = "condition",
                                     y= colname,
                                     fill = "condition",
                                     # group = "patient",
                                     # label="patient",
                                     # legend = "none",
                                     facet.by = "clusterID",
                                     repel = F,
                                     nrow=1,
                                     size=1, 
                                     # scales ="free_x",
                                     outlier.size = 0.01
                                     # ylim=c(0,100000)
      )+
        stat_compare_means(method = "t.test",
                           paired = F,
                           label = "p.signif", 
                           # ref.group = "Fresh Tissue",
                           hide.ns = F,
                           comparisons = my_comparisons,
                           bracket.size = 0.1,
                           label.y = temp_label_y,
                           symnum.args = symnum.args) +  # Pairwise comparison against all
        # stat_compare_means(method = "anova", 
        #                    label.y = temp_max_label_y) +
        xlab(paste0(i, " Cluster")) + ylab(temp_label) + theme(axis.text.x = element_blank()) +
        scale_fill_manual(values = temp_colours) +
        coord_trans(y = "log") +
        scale_y_continuous(breaks=temp_scale,
                           labels=temp_scale)
      
      temp_png_function(paste0("SUPP_03_nGenes_nUMIs/ggboxplot_sample_",i,"_", plottype,".png"))
      print(temp_gene_boxplot)
      dev.off()
      
      # temp_pdf_function(paste0("SUPP_03_nGenes_nUMIs/ggboxplot_sample_",i,"_", plottype,".pdf"))
      # print(temp_gene_boxplot)
      # dev.off()
      
    }
    
  }
}


# print metrics
for(plottype in c("Genes", "UMIs")) {
  print(plottype)
  if(plottype == "Genes") {colname <- "downsampled_nFeature_RNA"}
  if(plottype == "UMIs") {colname <- "downsampled_nCount_RNA"}
  
  temp_gene_data_df <- NULL
  temp_gene_data_df_summary <- NULL
  for(sampleID in temp_sampleIDs) {
    print(sampleID)
    temp_seurat_object <- 
      get(paste0("seurat_10X_downsampled_", sampleID))
    
    temp_df <- data.frame(temp_seurat_object@meta.data[,colname],
                          temp_seurat_object@meta.data[,"int_PCA_B_res.0.8"],
                          temp_seurat_object@meta.data[,"orig.ident"],
                          temp_seurat_object@meta.data[,"condition"])
    
    colnames(temp_df) <- c(colname, "clusterID", "orig.ident", "condition")
    
    temp_conditions <- as.vector(unique(temp_df$condition))
    temp_clusters <-  as.vector(unique(temp_df$clusterID))
    temp_clusters <- str_sort(temp_clusters,numeric=T)
    for(clusterID in temp_clusters){
      print(clusterID)
      temp_df_subset_subset <- temp_df[temp_df$clusterID == clusterID,]
      
      for(condition in temp_conditions) {
        print(condition)
        temp_df_subset_subset_cond <- temp_df_subset_subset[temp_df_subset_subset$condition == condition,]
        temp_df_summary_subset <- data.frame(sample = sampleID,
                                             condition = condition,
                                             cluster = clusterID,
                                             mean = mean(temp_df_subset_subset_cond[,colname])
        )
        
        temp_gene_data_df_summary <- rbind(temp_gene_data_df_summary, 
                                           temp_df_summary_subset)
      }
    }
    
    temp_gene_data_df <- rbind(temp_gene_data_df,
                               temp_df)
    
  }
  
  # patient
  temp_gene_data_df$patient <- NA
  for(i in unique(temp_gene_data_df$orig.ident)) {
    if(i %in% c("CID4471", "CID4471CELLSUS", "CID4471CHUNKS")){temp_pid = "BC P1"}
    if(i %in% c("CID44971", "CID44971CELLSUS", "CID44971CHUNKS")){temp_pid = "BC P2"}
    if(i %in% c("CID4513", "CID4513CELLSUS", "CID4513CHUNKS")){temp_pid = "BC P3"}
    if(i %in% c("CID4523", "CID4523CELLSUS", "CCID4523CHUNKS")){temp_pid = "BC P4"}
    if(i %in% c("PID17267", "PID17267CELLSUS", "PID17267CHUNKS")){temp_pid = "PC P1"}
    if(i %in% c("PID20033", "PID20033CELLSUS", "PID20033CHUNKS")){temp_pid = "PC P2"}
    if(i %in% c("SCC180161", "SCC180161SD", "SCC180161ND")){temp_pid = "M P1"}
    temp_gene_data_df$patient[temp_gene_data_df$orig.ident == paste0(i)] <- temp_pid
  }
  
  
  write.csv(temp_gene_data_df_summary,
            paste0("SUPP_03_nGenes_nUMIs/downsampled_mean_",plottype, ".csv"))
}



# 06 BULK CORRELATIONS DOWNSAMPLED ----------------------------------------------------

dir.create("05_Bulk_correlations_downsampled/")

# for(sampleID in temp_sampleIDs[!temp_sampleIDs %in% "PID20033"]) {
# for(sampleID in "PID20033") {
for(sampleID in temp_sampleIDs) {
  
  temp_seurat_object <- 
    get(paste0("seurat_10X_downsampledobject_", sampleID))
  
  Idents(temp_seurat_object) <- "condition"
  
  temp_cluster.averages <- 
    AverageExpression(temp_seurat_object,
                      assays = "RNA",
                      return.seurat = TRUE)
  
  temp_data_frame <- GetAssayData(object = temp_cluster.averages,
                                  assay = "RNA",
                                  slot = "data")
  
  temp_data_frame <- 
    as.data.frame(temp_data_frame)
  
  temp_Names <- colnames(temp_data_frame)
  
  if(sampleID == "PID20033"){
    temp_vector <- c(2)
  } else { temp_vector <- c(2,3) 
  }
  # get linear regression summary
  for(i in temp_vector) {
    temp_linearMod <- lm(get(temp_Names[1]) ~ get(temp_Names[i]),
                         data = temp_data_frame)
    temp <- summary(temp_linearMod)
    temp_Rsquared <- paste0("R = ", round(temp$adj.r.squared, digits=3))
    
    n <- paste0("temp_Rsquared_",sampleID, "_", i)
    assign(n, temp_Rsquared)
  }
  
  # temp_data_frame$gene <- row.names(temp_data_frame)
  # row.names(temp_data_frame) <- NULL
  # 
  # temp_data_frame_melted <- melt(temp_data_frame, 
  #                                id="gene")
  
  n <- paste0("temp_data_frame_",sampleID)
  assign(n, temp_data_frame)
}


# plot
# temp_pdf_function <-
#   function(x) {
#     pdf(
#       file = (x), 
#       width = 6, 
#       height = 6, 
#       useDingbats = F
#     )
#   }

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


# for(sampleID in temp_sampleIDs[!temp_sampleIDs %in% "PID20033"]) {
for(sampleID in "PID20033") {
  
  temp_data_frame <- 
    get(paste0("temp_data_frame_",sampleID))
  temp_Rsquared_2 <- get(paste0("temp_Rsquared_",sampleID, "_", 2))
  
  if(!sampleID == "PID20033"){
    temp_Rsquared_3 <- get(paste0("temp_Rsquared_",sampleID, "_", 3))
  }
  
  temp_label <- sampleID
  temp_Names <- c("FT", "CCS", "CT")
  if(sampleID == "SCC180161") {
    temp_Names <- c("FT", "CT", "CO")
  }
  if(sampleID == "PID20033") {
    temp_Names <- c("FT", "CT")
  }
  # overlay scatter plots of fresh vs cellsus & fresh vs chunks
  if(! sampleID == "PID20033"){
    temp_png_function(paste0("05_Bulk_correlations_downsampled/",sampleID,"_scatter_fresh_vs_cellsus_chunks_overlayed.png"))
    temp <- plot(x=temp_data_frame[,1], 
                 y=temp_data_frame[,3], 
                 main= temp_label, 
                 xlab = temp_Names[1],
                 ylab = paste0(temp_Names[2],
                               "(",temp_Rsquared_2,") ",
                               "/", 
                               temp_Names[3],
                               "(",temp_Rsquared_3,") "
                 ), 
                 col = "#d73027")
    points(x=temp_data_frame[,1], 
           y=temp_data_frame[,2], 
           main= temp_label, 
           xlab = temp_Names[1],
           ylab = paste0(temp_Names[2],
                         "(",temp_Rsquared_2,") ",
                         "/", 
                         temp_Names[3],
                         "(",temp_Rsquared_3,") "
           ), 
           col = "#4575b4", 
           pch = 6)
    lines(loess.smooth(x=temp_data_frame[,1], 
                       y=temp_data_frame[,2], 
                       main= temp_label, 
                       xlab = temp_Names[1],
                       ylab = temp_Names[2]),
          col = "#d73027",
          lwd = 4)
    lines(loess.smooth(x=temp_data_frame[,1], 
                       y=temp_data_frame[,3], 
                       main= temp_label, 
                       xlab = temp_Names[1],
                       ylab = temp_Names[3]),
          col = "#4575b4",
          lwd = 4)
    dev.off()
    
    temp_png_function(paste0("05_Bulk_correlations_downsampled/nolabel_",sampleID,"_scatter_fresh_vs_cellsus_chunks_overlayed.png"))
    temp <- plot(x=temp_data_frame[,1], 
                 y=temp_data_frame[,3], 
                 main= " ", 
                 xlab = " ",
                 ylab = " ", 
                 col = "#d73027")
    points(x=temp_data_frame[,1], 
           y=temp_data_frame[,2], 
           main= " ", 
           xlab = " ",
           ylab = " ", 
           col = "#4575b4", 
           pch = 6)
    lines(loess.smooth(x=temp_data_frame[,1], 
                       y=temp_data_frame[,2], 
                       main= temp_label, 
                       xlab = temp_Names[1],
                       ylab = temp_Names[2]),
          col = "#d73027",
          lwd = 4)
    lines(loess.smooth(x=temp_data_frame[,1], 
                       y=temp_data_frame[,3], 
                       main= temp_label, 
                       xlab = temp_Names[1],
                       ylab = temp_Names[3]),
          col = "#4575b4",
          lwd = 4)
    dev.off()
  }
  
  if(sampleID == "PID20033"){
    temp_png_function(paste0("05_Bulk_correlations_downsampled/nolabel_",sampleID,"_scatter_fresh_vs_cellsus_chunks_overlayed.png"))
    temp <- plot(x=temp_data_frame[,1],
                 y=temp_data_frame[,2],
                 main= " ",
                 xlab = " ",
                 ylab = " ",
                 col = "#4575b4")
    lines(loess.smooth(x=temp_data_frame[,1],
                       y=temp_data_frame[,2],
                       main= " ",
                       xlab = " ",
                       ylab = " "),
          col = "#4575b4",
          lwd = 4)
    dev.off()
    
    temp_png_function(paste0("05_Bulk_correlations_downsampled/",sampleID,"_scatter_fresh_vs_cellsus_chunks_overlayed.png"))
    temp <- plot(x=temp_data_frame[,1],
                 y=temp_data_frame[,2],
                 main= temp_label,
                 xlab = temp_Names[1],
                 ylab = paste0(temp_Names[2],
                               "(",temp_Rsquared_2,") "
                 ),
                 col = "#4575b4")
    lines(loess.smooth(x=temp_data_frame[,1],
                       y=temp_data_frame[,2],
                       main= temp_label,
                       xlab = temp_Names[1],
                       ylab = temp_Names[2]),
          col = "#4575b4",
          lwd = 4)
    dev.off()
    
  }
}





# 07 CLUSTER LEVEL CORRELATIONS FILTERED FOR CLUSTERS WITH ALL CELLS -------------------------------------------

dir.create("06_cluster_level_correlations/")

# generate dataframes for cluster level correlations
temp_rsquared_df <- NULL
for(sampleID in temp_sampleIDs) {
  print(sampleID)
  temp_seurat_object <-
    get(paste0("seurat_10X_downsampledobject_", sampleID))
  
  Idents(temp_seurat_object) <- "int_PCA_B_res.0.8"
  
  temp_cluster.averages <-
    AverageExpression(temp_seurat_object,
                      assays = "RNA",
                      return.seurat = TRUE,
                      add.ident = "condition")
  
  temp_data_frame <- GetAssayData(object = temp_cluster.averages,
                                  assay = "RNA",
                                  slot = "data")
  
  temp_data_frame <-
    as.data.frame(temp_data_frame)
  
  temp_Names <- colnames(temp_data_frame)
  
  # filter out all small noisey clusters
  temp_cluster_nums <- data.frame(table(temp_seurat_object@meta.data$int_PCA_B_res.0.8))
  temp_cluster_nums <- temp_cluster_nums[temp_cluster_nums$Freq > 100,]
  temp_cluster_nums <- as.vector(unique(temp_cluster_nums$Var1))
  
  # filter out all clusters that have 0 cells from any given condition
  temp_df_combined <- NULL
  for(i in unique(temp_seurat_object@meta.data$int_PCA_B_res.0.8)){
    # print(i)
    
    temp_df <-
      temp_seurat_object@meta.data[temp_seurat_object@meta.data$int_PCA_B_res.0.8 %in% i ,]
    
    temp_df <-
      as.data.frame((table(temp_df$condition)))
    
    temp_df$cluster <- i
    
    temp_df_combined <- rbind(temp_df_combined,
                              temp_df)
  }
  
  temp_df_combined <- arrange(temp_df_combined,
                              Freq)
  
  temp_df_zeros <- temp_df_combined[temp_df_combined$Freq < 19 ,]
  temp_zero_clusters <- unique(temp_df_zeros$cluster)
  
  # cluster IDs
  temp_cluster_ids <-
    temp_cluster_nums[! temp_cluster_nums %in% temp_zero_clusters]
  
  temp_conditions <-
    unique(temp_seurat_object@meta.data$condition)
  
  # get linear regression summary
  for(cluster in temp_cluster_ids) {
    print(cluster)
    temp_Names_subset <- c(paste0(cluster, "_", temp_conditions[1]),
                           paste0(cluster, "_", temp_conditions[2]),
                           paste0(cluster, "_", temp_conditions[3]))
    
    for(i in c(2,3)) {
      
      temp_linearMod <- lm(get(temp_Names_subset[1]) ~ get(temp_Names_subset[i]),
                           data = temp_data_frame)
      temp <- summary(temp_linearMod)
      # temp_Rsquared <- paste0("R = ", round(temp$adj.r.squared, digits=3))
      temp_Rsquared <- round(temp$adj.r.squared, digits=3)
      # print(temp_Rsquared)
      
      # n <- paste0("temp_Rsquared_",sampleID, "_cluster_", cluster, "_", i)
      # assign(n, temp_Rsquared)
      
      # cell type
      temp_df <- temp_seurat_object@meta.data[temp_seurat_object@meta.data$int_PCA_B_res.0.8 == cluster,]
      cluster_celltype <- unique(temp_df$int_cluster_major_PC_B_res.0.8)
      print(cluster_celltype)
      
      temp_rsquared_df_subset <- data.frame(sampleID=sampleID,
                                            cluster=cluster,
                                            condition=paste0("FT vs ", temp_conditions[i]),
                                            Rsquared=temp_Rsquared,
                                            celltype=cluster_celltype,
                                            cellnum = nrow(temp_df))
      
      temp_rsquared_df <- rbind(temp_rsquared_df,
                                temp_rsquared_df_subset)
      
    }
    
    
  }
  
  m <- paste0("temp_cluster_ids_",sampleID)
  assign(m,temp_cluster_ids)
  
  n <- paste0("temp_data_frame_bycluster_",sampleID)
  assign(n, temp_data_frame)
  
  # g <- paste0("temp_rsquared_df_",sampleID)
  # assign(n, temp_rsquared_df)
}

# temp_rsquared_df <- read.csv("../20191001_paper_figures/06_cluster_level_correlations/correlations_per_cluster_per_condition_filtered.csv", row.names = "X")

temp_rsquared_df$cluster <- paste0("c",temp_rsquared_df$cluster)
temp_rsquared_df$cluster <- factor(temp_rsquared_df$cluster,
                                   levels = str_sort(unique(temp_rsquared_df$cluster),
                                                     numeric=T)
)

temp_rsquared_df$condition <- factor(temp_rsquared_df$condition,
                                     levels=c("FT vs Cryopreserved Cell Suspension",
                                              "FT vs Cryopreserved Tissue",
                                              "FT vs Cryopreserved Overnight"))



# patient
temp_rsquared_df$patient <- NA
for(i in unique(temp_rsquared_df$sampleID)) {
  if(i %in% c("CID4471", "CID4471CELLSUS", "CID4471CHUNKS")){temp_pid = "BC P1"}
  if(i %in% c("CID44971", "CID44971CELLSUS", "CID44971CHUNKS")){temp_pid = "BC P2"}
  if(i %in% c("CID4513", "CID4513CELLSUS", "CID4513CHUNKS")){temp_pid = "BC P3"}
  # if(i %in% c("CID4523", "CID4523CELLSUS", "CCID4523CHUNKS")){temp_pid = "BC P4"}
  if(i %in% c("PID17267", "PID17267CELLSUS", "PID17267CHUNKS")){temp_pid = "PC P1"}
  if(i %in% c("PID20033", "PID20033CELLSUS", "PID20033CHUNKS")){temp_pid = "PC P2"}
  if(i %in% c("SCC180161", "SCC180161SD", "SCC180161ND")){temp_pid = "M P1"}
  temp_rsquared_df$patient[temp_rsquared_df$sampleID == paste0(i)] <- temp_pid
}
temp_rsquared_df$patient <- factor(temp_rsquared_df$patient,
                                   levels = c("BC P1",
                                              "BC P2",
                                              "BC P3",
                                              "PC P1",
                                              "M P1"
                                   ))


temp_rsquared_df_ordered <- 
  temp_rsquared_df[
    with(temp_rsquared_df, order(patient, cluster)),
    ]

write.csv(temp_rsquared_df_ordered,
          "06_cluster_level_correlations/correlations_per_cluster_per_condition_filtered.csv")


# plot
library(ggrepel)
temp_pdf_function <-
  function(x) {
    pdf(
      file = (x), 
      width = 20, 
      height = 6, 
      useDingbats = F
    )
  }

temp_ggplot <- ggplot(temp_rsquared_df, 
                      aes(x=condition, 
                          y=Rsquared, 
                          group=cluster)) +
  geom_line(aes(color=cluster),
            linetype="dashed")+
  geom_point(aes(color=cluster,
                 shape=condition),
             size=3)+
  facet_wrap(~patient, 
             scales="free_x", 
             nrow=1) +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        strip.text.x = element_text(size=20),
        legend.key = element_blank()) +
  scale_y_continuous(breaks = seq(0,1,0.1)) +
  geom_text_repel(aes(label=ifelse(Rsquared<0.9,
                                   as.character(cluster),'')
  ),
  hjust=0,
  vjust=0,
  size = 5,
  force=2,
  direction="y",
  nudge_x=0.25
  ) + # label all clusters with Rsquared less than 0.9 
  geom_hline(yintercept=0.9, 
             color = "red", 
             size=0.5) +
  xlab("Condition Comparisons")

temp_pdf_function(paste0("06_cluster_level_correlations/Correlations_by_cluster_correlations_by_cluster.pdf"))
print(temp_ggplot)
dev.off()





