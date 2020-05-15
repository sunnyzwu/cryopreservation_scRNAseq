# KEBT PLOTTING
#
#
# SETUP -------------------------------------------------------------------

temp_sampleIDs <- c("CID4471", "CID44971", "CID4513", "PID17267", "PID20033", "SCC180161")

library(Seurat)
library(ggplot2)
library(cowplot)
library(reshape2)
library(stringr)

# DIRECTORY -----------------------------------------------------------------

setwd("/share/ScratchGeneral/sunwu/projects/cryopreservation_paper/seurat_v3/03_output/")
dir.create("Results_summary_plots_v5")
setwd("Results_summary_plots_v5")

# dir.create("kBET_FIGURES")
# dir.create("kBET_FIGURES/PNGs/")
# dir.create("kBET_FIGURES/PDFs/")
# dir.create("kBET_FIGURES/SUMMARIES/")

# LOAD DATA ---------------------------------------------------------------

# BULK DATA RUNs
temp_wd <- "/share/ScratchGeneral/sunwu/projects/cryopreservation_paper/seurat_v3/03_output/v5_with_silhouette/run04_seuratmetrics/kBET_output/"
cellnum <- 250

plot_data_combined_allcells <- NULL
for(sampleID in temp_sampleIDs[!temp_sampleIDs %in% c("CID4523")]) {
  print(sampleID)
  {temp_condition_ids <- c("FT", "CCS", "CT") 
    
    if(sampleID == "SCC180161") {
      temp_condition_ids <- c("FT", "CT", "CO") } 
    
    if(sampleID == "PID20033") {
      temp_condition_ids <- c("FT", "CT") }
    
    if(sampleID == "CID4523") {
      temp_condition_ids <- c("FT", "CCS") }
    
  }
  
  
  
  temp_csv <- read.csv(paste0(temp_wd, 
                              sampleID,
                              "_", cellnum, 
                              "/08_silhouette_mixing_LSC_metrics.csv"),
                       row.names = "X")
  temp_csv$sample <- sampleID
  
  
  
  temp_csv$patient <- NA
    if(sampleID %in% c("CID4471")){temp_pid = "BC P1"}
    if(sampleID %in% c("CID44971")){temp_pid = "BC P2"}
    if(sampleID %in% c("CID4513")){temp_pid = "BC P3"}
    if(sampleID %in% c("PID17267")){temp_pid = "PC P1"}
    if(sampleID %in% c("PID20033")){temp_pid = "PC P2"}
    if(sampleID %in% c("SCC180161")){temp_pid = "M P1"}
    temp_csv$patient[temp_csv$sample == paste0(sampleID)] <- temp_pid
  
  
  plot_data_combined_allcells <- rbind(plot_data_combined_allcells, 
                                           temp_csv)
    }

  # factorise
# plot_data_combined_allcells$ID <- factor(plot_data_combined_allcells$ID,
#                                          levels=unique(plot_data_combined_allcells$ID))


# PLOT PER CELL TYPE ----------------------------------------------------

# plot
temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 8, 
      height = 5, 
      res = 300, 
      units = 'in'
    )
  }
# temp_pdf_function <-
#   function(x) {
#     pdf(
#       file = (x), 
#       width = 6, 
#       height = 4, 
#       useDingbats = F
#     )
#   }
temp_colours <- c("#1b9e77", "#d95f02",  "#7570b3", "#abd9e9")

plot_data_combined_allcells$condition <- factor(plot_data_combined_allcells$condition, 
                            levels=unique(plot_data_combined_allcells$condition))


for(metrics in c("silhouette","mixing.metric","local.struct")){
  
  
  print(metrics)
  temp_df_subset <- plot_data_combined_allcells[,c( metrics,"condition", "cluster", "patient")]
  colnames(temp_df_subset) <- c("value","condition", "cluster", "patient")
  
  if(metrics == "mixing.metric"){
    temp_range <- c(0,300)
  }
  if(metrics == "local.struct"){
    temp_range <- c(0,1)
  }
  if(metrics == "silhouette"){
    temp_range <- c(-1,1)  
    }

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
    scale_y_continuous(limits=temp_range) +
    theme(axis.text.x = element_text(angle=25, vjust = .25, size=8)) +
    scale_fill_manual(values = temp_colours) +
    facet_wrap(. ~ patient,
               scales="free_x",
               nrow =1) # free x removes white space
  
  
  temp_png_function(paste0("01_",metrics,".png"))
  print(temp_ggplot_clusters)
  dev.off()
  
  if(metrics == "silhouette"){
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
      # scale_y_continuous(limits=temp_range) +
      theme(axis.text.x = element_text(angle=25, vjust = .25, size=8)) +
      scale_fill_manual(values = temp_colours) +
      facet_wrap(. ~ patient,
                 scales="free_x",
                 nrow =1) # free x removes white space
    
    
    temp_png_function(paste0("01_",metrics,"_rescale.png"))
    print(temp_ggplot_clusters)
    dev.off()
    
    
  }
  
}


# PLOT FOR FIGURES BY CANCER --------------------------------------------------------

temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 10, 
      height = 4, 
      res = 300, 
      units = 'in'
    )
  }

temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 10,
      height = 4,
      useDingbats=F
    )
  }


for(cancertype in c("breast", "prostateMEL")) {
  print(cancertype)
  if(cancertype == "breast"){
    temp_pids <- c("BC P1", "BC P2", "BC P3")
  }
  if(cancertype == "prostateMEL"){
    temp_pids <- c("PC P1", "PC P2", "M P1")
  }
  plot_data_combined_allcells_subset <- 
    plot_data_combined_allcells[plot_data_combined_allcells$patient %in% temp_pids,]
  
  
  
  for(metrics in c("silhouette","mixing.metric","local.struct")){
    
    print(metrics)
    temp_df_subset <- plot_data_combined_allcells_subset[,c( metrics,"condition", "cluster", "patient")]
    colnames(temp_df_subset) <- c("value","condition", "cluster", "patient")
    
    temp_df_subset$patient <- factor(temp_df_subset$patient,
                                     levels=temp_pids)
    
    if(metrics == "mixing.metric"){
      temp_range <- c(0,300)
      
      temp_ggplot_clusters <- ggplot(temp_df_subset, 
                                     aes(condition, 
                                         value,
                                         label = condition,
                                         fill = condition)) + 
        geom_boxplot(outlier.shape = NA, 
                     size=0.3, 
                     alpha = 0.5) + 
        labs(y=paste0(metrics, ' Score'),
             title=metrics) +
        theme_bw() +
        scale_y_continuous(trans = "log",
                           breaks = seq(0,300,50),
                           labels = seq(0,300,50)) +
        theme(axis.text.x = element_blank(),
              # legend.position = "none",
              plot.margin = unit(c(0,1,0,0), "cm")) +
        scale_fill_manual(values = temp_colours) +
        facet_wrap(. ~ patient,
                   scales="free_x",
                   nrow =1) # free x removes white space
      
    }
    if(metrics == "local.struct"){
      temp_range <- c(0,1)
      
      temp_ggplot_clusters <- ggplot(temp_df_subset, 
                                     aes(condition, 
                                         value,
                                         label = condition,
                                         fill = condition)) + 
        geom_boxplot(outlier.shape = NA, 
                     size=0.3, 
                     alpha = 0.5) + 
        labs(y=paste0(metrics, ' Score'),
             title=metrics) +
        theme_bw() +
        scale_y_continuous(limits=temp_range) +
        theme(axis.text.x = element_blank(),
              plot.margin = unit(c(0,1,0,0), "cm")) +
        scale_fill_manual(values = temp_colours) +
        facet_wrap(. ~ patient,
                   scales="free_x",
                   nrow =1) # free x removes white space
      
    }
    if(metrics == "silhouette"){
      temp_range <- c(-1,1) 
      
      temp_ggplot_clusters <- ggplot(temp_df_subset, 
                                     aes(condition, 
                                         value,
                                         label = condition,
                                         fill = condition)) + 
        geom_boxplot(outlier.shape = NA, 
                     size=0.3, 
                     alpha = 0.5) + 
        labs(y=paste0(metrics, ' Score'),
             title=metrics) +
        theme_bw() +
        scale_y_continuous(limits=temp_range) +
        theme(axis.text.x = element_blank(),
              # legend.position = "none",
              plot.margin = unit(c(0,1,0,0), "cm") # T, r, b, l
              ) +
        scale_fill_manual(values = temp_colours) +
        facet_wrap(. ~ patient,
                   scales="free_x",
                   nrow =1) # free x removes white space
      
    }
    
    
    n <- paste0("temp_ggplot_clusters_",metrics)
    assign(n,
           temp_ggplot_clusters)
  }
  
  temp_grid <- 
    plot_grid(plotlist=mget(paste0("temp_ggplot_clusters_",c("silhouette","mixing.metric","local.struct"))),
              nrow = 1)
  
  temp_png_function(paste0("02_",cancertype,".png"))
  print(temp_grid)
  dev.off()
  
  temp_pdf_function(paste0("02_",cancertype,".pdf"))
  print(temp_grid)
  dev.off()
  
  }
  
  
  
  



# PLOT FOR FIGURES ALL --------------------------------------------------------

temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 13, 
      height = 4, 
      res = 300, 
      units = 'in'
    )
  }

temp_pdf_function <-
  function(x) {
    pdf(
      file = (x),
      width = 13,
      height = 4,
      useDingbats=F
    )
  }


# for(cancertype in c("breast", "prostateMEL")) {
#   print(cancertype)
#   if(cancertype == "breast"){
#     temp_pids <- c("BC P1", "BC P2", "BC P3")
#   }
#   if(cancertype == "prostateMEL"){
#     temp_pids <- c("PC P1", "PC P2", "M P1")
#   }
#   plot_data_combined_allcells_subset <- 
#     plot_data_combined_allcells[plot_data_combined_allcells$patient %in% temp_pids,]
#   
temp_strip_text_x <- 10
  for(metrics in c("silhouette","mixing.metric","local.struct")){
    
    print(metrics)
    temp_df_subset <- plot_data_combined_allcells[,c( metrics,"condition", "cluster", "patient")]
    colnames(temp_df_subset) <- c("value","condition", "cluster", "patient")
    
    temp_df_subset$patient <- factor(temp_df_subset$patient,
                                     levels=unique(temp_df_subset$patient))
    
    if(metrics == "mixing.metric"){
      temp_range <- c(0,300)
      
      temp_ggplot_clusters <- ggplot(temp_df_subset, 
                                     aes(condition, 
                                         value,
                                         label = condition,
                                         fill = condition)) + 
        geom_boxplot(outlier.shape = NA, 
                     size=0.3, 
                     alpha = 0.5) + 
        labs(y=paste0(metrics, ' Score'),
             title=metrics) +
        theme_bw() +
        scale_y_continuous(trans = "log",
                           breaks = seq(0,300,50),
                           labels = seq(0,300,50)) +
        theme(axis.text.x = element_blank(),
              legend.position = "none",
              plot.margin = unit(c(0,2,0,0), "cm"), 
              strip.text.x = element_text(size=temp_strip_text_x)) +
        scale_fill_manual(values = temp_colours) +
        facet_wrap(. ~ patient,
                   scales="free_x",
                   nrow =1) # free x removes white space
      
    }
    if(metrics == "local.struct"){
      temp_range <- c(0,1)
      
      temp_ggplot_clusters <- ggplot(temp_df_subset, 
                                     aes(condition, 
                                         value,
                                         label = condition,
                                         fill = condition)) + 
        geom_boxplot(outlier.shape = NA, 
                     size=0.3, 
                     alpha = 0.5) + 
        labs(y=paste0(metrics, ' Score'),
             title=metrics) +
        theme_bw() +
        scale_y_continuous(limits=temp_range) +
        theme(axis.text.x = element_blank(),
              plot.margin = unit(c(0,0,0,0), "cm"), 
              strip.text.x = element_text(size=temp_strip_text_x)) +
        scale_fill_manual(values = temp_colours) +
        facet_wrap(. ~ patient,
                   scales="free_x",
                   nrow =1) # free x removes white space
      
    }
    if(metrics == "silhouette"){
      temp_range <- c(-1,1) 
      
      temp_ggplot_clusters <- ggplot(temp_df_subset, 
                                     aes(condition, 
                                         value,
                                         label = condition,
                                         fill = condition)) + 
        geom_boxplot(outlier.shape = NA, 
                     size=0.3, 
                     alpha = 0.5) + 
        labs(y=paste0(metrics, ' Score'),
             title=metrics) +
        theme_bw() +
        scale_y_continuous(limits=temp_range) +
        theme(axis.text.x = element_blank(),
              legend.position = "none",
              strip.text.x = element_text(size=temp_strip_text_x),
              plot.margin = unit(c(0,2,0,0), "cm") # T, r, b, l
        ) +
        scale_fill_manual(values = temp_colours) +
        facet_wrap(. ~ patient,
                   scales="free_x",
                   nrow =1) # free x removes white space
      
    }
    
    
    n <- paste0("temp_ggplot_clusters_",metrics)
    assign(n,
           temp_ggplot_clusters)
  }
  
  temp_grid <- 
    plot_grid(plotlist=mget(paste0("temp_ggplot_clusters_",c("silhouette","mixing.metric","local.struct"))),
              nrow = 1)
  
  temp_png_function(paste0("03_combined.png"))
  print(temp_grid)
  dev.off()
  
  temp_pdf_function(paste0("03_combined.pdf"))
  print(temp_grid)
  dev.off()
  
# }







  
  
  
# DATA OUTPUT IN TABLE FORMAT ---------------------------------------------

  temp_df <- NULL
  for(sampleID in temp_sampleIDs[!temp_sampleIDs %in% "CID4523"]){
      print(sampleID)
    
      if(sampleID %in% c("CID4471")){temp_pid = "BC P1"}
      if(sampleID %in% c("CID44971")){temp_pid = "BC P2"}
      if(sampleID %in% c("CID4513")){temp_pid = "BC P3"}
      if(sampleID %in% c("PID17267")){temp_pid = "PC P1"}
      if(sampleID %in% c("PID20033")){temp_pid = "PC P2"}
      if(sampleID %in% c("SCC180161")){temp_pid = "M P1"}
      
      temp_df_subset_patient <- 
        plot_data_combined_allcells[plot_data_combined_allcells$patient == temp_pid,,drop=F]
      
      for(metrics in c("silhouette","mixing.metric","local.struct")){
        
        print(metrics)
        temp_df_subset <- temp_df_subset_patient[,c( metrics,"condition", "cluster", "patient")]
        colnames(temp_df_subset) <- c("value","condition", "cluster", "patient")
        
      {temp_condition_ids <- c("FT", "CCS", "CT") 
        
        if(sampleID == "SCC180161") {
          temp_condition_ids <- c("FT", "CT", "CO") } 
        
        if(sampleID == "PID20033") {
          temp_condition_ids <- c("FT", "CT") }

      }
      
      for(condition in temp_condition_ids){
        temp_df_subset_patient_condition <- 
          temp_df_subset[temp_df_subset$condition == condition,,drop=F]
        
        temp_df_subset_new <- data.frame(mean = mean(temp_df_subset_patient_condition$value),
                                     sd = sd(temp_df_subset_patient_condition$value))
        temp_df_subset_new$sampleID <- sampleID
        temp_df_subset_new$condition <- condition
        temp_df_subset_new$metrics <- metrics
        
        temp_df <- rbind(temp_df,temp_df_subset_new)
        
      }
      
      
    }
    
    
  }
  
  write.csv(temp_df,
            "04_mean_sd_per_condition.csv")