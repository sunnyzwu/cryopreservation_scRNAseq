# GO ENRICHMENT PER CRYO CONDITION
# SCRIPT A 
# R.v3.4.1 ENRICHPROFILER BUG IN R v3.5.0
# 01 PARSE ARGUMENTS  ------------------------------------------------------------

temp_start_time <- Sys.time()
print(temp_start_time)
temp_args <-
  commandArgs(trailingOnly = T)

# 01 CELL TYPE
temp_sampleIDs  <- 
  temp_args[1]

# 02 SETUP  ------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
 
# 03 LOAD DEG FILES -----------------------------------------------------

temp_conditions <- 
  c("Fresh Tissue", "Cryopreserved Cell Suspension", "Cryopreserved Tissue")

if(temp_sampleIDs == "SCC180161") {
  temp_conditions <- c("Fresh Tissue", "Cryopreserved Tissue", "Cryopreserved Overnight") } 

if(temp_sampleIDs == "PID20033") {
  temp_conditions <- c("Fresh Tissue", "Cryopreserved Tissue") 
  } 


for(sampleID in temp_sampleIDs) {
  
  temp_go_combined_df <- NULL
  for(condition in temp_conditions){
    print(condition)
  temp_cluster_allmarkers <- 
    read.csv(paste0("01_DGE_",sampleID,"_", condition, ".csv"))
  
  print("DEGs per cluster")
  print(table(temp_cluster_allmarkers$cluster))
  
  temp_clusterMarkers <- NULL
  for (cluster in as.character(unique(temp_cluster_allmarkers$cluster)))
  {
    temp_signatureEntrezIDs <- bitr(temp_cluster_allmarkers[temp_cluster_allmarkers$cluster == cluster, ]$gene,
                                    fromType = "SYMBOL",
                                    toType = "ENTREZID",
                                    OrgDb = "org.Hs.eg.db"
    )
    
    temp_cm <-
      data.frame(EntrezID = temp_signatureEntrezIDs$ENTREZID,
                 clusterID = cluster)
    if (!is.null(temp_clusterMarkers))
    {
      temp_clusterMarkers <- rbind(temp_clusterMarkers, 
                                   temp_cm)
    } else {
      temp_clusterMarkers <- temp_cm
    }}
  
  temp_go <- compareCluster(
    EntrezID ~ clusterID,
    data = temp_clusterMarkers,
    fun = "enrichGO",
    OrgDb = "org.Hs.eg.db",
    readable = T
  )
  
  temp_go <- 
    as.data.frame(temp_go)
  
  temp_go$condition <- paste0(condition)
  
  write.csv(temp_go,
            paste0("02_enrichGO_",condition,".csv"))
  
  temp_go_combined_df <- 
    rbind(temp_go_combined_df,
          temp_go)

  }
  
  temp_go_combined_df_sorted <- 
    temp_go_combined_df[ order( temp_go_combined_df[, "clusterID"],
                                # temp_go_combined_df[, "ONTOLOGY"],
                                temp_go_combined_df[, "Description"] ),]
  
  temp_go_combined_df_sorted$logp_val <- 
    (-(log10(temp_go_combined_df_sorted$p.adjust)))
  
  write.csv(temp_go_combined_df_sorted,
            paste0("03_enrichGO_combined_sorted.csv"))

}


# FINISH ------------------------------------------------------------------
  
  print("start time")
  print(temp_start_time)
  print("finish time")
  print(Sys.time())
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  