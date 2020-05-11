# cryopreservation_scRNAseq
Code related to the data processing and visualisation related to the Wu et al. manuscript "Cryopreservation of human cancers conserves tumour heterogeneity for single-cell multi-omics analysis". This study describes the benchmarking of cryopreservation across human cancers using scRNA-Seq on the 10X Chromium platform. 

### 01_cellranger_count_processing  
job submission script for single-cell RNA-Seq processing using cellranger v2.2.0, mapping to the GRCh38 reference genome.

### 02_seurat_v3_individual_replicate_processing
job submission script and R script for processing individual seurat objects

### 03_seurat_v3_integration_of_replicates
job submission script and R script for integrating seurat objects

### 04_cellranger_integration_of_replicates

### 05 clustering_metrics
silhouette scores, mixing metric and local structure scores

### 06 bulk_and_cluster_correlations 

### 07 pathway_enrichment
differential gene expresison and pathway enrichment

### 08 CITESeq_analysis


