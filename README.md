# cryopreservation_scRNAseq
Code related to the data processing and visualisation related to the Wu et al. manuscript "Cryopreservation of human cancers conserves tumour heterogeneity for single-cell multi-omics analysis". This study describes the benchmarking of cryopreservation across human cancers using scRNA-Seq on the 10X Chromium platform. 

### 01_cellranger_count_processing  
job submission shell script for single-cell RNA-Seq processing using cellranger v2.2.0 - count function. Reads are mapped to the GRCh38 reference genome. Input requires demultiplexed fastqs e.g. via cellranger mkfastq.

### 02_seurat_v3_individual_replicate_processing
job submission shell script and R script for processing individual seurat objects using the seurat v3 method.

### 03_seurat_v3_integration_of_replicates
job submission shell script and R script for integrating seurat objects using the seurat v3 anchoring method as described by Stuart et al. (2019).

### 04_cellranger_integration_of_replicates
job submission shell script and config files for the cellranger aggr function. Integration of matched replicates by normalising based on the number of mapped sequencing reads. 

### 05 clustering_metrics
job submission shell script and R scripts for the downsampling and computation of silhouette scores, mixing metric and local structure scores, as described by Stuart et al. (2019).

### 06 bulk_and_cluster_correlations 



### 07 pathway_enrichment
differential gene expresison and pathway enrichment

### 08 CITESeq_analysis


