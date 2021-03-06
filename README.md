# cryopreservation_scRNAseq
Code related to the data processing and visualisation related to the Wu et al. manuscript "Cryopreservation of human cancers conserves tumour heterogeneity for single-cell multi-omics analysis". This study describes the benchmarking of cryopreservation across human cancers using scRNA-Seq on the 10X Chromium platform. 

Wu, S.Z., Roden, D.L., Al-Eryani, G. et al. Cryopreservation of human cancers conserves tumour heterogeneity for single-cell multi-omics analysis. Genome Med 13, 81 (2021). https://doi.org/10.1186/s13073-021-00885-z

## contents of repository
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
Calculating metrics and correlations from cellranger downsampled data. Both are calculated on the bulk and integrated cluster level. 

### 07 pathway_enrichment
Differential gene expresison of integrated clusters from each cryopreserved replicate using the MAST method. Subseqeunt gene ontology (GO) pathway enrichment performed using enrichProfiler. 

### 08 CITESeq_analysis
Seurat processing and visualisation of CITESeq data of an independent breast cancer cryopreserved as solid tissue (CT).

## Data availability
Raw scRNA-Seq data from this study has been deposited in the European Genome Phenome (EGA) Archive under the accession code EGAS00001005115 (https://ega-archive.org/studies/EGAS00001005115) This depository includes the demultiplexed paired ended reads (R1 and R2), Illumina indices and bam files processed using the Cellranger software. These paired ended reads (R1 and R2) and Illumina indices can be used for input for data reprocessing from the '01 cellranger count processing' step.

Processed scRNA-Seq data, in the form of raw and normalised expression matrices, can be found through the interactive data portal: https://singlecell.broadinstitute.org/single_cell/study/SCP1415

## Contact
Please email s.wu@garvan.org.au for any additional questions about the analytical methods used in this paper. All other relevant data are available from the authors upon request.
