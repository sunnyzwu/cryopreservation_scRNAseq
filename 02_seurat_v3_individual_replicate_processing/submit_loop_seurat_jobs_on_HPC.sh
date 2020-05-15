# this script is for the analysis of all cryopreservation samples using the seurat_emptydrops R script
# specify input files in config/ folder
#     - config/sample_input_file.csv
#       sample_num,sampleID,path_to_matrix
#       1,CID4471,/share/path/matrix
#       2,CID4471CELLSUS,/share/path/matrix
#       3,CID4471CHUNKS,/share/path/matrix

R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.5.0/bin/R"
TEMPPWD=$(pwd)
mkdir output/

for samplename in $(cut -d ',' -f 2 ./config/input_sample_file.csv | tail -n +2); do
  mkdir output/seurat_$samplename/
done

for samplenum in $(cut -d ',' -f 1 $TEMPPWD/config/input_sample_file.csv | tail -n +2); do

  SAMPLENAME=$(cut -d ',' -f 2 $TEMPPWD/config/input_sample_file.csv | tail -n +2 | head -n $samplenum | tail -n 1)
  REFERENCEGENOME="human"
  SEURATGENEINPUTFILE="$TEMPPWD/config/seurat_gene_input_file.csv"
  SEURATPARAMFILE="$TEMPPWD/config/seurat_params_file.csv"
  SEURATSCRIPT="$TEMPPWD/config/seurat_HPC_processing.R"
  SEURATJOBNAME="s_$SAMPLENAME"
  MATRIXPATH=$(cut -d ',' -f 3 $TEMPPWD/config/input_sample_file.csv | tail -n +2 | head -n $samplenum | tail -n 1)
  TEMPWD="$TEMPPWD/output/seurat_$SAMPLENAME/"

  cd $TEMPPWD/output/seurat_$SAMPLENAME/

  qsub \
  -cwd \
  -pe smp 16 \
  -l h_vmem=300G \
  -P TumourProgression \
  -b y \
  -j y \
  -V \
  -N $SEURATJOBNAME\
  "${R} CMD BATCH \
  --no-save '--args \
  $SAMPLENAME \
  $REFERENCEGENOME \
  $TEMPWD \
  $MATRIXPATH \
  $SEURATGENEINPUTFILE \
  $SEURATPARAMFILE' \
  $SEURATSCRIPT"

  cd $TEMPPWD

done
