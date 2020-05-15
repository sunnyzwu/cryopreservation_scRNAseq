# Cryopreservation paper
# Submission of CCA jobs
#
#
# SAMPLE IDs
SAMPLEIDS="CID4471"
# SAMPLEIDS="CID4471 CID44971 CID4513 CID4523 PID17267 PID20033 SCC180161"
#
# GLOBS
R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R"
TEMPPWD=$(pwd)
# output directory structure
mkdir output
for SAMPLENAME in ${SAMPLEIDS}; do
  mkdir output/CCA_${SAMPLENAME}
done
#
# SUBMIT JOBS
for SAMPLENAME in ${SAMPLEIDS}; do

  echo ${SAMPLENAME}
  cd output/CCA_${SAMPLENAME}

  # subset sample sheet
  head ${TEMPPWD}/config/sample_input_file.csv -n 1 > ./sample_input_file_subset.csv
  grep "${SAMPLENAME}" ${TEMPPWD}/config/sample_input_file.csv >> ./sample_input_file_subset.csv
  SUBSETINPUTFILE="./sample_input_file_subset.csv"
  # reference genome
  SPECIES="human"
  # seurat gene input file
  SEURATGENEINPUTFILE="$TEMPPWD/config/seurat_gene_input_file.csv"
  # seurat params file
  SEURATPARAMFILE="$TEMPPWD/config/seurat_CCA_params_file.csv"
  # seurat script
  SEURATSCRIPT="$TEMPPWD/config/seurat_CCA_HPC_processing.R"
  # job name
  SEURATJOBNAME="sCCA_$SAMPLENAME"
  # path to individual seurat objects
  OBJECTSPATH="/share/ScratchGeneral/sunwu/projects/cryopreservation_paper/seurat_v3/individual_seurat_processing/output/"

  # cd to directory
  cd $TEMPPWD/output/CCA_${SAMPLENAME}/
  # submit job
  qsub \
  -cwd \
  -pe smp 32 \
  -l h_vmem=100G \
  -l mem_requested=10G \
  -b y \
  -j y \
  -V \
  -P TumourProgression \
  -N $SEURATJOBNAME\
  "${R} CMD BATCH \
  --no-save '--args \
  $SAMPLENAME \
  $SPECIES \
  $(pwd) \
  $SUBSETINPUTFILE \
  $SEURATGENEINPUTFILE \
  $SEURATPARAMFILE \
  $OBJECTSPATH' \
  $SEURATSCRIPT"

  cd ${TEMPPWD}

  done
