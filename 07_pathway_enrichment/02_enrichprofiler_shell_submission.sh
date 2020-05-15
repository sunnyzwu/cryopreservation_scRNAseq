# Cryopreservation paper
# Submission of DGE jobs across cryo conditions
# R 3.4.1
#
# SAMPLE IDs
SAMPLEIDS="CID4471 CID44971 CID4513 CID4523 PID17267 PID20033 SCC180161"
#
# GLOBS
R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R-3.4.1/bin/R"
TEMPPWD=$(pwd)
# output directory structure
# mkdir DGE_output
# for SAMPLENAME in ${SAMPLEIDS}; do
#   mkdir DGE_output/${SAMPLENAME}
# done
#
# SUBMIT JOBS
for SAMPLENAME in ${SAMPLEIDS}; do

  echo ${SAMPLENAME}
  cd DGE_output/${SAMPLENAME}/

  JOBNAME="enrichprofiler_${SAMPLENAME}"
  RSCRIPT="${TEMPPWD}/04_Rscript_processing.R"

  # submit job
  qsub \
  -cwd \
  -pe smp 4 \
  -l h_vmem=100G \
  -l mem_requested=10G \
  -b y \
  -j y \
  -V \
  -P TumourProgression \
  -N $JOBNAME \
  "${R} CMD BATCH \
  --no-save '--args \
  $SAMPLENAME' \
  $RSCRIPT" \
  ./02_Rlog.txt

  cd ${TEMPPWD}
  echo "  submitted job ${JOBNAME}"


done
