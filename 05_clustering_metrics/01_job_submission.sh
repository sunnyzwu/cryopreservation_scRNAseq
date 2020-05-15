# Cryopreservation paper
# Submission of cluster metrics jobs
#
# environment
source activate Renv
#
# SAMPLE IDs to run loop submission of jobs
SAMPLEIDS="CID4471 CID44971 CID4513 PID17267 PID20033 SCC180161"
# These correlate with BC-P1, BC-P2, BC-P3, PC-P1, PC-P2 and M-P1, respectively
#
# INTEGRATEDPARAMS="TRUE"
# RUNBYKNN="TRUE"
# GLOBS
R="/share/ClusterShare/software/contrib/CTP_single_cell/tools/R_developers/config_R-3.5.0/bin/R"
TEMPPWD=$(pwd)
# output directory structure
mkdir output
for SAMPLENAME in ${SAMPLEIDS}; do
    mkdir output/${SAMPLENAME}_${NUM}
done
#
# SUBMIT JOBS
for SAMPLENAME in ${SAMPLEIDS}; do
  echo ${SAMPLENAME}

    cd output/${SAMPLENAME}_${NUM}
    # JOB NAME
    JOBNAME="kBET_${SAMPLENAME}_${NUM}"
    # R SCRIPT
    RSCRIPT="${TEMPPWD}/02_Rscript_processing.R"

    # PARAMS
    ## cell number to downsample whole dataset to
    ### Set to 'lowest' if you want to downsample to the lowest replicate size divided by 2 for FTvsFT2 comparison
    CELLNUMBER="lowest" # irrelavent param if running by cluster RUNBYCLUSTER=T
    ## Run by cluster ? downsamples by cluster
    RUNBYCLUSTER="FALSE"
    ## cells to downsample per cluster
    CELLNUMPERCLUST=${NUM} # irrelavent param if running by cluster RUNBYCLUSTER=F
    ## Run using the integrated assay
    RUNBYINTEGRATED="FALSE"
    # RUN STUART ET AL METRICS
    RUNBYSEURATMETRICS="TRUE"

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
      $SAMPLENAME \
      $CELLNUMBER \
      $RUNBYCLUSTER \
      $CELLNUMPERCLUST \
      $RUNBYINTEGRATED \
      $RUNBYSEURATMETRICS' \
      $RSCRIPT" \
      ./01_Rlog.txt

      cd ${TEMPPWD}
      echo "  submitted job ${JOBNAME}"
    done
