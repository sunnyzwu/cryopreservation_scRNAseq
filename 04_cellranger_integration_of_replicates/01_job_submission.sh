# GLOBS
CELLRANGER="/share/ClusterShare/software/contrib/CTP_single_cell/cellranger/cellranger-2.2.0/cellranger"
TEMPPWD=$(pwd)

# SAMPLE ID NAMES TO RUN
SAMPLEIDS="CID44971 CID4471 CID4513 PID17267 PID20033 SCC180161"
TEMPWD=$(pwd)
# output directory structure
mkdir output
for SAMPLENAME in ${SAMPLEIDS}; do
  mkdir output/cellranger_aggr_${SAMPLENAME}
done

# SUBMIT JOBS
for SAMPLENAME in ${SAMPLEIDS}; do

  echo ${SAMPLENAME}
  cd "output/cellranger_aggr_${SAMPLENAME}/"

  JOBNAME="aggr_${SAMPLENAME}"
  CONFIGFILE="${TEMPWD}/config_files/${SAMPLENAME}.csv"

    qsub \
    -cwd \
    -pe smp 16 \
    -l h_vmem=100G \
    -l mem_requested=10G \
    -P TumourProgression \
    -b y \
    -j y \
    -V \
    -N ${JOBNAME} \
    "${CELLRANGER} aggr \
    --id=${SAMPLENAME} \
    --csv=${CONFIGFILE} \
    --normalize=mapped \
    --nosecondary"

    cd ${TEMPWD}

done
