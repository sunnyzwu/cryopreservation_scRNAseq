
# GLOBS
CELLRANGER="/share/ClusterShare/software/contrib/CTP_single_cell/cellranger/cellranger-2.2.0/cellranger"
TEMPPWD=$(pwd)

# SAMPLE ID NAMES TO RUN
SAMPLEIDS="CID4497_CELLSUS CID4497_CHUNKS CID4471_CELLSUS CID4471_CHUNKS CID4513_CELLSUS CID4513_CHUNKS CID4523_CELLSUS CID4523_CHUNKS"

# SUBMIT JOBS
for SAMPLENAME in ${SAMPLEIDS}; do
  echo ${SAMPLENAME}
  COUNTJOBNAME="count_${SAMPLENAME}"
  OUTPUT_ID_STRING="count_${SAMPLENAME}_GRCh38"
  TRANSCRIPTOME="/share/ClusterShare/software/contrib/CTP_single_cell/cellranger/refdata/refdata-cellranger-GRCh38-1.2.0/"

  if [[ "${SAMPLENAME}" == "CID4497_CELLSUS" || "${SAMPLENAME}" == "CID4497_CHUNKS" || "${SAMPLENAME}" == "CID4471_CELLSUS" || "${SAMPLENAME}" == "CID4471_CHUNKS" ]]; then
    FASTQPATH="/paella/TumourProgressionGroupTemp/projects/data/cellranger_mkfastq/10X_4497+4471-06062018/HVLNFBGX5_1234/outs/fastq_path/Clinical_10X/${SAMPLENAME}"
  elif [[ "${SAMPLENAME}" == "CID4513_CELLSUS" || "${SAMPLENAME}" == "CID4513_CHUNKS" ]]; then
    FASTQPATH="/paella/TumourProgressionGroupTemp/projects/data/cellranger_mkfastq/10X_17267+4513cryo-23072018/HVFYNBGX5_1234/outs/fastq_path/Clinical_10X/${SAMPLENAME}/"
  elif [[ "${SAMPLENAME}" == "CID4523_CELLSUS" || "${SAMPLENAME}" == "CID4523_CHUNKS" ]]; then
    FASTQPATH="/paella/TumourProgressionGroupTemp/projects/data/cellranger_mkfastq/10X_4523x3+CITE-27082018/HGK72BGX7_1234/outs/fastq_path/Clinical_10X/${SAMPLENAME}/"
  fi

  # echo ${FASTQPATH}
  # ls ${FASTQPATH}
  SAMPLENAMEPATH=${SAMPLENAME//"_"/}
  if [[ "${SAMPLENAME}" == "CID4497_CELLSUS" ]]; then
    SAMPLENAMEPATH="CID44971CELLSUS"
  elif [[ "${SAMPLENAME}" == "CID4497_CHUNKS" ]]; then
    SAMPLENAMEPATH="CID44971CHUNKS"
  fi
  echo ${SAMPLENAMEPATH}
  cd "/paella/TumourProgressionGroupTemp/projects/data/cellranger_count/cryopreservation_paper/${SAMPLENAMEPATH}/"

qsub \
-cwd \
-pe smp 32 \
-l h_vmem=100G \
-l mem_requested=8G \
-P TumourProgression \
-b y \
-j y \
-V \
-N $COUNTJOBNAME \
"${CELLRANGER} count \
--id=$OUTPUT_ID_STRING \
--sample=$SAMPLENAME \
--transcriptome=$TRANSCRIPTOME \
--fastqs=$FASTQPATH \
--expect-cells=9000" \

  done

cd ${TEMPPWD}
