#!/usr/bin/env bash
set -euo pipefail

# Usage:
# ./run_deepsomatic.sh <work_dir> <output_dir> <normal_dir> <ref_dir> <tumor_sample> <normal_sample>

WORK_DIR="${1:?Missing work_dir}"
OUTPUT_DIR="${2:?Missing output_dir}"
NORMAL_DIR="${3:?Missing normal_dir}"
REF_DIR="${4:?Missing ref_dir}"
TUMOR_SAMPLE="${5:?Missing tumor_sample}"
NORMAL_SAMPLE="${6:?Missing normal_sample}"

DOCKER_IMAGE="google/deepsomatic:1.8.0"
MODEL_TYPE="WGS"
NUM_SHARDS="30"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}/intermediate_results_dir"

docker run \
  -v "${WORK_DIR}:/work" \
  -v "${NORMAL_DIR}:/normal" \
  -v "${REF_DIR}:/ref" \
  -v "${OUTPUT_DIR}:/output" \
  "${DOCKER_IMAGE}" \
  run_deepsomatic \
  --model_type="${MODEL_TYPE}" \
  --ref=/ref/hg19.fa \
  --reads_tumor="/work/${TUMOR_SAMPLE}/BAMs/${TUMOR_SAMPLE}.bam" \
  --reads_normal="/normal/BAMs/${NORMAL_SAMPLE}.bam" \
  --output_vcf="/output/deepsomatic_${TUMOR_SAMPLE}.vcf.gz" \
  --output_gvcf="/output/deepsomatic_${TUMOR_SAMPLE}.g.vcf.gz" \
  --sample_name_tumor="tumor" \
  --sample_name_normal="normal" \
  --num_shards="${NUM_SHARDS}" \
  --logging_dir=/output/logs \
  --vcf_stats_report=true \
  --intermediate_results_dir /output/intermediate_results_dir
