#!/usr/bin/env bash
set -euo pipefail

# Reproducible bulk RNA-seq workflow for PRJNA1167170.
# It preserves SRA cache and FASTQ files and never downloads accessions outside
# the eight explicitly listed below.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
RAW_DIR="${PROJECT_ROOT}/data/raw/PRJNA1167170"
TEST_RAW_DIR="${RAW_DIR}/test"
RESULT_DIR="${PROJECT_ROOT}/results/PRJNA1167170_bulk"
QC_DIR="${RESULT_DIR}/QC"
SALMON_DIR="${RESULT_DIR}/salmon"
INDEX_DIR="${PROJECT_ROOT}/reference/salmon_mouse_GRCm39"
R_SCRIPT="${SCRIPT_DIR}/PRJNA1167170_tximport_edgeR_fgsea.R"
R_LIB="${PROJECT_ROOT}/reference/R_library"
THREADS=8

RUNS=(SRR30835575 SRR30835574 SRR30835573 SRR30835572 SRR30835571 SRR30835570 SRR30835569 SRR30835568)
SAMPLES=(con1 con2 con3 con4 SD1 SD2 SD3 SD4)
REMAINING_RUNS=(SRR30835574 SRR30835573 SRR30835572 SRR30835571 SRR30835570 SRR30835569 SRR30835568)

mkdir -p "${RAW_DIR}" "${RAW_DIR}/fasterq_tmp" "${QC_DIR}" "${SALMON_DIR}"

for tool in prefetch fasterq-dump fastqc multiqc salmon Rscript; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    printf 'ERROR: required tool not found: %s\n' "${tool}" >&2
    exit 1
  fi
done

if [[ "$(salmon --version)" != *"2.7.0"* ]]; then
  printf 'ERROR: Salmon 2.7.0 is required; found: %s\n' "$(salmon --version)" >&2
  exit 1
fi

if [[ ! -d "${INDEX_DIR}" ]]; then
  printf 'ERROR: Salmon index is missing: %s\n' "${INDEX_DIR}" >&2
  exit 1
fi

if [[ ! -f "${R_SCRIPT}" ]]; then
  printf 'ERROR: R analysis script is missing: %s\n' "${R_SCRIPT}" >&2
  exit 1
fi

fastq_path() {
  local run="$1"
  local mate="$2"
  if [[ "${run}" == "SRR30835575" ]]; then
    printf '%s/%s_%s.fastq\n' "${TEST_RAW_DIR}" "${run}" "${mate}"
  else
    printf '%s/%s_%s.fastq\n' "${RAW_DIR}" "${run}" "${mate}"
  fi
}

# Prevent a predictable disk-full failure. Based on the observed SRR30835575
# expansion, the seven missing runs need about 144 GiB for retained FASTQ+SRA.
# A further 35 GiB is reserved for fasterq-dump scratch space and filesystem
# safety. The gate is skipped once all seven FASTQ pairs already exist.
missing_fastq=0
for run in "${REMAINING_RUNS[@]}"; do
  r1="$(fastq_path "${run}" 1)"
  r2="$(fastq_path "${run}" 2)"
  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    missing_fastq=1
    break
  fi
done

if [[ "${missing_fastq}" -eq 1 ]]; then
  available_kib="$(df -Pk "${PROJECT_ROOT}" | awk 'NR==2 {print $4}')"
  required_kib=$((179 * 1024 * 1024))
  if (( available_kib < required_kib )); then
    available_gib="$(awk -v kb="${available_kib}" 'BEGIN {printf "%.1f", kb/1024/1024}')"
    printf 'ERROR: insufficient free space for retained SRA + uncompressed FASTQ + conversion scratch.\n' >&2
    printf 'Available: %s GiB; required before starting: at least 179 GiB.\n' "${available_gib}" >&2
    printf 'No new accession was downloaded. Add storage, then rerun this script.\n' >&2
    exit 2
  fi
fi

# Download and convert only the seven not-yet-downloaded runs. Existing complete
# FASTQ pairs are retained and skipped. A one-mate-only state is treated as an
# error rather than overwritten.
for run in "${REMAINING_RUNS[@]}"; do
  r1="$(fastq_path "${run}" 1)"
  r2="$(fastq_path "${run}" 2)"
  if [[ -s "${r1}" && -s "${r2}" ]]; then
    printf 'FASTQ pair exists; skipping conversion: %s\n' "${run}"
    continue
  fi
  if [[ -e "${r1}" || -e "${r2}" ]]; then
    printf 'ERROR: incomplete pre-existing FASTQ pair for %s; refusing to overwrite.\n' "${run}" >&2
    exit 1
  fi

  (
    cd "${PROJECT_ROOT}"
    prefetch --max-size u "${run}"
  )
  mkdir -p "${RAW_DIR}/fasterq_tmp/${run}"
  fasterq-dump \
    --split-files \
    --outdir "${RAW_DIR}" \
    --temp "${RAW_DIR}/fasterq_tmp/${run}" \
    -e "${THREADS}" \
    "${PROJECT_ROOT}/${run}"

  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    printf 'ERROR: fasterq-dump did not produce a complete FASTQ pair for %s.\n' "${run}" >&2
    exit 1
  fi
done

# Confirm that all eight paired inputs are present before QC or quantification.
FASTQS=()
for run in "${RUNS[@]}"; do
  r1="$(fastq_path "${run}" 1)"
  r2="$(fastq_path "${run}" 2)"
  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    printf 'ERROR: required FASTQ pair is missing for %s.\n' "${run}" >&2
    exit 1
  fi
  FASTQS+=("${r1}" "${r2}")
done

# FastQC is run only for files whose report is not already present.
for fastq in "${FASTQS[@]}"; do
  base="$(basename "${fastq}" .fastq)"
  if [[ -s "${QC_DIR}/${base}_fastqc.zip" && -s "${QC_DIR}/${base}_fastqc.html" ]]; then
    printf 'FastQC output exists; skipping: %s\n' "${base}"
  else
    fastqc --threads 2 --outdir "${QC_DIR}" "${fastq}"
  fi
done

if [[ ! -s "${QC_DIR}/multiqc_report.html" ]]; then
  multiqc "${QC_DIR}" --outdir "${QC_DIR}" --filename multiqc_report.html
else
  printf 'MultiQC report exists; skipping.\n'
fi

# Quantify every sample with the same index and Salmon settings. Existing
# complete quantifications are retained and skipped.
for i in "${!RUNS[@]}"; do
  run="${RUNS[$i]}"
  sample="${SAMPLES[$i]}"
  r1="$(fastq_path "${run}" 1)"
  r2="$(fastq_path "${run}" 2)"
  out="${SALMON_DIR}/${sample}"
  if [[ -s "${out}/quant.sf" && -s "${out}/aux_info/meta_info.json" && -s "${out}/lib_format_counts.json" ]]; then
    printf 'Salmon output exists; skipping: %s\n' "${sample}"
    continue
  fi
  if [[ -e "${out}" ]]; then
    printf 'ERROR: incomplete Salmon output exists for %s; refusing to overwrite: %s\n' "${sample}" "${out}" >&2
    exit 1
  fi
  salmon quant \
    -i "${INDEX_DIR}" \
    -l A \
    -1 "${r1}" \
    -2 "${r2}" \
    --validateMappings \
    -p "${THREADS}" \
    -o "${out}"
done

if [[ -e "${RESULT_DIR}/reproducibility_info.txt" ]]; then
  printf 'Analysis outputs already exist; preserving them and stopping before R analysis.\n' >&2
  exit 1
fi

R_LIBS_USER="${R_LIB}" Rscript "${R_SCRIPT}" "${PROJECT_ROOT}"

printf 'PRJNA1167170 bulk workflow completed successfully.\n'
