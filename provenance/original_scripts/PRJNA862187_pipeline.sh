#!/usr/bin/env bash
set -euo pipefail

# Stage-gated, reproducible FASTQ/QC/Salmon workflow for the 12 prespecified
# PRJNA862187 colon RNA-seq samples. Fecal metagenomic runs are never included.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RAW_DIR="${PROJECT_ROOT}/data/raw/PRJNA862187"
OUT_DIR="${PROJECT_ROOT}/results/PRJNA862187_validation"
QC_DIR="${OUT_DIR}/QC"
SALMON_DIR="${OUT_DIR}/salmon"
INDEX_DIR="${PROJECT_ROOT}/reference/salmon_mouse_GRCm39"
PLAN="${OUT_DIR}/PRJNA862187_analysis_plan_frozen.md"
METADATA="${OUT_DIR}/PRJNA862187_sample_metadata.csv"
PILOT_MARKER="${OUT_DIR}/PILOT_PASSED.txt"
THREADS="${THREADS:-8}"
MODE="${1:-pilot}"

RUNS=(SRR22944153 SRR22944152 SRR22944155 SRR22944148 SRR22944145 SRR22944156 SRR22944147 SRR22944146 SRR22944154 SRR22944151 SRR22944150 SRR22944149)
SAMPLES=(EC_SCD_1 EC_SCD_2 EC_SCD_3 SD_SCD_1 SD_SCD_2 SD_SCD_3 EC_HFD_1 EC_HFD_2 EC_HFD_3 SD_HFD_1 SD_HFD_2 SD_HFD_3)
VOLS=(053 052 055 048 045 056 047 046 054 051 050 049)
BYTES1=(1252359216 1300973311 1273821472 1276333064 1275789485 1331110739 1250287164 1250226907 1296023260 1217415325 1314777653 1287113440)
BYTES2=(1253440950 1301384346 1272839099 1277189821 1275771858 1330878762 1250817857 1250382200 1295797168 1216314765 1316030772 1287318883)
MD5_1=(31e3fd87536423d8977a84edd382d4bd bf5b68800716f4b808b08f535335316a 074482c241b7afd528d3c4fcfa5e3250 0a47b09130e50ecdc47e3f5ec372f28c e2746f5d3189eb4de2db7a26cdc1eb40 05a86f85ee8b30ca29190dce45131b7d 310bfcda97f8b68a0e9b0f5901251f4d 6fd24ed1f63c8556e3cd6614f50cedf4 871da90161c5f6829f56edfed12df9fa 34fad3c5aed6ce984edb456d2d865670 349713c10211f5e097658e034f2f38d4 1b13ef5e03d5acbf2657a2d8b97aff93)
MD5_2=(6d100ab4c8bd021dceb594018e2db7dd 0d422f9331abac22329e43f9efb6fef1 b1dc66ba49ab87b75bfece629775f84e 8f5c28c5d6b78f756ff440e848a95303 53095ce540c0b51517d627cbb01da063 37d824f87f72008cabd11b337aa74987 5339689288cbbdf1e421d27e55f46aef 7dee4c93242a5d65f5b3bcc6d687c02b 77b9b6d94d281854281144bd03c96e95 2ee1552c2823ee6efd12d5941e92427b f78479910cc4b35a68f0192815fe10c9 3e2dcd7d3bba7a1cda95e3f02119bc51)

if [[ "${MODE}" != "pilot" && "${MODE}" != "full" ]]; then
  printf 'Usage: %s [pilot|full]\n' "$0" >&2
  exit 2
fi

for tool in curl fastqc multiqc salmon jq unzip; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    printf 'ERROR: required tool not found: %s\n' "${tool}" >&2
    exit 1
  fi
done

if [[ ! -s "${PLAN}" || ! -s "${METADATA}" ]]; then
  printf 'ERROR: frozen plan or verified metadata is missing. Refusing data acquisition.\n' >&2
  exit 1
fi
if [[ "$(salmon --version)" != *"2.7.0"* ]]; then
  printf 'ERROR: Salmon 2.7.0 required; found %s\n' "$(salmon --version)" >&2
  exit 1
fi
if [[ ! -s "${INDEX_DIR}/info.json" ]]; then
  printf 'ERROR: Salmon GRCm39/M39 index is missing: %s\n' "${INDEX_DIR}" >&2
  exit 1
fi
if [[ "${MODE}" == "full" && ! -s "${PILOT_MARKER}" ]]; then
  printf 'ERROR: pilot gate has not passed; refusing the remaining 11 downloads.\n' >&2
  exit 1
fi

mkdir -p "${RAW_DIR}" "${QC_DIR}" "${SALMON_DIR}"

if [[ "${MODE}" == "full" ]]; then
  available_kib="$(df -Pk "${PROJECT_ROOT}" | awk 'NR==2 {print $4}')"
  required_kib=$((40 * 1024 * 1024))
  if (( available_kib < required_kib )); then
    printf 'ERROR: full mode requires at least 40 GiB free for compressed FASTQ and outputs; currently %.1f GiB.\n' "$(awk -v x="${available_kib}" 'BEGIN {print x/1024/1024}')" >&2
    exit 2
  fi
fi

download_mate() {
  local run="$1" vol="$2" mate="$3" expected="$4" expected_md5="$5"
  local out="${RAW_DIR}/${run}_${mate}.fastq.gz"
  local url="https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR229/${vol}/${run}/${run}_${mate}.fastq.gz"
  python3 "${SCRIPT_DIR}/http_range_download.py" \
    --url "${url}" --output "${out}" --size "${expected}" --md5 "${expected_md5}" \
    --workers "${DOWNLOAD_WORKERS:-128}" --chunk-mib 1
}

process_sample() {
  local i="$1"
  local run="${RUNS[$i]}" sample="${SAMPLES[$i]}" vol="${VOLS[$i]}"
  local r1="${RAW_DIR}/${run}_1.fastq.gz" r2="${RAW_DIR}/${run}_2.fastq.gz"
  download_mate "${run}" "${vol}" 1 "${BYTES1[$i]}" "${MD5_1[$i]}"
  download_mate "${run}" "${vol}" 2 "${BYTES2[$i]}" "${MD5_2[$i]}"

  for fq in "${r1}" "${r2}"; do
    base="$(basename "${fq}" .fastq.gz)"
    if [[ ! -s "${QC_DIR}/${base}_fastqc.zip" || ! -s "${QC_DIR}/${base}_fastqc.html" ]]; then
      fastqc --threads 2 --outdir "${QC_DIR}" "${fq}"
    fi
  done

  quant="${SALMON_DIR}/${sample}"
  if [[ ! -s "${quant}/quant.sf" || ! -s "${quant}/aux_info/meta_info.json" || ! -s "${quant}/lib_format_counts.json" ]]; then
    if [[ -e "${quant}" ]]; then
      printf 'ERROR: incomplete Salmon output exists; refusing overwrite: %s\n' "${quant}" >&2
      exit 1
    fi
    salmon quant \
      -i "${INDEX_DIR}" \
      -l A \
      -1 "${r1}" \
      -2 "${r2}" \
      --validateMappings \
      -p "${THREADS}" \
      -o "${quant}"
  fi
}

if [[ "${MODE}" == "pilot" ]]; then
  process_sample 0
  multiqc "${QC_DIR}" --outdir "${QC_DIR}" --filename pilot_multiqc_report.html --force

  pilot_quant="${SALMON_DIR}/EC_SCD_1"
  mapping_rate="$(jq -r '.percent_mapped' "${pilot_quant}/aux_info/meta_info.json")"
  processed="$(jq -r '.num_processed' "${pilot_quant}/aux_info/meta_info.json")"
  orphan="$(jq -r '.num_orphan' "${pilot_quant}/aux_info/meta_info.json")"
  library_type="$(jq -r '.detected_library_type' "${pilot_quant}/aux_info/meta_info.json")"
  compatible="$(jq -r '.compatible_fragment_ratio' "${pilot_quant}/lib_format_counts.json")"
  read_length_r1="$(unzip -p "${QC_DIR}/SRR22944153_1_fastqc.zip" '*/fastqc_data.txt' | awk -F '\t' '$1=="Sequence length" {print $2}')"
  read_length_r2="$(unzip -p "${QC_DIR}/SRR22944153_2_fastqc.zip" '*/fastqc_data.txt' | awk -F '\t' '$1=="Sequence length" {print $2}')"
  adapter_r1="$(unzip -p "${QC_DIR}/SRR22944153_1_fastqc.zip" '*/summary.txt' | awk -F '\t' '$2=="Adapter Content" {print $1}')"
  adapter_r2="$(unzip -p "${QC_DIR}/SRR22944153_2_fastqc.zip" '*/summary.txt' | awk -F '\t' '$2=="Adapter Content" {print $1}')"
  quality_r1="$(unzip -p "${QC_DIR}/SRR22944153_1_fastqc.zip" '*/summary.txt' | awk -F '\t' '$2=="Per base sequence quality" {print $1}')"
  quality_r2="$(unzip -p "${QC_DIR}/SRR22944153_2_fastqc.zip" '*/summary.txt' | awk -F '\t' '$2=="Per base sequence quality" {print $1}')"
  orphan_rate="$(awk -v o="${orphan}" -v p="${processed}" 'BEGIN {printf "%.6f", 100*o/p}')"

  printf 'sample,run,read_length_R1,read_length_R2,fastqc_per_base_quality_R1,fastqc_per_base_quality_R2,fastqc_adapter_R1,fastqc_adapter_R2,processed_fragments,mapping_rate_percent,library_type,orphan_fragments,orphan_rate_percent,compatible_ratio\n' > "${OUT_DIR}/pilot_qc_summary.csv"
  printf 'EC_SCD_1,SRR22944153,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
    "${read_length_r1}" "${read_length_r2}" "${quality_r1}" "${quality_r2}" "${adapter_r1}" "${adapter_r2}" \
    "${processed}" "${mapping_rate}" "${library_type}" "${orphan}" "${orphan_rate}" "${compatible}" >> "${OUT_DIR}/pilot_qc_summary.csv"

  if awk -v m="${mapping_rate}" -v c="${compatible}" 'BEGIN {exit !(m >= 60 && c >= 0.50)}' \
     && [[ "${read_length_r1}" == "151" && "${read_length_r2}" == "151" && "${library_type}" != "U" ]]; then
    printf 'Pilot passed on %s: paired 151-bp reads; mapping %.6f%%; library %s; compatible ratio %s.\n' \
      "$(date '+%Y-%m-%d %H:%M:%S %Z')" "${mapping_rate}" "${library_type}" "${compatible}" > "${PILOT_MARKER}"
    printf 'PILOT PASS\n'
  else
    printf 'PILOT STOP: technical gate failed; remaining samples were not downloaded.\n' >&2
    exit 3
  fi
else
  for i in "${!RUNS[@]}"; do
    process_sample "${i}"
  done
  multiqc "${QC_DIR}" --outdir "${QC_DIR}" --filename multiqc_report.html --force
  printf 'All 12 prespecified samples completed FastQC and Salmon.\n'
fi
