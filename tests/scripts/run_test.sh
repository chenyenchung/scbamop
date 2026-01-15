#!/usr/bin/env bash
set -euo pipefail

python=""
scbamop=""
mode="split"
seed=1
num_reads=200
work_dir=""
mapq_threshold=20
mapq_tie_rate=0.2
dup_rate=0.2
missing_cb_rate=0.05
missing_ub_rate=0.05
cb_unknown_rate=0.1
secondary_rate=0.05
supplementary_rate=0.05
unsafe_label="false"
unsafe_label_flag=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        --python)
            python="$2"
            shift 2
            ;;
        --scbamop)
            scbamop="$2"
            shift 2
            ;;
        --mode)
            mode="$2"
            shift 2
            ;;
        --seed)
            seed="$2"
            shift 2
            ;;
        --num-reads)
            num_reads="$2"
            shift 2
            ;;
        --work-dir)
            work_dir="$2"
            shift 2
            ;;
        --mapq-threshold)
            mapq_threshold="$2"
            shift 2
            ;;
        --mapq-tie-rate)
            mapq_tie_rate="$2"
            shift 2
            ;;
        --dup-rate)
            dup_rate="$2"
            shift 2
            ;;
        --missing-cb-rate)
            missing_cb_rate="$2"
            shift 2
            ;;
        --missing-ub-rate)
            missing_ub_rate="$2"
            shift 2
            ;;
        --cb-unknown-rate)
            cb_unknown_rate="$2"
            shift 2
            ;;
        --secondary-rate)
            secondary_rate="$2"
            shift 2
            ;;
        --supplementary-rate)
            supplementary_rate="$2"
            shift 2
            ;;
        --unsafe-label)
            unsafe_label="true"
            shift 1
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
 done

if [ -z "$python" ] || [ -z "$scbamop" ] || [ -z "$work_dir" ]; then
    echo "Missing required arguments" >&2
    exit 1
fi

if [ "$num_reads" -ge 1000000 ] && [ -z "${RUN_STRESS_TESTS:-}" ]; then
    echo "RUN_STRESS_TESTS not set; skipping stress test." >&2
    exit 0
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(CDPATH= cd -- "${script_dir}/../.." && pwd)

mkdir -p "$work_dir"
input_bam="${work_dir}/input.bam"
metadata="${work_dir}/metadata.csv"
output_dir="${work_dir}/output"

if [ -d "$output_dir" ]; then
    rm -rf "$output_dir"
fi
mkdir -p "$output_dir"

if [ "$unsafe_label" = "true" ]; then
    unsafe_label_flag="--unsafe-label"
fi

"$python" "${repo_dir}/tests/scripts/generate_sam_fixture.py" \
    --seed "$seed" \
    --output-bam "$input_bam" \
    --metadata "$metadata" \
    --num-reads "$num_reads" \
    --dup-rate "$dup_rate" \
    --missing-cb-rate "$missing_cb_rate" \
    --missing-ub-rate "$missing_ub_rate" \
    --cb-unknown-rate "$cb_unknown_rate" \
    --secondary-rate "$secondary_rate" \
    --supplementary-rate "$supplementary_rate" \
    --mapq-tie-rate "$mapq_tie_rate" \
    $unsafe_label_flag

if [ "$mode" = "dedup" ]; then
    "$scbamop" split -f "$input_bam" -m "$metadata" -o "$output_dir" -q "$mapq_threshold" -d
else
    "$scbamop" split -f "$input_bam" -m "$metadata" -o "$output_dir" -q "$mapq_threshold"
fi

"$python" "${repo_dir}/tests/scripts/validate_split.py" \
    --input-bam "$input_bam" \
    --metadata "$metadata" \
    --output-dir "$output_dir" \
    --mode "$mode" \
    --mapq-threshold "$mapq_threshold"
