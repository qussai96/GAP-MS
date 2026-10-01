#!/usr/bin/env bash

set -Eeuo pipefail

# Run the complete GAP-MS tutorial dataset using the image published on
# Docker Hub. Paths can be overridden with environment variables; for example:
#   REFRESH_IMAGE=1 OUTPUT_DIR=/scratch/$USER/gapms_test ./tests/run_gapms_image.sh

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"

IMAGE_URI="${IMAGE_URI:-docker://qussaiab96/gapms:latest}"
IMAGE_PATH="${IMAGE_PATH:-$SCRIPT_DIR/gapms_dockerhub_latest.sif}"
TUTORIAL_DIR="${TUTORIAL_DIR:-$REPO_DIR/tutorials}"
OUTPUT_DIR="${OUTPUT_DIR:-$SCRIPT_DIR/tutorial_image_output_$(date +%Y%m%d_%H%M%S)}"
REFRESH_IMAGE="${REFRESH_IMAGE:-0}"

required_files=(
    Helixer_prediction.gtf
    Helixer_prediction.faa
    NC_023759.2.test.bam
    NC_023759.2_genome_Sequence.fna
    Helixer_MS_identified_peptides.tsv
    Reference_annotation.gtf
    Reference_annotation.faa
)

if ! command -v apptainer >/dev/null 2>&1; then
    echo "Error: apptainer is not available in PATH." >&2
    exit 1
fi

for filename in "${required_files[@]}"; do
    if [[ ! -f "$TUTORIAL_DIR/$filename" ]]; then
        echo "Error: missing tutorial input: $TUTORIAL_DIR/$filename" >&2
        exit 1
    fi
done

mkdir -p "$(dirname -- "$IMAGE_PATH")"

if [[ "$REFRESH_IMAGE" == "1" ]]; then
    echo "[$(date '+%F %T')] Refreshing image from $IMAGE_URI"
    apptainer pull --force "$IMAGE_PATH" "$IMAGE_URI"
elif [[ ! -f "$IMAGE_PATH" ]]; then
    echo "[$(date '+%F %T')] Pulling image from $IMAGE_URI"
    apptainer pull "$IMAGE_PATH" "$IMAGE_URI"
else
    echo "[$(date '+%F %T')] Reusing image: $IMAGE_PATH"
    echo "Set REFRESH_IMAGE=1 to pull the latest image again."
fi

mkdir -p "$OUTPUT_DIR"
TUTORIAL_DIR="$(cd -- "$TUTORIAL_DIR" && pwd)"
OUTPUT_DIR="$(cd -- "$OUTPUT_DIR" && pwd)"
IMAGE_PATH="$(cd -- "$(dirname -- "$IMAGE_PATH")" && pwd)/$(basename -- "$IMAGE_PATH")"

echo "[$(date '+%F %T')] Starting GAP-MS tutorial run"
echo "Image:  $IMAGE_PATH"
echo "Inputs: $TUTORIAL_DIR"
echo "Output: $OUTPUT_DIR"

apptainer run --cleanenv \
    --bind "$TUTORIAL_DIR:/input:ro" \
    --bind "$OUTPUT_DIR:/output" \
    "$IMAGE_PATH" \
    -g /input/Helixer_prediction.gtf \
    -f /input/Helixer_prediction.faa \
    -b /input/NC_023759.2.test.bam \
    -a /input/NC_023759.2_genome_Sequence.fna \
    -p /input/Helixer_MS_identified_peptides.tsv \
    -rg /input/Reference_annotation.gtf \
    -rf /input/Reference_annotation.faa \
    -o /output

echo "[$(date '+%F %T')] GAP-MS tutorial run completed"
echo "Results: $OUTPUT_DIR"
