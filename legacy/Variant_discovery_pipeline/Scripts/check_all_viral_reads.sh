#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate lofreq-env

# Point BAM_DIR at your own indexed BAM directory (override via environment).
BAM_DIR="${BAM_DIR:-./input/${DATASET:-mouse_veev}/BAMs}"
CONTIG="KP282671.1"

echo "Sample | Viral reads | Mean depth"
echo "-------|-------------|-----------"
for bam in "$BAM_DIR"/*.bam; do
    s=$(basename "$bam" .bam)
    reads=$(samtools view -c -r "$CONTIG" "$bam" 2>/dev/null || echo 0)
    if [ "$reads" -gt 0 ]; then
        mean=$(samtools depth -a -r "$CONTIG" "$bam" 2>/dev/null | awk '{sum+=$3; n++} END{if(n>0) printf "%.0f", sum/n; else print 0}')
    else
        mean=0
    fi
    echo "$s | $reads | $mean"
done
