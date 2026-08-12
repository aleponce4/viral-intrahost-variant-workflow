#!/bin/bash
# Fast depth check using idxstats + targeted depth on viral contig only
source ~/miniconda3/etc/profile.d/conda.sh
conda activate lofreq-env

# Point BAM_DIR at your own indexed BAM directory (override via environment).
BAM_DIR="${BAM_DIR:-./input/${DATASET:-mouse_veev}/BAMs}"
CONTIG="KP282671.1"

# Check depth on a few high-viral-load samples
for s in ${SAMPLES:-sampleA sampleB}; do
    bam="$BAM_DIR/${s}.bam"
    echo "=== $s ==="
    samtools depth -a -r "$CONTIG" "$bam" | awk '
        BEGIN { min=999999; max=0; sum=0; n=0; a100=0; a1000=0; a5000=0 }
        { d=$3; sum+=d; n++; if(d<min)min=d; if(d>max)max=d; if(d>=100)a100++; if(d>=1000)a1000++; if(d>=5000)a5000++ }
        END { printf "  mean=%.0f  min=%d  max=%d  >=100x=%.1f%%  >=1000x=%.1f%%  >=5000x=%.1f%%\n", sum/n, min, max, a100/n*100, a1000/n*100, a5000/n*100 }
    '
done
