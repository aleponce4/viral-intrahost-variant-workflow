#!/usr/bin/env bash
set -euo pipefail

python3 tests/data/generate_fixtures.py

# Convert SAM to sorted BAM and index using docker samtools
CONTAINER_SAMTOOLS="quay.io/biocontainers/samtools:1.21--h50ea8bc_0"

for sample in sampleA sampleB; do
    samtools view -bS tests/data/${sample}.sam | samtools sort -o tests/data/${sample}.test.bam -
    samtools index tests/data/${sample}.test.bam
    rm -f tests/data/${sample}.sam
done

# Stock reference fixtures: index the stock, align it to the reference for the
# coordinate map, and transfer the annotation the way the pipeline does.
CONTAINER_MINIMAP2="quay.io/biocontainers/minimap2:2.30--h577a1d6_0"
CONTAINER_LIFTOFF="quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_2"

samtools faidx tests/data/stock.test.fasta

minimap2 -cx asm10 --cs tests/data/viral_ref.test.fasta tests/data/stock.test.fasta \
    > tests/data/stock.test.paf
minimap2 -cx asm10 --cs tests/data/viral_ref.test.fasta tests/data/stock.rc.test.fasta \
    > tests/data/stock.rc.test.paf

liftoff -g tests/data/viral_ref.test.gff3 -o tests/data/stock.test.gff3 \
    -u tests/data/stock.unmapped.txt -dir tests/data/liftoff_intermediate \
    tests/data/stock.test.fasta tests/data/viral_ref.test.fasta
rm -rf tests/data/liftoff_intermediate tests/data/stock.unmapped.txt \
       tests/data/viral_ref.test.gff3_db

# Outputs of earlier steps, committed so the liftover-BED and report modules can
# be tested without first running the steps that produce their inputs.
python3 bin/build_liftover_table.py \
    --paf tests/data/stock.test.paf \
    --output-tsv tests/data/stock.liftover.test.tsv \
    --summary-tsv tests/data/stock.liftover_summary.test.tsv

python3 bin/check_lifted_annotation.py \
    --fasta tests/data/stock.test.fasta \
    --gff tests/data/stock.test.gff3 \
    --reference-fasta tests/data/viral_ref.test.fasta \
    --reference-gff tests/data/viral_ref.test.gff3 \
    --output-tsv tests/data/stock.annotation_check.test.tsv

printf 'as assembled\n' > tests/data/stock.orientation.test.txt
: > tests/data/stock.unmapped.test.txt

echo "Fixtures successfully generated."
