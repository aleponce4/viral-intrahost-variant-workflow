# Test Fixtures

This directory contains test datasets used by `nf-test` and the `-profile test` execution profile:

- `viral_ref.test.fasta`: Synthetic VEEV-derived sequence (~11.5 kb) with contig `KP282671.1`.
- `viral_ref.test.gff3`: Matching GFF3 containing gene, mRNA, and CDS annotations.
- `sampleA.test_1.fastq.gz` / `sampleA.test_2.fastq.gz`: Paired-end FASTQ reads for sampleA.
- `sampleB.test_1.fastq.gz` / `sampleB.test_2.fastq.gz`: Paired-end FASTQ reads for sampleB.
- `sampleA.test.bam` / `sampleA.test.bam.bai`: Test BAM with seeded variants.
- `sampleB.test.bam` / `sampleB.test.bam.bai`: Test BAM with seeded variants.
- `primers.test.bed`: Synthetic primer BED file for amplicon protocol testing.
- `sampleA.capped.test.vcf`: Two LoFreq-style calls and a header that records `--max-depth 30`. Position 20 has 20 reads in `sampleA.test.depth.tsv`, under the cap. Position 100 has 60, so the cap lets LoFreq examine half of them. Input to `LOFREQ_DEPTH_CHECK`.
- `sampleA.test.ivar.tsv`: iVar `variants` output in iVar's real column layout, holding an SNV, an insertion inside a homopolymer run and an in-frame deletion.

Stock reference fixtures, used by `-entry BUILD_STOCK_REFERENCE`:

- `stock.test.fasta` / `.fai`: The test reference with one substitution (position 9000), a 2 bp insertion and a 3 bp deletion. All three sit after the CDS (50-7549), so the reading-frame checks test the annotation transfer rather than the edits.
- `stock.multi.test.fasta`: The same sequence split in two, so a test can show that a fragmented assembly is refused rather than truncated to its first contig.
- `stock.rc.test.fasta` / `stock.rc.test.paf`: The same stock written on the opposite strand. An assembler picks a strand per run, so a stock reference has to be oriented against the lab reference before anything reads it.
- `stock.test.paf`: `minimap2 -cx asm10 --cs` of the stock against the reference. Input to the coordinate map.
- `stock.test.gff3`: The reference annotation carried onto the stock by Liftoff.

- `stock.liftover.test.tsv`, `stock.liftover_summary.test.tsv`, `stock.annotation_check.test.tsv`, `stock.orientation.test.txt`, `stock.unmapped.test.txt`: Outputs of earlier steps, committed so the liftover-BED and report modules can be tested on their own.

Fixtures are generated deterministically using `tests/data/generate_fixtures.sh`.
