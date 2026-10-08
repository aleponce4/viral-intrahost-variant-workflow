# Test Fixtures

This directory contains test datasets used by `nf-test` and the `-profile test` execution profile:

- `viral_ref.test.fasta`: Synthetic VEEV-derived sequence (~11.5 kb) with contig `KP282671.1`.
- `viral_ref.test.gff3`: Matching GFF3 containing gene, mRNA, and CDS annotations.
- `sampleA.test_1.fastq.gz` / `sampleA.test_2.fastq.gz`: Paired-end FASTQ reads for sampleA.
- `sampleB.test_1.fastq.gz` / `sampleB.test_2.fastq.gz`: Paired-end FASTQ reads for sampleB.
- `sampleA.test.bam` / `sampleA.test.bam.bai`: Test BAM with seeded variants.
- `sampleB.test.bam` / `sampleB.test.bam.bai`: Test BAM with seeded variants.
- `primers.test.bed`: Synthetic primer BED file for amplicon protocol testing.
- `sampleA.test.ivar.tsv`: iVar `variants` output in iVar's real column layout, holding an SNV, an insertion inside a homopolymer run and an in-frame deletion.

Stock reference fixtures, used by `-entry BUILD_STOCK_REFERENCE`:

- `stock.test.fasta` / `.fai`: The test reference with one substitution (position 9000), a 2 bp insertion and a 3 bp deletion. All three sit after the CDS (50-7549), so the reading-frame checks test the annotation transfer rather than the edits.
- `stock.multi.test.fasta`: The same sequence split in two, so a test can show that a fragmented assembly is refused rather than truncated to its first contig.
- `stock.test.paf`: `minimap2 -cx asm10 --cs` of the stock against the reference. Input to the coordinate map.
- `stock.test.gff3`: The reference annotation carried onto the stock by Liftoff.

Fixtures are generated deterministically using `tests/data/generate_fixtures.sh`.
