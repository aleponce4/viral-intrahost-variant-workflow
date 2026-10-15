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
- `viral_ref.no_transcript.test.gff3`: `viral_ref.test.gff3` with the `mRNA` line removed and the CDS pointing at the gene. The shape `bcftools csq` cannot use, and the shape real annotations arrive in.
- `viral_ref.other_contig.test.gff3`: a sound annotation of a contig the test FASTA does not hold.
- `viral_ref.two_contigs.test.gff3`: annotations for the test contig and for one more, as a segmented genome's annotation looks. It has to be accepted.
- `purity.test.bam` / `.bai`: Reads mapped to two contigs, for `PURITY_BAM_STATS`. Contig `refA` has 4 reads at mapping quality 60 and 2 at 0, `refB` has 2, a secondary and a supplementary alignment are there to be left out, and 5 reads did not map: 3 plain, 1 poly-G and 1 with an Illumina adapter.
- `purity.stock_counts.test.tsv`, `purity.unmapped_kinds.test.tsv`, `purity.pool_counts.test.tsv`, `purity.kraken2_report.test.txt`: Counts and a Kraken2 report in the formats `summarize_purity.py` reads. The numbers are in `tests/unit/test_summarize_purity.py`.
- `stock.purity.test.tsv` / `.md`: `summarize_purity.py` run on the four files above. It flags a second pool genome and a species. Input to `STOCK_REPORT`.
- `kraken2_db.test/`: Not a database. A file named `hash.k2d`, which is what `KRAKEN2_CLASSIFY` looks for before it runs, and an `inspect.txt` naming three domains, which is the file a real Kraken2 download ships and the module reads to report what the database holds. Together they let the module be tested without a real database.
- `purity.kraken_db_info.test.tsv`: What `KRAKEN2_CLASSIFY` writes about the database it used. Input to `PURITY_SUMMARIZE`.
- `empty.test.fastq.gz`: A gzip file with no reads.
- `stock.offtarget.test.fasta` / `.fai`: The second contig of `stock.multi.test.fasta`. The `sampleA` reads come from the first contig, so none of them maps here. It stands for a stock whose reads are of another virus.
- `sampleA.test.ivar.tsv`: iVar `variants` output in iVar's real column layout, holding an SNV, an insertion inside a homopolymer run and an in-frame deletion.

Stock reference fixtures, used by `-entry BUILD_STOCK_REFERENCE`:

- `stock.test.fasta` / `.fai`: The test reference with one substitution (position 9000), a 2 bp insertion and a 3 bp deletion. All three sit after the CDS (50-7549), so the reading-frame checks test the annotation transfer rather than the edits.
- `stock.multi.test.fasta`: The same sequence split in two, so a test can show that a fragmented assembly is refused rather than truncated to its first contig.
- `stock.rc.test.fasta` / `stock.rc.test.paf`: The same stock written on the opposite strand. An assembler picks a strand per run, so a stock reference has to be oriented against the lab reference before anything reads it.
- `stock.test.paf`: `minimap2 -cx asm10 --cs` of the stock against the reference. Input to the coordinate map.
- `stock.test.gff3`: The reference annotation carried onto the stock by Liftoff.

- `stock.liftover.test.tsv`, `stock.liftover_summary.test.tsv`, `stock.annotation_check.test.tsv`, `stock.orientation.test.txt`, `stock.unmapped.test.txt`: Outputs of earlier steps, committed so the liftover-BED and report modules can be tested on their own.

Fixtures are generated deterministically using `tests/data/generate_fixtures.sh`.
