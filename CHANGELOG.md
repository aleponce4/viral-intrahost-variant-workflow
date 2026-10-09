# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Strand-bias rule for `lofreq filter` is configurable (`--lofreq_sb_thresh`, default 0)**:
  `LOFREQ_FILTER` passes `--sb-thresh` to `lofreq filter` when the value is above 0. That replaces
  LoFreq's built-in strand-bias rule, an FDR test, with a fixed phred cutoff. The default, `0`,
  keeps the built-in rule, so results are the same as before this parameter existed. Either rule
  drops a call only when about 85% of its alt reads also sit on one strand. On the TC-83 stock
  data (full depth, 2 million pairs and 200,000 pairs) neither the built-in rule nor any cutoff
  from 100 down to 40 removed a call, so a fixed cutoff is not needed there. High
  strand-bias values are common there and grow with depth (4 of 124 calls above 60 at 200,000 pairs,
  145 of 516 at 2 million, 1,156 of 2,369 at full depth), but none of those calls had 85% of its alt
  reads on one strand. The rule in use is written to `<sample>.qc_stats.txt`.
- **A primer scheme can be moved onto a stock (`--stock_primer_bed`, `bin/liftover_bed.py`)**:
  a scheme is designed once against the lab reference, but `ivar trim` needs the coordinates of
  whatever the reads were aligned to, and one indel between a stock and the reference shifts
  every primer after it. `BUILD_STOCK_REFERENCE` now takes a BED in reference coordinates and
  writes `<stock>.primers.bed` in the stock's, with `<stock>.primers_liftover.tsv` recording what
  happened to each interval. A primer whose site the stock has deleted is dropped, not moved
  somewhere plausible, because an amplicon that cannot be trimmed must not be trimmed wrongly.
  The script works in either direction, so a result called on a stock can also be put back on
  reference coordinates.
- **One report per stock (`bin/summarize_stock_reference.py`)**: `<stock>.stock_report.md` and
  `.tsv` gather the sequence, the distance from the lab reference, the per-CDS annotation check,
  the primer placement and a short list of anything worth a look. It reports rather than decides;
  the annotation check has already failed the run by this point if the transfer broke a frame.
- **A stock consensus is oriented against the lab reference before anything reads it**
  (`STOCK_PREPARE`): an assembler has no way to know which strand a genome is meant to be read
  on and picks one per run. Two runs of the same assembly pipeline over the same stock reads gave
  opposite orientations, so roughly half of all stocks would otherwise get a reference whose genes
  sit on the minus strand and whose coordinates run backwards against every other stock.
  `MINIMAP2_ALIGN` now runs first on the raw consensus, `STOCK_PREPARE` takes the strand from that
  alignment and reverse-complements when it is `-`, and the coordinate map is built from a second
  alignment of the oriented sequence. The chosen orientation is published as
  `<stock>.orientation.txt`. `MINIMAP2_ALIGN` no longer requires an index for its query, since the
  query is not indexed until after it has been oriented.
- **Per-stock references (`-entry BUILD_STOCK_REFERENCE`)**: turns a de novo consensus for one
  virus stock into the `--fasta` and `--gff` of a normal run, so calls describe the variation
  inside that stock instead of its fixed differences from a strain in GenBank.
  - `STOCK_PREPARE` names the contig after the stock and indexes it, and refuses a multi-contig
    assembly rather than silently using its first contig.
  - `LIFTOFF` carries the lab annotation onto the stock. `MINIMAP2_ALIGN` aligns the two with
    `-x asm10 --cs`.
  - `bin/check_lifted_annotation.py` checks each CDS for coordinates inside the contig, a length
    that is a multiple of three, and no stop codon before the end. The reference is passed as a
    control: a problem it already has is reported as `inherited` and does not stop the run, so
    the nsP3 opal (TGA) readthrough codon present in the lab's TC-83, VEEV INH-9813 and EEEV V105
    references is not mistaken for a broken transfer. Start and stop codons are reported but
    never fail, because alphavirus mature peptides are polyprotein cleavage products.
  - `bin/build_liftover_table.py` walks the PAF `cs` tag and writes one row per alignment column,
    so a variant called on a stock can be reported at its reference position and a
    reference-coordinate primer BED can be placed on the stock. It refuses a reverse-strand or
    multi-block alignment rather than mapping coordinates that cannot be mapped.
  - New parameters `--stock_consensus` and `--stock_name`.

- **Results on a stock can be reported on lab-reference coordinates (`--liftover_tsv`,
  `bin/append_ref_coords.py`)**: a call on a stock carries a stock position, which is not
  comparable across stocks and is not the position the lab's other datasets use. Pass
  `qc/<stock>.liftover.tsv` from `BUILD_STOCK_REFERENCE` and every annotated VCF is also
  written as `<sample>.csq.refcoords.vcf` with INFO `REF_CONTIG`, `REF_POS`, `REF_END` and
  `REF_STATUS`, and the variant summary table gains `ref_contig`, `ref_pos`, `ref_end` and
  `ref_status` columns. The call itself is not changed. A stock base the reference lacks has
  status `insertion` and no `REF_POS`, and a position outside the aligned block is
  `unmapped`. Nothing runs unless the parameter is set. On the TC-83 pilot, all 516 LoFreq and
  1,981 iVar calls that received a position carried the same REF base as the lab reference
  at that position; 6 iVar poly-A calls past the aligned block were `unmapped`.

- **The Stage A recipe is in the repository (`assets/stage_a_viralmetagenome.params.yaml`,
  README "Stage A")**: the stock workflow starts from a de novo consensus built by
  nf-core/viralmetagenome 1.2.0, which this repository does not contain. The settings that
  decide whether that consensus is usable (`deduplicate: false`, `normalise_reads: true`, the
  `skip_*` flags) existed only in a local run script. They are now a params file with the
  reason for each setting in the README. Documentation and one data file; no code changes.

- **One command for the whole stock workflow (`scripts/run_stock_workflow.sh`)**: chains
  Stage A (nf-core/viralmetagenome 1.2.0), `BUILD_STOCK_REFERENCE` and the main workflow for
  one stock. They stay three separate Nextflow runs, each in its own directory, joined by
  files, so each keeps its own work directory and `-resume`. The script finds the consensus
  Stage A wrote, hands it on, and passes the stock reference and its liftover table to the
  variant run. Stage A is optional: `--consensus` skips it, and `--step` runs one step. The
  other options are `--primer-bed`, `--profile`, `--config` and `--dry-run`. It writes the
  commands and the Nextflow and Java versions to `stock_run.log`. `NXF_STAGE_A` and `NXF_STAGE_B`
  pick the Nextflow binary per stage, because the pilot ran the two under different versions.
  `NXF_STAGE_A_JAVA_HOME` and `NXF_STAGE_B_JAVA_HOME` do the same for Java, because Stage A
  needs Java 17 or newer and the rest of the pipeline does not.
- **`apptainer` and `singularity` profiles**: Apptainer used to need a local override config
  on top of `-profile docker`. `-profile apptainer` and `-profile singularity` now work for
  every entry point, and CI parses both.

- **The stock workflow is in the diagrams**: the README overview and `ARCHITECTURE.md` now show
  Stage A, `BUILD_STOCK_REFERENCE` with its steps, the launcher, and where the reference
  coordinates and the primer BED attach. Documentation only. The diagrams were checked with
  Mermaid's own parser.

- **The launcher can cut Stage A's reads (`--stage-a-pairs N`) and takes a config per stage
  (`--stage-a-config`, `--variants-config`)**: Stage A maps the reads again in three rounds, so
  at full depth it is estimated to add several hours, and the consensus does not change (200,000
  and 2,000,000 pairs of one stock gave the identical consensus). `--stage-a-pairs` keeps the first
  N pairs of `--stock-reads` for Stage A and leaves `--samples` at full depth. It refuses a value
  that is not a positive whole number, a file that is not a gzip FASTQ, and an R1 and R2 of
  different length. The two config options exist because `--config` reaches every step and Stage A
  has processes named like ours (`IVAR_VARIANTS`, `SAMTOOLS_SORT`), so a setting aimed at the
  variant run would also change Stage A.

### Fixed
- **iVar indels are valid VCF and left-aligned before annotation (`bin/ivar_variants_to_vcf.py`
  1.2.0, `BCFTOOLS_CSQ`)**:
  iVar writes an insertion as ALT `+T` and a deletion as ALT `-AC`. The converter copied that
  into the VCF ALT column, which is not valid VCF, so `bcftools csq` could not annotate the
  record. Both now use VCF alleles anchored on the reference base (`T` to `TT`, `AC` to `A`).
  - iVar reports an insertion inside a homopolymer at the last base of the run, not the
    first. On 22 of 22 real files with an insertion, `bcftools norm` moved the record. `BCFTOOLS_CSQ`
    now adds the `##contig` lines from the FASTA index (iVar and LoFreq VCFs have none), left-aligns
    with `bcftools norm`, then annotates. SNV records, including every LoFreq SNV, are unchanged.
  - `tests/subworkflows/annotation` now covers an SNV, a homopolymer insertion and an in-frame
    deletion, and asserts the left-aligned position and the consequence.
- **LoFreq `INFO/AF` is not a VAF; iVar and LoFreq are given the same reads**:
  - LoFreq's `AF` divides an alt count filtered at `--min-bq` by a depth with no base-quality
    filter, so it reads low and tiered variants could land in the wrong band. Reports and
    plots now compute the VAF from `DP4`. The columns are renamed `af_*` to `vaf_*`, and the
    caller's own value is kept as `lofreq_info_af`. VCFs without `DP4` (the iVar route) keep `AF`,
    which is a real base-count ratio there.
  - iVar's `samtools mpileup` had no mapping-quality floor and left BAQ on. It now takes
    `-B` and `-q ${params.ivar_min_mq}` (default 20, the same as `lofreq_min_mq`).
- **iVar TSV parsed by column name, not position (`bin/ivar_variants_to_vcf.py`, 1.1.0)**:
  the parser read `ALT_DP` from the 6th column and `ALT_FREQ` from the 10th. Real iVar output
  has `REF_DP REF_RV REF_QUAL ALT_DP ALT_RV ALT_QUAL ALT_FREQ TOTAL_DP PVAL PASS`, so `ALT_DP`
  is the 8th column and `ALT_FREQ` the 11th. Every record carried the wrong values: `ALT_QUAL`
  was reported as `AF`, and the `ALT_FREQ` text was read as `TOTAL_DP`, which `isdigit()`
  turned into depth 0. The columns are now looked up by header name.
  - Depths of a million or more, which iVar prints as `1.16389e+06`, were also turned into 0.
    They are now parsed as numbers.
  - A missing column, a short row or an unparseable value now stops the run with the file and
    line, instead of becoming 0 or 1.0 without a message.
  - `tests/data/sampleA.test.ivar.tsv` used an invented column layout, which is why the old
    indices passed CI. It now uses iVar's real header, and `tests/unit/` pins the layout. CI
    runs the new unit tests.
- **Fabricated report output removed (reporting integrity)**:
  - `bin/generate_run_summary.py` no longer emits a hardcoded two-row `sampleA`/`sampleB`
    `COMPLETED` table. It now derives one row per sample from the QC artefacts actually
    staged into the task directory (`*.qc_stats.txt` from `LOFREQ_FILTER`, merged with
    `*coverage_summary.tsv` from `COVERAGE_SUMMARIZE`). With no staged inputs it writes a
    header-only TSV and warns on stderr instead of inventing samples.
  - `bin/generate_coverage_plots.py` no longer writes a hardcoded 1×1-pixel PNG. It parses
    the staged `samtools depth` TSVs and renders a real matplotlib figure: per-sample
    log-scaled depth profile with 100×/1,000×/5,000× tier reference lines, plus a mean-depth
    bar panel. Per-sample statistics are computed with `summarize_coverage.summarize_depth`
    so figure and TSV numbers agree. An empty run yields an explicitly labelled
    "no coverage data available" figure rather than a placeholder image.
  - `bin/plot_coverage.py` is no longer a no-op stub; it is a thin wrapper over the single
    coverage-plot implementation, exposing the `--coverage-dir` / `--output-dir` spelling.
  - `modules/local/report/run_summary` and `modules/local/report/coverage_plots` were passing
    flags the scripts do not accept (`--output-tsv`) or omitting `--dataset`; their commands
    now match the real CLIs.

### Changed
- **Documentation**: Replaced the internal `docs/history/` build specifications with
  `docs/migration.md`, a human-facing account of the legacy Bash → Nextflow DSL2 migration,
  linked from the README.
- **Third-party attribution**: Added `THIRD_PARTY_NOTICES.md` recording that
  `assets/snpgenie/snpgenie.pl` is SNPGenie by Chase W. Nelson, retained under GPLv3 under
  its own terms (with upstream URL, pinned commit, and citation). The project's own code
  remains MIT-licensed. Added a "Third-party Components" section to the README.

### Removed
- **Study data**: Deleted all committed real study outputs under
  `legacy/Variant_discovery_pipeline/Analysis_output/` (per-position intra-host allele
  frequency CSVs and the derived TIFF figure). Cleared all notebook outputs from
  `legacy/Variant_discovery_pipeline/Scripts/analysis.ipynb` and replaced leaked absolute
  paths and internal specimen identifiers throughout `legacy/` with relative or clearly
  invented placeholder values. Legacy docs now reference the synthetic fixtures in
  `tests/data/`.

## [1.3.0] - 2026-07-29

### Fixed
- **Scientific Audit Adjustments (RA-1 to RA-12)**:
  - **iVar Frequency Floor (RA-1)**: Lowered `ivar_min_freq` default from `0.01` (1%) to `0.001` (0.1%) to align reporting with 0.1%–1% iSNV discovery targets.
  - **Coverage Inclusion Floors (RA-2)**: Lowered inclusion depth floors `ivar_min_depth` and `lofreq_min_depth` from `1000` to `10` to eliminate artificial dead zones at amplicon shoulders and prevent downward bias in SNPGenie diversity denominators.
  - **Dead Parameter Cleanup (RA-3)**: Removed dead `lofreq_min_freq` parameter from `nextflow.config` and schema.
  - **LoFreq Mapping Quality Floor (RA-4)**: Relaxed `lofreq_min_mq` from `60` to `20` to avoid stripping informative reads from naturally variable/repeated viral regions.
  - **Haplotype Parameter Explicit Geometry (RA-5 & RA-6)**: Added explicit parameters `viloca_window` (`150`), `viloca_shift` (`50`), and `cliquesnv_min_freq` (`0.001`), wiring them into process invocations and recording them in `methods_key_parameters.tsv`.
  - **Reporting Container & Hard-Fail Statistical Guard (RA-7)**: Updated `container_python_reporting` to a multi-package data science biocontainer with pandas and scipy, and updated `analyze_delta_selection.py` to hard-fail (`exit 1`) if scipy is unavailable.
  - **iVar VCF Output Integrity (RA-8)**: Corrected fabricated fields in `ivar_variants_to_vcf.py` by emitting haploid GT (`1`) and converting binomial p-values to Phred-scaled QUAL (`-10*log10(pval)`).
  - **Proportion to Percent Reporting Contract (RA-9)**: Standardized raw machine files to proportions (0–1) and human-facing summary tables/plots to percentages (`100 × proportion`) with explicit `%` labels.
  - **Frequency Tiering Policy (RA-10)**: Documented 3-tier iSNV classification policy (`≥1%` high-confidence candidate, `0.1–1%` low-frequency candidate, `<0.1%` exploratory) in documentation and executive reporting.
  - **Explicit SNPGenie Frequency Floor (RA-11)**: Passed `--minfreq=0` explicitly to `snpgenie.pl` in `SNPGENIE_RUN`.
  - **Indel Scope Documentation (RA-12)**: Documented SNV-only indel scope under default `lofreq_enable_indelqual = false`.
  - **VILOCA Module Fix & Haplotype Reporting Integration**: Updated VILOCA module to execute real `viloca run` using pinned `quay.io/biocontainers/viloca:1.1.1--py310h563914a_0` container; split HAPLOTYPE subworkflow emits into dedicated CliqueSNV and VILOCA channels; added `build_haplotype_tables.py`, `summarize_linked_mutations.py`, and `plot_haplotypes.py` (composition stacked bar plot + Minimum Spanning Network with treatment pie nodes); integrated reporting into `REPORTING` subworkflow and `assets/executive_report.qmd`.

## [1.2.0] - 2026-07-29

### Added
- **Executive HTML Reporting**: Added self-contained HTML executive summary report module (`EXECUTIVE_REPORT`) rendering key metrics, depth profiles, variant density, selection statistics, and software provenance in a single report artifact.
- **Public & Hermetic Test Datasets**: Standardized `-profile test` for zero-friction demo execution, while retaining offline hermetic `nf-test` suite.

## [1.1.0] - 2026-07-29

### Added
- **Production Hardening & Schema Validation**: Integrated `plugin/nf-schema@2.1.1` for JSON Schema Draft 2020-12 input validation and fast parameter checking.
- **FASTQ Ingress & Read Preprocessing**: Added `READ_PREPROCESSING` subworkflow (`FASTQC`, `BWA_INDEX`, `BWA_MEM`, `SAMTOOLS_STATS`, `IVAR_TRIM`) for raw FASTQ input handling.
- **Subworkflow Architecture & Modularization**: Relocated `selection` and `haplotype` subworkflows to `subworkflows/local/` and wired explicit treatment group channels.
- **Software Provenance & MultiQC**: Added `DUMP_SOFTWARE_VERSIONS` module and `MULTIQC` terminal dashboard process.

## [1.0.0] - 2026-07-29

### Added
- Initial release of the `alphavirus-variant-analysis` Nextflow DSL2 workflow.
