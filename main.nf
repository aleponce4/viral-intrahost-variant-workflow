#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    alphavirus-variant-analysis Main Workflow Entrypoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { INPUT_CHECK            } from './subworkflows/input_check/main'
include { READ_PREPROCESSING     } from './subworkflows/local/read_preprocessing'
include { VARIANT_CALLING        } from './subworkflows/variant_calling/main'
include { ANNOTATION             } from './subworkflows/annotation/main'
include { COVERAGE_QC            } from './subworkflows/coverage_qc/main'
include { LOFREQ_DEPTH_CHECK     } from './modules/local/lofreq/depth_check/main'
include { SELECTION              } from './subworkflows/local/selection'
include { HAPLOTYPE              } from './subworkflows/local/haplotype'
include { REPORTING              } from './subworkflows/reporting/main'
include { DUMP_SOFTWARE_VERSIONS } from './modules/local/dumpsoftwareversions/main'
include { MULTIQC                } from './modules/local/multiqc/main'
include { SAMTOOLS_FAIDX        } from './modules/local/samtools/faidx/main'
include { STOCK_REFERENCE       } from './subworkflows/local/stock_reference'

/*
 * Build a per-stock reference from a de novo consensus:
 *
 *   nextflow run . -entry BUILD_STOCK_REFERENCE -profile docker \
 *       --stock_consensus stock.fasta --stock_name TC83-stock \
 *       --fasta lab_reference.fasta --gff lab_reference.gff3 --outdir refs/
 *
 * The outputs are the --fasta and --gff for a normal run of this pipeline over
 * that stock's samples. This entry does not read a samplesheet, so it checks its
 * own parameters instead of calling validateParameters(), whose required list
 * covers the main workflow.
 */
workflow BUILD_STOCK_REFERENCE {
    if (!params.stock_consensus) {
        error "Parameter --stock_consensus must be specified (the de novo consensus FASTA for one stock)."
    }
    if (!params.fasta) {
        error "Parameter --fasta must be specified (the lab reference to transfer the annotation from)."
    }
    if (!params.gff) {
        error "Parameter --gff must be specified (the lab reference annotation, with CDS features)."
    }

    def stock_id = params.stock_name ?: file(params.stock_consensus).getBaseName()

    ch_consensus = Channel.of([ [ id: stock_id ], file(params.stock_consensus, checkIfExists: true) ])
    ch_ref_gff   = Channel.of([ [ id: 'reference' ], file(params.gff, checkIfExists: true) ])

    SAMTOOLS_FAIDX(
        Channel.of([ [ id: 'reference' ], file(params.fasta, checkIfExists: true) ])
    )

    // Optional: a primer scheme in lab-reference coordinates, moved onto the stock.
    ch_primer_bed = params.stock_primer_bed
        ? Channel.value(file(params.stock_primer_bed, checkIfExists: true))
        : Channel.value(file("${projectDir}/assets/NO_FILE"))

    STOCK_REFERENCE(ch_consensus, SAMTOOLS_FAIDX.out.fai, ch_ref_gff, ch_primer_bed)

    DUMP_SOFTWARE_VERSIONS(
        STOCK_REFERENCE.out.versions
            .mix(SAMTOOLS_FAIDX.out.versions)
            .unique()
            .collectFile(name: 'collated_versions.yml', newLine: true)
    )
}

workflow {
    validateParameters()
    log.info paramsSummaryLog(workflow)

    log.info """
    ================================================================
    A L P H A V I R U S   V A R I A N T   A N A L Y S I S
    ================================================================
    Dataset         : ${params.dataset}
    Input           : ${params.input}
    FASTA           : ${params.fasta}
    GFF             : ${params.gff}
    Contig          : ${params.viral_contig}
    Protocol        : ${params.protocol}
    Primer BED      : ${params.primer_bed}
    Outdir          : ${params.outdir}
    Run iVar        : ${params.run_ivar}
    Run LoFreq      : ${params.run_lofreq}
    Run Annotation  : ${params.run_annotation}
    Run Coverage    : ${params.run_coverage}
    Run SNPGenie    : ${params.run_snpgenie}
    Run Haplotype   : ${params.run_haplotype}
    ================================================================
    """

    // 1. Input Validation
    INPUT_CHECK(file(params.input, checkIfExists: true))

    // 2. Read Preprocessing & Alignment
    READ_PREPROCESSING(INPUT_CHECK.out.fastqs, INPUT_CHECK.out.fasta)

    // 3. Variant Calling
    VARIANT_CALLING(READ_PREPROCESSING.out.bams, INPUT_CHECK.out.fasta, INPUT_CHECK.out.gff)

    // 4. Annotation
    ANNOTATION(VARIANT_CALLING.out.ivar_tsv, VARIANT_CALLING.out.lofreq_vcf, INPUT_CHECK.out.fasta, INPUT_CHECK.out.gff)

    // 5. Coverage QC
    COVERAGE_QC(VARIANT_CALLING.out.viral_bams)

    // 5b. Flag LoFreq calls made where LoFreq counted only some of the reads. Runs
    // when both the LoFreq calls and the depth table exist.
    ch_lofreq_depth = VARIANT_CALLING.out.lofreq_vcf
        .map { meta, vcf, tbi -> [ meta.id, meta, vcf ] }
        .join(COVERAGE_QC.out.depth.map { meta, depth -> [ meta.id, depth ] })
        .map { id, meta, vcf, depth -> [ meta, vcf, depth ] }
    LOFREQ_DEPTH_CHECK(ch_lofreq_depth)

    // 6. Selection Analysis (optional, --run_snpgenie)
    ch_snpgenie_vcfs = VARIANT_CALLING.out.lofreq_vcf
        .map { meta, vcf, tbi -> [ meta, vcf ] }

    ch_treatment_groups = INPUT_CHECK.out.fastqs
        .map { meta, fastq_1, fastq_2 -> [ meta.treatment, meta.id ] }
        .groupTuple()
        .map { treatment, ids -> [ treatment: treatment, n_replicates: ids.size() ] }

    SELECTION(
        ch_snpgenie_vcfs,
        INPUT_CHECK.out.fasta,
        INPUT_CHECK.out.gff,
        file(params.input),
        ch_treatment_groups
    )

    // 7. Haplotype Reconstruction (optional, --run_haplotype)
    HAPLOTYPE(VARIANT_CALLING.out.viral_bams, INPUT_CHECK.out.fasta)

    // 8. Software Versions Aggregation
    ch_versions = Channel.empty()
        .mix(INPUT_CHECK.out.versions)
        .mix(READ_PREPROCESSING.out.versions)
        .mix(VARIANT_CALLING.out.versions)
        .mix(ANNOTATION.out.versions)
        .mix(COVERAGE_QC.out.versions)
        .mix(LOFREQ_DEPTH_CHECK.out.versions)
        .mix(SELECTION.out.versions)
        .mix(HAPLOTYPE.out.versions)

    DUMP_SOFTWARE_VERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml', newLine: true))

    // 9. Executive HTML & Biological Reporting
    REPORTING(
        VARIANT_CALLING.out.lofreq_qc.map { meta, qc -> qc },
        COVERAGE_QC.out.coverage_summary.map { meta, summary -> summary },
        VARIANT_CALLING.out.lofreq_vcf.map { meta, vcf, tbi -> vcf },
        SELECTION.out.selection_tables,
        HAPLOTYPE.out.cliquesnv_fasta,
        HAPLOTYPE.out.viloca_cooccurrence,
        INPUT_CHECK.out.fasta,
        INPUT_CHECK.out.gff,
        file(params.input),
        DUMP_SOFTWARE_VERSIONS.out.yml
    )

    // 10. MultiQC Terminal Dashboard
    ch_multiqc_files = Channel.empty()
        .mix(READ_PREPROCESSING.out.fastqc_zip.map { meta, zip -> zip })
        .mix(READ_PREPROCESSING.out.bwa_log.map    { meta, log -> log })
        .mix(READ_PREPROCESSING.out.stats.map      { meta, stats, flagstat, idxstats -> [ stats, flagstat, idxstats ] })
        .mix(VARIANT_CALLING.out.lofreq_qc.map     { meta, qc -> qc })
        .mix(DUMP_SOFTWARE_VERSIONS.out.mqc_yml)
        .collect()

    MULTIQC(ch_multiqc_files, file("${projectDir}/assets/multiqc_config.yml"))
}
