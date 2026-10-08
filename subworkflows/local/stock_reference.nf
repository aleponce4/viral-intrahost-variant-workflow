/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: STOCK_REFERENCE

    Turn a de novo consensus for one virus stock into a reference this pipeline
    can be pointed at: a named FASTA with an index, the lab annotation
    transferred onto it, a check that the transfer did not break a reading
    frame, and a coordinate map back to the lab reference.

    The main workflow maps every sample to one reference. Characterising a stock
    needs the opposite: each stock gets its own reference, so that calls against
    it are the variation inside that stock rather than its fixed differences
    from a strain in GenBank. Run this once per stock, then run the main
    workflow with --fasta and --gff pointing at what it produced.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { STOCK_PREPARE                          } from '../../modules/local/stock/prepare/main'
include { LIFTOFF                                } from '../../modules/local/liftoff/main'
include { MINIMAP2_ALIGN as MINIMAP2_ORIENT      } from '../../modules/local/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_LIFTOVER    } from '../../modules/local/minimap2/align/main'
include { STOCK_LIFTOVER_TABLE                   } from '../../modules/local/stock/liftover_table/main'
include { STOCK_CHECK_ANNOTATION                 } from '../../modules/local/stock/check_annotation/main'
include { STOCK_LIFTOVER_BED                     } from '../../modules/local/stock/liftover_bed/main'
include { STOCK_REPORT                           } from '../../modules/local/stock/report/main'

workflow STOCK_REFERENCE {
    take:
    ch_consensus  // channel: [ val(meta), path(consensus_fasta) ]
    ch_ref_fasta  // channel: [ val(meta), path(fasta), path(fai) ]
    ch_ref_gff    // channel: [ val(meta), path(gff) ]
    ch_primer_bed // value:   path(primer bed in reference coordinates), or a NO_FILE placeholder

    main:
    ch_versions = Channel.empty()

    // Which strand did the assembler happen to write the genome on? The answer
    // decides whether the consensus is used as it came or reverse-complemented,
    // so it has to be settled before anything else reads the sequence.
    MINIMAP2_ORIENT(ch_consensus, ch_ref_fasta.first())

    // Orient it, name the contig after the stock, and index it.
    STOCK_PREPARE(ch_consensus, MINIMAP2_ORIENT.out.paf)
    ch_stock_fasta = STOCK_PREPARE.out.fasta

    // Carry the lab annotation across to the stock's own coordinates.
    LIFTOFF(ch_stock_fasta, ch_ref_fasta.first(), ch_ref_gff.first())

    // Refuse to hand on an annotation whose reading frames the transfer broke.
    STOCK_CHECK_ANNOTATION(
        ch_stock_fasta,
        LIFTOFF.out.gff,
        ch_ref_fasta.first(),
        ch_ref_gff.first()
    )

    // Map every position to the lab reference, so results from different stocks
    // can be put on one coordinate system and reference-coordinate intervals
    // such as a primer BED can be placed on the stock. This realigns the
    // oriented sequence, because the first alignment was of the raw consensus.
    MINIMAP2_LIFTOVER(ch_stock_fasta.map { meta, fasta, fai -> [ meta, fasta ] }, ch_ref_fasta.first())
    STOCK_LIFTOVER_TABLE(MINIMAP2_LIFTOVER.out.paf)

    // Amplicon stocks: a primer scheme designed against the lab reference has to
    // be moved onto this stock before ivar trim can use it.
    ch_primer_out = Channel.empty()
    ch_primer_report = Channel.empty()
    ch_primer_for_report = Channel.value(file("${projectDir}/assets/NO_FILE"))
    if (params.stock_primer_bed) {
        STOCK_LIFTOVER_BED(STOCK_LIFTOVER_TABLE.out.table, ch_primer_bed)
        ch_primer_out = STOCK_LIFTOVER_BED.out.bed
        ch_primer_report = STOCK_LIFTOVER_BED.out.report
        ch_primer_for_report = STOCK_LIFTOVER_BED.out.report.map { meta, tsv -> tsv }
        ch_versions = ch_versions.mix(STOCK_LIFTOVER_BED.out.versions)
    }

    // One page per stock: how far it sits from the lab reference, whether the
    // annotation came across intact, and anything a reader should look at.
    STOCK_REPORT(
        ch_stock_fasta,
        STOCK_CHECK_ANNOTATION.out.report,
        STOCK_LIFTOVER_TABLE.out.summary,
        STOCK_PREPARE.out.orientation,
        LIFTOFF.out.unmapped,
        ch_primer_for_report
    )

    ch_versions = ch_versions
        .mix(MINIMAP2_ORIENT.out.versions)
        .mix(STOCK_PREPARE.out.versions)
        .mix(LIFTOFF.out.versions)
        .mix(STOCK_CHECK_ANNOTATION.out.versions)
        .mix(MINIMAP2_LIFTOVER.out.versions)
        .mix(STOCK_LIFTOVER_TABLE.out.versions)
        .mix(STOCK_REPORT.out.versions)

    emit:
    fasta            = ch_stock_fasta                 // [ meta, fasta, fai ]
    gff              = LIFTOFF.out.gff                // [ meta, gff ]
    unmapped         = LIFTOFF.out.unmapped           // [ meta, txt ]
    orientation      = STOCK_PREPARE.out.orientation
    annotation_check = STOCK_CHECK_ANNOTATION.out.report
    liftover         = STOCK_LIFTOVER_TABLE.out.table
    liftover_summary = STOCK_LIFTOVER_TABLE.out.summary
    paf              = MINIMAP2_LIFTOVER.out.paf
    primer_bed       = ch_primer_out       // only when --stock_primer_bed was given
    primer_report    = ch_primer_report
    report_md        = STOCK_REPORT.out.markdown
    report_tsv       = STOCK_REPORT.out.tsv
    versions         = ch_versions
}
