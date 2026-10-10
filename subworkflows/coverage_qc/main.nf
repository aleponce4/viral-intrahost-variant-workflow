/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: COVERAGE_QC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COVERAGE_DEPTH     } from '../../modules/local/coverage/depth/main'
include { COVERAGE_SUMMARIZE } from '../../modules/local/coverage/summarize/main'

workflow COVERAGE_QC {
    take:
    ch_viral_bams // channel: [ val(meta), path(bam), path(bai) ]

    main:
    ch_versions = Channel.empty()
    ch_coverage_summary = Channel.empty()
    ch_depth            = Channel.empty()

    if (params.run_coverage) {
        COVERAGE_DEPTH(ch_viral_bams)
        ch_depth    = COVERAGE_DEPTH.out.depth
        ch_versions = ch_versions.mix(COVERAGE_DEPTH.out.versions)

        COVERAGE_SUMMARIZE(COVERAGE_DEPTH.out.depth)
        ch_coverage_summary = COVERAGE_SUMMARIZE.out.summary
        ch_versions         = ch_versions.mix(COVERAGE_SUMMARIZE.out.versions)
    }

    emit:
    coverage_summary = ch_coverage_summary
    depth            = ch_depth // channel: [ val(meta), path(depth.tsv) ]
    versions         = ch_versions
}
