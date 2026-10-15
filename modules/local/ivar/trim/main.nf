process IVAR_TRIM {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_ivar}"

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed

    output:
    tuple val(meta), path("*.trimmed.sorted.bam"), path("*.trimmed.sorted.bam.bai"), emit: bam
    tuple val(meta), path("*.ivar_trim.log")                                       , emit: log
    path "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // -e keeps reads that do not begin at a primer. Without it ivar trim discards them, which
    // is right for a library whose reads start at primer sites and wrong for one that was
    // fragmented after amplification. On a real fragmented amplicon library 61% of reads began
    // outside a primer, so dropping them would have thrown away most of the data.
    //
    // ivar trim counts how many reads it trimmed, how many began outside a primer, and how many
    // it dropped on quality. Those numbers say whether the scheme matches the library, so the
    // log is kept rather than left in the work directory.
    """
    set -o pipefail
    ivar trim -i ${bam} -b ${primer_bed} -p ${prefix}.trimmed -q ${params.ivar_min_bq} -m 30 -e 2>&1 \\
        | tee ${prefix}.ivar_trim.log
    samtools sort -@ ${task.cpus} -o ${prefix}.trimmed.sorted.bam ${prefix}.trimmed.bam
    samtools index ${prefix}.trimmed.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | grep -i version | sed 's/iVar version //' | sed 's/iVar //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.trimmed.sorted.bam
    touch ${prefix}.trimmed.sorted.bam.bai
    touch ${prefix}.ivar_trim.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ivar: \$(ivar version 2>&1 | grep -i version | sed 's/iVar version //' | sed 's/iVar //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
