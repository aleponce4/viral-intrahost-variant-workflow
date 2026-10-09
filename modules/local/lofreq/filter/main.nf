process LOFREQ_FILTER {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_lofreq}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.variants.filtered.vcf.gz"), path("*.variants.filtered.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.qc_stats.txt")                                                    , emit: qc_stats
    path "versions.yml"                                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // LoFreq reads 0 as a real phred cutoff rather than "off", so the flag is
    // dropped entirely to disable the filter. It also conflicts with --sb-mtc.
    def sb_arg = params.lofreq_sb_thresh ? "--sb-thresh ${params.lofreq_sb_thresh}" : ''
    """
    lofreq filter -i ${vcf} -o ${prefix}.variants.filtered.vcf --snvqual-thresh 20 --indelqual-thresh 20 ${sb_arg}

    RAW_COUNT=\$(grep -v "^#" ${vcf} | wc -l || echo 0)
    FILTERED_COUNT=\$(grep -v "^#" ${prefix}.variants.filtered.vcf | wc -l || echo 0)

    bgzip -f ${prefix}.variants.filtered.vcf
    tabix -f -p vcf ${prefix}.variants.filtered.vcf.gz

    echo "Sample: ${prefix}" > ${prefix}.qc_stats.txt
    echo "Contig: ${params.viral_contig}" >> ${prefix}.qc_stats.txt
    echo "Strand-bias phred threshold: ${params.lofreq_sb_thresh ?: 'disabled'}" >> ${prefix}.qc_stats.txt
    echo "Raw variants: \${RAW_COUNT}" >> ${prefix}.qc_stats.txt
    echo "Filtered variants: \${FILTERED_COUNT}" >> ${prefix}.qc_stats.txt
    echo "Status: COMPLETE" >> ${prefix}.qc_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(lofreq version | grep version | sed 's/version: //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variants.filtered.vcf.gz
    touch ${prefix}.variants.filtered.vcf.gz.tbi
    touch ${prefix}.qc_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(lofreq version | grep version | sed 's/version: //')
    END_VERSIONS
    """
}
