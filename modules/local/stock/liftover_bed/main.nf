process STOCK_LIFTOVER_BED {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    tuple val(meta), path(liftover_tsv)
    path bed

    output:
    tuple val(meta), path("${meta.id}.primers.bed")            , emit: bed
    tuple val(meta), path("${meta.id}.primers_liftover.tsv")   , emit: report
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # A primer scheme is designed once against the lab reference, but ivar trim
    # needs the positions of whatever the reads were aligned to. One indel
    # between the stock and the reference shifts every primer after it.
    python3 \$(which liftover_bed.py || echo ${projectDir}/bin/liftover_bed.py) \\
        --liftover-tsv ${liftover_tsv} \\
        --input-bed ${bed} \\
        --output-bed ${prefix}.primers.bed \\
        --report-tsv ${prefix}.primers_liftover.tsv \\
        --direction ref-to-stock

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.primers.bed
    touch ${prefix}.primers_liftover.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
