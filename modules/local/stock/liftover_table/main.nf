process STOCK_LIFTOVER_TABLE {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    tuple val(meta), path(paf)

    output:
    tuple val(meta), path("${meta.id}.liftover.tsv"), emit: table
    tuple val(meta), path("${meta.id}.liftover_summary.tsv"), emit: summary
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 \$(which build_liftover_table.py || echo ${projectDir}/bin/build_liftover_table.py) \\
        --paf ${paf} \\
        --output-tsv ${prefix}.liftover.tsv \\
        --summary-tsv ${prefix}.liftover_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.liftover.tsv
    touch ${prefix}.liftover_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
