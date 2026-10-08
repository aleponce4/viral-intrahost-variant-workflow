process STOCK_REPORT {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    tuple val(meta), path(fasta), path(fai)
    tuple val(meta_check), path(annotation_check)
    tuple val(meta_sum), path(liftover_summary)
    tuple val(meta_or), path(orientation)
    tuple val(meta_un), path(unmapped)
    path primer_report

    output:
    tuple val(meta), path("${meta.id}.stock_report.md") , emit: markdown
    tuple val(meta), path("${meta.id}.stock_report.tsv"), emit: tsv
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def primers = primer_report.name != 'NO_FILE' ? "--primer-report ${primer_report}" : ''
    """
    python3 \$(which summarize_stock_reference.py || echo ${projectDir}/bin/summarize_stock_reference.py) \\
        --stock-name ${prefix} \\
        --fasta ${fasta} \\
        --annotation-check ${annotation_check} \\
        --liftover-summary ${liftover_summary} \\
        --orientation ${orientation} \\
        --unmapped ${unmapped} \\
        ${primers} \\
        --output-md ${prefix}.stock_report.md \\
        --output-tsv ${prefix}.stock_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stock_report.md
    touch ${prefix}.stock_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
