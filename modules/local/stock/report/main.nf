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
    // Each optional input is staged in its own directory, because the placeholder for
    // "not given" has the same name every time and two of them would collide.
    path primer_report, stageAs: 'primer/*'
    path purity_tsv, stageAs: 'purity_tsv/*'
    path purity_md, stageAs: 'purity_md/*'

    output:
    tuple val(meta), path("${meta.id}.stock_report.md") , emit: markdown
    tuple val(meta), path("${meta.id}.stock_report.tsv"), emit: tsv
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // A staged placeholder is named primer/NO_FILE, so the check looks at the last segment.
    def given   = { f -> f.toString().tokenize('/').last() != 'NO_FILE' }
    def primers = given(primer_report) ? "--primer-report ${primer_report}" : ''
    def purity  = given(purity_tsv) ? "--purity-tsv ${purity_tsv} --purity-md ${purity_md}" : ''
    """
    python3 \$(which summarize_stock_reference.py || echo ${projectDir}/bin/summarize_stock_reference.py) \\
        --stock-name ${prefix} \\
        --fasta ${fasta} \\
        --annotation-check ${annotation_check} \\
        --liftover-summary ${liftover_summary} \\
        --orientation ${orientation} \\
        --unmapped ${unmapped} \\
        ${primers} \\
        ${purity} \\
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
