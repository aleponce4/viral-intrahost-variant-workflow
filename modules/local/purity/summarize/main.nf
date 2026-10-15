process PURITY_SUMMARIZE {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    tuple val(meta), path(stock_counts), path(unmapped_kinds)
    // Each is the NO_FILE placeholder when not given. They are staged in their own directories,
    // because the placeholder has the same name every time and two of them would collide.
    path pool_counts, stageAs: 'pool/*'
    path kraken_report, stageAs: 'kraken/*'
    path kraken_db_info, stageAs: 'kraken_db/*'

    output:
    tuple val(meta), path("${meta.id}.purity.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}.purity.md") , emit: markdown
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ''
    // A staged placeholder is named pool/NO_FILE, so the check looks at the last segment.
    def given  = { f -> f.toString().tokenize('/').last() != 'NO_FILE' }
    def pool   = given(pool_counts) ? "--pool-counts ${pool_counts}" : ''
    def kraken = given(kraken_report) ? "--kraken-report ${kraken_report}" : ''
    def dbinfo = given(kraken_db_info) ? "--kraken-db-info ${kraken_db_info}" : ''
    """
    python3 \$(which summarize_purity.py || echo ${projectDir}/bin/summarize_purity.py) \\
        --stock-name ${prefix} \\
        --stock-counts ${stock_counts} \\
        --unmapped-kinds ${unmapped_kinds} \\
        ${pool} \\
        ${kraken} \\
        ${dbinfo} \\
        --output-tsv ${prefix}.purity.tsv \\
        --output-md ${prefix}.purity.md \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.purity.tsv
    touch ${prefix}.purity.md

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
