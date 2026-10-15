process PRIMER_CALL_CHECK {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    tuple val(meta), path(vcf), path(depth)
    path primer_bed

    output:
    tuple val(meta), path("*.primer_check.tsv"), emit: tsv
    tuple val(meta), path("*.primer_check.txt"), emit: summary
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ''
    """
    python3 \$(which check_primer_calls.py || echo ${projectDir}/bin/check_primer_calls.py) \\
        --vcf ${vcf} \\
        --primer-bed ${primer_bed} \\
        --depth ${depth} \\
        --output ${prefix}.primer_check.tsv \\
        ${args} > ${prefix}.primer_check.txt
    cat ${prefix}.primer_check.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.primer_check.tsv ${prefix}.primer_check.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
