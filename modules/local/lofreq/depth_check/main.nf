process LOFREQ_DEPTH_CHECK {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    tuple val(meta), path(vcf), path(depth)

    output:
    tuple val(meta), path("*.depth_check.tsv"), emit: tsv
    tuple val(meta), path("*.depth_check.txt"), emit: summary
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ''
    """
    python3 \$(which check_lofreq_depth.py || echo ${projectDir}/bin/check_lofreq_depth.py) \\
        --vcf ${vcf} \\
        --depth ${depth} \\
        --output ${prefix}.depth_check.tsv \\
        ${args} > ${prefix}.depth_check.txt
    cat ${prefix}.depth_check.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.depth_check.tsv ${prefix}.depth_check.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
