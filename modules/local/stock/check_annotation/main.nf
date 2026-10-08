process STOCK_CHECK_ANNOTATION {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    tuple val(meta), path(stock_fasta), path(stock_fai)
    tuple val(meta_gff), path(stock_gff)
    tuple val(meta_ref), path(ref_fasta), path(ref_fai)
    tuple val(meta_refgff), path(ref_gff)

    output:
    tuple val(meta), path("${meta.id}.annotation_check.tsv"), emit: report
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # The reference is passed as a control: a CDS problem that the reference
    # already has (the nsP3 opal readthrough codon, for one) is inherited, not
    # caused by the transfer, and must not stop the run.
    python3 \$(which check_lifted_annotation.py || echo ${projectDir}/bin/check_lifted_annotation.py) \\
        --fasta ${stock_fasta} \\
        --gff ${stock_gff} \\
        --reference-fasta ${ref_fasta} \\
        --reference-gff ${ref_gff} \\
        --output-tsv ${prefix}.annotation_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.annotation_check.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
