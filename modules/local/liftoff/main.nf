process LIFTOFF {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_liftoff}"

    input:
    tuple val(meta), path(stock_fasta), path(stock_fai)
    tuple val(meta_ref), path(ref_fasta), path(ref_fai)
    tuple val(meta_gff), path(ref_gff)

    output:
    tuple val(meta), path("${meta.id}.gff3")   , emit: gff
    tuple val(meta), path("${meta.id}.unmapped.txt"), emit: unmapped
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail

    # Liftoff writes its gffutils database next to the input GFF, so give it a
    # writable copy rather than the staged read-only one.
    cp ${ref_gff} reference_input.gff3

    liftoff \\
        -g reference_input.gff3 \\
        -o ${prefix}.gff3 \\
        -u ${prefix}.unmapped.txt \\
        -dir intermediate_files \\
        -p ${task.cpus} \\
        ${stock_fasta} \\
        ${ref_fasta}

    # Liftoff writes the unmapped file only when it has something to say.
    touch ${prefix}.unmapped.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        liftoff: \$(liftoff --version 2>&1 | tail -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff3
    touch ${prefix}.unmapped.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        liftoff: \$(liftoff --version 2>&1 | tail -n1)
    END_VERSIONS
    """
}
