process STOCK_PREPARE {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_samtools}"

    input:
    tuple val(meta), path(consensus)

    output:
    tuple val(meta), path("${meta.id}.fasta"), path("${meta.id}.fasta.fai"), emit: fasta
    path "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail

    # An assembler names its contigs after itself, so every stock would arrive
    # carrying a header like NODE_1_length_11446_cov_902. Rename it to the stock,
    # because that name becomes the CHROM of every variant called against it.
    if [ "\$(grep -c '^>' ${consensus})" -ne 1 ]; then
        echo "ERROR: ${consensus} holds \$(grep -c '^>' ${consensus}) sequences; a stock reference must be a single contig" >&2
        exit 1
    fi

    awk -v name="${prefix}" 'NR==1 { print ">" name; next } { print }' ${consensus} > ${prefix}.fasta
    samtools faidx ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    touch ${prefix}.fasta.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
