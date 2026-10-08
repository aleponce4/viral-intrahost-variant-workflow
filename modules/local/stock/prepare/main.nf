process STOCK_PREPARE {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_samtools}"

    input:
    tuple val(meta), path(consensus)
    tuple val(meta_paf), path(paf)

    output:
    tuple val(meta), path("${meta.id}.fasta"), path("${meta.id}.fasta.fai"), emit: fasta
    tuple val(meta), path("${meta.id}.orientation.txt")                    , emit: orientation
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

    awk -v name="${prefix}" 'NR==1 { print ">" name; next } { print }' ${consensus} > renamed.fasta
    samtools faidx renamed.fasta

    # An assembler has no way to know which strand a genome is meant to be read
    # on, and picks one per run: on the same reads, two runs of the same assembly
    # pipeline gave opposite orientations. Left alone, half the stocks would get a
    # reference whose genes sit on the minus strand and whose coordinates run
    # backwards against every other stock. Take the orientation from the best
    # alignment to the lab reference.
    STRAND=\$(sort -k10,10nr ${paf} | head -n 1 | cut -f5)
    if [ -z "\$STRAND" ]; then
        echo "ERROR: ${paf} is empty; cannot tell which strand ${consensus} is on" >&2
        exit 1
    fi

    if [ "\$STRAND" = "-" ]; then
        samtools faidx -i renamed.fasta ${prefix} \\
            | awk -v name="${prefix}" 'NR==1 { print ">" name; next } { print }' > ${prefix}.fasta
        echo "reverse-complemented" > ${prefix}.orientation.txt
    else
        mv renamed.fasta ${prefix}.fasta
        echo "as assembled" > ${prefix}.orientation.txt
    fi

    rm -f renamed.fasta.fai
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
    touch ${prefix}.orientation.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
