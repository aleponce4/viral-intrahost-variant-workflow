process PURITY_SUBSAMPLE {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_bwa}"

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)

    output:
    tuple val(meta), path("*.sub_1.fastq.gz"), path("*.sub_2.fastq.gz"), emit: reads
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The first N read pairs. A purity check needs a sample of the library, not all of it.
    // At a few hundred thousand pairs a contaminant at 0.1% still has hundreds of reads.
    def lines = (params.purity_pairs as long) * 4
    """
    gzip -dcf ${fastq_1} | head -n ${lines} | gzip -1 > ${prefix}.sub_1.fastq.gz
    gzip -dcf ${fastq_2} | head -n ${lines} | gzip -1 > ${prefix}.sub_2.fastq.gz

    n1=\$(gzip -dc ${prefix}.sub_1.fastq.gz | wc -l)
    n2=\$(gzip -dc ${prefix}.sub_2.fastq.gz | wc -l)
    if [ "\$n1" -eq 0 ]; then
        echo "ERROR: no reads could be read from ${fastq_1}" >&2
        exit 1
    fi
    if [ "\$n1" -ne "\$n2" ]; then
        echo "ERROR: the two read files hold different numbers of reads (\$((n1 / 4)) and \$((n2 / 4)))" >&2
        exit 1
    fi
    echo "Kept \$((n1 / 4)) read pairs" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version 2>&1 | head -n 1 | sed 's/^gzip //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sub_1.fastq.gz
    touch ${prefix}.sub_2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: stub
    END_VERSIONS
    """
}
