process PURITY_BAM_STATS {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_samtools}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.counts.tsv")        , emit: counts
    tuple val(meta), path("*.unmapped_kinds.tsv"), emit: kinds
    tuple val(meta), path("*.flagstat.txt")      , emit: flagstat
    tuple val(meta), path("*.unmapped.fastq.gz") , emit: unmapped
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Counts are of primary alignments. -F 2308 drops unmapped (4), secondary (256) and
    // supplementary (2048) records, so a read is counted once. The last row, `*`, holds the
    // reads that did not map.
    //
    // Most reads that do not map are not another virus. A read with a run of 20 Gs is a
    // poly-G read, which two-colour instruments write where they saw no signal. A read with
    // an Illumina adapter (TruSeq or Nextera) is an insert shorter than the read. Both are
    // counted in the kinds file and left out of the FASTQ that goes on to Kraken2, so they
    // are not mistaken for contamination. Left in, poly-G reads were classified as a phage in
    // a real run.
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${prefix}.flagstat.txt

    samtools view -F 2308 ${bam} \\
        | awk -F'\\t' '{n[\$3]++; if (\$5 == 0) z[\$3]++} END {for (c in n) printf "%s\\t%d\\t%d\\n", c, n[c], z[c] + 0}' \\
        > ${prefix}.mapped.tmp

    printf 'contig\\tlength\\tmapped\\tmapq0\\n' > ${prefix}.counts.tsv
    samtools idxstats ${bam} \\
        | awk -F'\\t' 'BEGIN {OFS = "\\t"} NR == FNR {n[\$1] = \$2; z[\$1] = \$3; next} \$1 != "*" {print \$1, \$2, n[\$1] + 0, z[\$1] + 0}' ${prefix}.mapped.tmp - \\
        >> ${prefix}.counts.tsv
    unmapped=\$(samtools view -c -f 4 -F 2304 ${bam})
    printf '*\\t0\\t%s\\t0\\n' "\$unmapped" >> ${prefix}.counts.tsv
    rm ${prefix}.mapped.tmp

    samtools fastq -@ ${task.cpus} -f 4 -F 2304 ${bam} 2> ${prefix}.fastq.log \\
        | awk 'BEGIN {polyg = "GGGGGGGGGGGGGGGGGGGG"; a1 = "AGATCGGAAGAGC"; a2 = "CTGTCTCTTATACACATCT"}
               {h = \$0; getline s; getline p; getline q; total++
                if (index(s, polyg) > 0) g++
                else if (index(s, a1) > 0 || index(s, a2) > 0) a++
                else print h "\\n" s "\\n" p "\\n" q}
               END {printf "kind\\treads\\ntotal\\t%d\\npolyg\\t%d\\nadapter\\t%d\\n", total + 0, g + 0, a + 0 > "${prefix}.unmapped_kinds.tsv"}' \\
        | gzip -1 > ${prefix}.unmapped.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'contig\\tlength\\tmapped\\tmapq0\\n*\\t0\\t0\\t0\\n' > ${prefix}.counts.tsv
    printf 'kind\\treads\\ntotal\\t0\\npolyg\\t0\\nadapter\\t0\\n' > ${prefix}.unmapped_kinds.tsv
    touch ${prefix}.flagstat.txt
    touch ${prefix}.unmapped.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
