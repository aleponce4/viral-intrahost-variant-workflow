process KRAKEN2_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_kraken2}"

    input:
    tuple val(meta), path(reads)
    path db   // a Kraken2 database directory, or a .tar.gz of one

    output:
    tuple val(meta), path("*.kraken2.report.txt"), emit: report
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ''
    // Kraken2 loads the whole database into memory, so this process needs at least as much
    // memory as the database is large. The viral database is about 0.5 GB.
    """
    db="${db}"
    case "\$db" in
        *.tar.gz|*.tgz)
            mkdir kraken2_db
            tar -xzf "\$db" -C kraken2_db
            hash=\$(find kraken2_db -name hash.k2d | head -n 1)
            if [ -z "\$hash" ]; then
                echo "ERROR: no hash.k2d in ${db}, so it is not a Kraken2 database" >&2
                exit 1
            fi
            db=\$(dirname "\$hash")
            ;;
    esac
    if [ ! -s "\$db/hash.k2d" ]; then
        echo "ERROR: no hash.k2d in ${db}, so it is not a Kraken2 database" >&2
        exit 1
    fi

    if [ -z "\$(gzip -dc ${reads} | head -c 1)" ]; then
        # Kraken2 stops on an empty file, and no unmapped reads is a valid answer.
        printf ' 100.00\\t0\\t0\\tU\\t0\\tunclassified\\n' > ${prefix}.kraken2.report.txt
    else
        kraken2 \\
            --db "\$db" \\
            --threads ${task.cpus} \\
            --gzip-compressed \\
            --report ${prefix}.kraken2.report.txt \\
            --output /dev/null \\
            ${args} \\
            ${reads}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(kraken2 --version 2>&1 | head -n 1 | sed 's/^.*version //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf ' 100.00\\t0\\t0\\tU\\t0\\tunclassified\\n' > ${prefix}.kraken2.report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: stub
    END_VERSIONS
    """
}
