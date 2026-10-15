process KRAKEN2_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_kraken2}"

    input:
    tuple val(meta), path(reads)
    path db   // a Kraken2 database directory, or a .tar.gz of one

    output:
    tuple val(meta), path("*.kraken2.report.txt") , emit: report
    tuple val(meta), path("*.kraken2.db_info.txt"), emit: db_info
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: ''
    // Kraken2 loads hash.k2d into memory, so this process needs at least that much. The viral
    // database is about 0.5 GB, PlusPF-16 about 15 GB, full PlusPF about 104 GB. The check below
    // says so before Kraken2 starts, instead of letting it die partway.
    def mem_bytes = task.memory ? task.memory.toBytes() : 0
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

    need=\$(stat -c %s "\$db/hash.k2d")
    have=${mem_bytes}
    if [ "\$have" -gt 0 ] && [ "\$need" -gt "\$have" ]; then
        echo "ERROR: this Kraken2 database needs \$(( need / 1024 / 1024 / 1024 + 1 )) GB of memory and this step was given \$(( have / 1024 / 1024 / 1024 )) GB." >&2
        echo "       Kraken2 loads hash.k2d into memory in full. Raise the memory for this step, for example:" >&2
        echo "         process { withName: KRAKEN2_CLASSIFY { memory = \$(( need / 1024 / 1024 / 1024 + 2 )).GB } }" >&2
        echo "       conf/kraken2_large_db.example.config has a worked example." >&2
        exit 1
    fi

    # Every Kraken2 download ships inspect.txt, which lists what is in the database. Its domain
    # rows say what this run could possibly have named. Without them a large unclassified
    # fraction is unreadable: it may mean nothing is there, or that the database never held it.
    printf 'database\\t%s\\n' "\$(basename "${db}")" > ${prefix}.kraken2.db_info.txt
    if [ -s "\$db/inspect.txt" ]; then
        awk -F'\\t' '\$4 == "D" { gsub(/^ +/, "", \$6); printf "domain\\t%s\\n", \$6 }' "\$db/inspect.txt" \
            >> ${prefix}.kraken2.db_info.txt
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
    printf 'database\tstub\ndomain\tViruses\n' > ${prefix}.kraken2.db_info.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: stub
    END_VERSIONS
    """
}
