process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    container "${params.container_minimap2}"

    input:
    tuple val(meta), path(query), path(query_fai)
    tuple val(meta_ref), path(target), path(target_fai)

    output:
    tuple val(meta), path("*.paf"), emit: paf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -o pipefail

    # asm10 is the preset for assemblies up to a few percent divergent, which is
    # where two strains of the same alphavirus sit. --cs records the per-base
    # difference string that build_liftover_table.py walks, so it is required.
    minimap2 -cx asm10 --cs -t ${task.cpus} ${target} ${query} > ${prefix}.paf

    if [ ! -s ${prefix}.paf ]; then
        echo "ERROR: ${query} did not align to ${target}. Check that the stock and the reference are the same virus." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
