process APPEND_REF_COORDS {
    tag "$meta.id"
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    tuple val(meta), path(variants)
    path liftover_tsv

    output:
    tuple val(meta), path("*.refcoords.*"), emit: annotated
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // sample.csq.vcf -> sample.csq.refcoords.vcf, sample.vcf.gz -> sample.refcoords.vcf
    def stem = variants.name.replaceFirst(/\.gz$/, '')
    def ext  = stem.tokenize('.').last()
    def base = stem.substring(0, stem.length() - ext.length() - 1)
    """
    python3 \$(which append_ref_coords.py || echo ${projectDir}/bin/append_ref_coords.py) \\
        --liftover-tsv ${liftover_tsv} \\
        --input ${variants} \\
        --output ${base}.refcoords.${ext} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def stem = variants.name.replaceFirst(/\.gz$/, '')
    def ext  = stem.tokenize('.').last()
    def base = stem.substring(0, stem.length() - ext.length() - 1)
    """
    touch ${base}.refcoords.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
