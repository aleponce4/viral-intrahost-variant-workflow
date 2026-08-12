process REPORT_COVERAGE_PLOTS {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path coverage_dir

    output:
    path "*.png"        , emit: plots, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_coverage.py --coverage-dir ${coverage_dir} --output-dir . --dataset ${params.dataset}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${params.dataset}_coverage_summary.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
