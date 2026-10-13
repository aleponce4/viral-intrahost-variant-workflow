/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: INPUT_CHECK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList } from 'plugin/nf-schema'
include { SAMTOOLS_FAIDX    } from '../../modules/local/samtools/faidx/main'

workflow INPUT_CHECK {
    take:
    samplesheet // path: samplesheet.csv

    main:
    if (!params.fasta) {
        error "Parameter 'params.fasta' must be specified."
    }
    if (!params.gff) {
        error "Parameter 'params.gff' must be specified."
    }
    if (params.protocol == 'amplicon' && !params.primer_bed) {
        error "Parameter --primer_bed is required when --protocol amplicon."
    }
    if (params.protocol == 'metagenomic' && params.primer_bed) {
        log.warn "--primer_bed supplied but --protocol is 'metagenomic'; primer trimming is OFF."
    }

    // Read the annotation once and check the three things that otherwise surface late or not
    // at all: no CDS, a CDS that bcftools csq cannot reach, and names that do not match the
    // FASTA. Each of these used to cost a full alignment and variant-calling run first.
    def gff_file = file(params.gff, checkIfExists: true)
    def feature_type = [:]   // feature ID -> its type, for resolving a CDS Parent
    def cds_parents  = []    // the Parent of every CDS, in file order
    def cds_contigs  = [] as Set
    gff_file.eachLine { line ->
        if (line.startsWith('#')) return
        def field = line.split('\t')
        if (field.size() < 9) return
        def attrs = [:]
        field[8].split(';').each { pair ->
            def sep = pair.indexOf('=')
            if (sep > 0) attrs[pair.substring(0, sep).trim()] = pair.substring(sep + 1).trim()
        }
        if (attrs.ID) feature_type[attrs.ID] = field[2]
        if (field[2] == 'CDS') {
            cds_contigs << field[0]
            cds_parents << (attrs.Parent ? attrs.Parent.split(',')[0].trim() : null)
        }
    }

    if (cds_parents.isEmpty()) {
        error "Provided GFF file ${params.gff} does not contain CDS features required for annotation and selection analysis."
    }

    // bcftools csq walks gene -> transcript -> CDS. A GFF that attaches a CDS straight to its
    // gene indexes zero transcripts, and the run dies at BCFTOOLS_CSQ with "No usable
    // transcripts found" after everything upstream has already run.
    def reachable = cds_parents.count { parent ->
        parent != null && feature_type[parent] != null && feature_type[parent] != 'gene'
    }
    if (reachable == 0) {
        def example = cds_parents.find { it != null }
        def why = example == null ? "its CDS features have no Parent attribute"
                : feature_type[example] == 'gene' ? "its CDS features name a gene (ID=${example}) as their Parent"
                : "the Parent of its CDS features (${example}) is not defined in the file"
        error """Provided GFF file ${params.gff} has no CDS that bcftools csq can use, because ${why}.
  bcftools csq reads gene -> transcript -> CDS, so a CDS needs a transcript (mRNA) between it and
  its gene. Without one it indexes zero transcripts and annotation fails after alignment and
  variant calling have already run.
  Add a transcript line and point the CDS at it, for example:
    <contig>  <source>  mRNA  <start>  <end>  .  <strand>  .  ID=rna-X;Parent=gene-X
    <contig>  <source>  CDS   <start>  <end>  .  <strand>  0  ID=cds-X;Parent=rna-X
  Or rebuild the annotation from its GenBank record with bin/parse_gb_to_gff3.py, which writes
  the transcript line for you."""
    }
    if (reachable < cds_parents.size()) {
        log.warn "${cds_parents.size() - reachable} of ${cds_parents.size()} CDS features in ${params.gff} are not reachable through a transcript; bcftools csq will skip those."
    }

    // bcftools csq annotates nothing, and says nothing, when the GFF names a contig the FASTA
    // does not have. One annotation file covering several segments is fine, so this only fails
    // when no CDS at all belongs to the reference being used.
    def fasta_contigs = [] as Set
    file(params.fasta, checkIfExists: true).eachLine { line ->
        if (line.startsWith('>')) fasta_contigs << line.substring(1).trim().split(/\s/)[0]
    }
    def shared = cds_contigs.intersect(fasta_contigs)
    if (shared.isEmpty()) {
        error """The GFF and the FASTA do not name the same contig, so bcftools csq would annotate
  nothing and report no error.
    ${params.gff} has CDS on: ${cds_contigs.sort().join(', ')}
    ${params.fasta} holds:    ${fasta_contigs.sort().join(', ')}
  Use an annotation built against this reference, or rename the contig in one of them."""
    }
    if (cds_contigs.size() > shared.size()) {
        log.warn "${params.gff} also annotates ${(cds_contigs - fasta_contigs).sort().join(', ')}, which ${params.fasta} does not hold. Those features are ignored."
    }

    // --viral_contig picks the contig the BAM is cut down to and the one ivar reads. A name that
    // is not in the reference yields an empty BAM and no calls, with nothing said about why.
    if (params.viral_contig && !fasta_contigs.contains(params.viral_contig as String)) {
        error "--viral_contig '${params.viral_contig}' is not in ${params.fasta}, which holds: ${fasta_contigs.sort().join(', ')}"
    }

    def root_dir = workflow.projectDir.toString().split('/.nf-test')[0].split('\\.nf-test')[0]

    ch_samplesheet = Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map { row ->
            def meta = [:]
            meta.id        = row[0]
            meta.treatment = row[3]
            meta.condition = row[3]
            def f1 = row[1].toString()
            def f2 = row[2].toString()
            def file1 = f1.startsWith('/') ? file(f1, checkIfExists: true) : (file(f1).exists() ? file(f1) : file("${root_dir}/${f1}", checkIfExists: true))
            def file2 = f2.startsWith('/') ? file(f2, checkIfExists: true) : (file(f2).exists() ? file(f2) : file("${root_dir}/${f2}", checkIfExists: true))
            return [ meta, file1, file2 ]
        }

    ch_raw_fasta = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> [ [id: 'reference'], fasta ] }

    SAMTOOLS_FAIDX(ch_raw_fasta)
    ch_indexed_fasta = SAMTOOLS_FAIDX.out.fai

    ch_gff = Channel
        .fromPath(params.gff, checkIfExists: true)
        .map { gff -> [ [id: 'annotation'], gff ] }

    emit:
    fastqs   = ch_samplesheet // channel: [ val(meta), path(fastq_1), path(fastq_2) ]
    fasta    = ch_indexed_fasta// channel: [ val(meta), path(fasta), path(fai) ]
    gff      = ch_gff         // channel: [ val(meta), path(gff) ]
    versions = SAMTOOLS_FAIDX.out.versions
}
