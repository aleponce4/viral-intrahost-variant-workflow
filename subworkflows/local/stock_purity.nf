/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: STOCK_PURITY

    A stock is meant to hold one virus. This takes a sample of its reads and asks
    three questions of them:

      1. How many map to the stock's own reference? A pure stock puts nearly all there.
      2. Where do they go when they can choose between the lab's other genomes (the
         reference pool)? A second genome with a real share means a second strain.
         Kraken2 cannot tell close strains apart, so this uses mapping.
      3. What are the reads that map to nothing? Most are sequencing artefacts (poly-G and
         adapter reads), which are counted and set aside. Kraken2 on the rest finds other
         viruses, and bacteria or host when the database holds them. Only when a database
         is given.

    It reports. A flag says where to look.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PURITY_SUBSAMPLE                        } from '../../modules/local/purity/subsample/main'
include { BWA_INDEX as PURITY_STOCK_INDEX         } from '../../modules/local/bwa/index/main'
include { BWA_INDEX as PURITY_POOL_INDEX          } from '../../modules/local/bwa/index/main'
include { BWA_MEM as PURITY_STOCK_MEM             } from '../../modules/local/bwa/mem/main'
include { BWA_MEM as PURITY_POOL_MEM              } from '../../modules/local/bwa/mem/main'
include { SAMTOOLS_SORT as PURITY_STOCK_SORT      } from '../../modules/local/samtools/sort/main'
include { SAMTOOLS_SORT as PURITY_POOL_SORT       } from '../../modules/local/samtools/sort/main'
include { PURITY_BAM_STATS as PURITY_STOCK_STATS  } from '../../modules/local/purity/bam_stats/main'
include { PURITY_BAM_STATS as PURITY_POOL_STATS   } from '../../modules/local/purity/bam_stats/main'
include { KRAKEN2_CLASSIFY                        } from '../../modules/local/kraken2/classify/main'
include { PURITY_SUMMARIZE                        } from '../../modules/local/purity/summarize/main'

workflow STOCK_PURITY {
    take:
    ch_stock_fasta // channel: [ val(meta), path(fasta), path(fai) ]
    ch_reads       // channel: [ val(meta), path(fastq_1), path(fastq_2) ]
    ch_pool        // value:   [ val(meta), path(fasta) ] when --purity_pool is set, else not read
    ch_kraken_db   // value:   path(database) when --purity_kraken2_db is set, else not read

    main:
    ch_versions = Channel.empty()

    // A sample of the library. The check does not need all of it.
    PURITY_SUBSAMPLE(ch_reads)
    ch_sample = PURITY_SUBSAMPLE.out.reads
    ch_versions = ch_versions.mix(PURITY_SUBSAMPLE.out.versions)

    // 1. The reads against the stock's own reference. The reads that do not map go on to
    //    Kraken2. The names carry the step, because both mappings write a counts file
    //    into the same directory.
    PURITY_STOCK_INDEX(ch_stock_fasta.map { meta, fasta, fai -> [ [ id: "${meta.id}.vs_stock" ], fasta ] })
    PURITY_STOCK_MEM(
        ch_sample.map { meta, r1, r2 -> [ [ id: "${meta.id}.vs_stock" ], r1, r2 ] },
        PURITY_STOCK_INDEX.out.index.first()
    )
    PURITY_STOCK_SORT(PURITY_STOCK_MEM.out.sam)
    PURITY_STOCK_STATS(PURITY_STOCK_SORT.out.bam)
    ch_versions = ch_versions
        .mix(PURITY_STOCK_INDEX.out.versions)
        .mix(PURITY_STOCK_MEM.out.versions)
        .mix(PURITY_STOCK_SORT.out.versions)
        .mix(PURITY_STOCK_STATS.out.versions)

    // 2. The same reads against the reference pool.
    ch_pool_counts = Channel.value(file("${projectDir}/assets/NO_FILE"))
    if (params.purity_pool) {
        PURITY_POOL_INDEX(ch_pool.map { meta, fasta -> [ [ id: "${meta.id}.vs_pool" ], fasta ] })
        PURITY_POOL_MEM(
            ch_sample.map { meta, r1, r2 -> [ [ id: "${meta.id}.vs_pool" ], r1, r2 ] },
            PURITY_POOL_INDEX.out.index.first()
        )
        PURITY_POOL_SORT(PURITY_POOL_MEM.out.sam)
        PURITY_POOL_STATS(PURITY_POOL_SORT.out.bam)
        ch_pool_counts = PURITY_POOL_STATS.out.counts.map { meta, counts -> counts }
        ch_versions = ch_versions
            .mix(PURITY_POOL_INDEX.out.versions)
            .mix(PURITY_POOL_MEM.out.versions)
            .mix(PURITY_POOL_SORT.out.versions)
            .mix(PURITY_POOL_STATS.out.versions)
    }

    // 3. What the reads that did not map are. Poly-G and adapter reads were set aside in the
    //    stats step, so Kraken2 sees only reads that nothing explains yet.
    ch_kraken_report = Channel.value(file("${projectDir}/assets/NO_FILE"))
    if (params.purity_kraken2_db) {
        KRAKEN2_CLASSIFY(PURITY_STOCK_STATS.out.unmapped, ch_kraken_db)
        ch_kraken_report = KRAKEN2_CLASSIFY.out.report.map { meta, report -> report }
        ch_versions = ch_versions.mix(KRAKEN2_CLASSIFY.out.versions)
    }

    // The summary is named after the stock, not after the step.
    PURITY_SUMMARIZE(
        PURITY_STOCK_STATS.out.counts
            .join(PURITY_STOCK_STATS.out.kinds)
            .map { meta, counts, kinds -> [ [ id: meta.id.replaceFirst(/\.vs_stock$/, '') ], counts, kinds ] },
        ch_pool_counts,
        ch_kraken_report
    )
    ch_versions = ch_versions.mix(PURITY_SUMMARIZE.out.versions)

    emit:
    tsv      = PURITY_SUMMARIZE.out.tsv       // [ meta, tsv ]
    markdown = PURITY_SUMMARIZE.out.markdown  // [ meta, md ]
    versions = ch_versions
}
