## Purity

The first 500,000 read pairs (1,000,000 reads) of the stock's sequencing run. The counts are of reads, so a pair counts twice.

### Reads against the stock reference

| | Reads | Share |
| :--- | ---: | ---: |
| Map to the stock reference | 950,000 | 95.00% |
| Do not map | 50,000 | 5.00% |
| of which poly-G reads | 6,000 | 0.60% |
| of which carry an Illumina adapter | 4,000 | 0.40% |
| of which neither | 40,000 | 4.00% |

Poly-G reads are where the instrument saw no signal and read G. Adapter reads are inserts shorter than the read. Neither is another virus, so they are counted and set aside before the reads are classified.

### Reads against the reference pool

Each read goes to the genome it fits best. A stock of one virus puts nearly all of its reads on the genome closest to it. The reads with mapping quality 0 fit two genomes equally well, so shares are of the other reads.

981,500 reads map to the pool in all.

| Genome | Reads | Mapping quality 0 | Share of the rest |
| :--- | ---: | ---: | ---: |
| `strainA` | 900,000 | 10,000 | 92.61% |
| `strainB` | 80,000 | 10,000 | 7.28% |
| `strainC` | 1,500 | 500 | 0.10% |

### What the unmapped reads are

Kraken2 classified 30,000 of the 50,000 reads left after the poly-G and adapter reads. 40.00% are unclassified.

| Domain | Reads |
| :--- | ---: |
| Viruses | 25,000 |
| Bacteria | 5,000 |

| Top species | Reads | Share of all reads |
| :--- | ---: | ---: |
| Venezuelan equine encephalitis virus | 20,000 | 2.00% |
| Escherichia coli | 5,000 | 0.50% |
| Homo sapiens | 100 | 0.01% |

### Flags

- 70,000 reads (7.3% of the reads placed uniquely on the pool) fit strainB best, not strainA, which is the genome closest to the stock
- 20,000 reads (2.0% of all) classify as Venezuelan equine encephalitis virus
