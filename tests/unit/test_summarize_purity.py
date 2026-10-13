"""Unit tests for bin/summarize_purity.py.

Run from the repository root:  python -m unittest discover -s tests/unit -v

The script says whether a stock's reads look like one virus. A wrong flag either hides a
second strain in the tube or sends someone looking for one that is not there, so each
test states the expected number or flag explicitly.

Stock counts used below (reads, so a pair counts twice):

  mapped to the stock   950000
  not mapped             50000     -> 95.00% mapped, not flagged at the 90% default

Pool counts (reads, and how many of them have mapping quality 0):

  genome     mapped    mapq0    placed uniquely    share
  strainA    900000    10000        890000         0.9261
  strainB     80000    10000         70000         0.0728  (flagged: over 1% and over 100 reads)
  strainC      1500      500          1000         0.0010  (not flagged: under 1%)
"""
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SCRIPT = os.path.join(REPO, 'bin', 'summarize_purity.py')
sys.path.insert(0, os.path.join(REPO, 'bin'))

import summarize_purity as sp  # noqa: E402

HEADER = "contig\tlength\tmapped\tmapq0\n"

STOCK = HEADER + "stock\t11446\t950000\t0\n*\t0\t50000\t0\n"

POOL = HEADER + (
    "strainA\t11446\t900000\t10000\n"
    "strainB\t11417\t80000\t10000\n"
    "strainC\t11400\t1500\t500\n"
)

# Of the 50000 reads that did not map: 6000 poly-G, 4000 adapter, 40000 neither.
KINDS = "kind\treads\ntotal\t50000\npolyg\t6000\nadapter\t4000\n"

# A Kraken2 report. Columns: percent, reads in clade, reads at taxon, rank, taxid, name.
KRAKEN = (
    " 40.00\t20000\t20000\tU\t0\tunclassified\n"
    " 60.00\t30000\t0\tR\t1\troot\n"
    " 50.00\t25000\t0\tD\t10239\t  Viruses\n"
    " 10.00\t5000\t0\tD\t2\t  Bacteria\n"
    " 40.00\t20000\t20000\tS\t11036\t      Venezuelan equine encephalitis virus\n"
    " 10.00\t5000\t5000\tS\t562\t        Escherichia coli\n"
    "  0.20\t100\t100\tS\t9606\t        Homo sapiens\n"
)


class Workdir(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = self._tmp.name

    def path(self, name):
        return os.path.join(self.dir, name)

    def write(self, name, text):
        with open(self.path(name), 'w', encoding='utf-8', newline='') as f:
            f.write(text)
        return self.path(name)

    def run_script(self, *args):
        return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)

    def summarise(self, stock=STOCK, pool=None, kraken=None, kinds=None, extra=()):
        args = ['--stock-name', 'S', '--stock-counts', self.write('stock.tsv', stock),
                '--output-tsv', self.path('out.tsv'), '--output-md', self.path('out.md')]
        if kinds is not None:
            args += ['--unmapped-kinds', self.write('kinds.tsv', kinds)]
        if pool is not None:
            args += ['--pool-counts', self.write('pool.tsv', pool)]
        if kraken is not None:
            args += ['--kraken-report', self.write('k.txt', kraken)]
        result = self.run_script(*args, *extra)
        row, md = {}, ''
        if os.path.exists(self.path('out.tsv')):
            with open(self.path('out.tsv'), encoding='utf-8') as f:
                lines = f.read().splitlines()
            row = dict(zip(lines[0].split('\t'), lines[1].split('\t')))
        if os.path.exists(self.path('out.md')):
            with open(self.path('out.md'), encoding='utf-8') as f:
                md = f.read()
        return result, row, md


class TestStockMapping(Workdir):
    def test_output_has_the_documented_columns(self):
        result, row, _ = self.summarise()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(list(row), sp.COLUMNS)

    def test_fraction_mapped_is_mapped_over_all_reads(self):
        _, row, _ = self.summarise()
        self.assertEqual(row['reads_sampled'], '1000000')
        self.assertEqual(row['reads_mapped_to_stock'], '950000')
        self.assertEqual(row['reads_unmapped_to_stock'], '50000')
        self.assertEqual(row['fraction_mapped_to_stock'], '0.9500')

    def test_a_pure_stock_has_no_flag(self):
        _, row, md = self.summarise()
        self.assertEqual(row['flags'], 'none')
        self.assertIn('Nothing flagged.', md)

    def test_low_mapped_fraction_is_flagged(self):
        stock = HEADER + "stock\t11446\t500000\t0\n*\t0\t500000\t0\n"
        _, row, md = self.summarise(stock=stock)
        self.assertEqual(row['flags'], 'only 50.0% of the reads map to the stock reference')
        self.assertIn('only 50.0% of the reads map to the stock reference', md)

    def test_min_mapped_is_a_strict_lower_bound(self):
        # 95% mapped is not under 0.95, but is under 0.96.
        _, row, _ = self.summarise(extra=('--min-mapped', '0.95'))
        self.assertEqual(row['flags'], 'none')
        _, row, _ = self.summarise(extra=('--min-mapped', '0.96'))
        self.assertIn('only 95.0%', row['flags'])

    def test_a_stock_with_several_contigs_sums_them(self):
        stock = HEADER + "a\t100\t300000\t0\nb\t100\t650000\t0\n*\t0\t50000\t0\n"
        _, row, _ = self.summarise(stock=stock)
        self.assertEqual(row['reads_mapped_to_stock'], '950000')

    def test_markdown_states_the_read_counts(self):
        _, _, md = self.summarise()
        self.assertIn('The first 500,000 read pairs (1,000,000 reads)', md)
        self.assertIn('| Map to the stock reference | 950,000 | 95.00% |', md)
        self.assertIn('| Do not map | 50,000 | 5.00% |', md)


class TestKinds(Workdir):
    def test_poly_g_and_adapter_reads_are_counted(self):
        _, row, _ = self.summarise(kinds=KINDS)
        self.assertEqual(row['unmapped_polyg'], '6000')
        self.assertEqual(row['unmapped_adapter'], '4000')
        # 40000 of the 1,000,000 reads are neither mapped nor an artefact.
        self.assertEqual(row['fraction_unexplained'], '0.0400')

    def test_markdown_lists_them_and_says_why_they_are_set_aside(self):
        _, _, md = self.summarise(kinds=KINDS)
        self.assertIn('| of which poly-G reads | 6,000 | 0.60% |', md)
        self.assertIn('| of which carry an Illumina adapter | 4,000 | 0.40% |', md)
        self.assertIn('| of which neither | 40,000 | 4.00% |', md)
        self.assertIn('Neither is another virus', md)

    def test_without_the_kinds_file_the_columns_are_NA_and_everything_unmapped_is_unexplained(self):
        _, row, md = self.summarise()
        self.assertEqual(row['unmapped_polyg'], 'NA')
        self.assertEqual(row['unmapped_adapter'], 'NA')
        self.assertEqual(row['fraction_unexplained'], '0.0500')
        self.assertNotIn('poly-G', md)

    def test_artefact_reads_count_as_explained_for_the_low_mapped_flag(self):
        # Half the reads map, and 45% more are artefacts: 95% is explained, so no flag.
        stock = HEADER + "stock\t11446\t500000\t0\n*\t0\t500000\t0\n"
        kinds = "kind\treads\ntotal\t500000\npolyg\t300000\nadapter\t150000\n"
        _, row, _ = self.summarise(stock=stock, kinds=kinds)
        self.assertEqual(row['flags'], 'none')
        # Without the kinds the same counts are flagged.
        _, row, _ = self.summarise(stock=stock)
        self.assertIn('only 50.0%', row['flags'])

    def test_the_low_mapped_flag_names_the_artefacts_when_a_lot_is_still_unexplained(self):
        stock = HEADER + "stock\t11446\t500000\t0\n*\t0\t500000\t0\n"
        kinds = "kind\treads\ntotal\t500000\npolyg\t100000\nadapter\t50000\n"
        _, row, _ = self.summarise(stock=stock, kinds=kinds)
        self.assertEqual(row['flags'],
                         'only 50.0% of the reads map to the stock reference, '
                         'and 15.0% more are poly-G or adapter reads')

    def test_kraken_text_says_the_reads_were_left_after_the_artefacts(self):
        _, _, md = self.summarise(kinds=KINDS, kraken=KRAKEN)
        self.assertIn('of the 50,000 reads left after the poly-G and adapter reads', md)
        _, _, md = self.summarise(kraken=KRAKEN)
        self.assertIn('of the 50,000 reads that did not map to the stock', md)

    def test_a_poly_g_read_is_not_also_counted_as_adapter(self):
        # The file gives each read one kind, so the two counts never overlap.
        kinds = "kind\treads\ntotal\t50000\npolyg\t50000\nadapter\t0\n"
        _, row, _ = self.summarise(kinds=kinds)
        self.assertEqual(row['fraction_unexplained'], '0.0000')


class TestPool(Workdir):
    def test_top_genome_and_its_share(self):
        _, row, _ = self.summarise(pool=POOL)
        self.assertEqual(row['pool_top_genome'], 'strainA')
        # 890000 of 961000 uniquely placed reads.
        self.assertEqual(row['pool_top_share'], '0.9261')

    def test_the_best_other_genome_is_reported_and_flagged(self):
        _, row, md = self.summarise(pool=POOL)
        self.assertEqual(row['pool_other_max_genome'], 'strainB')
        self.assertEqual(row['pool_other_max_share'], '0.0728')
        self.assertIn('70,000 reads (7.3% of the reads placed uniquely on the pool) fit strainB best',
                      row['flags'])
        self.assertIn('not strainA, which is the genome closest to the stock', row['flags'])
        self.assertIn('fit strainB best', md)

    def test_a_genome_under_the_flag_fraction_is_not_flagged(self):
        _, row, _ = self.summarise(pool=POOL)
        self.assertNotIn('strainC', row['flags'])

    def test_flag_fraction_can_be_raised(self):
        _, row, _ = self.summarise(pool=POOL, extra=('--flag-fraction', '0.1'))
        self.assertEqual(row['flags'], 'none')

    def test_min_reads_keeps_a_handful_of_reads_from_raising_a_flag(self):
        pool = HEADER + "strainA\t11446\t100000\t0\nstrainB\t11417\t50\t0\n"
        # 50 of 100050 is 0.05%, but also under 100 reads. Lower the fraction to isolate min-reads.
        _, row, _ = self.summarise(pool=pool, extra=('--flag-fraction', '0.0001'))
        self.assertEqual(row['flags'], 'none')
        _, row, _ = self.summarise(pool=pool, extra=('--flag-fraction', '0.0001', '--min-reads', '10'))
        self.assertIn('fit strainB best', row['flags'])

    def test_mapping_quality_zero_reads_are_left_out_of_the_shares(self):
        # Every read on the second genome also fits the first equally well: not a second strain.
        pool = HEADER + "strainA\t11446\t900000\t0\nstrainB\t11417\t80000\t80000\n"
        _, row, _ = self.summarise(pool=pool)
        self.assertEqual(row['pool_other_max_share'], '0.0000')
        self.assertEqual(row['flags'], 'none')

    def test_a_single_genome_pool_has_no_other(self):
        pool = HEADER + "only\t11446\t900000\t0\n"
        _, row, _ = self.summarise(pool=pool)
        self.assertEqual(row['pool_top_genome'], 'only')
        self.assertEqual(row['pool_top_share'], '1.0000')
        self.assertEqual(row['pool_other_max_genome'], 'NA')
        self.assertEqual(row['flags'], 'none')

    def test_a_pool_with_no_unique_reads_is_flagged(self):
        pool = HEADER + "a\t100\t500\t500\nb\t100\t500\t500\n"
        _, row, _ = self.summarise(pool=pool)
        self.assertEqual(row['pool_top_genome'], 'NA')
        self.assertEqual(row['flags'], 'no read maps uniquely to the reference pool')

    def test_without_a_pool_the_pool_columns_are_NA(self):
        _, row, md = self.summarise()
        for col in ('pool_top_genome', 'pool_top_share', 'pool_other_max_genome', 'pool_other_max_share'):
            self.assertEqual(row[col], 'NA')
        self.assertNotIn('reference pool', md)

    def test_markdown_lists_genomes_by_reads_and_caps_the_table(self):
        pool = HEADER + ''.join(f"g{i}\t100\t{1000 - i}\t0\n" for i in range(12))
        _, _, md = self.summarise(pool=pool)
        self.assertIn('| `g0` |', md)
        self.assertNotIn('| `g8` |', md)
        self.assertIn('| 4 more genome(s) | | | |', md)


class TestKraken(Workdir):
    def test_unclassified_fraction_and_top_species(self):
        _, row, _ = self.summarise(kraken=KRAKEN)
        self.assertEqual(row['kraken_reads'], '50000')
        self.assertEqual(row['kraken_unclassified_fraction'], '0.4000')
        self.assertEqual(row['kraken_top_species'], 'Venezuelan equine encephalitis virus')
        self.assertEqual(row['kraken_top_species_reads'], '20000')

    def test_a_species_with_enough_reads_is_flagged(self):
        _, row, _ = self.summarise(kraken=KRAKEN)
        # 20000 of 1,000,000 reads is 2.0%, 5000 is 0.5%, 100 is 0.01%.
        self.assertIn('20,000 reads (2.0% of all) classify as Venezuelan equine encephalitis virus', row['flags'])
        self.assertNotIn('Escherichia coli', row['flags'])
        self.assertNotIn('Homo sapiens', row['flags'])

    def test_a_lower_flag_fraction_flags_more_species(self):
        _, row, _ = self.summarise(kraken=KRAKEN, extra=('--flag-fraction', '0.005'))
        self.assertIn('Escherichia coli', row['flags'])
        self.assertNotIn('Homo sapiens', row['flags'])

    def test_markdown_has_domains_and_top_species(self):
        _, _, md = self.summarise(kraken=KRAKEN)
        self.assertIn('Kraken2 classified 30,000 of the 50,000 reads that did not map to the stock.', md)
        self.assertIn('| Viruses | 25,000 |', md)
        self.assertIn('| Bacteria | 5,000 |', md)
        self.assertIn('| Escherichia coli | 5,000 | 0.50% |', md)

    def test_only_domains_are_listed_not_kingdoms_or_ranks_between(self):
        # In the viral database, Viruses (D) holds kingdoms (K) and ranks such as D1.
        report = (
            " 40.00\t20000\t20000\tU\t0\tunclassified\n"
            " 60.00\t30000\t0\tR\t1\troot\n"
            " 60.00\t30000\t1\tD\t10239\t  Viruses\n"
            " 50.00\t25000\t0\tD1\t2731341\t    Duplodnaviria\n"
            " 50.00\t25000\t0\tK\t2731360\t      Heunggongvirae\n"
        )
        _, _, md = self.summarise(kraken=report)
        self.assertIn('| Viruses | 30,000 |', md)
        self.assertNotIn('Duplodnaviria', md)
        self.assertNotIn('Heunggongvirae', md)

    def test_an_empty_report_means_no_reads_to_classify(self):
        report = " 100.00\t0\t0\tU\t0\tunclassified\n"
        _, row, md = self.summarise(kraken=report)
        self.assertEqual(row['kraken_reads'], '0')
        self.assertEqual(row['kraken_unclassified_fraction'], 'NA')
        self.assertIn('There were no reads left to classify.', md)

    def test_without_a_report_the_kraken_columns_are_NA(self):
        _, row, md = self.summarise()
        self.assertEqual(row['kraken_reads'], 'NA')
        self.assertEqual(row['kraken_top_species'], 'NA')
        self.assertNotIn('Kraken2', md)

    def test_malformed_report_lines_are_skipped(self):
        report = "not a report line\n" + KRAKEN + "bad\tline\n"
        _, row, _ = self.summarise(kraken=report)
        self.assertEqual(row['kraken_reads'], '50000')


class TestInputs(Workdir):
    def test_crlf_line_endings_are_read(self):
        _, row, _ = self.summarise(stock=STOCK.replace('\n', '\r\n'), pool=POOL.replace('\n', '\r\n'),
                                   kraken=KRAKEN.replace('\n', '\r\n'))
        self.assertEqual(row['reads_sampled'], '1000000')
        self.assertEqual(row['pool_top_genome'], 'strainA')
        self.assertEqual(row['kraken_reads'], '50000')

    def test_the_unmapped_row_may_come_first_or_last(self):
        stock = HEADER + "*\t0\t50000\t0\nstock\t11446\t950000\t0\n"
        _, row, _ = self.summarise(stock=stock)
        self.assertEqual(row['reads_unmapped_to_stock'], '50000')


class TestErrors(Workdir):
    def test_missing_stock_counts_fails(self):
        result = self.run_script('--stock-name', 'S', '--stock-counts', self.path('none.tsv'),
                                 '--output-tsv', self.path('o.tsv'), '--output-md', self.path('o.md'))
        self.assertEqual(result.returncode, 1)
        self.assertIn('file not found', result.stderr)

    def test_missing_optional_file_fails(self):
        result = self.run_script('--stock-name', 'S', '--stock-counts', self.write('s.tsv', STOCK),
                                 '--pool-counts', self.path('none.tsv'),
                                 '--output-tsv', self.path('o.tsv'), '--output-md', self.path('o.md'))
        self.assertEqual(result.returncode, 1)
        self.assertIn('file not found', result.stderr)

    def test_wrong_header_is_rejected(self):
        result, _, _ = self.summarise(stock="name\tlen\treads\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('expected a header of contig, length, mapped, mapq0', result.stderr)

    def test_empty_counts_file_is_rejected(self):
        result, _, _ = self.summarise(stock="")
        self.assertEqual(result.returncode, 1)
        self.assertIn('expected a header', result.stderr)

    def test_short_row_is_rejected(self):
        result, _, _ = self.summarise(stock=HEADER + "stock\t11446\t950000\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('expected 4 columns', result.stderr)

    def test_non_numeric_count_is_rejected(self):
        result, _, _ = self.summarise(stock=HEADER + "stock\t11446\tmany\t0\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('must be whole numbers', result.stderr)

    def test_mapq0_above_mapped_is_rejected(self):
        result, _, _ = self.summarise(stock=HEADER + "stock\t11446\t10\t20\n*\t0\t5\t0\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('mapq0 (20) is larger than mapped (10)', result.stderr)

    def test_no_reads_at_all_is_rejected(self):
        result, _, _ = self.summarise(stock=HEADER + "stock\t11446\t0\t0\n*\t0\t0\t0\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('hold no reads', result.stderr)

    def test_kinds_total_must_match_the_unmapped_count(self):
        result, _, _ = self.summarise(kinds="kind\treads\ntotal\t49999\npolyg\t0\nadapter\t0\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('counts 49999 unmapped reads, the stock counts 50000', result.stderr)

    def test_kinds_that_add_up_to_more_than_the_total_are_rejected(self):
        result, _, _ = self.summarise(kinds="kind\treads\ntotal\t50000\npolyg\t40000\nadapter\t20000\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('add up to more than the total', result.stderr)

    def test_kinds_with_a_wrong_header_are_rejected(self):
        result, _, _ = self.summarise(kinds="what\tn\ntotal\t50000\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('expected a header of kind, reads', result.stderr)

    def test_kinds_missing_a_row_are_rejected(self):
        result, _, _ = self.summarise(kinds="kind\treads\ntotal\t50000\npolyg\t0\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('no row for adapter', result.stderr)

    def test_kinds_with_a_non_numeric_count_are_rejected(self):
        result, _, _ = self.summarise(kinds="kind\treads\ntotal\tmany\npolyg\t0\nadapter\t0\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('reads must be a whole number', result.stderr)

    def test_a_missing_kinds_file_fails(self):
        result = self.run_script('--stock-name', 'S', '--stock-counts', self.write('s.tsv', STOCK),
                                 '--unmapped-kinds', self.path('none.tsv'),
                                 '--output-tsv', self.path('o.tsv'), '--output-md', self.path('o.md'))
        self.assertEqual(result.returncode, 1)
        self.assertIn('file not found', result.stderr)

    def test_bad_thresholds_are_rejected(self):
        for flag, value in (('--min-mapped', '0'), ('--min-mapped', '1.5'),
                            ('--flag-fraction', '0'), ('--flag-fraction', '2'),
                            ('--min-reads', '0')):
            result, _, _ = self.summarise(extra=(flag, value))
            self.assertEqual(result.returncode, 1, (flag, value))
            self.assertIn(flag, result.stderr)

    def test_failure_leaves_no_output_files(self):
        self.summarise(stock="bad")
        self.assertFalse(os.path.exists(self.path('out.tsv')))
        self.assertFalse(os.path.exists(self.path('out.md')))


if __name__ == '__main__':
    unittest.main()
