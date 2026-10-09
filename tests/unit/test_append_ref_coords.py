"""Unit tests for bin/append_ref_coords.py.

Run from the repository root:  python -m unittest discover -s tests/unit -v

This script decides which lab-reference position a stock call is reported at, so
a wrong lookup moves a variant onto the wrong base without any error. Each test
states the expected reference position explicitly rather than deriving it.

The table below is a stock aligned to a reference that starts 8 bases earlier,
with a substitution at stock 3, a 1 bp insertion at stock 4 and a 1 bp deletion
between stock 5 and 6:

  stock    1   2   3   4    5   (del)  6   7
  ref      9  10  11  NA   12    13   14  15
  status   m   m  mis  ins   m   del   m   m
"""
import gzip
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SCRIPT = os.path.join(REPO, 'bin', 'append_ref_coords.py')
sys.path.insert(0, os.path.join(REPO, 'bin'))

import append_ref_coords as arc  # noqa: E402

TABLE = (
    "stock_contig\tstock_pos\tref_contig\tref_pos\tstatus\tstock_base\tref_base\n"
    "stock\t1\tref\t9\tmatch\t.\t.\n"
    "stock\t2\tref\t10\tmatch\t.\t.\n"
    "stock\t3\tref\t11\tmismatch\tC\tT\n"
    "stock\t4\tref\tNA\tinsertion\tG\tNA\n"
    "stock\t5\tref\t12\tmatch\t.\t.\n"
    "stock\tNA\tref\t13\tdeletion\tNA\tA\n"
    "stock\t6\tref\t14\tmatch\t.\t.\n"
    "stock\t7\tref\t15\tmatch\t.\t.\n"
)

VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=stock,length=7>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "stock\t1\t.\tA\tG\t.\tPASS\tDP=100;AF=0.5\n"
    "stock\t3\t.\tC\tT\t.\tPASS\tDP=100;AF=0.1\n"
    "stock\t4\t.\tG\tA\t.\tPASS\t.\n"
    "stock\t5\t.\tAC\tA\t.\tPASS\tDP=90\n"
    "stock\t7\t.\tA\tC\t.\tPASS\tDP=80\n"
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
        return subprocess.run([sys.executable, SCRIPT, *args],
                              capture_output=True, text=True)

    def read(self, path):
        with open(path, encoding='utf-8') as f:
            return f.read()

    def info(self, vcf_text):
        """Return {pos: INFO string} for the records of a VCF."""
        out = {}
        for line in vcf_text.splitlines():
            if line and not line.startswith('#'):
                parts = line.split('\t')
                out[int(parts[1])] = parts[7]
        return out


class ReadLiftoverTests(Workdir):
    def test_deletion_rows_are_not_keyed_by_stock_position(self):
        table = arc.read_liftover(self.write('t.tsv', TABLE))
        self.assertEqual(sorted(table['stock']), [1, 2, 3, 4, 5, 6, 7])

    def test_insertion_has_no_reference_position(self):
        table = arc.read_liftover(self.write('t.tsv', TABLE))
        self.assertEqual(table['stock'][4], ('ref', None, 'insertion'))

    def test_missing_column_is_named(self):
        bad = self.write('bad.tsv', "stock_contig\tstock_pos\nstock\t1\n")
        with self.assertRaisesRegex(ValueError, 'ref_contig'):
            arc.read_liftover(bad)


class LiftTests(Workdir):
    def setUp(self):
        super().setUp()
        self.cmap = arc.read_liftover(self.write('t.tsv', TABLE))['stock']

    def test_first_base_shifts_by_the_five_prime_gap(self):
        self.assertEqual(arc.lift(self.cmap, 1, 1), ('ref', 9, None, 'match'))

    def test_substitution_keeps_its_reference_position(self):
        self.assertEqual(arc.lift(self.cmap, 3, 1), ('ref', 11, None, 'mismatch'))

    def test_base_after_a_deletion_skips_the_deleted_reference_base(self):
        # stock 6 is ref 14: ref 13 was deleted from the stock.
        self.assertEqual(arc.lift(self.cmap, 6, 1), ('ref', 14, None, 'match'))

    def test_base_after_an_insertion_does_not_shift(self):
        self.assertEqual(arc.lift(self.cmap, 5, 1), ('ref', 12, None, 'match'))

    def test_inserted_base_has_status_but_no_position(self):
        self.assertEqual(arc.lift(self.cmap, 4, 1), ('ref', None, None, 'insertion'))

    def test_multibase_ref_reports_the_last_base(self):
        # REF spans stock 5-6, which are ref 12 and 14.
        self.assertEqual(arc.lift(self.cmap, 5, 2), ('ref', 12, 14, 'match'))

    def test_position_outside_the_block_is_unmapped(self):
        self.assertEqual(arc.lift(self.cmap, 99, 1), (None, None, None, 'unmapped'))


class VcfTests(Workdir):
    def annotate(self, vcf_text=VCF, name='in.vcf'):
        src = self.write(name, vcf_text)
        dst = self.path('out.vcf')
        result = self.run_script('--liftover-tsv', self.write('t.tsv', TABLE),
                                 '--input', src, '--output', dst)
        return result, dst

    def test_each_call_gets_its_reference_position(self):
        result, dst = self.annotate()
        self.assertEqual(result.returncode, 0, result.stderr)
        info = self.info(self.read(dst))
        self.assertEqual(info[1], 'DP=100;AF=0.5;REF_CONTIG=ref;REF_POS=9;REF_STATUS=match')
        self.assertEqual(info[3], 'DP=100;AF=0.1;REF_CONTIG=ref;REF_POS=11;REF_STATUS=mismatch')
        self.assertEqual(info[7], 'DP=80;REF_CONTIG=ref;REF_POS=15;REF_STATUS=match')

    def test_inserted_base_has_no_ref_pos_and_an_empty_info_is_replaced(self):
        _, dst = self.annotate()
        info = self.info(self.read(dst))
        self.assertEqual(info[4], 'REF_CONTIG=ref;REF_STATUS=insertion')

    def test_deletion_record_gets_ref_end(self):
        _, dst = self.annotate()
        info = self.info(self.read(dst))
        self.assertEqual(info[5], 'DP=90;REF_CONTIG=ref;REF_POS=12;REF_END=14;REF_STATUS=match')

    def test_header_declares_the_new_fields_before_the_chrom_line(self):
        _, dst = self.annotate()
        lines = self.read(dst).splitlines()
        chrom_at = next(i for i, l in enumerate(lines) if l.startswith('#CHROM'))
        declared = [l for l in lines[:chrom_at] if l.startswith('##INFO=<ID=REF_')]
        self.assertEqual([l.split(',')[0] for l in declared],
                         ['##INFO=<ID=REF_CONTIG', '##INFO=<ID=REF_POS',
                          '##INFO=<ID=REF_END', '##INFO=<ID=REF_STATUS'])

    def test_calls_are_not_altered(self):
        _, dst = self.annotate()
        before = [l.split('\t')[:7] for l in VCF.splitlines() if not l.startswith('#')]
        after = [l.split('\t')[:7] for l in self.read(dst).splitlines()
                 if not l.startswith('#')]
        self.assertEqual(before, after)

    def test_gzip_input_is_read(self):
        src = self.path('in.vcf.gz')
        with gzip.open(src, 'wt', encoding='utf-8') as f:
            f.write(VCF)
        dst = self.path('out.vcf')
        result = self.run_script('--liftover-tsv', self.write('t.tsv', TABLE),
                                 '--input', src, '--output', dst)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('REF_POS=9', self.read(dst))

    def test_summary_counts_each_status(self):
        result, _ = self.annotate()
        self.assertIn('Annotated 5 record(s): 1 insertion, 3 match, 1 mismatch', result.stdout)

    def test_wrong_table_fails_and_names_the_contigs(self):
        vcf = VCF.replace('stock\t1\t', 'other\t1\t')
        result, _ = self.annotate(vcf)
        self.assertEqual(result.returncode, 1)
        self.assertIn("'other'", result.stderr)
        self.assertIn('stock', result.stderr)

    def test_annotating_twice_is_refused(self):
        _, dst = self.annotate()
        again = self.run_script('--liftover-tsv', self.path('t.tsv'),
                                '--input', dst, '--output', self.path('out2.vcf'))
        self.assertEqual(again.returncode, 1)
        self.assertIn('annotated before', again.stderr)


class TableTests(Workdir):
    SUMMARY = ("sample\tchrom\tpos\tref\talt\tvaf_percent\n"
               "s1\tstock\t1\tA\tG\t5.00\n"
               "s1\tstock\t4\tG\tA\t2.00\n"
               "s2\tstock\t5\tAC\tA\t1.00\n"
               "s2\tstock\t50\tT\tC\t0.50\n")

    def annotate(self, text=None, *extra):
        src = self.write('summary.tsv', text or self.SUMMARY)
        dst = self.path('out.tsv')
        result = self.run_script('--liftover-tsv', self.write('t.tsv', TABLE),
                                 '--input', src, '--output', dst, *extra)
        return result, dst

    def test_columns_are_appended_with_na_for_missing_values(self):
        result, dst = self.annotate()
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = [l.split('\t') for l in self.read(dst).splitlines()]
        self.assertEqual(rows[0][-4:], ['ref_contig', 'ref_pos', 'ref_end', 'ref_status'])
        self.assertEqual(rows[1][-4:], ['ref', '9', 'NA', 'match'])
        self.assertEqual(rows[2][-4:], ['ref', 'NA', 'NA', 'insertion'])
        self.assertEqual(rows[3][-4:], ['ref', '12', '14', 'match'])
        self.assertEqual(rows[4][-4:], ['NA', 'NA', 'NA', 'unmapped'])

    def test_original_columns_are_untouched(self):
        _, dst = self.annotate()
        rows = [l.split('\t') for l in self.read(dst).splitlines()]
        self.assertEqual([r[:6] for r in rows],
                         [l.split('\t') for l in self.SUMMARY.splitlines()])

    def test_custom_column_names(self):
        text = "id\tcontig\tposition\n1\tstock\t3\n"
        result, dst = self.annotate(text, '--chrom-col', 'contig', '--pos-col', 'position')
        self.assertEqual(result.returncode, 0, result.stderr)
        last = self.read(dst).splitlines()[1].split('\t')
        self.assertEqual(last[-4:], ['ref', '11', 'NA', 'mismatch'])

    def test_missing_position_column_is_named(self):
        result, _ = self.annotate("sample\tchrom\nx\tstock\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn("'pos'", result.stderr)

    def test_unknown_extension_needs_an_explicit_format(self):
        src = self.write('summary.dat', self.SUMMARY)
        result = self.run_script('--liftover-tsv', self.write('t.tsv', TABLE),
                                 '--input', src, '--output', self.path('o.tsv'))
        self.assertEqual(result.returncode, 1)
        self.assertIn('--format', result.stderr)


if __name__ == '__main__':
    unittest.main()
