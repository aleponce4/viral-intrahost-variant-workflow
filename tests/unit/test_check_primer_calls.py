"""Unit tests for bin/check_primer_calls.py.

Run from the repository root:  python -m unittest discover -s tests/unit -v

The script decides whether a call sits under a primer. A wrong answer either hides a call that
needs checking or sends someone chasing a sound one, so each test states the expected answer.

A primer BED is 0-based and half-open, a VCF POS is 1-based. The fixture scheme below is built
so the boundaries are explicit:

  BED interval 10 20  covers 1-based positions 11 to 20
  BED interval 50 60  covers 1-based positions 51 to 60

So position 10 is outside, 11 is the first inside, 20 is the last inside, and 21 is outside.
"""
import gzip
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SCRIPT = os.path.join(REPO, 'bin', 'check_primer_calls.py')
sys.path.insert(0, os.path.join(REPO, 'bin'))

import check_primer_calls as cpc  # noqa: E402

BED = (
    "c\t10\t20\tP1_LEFT\t0\t+\n"
    "c\t50\t60\tP1_RIGHT\t0\t-\n"
)

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

VCF = HEADER + (
    "c\t10\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.05\n"      # just outside
    "c\t11\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.05\n"      # first inside
    "c\t20\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.05\n"      # last inside
    "c\t21\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.05\n"      # just outside
    "c\t55\t.\tA\tG\t100\tPASS\tDP=200;AF=0.02\n"       # inside the second primer
    "c\t99\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.30\n"      # far from any primer
)

# Median of these is 1000, so a call at DP 200 is 0.20x the median.
DEPTH = "".join("c\t%d\t1000\n" % p for p in range(1, 100))


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

    def check(self, vcf=VCF, bed=BED, depth=None, extra=()):
        args = ['--vcf', self.write('in.vcf', vcf),
                '--primer-bed', self.write('primers.bed', bed),
                '--output', self.path('out.tsv')]
        if depth is not None:
            args += ['--depth', self.write('depth.tsv', depth)]
        result = self.run_script(*args, *extra)
        rows = []
        if os.path.exists(self.path('out.tsv')):
            with open(self.path('out.tsv'), encoding='utf-8') as f:
                lines = f.read().splitlines()
            header = lines[0].split('\t')
            rows = [dict(zip(header, line.split('\t'))) for line in lines[1:]]
        return result, rows


def by_pos(rows):
    return {int(r['pos']): r for r in rows}


class TestBoundaries(Workdir):
    def test_output_has_the_documented_columns(self):
        result, rows = self.check()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(list(rows[0]), cpc.COLUMNS)
        self.assertEqual(len(rows), 6)

    def test_bed_half_open_coordinates_are_respected(self):
        _, rows = self.check()
        flags = {pos: r['under_primer'] for pos, r in by_pos(rows).items()}
        self.assertEqual(flags[10], 'no')
        self.assertEqual(flags[11], 'yes')
        self.assertEqual(flags[20], 'yes')
        self.assertEqual(flags[21], 'no')

    def test_the_covering_primer_is_named(self):
        _, rows = self.check()
        r = by_pos(rows)
        self.assertEqual(r[11]['primers'], 'P1_LEFT')
        self.assertEqual(r[55]['primers'], 'P1_RIGHT')
        self.assertEqual(r[99]['primers'], 'NA')

    def test_overlapping_primers_are_all_named(self):
        bed = BED + "c\t15\t25\tP2_LEFT\t0\t+\n"
        _, rows = self.check(bed=bed)
        self.assertEqual(by_pos(rows)[20]['primers'], 'P1_LEFT,P2_LEFT')

    def test_a_call_on_another_contig_is_not_marked(self):
        vcf = HEADER + "other\t15\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.05\n"
        result, rows = self.check(vcf=vcf)
        self.assertEqual(rows[0]['under_primer'], 'no')
        self.assertIn('so none of those could be marked', result.stdout)


class TestDepth(Workdir):
    def test_depth_is_compared_with_the_median(self):
        _, rows = self.check(depth=DEPTH)
        r = by_pos(rows)
        self.assertEqual(r[11]['depth_vs_median'], '1.00')
        self.assertEqual(r[55]['depth_vs_median'], '0.20')

    def test_without_a_depth_table_the_ratio_is_NA(self):
        _, rows = self.check()
        self.assertEqual(by_pos(rows)[11]['depth_vs_median'], 'NA')

    def test_the_summary_gives_the_ratio_when_it_has_one(self):
        result, _ = self.check(depth=DEPTH, extra=('--min-vaf', '0.01'))
        self.assertIn('0.20x the median depth', result.stdout)


class TestSummary(Workdir):
    def test_counts_are_reported(self):
        result, _ = self.check()
        self.assertIn('6 call(s): 3 under a primer.', result.stdout)

    def test_calls_under_a_primer_are_listed_above_min_vaf(self):
        result, _ = self.check()
        self.assertIn('Under a primer at a VAF of 1.00% or more: 3', result.stdout)
        self.assertIn('c:11 A>G', result.stdout)
        self.assertIn('c:55 A>G', result.stdout)
        self.assertNotIn('c:99', result.stdout)

    def test_min_vaf_narrows_the_listing(self):
        result, _ = self.check(extra=('--min-vaf', '0.03'))
        self.assertIn('Under a primer at a VAF of 3.00% or more: 2', result.stdout)
        self.assertNotIn('c:55 A>G', result.stdout)

    def test_no_calls_under_a_primer_gives_no_listing(self):
        vcf = HEADER + "c\t99\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.30\n"
        result, rows = self.check(vcf=vcf)
        self.assertEqual(rows[0]['under_primer'], 'no')
        self.assertIn('1 call(s): 0 under a primer.', result.stdout)
        self.assertNotIn('check these against', result.stdout)

    def test_the_advice_names_what_to_do(self):
        result, _ = self.check()
        self.assertIn('check these against an untrimmed alignment', result.stdout)


class TestVaf(Workdir):
    def test_vaf_comes_from_dp4_when_present(self):
        vcf = HEADER + "c\t11\t.\tA\tG\t100\tPASS\tDP=1000;AF=0.5;DP4=400,400,100,100\n"
        _, rows = self.check(vcf=vcf)
        self.assertEqual(rows[0]['vaf'], '0.200000')

    def test_vaf_falls_back_to_af(self):
        _, rows = self.check()
        self.assertEqual(by_pos(rows)[11]['vaf'], '0.050000')

    def test_a_record_without_dp_still_works(self):
        vcf = HEADER + "c\t11\t.\tA\tG\t100\tPASS\tAF=0.1\n"
        result, rows = self.check(vcf=vcf)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows[0]['depth'], '0')


class TestInputs(Workdir):
    def test_gzip_inputs_are_read(self):
        vcf, bed = self.path('in.vcf.gz'), self.path('p.bed.gz')
        with gzip.open(vcf, 'wt', encoding='utf-8') as f:
            f.write(VCF)
        with gzip.open(bed, 'wt', encoding='utf-8') as f:
            f.write(BED)
        result = self.run_script('--vcf', vcf, '--primer-bed', bed, '--output', self.path('out.tsv'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('6 call(s)', result.stdout)

    def test_bed_track_and_comment_lines_are_skipped(self):
        bed = "# a comment\ntrack name=primers\n" + BED
        _, rows = self.check(bed=bed)
        self.assertEqual(by_pos(rows)[11]['under_primer'], 'yes')

    def test_crlf_line_endings_are_read(self):
        _, rows = self.check(vcf=VCF.replace('\n', '\r\n'), bed=BED.replace('\n', '\r\n'))
        self.assertEqual(len(rows), 6)
        self.assertEqual(by_pos(rows)[11]['under_primer'], 'yes')

    def test_an_empty_vcf_gives_a_header_only_file(self):
        result, rows = self.check(vcf=HEADER)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows, [])
        self.assertIn('0 call(s)', result.stdout)


class TestErrors(Workdir):
    def test_missing_vcf_fails(self):
        result = self.run_script('--vcf', self.path('none.vcf'),
                                 '--primer-bed', self.write('p.bed', BED),
                                 '--output', self.path('o.tsv'))
        self.assertEqual(result.returncode, 1)
        self.assertIn('file not found', result.stderr)

    def test_missing_bed_fails(self):
        result = self.run_script('--vcf', self.write('in.vcf', VCF),
                                 '--primer-bed', self.path('none.bed'),
                                 '--output', self.path('o.tsv'))
        self.assertEqual(result.returncode, 1)
        self.assertIn('file not found', result.stderr)

    def test_a_short_bed_line_is_rejected(self):
        result, _ = self.check(bed="c\t10\t20\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('at least 4 columns', result.stderr)

    def test_a_non_numeric_bed_coordinate_is_rejected(self):
        result, _ = self.check(bed="c\tten\t20\tP1\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('must be whole numbers', result.stderr)

    def test_an_inverted_bed_interval_is_rejected(self):
        result, _ = self.check(bed="c\t20\t10\tP1\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('must be greater than start', result.stderr)

    def test_an_empty_bed_is_rejected(self):
        result, _ = self.check(bed="# only a comment\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('no primer intervals found', result.stderr)

    def test_a_short_vcf_record_is_rejected(self):
        result, _ = self.check(vcf=HEADER + "c\t11\t.\tA\tG\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('8 columns', result.stderr)

    def test_failure_leaves_no_output_file(self):
        self.check(bed="c\t20\t10\tP1\n")
        self.assertFalse(os.path.exists(self.path('out.tsv')))


if __name__ == '__main__':
    unittest.main()
