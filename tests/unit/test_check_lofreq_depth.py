"""Unit tests for bin/check_lofreq_depth.py.

Run from the repository root:  python -m unittest discover -s tests/unit -v

The script decides whether a LoFreq VAF came from all the reads at a position or from a
subset. A wrong flag either hides an inflated VAF or marks a sound one, so each test
states the expected flag and the expected fraction explicitly.

LoFreq can examine at most max_depth reads, so fraction_examined = min(1, max_depth / real
depth). With the cap at 1,000,000 (contig "c"):

  pos   LoFreq DP   real depth    fraction     expected
  10      900000      4400000      0.2273      capped
  20      990000      1000000      1.0         not capped, depth is at the cap
  30      800000      1000000      1.0         not capped, LoFreq's own DP is lower only
                                               because it drops reads on quality
  40       99999      1100000      0.9091      not capped, 10% over the cap
  45      900000      1250000      0.8         not capped, 0.8 is not under the 0.8 ratio
  46      900000      1250001      0.7999994   capped
  50         500            0      NA          unknown, real depth 0
  60         500      (missing)    NA          unknown, position not in the depth table
"""
import gzip
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SCRIPT = os.path.join(REPO, 'bin', 'check_lofreq_depth.py')
sys.path.insert(0, os.path.join(REPO, 'bin'))

import check_lofreq_depth as cld  # noqa: E402

SOURCE = ("##source=lofreq call-parallel -f ref.fasta --min-cov 10 --min-bq 30 --sig 0.01 "
          "--max-depth {cap} -o out.vcf in.bam\n")
COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
HEADER = "##fileformat=VCFv4.2\n" + SOURCE.format(cap=1000000) + COLS
HEADER_NO_CAP = "##fileformat=VCFv4.2\n##source=lofreq call -f ref.fasta --sig 0.01 -o out.vcf in.bam\n" + COLS

VCF = HEADER + (
    "c\t10\t.\tA\tG\t100\tPASS\tDP=900000;AF=0.018;DP4=440000,440000,8000,12000\n"
    "c\t20\t.\tC\tT\t100\tPASS\tDP=990000;AF=0.002;DP4=500000,488000,1000,1000\n"
    "c\t30\t.\tG\tA\t100\tPASS\tDP=800000;AF=0.050;DP4=400000,360000,20000,20000\n"
    "c\t40\t.\tT\tC\t100\tPASS\tDP=99999;AF=0.100\n"
    "c\t45\t.\tT\tC\t100\tPASS\tDP=900000;AF=0.100\n"
    "c\t46\t.\tT\tC\t100\tPASS\tDP=900000;AF=0.100\n"
    "c\t50\t.\tA\tT\t100\tPASS\tDP=500;AF=0.200\n"
    "c\t60\t.\tA\tT\t100\tPASS\tDP=500;AF=0.300\n"
)

DEPTH = (
    "c\t10\t4400000\n"
    "c\t20\t1000000\n"
    "c\t30\t1000000\n"
    "c\t40\t1100000\n"
    "c\t45\t1250000\n"
    "c\t46\t1250001\n"
    "c\t50\t0\n"
)

N_CALLS = 8


def by_pos(rows):
    return {int(r['pos']): r for r in rows}


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

    def check(self, vcf=VCF, depth=DEPTH, extra=()):
        result = self.run_script('--vcf', self.write('in.vcf', vcf),
                                 '--depth', self.write('depth.tsv', depth),
                                 '--output', self.path('out.tsv'), *extra)
        rows = []
        if os.path.exists(self.path('out.tsv')):
            with open(self.path('out.tsv'), encoding='utf-8') as f:
                lines = f.read().splitlines()
            header = lines[0].split('\t')
            rows = [dict(zip(header, line.split('\t'))) for line in lines[1:]]
        return result, rows


class TestFlags(Workdir):
    def test_one_row_per_call_with_the_documented_columns(self):
        result, rows = self.check()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(rows), N_CALLS)
        self.assertEqual(list(rows[0]), cld.COLUMNS)

    def test_flags_follow_the_table(self):
        _, rows = self.check()
        flags = {pos: r['depth_capped'] for pos, r in by_pos(rows).items()}
        self.assertEqual(flags, {10: 'yes', 20: 'no', 30: 'no', 40: 'no', 45: 'no',
                                 46: 'yes', 50: 'unknown', 60: 'unknown'})

    def test_fraction_examined_is_the_cap_over_the_real_depth(self):
        _, rows = self.check()
        r = by_pos(rows)
        self.assertEqual(r[10]['fraction_examined'], '0.2273')
        self.assertEqual(r[10]['real_depth'], '4400000')
        self.assertEqual(r[10]['max_depth'], '1000000')
        self.assertEqual(r[10]['lofreq_depth'], '900000')
        self.assertEqual(r[40]['fraction_examined'], '0.9091')

    def test_fraction_is_one_when_the_depth_is_under_the_cap(self):
        _, rows = self.check()
        r = by_pos(rows)
        self.assertEqual(r[20]['fraction_examined'], '1.0000')
        self.assertEqual(r[30]['fraction_examined'], '1.0000')

    def test_low_lofreq_dp_alone_does_not_flag_a_call(self):
        # LoFreq drops reads on base and mapping quality, so its DP sits under the real depth
        # at positions well below the cap. That is not the cap at work.
        vcf = HEADER + "c\t10\t.\tA\tG\t100\tPASS\tDP=470000;AF=0.1\n"
        _, rows = self.check(vcf=vcf, depth="c\t10\t688000\n")
        self.assertEqual(rows[0]['depth_capped'], 'no')

    def test_unknown_rows_have_no_fraction(self):
        _, rows = self.check()
        r = by_pos(rows)
        self.assertEqual(r[50]['fraction_examined'], 'NA')
        self.assertEqual(r[50]['real_depth'], '0')
        self.assertEqual(r[60]['fraction_examined'], 'NA')
        self.assertEqual(r[60]['real_depth'], 'NA')

    def test_ratio_is_a_strict_lower_bound(self):
        # 1000000 / 1250000 is exactly 0.8, so it is not under the ratio.
        _, rows = self.check()
        self.assertEqual(by_pos(rows)[45]['depth_capped'], 'no')
        _, rows = self.check(extra=('--ratio', '0.81'))
        self.assertEqual(by_pos(rows)[45]['depth_capped'], 'yes')

    def test_ratio_of_one_flags_any_position_over_the_cap(self):
        _, rows = self.check(extra=('--ratio', '1'))
        flags = {pos: r['depth_capped'] for pos, r in by_pos(rows).items()}
        self.assertEqual(flags[40], 'yes')
        self.assertEqual(flags[20], 'no')


class TestCap(Workdir):
    def test_cap_is_read_from_the_vcf_header(self):
        vcf = "##fileformat=VCFv4.2\n" + SOURCE.format(cap=30) + COLS + "c\t10\t.\tA\tG\t1\tPASS\tDP=20;AF=0.1\n"
        result, rows = self.check(vcf=vcf, depth="c\t10\t60\n")
        self.assertEqual(rows[0]['max_depth'], '30')
        self.assertEqual(rows[0]['fraction_examined'], '0.5000')
        self.assertEqual(rows[0]['depth_capped'], 'yes')
        self.assertIn('LoFreq depth cap: 30 (from the VCF header).', result.stdout)

    def test_header_with_equals_sign_is_read(self):
        vcf = ("##fileformat=VCFv4.2\n##source=lofreq call --max-depth=30 -o x.vcf y.bam\n" + COLS +
               "c\t10\t.\tA\tG\t1\tPASS\tDP=20;AF=0.1\n")
        _, rows = self.check(vcf=vcf, depth="c\t10\t60\n")
        self.assertEqual(rows[0]['max_depth'], '30')

    def test_flag_overrides_the_header(self):
        result, rows = self.check(extra=('--max-depth', '5000000'))
        self.assertEqual(by_pos(rows)[10]['max_depth'], '5000000')
        self.assertEqual(by_pos(rows)[10]['depth_capped'], 'no')
        self.assertIn('LoFreq depth cap: 5,000,000 (from --max-depth).', result.stdout)

    def test_lofreq_default_is_used_when_the_header_has_no_cap(self):
        # lofreq call-parallel records --max-depth only when it was given.
        vcf = HEADER_NO_CAP + "c\t10\t.\tA\tG\t1\tPASS\tDP=900000;AF=0.1\n"
        result, rows = self.check(vcf=vcf, depth="c\t10\t4000000\n")
        self.assertEqual(rows[0]['max_depth'], '1000000')
        self.assertEqual(rows[0]['depth_capped'], 'yes')
        self.assertIn("LoFreq's default, because the VCF header has none", result.stdout)

    def test_only_a_source_line_counts_as_the_header_cap(self):
        vcf = ("##fileformat=VCFv4.2\n##note=ran with --max-depth 5\n" + HEADER_NO_CAP[len("##fileformat=VCFv4.2\n"):] +
               "c\t10\t.\tA\tG\t1\tPASS\tDP=900000;AF=0.1\n")
        _, rows = self.check(vcf=vcf, depth="c\t10\t4000000\n")
        self.assertEqual(rows[0]['max_depth'], '1000000')


class TestVaf(Workdir):
    def test_vaf_comes_from_dp4_when_present(self):
        # DP4 counts 8000 + 12000 alt reads out of 900000 reads in all.
        _, rows = self.check()
        self.assertEqual(by_pos(rows)[10]['lofreq_vaf'], '0.022222')

    def test_vaf_falls_back_to_af_without_dp4(self):
        _, rows = self.check()
        self.assertEqual(by_pos(rows)[40]['lofreq_vaf'], '0.100000')

    def test_vaf_falls_back_to_af_when_dp4_is_all_zero(self):
        vcf = HEADER + "c\t10\t.\tA\tG\t1\tPASS\tDP=10;AF=0.25;DP4=0,0,0,0\n"
        _, rows = self.check(vcf=vcf, depth="c\t10\t10\n")
        self.assertEqual(rows[0]['lofreq_vaf'], '0.250000')


class TestSummary(Workdir):
    def test_summary_counts_and_lists_capped_calls_at_or_above_min_vaf(self):
        result, _ = self.check()
        self.assertIn('LoFreq depth cap: 1,000,000 (from the VCF header).', result.stdout)
        self.assertIn('8 call(s): 2 at positions where the cap let LoFreq examine under 80% of the reads, '
                      '2 with no depth to compare.', result.stdout)
        self.assertIn('c:10 A>G', result.stdout)
        self.assertIn('4,400,000 reads, LoFreq examined at most 1,000,000 (23%)', result.stdout)
        self.assertIn('c:46 T>C', result.stdout)

    def test_min_vaf_hides_low_vaf_calls_from_the_listing(self):
        # Position 10 has a VAF of 2.2% and position 46 has 10%.
        result, _ = self.check(extra=('--min-vaf', '0.05'))
        self.assertIn('c:46 T>C', result.stdout)
        self.assertNotIn('c:10 A>G', result.stdout)

    def test_calls_that_are_not_capped_are_never_listed(self):
        result, _ = self.check(extra=('--min-vaf', '0'))
        self.assertNotIn('c:20', result.stdout)
        self.assertNotIn('c:30', result.stdout)

    def test_no_capped_calls_gives_no_warning(self):
        vcf = HEADER + "c\t20\t.\tC\tT\t100\tPASS\tDP=990000;AF=0.002\n"
        result, rows = self.check(vcf=vcf)
        self.assertEqual(rows[0]['depth_capped'], 'no')
        self.assertNotIn('subset of the reads', result.stdout)

    def test_warning_names_the_parameter_to_raise(self):
        result, _ = self.check()
        self.assertIn('--lofreq_max_depth', result.stdout)


class TestInputs(Workdir):
    def test_empty_vcf_gives_header_only_output(self):
        result, rows = self.check(vcf=HEADER)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(rows, [])
        with open(self.path('out.tsv'), encoding='utf-8') as f:
            self.assertEqual(f.read(), '\t'.join(cld.COLUMNS) + '\n')
        self.assertIn('0 call(s)', result.stdout)

    def test_gzip_inputs_are_read(self):
        vcf = self.path('in.vcf.gz')
        depth = self.path('depth.tsv.gz')
        with gzip.open(vcf, 'wt', encoding='utf-8') as f:
            f.write(VCF)
        with gzip.open(depth, 'wt', encoding='utf-8') as f:
            f.write(DEPTH)
        result = self.run_script('--vcf', vcf, '--depth', depth, '--output', self.path('out.tsv'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(f'{N_CALLS} call(s)', result.stdout)

    def test_scientific_notation_depth_is_read(self):
        # samtools can print a large depth as 1.16389e+06.
        vcf = HEADER + "c\t10\t.\tA\tG\t100\tPASS\tDP=900000;AF=0.1\n"
        _, rows = self.check(vcf=vcf, depth="c\t10\t4.2e+06\n")
        self.assertEqual(rows[0]['real_depth'], '4200000')
        self.assertEqual(rows[0]['depth_capped'], 'yes')

    def test_crlf_line_endings_are_read(self):
        result, rows = self.check(vcf=VCF.replace('\n', '\r\n'), depth=DEPTH.replace('\n', '\r\n'))
        self.assertEqual(len(rows), N_CALLS)
        self.assertEqual(by_pos(rows)[10]['depth_capped'], 'yes')
        self.assertIn('1,000,000 (from the VCF header)', result.stdout)

    def test_blank_lines_in_the_depth_table_are_ignored(self):
        _, rows = self.check(depth=DEPTH + "\n\n")
        self.assertEqual(len(rows), N_CALLS)


class TestErrors(Workdir):
    def test_missing_vcf_fails(self):
        result = self.run_script('--vcf', self.path('none.vcf'),
                                 '--depth', self.write('depth.tsv', DEPTH),
                                 '--output', self.path('out.tsv'))
        self.assertEqual(result.returncode, 1)
        self.assertIn('file not found', result.stderr)

    def test_missing_depth_fails(self):
        result = self.run_script('--vcf', self.write('in.vcf', VCF),
                                 '--depth', self.path('none.tsv'),
                                 '--output', self.path('out.tsv'))
        self.assertEqual(result.returncode, 1)
        self.assertIn('file not found', result.stderr)

    def test_vcf_without_dp_is_rejected(self):
        result, _ = self.check(vcf=HEADER + "c\t10\t.\tA\tG\t100\tPASS\tAF=0.1\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('no INFO/DP', result.stderr)

    def test_short_vcf_record_is_rejected(self):
        result, _ = self.check(vcf=HEADER + "c\t10\t.\tA\tG\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('8 columns', result.stderr)

    def test_short_depth_line_is_rejected(self):
        result, _ = self.check(depth="c\t10\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('contig, position and depth', result.stderr)

    def test_non_numeric_depth_is_rejected(self):
        result, _ = self.check(depth="c\t10\tmany\n")
        self.assertEqual(result.returncode, 1)
        self.assertIn('must be numbers', result.stderr)

    def test_bad_ratio_is_rejected(self):
        for ratio in ('0', '-1', '1.5'):
            result, _ = self.check(extra=('--ratio', ratio))
            self.assertEqual(result.returncode, 1, ratio)
            self.assertIn('--ratio', result.stderr)

    def test_bad_max_depth_is_rejected(self):
        for cap in ('0', '-5'):
            result, _ = self.check(extra=('--max-depth', cap))
            self.assertEqual(result.returncode, 1, cap)
            self.assertIn('--max-depth', result.stderr)

    def test_failure_leaves_no_output_file(self):
        self.check(vcf=HEADER + "c\t10\t.\tA\tG\t100\tPASS\tAF=0.1\n")
        self.assertFalse(os.path.exists(self.path('out.tsv')))


if __name__ == '__main__':
    unittest.main()
