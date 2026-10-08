"""Unit tests for bin/ivar_variants_to_vcf.py.

Run from the repository root:  python -m unittest discover -s tests/unit -v

Every TSV here uses iVar's real column layout (ALT_DP after REF_QUAL, ALT_FREQ
after ALT_QUAL). Values are synthetic. Each numeric column gets a distinct value,
so a parser that reads the wrong column fails an assertion.
"""
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SCRIPT = os.path.join(REPO, 'bin', 'ivar_variants_to_vcf.py')
sys.path.insert(0, os.path.join(REPO, 'bin'))

import ivar_variants_to_vcf as conv  # noqa: E402

HEADER = ('REGION\tPOS\tREF\tALT\tREF_DP\tREF_RV\tREF_QUAL\tALT_DP\tALT_RV\tALT_QUAL\t'
          'ALT_FREQ\tTOTAL_DP\tPVAL\tPASS\tGFF_FEATURE\tREF_CODON\tREF_AA\tALT_CODON\t'
          'ALT_AA\tPOS_AA')


def row(pos=100, ref='T', alt='C', ref_dp='9000', alt_dp='1000', alt_freq='0.1',
        total_dp='10000', pval='0.001', passed='TRUE'):
    # REF_RV=11, REF_QUAL=22, ALT_RV=33, ALT_QUAL=44 are decoys for the wrong-column bug.
    return (f'KP282671.1\t{pos}\t{ref}\t{alt}\t{ref_dp}\t11\t22\t{alt_dp}\t33\t44\t'
            f'{alt_freq}\t{total_dp}\t{pval}\t{passed}\tCDS\tAAA\tK\tAGA\tR\t1')


class ParseTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def tsv(self, *lines, header=HEADER):
        path = os.path.join(self.tmp.name, 'sample.tsv')
        with open(path, 'w', encoding='utf-8', newline='') as f:
            f.write('\n'.join([header, *lines]) + '\n')
        return path

    def test_real_layout_reads_the_right_columns(self):
        v, = conv.parse_ivar_tsv(self.tsv(row()))
        self.assertEqual(v['REF_DP'], 9000)
        self.assertEqual(v['ALT_DP'], 1000)       # not REF_RV (11)
        self.assertEqual(v['ALT_FREQ'], 0.1)      # not ALT_QUAL (44)
        self.assertEqual(v['TOTAL_DP'], 10000)    # not ALT_FREQ
        self.assertEqual(v['PVAL'], 0.001)        # not TOTAL_DP
        self.assertEqual(v['PASS'], 'TRUE')       # not PVAL

    def test_scientific_notation_depth_is_kept(self):
        v, = conv.parse_ivar_tsv(self.tsv(row(
            ref_dp='1.09161e+06', alt_dp='72140', alt_freq='0.061982', total_dp='1.16389e+06')))
        self.assertEqual(v['REF_DP'], 1091610)
        self.assertEqual(v['TOTAL_DP'], 1163890)

    def test_column_order_does_not_matter(self):
        names = HEADER.split('\t')
        values = row().split('\t')
        order = list(reversed(range(len(names))))
        header = '\t'.join(names[i] for i in order)
        line = '\t'.join(values[i] for i in order)
        v, = conv.parse_ivar_tsv(self.tsv(line, header=header))
        self.assertEqual((v['ALT_DP'], v['ALT_FREQ'], v['TOTAL_DP']), (1000, 0.1, 10000))

    def test_missing_column_is_an_error(self):
        header = HEADER.replace('\tALT_FREQ', '\tFREQ_RENAMED')
        with self.assertRaisesRegex(ValueError, 'ALT_FREQ'):
            conv.parse_ivar_tsv(self.tsv(row(), header=header))

    def test_unparseable_value_is_an_error_not_zero(self):
        with self.assertRaisesRegex(ValueError, r'sample\.tsv:2'):
            conv.parse_ivar_tsv(self.tsv(row(alt_freq='NA')))

    def test_short_row_is_an_error(self):
        with self.assertRaisesRegex(ValueError, 'expected at least'):
            conv.parse_ivar_tsv(self.tsv('KP282671.1\t100\tT\tC'))

    def test_header_only_file_gives_no_variants(self):
        self.assertEqual(conv.parse_ivar_tsv(self.tsv()), [])

    def test_crlf_line_endings(self):
        path = os.path.join(self.tmp.name, 'crlf.tsv')
        with open(path, 'w', encoding='utf-8', newline='') as f:
            f.write(HEADER + '\r\n' + row() + '\r\n')
        v, = conv.parse_ivar_tsv(path)
        self.assertEqual(v['PASS'], 'TRUE')


class CliTests(unittest.TestCase):
    def run_cli(self, tsv_text):
        with tempfile.TemporaryDirectory() as d:
            tsv = os.path.join(d, 's1.tsv')
            ref = os.path.join(d, 'ref.fa')
            out = os.path.join(d, 's1.vcf')
            open(tsv, 'w', encoding='utf-8').write(tsv_text)
            open(ref, 'w', encoding='utf-8').write('>KP282671.1\nACGT\n')
            proc = subprocess.run(
                [sys.executable, SCRIPT, '--input-tsv', tsv, '--output-vcf', out,
                 '--reference-fasta', ref], capture_output=True, text=True)
            vcf = None
            if os.path.exists(out):
                with open(out, encoding='utf-8') as f:
                    vcf = f.read()
            return proc, vcf

    def test_vcf_record_carries_the_right_values(self):
        proc, vcf = self.run_cli(HEADER + '\n' + row(alt_freq='0.0125', passed='FALSE') + '\n')
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rec = [l for l in vcf.splitlines() if not l.startswith('#')][0].split('\t')
        self.assertEqual(rec[1:5], ['100', '.', 'T', 'C'])
        self.assertEqual(rec[6], 'FAIL')
        self.assertEqual(rec[7], 'DP=10000;AF=0.0125')
        self.assertEqual(rec[9], '1:10000:9000,1000:0.0125')

    def test_bad_input_exits_nonzero_without_a_vcf(self):
        proc, vcf = self.run_cli(HEADER + '\n' + row(total_dp='x') + '\n')
        self.assertEqual(proc.returncode, 1)
        self.assertIn('ERROR', proc.stderr)
        self.assertIsNone(vcf)

    def test_bad_flag_exits_2(self):
        proc = subprocess.run([sys.executable, SCRIPT, '--badflag'], capture_output=True, text=True)
        self.assertEqual(proc.returncode, 2)


if __name__ == '__main__':
    unittest.main()
