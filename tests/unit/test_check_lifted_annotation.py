"""Unit tests for bin/check_lifted_annotation.py.

Run from the repository root:  python -m unittest discover -s tests/unit -v

The check exists to stop a broken reference reaching variant calling, so the
tests cover both halves of that: a real break must fail, and a problem the
reference already had must not.
"""
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SCRIPT = os.path.join(REPO, 'bin', 'check_lifted_annotation.py')
sys.path.insert(0, os.path.join(REPO, 'bin'))

import check_lifted_annotation as cla  # noqa: E402

# A 60 bp CDS with no internal stop, then the same with an opal TGA at codon 3.
CLEAN_CDS = 'ATG' + 'AAG' * 18 + 'TAA'
OPAL_CDS = 'ATG' + 'AAG' * 2 + 'TGA' + 'AAG' * 15 + 'TAA'
FLANK = 'CCCCCCCCCC'


def gff(seqid='stock', start=11, end=70, name='nsP1', strand='+'):
    return ('##gff-version 3\n'
            f'{seqid}\tTest\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=gene:{name};Name={name}\n'
            f'{seqid}\tTest\tCDS\t{start}\t{end}\t.\t{strand}\t0\t'
            f'ID=cds:{name};Parent=gene:{name};Name={name}\n')


class Case:
    """A temporary stock/reference pair, written to disk."""

    def __init__(self, tmp, stock_cds, ref_cds, seqid='stock', stock_kw=None, **gff_kw):
        # stock_kw changes the stock's GFF only, so a test can break the
        # transferred annotation while the reference control stays sound.
        self.dir = tmp
        self.stock_fasta = self._fasta('stock.fasta', seqid, FLANK + stock_cds + FLANK)
        self.ref_fasta = self._fasta('ref.fasta', 'ref', FLANK + ref_cds + FLANK)
        self.stock_gff = self._gff('stock.gff3', gff(seqid=seqid, **{**gff_kw, **(stock_kw or {})}))
        self.ref_gff = self._gff('ref.gff3', gff(seqid='ref', **gff_kw))
        self.out = os.path.join(tmp, 'check.tsv')

    def _fasta(self, name, seqid, seq):
        p = os.path.join(self.dir, name)
        with open(p, 'w', encoding='utf-8', newline='') as f:
            f.write(f'>{seqid} test\n{seq}\n')
        return p

    def _gff(self, name, text):
        p = os.path.join(self.dir, name)
        with open(p, 'w', encoding='utf-8', newline='') as f:
            f.write(text)
        return p

    def run(self, *extra, control=True):
        cmd = [sys.executable, SCRIPT, '--fasta', self.stock_fasta, '--gff', self.stock_gff,
               '--output-tsv', self.out]
        if control:
            cmd += ['--reference-fasta', self.ref_fasta, '--reference-gff', self.ref_gff]
        cmd += list(extra)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        rows = []
        if os.path.exists(self.out):
            with open(self.out, encoding='utf-8') as f:
                lines = f.read().splitlines()
            keys = lines[0].split('\t')
            rows = [dict(zip(keys, l.split('\t'))) for l in lines[1:]]
        return proc, rows


class CheckTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def case(self, stock_cds, ref_cds, **kw):
        return Case(self.tmp.name, stock_cds, ref_cds, **kw)

    def test_clean_transfer_passes(self):
        proc, rows = self.case(CLEAN_CDS, CLEAN_CDS).run()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rows[0]['in_contig'], 'pass')
        self.assertEqual(rows[0]['length_mod3'], 'pass')
        self.assertEqual(rows[0]['no_internal_stop'], 'pass')

    def test_a_stop_the_reference_also_has_is_inherited_not_a_failure(self):
        # The alphavirus nsP3 opal codon in miniature: a real stop inside a real
        # CDS, carried over from the reference rather than caused by the transfer.
        proc, rows = self.case(OPAL_CDS, OPAL_CDS).run()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rows[0]['no_internal_stop'], 'inherited')
        self.assertIn('codon 4', rows[0]['internal_stop_codon'])
        self.assertIn('also present in the reference', proc.stdout)

    def test_the_same_stop_without_a_control_is_reported_as_a_failure(self):
        proc, rows = self.case(OPAL_CDS, OPAL_CDS).run(control=False)
        self.assertEqual(proc.returncode, 1)
        self.assertEqual(rows[0]['no_internal_stop'], 'FAIL')

    def test_a_new_stop_fails_even_when_the_reference_is_clean(self):
        proc, rows = self.case(OPAL_CDS, CLEAN_CDS).run()
        self.assertEqual(proc.returncode, 1)
        self.assertEqual(rows[0]['no_internal_stop'], 'FAIL')
        self.assertIn('stop codon TGA', proc.stderr)

    def test_broken_reading_frame_fails(self):
        # 59 bp on the stock, 60 bp on the reference: the transfer lost a base.
        proc, rows = self.case(CLEAN_CDS, CLEAN_CDS, stock_kw={'end': 69}).run()
        self.assertEqual(proc.returncode, 1)
        self.assertEqual(rows[0]['length_mod3'], 'FAIL')
        self.assertIn('multiple of three', proc.stderr)

    def test_coordinates_past_the_end_of_the_contig_fail(self):
        proc, rows = self.case(CLEAN_CDS, CLEAN_CDS, stock_kw={'end': 9000}).run()
        self.assertEqual(proc.returncode, 1)
        self.assertEqual(rows[0]['in_contig'], 'FAIL')
        self.assertIn('outside', proc.stderr)

    def test_a_frame_break_the_reference_already_had_is_inherited(self):
        # Both sides 59 bp: the lab annotation was already out of frame, which
        # the transfer did not cause and must not be blamed for.
        proc, rows = self.case(CLEAN_CDS, CLEAN_CDS, end=69).run()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rows[0]['length_mod3'], 'inherited')

    def test_a_cds_lost_in_transfer_fails(self):
        case = self.case(CLEAN_CDS, CLEAN_CDS)
        with open(case.ref_gff, 'a', encoding='utf-8') as f:
            f.write('ref\tTest\tCDS\t11\t70\t.\t+\t0\tID=cds:E2;Parent=gene:E2;Name=E2\n')
        proc, _ = case.run()
        self.assertEqual(proc.returncode, 1)
        self.assertIn('not transferred: E2', proc.stderr)

    def test_start_and_stop_are_informational_by_default(self):
        # A polyprotein cleavage product has neither; that must not fail.
        headless = 'AAG' * 20
        proc, rows = self.case(headless, headless).run()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual((rows[0]['starts_atg'], rows[0]['ends_stop']), ('no', 'no'))

    def test_require_start_stop_turns_them_into_failures(self):
        headless = 'AAG' * 20
        proc, _ = self.case(headless, headless).run('--require-start-stop')
        self.assertEqual(proc.returncode, 1)
        self.assertIn('does not start with ATG', proc.stderr)

    def test_reverse_strand_cds_is_read_on_the_right_strand(self):
        rc = CLEAN_CDS.translate(cla.COMPLEMENT)[::-1]
        proc, rows = self.case(rc, rc, strand='-').run()
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(rows[0]['starts_atg'], 'yes')
        self.assertEqual(rows[0]['ends_stop'], 'yes')

    def test_half_a_control_is_refused(self):
        case = self.case(CLEAN_CDS, CLEAN_CDS)
        proc = subprocess.run(
            [sys.executable, SCRIPT, '--fasta', case.stock_fasta, '--gff', case.stock_gff,
             '--reference-gff', case.ref_gff, '--output-tsv', case.out],
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, 1)
        self.assertIn('must be given together', proc.stderr)

    def test_gff_without_cds_is_an_error(self):
        case = self.case(CLEAN_CDS, CLEAN_CDS)
        with open(case.stock_gff, 'w', encoding='utf-8') as f:
            f.write('##gff-version 3\nstock\tTest\tgene\t11\t70\t.\t+\t.\tID=g;Name=g\n')
        proc, _ = case.run()
        self.assertEqual(proc.returncode, 1)
        self.assertIn('no CDS features', proc.stderr)


class FixtureTests(unittest.TestCase):
    def test_the_committed_stock_fixture_passes_against_its_reference(self):
        data = os.path.join(REPO, 'tests', 'data')
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, 'check.tsv')
            proc = subprocess.run(
                [sys.executable, SCRIPT,
                 '--fasta', os.path.join(data, 'stock.test.fasta'),
                 '--gff', os.path.join(data, 'stock.test.gff3'),
                 '--reference-fasta', os.path.join(data, 'viral_ref.test.fasta'),
                 '--reference-gff', os.path.join(data, 'viral_ref.test.gff3'),
                 '--output-tsv', out], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with open(out, encoding='utf-8') as f:
                lines = f.read().splitlines()
            self.assertEqual(len(lines), 2)  # header plus the one CDS

    def test_bad_flag_exits_2(self):
        proc = subprocess.run([sys.executable, SCRIPT, '--badflag'],
                              capture_output=True, text=True)
        self.assertEqual(proc.returncode, 2)


if __name__ == '__main__':
    unittest.main()
