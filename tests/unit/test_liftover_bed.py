"""Unit tests for bin/liftover_bed.py.

Run from the repository root:  python -m unittest discover -s tests/unit -v

A primer placed one base off trims one base too many or too few from every read
in its amplicon, so these tests state the expected coordinates outright. BED is
half-open and 0-based throughout.
"""
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SCRIPT = os.path.join(REPO, 'bin', 'liftover_bed.py')
sys.path.insert(0, os.path.join(REPO, 'bin'))

import liftover_bed as lb  # noqa: E402

HEADER = ['stock_contig', 'stock_pos', 'ref_contig', 'ref_pos', 'status', 'stock_base', 'ref_base']


def table(rows, stock='stockX', ref='refX'):
    """rows: list of (stock_pos, ref_pos, status) with None meaning NA."""
    out = ['\t'.join(HEADER)]
    for s, r, status in rows:
        out.append('\t'.join([stock, 'NA' if s is None else str(s), ref,
                              'NA' if r is None else str(r), status, '.', '.']))
    return '\n'.join(out) + '\n'


def identity_rows(n, offset=0):
    """n positions where stock i maps to ref i+offset."""
    return [(i, i + offset, 'match') for i in range(1, n + 1)]


class LiftIntervalTests(unittest.TestCase):
    def test_one_to_one_interval_keeps_its_coordinates(self):
        m = {p: p for p in range(1, 101)}
        self.assertEqual(lb.lift_interval(10, 20, m), (10, 20, 10))

    def test_a_constant_shift_moves_the_interval(self):
        # every reference base sits 5 further along the stock
        m = {p: p + 5 for p in range(1, 101)}
        self.assertEqual(lb.lift_interval(10, 20, m), (15, 25, 10))

    def test_an_interval_over_a_deletion_comes_back_shorter(self):
        # reference 13,14,15 are absent from the stock
        m = {p: p for p in range(1, 13)}
        m.update({p: p - 3 for p in range(16, 101)})
        start, end, kept = lb.lift_interval(10, 20, m)
        self.assertEqual(kept, 7)              # 10 reference bases, 3 of them gone
        self.assertEqual((start, end), (10, 17))

    def test_an_interval_entirely_inside_a_deletion_is_dropped(self):
        m = {p: p for p in range(1, 11)}       # reference 11+ absent
        self.assertEqual(lb.lift_interval(20, 30, m), (None, None, 0))

    def test_only_the_overlapping_part_is_used(self):
        m = {p: p for p in range(1, 16)}       # alignment stops at 15
        start, end, kept = lb.lift_interval(10, 20, m)
        self.assertEqual((start, end, kept), (10, 15, 5))


class CliTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_cli(self, table_text, bed_text, *extra):
        d = self.tmp.name
        tsv, bed, out, rep = (os.path.join(d, n) for n in
                              ('map.tsv', 'in.bed', 'out.bed', 'report.tsv'))
        for path, text in ((tsv, table_text), (bed, bed_text)):
            with open(path, 'w', encoding='utf-8', newline='') as f:
                f.write(text)
        proc = subprocess.run(
            [sys.executable, SCRIPT, '--liftover-tsv', tsv, '--input-bed', bed,
             '--output-bed', out, '--report-tsv', rep, *extra],
            capture_output=True, text=True)
        lines = []
        if os.path.exists(out):
            with open(out, encoding='utf-8') as f:
                lines = [l.split('\t') for l in f.read().splitlines() if l.strip()]
        report = []
        if os.path.exists(rep):
            with open(rep, encoding='utf-8') as f:
                rows = f.read().splitlines()
                keys = rows[0].split('\t')
                report = [dict(zip(keys, r.split('\t'))) for r in rows[1:]]
        return proc, lines, report

    def test_primer_bed_moves_onto_the_stock_and_keeps_its_columns(self):
        # The stock carries a 3 bp deletion of reference 13-15, so every primer
        # after it sits 3 bases earlier.
        rows = [(p, p, 'match') for p in range(1, 13)]
        rows += [(None, p, 'deletion') for p in (13, 14, 15)]
        rows += [(p - 3, p, 'match') for p in range(16, 101)]
        bed = ('refX\t5\t10\tprimer_1_LEFT\t1\t+\n'
               'refX\t40\t50\tprimer_1_RIGHT\t1\t-\n')
        proc, out, report = self.run_cli(table(rows), bed)
        self.assertEqual(proc.returncode, 0, proc.stderr)

        # Before the deletion: unchanged. After it: shifted by exactly 3.
        self.assertEqual(out[0][:4], ['stockX', '5', '10', 'primer_1_LEFT'])
        self.assertEqual(out[1][:4], ['stockX', '37', '47', 'primer_1_RIGHT'])
        # ivar trim reads the score and strand columns, so they must survive.
        self.assertEqual(out[0][4:], ['1', '+'])
        self.assertEqual(out[1][4:], ['1', '-'])
        self.assertTrue(all(r['status'] == 'exact' for r in report))

    def test_a_primer_over_a_deletion_is_reported_as_shortened(self):
        rows = [(p, p, 'match') for p in range(1, 13)]
        rows += [(None, p, 'deletion') for p in (13, 14, 15)]
        rows += [(p - 3, p, 'match') for p in range(16, 101)]
        proc, out, report = self.run_cli(table(rows), 'refX\t10\t20\tspans_del\t1\t+\n')
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(report[0]['status'], 'shortened')
        self.assertEqual((report[0]['source_width'], report[0]['target_width']), ('10', '7'))
        self.assertIn('shortened', proc.stderr)

    def test_a_primer_the_stock_has_deleted_is_dropped_not_guessed(self):
        # Reference 20-30 is gone from the stock: that amplicon cannot be trimmed.
        rows = [(p, p, 'match') for p in range(1, 20)]
        rows += [(None, p, 'deletion') for p in range(20, 31)]
        rows += [(p - 11, p, 'match') for p in range(31, 101)]
        bed = 'refX\t19\t30\tlost_primer\t1\t+\nrefX\t40\t50\tkept_primer\t1\t-\n'
        proc, out, report = self.run_cli(table(rows), bed)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0][3], 'kept_primer')
        self.assertEqual([r['status'] for r in report], ['dropped', 'exact'])
        self.assertIn('dropped', proc.stderr)

    def test_strict_fails_when_anything_was_dropped(self):
        rows = [(p, p, 'match') for p in range(1, 20)]
        rows += [(None, p, 'deletion') for p in range(20, 31)]
        rows += [(p - 11, p, 'match') for p in range(31, 101)]
        bed = 'refX\t19\t30\tlost\t1\t+\nrefX\t40\t50\tkept\t1\t-\n'
        proc, _, _ = self.run_cli(table(rows), bed, '--strict')
        self.assertEqual(proc.returncode, 1)
        self.assertIn('--strict', proc.stderr)

    def test_the_other_direction(self):
        # stock-to-ref on a stock that gained 2 bases after position 10
        rows = [(p, p, 'match') for p in range(1, 11)]
        rows += [(11, None, 'insertion'), (12, None, 'insertion')]
        rows += [(p, p - 2, 'match') for p in range(13, 101)]
        proc, out, _ = self.run_cli(table(rows), 'stockX\t20\t30\treg\t0\t+\n',
                                    '--direction', 'stock-to-ref')
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out[0][:3], ['refX', '18', '28'])

    def test_inserted_stock_bases_have_no_reference_position(self):
        rows = [(p, p, 'match') for p in range(1, 11)]
        rows += [(11, None, 'insertion'), (12, None, 'insertion')]
        rows += [(p, p - 2, 'match') for p in range(13, 101)]
        # the interval covers only the two inserted bases
        proc, out, report = self.run_cli(table(rows), 'stockX\t10\t12\tins_only\t0\t+\n',
                                         '--direction', 'stock-to-ref')
        self.assertEqual(report[0]['status'], 'dropped')
        self.assertEqual(len(out), 0)
        self.assertEqual(proc.returncode, 1)  # nothing placed at all

    def test_comments_and_track_lines_are_skipped(self):
        proc, out, _ = self.run_cli(table(identity_rows(100)),
                                    '# a comment\ntrack name=x\nrefX\t5\t10\tp\t1\t+\n')
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(len(out), 1)

    def test_a_three_column_bed_works(self):
        proc, out, _ = self.run_cli(table(identity_rows(100)), 'refX\t5\t10\n')
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(out[0], ['refX'.replace('refX', 'stockX'), '5', '10'])

    def test_a_malformed_bed_is_an_error(self):
        proc, _, _ = self.run_cli(table(identity_rows(100)), 'refX\t5\n')
        self.assertEqual(proc.returncode, 1)
        self.assertIn('at least 3 columns', proc.stderr)

    def test_a_backwards_interval_is_an_error(self):
        proc, _, _ = self.run_cli(table(identity_rows(100)), 'refX\t20\t10\tbad\n')
        self.assertEqual(proc.returncode, 1)
        self.assertIn('greater than start', proc.stderr)

    def test_bad_flag_exits_2(self):
        proc = subprocess.run([sys.executable, SCRIPT, '--badflag'],
                              capture_output=True, text=True)
        self.assertEqual(proc.returncode, 2)


class FixtureTests(unittest.TestCase):
    def test_the_committed_primer_bed_moves_onto_the_stock_fixture(self):
        """The real chain: build the map from the fixture PAF, then move the primers."""
        data = os.path.join(REPO, 'tests', 'data')
        with tempfile.TemporaryDirectory() as d:
            tsv = os.path.join(d, 'map.tsv')
            proc = subprocess.run(
                [sys.executable, os.path.join(REPO, 'bin', 'build_liftover_table.py'),
                 '--paf', os.path.join(data, 'stock.test.paf'), '--output-tsv', tsv],
                capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)

            out = os.path.join(d, 'out.bed')
            proc = subprocess.run(
                [sys.executable, SCRIPT, '--liftover-tsv', tsv,
                 '--input-bed', os.path.join(data, 'primers.test.bed'),
                 '--output-bed', out, '--direction', 'ref-to-stock', '--strict'],
                capture_output=True, text=True)
            # Every primer in the fixture sits well before the stock's first
            # difference, so all four must come across unchanged and --strict
            # must be satisfied.
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with open(out, encoding='utf-8') as f:
                rows = [l.split('\t') for l in f.read().splitlines() if l.strip()]
            self.assertEqual(len(rows), 4)
            self.assertTrue(all(r[0] == 'stock-test' for r in rows))
            self.assertEqual([r[1:4] for r in rows],
                             [['10', '30', 'primer_1_F'], ['120', '140', 'primer_1_R'],
                              ['200', '220', 'primer_2_F'], ['280', '300', 'primer_2_R']])


if __name__ == '__main__':
    unittest.main()
