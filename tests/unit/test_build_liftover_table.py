"""Unit tests for bin/build_liftover_table.py.

Run from the repository root:  python -m unittest discover -s tests/unit -v

The coordinate map is what lets a variant called on a stock be reported at its
reference position, so an off-by-one here would silently move every call. Each
test states the expected position explicitly rather than deriving it.
"""
import os
import subprocess
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
SCRIPT = os.path.join(REPO, 'bin', 'build_liftover_table.py')
sys.path.insert(0, os.path.join(REPO, 'bin'))

import build_liftover_table as blt  # noqa: E402


def paf(cs, qlen=10, qstart=0, qend=10, tlen=10, tstart=0, tend=10,
        query='stock', target='ref', strand='+'):
    cols = [query, qlen, qstart, qend, strand, target, tlen, tstart, tend,
            qend - qstart, qend - qstart, 60, f'cs:Z:{cs}']
    return '\t'.join(str(c) for c in cols) + '\n'


class ParseCsTests(unittest.TestCase):
    def test_each_operation(self):
        ops = list(blt.parse_cs(':5*ac+gg-tt:3'))
        self.assertEqual(ops, [(':', 5), ('*', 'ac'), ('+', 'gg'), ('-', 'tt'), (':', 3)])

    def test_intron_operation_is_rejected(self):
        # A spliced alignment would silently skip reference bases.
        with self.assertRaisesRegex(ValueError, 'intron'):
            list(blt.parse_cs(':5~gt100ag:5'))

    def test_junk_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'cannot parse'):
            list(blt.parse_cs(':5?zz'))


class BuildRowsTests(unittest.TestCase):
    def rows(self, cs, **kw):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'a.paf')
            with open(p, 'w', encoding='utf-8') as f:
                f.write(paf(cs, **kw))
            aln = blt.pick_alignment(p)
            return list(blt.build_rows(aln))

    def test_all_match_is_one_to_one_and_one_based(self):
        rows = self.rows(':3', qlen=3, qend=3, tlen=3, tend=3)
        self.assertEqual([(r['stock_pos'], r['ref_pos']) for r in rows],
                         [(1, 1), (2, 2), (3, 3)])
        self.assertTrue(all(r['status'] == 'match' for r in rows))

    def test_substitution_records_both_bases(self):
        rows = self.rows(':1*ac:1', qlen=3, qend=3, tlen=3, tend=3)
        mid = rows[1]
        self.assertEqual((mid['stock_pos'], mid['ref_pos'], mid['status']), (2, 2, 'mismatch'))
        self.assertEqual((mid['ref_base'], mid['stock_base']), ('A', 'C'))

    def test_insertion_has_no_reference_position_and_shifts_later_bases(self):
        # stock gains 2 bases: stock 2,3 are new; stock 4 lines up with ref 2.
        rows = self.rows(':1+gg:1', qlen=4, qend=4, tlen=2, tend=2)
        self.assertEqual([(r['stock_pos'], r['ref_pos'], r['status']) for r in rows],
                         [(1, 1, 'match'), (2, 'NA', 'insertion'),
                          (3, 'NA', 'insertion'), (4, 2, 'match')])

    def test_deletion_has_no_stock_position_and_shifts_later_bases(self):
        # ref 2,3 are absent from the stock; stock 2 lines up with ref 4.
        rows = self.rows(':1-gg:1', qlen=2, qend=2, tlen=4, tend=4)
        self.assertEqual([(r['stock_pos'], r['ref_pos'], r['status']) for r in rows],
                         [(1, 1, 'match'), ('NA', 2, 'deletion'),
                          ('NA', 3, 'deletion'), (2, 4, 'match')])

    def test_alignment_not_starting_at_zero_keeps_its_offset(self):
        rows = self.rows(':2', qlen=10, qstart=4, qend=6, tlen=10, tstart=7, tend=9)
        self.assertEqual([(r['stock_pos'], r['ref_pos']) for r in rows], [(5, 8), (6, 9)])

    def test_cs_disagreeing_with_the_paf_record_is_an_error(self):
        with self.assertRaisesRegex(ValueError, 'consumed'):
            self.rows(':2', qlen=10, qend=10, tlen=10, tend=10)


class PickAlignmentTests(unittest.TestCase):
    def write(self, text):
        d = tempfile.mkdtemp()
        self.addCleanup(lambda: None)
        p = os.path.join(d, 'a.paf')
        with open(p, 'w', encoding='utf-8') as f:
            f.write(text)
        return p

    def test_missing_cs_tag_is_an_error(self):
        line = '\t'.join(['s', '10', '0', '10', '+', 'r', '10', '0', '10', '10', '10', '60']) + '\n'
        with self.assertRaisesRegex(ValueError, 'no cs tag'):
            blt.pick_alignment(self.write(line))

    def test_no_alignment_is_an_error(self):
        with self.assertRaisesRegex(ValueError, 'no alignment'):
            blt.pick_alignment(self.write(''))

    def test_several_blocks_are_an_error(self):
        # Two blocks mean a rearrangement or the wrong reference; the table
        # cannot express either, so it must not silently use the first one.
        with self.assertRaisesRegex(ValueError, 'one alignment block'):
            blt.pick_alignment(self.write(paf(':10') + paf(':10', qstart=10, qend=20)))

    def test_reverse_strand_is_an_error(self):
        with self.assertRaisesRegex(ValueError, "'-' strand"):
            blt.pick_alignment(self.write(paf(':10', strand='-')))


class FixtureTests(unittest.TestCase):
    """The committed fixture: an SNV, a 2 bp insertion and a 3 bp deletion."""

    def test_cli_on_the_test_fixture(self):
        fixture = os.path.join(REPO, 'tests', 'data', 'stock.test.paf')
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, 'map.tsv')
            summary = os.path.join(d, 'summary.tsv')
            proc = subprocess.run(
                [sys.executable, SCRIPT, '--paf', fixture, '--output-tsv', out,
                 '--summary-tsv', summary], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with open(out, encoding='utf-8') as f:
                lines = f.read().splitlines()
            header, rows = lines[0].split('\t'), [l.split('\t') for l in lines[1:]]

            self.assertEqual(header, ['stock_contig', 'stock_pos', 'ref_contig', 'ref_pos',
                                      'status', 'stock_base', 'ref_base'])
            # Every reference base and every stock base appears exactly once.
            self.assertEqual(sum(1 for r in rows if r[3] != 'NA'), 11444)
            self.assertEqual(sum(1 for r in rows if r[1] != 'NA'), 11443)

            by_status = {}
            for r in rows:
                by_status.setdefault(r[4], []).append(r)
            self.assertEqual(len(by_status['mismatch']), 1)
            self.assertEqual(len(by_status['insertion']), 2)
            self.assertEqual(len(by_status['deletion']), 3)
            self.assertEqual(by_status['mismatch'][0][1], '9000')
            self.assertEqual(by_status['mismatch'][0][3], '9000')

            with open(summary, encoding='utf-8') as f:
                keys, values = [l.split('\t') for l in f.read().splitlines()]
            s = dict(zip(keys, values))
            self.assertEqual(s['stock_contig'], 'stock-test')
            self.assertEqual(s['ref_contig'], 'KP282671.1')
            self.assertEqual((s['match'], s['mismatch'], s['insertion'], s['deletion']),
                             ('11440', '1', '2', '3'))

    def test_missing_paf_exits_nonzero(self):
        proc = subprocess.run(
            [sys.executable, SCRIPT, '--paf', 'does_not_exist.paf', '--output-tsv', 'x.tsv'],
            capture_output=True, text=True)
        self.assertEqual(proc.returncode, 1)
        self.assertIn('ERROR', proc.stderr)

    def test_bad_flag_exits_2(self):
        proc = subprocess.run([sys.executable, SCRIPT, '--badflag'],
                              capture_output=True, text=True)
        self.assertEqual(proc.returncode, 2)


if __name__ == '__main__':
    unittest.main()
