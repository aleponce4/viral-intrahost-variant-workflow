#!/usr/bin/env python3
"""Build a per-position coordinate map between a stock consensus and a reference.

Input is a PAF alignment carrying a short `cs` tag, produced by
`minimap2 -cx asm10 --cs <reference.fa> <stock.fa>`, so the stock is the PAF
query and the reference is the PAF target.

Output is one row per alignment column, which lets coordinates move in either
direction: a variant called on the stock can be reported at its reference
position, and a reference-coordinate interval such as a primer BED can be
placed on the stock. Positions are 1-based.

  status      stock_pos   ref_pos   meaning
  match       set         set       same base in both
  mismatch    set         set       substitution
  insertion   set         NA        base present only in the stock
  deletion    NA          set       base present only in the reference
"""
import argparse
import os
import re
import sys

__version__ = "1.0.0"

# Short cs operations: `:N` identical run, `*xy` substitution, `+seq` insertion
# into the query, `-seq` deletion from the query. `~` (intron) cannot occur in a
# genome-to-genome asm alignment and is rejected rather than guessed at.
CS_TOKEN = re.compile(r'([:*+\-~])([A-Za-z0-9]+)')


def parse_cs(cs):
    """Yield (op, value) pairs from a short cs string."""
    pos = 0
    while pos < len(cs):
        m = CS_TOKEN.match(cs, pos)
        if not m:
            raise ValueError(f"cannot parse cs tag at offset {pos}: {cs[pos:pos + 20]!r}")
        op, val = m.group(1), m.group(2)
        if op == ':':
            if not val.isdigit():
                raise ValueError(f"cs run length is not a number: {val!r}")
            yield op, int(val)
        elif op == '*':
            if len(val) != 2:
                raise ValueError(f"cs substitution needs two bases, found {val!r}")
            yield op, val
        elif op == '~':
            raise ValueError("cs tag contains an intron operation; this script expects "
                             "a genome-to-genome alignment (minimap2 -x asm*)")
        else:
            yield op, val
        pos = m.end()


def pick_alignment(paf_file):
    """Return the single best PAF line, as a dict.

    A stock consensus should align to its reference as one block. Several lines
    mean a rearrangement, a multi-contig assembly or the wrong reference, none of
    which this table can express, so they stop the run.
    """
    rows = []
    with open(paf_file, 'r', encoding='utf-8') as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 12:
                raise ValueError(f"{paf_file}:{line_no}: expected at least 12 PAF columns, "
                                 f"found {len(fields)}")
            cs = None
            for tag in fields[12:]:
                if tag.startswith('cs:Z:'):
                    cs = tag[5:]
            if cs is None:
                raise ValueError(f"{paf_file}:{line_no}: no cs tag. Run minimap2 with --cs.")
            rows.append({
                'query': fields[0], 'query_start': int(fields[2]), 'query_end': int(fields[3]),
                'strand': fields[4],
                'target': fields[5], 'target_start': int(fields[7]), 'target_end': int(fields[8]),
                'matches': int(fields[9]), 'block': int(fields[10]),
                'cs': cs, 'line_no': line_no,
            })
    if not rows:
        raise ValueError(f"{paf_file}: no alignment. The stock and the reference did not align.")
    if len(rows) > 1:
        where = ', '.join(f"{r['query']}:{r['query_start']}-{r['query_end']}" for r in rows)
        raise ValueError(
            f"{paf_file}: expected one alignment block, found {len(rows)} ({where}). "
            "A stock consensus must align to its reference end to end."
        )
    row = rows[0]
    if row['strand'] != '+':
        raise ValueError(f"{paf_file}: alignment is on the '-' strand. Reverse-complement the "
                         "stock consensus before building the table.")
    return row


def build_rows(aln):
    """Walk the cs tag and yield one dict per alignment column."""
    q = aln['query_start']   # 0-based cursor into the stock
    t = aln['target_start']  # 0-based cursor into the reference
    for op, val in parse_cs(aln['cs']):
        if op == ':':
            for _ in range(val):
                yield {'stock_pos': q + 1, 'ref_pos': t + 1, 'status': 'match',
                       'stock_base': '.', 'ref_base': '.'}
                q += 1
                t += 1
        elif op == '*':
            ref_base, stock_base = val[0].upper(), val[1].upper()
            yield {'stock_pos': q + 1, 'ref_pos': t + 1, 'status': 'mismatch',
                   'stock_base': stock_base, 'ref_base': ref_base}
            q += 1
            t += 1
        elif op == '+':
            for base in val.upper():
                yield {'stock_pos': q + 1, 'ref_pos': 'NA', 'status': 'insertion',
                       'stock_base': base, 'ref_base': 'NA'}
                q += 1
        elif op == '-':
            for base in val.upper():
                yield {'stock_pos': 'NA', 'ref_pos': t + 1, 'status': 'deletion',
                       'stock_base': 'NA', 'ref_base': base}
                t += 1
    if q != aln['query_end']:
        raise ValueError(f"cs tag consumed {q} stock bases but the PAF record ends at "
                         f"{aln['query_end']}")
    if t != aln['target_end']:
        raise ValueError(f"cs tag consumed {t} reference bases but the PAF record ends at "
                         f"{aln['target_end']}")


def summarize(rows, aln):
    counts = {'match': 0, 'mismatch': 0, 'insertion': 0, 'deletion': 0}
    for r in rows:
        counts[r['status']] += 1
    aligned = counts['match'] + counts['mismatch']
    identity = (counts['match'] / aligned) if aligned else 0.0
    return {
        'stock_contig': aln['query'], 'ref_contig': aln['target'],
        'stock_aligned_from': aln['query_start'] + 1, 'stock_aligned_to': aln['query_end'],
        'ref_aligned_from': aln['target_start'] + 1, 'ref_aligned_to': aln['target_end'],
        'identity': f"{identity:.6f}", **counts,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Build a per-position coordinate map between a stock consensus and a reference.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--paf", required=True,
                        help="PAF from `minimap2 -cx asm10 --cs <reference> <stock>`")
    parser.add_argument("--output-tsv", required=True, help="Per-position coordinate map")
    parser.add_argument("--summary-tsv", help="One-line alignment summary")
    args = parser.parse_args()

    if not os.path.exists(args.paf):
        print(f"ERROR: PAF file not found: {args.paf}", file=sys.stderr)
        sys.exit(1)

    try:
        aln = pick_alignment(args.paf)
        rows = list(build_rows(aln))
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    cols = ['stock_contig', 'stock_pos', 'ref_contig', 'ref_pos', 'status', 'stock_base', 'ref_base']
    with open(args.output_tsv, 'w', encoding='utf-8', newline='') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in [
                aln['query'], r['stock_pos'], aln['target'], r['ref_pos'],
                r['status'], r['stock_base'], r['ref_base']]) + '\n')

    summary = summarize(rows, aln)
    if args.summary_tsv:
        with open(args.summary_tsv, 'w', encoding='utf-8', newline='') as f:
            f.write('\t'.join(summary) + '\n')
            f.write('\t'.join(str(summary[k]) for k in summary) + '\n')

    print(f"{aln['query']} vs {aln['target']}: {summary['match']} match, "
          f"{summary['mismatch']} mismatch, {summary['insertion']} inserted, "
          f"{summary['deletion']} deleted, identity {summary['identity']}")


if __name__ == "__main__":
    main()
