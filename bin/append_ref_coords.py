#!/usr/bin/env python3
"""Add lab-reference coordinates to variants that were called on a stock.

A stock is analysed against its own de novo reference, so every call carries a
stock position. That position is not comparable across stocks, and it is not the
one the lab's other datasets use. This script looks each call up in the table
written by `build_liftover_table.py` and records where it falls on the lab
reference.

It reads a VCF (plain or gzip) or a tab-separated table and writes the same
records with these additions. Positions are 1-based.

  REF_CONTIG   lab reference contig the position maps to
  REF_POS      lab reference position of the first REF base
  REF_END      lab reference position of the last REF base, only when REF spans
               more than one base
  REF_STATUS   match | mismatch | insertion | unmapped, for the first REF base

For a table the same values go in columns named ref_contig, ref_pos, ref_end and
ref_status, with NA for a missing value. `insertion` means the stock has a base
the reference lacks, so there is no reference position. `unmapped` means the
position is outside the aligned block, for example soft-clipped ends.

The call itself is not changed. A `mismatch` row is a position where the stock
differs from the lab reference, and a variant called there is still a variant
within the stock.
"""
import argparse
import gzip
import os
import sys

__version__ = "1.0.0"

VCF_HEADER_LINES = [
    '##INFO=<ID=REF_CONTIG,Number=1,Type=String,Description="Lab reference contig this position maps to">',
    '##INFO=<ID=REF_POS,Number=1,Type=Integer,Description="Lab reference position of the first REF base (1-based)">',
    '##INFO=<ID=REF_END,Number=1,Type=Integer,Description="Lab reference position of the last REF base, when REF spans more than one base">',
    '##INFO=<ID=REF_STATUS,Number=1,Type=String,Description="Alignment status of the first REF base against the lab reference: match, mismatch, insertion or unmapped">',
]


def open_text(path, mode='rt'):
    if str(path).endswith('.gz'):
        return gzip.open(path, mode, encoding='utf-8')
    return open(path, mode, encoding='utf-8', newline='' if 'w' in mode else None)


def read_liftover(path):
    """Return {stock_contig: {stock_pos: (ref_contig, ref_pos or None, status)}}."""
    needed = ('stock_contig', 'stock_pos', 'ref_contig', 'ref_pos', 'status')
    table = {}
    with open(path, 'r', encoding='utf-8') as f:
        header = f.readline().rstrip('\r\n').split('\t')
        missing = [c for c in needed if c not in header]
        if missing:
            raise ValueError(f"{path}: missing column(s): {', '.join(missing)}")
        idx = {name: header.index(name) for name in needed}
        for line_no, line in enumerate(f, start=2):
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            fields = line.split('\t')
            stock_pos = fields[idx['stock_pos']]
            if stock_pos == 'NA':
                continue  # a deletion row: the reference has a base the stock lacks
            ref_pos = fields[idx['ref_pos']]
            try:
                stock_pos = int(stock_pos)
                ref_pos = None if ref_pos == 'NA' else int(ref_pos)
            except ValueError:
                raise ValueError(f"{path}:{line_no}: position is not a number") from None
            table.setdefault(fields[idx['stock_contig']], {})[stock_pos] = (
                fields[idx['ref_contig']], ref_pos, fields[idx['status']])
    if not table:
        raise ValueError(f"{path}: no positions found")
    return table


def lift(contig_map, pos, ref_len):
    """Return (ref_contig, ref_pos, ref_end, status) for a call at `pos`.

    Values are None where there is nothing to report.
    """
    first = contig_map.get(pos)
    if first is None:
        return None, None, None, 'unmapped'
    ref_contig, ref_pos, status = first
    ref_end = None
    if ref_len > 1:
        last = contig_map.get(pos + ref_len - 1)
        ref_end = last[1] if last else None
    return ref_contig, ref_pos, ref_end, status


def contig_map_for(table, contig):
    if contig not in table:
        raise ValueError(
            f"contig {contig!r} is not in the liftover table. The table lists: "
            f"{', '.join(sorted(table))}. Use the table built for this stock.")
    return table[contig]


def annotate_vcf(table, src, dst, counts):
    with open_text(src) as fin, open_text(dst, 'wt') as fout:
        header_written = False
        for line_no, line in enumerate(fin, start=1):
            line = line.rstrip('\r\n')
            if line.startswith('##'):
                if line.startswith('##INFO=<ID=REF_POS,'):
                    raise ValueError(f"{src}: already carries REF_POS; it was annotated before")
                fout.write(line + '\n')
                continue
            if line.startswith('#'):
                for extra in VCF_HEADER_LINES:
                    fout.write(extra + '\n')
                fout.write(line + '\n')
                header_written = True
                continue
            if not line.strip():
                continue
            if not header_written:
                raise ValueError(f"{src}:{line_no}: record before the #CHROM line")
            parts = line.split('\t')
            if len(parts) < 8:
                raise ValueError(f"{src}:{line_no}: a VCF record needs 8 columns, found {len(parts)}")
            try:
                pos = int(parts[1])
            except ValueError:
                raise ValueError(f"{src}:{line_no}: POS is not a number") from None
            cmap = contig_map_for(table, parts[0])
            ref_contig, ref_pos, ref_end, status = lift(cmap, pos, len(parts[3]))
            extra = []
            if ref_contig is not None:
                extra.append(f"REF_CONTIG={ref_contig}")
            if ref_pos is not None:
                extra.append(f"REF_POS={ref_pos}")
            if ref_end is not None:
                extra.append(f"REF_END={ref_end}")
            extra.append(f"REF_STATUS={status}")
            parts[7] = ';'.join(extra) if parts[7] in ('', '.') else parts[7] + ';' + ';'.join(extra)
            fout.write('\t'.join(parts) + '\n')
            counts[status] = counts.get(status, 0) + 1


def annotate_table(table, src, dst, chrom_col, pos_col, ref_col, counts):
    with open_text(src) as fin, open_text(dst, 'wt') as fout:
        header = fin.readline().rstrip('\r\n').split('\t')
        for col in (chrom_col, pos_col):
            if col not in header:
                raise ValueError(f"{src}: no column named {col!r}; found {', '.join(header)}")
        added = ['ref_contig', 'ref_pos', 'ref_end', 'ref_status']
        clash = [c for c in added if c in header]
        if clash:
            raise ValueError(f"{src}: already has column(s) {', '.join(clash)}; it was annotated before")
        ci, pi = header.index(chrom_col), header.index(pos_col)
        ri = header.index(ref_col) if ref_col in header else None
        fout.write('\t'.join(header + added) + '\n')
        for line_no, line in enumerate(fin, start=2):
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) != len(header):
                raise ValueError(f"{src}:{line_no}: expected {len(header)} columns, found {len(fields)}")
            try:
                pos = int(fields[pi])
            except ValueError:
                raise ValueError(f"{src}:{line_no}: {pos_col} is not a number") from None
            ref_len = len(fields[ri]) if ri is not None else 1
            cmap = contig_map_for(table, fields[ci])
            ref_contig, ref_pos, ref_end, status = lift(cmap, pos, ref_len)
            fout.write('\t'.join(fields + [
                'NA' if ref_contig is None else ref_contig,
                'NA' if ref_pos is None else str(ref_pos),
                'NA' if ref_end is None else str(ref_end),
                status]) + '\n')
            counts[status] = counts.get(status, 0) + 1


def detect_format(path):
    name = str(path)
    if name.endswith(('.vcf', '.vcf.gz')):
        return 'vcf'
    if name.endswith(('.tsv', '.tsv.gz', '.txt')):
        return 'tsv'
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Add lab-reference coordinates to variants called on a stock.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--liftover-tsv", required=True, help="Table from build_liftover_table.py")
    parser.add_argument("--input", required=True, help="VCF (plain or gzip) or TSV of variants on the stock")
    parser.add_argument("--output", required=True, help="Annotated file, written uncompressed")
    parser.add_argument("--format", choices=['auto', 'vcf', 'tsv'], default='auto',
                        help="Input format. Default: from the file name")
    parser.add_argument("--chrom-col", default="chrom", help="TSV column holding the contig. Default: chrom")
    parser.add_argument("--pos-col", default="pos", help="TSV column holding the position. Default: pos")
    parser.add_argument("--ref-col", default="ref",
                        help="TSV column holding the REF allele, used to find REF_END. Default: ref")
    args = parser.parse_args()

    for path in (args.liftover_tsv, args.input):
        if not os.path.exists(path):
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    fmt = args.format if args.format != 'auto' else detect_format(args.input)
    if fmt is None:
        print(f"ERROR: cannot tell the format of {args.input}; pass --format vcf or --format tsv",
              file=sys.stderr)
        sys.exit(1)

    counts = {}
    try:
        table = read_liftover(args.liftover_tsv)
        if fmt == 'vcf':
            annotate_vcf(table, args.input, args.output, counts)
        else:
            annotate_table(table, args.input, args.output,
                           args.chrom_col, args.pos_col, args.ref_col, counts)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    total = sum(counts.values())
    detail = ', '.join(f"{n} {s}" for s, n in sorted(counts.items())) or 'none'
    print(f"Annotated {total} record(s): {detail}")


if __name__ == "__main__":
    main()
