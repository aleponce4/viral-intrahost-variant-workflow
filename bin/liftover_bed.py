#!/usr/bin/env python3
"""Move a BED file between a lab reference and a stock's own coordinates.

Amplicon work needs this. A primer scheme is designed once, against the lab
reference, but `ivar trim` has to be given primer positions in the coordinates
of whatever the reads were aligned to. Map reads to a per-stock reference and
those coordinates no longer agree: one indel between the stock and the
reference shifts every primer after it.

Input is the table written by `build_liftover_table.py`, which holds one row per
alignment column, and a BED in the source coordinates. Output is the same BED in
the target coordinates, with columns 4 onward carried through untouched.

BED is half-open and 0-based; the liftover table is 1-based. This script does the
conversion, so give it ordinary BED and you get ordinary BED back.

An interval is dropped, not guessed at, when none of its bases exist in the
target. That matters for primers: an amplicon whose primer site the stock has
deleted cannot be trimmed, and silently moving it somewhere plausible would
corrupt every read in that amplicon. Dropped and shortened intervals are both
reported, and --strict turns any of either into a non-zero exit.
"""
import argparse
import os
import sys

__version__ = "1.0.0"


def read_liftover(path, direction):
    """Return {source 1-based position: target 1-based position} and the contig names."""
    src_col, dst_col = ('ref_pos', 'stock_pos') if direction == 'ref-to-stock' else ('stock_pos', 'ref_pos')
    src_name, dst_name = ('ref_contig', 'stock_contig') if direction == 'ref-to-stock' else ('stock_contig', 'ref_contig')
    mapping, source_contig, target_contig = {}, None, None
    with open(path, 'r', encoding='utf-8') as f:
        header = f.readline().rstrip('\r\n').split('\t')
        missing = [c for c in (src_col, dst_col, src_name, dst_name) if c not in header]
        if missing:
            raise ValueError(f"{path}: missing column(s): {', '.join(missing)}")
        idx = {name: header.index(name) for name in header}
        for line_no, line in enumerate(f, start=2):
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            fields = line.split('\t')
            source_contig = source_contig or fields[idx[src_name]]
            target_contig = target_contig or fields[idx[dst_name]]
            s, d = fields[idx[src_col]], fields[idx[dst_col]]
            if s == 'NA' or d == 'NA':
                continue
            try:
                mapping[int(s)] = int(d)
            except ValueError:
                raise ValueError(f"{path}:{line_no}: position is not a number") from None
    if not mapping:
        raise ValueError(f"{path}: no position maps in the {direction} direction")
    return mapping, source_contig, target_contig


def lift_interval(start0, end0, mapping):
    """Map a half-open 0-based interval. Returns (start0, end0, n_source_bases_kept)."""
    # Bases of the interval that exist in the target, as 1-based source positions.
    present = [mapping[p] for p in range(start0 + 1, end0 + 1) if p in mapping]
    if not present:
        return None, None, 0
    return min(present) - 1, max(present), len(present)


def main():
    parser = argparse.ArgumentParser(
        description="Move a BED file between a lab reference and a stock's own coordinates.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--liftover-tsv", required=True,
                        help="Table from build_liftover_table.py")
    parser.add_argument("--input-bed", required=True, help="BED in the source coordinates")
    parser.add_argument("--output-bed", required=True, help="BED in the target coordinates")
    parser.add_argument("--direction", choices=['ref-to-stock', 'stock-to-ref'],
                        default='ref-to-stock', help="Default: ref-to-stock")
    parser.add_argument("--report-tsv", help="One row per input interval, with what happened to it")
    parser.add_argument("--strict", action="store_true",
                        help="Exit non-zero if any interval was dropped or changed length")
    args = parser.parse_args()

    for path in (args.liftover_tsv, args.input_bed):
        if not os.path.exists(path):
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    try:
        mapping, source_contig, target_contig = read_liftover(args.liftover_tsv, args.direction)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    rows, dropped, shortened, kept = [], 0, 0, 0
    with open(args.input_bed, 'r', encoding='utf-8') as fin, \
            open(args.output_bed, 'w', encoding='utf-8', newline='') as fout:
        for line_no, line in enumerate(fin, start=1):
            stripped = line.rstrip('\r\n')
            if not stripped.strip() or stripped.startswith(('#', 'track', 'browser')):
                continue
            fields = stripped.split('\t')
            if len(fields) < 3:
                print(f"ERROR: {args.input_bed}:{line_no}: BED needs at least 3 columns, "
                      f"found {len(fields)}", file=sys.stderr)
                sys.exit(1)
            try:
                start0, end0 = int(fields[1]), int(fields[2])
            except ValueError:
                print(f"ERROR: {args.input_bed}:{line_no}: start and end must be numbers",
                      file=sys.stderr)
                sys.exit(1)
            if end0 <= start0:
                print(f"ERROR: {args.input_bed}:{line_no}: end must be greater than start",
                      file=sys.stderr)
                sys.exit(1)

            name = fields[3] if len(fields) > 3 else f"interval_{line_no}"
            new_start, new_end, n_kept = lift_interval(start0, end0, mapping)
            width_in = end0 - start0

            if new_start is None:
                dropped += 1
                status = 'dropped'
                rows.append((name, fields[0], start0, end0, width_in, 'NA', 'NA', 'NA', 0, status))
                continue

            width_out = new_end - new_start
            if n_kept < width_in or width_out != width_in:
                shortened += 1
                status = 'shortened' if n_kept < width_in else 'resized'
            else:
                status = 'exact'
            kept += 1
            out = [target_contig, str(new_start), str(new_end)] + fields[3:]
            fout.write('\t'.join(out) + '\n')
            rows.append((name, fields[0], start0, end0, width_in,
                         target_contig, new_start, new_end, width_out, status))

    if args.report_tsv:
        with open(args.report_tsv, 'w', encoding='utf-8', newline='') as f:
            f.write('\t'.join(['name', 'source_contig', 'source_start', 'source_end',
                               'source_width', 'target_contig', 'target_start', 'target_end',
                               'target_width', 'status']) + '\n')
            for r in rows:
                f.write('\t'.join(str(x) for x in r) + '\n')

    print(f"{kept} interval(s) placed on {target_contig}, {dropped} dropped, "
          f"{shortened} changed width")
    for r in rows:
        if r[-1] == 'dropped':
            print(f"  dropped: {r[0]} at {r[1]}:{r[2]}-{r[3]} has no bases in the target",
                  file=sys.stderr)
        elif r[-1] in ('shortened', 'resized'):
            print(f"  {r[-1]}: {r[0]} {r[4]} bp -> {r[8]} bp", file=sys.stderr)

    if args.strict and (dropped or shortened):
        print(f"FAIL: {dropped} dropped and {shortened} resized under --strict", file=sys.stderr)
        sys.exit(1)
    if kept == 0:
        print("ERROR: no interval could be placed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
