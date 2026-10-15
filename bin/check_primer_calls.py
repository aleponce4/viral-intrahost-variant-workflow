#!/usr/bin/env python3
"""Flag variant calls that sit under a primer.

In an amplicon run, the bases under a primer come from the primer's own sequence as much as
from the template, so a mismatch between a primer and the sample can surface as a low-frequency
variant. `ivar trim` removes primer sequence from reads that carry it, which is why those
positions usually survive on fewer reads than the rest of the genome. A call made on what is
left deserves a second look before it is believed.

On one real amplicon library, 12 of 64 calls sat under a primer. The four at 1% or more were all
at positions covered far below the genome's mean depth, and three of them clustered within 11
bases of one primer junction.

This does not decide anything. It marks the calls a reader should check, and gives the two
numbers that decide it: which primer covers the call, and how the depth there compares with the
rest of the genome.

Input:
  --vcf          A VCF, plain or gzip. Depth comes from INFO/DP, VAF from INFO/DP4 or INFO/AF.
  --primer-bed   BED6 of the primer scheme, in the coordinates the calls were made against.
  --depth        Optional samtools depth table, used for the genome's median depth.

Output: one row per call, tab-separated, and a summary on stdout.
  under_primer   yes  the position is inside a primer interval
                 no   it is not
"""
import argparse
import gzip
import os
import statistics
import sys

__version__ = "1.0.0"

COLUMNS = ["chrom", "pos", "ref", "alt", "vaf", "depth", "under_primer", "primers",
           "depth_vs_median"]


def opener(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def read_primers(path):
    """Return {contig: [(start_1based, end_1based, name), ...]}."""
    primers = {}
    with opener(path) as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_no}: a primer BED needs at least 4 columns "
                                 "(contig, start, end, name)")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                raise ValueError(f"{path}:{line_no}: start and end must be whole numbers") from None
            if end <= start:
                raise ValueError(f"{path}:{line_no}: end ({end}) must be greater than start ({start})")
            # BED is 0-based half-open. A call's POS is 1-based.
            primers.setdefault(fields[0], []).append((start + 1, end, fields[3]))
    if not primers:
        raise ValueError(f"{path}: no primer intervals found")
    return primers


def read_depth(path):
    with opener(path) as f:
        values = []
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 3:
                try:
                    values.append(int(float(fields[2])))
                except ValueError:
                    continue
    return values


def read_calls(path):
    with opener(path) as f:
        for line_no, line in enumerate(f, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_no}: a VCF record needs 8 columns, found {len(fields)}")
            info = dict(kv.split("=", 1) for kv in fields[7].split(";") if "=" in kv)
            try:
                pos = int(fields[1])
                dp = int(info["DP"]) if "DP" in info else 0
                dp4 = [int(x) for x in info["DP4"].split(",")] if "DP4" in info else None
                if dp4 and len(dp4) == 4 and sum(dp4) > 0:
                    vaf = (dp4[2] + dp4[3]) / sum(dp4)
                else:
                    vaf = float(info.get("AF", "0").split(",")[0])
            except ValueError:
                raise ValueError(f"{path}:{line_no}: could not read POS, DP, DP4 or AF") from None
            yield fields[0], pos, fields[3], fields[4], vaf, dp


def main():
    parser = argparse.ArgumentParser(description="Flag variant calls that sit under a primer.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--primer-bed", required=True)
    parser.add_argument("--depth", help="samtools depth table, for the genome's median depth")
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-vaf", type=float, default=0.01,
                        help="Calls under a primer at or above this VAF are listed. Default: 0.01")
    args = parser.parse_args()

    for path in (args.vcf, args.primer_bed, args.depth):
        if path and not os.path.exists(path):
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    try:
        primers = read_primers(args.primer_bed)
        depths = read_depth(args.depth) if args.depth else []
        calls = list(read_calls(args.vcf))
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    median = statistics.median(depths) if depths else 0

    rows = []
    for chrom, pos, ref, alt, vaf, dp in calls:
        hits = [name for start, end, name in primers.get(chrom, []) if start <= pos <= end]
        ratio = f"{dp / median:.2f}" if median else "NA"
        rows.append((chrom, pos, ref, alt, vaf, dp, "yes" if hits else "no",
                     ",".join(hits) if hits else "NA", ratio))

    with open(args.output, "w", encoding="utf-8", newline="") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for chrom, pos, ref, alt, vaf, dp, flag, names, ratio in rows:
            out.write("\t".join([chrom, str(pos), ref, alt, f"{vaf:.6f}", str(dp),
                                 flag, names, ratio]) + "\n")

    under = [r for r in rows if r[6] == "yes"]
    listed = [r for r in under if r[4] >= args.min_vaf]
    contigs = sorted(primers)
    unseen = sorted({r[0] for r in rows} - set(contigs))
    print(f"{len(rows)} call(s): {len(under)} under a primer.")
    if unseen:
        print(f"The primer scheme covers {', '.join(contigs)}, and these calls are on "
              f"{', '.join(unseen)}, so none of those could be marked.")
    if listed:
        print(f"Under a primer at a VAF of {args.min_vaf:.2%} or more: {len(listed)}")
        for chrom, pos, ref, alt, vaf, dp, _, names, ratio in sorted(listed, key=lambda r: r[1]):
            tail = f", {ratio}x the median depth" if ratio != "NA" else ""
            print(f"  {chrom}:{pos} {ref}>{alt}  VAF {vaf:.2%}  {dp:,} reads{tail}  under {names}")
        print("The bases under a primer come partly from the primer's own sequence, so check "
              "these against an untrimmed alignment or a second primer scheme before believing them.")


if __name__ == "__main__":
    main()
