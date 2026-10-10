#!/usr/bin/env python3
"""Flag LoFreq calls made where the depth cap let LoFreq examine only some of the reads.

`lofreq call` stops counting reads at --max-depth, 1,000,000 by default. At a position
deeper than that, the allele counts, and so the VAF, come from a subset of the reads that
cover it, and the VAF can read well above the truth. On a library sequenced to about 1.2
million times coverage, a position at 4.4 million times gave 1.8% from LoFreq's capped
subset and 0.5% from iVar at full depth. With the cap raised to 10 million LoFreq gave 0.56%.

This script compares the cap with the real depth from `samtools depth`. LoFreq can examine at
most max_depth reads, so fraction_examined = min(1, max_depth / real_depth). A call is marked
capped when that is under --ratio. The comparison does not use LoFreq's own INFO/DP, because
LoFreq also drops reads on base and mapping quality, which lowers DP at depths well under the cap.

Input:
  --vcf        LoFreq VCF, plain or gzip. VAF comes from INFO/DP4 when present, else INFO/AF.
  --depth      samtools depth table: contig, 1-based position, depth.
  --max-depth  The cap LoFreq ran with. Default: the value in the VCF's ##source line, and
               LoFreq's own default, 1000000, when the line has none.

Output: one row per call, tab-separated, and a summary on stdout.
  depth_capped  yes     the cap let LoFreq examine under --ratio of the reads
                no      it did not
                unknown the position is not in the depth table, or its real depth is 0
"""
import argparse
import gzip
import os
import re
import sys

__version__ = "2.0.0"

LOFREQ_DEFAULT_MAX_DEPTH = 1000000

COLUMNS = ["chrom", "pos", "ref", "alt", "lofreq_vaf", "lofreq_depth", "real_depth",
           "max_depth", "fraction_examined", "depth_capped"]

MAX_DEPTH_RE = re.compile(r"--max-depth[ =](\d+)")


def opener(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def read_depth(path):
    """Return {(contig, pos): depth}."""
    depth = {}
    with opener(path) as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_no}: expected contig, position and depth")
            try:
                depth[(fields[0], int(fields[1]))] = int(float(fields[2]))
            except ValueError:
                raise ValueError(f"{path}:{line_no}: position and depth must be numbers") from None
    return depth


def read_vcf(path):
    """Return (max_depth from the header or None, [(chrom, pos, ref, alt, vaf, lofreq_depth), ...])."""
    header_max_depth = None
    calls = []
    with opener(path) as f:
        for line_no, line in enumerate(f, start=1):
            if line.startswith("##"):
                if header_max_depth is None and line.startswith("##source="):
                    match = MAX_DEPTH_RE.search(line)
                    if match:
                        header_max_depth = int(match.group(1))
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_no}: a VCF record needs 8 columns, found {len(fields)}")
            info = dict(kv.split("=", 1) for kv in fields[7].split(";") if "=" in kv)
            if "DP" not in info:
                raise ValueError(f"{path}:{line_no}: no INFO/DP, so this is not a LoFreq VCF")
            try:
                pos, dp = int(fields[1]), int(info["DP"])
                dp4 = [int(x) for x in info["DP4"].split(",")] if "DP4" in info else None
                if dp4 and len(dp4) == 4 and sum(dp4) > 0:
                    vaf = (dp4[2] + dp4[3]) / sum(dp4)
                else:
                    vaf = float(info.get("AF", "0").split(",")[0])
            except ValueError:
                raise ValueError(f"{path}:{line_no}: could not read POS, DP, DP4 or AF") from None
            calls.append((fields[0], pos, fields[3], fields[4], vaf, dp))
    return header_max_depth, calls


def main():
    parser = argparse.ArgumentParser(
        description="Flag LoFreq calls made where the depth cap let LoFreq examine only some of the reads.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--vcf", required=True, help="LoFreq VCF, plain or gzip")
    parser.add_argument("--depth", required=True, help="samtools depth table: contig, position, depth")
    parser.add_argument("--output", required=True, help="Output TSV, one row per call")
    parser.add_argument("--max-depth", type=int, default=None,
                        help="The cap LoFreq ran with. Default: read from the VCF header, else 1000000")
    parser.add_argument("--ratio", type=float, default=0.8,
                        help="A call is capped when the cap lets LoFreq examine under this fraction "
                             "of the reads. Default: 0.8")
    parser.add_argument("--min-vaf", type=float, default=0.01,
                        help="Capped calls at or above this VAF are listed in the summary. Default: 0.01")
    args = parser.parse_args()

    if not 0 < args.ratio <= 1:
        print("ERROR: --ratio must be above 0 and at most 1", file=sys.stderr)
        sys.exit(1)
    if args.max_depth is not None and args.max_depth < 1:
        print("ERROR: --max-depth must be at least 1", file=sys.stderr)
        sys.exit(1)
    for path in (args.vcf, args.depth):
        if not os.path.exists(path):
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    try:
        depth = read_depth(args.depth)
        header_max_depth, calls = read_vcf(args.vcf)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    if args.max_depth is not None:
        max_depth, where = args.max_depth, "from --max-depth"
    elif header_max_depth is not None:
        max_depth, where = header_max_depth, "from the VCF header"
    else:
        max_depth, where = LOFREQ_DEFAULT_MAX_DEPTH, "LoFreq's default, because the VCF header has none"

    rows = []
    for chrom, pos, ref, alt, vaf, dp in calls:
        real = depth.get((chrom, pos))
        if not real:
            rows.append((chrom, pos, ref, alt, vaf, dp, "NA" if real is None else real, max_depth, "NA", "unknown"))
            continue
        fraction = min(1.0, max_depth / real)
        rows.append((chrom, pos, ref, alt, vaf, dp, real, max_depth, fraction,
                     "yes" if fraction < args.ratio else "no"))

    with open(args.output, "w", encoding="utf-8", newline="") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for chrom, pos, ref, alt, vaf, dp, real, cap, fraction, flag in rows:
            out.write("\t".join([
                chrom, str(pos), ref, alt, f"{vaf:.6f}", str(dp), str(real), str(cap),
                fraction if isinstance(fraction, str) else f"{fraction:.4f}", flag]) + "\n")

    capped = [r for r in rows if r[9] == "yes"]
    unknown = [r for r in rows if r[9] == "unknown"]
    listed = [r for r in capped if r[4] >= args.min_vaf]
    print(f"LoFreq depth cap: {max_depth:,} ({where}).")
    print(f"{len(rows)} call(s): {len(capped)} at positions where the cap let LoFreq examine under "
          f"{args.ratio:.0%} of the reads, {len(unknown)} with no depth to compare.")
    if listed:
        print(f"Capped calls with a VAF of {args.min_vaf:.2%} or more: {len(listed)}")
        for chrom, pos, ref, alt, vaf, dp, real, cap, fraction, _ in sorted(listed, key=lambda r: r[1]):
            print(f"  {chrom}:{pos} {ref}>{alt}  LoFreq VAF {vaf:.2%}  "
                  f"{real:,} reads, LoFreq examined at most {cap:,} ({fraction:.0%})")
    if capped:
        print("A capped VAF comes from a subset of the reads and can read high. "
              "Check these calls against iVar, or raise --lofreq_max_depth to the deepest position.")


if __name__ == "__main__":
    main()
