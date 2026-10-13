#!/usr/bin/env python3
"""Summarise a purity check of one virus stock's reads.

A stock is meant to hold one virus. This puts three numbers beside each other:

  1. How many reads map to the stock's own reference. A pure stock puts nearly all of
     them there. Most of the rest are usually sequencing artefacts, not another virus:
     poly-G reads, where the instrument saw no signal, and reads that are adapter. Both
     are counted and set aside, so they are not mistaken for contamination.
  2. Where the reads go when they can choose between the lab's other genomes (the
     reference pool). A stock of one virus puts them on the genome closest to it, and
     almost none on the rest. A second genome with a real share points to another
     strain in the tube. Kraken2 cannot tell close strains apart. Mapping can.
  3. What the reads left over are, from a Kraken2 report. This finds other viruses, and
     bacteria or host when the database holds them.

It reports. It does not decide. A flag says where to look. A stock can be sound and
still carry a flag, for example when the pool holds two strains that are close to each
other, so read the numbers.

Input counts files have the columns contig, length, mapped, mapq0, one row per reference
contig, and a last row `*` whose `mapped` column holds the reads that did not map. `mapped`
counts primary alignments. `mapq0` is how many of those have mapping quality 0, which means
the read fits two places equally well. The unmapped-kinds file has the columns kind and reads,
with rows total, polyg and adapter, counted among the reads that did not map.
"""
import argparse
import os
import sys

__version__ = "2.0.0"

COLUMNS = ["stock", "reads_sampled", "reads_mapped_to_stock", "fraction_mapped_to_stock",
           "reads_unmapped_to_stock", "unmapped_polyg", "unmapped_adapter",
           "fraction_unexplained", "pool_top_genome", "pool_top_share",
           "pool_other_max_genome", "pool_other_max_share", "kraken_reads",
           "kraken_unclassified_fraction", "kraken_top_species", "kraken_top_species_reads",
           "flags"]


def read_counts(path):
    """Return ([(contig, length, mapped, mapq0), ...], unmapped)."""
    rows, unmapped = [], 0
    with open(path, "r", encoding="utf-8") as f:
        lines = [l.rstrip("\r\n") for l in f if l.strip()]
    if not lines or lines[0].split("\t")[:4] != ["contig", "length", "mapped", "mapq0"]:
        raise ValueError(f"{path}: expected a header of contig, length, mapped, mapq0")
    for line_no, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) < 4:
            raise ValueError(f"{path}:{line_no}: expected 4 columns, found {len(fields)}")
        try:
            length, mapped, mapq0 = int(fields[1]), int(fields[2]), int(fields[3])
        except ValueError:
            raise ValueError(f"{path}:{line_no}: length, mapped and mapq0 must be whole numbers") from None
        if mapq0 > mapped:
            raise ValueError(f"{path}:{line_no}: mapq0 ({mapq0}) is larger than mapped ({mapped})")
        if fields[0] == "*":
            unmapped = mapped
        else:
            rows.append((fields[0], length, mapped, mapq0))
    return rows, unmapped


def read_kinds(path):
    """Return {'total': n, 'polyg': n, 'adapter': n} from the unmapped-kinds file."""
    kinds = {}
    with open(path, "r", encoding="utf-8") as f:
        lines = [l.rstrip("\r\n") for l in f if l.strip()]
    if not lines or lines[0].split("\t")[:2] != ["kind", "reads"]:
        raise ValueError(f"{path}: expected a header of kind, reads")
    for line_no, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) < 2:
            raise ValueError(f"{path}:{line_no}: expected 2 columns, found {len(fields)}")
        try:
            kinds[fields[0]] = int(fields[1])
        except ValueError:
            raise ValueError(f"{path}:{line_no}: reads must be a whole number") from None
    for key in ("total", "polyg", "adapter"):
        if key not in kinds:
            raise ValueError(f"{path}: no row for {key}")
    if kinds["polyg"] + kinds["adapter"] > kinds["total"]:
        raise ValueError(f"{path}: poly-G and adapter reads add up to more than the total")
    return kinds


def read_kraken(path):
    """Parse a Kraken2 report. Return dict with total, unclassified, domains, species."""
    unclassified, root, domains, species = 0, 0, [], []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 6:
                continue
            try:
                clade = int(fields[1])
            except ValueError:
                continue
            rank, name = fields[3], fields[5].strip()
            if rank == "U":
                unclassified = clade
            elif rank == "R":
                root = clade
            elif rank == "D":
                # Kraken2 rank codes: D is a domain, D1 a rank between domain and kingdom,
                # K a kingdom. The viral database has kingdoms under Viruses, and they are
                # not domains.
                domains.append((name, clade))
            elif rank == "S":
                species.append((name, clade))
    species.sort(key=lambda s: -s[1])
    domains.sort(key=lambda s: -s[1])
    return {"total": unclassified + root, "unclassified": unclassified,
            "domains": domains, "species": species}


def pct(part, whole):
    return f"{100 * part / whole:.2f}%" if whole else "NA"


def main():
    parser = argparse.ArgumentParser(description="Summarise a purity check of one virus stock's reads.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--stock-name", required=True)
    parser.add_argument("--stock-counts", required=True,
                        help="Counts of the reads mapped to the stock reference")
    parser.add_argument("--unmapped-kinds",
                        help="How many of the reads that did not map are poly-G or adapter reads")
    parser.add_argument("--pool-counts", help="Counts of the same reads mapped to the reference pool")
    parser.add_argument("--kraken-report", help="Kraken2 report of the reads left after the artefacts")
    parser.add_argument("--min-mapped", type=float, default=0.90,
                        help="Flag when a smaller fraction of the reads maps to the stock or is a "
                             "sequencing artefact. Default: 0.90")
    parser.add_argument("--flag-fraction", type=float, default=0.01,
                        help="Flag a pool genome other than the closest one, or a Kraken2 species, "
                             "when it holds at least this fraction of the reads it is measured "
                             "against. Default: 0.01")
    parser.add_argument("--min-reads", type=int, default=100,
                        help="A genome or species needs at least this many reads to be flagged, "
                             "so a handful of reads cannot raise a flag. Default: 100")
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-md", required=True)
    args = parser.parse_args()

    if not 0 < args.min_mapped <= 1:
        print("ERROR: --min-mapped must be above 0 and at most 1", file=sys.stderr)
        sys.exit(1)
    if not 0 < args.flag_fraction <= 1:
        print("ERROR: --flag-fraction must be above 0 and at most 1", file=sys.stderr)
        sys.exit(1)
    if args.min_reads < 1:
        print("ERROR: --min-reads must be at least 1", file=sys.stderr)
        sys.exit(1)
    for path in (args.stock_counts, args.unmapped_kinds, args.pool_counts, args.kraken_report):
        if path and not os.path.exists(path):
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    try:
        stock_rows, unmapped = read_counts(args.stock_counts)
        kinds = read_kinds(args.unmapped_kinds) if args.unmapped_kinds else None
        pool_rows = read_counts(args.pool_counts)[0] if args.pool_counts else None
        kraken = read_kraken(args.kraken_report) if args.kraken_report else None
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    mapped = sum(r[2] for r in stock_rows)
    total = mapped + unmapped
    if total == 0:
        print("ERROR: the stock counts hold no reads", file=sys.stderr)
        sys.exit(1)
    if kinds is not None and kinds["total"] != unmapped:
        print(f"ERROR: the unmapped-kinds file counts {kinds['total']} unmapped reads, "
              f"the stock counts {unmapped}", file=sys.stderr)
        sys.exit(1)

    fraction = mapped / total
    polyg = kinds["polyg"] if kinds else 0
    adapter = kinds["adapter"] if kinds else 0
    artefacts = polyg + adapter
    unexplained = unmapped - artefacts
    explained = (mapped + artefacts) / total
    flags = []
    if explained < args.min_mapped:
        text = f"only {fraction:.1%} of the reads map to the stock reference"
        if artefacts:
            text += f", and {artefacts / total:.1%} more are poly-G or adapter reads"
        flags.append(text)

    # Reads on the pool. A read that fits two genomes equally well has mapping quality 0
    # and says nothing about which one it came from, so shares use the other reads.
    pool_view, top_name, top_share = [], "NA", "NA"
    other_name, other_share = "NA", "NA"
    if pool_rows is not None:
        unique_total = sum(r[2] - r[3] for r in pool_rows)
        pool_view = sorted(((name, length, mapped_, q0, mapped_ - q0)
                            for name, length, mapped_, q0 in pool_rows), key=lambda r: -r[4])
        if pool_view and unique_total:
            top_name = pool_view[0][0]
            top_share = f"{pool_view[0][4] / unique_total:.4f}"
            others = pool_view[1:]
            if others:
                best = others[0]
                other_name, other_share = best[0], f"{best[4] / unique_total:.4f}"
                for name, _, _, _, unique in others:
                    if unique >= args.min_reads and unique / unique_total >= args.flag_fraction:
                        flags.append(f"{unique:,} reads ({unique / unique_total:.1%} of the reads placed "
                                     f"uniquely on the pool) fit {name} best, not {top_name}, which is "
                                     f"the genome closest to the stock")
        elif pool_view:
            flags.append("no read maps uniquely to the reference pool")

    k_unclassified, k_top, k_top_reads = "NA", "NA", "NA"
    if kraken is not None:
        if kraken["total"]:
            k_unclassified = f"{kraken['unclassified'] / kraken['total']:.4f}"
        if kraken["species"]:
            k_top, k_top_reads = kraken["species"][0][0], str(kraken["species"][0][1])
        for name, reads in kraken["species"]:
            if reads >= args.min_reads and reads / total >= args.flag_fraction:
                flags.append(f"{reads:,} reads ({reads / total:.1%} of all) classify as {name}")

    stats = {
        "stock": args.stock_name,
        "reads_sampled": total,
        "reads_mapped_to_stock": mapped,
        "fraction_mapped_to_stock": f"{fraction:.4f}",
        "reads_unmapped_to_stock": unmapped,
        "unmapped_polyg": polyg if kinds else "NA",
        "unmapped_adapter": adapter if kinds else "NA",
        "fraction_unexplained": f"{unexplained / total:.4f}",
        "pool_top_genome": top_name,
        "pool_top_share": top_share,
        "pool_other_max_genome": other_name,
        "pool_other_max_share": other_share,
        "kraken_reads": kraken["total"] if kraken is not None else "NA",
        "kraken_unclassified_fraction": k_unclassified,
        "kraken_top_species": k_top,
        "kraken_top_species_reads": k_top_reads,
        "flags": "; ".join(flags) if flags else "none",
    }
    with open(args.output_tsv, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(COLUMNS) + "\n")
        f.write("\t".join(str(stats[c]) for c in COLUMNS) + "\n")

    md = ["## Purity", "",
          f"The first {total // 2:,} read pairs ({total:,} reads) of the stock's sequencing run. "
          "The counts are of reads, so a pair counts twice.", "",
          "### Reads against the stock reference", "",
          "| | Reads | Share |", "| :--- | ---: | ---: |",
          f"| Map to the stock reference | {mapped:,} | {pct(mapped, total)} |",
          f"| Do not map | {unmapped:,} | {pct(unmapped, total)} |"]
    if kinds is not None:
        md += [f"| of which poly-G reads | {polyg:,} | {pct(polyg, total)} |",
               f"| of which carry an Illumina adapter | {adapter:,} | {pct(adapter, total)} |",
               f"| of which neither | {unexplained:,} | {pct(unexplained, total)} |"]
    md.append("")
    if kinds is not None:
        md += ["Poly-G reads are where the instrument saw no signal and read G. Adapter reads are "
               "inserts shorter than the read. Neither is another virus, so they are counted and set "
               "aside before the reads are classified.", ""]

    if pool_rows is not None:
        pool_total = sum(r[2] for r in pool_rows)
        md += ["### Reads against the reference pool", "",
               "Each read goes to the genome it fits best. A stock of one virus puts nearly all "
               "of its reads on the genome closest to it. The reads with mapping quality 0 fit "
               "two genomes equally well, so shares are of the other reads.", "",
               f"{pool_total:,} reads map to the pool in all.", "",
               "| Genome | Reads | Mapping quality 0 | Share of the rest |",
               "| :--- | ---: | ---: | ---: |"]
        unique_total = sum(r[4] for r in pool_view)
        for name, _, mapped_, q0, unique in pool_view[:8]:
            md.append(f"| `{name}` | {mapped_:,} | {q0:,} | {pct(unique, unique_total)} |")
        if len(pool_view) > 8:
            md.append(f"| {len(pool_view) - 8} more genome(s) | | | |")
        md.append("")

    if kraken is not None:
        md += ["### What the unmapped reads are", ""]
        if not kraken["total"]:
            md += ["There were no reads left to classify.", ""]
        else:
            left = "reads left after the poly-G and adapter reads" if kinds is not None \
                else "reads that did not map to the stock"
            md += [f"Kraken2 classified {kraken['total'] - kraken['unclassified']:,} of the "
                   f"{kraken['total']:,} {left}. "
                   f"{pct(kraken['unclassified'], kraken['total'])} are unclassified.", ""]
            if kraken["domains"]:
                md += ["| Domain | Reads |", "| :--- | ---: |"]
                md += [f"| {n} | {c:,} |" for n, c in kraken["domains"] if c]
                md.append("")
            if kraken["species"]:
                md += ["| Top species | Reads | Share of all reads |", "| :--- | ---: | ---: |"]
                md += [f"| {n} | {c:,} | {pct(c, total)} |" for n, c in kraken["species"][:5]]
                md.append("")

    md += ["### Flags", ""]
    md += ([f"- {f}" for f in flags] if flags else ["Nothing flagged."])
    md.append("")

    with open(args.output_md, "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(md))

    print(f"{args.stock_name}: {fraction:.1%} of {total:,} reads map to the stock; "
          f"flags: {stats['flags']}")


if __name__ == "__main__":
    main()
