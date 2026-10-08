#!/usr/bin/env python3
"""Summarise one stock's reference in a page someone can actually check.

Building a per-stock reference produces several small files, and the question a
reader has is a single one: can this reference be trusted, and how far is this
stock from the lab's own? This gathers the answers into one Markdown page and
one TSV.

It reports. It does not decide: the annotation check already fails the run on a
reading frame the transfer broke, so by the time this runs the reference has
passed. What is left is to say plainly how the stock differs from the reference
it was built against, and to flag anything a reader should look at.
"""
import argparse
import os
import sys

__version__ = "1.0.0"


def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, 'r', encoding='utf-8') as f:
        lines = [l.rstrip('\r\n') for l in f if l.strip()]
    if not lines:
        return []
    keys = lines[0].split('\t')
    return [dict(zip(keys, l.split('\t'))) for l in lines[1:]]


def read_text(path, default='unknown'):
    if not path or not os.path.exists(path):
        return default
    with open(path, 'r', encoding='utf-8') as f:
        return f.read().strip() or default


def count_lines(path):
    if not path or not os.path.exists(path):
        return 0
    with open(path, 'r', encoding='utf-8') as f:
        return sum(1 for l in f if l.strip())


def fasta_stats(path):
    if not path or not os.path.exists(path):
        return None, 0, 0
    name, seq = None, []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
    s = ''.join(seq).upper()
    return name, len(s), s.count('N')


def main():
    parser = argparse.ArgumentParser(description="Summarise one stock's reference.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--stock-name", required=True)
    parser.add_argument("--fasta", required=True, help="The stock reference FASTA")
    parser.add_argument("--annotation-check", required=True, help="Per-CDS check TSV")
    parser.add_argument("--liftover-summary", required=True, help="One-line alignment summary")
    parser.add_argument("--orientation", help="Whether the consensus was reverse-complemented")
    parser.add_argument("--unmapped", help="Features Liftoff could not place")
    parser.add_argument("--primer-report", help="Primer placement report, if a BED was given")
    parser.add_argument("--output-md", required=True)
    parser.add_argument("--output-tsv", required=True)
    args = parser.parse_args()

    for path in (args.fasta, args.annotation_check, args.liftover_summary):
        if not os.path.exists(path):
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)

    contig, length, n_count = fasta_stats(args.fasta)
    cds = read_tsv(args.annotation_check)
    summary = (read_tsv(args.liftover_summary) or [{}])[0]
    orientation = read_text(args.orientation)
    unmapped = count_lines(args.unmapped)
    primers = read_tsv(args.primer_report)

    inherited = [c for c in cds if 'inherited' in (c.get('in_contig', ''), c.get('length_mod3', ''),
                                                   c.get('no_internal_stop', ''))]
    mismatch = int(summary.get('mismatch', 0) or 0)
    insertion = int(summary.get('insertion', 0) or 0)
    deletion = int(summary.get('deletion', 0) or 0)
    identity = summary.get('identity', 'NA')
    ref_contig = summary.get('ref_contig', 'NA')

    flags = []
    if n_count:
        flags.append(f"{n_count} ambiguous base(s) in the consensus")
    if unmapped:
        flags.append(f"{unmapped} annotation feature(s) could not be placed")
    aligned_from = int(summary.get('ref_aligned_from', 1) or 1)
    aligned_to = int(summary.get('ref_aligned_to', 0) or 0)
    if aligned_from > 1:
        flags.append(f"the first {aligned_from - 1} base(s) of {ref_contig} are not covered")
    dropped = [p for p in primers if p.get('status') == 'dropped']
    resized = [p for p in primers if p.get('status') in ('shortened', 'resized')]
    if dropped:
        flags.append(f"{len(dropped)} primer interval(s) could not be placed on this stock")
    if resized:
        flags.append(f"{len(resized)} primer interval(s) changed width")

    stats = {
        'stock': args.stock_name,
        'contig': contig or 'NA',
        'length': length,
        'ambiguous_bases': n_count,
        'orientation': orientation,
        'reference': ref_contig,
        'identity_to_reference': identity,
        'substitutions': mismatch,
        'inserted_vs_reference': insertion,
        'deleted_vs_reference': deletion,
        'cds_transferred': len(cds),
        'cds_unmapped': unmapped,
        'cds_inherited_findings': len(inherited),
        'primers_placed': len(primers) - len(dropped) if primers else 'NA',
        'primers_dropped': len(dropped) if primers else 'NA',
        'flags': '; '.join(flags) if flags else 'none',
    }

    with open(args.output_tsv, 'w', encoding='utf-8', newline='') as f:
        f.write('\t'.join(stats) + '\n')
        f.write('\t'.join(str(v) for v in stats.values()) + '\n')

    md = [f"# Stock reference: {args.stock_name}", ""]
    md.append(f"Built from a de novo consensus and the annotation of `{ref_contig}`. "
              f"Use `{args.stock_name}.fasta` and `{args.stock_name}.gff3` as `--fasta` and "
              f"`--gff` when calling variants in this stock's samples.")
    md += ["", "## The sequence", "",
           "| | |", "| :--- | :--- |",
           f"| Contig | `{contig}` |",
           f"| Length | {length:,} bp |",
           f"| Ambiguous bases | {n_count} |",
           f"| Consensus used | {orientation} |", ""]

    md += ["## Distance from the lab reference", "",
           f"Against `{ref_contig}`, identity {identity} over the aligned region "
           f"({aligned_from:,}-{aligned_to:,}).", "",
           "| | |", "| :--- | :--- |",
           f"| Substitutions | {mismatch} |",
           f"| Bases only in the stock | {insertion} |",
           f"| Bases only in the reference | {deletion} |", "",
           "These are the stock's fixed differences from the lab reference. Calling this "
           "stock against the reference would report every one of them as a variant. "
           "Calling against this file does not, which is the point of building it.", ""]

    md += ["## Annotation", "",
           f"{len(cds)} CDS transferred, {unmapped} feature(s) unplaced.", ""]
    if cds:
        md += ["| CDS | Position | In contig | In frame | No internal stop |",
               "| :--- | :--- | :--- | :--- | :--- |"]
        for c in cds:
            md.append(f"| {c.get('cds')} | {c.get('start')}-{c.get('end')} | "
                      f"{c.get('in_contig')} | {c.get('length_mod3')} | "
                      f"{c.get('no_internal_stop')} |")
        md.append("")
    if inherited:
        md += ["`inherited` means the lab reference has the same finding, so the transfer did "
               "not cause it. For these genomes that is expected: nsP3 ends in an opal (TGA) "
               "codon the virus reads through.", ""]

    if primers:
        md += ["## Primers", "",
               f"{len(primers) - len(dropped)} of {len(primers)} interval(s) placed on this stock.", ""]
        if dropped:
            md += ["Could not be placed, so those amplicons cannot be trimmed:", ""]
            md += [f"- `{p['name']}` at {p['source_contig']}:{p['source_start']}-{p['source_end']}"
                   for p in dropped]
            md.append("")
        if resized:
            md += ["Changed width. Some of the reference bases under these primers are not in "
                   "this stock, either deleted or outside the assembled region, so the primer "
                   "covers fewer bases here than it was designed to:", ""]
            md += [f"- `{p['name']}` {p['source_width']} bp to {p['target_width']} bp"
                   for p in resized]
            md.append("")

    md += ["## Worth a look", ""]
    md += ([f"- {f}" for f in flags] if flags else ["Nothing flagged."])
    md.append("")

    with open(args.output_md, 'w', encoding='utf-8', newline='') as f:
        f.write('\n'.join(md))

    print(f"{args.stock_name}: {length} bp, {mismatch} substitution(s), "
          f"{insertion} inserted, {deletion} deleted vs {ref_contig}; "
          f"{len(cds)} CDS; flags: {stats['flags']}")


if __name__ == "__main__":
    main()
