#!/usr/bin/env python3
"""Check an annotation transferred onto a stock consensus.

A per-stock reference is only usable if its CDS features still describe the same
reading frames as the reference they came from. This compares the transferred
GFF3 against the stock FASTA, and writes one row per CDS.

Checks applied to every CDS:

  in_contig        coordinates lie inside the stock contig
  length_mod3      length is a multiple of three, so the frame is intact
  no_internal_stop no stop codon before the last codon

The question worth answering is whether the transfer broke anything, not whether
the annotation was perfect to begin with. Give --reference-fasta with
--reference-gff and the same checks run on the reference first. A result that
already failed there is reported as `inherited` and does not fail the run; only
a new failure does.

Two points of alphavirus biology matter here, and both are why the reference
control is not optional:

  * The mature peptides (nsP1..nsP4, Capsid, E3, E2, 6K, E1) are cleavage
    products of two polyproteins. An individual peptide does not begin with ATG
    or end in a stop codon, so start and stop are reported for information and
    never fail. Liftoff's own `valid_ORF=False` is ignored for the same reason.
  * nsP3 ends in an opal (TGA) codon that the virus reads through to make
    nsP1234. It is a real stop codon inside a real CDS, present in the lab's
    TC-83, VEEV INH-9813 and EEEV V105 references. Against a reference control it
    is `inherited`; without one it would be a false alarm.

Exit status is 1 on a new failure, so a pipeline stops before a broken reference
reaches variant calling.
"""
import argparse
import os
import sys

__version__ = "1.1.0"

STOPS = {'TAA', 'TAG', 'TGA'}
COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CHECKS = ('in_contig', 'length_mod3', 'no_internal_stop')


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks)
                name = line[1:].split()[0] if len(line) > 1 else ''
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(chunks)
    return seqs


def read_cds(path):
    """Return CDS features as dicts, in file order."""
    feats = []
    with open(path, 'r', encoding='utf-8') as f:
        for line_no, line in enumerate(f, start=1):
            if line.startswith('#'):
                continue
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            attrs = {}
            for item in fields[8].split(';'):
                if '=' in item:
                    k, v = item.split('=', 1)
                    attrs[k.strip()] = v.strip()
            feats.append({
                'line_no': line_no, 'seqid': fields[0],
                'start': int(fields[3]), 'end': int(fields[4]),
                'strand': fields[6],
                'name': attrs.get('Name') or attrs.get('ID') or f"CDS@{fields[3]}",
            })
    return feats


def cds_sequence(seq, feat):
    sub = seq[feat['start'] - 1:feat['end']]
    if feat['strand'] == '-':
        sub = sub.translate(COMPLEMENT)[::-1]
    return sub.upper()


def check_feature(feat, seqs):
    """Return a result dict for one CDS. Each name in CHECKS maps to True/False."""
    res = {
        'cds': feat['name'], 'seqid': feat['seqid'],
        'start': feat['start'], 'end': feat['end'], 'strand': feat['strand'],
        'length': feat['end'] - feat['start'] + 1,
        'starts_atg': 'NA', 'ends_stop': 'NA', 'internal_stop_codon': 'NA',
        'detail': {},
    }
    for name in CHECKS:
        res[name] = None  # not evaluated

    seq = seqs.get(feat['seqid'])
    if seq is None:
        res['in_contig'] = False
        res['detail']['in_contig'] = f"sequence {feat['seqid']} is not in the FASTA"
        return res

    in_contig = feat['start'] >= 1 and feat['end'] <= len(seq) and feat['start'] <= feat['end']
    res['in_contig'] = in_contig
    if not in_contig:
        res['detail']['in_contig'] = (f"coordinates {feat['start']}-{feat['end']} fall outside "
                                      f"{feat['seqid']} (length {len(seq)})")
        return res

    length_ok = res['length'] % 3 == 0
    res['length_mod3'] = length_ok
    if not length_ok:
        res['detail']['length_mod3'] = (f"length {res['length']} is not a multiple of three, so "
                                        "the reading frame is broken")

    sub = cds_sequence(seq, feat)
    res['starts_atg'] = 'yes' if sub[:3] == 'ATG' else 'no'
    res['ends_stop'] = 'yes' if sub[-3:] in STOPS else 'no'

    if length_ok:
        internal = [i for i in range(0, len(sub) - 3, 3) if sub[i:i + 3] in STOPS]
        res['no_internal_stop'] = not internal
        if internal:
            first = internal[0]
            res['internal_stop_codon'] = f"codon {first // 3 + 1} ({sub[first:first + 3]})"
            res['detail']['no_internal_stop'] = (
                f"stop codon {sub[first:first + 3]} at codon {first // 3 + 1} of "
                f"{len(sub) // 3}, before the end of the feature")
        else:
            res['internal_stop_codon'] = 'none'
    return res


def run_checks(fasta, gff):
    seqs = read_fasta(fasta)
    feats = read_cds(gff)
    return [check_feature(f, seqs) for f in feats], feats


def main():
    parser = argparse.ArgumentParser(
        description="Check a GFF3 transferred onto a stock consensus against that consensus.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--fasta", required=True, help="Stock consensus FASTA")
    parser.add_argument("--gff", required=True, help="Transferred GFF3")
    parser.add_argument("--reference-fasta", help="Source FASTA, used as a control")
    parser.add_argument("--reference-gff", help="Source GFF3, used as a control and to confirm "
                                                "no CDS was lost")
    parser.add_argument("--output-tsv", required=True, help="One row per CDS")
    parser.add_argument("--require-start-stop", action="store_true",
                        help="Also require each CDS to start with ATG and end in a stop codon. "
                             "Off by default: alphavirus mature peptides are polyprotein "
                             "cleavage products and have neither.")
    args = parser.parse_args()

    needed = [args.fasta, args.gff]
    if args.reference_fasta:
        needed.append(args.reference_fasta)
    if args.reference_gff:
        needed.append(args.reference_gff)
    for path in needed:
        if not os.path.exists(path):
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)
    if bool(args.reference_fasta) != bool(args.reference_gff):
        print("ERROR: --reference-fasta and --reference-gff must be given together",
              file=sys.stderr)
        sys.exit(1)

    results, feats = run_checks(args.fasta, args.gff)
    if not feats:
        print(f"ERROR: {args.gff} has no CDS features", file=sys.stderr)
        sys.exit(1)

    # Control: the same checks on the reference. A failure already present there
    # came with the annotation, not with the transfer.
    inherited = set()
    if args.reference_fasta:
        ref_results, _ = run_checks(args.reference_fasta, args.reference_gff)
        for r in ref_results:
            for name in CHECKS:
                if r[name] is False:
                    inherited.add((r['cds'], name))

    failures, notes = [], []
    for res in results:
        for name in CHECKS:
            if res[name] is False:
                msg = f"{res['cds']}: {res['detail'].get(name, name + ' failed')}"
                if (res['cds'], name) in inherited:
                    res[name] = 'inherited'
                    notes.append(msg)
                else:
                    res[name] = 'FAIL'
                    failures.append(msg)
            elif res[name] is True:
                res[name] = 'pass'
            else:
                res[name] = 'NA'
        if args.require_start_stop:
            if res['starts_atg'] == 'no':
                failures.append(f"{res['cds']}: does not start with ATG")
            if res['ends_stop'] == 'no':
                failures.append(f"{res['cds']}: does not end with a stop codon")

    if args.reference_gff:
        ref_names = {f['name'] for f in read_cds(args.reference_gff)}
        lost = sorted(ref_names - {f['name'] for f in feats})
        if lost:
            failures.append(f"CDS present in the reference but not transferred: {', '.join(lost)}")

    cols = ['cds', 'seqid', 'start', 'end', 'strand', 'length', 'in_contig', 'length_mod3',
            'no_internal_stop', 'internal_stop_codon', 'starts_atg', 'ends_stop']
    with open(args.output_tsv, 'w', encoding='utf-8', newline='') as f:
        f.write('\t'.join(cols) + '\n')
        for res in results:
            f.write('\t'.join(str(res.get(c, 'NA')) for c in cols) + '\n')

    for msg in notes:
        print(f"note (also present in the reference): {msg}")

    if failures:
        print(f"FAIL: {len(failures)} new problem(s) in {args.gff}", file=sys.stderr)
        for msg in failures:
            print(f"  - {msg}", file=sys.stderr)
        sys.exit(1)

    carried = f", {len(notes)} inherited from the reference" if notes else ""
    print(f"OK: {len(results)} CDS features check out against {args.fasta}{carried}")


if __name__ == "__main__":
    main()
