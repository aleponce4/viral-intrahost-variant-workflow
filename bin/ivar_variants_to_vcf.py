#!/usr/bin/env python3
import sys
import os
import argparse
import math
from datetime import datetime

__version__ = "1.1.0"

# iVar variants columns this script reads, looked up by name. iVar's layout is
# REGION POS REF ALT REF_DP REF_RV REF_QUAL ALT_DP ALT_RV ALT_QUAL ALT_FREQ
# TOTAL_DP PVAL PASS, then GFF columns when a GFF is given. Position-based
# indexing read ALT_QUAL as ALT_FREQ and REF_RV as ALT_DP.
REQUIRED_COLUMNS = (
    'REGION', 'POS', 'REF', 'ALT', 'REF_DP', 'ALT_DP',
    'ALT_FREQ', 'TOTAL_DP', 'PVAL', 'PASS',
)

def _count(value):
    """A depth from an iVar column. iVar prints counts of a million or more in
    scientific notation ("1.16389e+06"), which str.isdigit() rejects."""
    return int(round(float(value)))

def parse_ivar_tsv(tsv_file):
    """Parse iVar TSV file by column name and extract variant information.

    Raises ValueError on a missing column or an unparseable value, so a bad
    file stops the run instead of producing zeros.
    """
    variants = []
    with open(tsv_file, 'r', encoding='utf-8') as f:
        header = f.readline().rstrip('\r\n').split('\t')
        missing = [c for c in REQUIRED_COLUMNS if c not in header]
        if missing:
            raise ValueError(
                f"{tsv_file}: missing required iVar column(s): {', '.join(missing)}"
            )
        col = {name: header.index(name) for name in REQUIRED_COLUMNS}
        width = max(col.values()) + 1
        for line_no, line in enumerate(f, start=2):
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            fields = line.split('\t')
            if len(fields) < width:
                raise ValueError(
                    f"{tsv_file}:{line_no}: expected at least {width} columns, found {len(fields)}"
                )
            try:
                variants.append({
                    'CHROM': fields[col['REGION']],
                    'POS': int(fields[col['POS']]),
                    'REF': fields[col['REF']],
                    'ALT': fields[col['ALT']],
                    'REF_DP': _count(fields[col['REF_DP']]),
                    'ALT_DP': _count(fields[col['ALT_DP']]),
                    'ALT_FREQ': float(fields[col['ALT_FREQ']]),
                    'TOTAL_DP': _count(fields[col['TOTAL_DP']]),
                    'PVAL': float(fields[col['PVAL']]),
                    'PASS': fields[col['PASS']],
                })
            except ValueError as exc:
                raise ValueError(f"{tsv_file}:{line_no}: {exc}") from None
    return variants

def write_vcf_header(output_file, reference_file, sample_name):
    """Write VCF header."""
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        f.write("##source=ivar_variants_to_vcf.py\n")
        f.write(f"##reference={reference_file}\n")
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n")
        f.write("##FORMAT=<ID=ALT_FREQ,Number=1,Type=Float,Description=\"Alternative allele frequency\">\n")
        f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        f.write("##FILTER=<ID=FAIL,Description=\"Failed quality filters\">\n")
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

def convert_to_vcf(variants, output_file, reference_file, sample_name):
    """Convert parsed variants to VCF format."""
    write_vcf_header(output_file, reference_file, sample_name)
    with open(output_file, 'a', encoding='utf-8') as f:
        for variant in variants:
            chrom = variant['CHROM']
            pos = variant['POS']
            ref = variant['REF']
            alt = variant['ALT']
            pval = variant.get('PVAL', 1.0)
            qual = -10.0 * math.log10(max(pval, 1e-300))
            qual_str = f"{qual:.2f}"
            filter_field = "PASS" if variant['PASS'] == 'TRUE' else "FAIL"
            total_dp = variant['TOTAL_DP']
            alt_freq = variant['ALT_FREQ']
            info = f"DP={total_dp};AF={alt_freq}"
            format_field = "GT:DP:AD:ALT_FREQ"
            genotype = "1"
            ref_dp = variant['REF_DP']
            alt_dp = variant['ALT_DP']
            sample_data = f"{genotype}:{total_dp}:{ref_dp},{alt_dp}:{alt_freq}"
            vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual_str}\t{filter_field}\t{info}\t{format_field}\t{sample_data}\n"
            f.write(vcf_line)

def main():
    parser = argparse.ArgumentParser(description="Convert iVar TSV output to VCF format.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--input-tsv", required=True, help="Input iVar variants TSV file")
    parser.add_argument("--output-vcf", required=True, help="Output VCF file")
    parser.add_argument("--reference-fasta", required=True, help="Reference FASTA file")

    args = parser.parse_args()

    if not os.path.exists(args.input_tsv):
        print(f"ERROR: Input TSV file not found: {args.input_tsv}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.reference_fasta):
        print(f"ERROR: Reference file not found: {args.reference_fasta}", file=sys.stderr)
        sys.exit(1)

    sample_name = os.path.splitext(os.path.basename(args.input_tsv))[0]
    try:
        variants = parse_ivar_tsv(args.input_tsv)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
    convert_to_vcf(variants, args.output_vcf, args.reference_fasta, sample_name)
    print(f"Successfully converted {len(variants)} variants from {args.input_tsv} to {args.output_vcf}")

if __name__ == "__main__":
    main()
