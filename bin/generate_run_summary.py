#!/usr/bin/env python3
"""
generate_run_summary.py — Consolidated per-sample run summary.

Derives one row per sample from the QC artefacts staged into the task directory by the
REPORTING subworkflow:

  * ``*.qc_stats.txt``        — emitted by LOFREQ_FILTER (sample, contig, variant counts,
                                status)
  * ``*coverage_summary.tsv`` — emitted by COVERAGE_SUMMARIZE (mean/min/max depth and
                                depth-breadth fractions), merged on sample id when present

Samples are discovered from the staged files. Nothing is invented: if no QC artefacts are
staged, a header-only TSV is written and a warning is emitted on stderr.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

__version__ = "2.0.0"

QC_FIELDS = {
    "sample": "sample",
    "contig": "contig",
    "raw variants": "raw_variants",
    "filtered variants": "filtered_variants",
    "status": "status",
}

COVERAGE_FIELDS = [
    "mean_depth",
    "min_depth",
    "max_depth",
    "percent_above_100x",
    "percent_above_1000x",
    "percent_above_5000x",
]

COLUMNS = [
    "sample",
    "dataset",
    "contig",
    "raw_variants",
    "filtered_variants",
] + COVERAGE_FIELDS + ["status"]


def _sample_from_filename(path: Path, *suffixes: str) -> str:
    name = path.name
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def parse_qc_stats(path: Path) -> dict:
    """Parse a LOFREQ_FILTER ``qc_stats.txt`` file of ``Key: value`` lines."""
    record = {"sample": _sample_from_filename(path, ".qc_stats.txt", "_qc_stats.txt")}
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if ":" not in line:
                    continue
                key, _, value = line.partition(":")
                column = QC_FIELDS.get(key.strip().lower())
                if column:
                    record[column] = value.strip()
    except OSError as exc:
        print(f"Warning: could not read {path}: {exc}", file=sys.stderr)
    if not record.get("sample"):
        record["sample"] = _sample_from_filename(path, ".qc_stats.txt")
    return record


def parse_coverage_summary(path: Path) -> list[dict]:
    """Parse a COVERAGE_SUMMARIZE TSV (header row plus one row per sample)."""
    records = []
    try:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                sample = (row.get("sample") or "").strip()
                if not sample:
                    continue
                record = {"sample": sample}
                for field in COVERAGE_FIELDS:
                    if row.get(field) not in (None, ""):
                        record[field] = row[field].strip()
                records.append(record)
    except OSError as exc:
        print(f"Warning: could not read {path}: {exc}", file=sys.stderr)
    return records


def collect(results_dir: Path) -> list[dict]:
    """Walk ``results_dir`` and merge every QC artefact found into per-sample rows."""
    samples: dict[str, dict] = {}

    def slot(sample_id: str) -> dict:
        return samples.setdefault(sample_id, {"sample": sample_id})

    qc_files = sorted(results_dir.rglob("*qc_stats.txt"))
    for qc_file in qc_files:
        record = parse_qc_stats(qc_file)
        slot(record["sample"]).update({k: v for k, v in record.items() if v != ""})

    cov_files = sorted(results_dir.rglob("*coverage_summary.tsv"))
    for cov_file in cov_files:
        for record in parse_coverage_summary(cov_file):
            slot(record["sample"]).update(record)

    if not qc_files and not cov_files:
        print(
            f"Warning: no *qc_stats.txt or *coverage_summary.tsv files found under "
            f"{results_dir}; writing header-only summary.",
            file=sys.stderr,
        )

    return [samples[key] for key in sorted(samples)]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate a consolidated per-sample run summary from staged QC files."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--results-dir", required=True, help="Directory holding staged QC files")
    parser.add_argument("--outdir", required=True, help="Output directory for the summary TSV")
    parser.add_argument("--dataset", default="viral_analysis", help="Dataset label")

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.is_dir():
        print(f"Error: --results-dir {results_dir} is not a directory", file=sys.stderr)
        return 1

    outdir = Path(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    summary_path = outdir / "consolidated_sample_summary.tsv"

    rows = collect(results_dir)

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore", restval="NA"
        )
        writer.writeheader()
        for row in rows:
            row.setdefault("dataset", args.dataset)
            row["dataset"] = args.dataset
            writer.writerow(row)

    print(f"Run summary for {len(rows)} sample(s) written to {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
