#!/usr/bin/env python3
"""
plot_coverage.py — Depth/coverage figure generator (``--output-dir`` flavoured CLI).

Thin wrapper around :mod:`generate_coverage_plots` so that callers using the
``--coverage-dir`` / ``--output-dir`` argument spelling get the same real figure. All
plotting logic lives in ``generate_coverage_plots.py``; there is exactly one implementation.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from generate_coverage_plots import (  # noqa: E402
    find_depth_files,
    plot_coverage_figure,
    read_depth_profile,
    sample_id_from_path,
)
from summarize_coverage import summarize_depth  # noqa: E402

__version__ = "2.0.0"


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate coverage distribution plots.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--coverage-dir", required=True, help="Input directory containing samtools depth TSVs"
    )
    parser.add_argument("--output-dir", required=True, help="Output directory for plots")
    parser.add_argument("--dataset", default="viral_analysis", help="Dataset label")
    parser.add_argument(
        "--windowsize", type=int, default=50, help="Rolling-mean window for the depth profile"
    )

    args = parser.parse_args()

    coverage_dir = Path(args.coverage_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not coverage_dir.is_dir():
        print(f"Error: --coverage-dir {coverage_dir} is not a directory", file=sys.stderr)
        return 1

    profiles = {}
    stats = {}
    for depth_file in find_depth_files(coverage_dir):
        sample = sample_id_from_path(depth_file)
        positions, depths = read_depth_profile(depth_file)
        if not depths:
            continue
        profiles[sample] = (positions, depths)
        stats[sample] = summarize_depth(str(depth_file), sample)

    if not profiles:
        print(
            f"Warning: no usable depth records found under {coverage_dir}; "
            "emitting an explicit 'no coverage data' figure.",
            file=sys.stderr,
        )

    out_path = plot_coverage_figure(profiles, stats, outdir, args.dataset, args.windowsize)
    print(f"Coverage figure for {len(profiles)} sample(s) written to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
