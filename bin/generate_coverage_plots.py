#!/usr/bin/env python3
"""
generate_coverage_plots.py — Real depth/coverage figures from samtools depth TSVs.

Reads the per-base depth tables staged by the COVERAGE_QC subworkflow
(``<contig>\\t<position>\\t<depth>``) and produces:

  * ``<dataset>_coverage_summary.png`` — per-sample depth profile along the genome
    (log-scaled), with the 100x/1,000x/5,000x reference lines used by the reporting tiers,
    plus a companion mean-depth bar panel.

Per-sample statistics are computed with ``summarize_coverage.summarize_depth`` so the
numbers in the figure are the same ones written to the coverage summary TSVs.

If no depth records are staged, an explicitly labelled "no coverage data" figure is
produced. No synthetic or placeholder depth values are ever plotted.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from summarize_coverage import summarize_depth  # noqa: E402

__version__ = "2.0.0"

DEPTH_SUFFIXES = (".depth.tsv", "_depth.tsv", ".depth.txt", ".tsv", ".txt")
TIER_LINES = ((100, "#8c8c8c", "100x"), (1000, "#e06d53", "1,000x"), (5000, "#2b5c8f", "5,000x"))


def sample_id_from_path(path: Path) -> str:
    name = path.name
    for suffix in DEPTH_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def read_depth_profile(path: Path) -> tuple[list[int], list[int]]:
    """Return (positions, depths) from a samtools depth TSV."""
    positions: list[int] = []
    depths: list[int] = []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                try:
                    positions.append(int(parts[1]))
                    depths.append(int(parts[2]))
                except ValueError:
                    continue
    except OSError as exc:
        print(f"Warning: could not read {path}: {exc}", file=sys.stderr)
    return positions, depths


def rolling_mean(values: list[int], window: int) -> list[float]:
    """Simple centred rolling mean; returns the input as floats when window <= 1."""
    if window <= 1 or len(values) <= window:
        return [float(v) for v in values]
    out: list[float] = []
    running = 0.0
    for i, value in enumerate(values):
        running += value
        if i >= window:
            running -= values[i - window]
        out.append(running / min(i + 1, window))
    return out


def find_depth_files(coverage_dir: Path) -> list[Path]:
    """Candidate per-base depth tables, excluding already-summarised coverage tables."""
    return [
        p
        for p in sorted(coverage_dir.rglob("*"))
        if p.is_file()
        and p.suffix.lower() in (".tsv", ".txt")
        and not p.name.endswith("coverage_summary.tsv")
    ]


def plot_coverage_figure(profiles: dict, stats: dict, outdir: Path, dataset: str, window: int) -> Path:
    """Render the coverage figure. Imported lazily so --help/--version stay dependency-free."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_path = outdir / f"{dataset}_coverage_summary.png"

    if not profiles:
        fig, ax = plt.subplots(figsize=(9, 4), dpi=200)
        ax.text(
            0.5,
            0.5,
            "No coverage data available for this run",
            ha="center",
            va="center",
            fontsize=13,
        )
        ax.set_axis_off()
        fig.suptitle(f"Depth & Coverage QC ({dataset})", fontsize=13, fontweight="bold")
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        return out_path

    fig, (ax_profile, ax_mean) = plt.subplots(
        2, 1, figsize=(10, 8), dpi=200, gridspec_kw={"height_ratios": [2, 1]}
    )

    for sample in sorted(profiles):
        positions, depths = profiles[sample]
        ax_profile.plot(
            positions,
            rolling_mean(depths, window),
            linewidth=1.2,
            alpha=0.9,
            label=sample,
        )

    for level, colour, label in TIER_LINES:
        ax_profile.axhline(level, color=colour, linestyle="--", linewidth=0.9, alpha=0.7)
        ax_profile.annotate(
            label,
            xy=(1.002, level),
            xycoords=("axes fraction", "data"),
            color=colour,
            fontsize=8,
            va="center",
        )

    ax_profile.set_yscale("symlog", linthresh=10)
    ax_profile.set_xlabel("Genomic Position (bp)", fontsize=11)
    ax_profile.set_ylabel(f"Depth (x, {window}-bp rolling mean)", fontsize=11)
    ax_profile.set_title(
        f"Per-base Sequencing Depth ({dataset})", fontsize=13, fontweight="bold", pad=10
    )
    ax_profile.grid(True, which="major", alpha=0.25)
    ax_profile.legend(fontsize=8, ncol=2, loc="upper right", framealpha=0.9)

    samples = sorted(stats)
    means = [stats[s]["mean_depth"] for s in samples]
    ax_mean.bar(samples, means, color="#2b5c8f", alpha=0.85, edgecolor="black", linewidth=0.5)
    ax_mean.set_ylabel("Mean Depth (x)", fontsize=11)
    ax_mean.set_title("Mean Depth per Sample", fontsize=12, pad=8)
    ax_mean.grid(True, axis="y", alpha=0.25)
    ax_mean.tick_params(axis="x", labelrotation=45)
    for label in ax_mean.get_xticklabels():
        label.set_horizontalalignment("right")

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate depth/coverage figures from samtools depth TSVs."
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--coverage-dir", required=True, help="Directory of samtools depth TSVs")
    parser.add_argument("--outdir", required=True, help="Output directory for plots")
    parser.add_argument("--dataset", default="viral_analysis", help="Dataset label")
    parser.add_argument(
        "--windowsize", type=int, default=50, help="Rolling-mean window for the depth profile"
    )

    args = parser.parse_args()

    coverage_dir = Path(args.coverage_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not coverage_dir.is_dir():
        print(f"Error: --coverage-dir {coverage_dir} is not a directory", file=sys.stderr)
        return 1

    profiles: dict[str, tuple[list[int], list[int]]] = {}
    stats: dict[str, dict] = {}

    for depth_file in find_depth_files(coverage_dir):
        sample = sample_id_from_path(depth_file)
        positions, depths = read_depth_profile(depth_file)
        if not depths:
            print(f"Note: {depth_file} contained no depth records; skipping.", file=sys.stderr)
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

    for sample in sorted(stats):
        s = stats[sample]
        print(
            f"{sample}: mean={s['mean_depth']}x min={s['min_depth']}x max={s['max_depth']}x "
            f">=100x={s['percent_above_100x']}% >=1000x={s['percent_above_1000x']}% "
            f">=5000x={s['percent_above_5000x']}%"
        )
    print(f"Coverage figure for {len(profiles)} sample(s) written to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
