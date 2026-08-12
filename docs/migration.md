# From Bash scripts to Nextflow DSL2

This pipeline did not start as a workflow-manager project. It started as a directory of
Bash scripts that grew around one real analysis problem: calling low-frequency intra-host
variants (iSNVs) in alphavirus RNA-seq data deeply enough to say something about selection.
That original code is preserved under [`legacy/`](../legacy) so the lineage is visible. This
note explains what it did, why it was rewritten, and what actually changed.

## What the legacy pipeline did

`legacy/run_full_pipeline.sh` was the entry point. It sourced a single
`Variant_discovery_pipeline/config.sh` holding every threshold as a shell variable, then ran
six phases in sequence, each gated by a `RUN_PHASE_*` environment flag:

1. **Extract viral BAMs** — `samtools view`/`sort`/`index` to pull one viral contig out of
   host-dominated alignments.
2. **Variant calling** — `lofreq call-parallel` + `lofreq filter` for statistically modelled
   calls, and `samtools mpileup | ivar variants` + `ivar consensus` for frequency-based
   calls, so the two callers could be compared.
3. **Annotation** — `ivar_variants_to_vcf.py` to turn iVar TSVs into valid VCFs, then
   `bcftools csq --local-csq` for consequence prediction on a haploid viral CDS.
4. **Coverage QC** — `samtools depth` plus an awk/Python summary of mean depth and the
   fraction of the genome above 100×/1,000×/5,000×.
5. **Selection analysis** — SNPGenie (Perl) per sample, then Python and R helpers for
   πN/πS deltas, Kruskal–Wallis tests, and limma models.
6. **Haplotype reconstruction** — CliqueSNV and VILOCA/ShoRAH, with helper scripts to build
   haplotype tables and PopART/Datamonkey inputs.

It worked, and it produced real results. But the operational shape of it was the problem.

## Why it was rewritten

Every complaint below is something the shell version actually cost us, not a theoretical
concern:

- **Parallelism was hand-rolled.** Each script had its own `MAX_JOBS` loop backgrounding
  jobs with `&` and reaping them with `wait`. Tuning meant editing `THREADS` and `MAX_JOBS`
  in `config.sh` per machine, and a crashed background job could be silently swallowed
  because `set -e` had to stay off (a non-zero background exit would kill the whole run).
- **Resume was a `[ -f ... ]` check.** Restartability was implemented as "skip the sample
  if the output file exists," controlled by a `FORCE_RECALL` flag. A half-written output
  from an interrupted run looked exactly like a finished one.
- **Environments were the developer's, not the pipeline's.** Six Conda environments
  (`lofreq-env`, `ivar_env`, `annotation-env`, `env_cliquesnv`, `env_viloca`, plus R
  packages) were activated by name mid-script. Nothing recorded which versions ran, and
  moving to a cluster meant recreating all six by hand.
- **Paths were absolute and machine-specific.** The scripts carried hard-coded WSL and
  `/mnt/d` paths, and the last phase copied results across the WSL/Windows boundary. That
  is not portable and it leaks the author's filesystem into the repository.
- **No test surface.** There was no way to exercise the pipeline without a full real
  dataset, so a change to an annotation helper could only be validated by rerunning hours
  of variant calling.

Nextflow DSL2 addresses all five directly: the executor owns scheduling, the work-directory
hash owns resume, containers own the environment, channels own the paths, and `nf-test`
gives per-module assertions.

## What changed structurally

The one-to-one mapping from legacy script to DSL2 component is in the README's migration
table. The architectural differences are these:

| Concern | Legacy | Now |
|---|---|---|
| Orchestration | Sequential phases, `RUN_PHASE_*` env flags | DAG of subworkflows, `--run_*` params |
| Scheduling | `&` / `wait` with `MAX_JOBS` | Executor (`local`, `slurm`, `awsbatch`) |
| Restart | `[ -f output ]` + `FORCE_RECALL` | Content-hashed task cache, `-resume` |
| Software | 6 named Conda envs | Pinned Biocontainers per process, lint-enforced |
| Config | One `config.sh` of shell vars | `params` + `conf/{base,modules,containers,test}.config` |
| Output paths | Hard-coded absolute paths + a copy step | `publishDir` with `--outdir` |
| Provenance | None | `versions.yml` per process → MultiQC + executive report |
| Validation | Manual eyeballing | `nf-test` suite, schema validation, CI lint |
| Helpers | Loose scripts sourcing shell globals | `bin/` executables with strict `argparse`/`optparse` CLIs |

Three things deliberately did **not** survive: the WSL→Windows result-sync step (replaced by
`publishDir`), the skip-if-exists logic (replaced by Nextflow's cache), and the
`THREADS`/`MAX_JOBS` variables (replaced by per-process `cpus` with `check_max()` capping).
Legacy caller thresholds were carried over as `params` so behaviour stayed comparable across
the cut, and were later revised on their own merits (see `CHANGELOG.md` 1.3.0 and
`SCIENTIFIC_AUDIT.md`).

## What was gained

- **Reproducibility.** A run is defined by a samplesheet, a reference, a params set, and a
  set of container digests. Nothing depends on what is installed on the host.
- **Portability.** The same code runs as `-profile test` on a laptop, `-profile slurm` on
  an HPC cluster, and `-profile awsbatch` in the cloud, with no edits to process logic.
- **Cheap iteration.** Synthetic fixtures under `tests/data/` plus `-stub-run` mean a
  structural change can be validated in seconds instead of hours, and CI does it on every
  push.
- **Auditable provenance.** Every process emits its tool versions; MultiQC and the
  executive HTML report collect them alongside the results, so a figure can be traced back
  to the exact software that produced it.
- **A reviewable unit of change.** Modules are small, individually snapshot-tested files.
  Reviewing "did the LoFreq filter change?" is now reading twenty lines, not diffing a
  three-hundred-line script that also does logging, parallelism, and path management.

## What the legacy tree is still good for

`legacy/` is reference material, not a supported code path. It documents the parameter
choices the DSL2 version inherited and the exploratory notebook analysis
(`legacy/Variant_discovery_pipeline/Scripts/analysis.ipynb`) that motivated several of the
reporting tiers. It is excluded from language statistics via `.gitattributes` and contains
no study data — the notebook is committed with outputs cleared, and its paths are relative
or placeholders.
