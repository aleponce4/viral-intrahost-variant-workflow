#!/usr/bin/env bash
# Run the whole stock workflow for one stock, in three steps:
#
#   stage-a    de novo consensus, with nf-core/viralmetagenome (an external pipeline)
#   reference  turn that consensus into a per-stock reference (-entry BUILD_STOCK_REFERENCE)
#   variants   call variants against that reference (the main workflow)
#
# Each step is its own Nextflow run in its own directory under --outdir, so each keeps
# its own work directory and `-resume` works per step. Run `--help` for the options.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Stage A is pinned. Its output layout (consensus/seq/variant-calling/<sample>/) is
# read below, so a new revision needs that path checked before the pin moves.
VMG_REVISION="1.2.0"
STAGE_A_PARAMS="${REPO_ROOT}/assets/stage_a_viralmetagenome.params.yaml"

NXF_STAGE_A="${NXF_STAGE_A:-nextflow}"
NXF_STAGE_B="${NXF_STAGE_B:-nextflow}"

usage() {
    cat <<'EOF'
Usage: run_stock_workflow.sh --stock NAME --outdir DIR [options]

Required
  --stock NAME            Name of the stock. It becomes the contig name and the sample
                          name in Stage A.
  --outdir DIR            Where everything is written (one subdirectory per step).

Stage A (de novo consensus)
  --stock-reads R1 R2     The stock's paired FASTQ files.
  --reference-pool FASTA  Related genomes the assembly is scaffolded against.
  --consensus FILE        Use this consensus and skip Stage A.

Stock reference
  --lab-fasta FASTA       The lab reference genome.
  --lab-gff GFF3          The lab reference annotation (with CDS features).
  --primer-bed BED        Optional. Primer scheme in lab-reference coordinates. It is
                          moved onto the stock and the variant run switches to
                          --protocol amplicon.

Variant calling
  --samples CSV           Samplesheet (sample,fastq_1,fastq_2,treatment) of the
                          samples that belong to this stock.

Run control
  --step STEP             all (default), stage-a, reference or variants.
  --profile NAME          Nextflow profile for every step: docker, apptainer or singularity.
                          Default: docker, or the value of STOCK_WORKFLOW_PROFILE.
  --config FILE           Extra Nextflow config for every step, for example resource limits.
  --stage-a-extra "ARGS"  Extra arguments for the Stage A run.
  --variants-extra "ARGS" Extra arguments for the variant run.
  --dry-run               Print the commands. Run nothing and write nothing.
  -h, --help              This text.

Environment
  STOCK_WORKFLOW_PROFILE     Default for --profile, so a site sets it once.
  NXF_STAGE_A, NXF_STAGE_B   The Nextflow binary for Stage A and for the other two steps.
                             Default: nextflow. Use these when the two need different versions.

Steps that finished are reused, because every run passes -resume. Re-running the same
command after a failure picks up where it stopped.
EOF
}

die() { echo "ERROR: $*" >&2; exit 2; }

abs() {
    case "$1" in
        /*) printf '%s\n' "$1" ;;
        *)  printf '%s/%s\n' "$PWD" "${1#./}" ;;
    esac
}

STOCK="" OUTDIR="" READS_1="" READS_2="" POOL="" CONSENSUS_ARG=""
LAB_FASTA="" LAB_GFF="" PRIMER_BED="" SAMPLES=""
STEP="all" PROFILE="${STOCK_WORKFLOW_PROFILE:-docker}" CONFIG="" STAGE_A_EXTRA="" VARIANTS_EXTRA="" DRY_RUN=0

while [ $# -gt 0 ]; do
    case "$1" in
        --stock)           [ $# -ge 2 ] || die "--stock needs a value";           STOCK="$2"; shift 2 ;;
        --outdir)          [ $# -ge 2 ] || die "--outdir needs a value";          OUTDIR="$2"; shift 2 ;;
        --stock-reads)     [ $# -ge 3 ] || die "--stock-reads needs two files";   READS_1="$2"; READS_2="$3"; shift 3 ;;
        --reference-pool)  [ $# -ge 2 ] || die "--reference-pool needs a value";  POOL="$2"; shift 2 ;;
        --consensus)       [ $# -ge 2 ] || die "--consensus needs a value";       CONSENSUS_ARG="$2"; shift 2 ;;
        --lab-fasta)       [ $# -ge 2 ] || die "--lab-fasta needs a value";       LAB_FASTA="$2"; shift 2 ;;
        --lab-gff)         [ $# -ge 2 ] || die "--lab-gff needs a value";         LAB_GFF="$2"; shift 2 ;;
        --primer-bed)      [ $# -ge 2 ] || die "--primer-bed needs a value";      PRIMER_BED="$2"; shift 2 ;;
        --samples)         [ $# -ge 2 ] || die "--samples needs a value";         SAMPLES="$2"; shift 2 ;;
        --step)            [ $# -ge 2 ] || die "--step needs a value";            STEP="$2"; shift 2 ;;
        --profile)         [ $# -ge 2 ] || die "--profile needs a value";         PROFILE="$2"; shift 2 ;;
        --config)          [ $# -ge 2 ] || die "--config needs a value";          CONFIG="$2"; shift 2 ;;
        --stage-a-extra)   [ $# -ge 2 ] || die "--stage-a-extra needs a value";   STAGE_A_EXTRA="$2"; shift 2 ;;
        --variants-extra)  [ $# -ge 2 ] || die "--variants-extra needs a value";  VARIANTS_EXTRA="$2"; shift 2 ;;
        --dry-run)         DRY_RUN=1; shift ;;
        -h|--help)         usage; exit 0 ;;
        *)                 usage >&2; die "unknown option: $1" ;;
    esac
done

[ -n "$STOCK" ]  || { usage >&2; die "--stock is required"; }
[ -n "$OUTDIR" ] || { usage >&2; die "--outdir is required"; }
case "$STOCK" in
    *[!A-Za-z0-9._-]*) die "--stock may hold only letters, digits, '.', '_' and '-'. It becomes a contig name and a sample name." ;;
esac
case "$STEP" in all|stage-a|reference|variants) ;; *) die "--step must be all, stage-a, reference or variants" ;; esac

# --- which steps run ---------------------------------------------------------
RUN_A=0; RUN_REF=0; RUN_VAR=0
case "$STEP" in
    all)       RUN_A=1; RUN_REF=1; RUN_VAR=1 ;;
    stage-a)   RUN_A=1 ;;
    reference) RUN_REF=1 ;;
    variants)  RUN_VAR=1 ;;
esac
# A supplied consensus replaces Stage A.
if [ -n "$CONSENSUS_ARG" ] && [ "$RUN_A" -eq 1 ]; then
    if [ "$STEP" = "stage-a" ]; then die "--consensus replaces Stage A, so --step stage-a has nothing to do"; fi
    RUN_A=0
    echo "Note: --consensus given, so Stage A is skipped."
fi

# --- inputs for the steps that will run ---------------------------------------
need_file() { [ -n "$2" ] || die "$1 is required for this step"; [ -f "$2" ] || die "$1: file not found: $2"; }

if [ "$RUN_A" -eq 1 ]; then
    [ -n "$READS_1" ] && [ -n "$READS_2" ] || die "--stock-reads R1 R2 is required for Stage A"
    need_file "--stock-reads (R1)" "$READS_1"
    need_file "--stock-reads (R2)" "$READS_2"
    need_file "--reference-pool" "$POOL"
fi
if [ "$RUN_REF" -eq 1 ]; then
    need_file "--lab-fasta" "$LAB_FASTA"
    need_file "--lab-gff" "$LAB_GFF"
    [ -z "$PRIMER_BED" ] || need_file "--primer-bed" "$PRIMER_BED"
    [ -z "$CONSENSUS_ARG" ] || need_file "--consensus" "$CONSENSUS_ARG"
fi
if [ "$RUN_VAR" -eq 1 ]; then
    need_file "--samples" "$SAMPLES"
fi
[ -z "$CONFIG" ] || need_file "--config" "$CONFIG"
[ -f "$STAGE_A_PARAMS" ] || die "missing Stage A params file: $STAGE_A_PARAMS"

# --- paths -------------------------------------------------------------------
OUTDIR="$(abs "$OUTDIR")"
[ -z "$READS_1" ] || { READS_1="$(abs "$READS_1")"; READS_2="$(abs "$READS_2")"; }
[ -z "$POOL" ] || POOL="$(abs "$POOL")"
[ -z "$CONSENSUS_ARG" ] || CONSENSUS_ARG="$(abs "$CONSENSUS_ARG")"
[ -z "$LAB_FASTA" ] || LAB_FASTA="$(abs "$LAB_FASTA")"
[ -z "$LAB_GFF" ] || LAB_GFF="$(abs "$LAB_GFF")"
[ -z "$PRIMER_BED" ] || PRIMER_BED="$(abs "$PRIMER_BED")"
[ -z "$SAMPLES" ] || SAMPLES="$(abs "$SAMPLES")"
[ -z "$CONFIG" ] || CONFIG="$(abs "$CONFIG")"

A_DIR="${OUTDIR}/stage_a"
REF_DIR="${OUTDIR}/reference"
VAR_DIR="${OUTDIR}/variants"
A_OUT="${A_DIR}/results"
STOCK_REF="${REF_DIR}/results/StockReference/${STOCK}"
REF_FASTA="${STOCK_REF}/${STOCK}.fasta"
REF_GFF="${STOCK_REF}/${STOCK}.gff3"
REF_LIFTOVER="${STOCK_REF}/qc/${STOCK}.liftover.tsv"
REF_PRIMERS="${STOCK_REF}/${STOCK}.primers.bed"
LOG="${OUTDIR}/stock_run.log"

CONFIG_ARGS=()
[ -z "$CONFIG" ] || CONFIG_ARGS=(-c "$CONFIG")

# --- helpers -----------------------------------------------------------------
quote_cmd() { printf '%q ' "$@"; printf '\n'; }

log() {
    [ "$DRY_RUN" -eq 1 ] && return 0
    mkdir -p "$OUTDIR"
    printf '%s\n' "$*" >> "$LOG"
}

nxf_version() {
    "$1" -version 2>/dev/null | grep -i 'version' | head -n 1 | sed 's/^ *//' || true
}

# The one consensus Stage A wrote for this stock. Stage A names its output by the sample
# name, which is --stock. More than one file means more than one assembled cluster.
find_consensus() {
    local dir="${A_OUT}/consensus/seq/variant-calling/${STOCK}"
    local hits=()
    local f
    if [ -d "$dir" ]; then
        for f in "$dir"/*.consensus.fasta; do
            [ -e "$f" ] && hits+=("$f")
        done
    fi
    if [ "${#hits[@]}" -eq 0 ]; then
        echo "ERROR: no consensus found in ${dir}" >&2
        echo "       Run Stage A first (--step stage-a), or pass --consensus FILE." >&2
        return 1
    fi
    if [ "${#hits[@]}" -gt 1 ]; then
        echo "ERROR: more than one consensus in ${dir}:" >&2
        printf '         %s\n' "${hits[@]}" >&2
        echo "       Pick the right one and pass it with --consensus FILE." >&2
        return 1
    fi
    printf '%s\n' "${hits[0]}"
}

# run_step NAME DIR CMD...  : run CMD inside DIR, or print it on a dry run.
run_step() {
    local name="$1" dir="$2"
    shift 2
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "[dry-run] ${name}: cd ${dir} && $(quote_cmd "$@")"
        return 0
    fi
    mkdir -p "$dir"
    log "== ${name} =="
    log "cd ${dir} && $(quote_cmd "$@")"
    echo "== ${name} =="
    ( cd "$dir" && "$@" )
}

if [ "$DRY_RUN" -eq 0 ]; then
    log ""
    log "run started $(date -u +%Y-%m-%dT%H:%M:%SZ), stock ${STOCK}, step ${STEP}"
    log "pipeline commit: $(git -C "$REPO_ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
fi

# --- Stage A -----------------------------------------------------------------
if [ "$RUN_A" -eq 1 ]; then
    if [ "$DRY_RUN" -eq 1 ]; then
        echo "[dry-run] would write ${A_DIR}/samplesheet.csv:"
        printf '            sample,fastq_1,fastq_2\n            %s,%s,%s\n' "$STOCK" "$READS_1" "$READS_2"
    else
        mkdir -p "$A_DIR"
        printf 'sample,fastq_1,fastq_2\n%s,%s,%s\n' "$STOCK" "$READS_1" "$READS_2" > "${A_DIR}/samplesheet.csv"
    fi

    stage_a_cmd=("$NXF_STAGE_A" run nf-core/viralmetagenome -r "$VMG_REVISION" -profile "$PROFILE")
    stage_a_cmd+=("${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"}")
    stage_a_cmd+=(-params-file "$STAGE_A_PARAMS" --input "${A_DIR}/samplesheet.csv"
                  --reference_pool "$POOL" --outdir "$A_OUT" -resume)
    if [ -n "$STAGE_A_EXTRA" ]; then
        read -r -a extra_a <<< "$STAGE_A_EXTRA"
        stage_a_cmd+=("${extra_a[@]}")
    fi
    if [ "$DRY_RUN" -eq 0 ]; then
        log "Stage A nextflow: $(nxf_version "$NXF_STAGE_A")"
    fi
    run_step "Stage A: de novo consensus" "$A_DIR" "${stage_a_cmd[@]}"
fi

# --- reference ---------------------------------------------------------------
if [ "$RUN_REF" -eq 1 ]; then
    if [ -n "$CONSENSUS_ARG" ]; then
        CONSENSUS="$CONSENSUS_ARG"
    elif [ "$DRY_RUN" -eq 1 ] && [ "$RUN_A" -eq 1 ]; then
        CONSENSUS="CONSENSUS_FROM_STAGE_A"
    else
        CONSENSUS="$(find_consensus)" || exit 2
    fi

    ref_cmd=("$NXF_STAGE_B" run "$REPO_ROOT" -entry BUILD_STOCK_REFERENCE -profile "$PROFILE")
    ref_cmd+=("${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"}")
    ref_cmd+=(--stock_consensus "$CONSENSUS" --stock_name "$STOCK" --fasta "$LAB_FASTA" --gff "$LAB_GFF")
    [ -z "$PRIMER_BED" ] || ref_cmd+=(--stock_primer_bed "$PRIMER_BED")
    ref_cmd+=(--outdir "${REF_DIR}/results" -resume)
    if [ "$DRY_RUN" -eq 0 ]; then
        log "Stage B nextflow: $(nxf_version "$NXF_STAGE_B")"
    fi
    run_step "Stock reference" "$REF_DIR" "${ref_cmd[@]}"
fi

# --- variants ----------------------------------------------------------------
if [ "$RUN_VAR" -eq 1 ]; then
    if [ "$DRY_RUN" -eq 0 ]; then
        for f in "$REF_FASTA" "$REF_GFF" "$REF_LIFTOVER"; do
            [ -f "$f" ] || die "missing ${f}. Run the reference step first (--step reference)."
        done
        if [ -n "$PRIMER_BED" ] && [ ! -f "$REF_PRIMERS" ]; then
            die "missing ${REF_PRIMERS}. The reference was built without --primer-bed. Re-run --step reference with it."
        fi
    fi

    var_cmd=("$NXF_STAGE_B" run "$REPO_ROOT" -profile "$PROFILE")
    var_cmd+=("${CONFIG_ARGS[@]+"${CONFIG_ARGS[@]}"}")
    var_cmd+=(--input "$SAMPLES" --fasta "$REF_FASTA" --gff "$REF_GFF"
              --viral_contig "$STOCK" --dataset "$STOCK" --liftover_tsv "$REF_LIFTOVER")
    [ -z "$PRIMER_BED" ] || var_cmd+=(--protocol amplicon --primer_bed "$REF_PRIMERS")
    var_cmd+=(--outdir "${VAR_DIR}/results" -resume)
    if [ -n "$VARIANTS_EXTRA" ]; then
        read -r -a extra_v <<< "$VARIANTS_EXTRA"
        var_cmd+=("${extra_v[@]}")
    fi
    run_step "Variant calling" "$VAR_DIR" "${var_cmd[@]}"
fi

# --- summary -----------------------------------------------------------------
if [ "$DRY_RUN" -eq 1 ]; then
    echo "[dry-run] nothing was run."
else
    echo
    echo "Done. Outputs under ${OUTDIR}:"
    [ "$RUN_A" -eq 0 ]   || echo "  Stage A      ${A_OUT}"
    [ "$RUN_REF" -eq 0 ] || echo "  Reference    ${STOCK_REF}"
    [ "$RUN_VAR" -eq 0 ] || echo "  Variants     ${VAR_DIR}/results"
    echo "  Commands     ${LOG}"
fi
