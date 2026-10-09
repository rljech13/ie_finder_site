#!/bin/bash
# FASTA in, genome-coordinate GFF3 and GenBank out.
# This directory is a reduced copy of the MGE_finder discovery steps and runs on its own.
#
#   ./run.sh strain.fasta outdir
#   ./run.sh genomes_dir outdir
#   ANNOTATE_ALL=1 ./run.sh strain.fasta outdir
#   KEEP_WORK=1 ./run.sh strain.fasta outdir
#   SEARCH_PARAMS=strict.yaml ./run.sh strain.fasta outdir
#
# The first run creates the conda environment named in envs/IE_finder.yaml
# and puts it on PATH. Later runs reuse it. Search thresholds are in
# search_params.yaml. SEARCH_PARAMS names a file whose keys override them.
#
# Extra arguments after the outdir are passed to Snakemake.

set -euo pipefail

CALL_DIR="$(pwd)"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Relative paths are taken from the directory where the script was started,
# not from this repository after the cd below.
abspath() {
  case "$1" in
    /*) realpath -m "$1" ;;
    *) realpath -m "$CALL_DIR/$1" ;;
  esac
}

cd "$SCRIPT_DIR"

if [[ $# -lt 1 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  echo "Usage: ./run.sh genome.fasta|genomes_dir [outdir] [snakemake args...]" >&2
  exit 2
fi

INPUT="$(abspath "$1")"
shift
if [[ $# -ge 1 && "$1" != -* ]]; then
  OUTDIR="$(abspath "$1")"
  shift
else
  OUTDIR="$(realpath -m "$SCRIPT_DIR/outdir")"
fi
mkdir -p "$OUTDIR"

if [[ ! -e "$INPUT" ]]; then
  echo "Input not found: $INPUT" >&2
  exit 1
fi

BASE_PARAMS="$SCRIPT_DIR/search_params.yaml"
USER_PARAMS=""
if [[ -n "${SEARCH_PARAMS:-}" ]]; then
  USER_PARAMS="$(abspath "$SEARCH_PARAMS")"
  if [[ ! -f "$USER_PARAMS" ]]; then
    echo "Search parameters not found: $USER_PARAMS" >&2
    exit 1
  fi
fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ie_finder_site.XXXXXX")"
# The values this run actually used: search_params.yaml plus SEARCH_PARAMS.
PARAMS_USED="$WORK/search_params.yaml"
UPSTREAM_COMMIT="$(awk -F': ' '/^commit:/{print $2; exit}' "$SCRIPT_DIR/UPSTREAM")"

write_reports() {
  shopt -s nullglob
  for sample_dir in "$RESULTS"/*/; do
    sample="$(basename "$sample_dir")"
    if [[ ! -f "$sample_dir/integrase_hits_summary.tsv" && ! -f "$sample_dir/ie_filter_audit.tsv" ]]; then
      continue
    fi
    python3 "$SCRIPT_DIR/scripts/write_report.py" \
      --sample-dir "$sample_dir" \
      --sample "$sample" \
      --mge-finder "$SCRIPT_DIR" \
      --upstream-commit "$UPSTREAM_COMMIT" \
      --params "${PARAMS_USED:-}" \
      --out "$OUTDIR/${sample}.ie.report.txt" \
      || echo "Could not write a report for $sample" >&2
  done
  shopt -u nullglob
}

cleanup() {
  local status=$?
  if [[ $status -ne 0 && -d "${RESULTS:-}" ]]; then
    write_reports || true
  fi
  if [[ $status -eq 0 && "${KEEP_WORK:-0}" != "1" ]]; then
    rm -rf "$WORK"
  elif [[ -d "$WORK" ]]; then
    echo "Work left at $WORK" >&2
  fi
  exit "$status"
}
trap cleanup EXIT

GENOMES="$WORK/genomes"
RESULTS="$WORK/results"
mkdir -p "$GENOMES" "$RESULTS"

link_fasta() {
  local src="$1"
  local base sample dest
  base="$(basename "$src")"
  sample="${base%.*}"
  if [[ -z "$sample" || "$sample" == .* ]]; then
    echo "Could not derive a sample name from $src" >&2
    exit 1
  fi
  dest="$GENOMES/${sample}.fna"
  if [[ -e "$dest" ]]; then
    echo "Two inputs share the sample name $sample" >&2
    exit 1
  fi
  ln -s "$src" "$dest"
}

shopt -s nullglob
if [[ -d "$INPUT" ]]; then
  files=("$INPUT"/*.fna "$INPUT"/*.fa "$INPUT"/*.fasta "$INPUT"/*.fsa)
  if [[ ${#files[@]} -eq 0 ]]; then
    echo "No FASTA files in $INPUT" >&2
    exit 1
  fi
  for src in "${files[@]}"; do
    link_fasta "$src"
  done
else
  if [[ ! -s "$INPUT" ]]; then
    echo "FASTA is empty: $INPUT" >&2
    exit 1
  fi
  link_fasta "$INPUT"
fi
shopt -u nullglob

find_conda() {
  local candidate
  if [[ -n "${CONDA_EXE:-}" && -x "$CONDA_EXE" ]]; then
    echo "$CONDA_EXE"
    return 0
  fi
  candidate="$(command -v conda 2>/dev/null || true)"
  if [[ -n "$candidate" && -x "$candidate" ]]; then
    echo "$candidate"
    return 0
  fi
  for candidate in \
    "$HOME/miniforge3/bin/conda" \
    "$HOME/mambaforge/bin/conda" \
    "$HOME/miniconda3/bin/conda" \
    "$HOME/anaconda3/bin/conda" \
    "/opt/miniforge3/bin/conda" \
    "/opt/conda/bin/conda"
  do
    if [[ -x "$candidate" ]]; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

install_miniforge() {
  local dest arch url installer
  dest="${MINIFORGE_ROOT:-$HOME/miniforge3}"
  arch="$(uname -m)"
  case "$arch" in
    x86_64|aarch64) ;;
    *)
      echo "conda is not installed, and Miniforge has no build for $arch" >&2
      exit 1
      ;;
  esac
  if [[ -x "$dest/bin/conda" ]]; then
    echo "$dest/bin/conda"
    return 0
  fi
  url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname -s)-${arch}.sh"
  installer="$(mktemp "${TMPDIR:-/tmp}/Miniforge3.XXXXXX.sh")"
  echo "Installing Miniforge to $dest" >&2
  curl -fL --retry 3 -o "$installer" "$url"
  bash "$installer" -b -p "$dest"
  rm -f "$installer"
  echo "$dest/bin/conda"
}

ensure_runtime() {
  local conda_bin base name prefix
  if ! conda_bin="$(find_conda)"; then
    conda_bin="$(install_miniforge)"
  fi
  name="$(awk '/^name:/{print $2; exit}' "$SCRIPT_DIR/envs/IE_finder.yaml")"
  base="$("$conda_bin" info --base)"
  prefix="$base/envs/$name"
  if [[ ! -x "$prefix/bin/snakemake" || ! -x "$prefix/bin/prodigal" || ! -x "$prefix/bin/hmmscan" || ! -x "$prefix/bin/aragorn" || ! -x "$prefix/bin/blastn" ]]; then
    echo "Building conda environment $name from envs/IE_finder.yaml" >&2
    if [[ -d "$prefix" ]]; then
      "$conda_bin" env update -p "$prefix" -f "$SCRIPT_DIR/envs/IE_finder.yaml"
    else
      "$conda_bin" env create -f "$SCRIPT_DIR/envs/IE_finder.yaml"
    fi
  fi
  if [[ ! -x "$prefix/bin/snakemake" || ! -x "$prefix/bin/python3" ]]; then
    echo "Environment $name has no snakemake after creation: $prefix" >&2
    exit 1
  fi
  export PATH="$prefix/bin:${PATH:-}"
  hash -r
}

ensure_runtime

CONFIG="$WORK/config.yaml"
python3 - "$SCRIPT_DIR/ie_finder_config.yaml" "$CONFIG" "$SCRIPT_DIR" "$GENOMES" "$RESULTS" "${ANNOTATE_ALL:-0}" \
  "$BASE_PARAMS" "$USER_PARAMS" "$PARAMS_USED" << 'PY'
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("PyYAML is required")

src, dst, root, genomes, results, annotate_all, base_params, user_params, params_used = sys.argv[1:10]
root = Path(root)
sys.path.insert(0, str(root / "scripts"))
from v3ps_filters import read_params_file, thresholds_from_params, thresholds_to_params

# Each file is checked on its own so an error names the file it came from.
try:
    params = read_params_file(base_params)
    thresholds_from_params(params, base_params)
    if user_params:
        override = read_params_file(user_params)
        thresholds_from_params(override, user_params)
        params.update(override)
    thresholds = thresholds_from_params(params)
except (OSError, ValueError, yaml.YAMLError) as exc:
    sys.exit(f"Search parameters: {exc}")
Path(params_used).write_text(yaml.safe_dump(thresholds_to_params(thresholds), sort_keys=False))

cfg = yaml.safe_load(Path(src).read_text()) or {}
cfg.setdefault("paths", {})
cfg["paths"].pop("config_file", None)
cfg["paths"]["search_params"] = params_used
cfg["paths"]["genomes_dir"] = genomes
cfg["paths"]["results_dir"] = results
cfg.setdefault("execution", {})
cfg["execution"]["conda_env"] = str(root / "envs" / "IE_finder.yaml")
cfg["annotate"] = {"all_candidates": annotate_all == "1"}
cfg["pfam_profiles"] = [
    str(root / "pfam" / "PF00589.hmm"),
    str(root / "pfam" / "PF22022.hmm"),
]
Path(dst).write_text(yaml.safe_dump(cfg, sort_keys=False))
PY

SM="${SNAKEMAKE:-snakemake}"
if ! command -v "$SM" >/dev/null 2>&1; then
  echo "snakemake is not on PATH after building the environment" >&2
  exit 1
fi

PYTHONWARNINGS="${PYTHONWARNINGS:-ignore::FutureWarning}" "$SM" \
  --snakefile "$SCRIPT_DIR/Snakefile" \
  --configfile "$CONFIG" \
  --cores "$(nproc)" \
  --printshellcmds \
  --show-failed-logs \
  --rerun-incomplete \
  "$@"

shopt -s nullglob
found=0
for gff in "$RESULTS"/*/*.ie.gff3; do
  sample="$(basename "$(dirname "$gff")")"
  cp -f "$gff" "$OUTDIR/${sample}.ie.gff3"
  cp -f "$RESULTS/$sample/${sample}.ie.gbk" "$OUTDIR/${sample}.ie.gbk"
  echo "GFF3:    $OUTDIR/${sample}.ie.gff3"
  echo "GenBank: $OUTDIR/${sample}.ie.gbk"
  echo "Report:  $OUTDIR/${sample}.ie.report.txt"
  found=1
done
shopt -u nullglob
write_reports
if [[ $found -eq 0 ]]; then
  echo "Snakemake finished without a GFF3" >&2
  exit 1
fi
