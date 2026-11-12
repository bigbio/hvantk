#!/usr/bin/env bash
# Stage the ml-features eval, once per cluster. Login-node orchestrator: downloads go
# through the image (each plugin's own downloader, or curl for the INSIDER BED, which
# has none), the three variant-level files are cut to chromosome 22 by slice_text.py as
# short jobs, the gene index and the HGNC table are copied from the newest gene-index
# control that fully passed (it has control.passed), and the eval spec and gene list
# are copied next to them. Writes staged.sha256.
#
#   HVANTK_SITE_PROFILE=... DBNSFP_RAW=/abs/path/dbNSFP4.9a_variant... \
#     [CLINVAR_VERSION=20261108] bash skills/ml-features/eval/stage.sh
#
# DBNSFP_RAW: dbNSFP is bring-your-own (academic registration with the dbNSFP maintainers; there is
# no download URL), so name a dbNSFP 4.9a variant file you already have, gzipped or
# block-gzipped: the all-chromosomes file or the chr22 one. Only its chr22 rows are
# staged. CLINVAR_VERSION pins the ClinVar release (a dated weekly archive); the
# reference values were produced from the default.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
bash "$HERE/../../_site/check-site-profile.sh" >/dev/null
# shellcheck disable=SC1090
source "$HVANTK_SITE_PROFILE"
R="$(cd "$HVANTK_REPO" && pwd -P)"
[[ "$HERE" == "$R"/* ]] || { echo "ERROR: this eval folder ($HERE) is not inside HVANTK_REPO ($R)" >&2; exit 1; }
# Paths handed to a job use the profile's spelling of the checkout (what the job binds).
REL="${HERE#"$R"/}"
SLICE="$HVANTK_REPO/$REL/slice_text.py"
RUN_IN_IMAGE="$HVANTK_REPO/skills/hpc-hail-job/templates/run-in-image.sh"
SUBMIT="$HVANTK_REPO/skills/hpc-hail-job/templates/submit-hail-job.sh"

: "${DBNSFP_RAW:?set DBNSFP_RAW to a dbNSFP 4.9a variant file (see the header of this script)}"
[[ -s "$DBNSFP_RAW" ]] || { echo "ERROR: DBNSFP_RAW ($DBNSFP_RAW) is missing or empty" >&2; exit 1; }
CLINVAR_VERSION="${CLINVAR_VERSION:-20261108}"
[[ "$CLINVAR_VERSION" =~ ^[0-9]{8}$ ]] || { echo "ERROR: CLINVAR_VERSION must be YYYYMMDD" >&2; exit 1; }

EVAL_ROOT="${HVANTK_EVAL_DATA:-$HVANTK_WORK/hvantk-eval}"
STAGE="$EVAL_ROOT/ml-features"
GI="$EVAL_ROOT/gene-index"
mkdir -p "$STAGE/download/clinvar" "$STAGE/download/insider" "$STAGE/upstream"
for d in clinvar dbnsfp insider_variants insider_interfaces uniprot_ptm gnomad_metrics gevir; do
  mkdir -p "$STAGE/raw/$d"
done

# Every raw dir must hold exactly one file matching *pattern* and nothing else (copied
# from rerank-evaluate/eval/stage.sh: a partial download leaves a <name>.<random> temp
# file behind rather than disappearing).
validate_raw_dir() {
  local dir="$1" pattern="$2" label="$3"
  local -a matching=() extra=()
  local f
  while IFS= read -r -d '' f; do
    case "$(basename "$f")" in
      $pattern) matching+=("$f") ;;
      *) extra+=("$f") ;;
    esac
  done < <(find "$dir" -maxdepth 1 -type f ! -name '.*' -print0)
  if (( ${#extra[@]} > 0 )); then
    echo "ERROR: $label raw dir ($dir) holds an unexpected file: ${extra[0]} (a partial download looks like <name>.<random>; remove it and re-stage)" >&2
    exit 1
  fi
  if (( ${#matching[@]} != 1 )); then
    echo "ERROR: $label raw dir ($dir) must hold exactly one file matching '$pattern'; found ${#matching[@]}: ${matching[*]:-(none)}" >&2
    exit 1
  fi
}

# slice IN OUT CONTIG... : run slice_text.py as a short job unless OUT already exists.
slice() {
  local in="$1" out="$2" name="$3" time="$4"; shift 4
  if [[ -s "$out" ]]; then echo "already staged: $out"; return; fi
  local -a contig_args=()
  local c
  for c in "$@"; do contig_args+=(--contig "$c"); done
  # A hidden temp name in the same directory: it keeps the final name's suffix
  # (slice_text.py writes BGZF only for *.bgz), validate_raw_dir ignores it, and a
  # failed run never looks staged.
  local tmp
  tmp="$(dirname "$out")/.part-$(basename "$out")"
  rm -f "$tmp"
  bash "$SUBMIT" --name "stage-mlf-$name" --cpus 2 --mem-gb 8 --time "$time" --wait -- \
    python "$SLICE" --input "$in" --output "$tmp" "${contig_args[@]}" >/dev/null
  mv "$tmp" "$out"
  # Hadoop's local filesystem writes a checksum next to a file Hail writes; it is named
  # after the temp name, so drop it rather than leave it behind.
  rm -f "$(dirname "$out")/.$(basename "$tmp").crc"
}

# --- upstream: the newest gene-index control that fully passed ---
PICK="$(python3 -c '
import os, sys
root = sys.argv[1]
if os.path.isdir(root):
    for name in sorted(os.listdir(root), reverse=True):
        if os.path.isfile(os.path.join(root, name, "control.passed")):
            print(os.path.join(root, name))
            sys.exit(0)
sys.exit(1)
' "$GI/control")" || {
  echo "ERROR: no gene-index control directory under $GI/control has control.passed; run skills/gene-index/eval/stage.sh and control.sh first" >&2
  exit 1
}
echo "using gene-index control: $PICK"
for t in gene-index.ht hgnc-lookup.ht; do
  [[ -e "$PICK/pos/$t/_SUCCESS" ]] || { echo "ERROR: $PICK/pos/$t has no _SUCCESS marker" >&2; exit 1; }
  rsync -a --delete "$PICK/pos/$t/" "$STAGE/upstream/$t/"
done
cp "$PICK/pos/hgnc-lookup.ht.provenance.json" "$STAGE/upstream/hgnc-lookup.ht.provenance.json"
python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['gene_index_content_sha256'])" \
  "$PICK/control.passed" > "$STAGE/upstream/gene-index.content_sha256"

# --- ClinVar: the pinned release, cut to chr22 ---
shopt -s nullglob
clinvar_dl=("$STAGE"/download/clinvar/*.vcf.gz)
shopt -u nullglob
if (( ${#clinvar_dl[@]} == 0 )); then
  bash "$RUN_IN_IMAGE" -- hvantk download clinvar --version "$CLINVAR_VERSION" \
    --output-dir "$STAGE/download/clinvar"
  clinvar_dl=("$STAGE"/download/clinvar/*.vcf.gz)
fi
(( ${#clinvar_dl[@]} == 1 )) || { echo "ERROR: expected one ClinVar VCF under $STAGE/download/clinvar" >&2; exit 1; }
slice "${clinvar_dl[0]}" "$STAGE/raw/clinvar/clinvar_${CLINVAR_VERSION}.chr22.vcf.bgz" clinvar 0-00:30:00 22 chr22
validate_raw_dir "$STAGE/raw/clinvar" "*.vcf.bgz" clinvar

# --- dbNSFP: bring-your-own, cut to chr22 (a full pass over the file) ---
slice "$DBNSFP_RAW" "$STAGE/raw/dbnsfp/dbNSFP4.9a_variant.chr22.bgz" dbnsfp 0-03:00:00 22 chr22
validate_raw_dir "$STAGE/raw/dbnsfp" "*.bgz" dbnsfp

# --- INSIDER variants: the BED has no downloader; its URL is in the plugin catalog ---
INSIDER_BED="$STAGE/download/insider/Whole_Human_Interactome_Interface_hg38.bed"
if [[ ! -s "$INSIDER_BED" ]]; then
  URL="$(python3 -c '
import json, sys
urls = [f.get("url", "") for d in json.load(open(sys.argv[1])) for f in d.get("files", [])]
print(next((u for u in urls if u.endswith(".bed")), ""))
' "$HVANTK_REPO/hvantk/skills/insider/catalog/datasets.json")"
  [[ -n "$URL" ]] || { echo "ERROR: no .bed URL in the INSIDER catalog" >&2; exit 1; }
  curl -fsSL -o "$INSIDER_BED.part" "$URL" && mv "$INSIDER_BED.part" "$INSIDER_BED"
fi
slice "$INSIDER_BED" "$STAGE/raw/insider_variants/insider.chr22.bed" insider 0-00:30:00 chr22 22
validate_raw_dir "$STAGE/raw/insider_variants" "*.bed" insider_variants

# --- small whole-genome inputs, each through its plugin's downloader ---
stage_download() {
  local dir="$1" pattern="$2" label="$3"; shift 3
  shopt -s nullglob
  local -a have=("$dir"/$pattern)
  shopt -u nullglob
  if (( ${#have[@]} > 0 )); then echo "already staged: $dir"; else bash "$RUN_IN_IMAGE" -- "$@"; fi
  validate_raw_dir "$dir" "$pattern" "$label"
}
stage_download "$STAGE/raw/insider_interfaces" "H_sapiens_interfacesALL.txt" insider_interfaces \
  hvantk download insider-interfaces --output "$STAGE/raw/insider_interfaces/H_sapiens_interfacesALL.txt"
stage_download "$STAGE/raw/uniprot_ptm" "*.tsv" uniprot_ptm \
  hvantk download uniprot-ptm --output-dir "$STAGE/raw/uniprot_ptm"
stage_download "$STAGE/raw/gnomad_metrics" "*.bgz" gnomad_metrics \
  hvantk download gnomad-metrics --output "$STAGE/raw/gnomad_metrics/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
stage_download "$STAGE/raw/gevir" "*.tsv.bgz" gevir \
  hvantk download gevir --output "$STAGE/raw/gevir/gevir.tsv.bgz"

# --- the spec and gene list the agent is handed ---
cp "$HERE/spec.yaml" "$STAGE/spec.yaml"
cp "$HERE/genes.txt" "$STAGE/genes.txt"

( cd "$STAGE" && sha256sum raw/*/* upstream/gene-index.content_sha256 \
    upstream/hgnc-lookup.ht.provenance.json spec.yaml genes.txt > staged.sha256 )
cat "$STAGE/staged.sha256"
echo "staged under $STAGE (gene index from $PICK)"
