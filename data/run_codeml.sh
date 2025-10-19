#!/usr/bin/env bash
set -euo pipefail

# ========= CONFIG =========
: "${OUT:?Set OUT to your output directory (e.g., export OUT=/path/to/output)}"
# Branch-site model with human as foreground
SEQTYPE=1           # 1: codons
CLEANDATA=1         # 1: remove codons with gaps/ambiguous
CODONFREQ=2
RUNMODE=0
CLOCK=0
MODEL=2             # branch-site
NSSITES=2
FIX_OMEGA=0
OMEGA=1.5
FIX_KAPPA=0
KAPPA=2
NOISY=9
# =========================

# Folders
mkdir -p "$OUT/logs" "$OUT/tree" "$OUT/codeml"

# Helpful logger
log(){ printf '%s\n' "$*" >&2; }

# Sanity: codeml present?
if ! command -v codeml >/dev/null 2>&1; then
  log "[ERR] 'codeml' not found in PATH. Install PAML or add it to PATH."
  exit 127
fi

# Make loop skip if no matches
shopt -s nullglob

# Counters
ok=0; fail=0; skip=0

for CODON in "$OUT"/codeml/*.codon.phy; do
  G=$(basename "$CODON" .codon.phy)

  TREE_RAW="$OUT/tree/${G}.nwk"
  TREE_FG="$OUT/tree/${G}.fg.nwk"
  CTL="$OUT/codeml/${G}.ctl"
  OUTFILE="$OUT/codeml/${G}.codeml.out"
  LOGF="$OUT/logs/codeml_${G}.log"

  # --- Input checks ---
  [[ -s "$CODON"    ]] || { log "[WARN] $G: missing codon file $CODON"; ((skip++)); continue; }
  [[ -s "$TREE_RAW" ]] || { log "[WARN] $G: missing tree $TREE_RAW";   ((skip++)); continue; }

  # PHYLIP header: expect 3 seqs
  NSEQ=$(awk 'NR==1{print $1}' "$CODON")
  if [[ "$NSEQ" != "3" ]]; then
    log "[SKIP] $G: PHYLIP header reports $NSEQ sequences (expected 3)"
    echo "$G $NSEQ" >> "$OUT/logs/bad_phylip_header.log"
    ((skip++)); continue
  fi

  # Tree must include 'human'
  if ! grep -q 'human' "$TREE_RAW"; then
    log "[SKIP] $G: 'human' not found in $TREE_RAW"
    echo "$G" >> "$OUT/logs/tree_name_mismatch.log"
    ((skip++)); continue
  fi

  # --- Foreground tree ---
  # If already tagged, copy; else tag 'human' -> 'human#1'
  if grep -q 'human#1' "$TREE_RAW"; then
    cp -f "$TREE_RAW" "$TREE_FG"
  else
    # Word-boundary match for 'human' (GNU sed supports \< \>)
    sed 's/\<human\>/human#1/g' "$TREE_RAW" > "$TREE_FG"
  fi

  # --- Control file ---
  # codeml is run from $OUT/codeml, so use basenames/relative paths
  cat > "$CTL" <<EOF
seqfile = $(basename "$CODON")
treefile = ../tree/$(basename "$TREE_FG")
outfile  = $(basename "$OUTFILE")
noisy    = $NOISY
runmode  = $RUNMODE
seqtype  = $SEQTYPE
CodonFreq= $CODONFREQ
clock    = $CLOCK
aaDist   = 0
model    = $MODEL
NSsites  = $NSSITES
fix_omega= $FIX_OMEGA
omega    = $OMEGA
cleandata= $CLEANDATA
fix_kappa= $FIX_KAPPA
kappa    = $KAPPA
EOF

  # --- Run codeml ---
  (
    set +e
    cd "$OUT/codeml" || exit 1
    rm -f "$(basename "$OUTFILE")"
    codeml "$(basename "$CTL")" >"$LOGF" 2>&1
    status=$?
    # success requires exit status 0 AND a nonempty outfile
    if [[ $status -ne 0 ]]; then
      echo "$G" >> "$OUT/logs/codeml_failed.log"
      exit 25
    fi
    if [[ ! -s "$(basename "$OUTFILE")" ]]; then
      echo "$G" >> "$OUT/logs/codeml_empty_outfile.log"
      exit 26
    fi
    exit 0
  )
  case $? in
    0)  log "[OK]   $G: $(basename "$OUTFILE")"; ((ok++)) ;;
    25) log "[ERR]  $G: codeml exited non-zero (see $LOGF)"; ((fail++));;
    26) log "[ERR]  $G: codeml produced no outfile (see $LOGF)"; ((fail++));;
    *)  log "[ERR]  $G: unknown error (see $LOGF)"; ((fail++));;
  esac

done

# Restore if you care
shopt -u nullglob || true

log "==== SUMMARY ===="
log "OK: $ok   FAIL: $fail   SKIP: $skip"
[[ $fail -eq 0 ]] || exit 1
