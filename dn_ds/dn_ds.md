```bash
set -euo pipefail

# =========[ 0. environment & config ]=========
# activate your env (contains bcftools, samtools, seqkit, mafft). gffread & paml need manual install.
micromamba activate bioinfo || true

# working dirs
ROOT=~/work/language_gwas
DAT=$ROOT/data
mkdir -p "$DAT"
cd "$DAT"

# input lists/maps you must provide
GENE_LIST=$DAT/GWAS_Suggestive_Genes_ID.txt          # one ENSG per line; header 'gene_id' allowed
ORTHMAP=$DAT/orth_hum_chimp.tsv                      # BioMart TSV (7 cols; see notes below)

# output area
OUT=$DAT/output
mkdir -p "${OUT}"/{logs,fa,cds,human,nean,chimp,aa,aln,tree,codeml,tmp}

# tools (adjust if needed)
# get the latest (v14 is the common one)
cd ~/bin
mkdir pal2nal.v14
cd pal2nal.v14
curl -LO http://www.bork.embl.de/pal2nal/distribution/pal2nal.v14.tar.gz
tar xzf pal2nal.v14.tar.gz

cd "$DAT"

PAL2NAL=~/bin/pal2nal.v14/pal2nal.v14/pal2nal.pl        # full path to pal2nal.pl
which gffread >/dev/null 2>&1 || { echo "[FATAL] gffread not found"; exit 1; }
which seqkit  >/dev/null 2>&1 || { echo "[FATAL] seqkit not found";  exit 1; }
which mafft   >/dev/null 2>&1 || { echo "[FATAL] mafft not found";   exit 1; }
which samtools>/dev/null 2>&1 || { echo "[FATAL] samtools not found";exit 1; }
which bcftools>/dev/null 2>&1 || { echo "[FATAL] bcftools not found";exit 1; }
which codeml  >/dev/null 2>&1 || { echo "[FATAL] codeml not found";  exit 1; }
[ -x "$PAL2NAL" ] || { echo "[FATAL] PAL2NAL not executable at $PAL2NAL"; exit 1; }

# =========[ 1. references ]=========
# Human GRCh37 (to match Vindija Neandertal)
REF=homo_sapiens.GRCh37.dna.fa
GTF=homo_sapiens.GRCh37.dna.gtf
if [ ! -s "$REF" ]; then
  wget -c ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
  gunzip -c Homo_sapiens.GRCh37.dna.toplevel.fa.gz > "$REF"
fi
if [ ! -s "$GTF" ]; then
  wget -c ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
  gunzip -c Homo_sapiens.GRCh37.87.gtf.gz > "$GTF"
fi
samtools faidx "$REF"

# Neandertal VCF (Vindija33.19; concatenated, indexed)
NEAN_VCF=Vindija33_19.merged.vcf.gz
NEAN_SAMPLE=Vindija33.19
if [ ! -s "$NEAN_VCF" ]; then
  base_url="http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19"
  for chr in {1..22} X; do
    wget -c "${base_url}/chr${chr}_mq25_mapab100.vcf.gz"
    wget -c "${base_url}/chr${chr}_mq25_mapab100.vcf.gz.tbi"
  done
  bcftools concat -a -Oz chr{1..22}_mq25_mapab100.vcf.gz chrX_mq25_mapab100.vcf.gz -o "$NEAN_VCF"
  bcftools index -t "$NEAN_VCF"
fi

# Chimp references
CHIMP_CDS="$DAT/Pan_troglodytes.Pan_tro_3.0.cds.all.fa"
CHIMP_GTF=Pan_troglodytes.gtf
if [ ! -s "$CHIMP_CDS" ]; then
  # Replace your CHIMP_CDS file with CDS
wget https://ftp.ensembl.org/pub/release-109/fasta/pan_troglodytes/cds/Pan_troglodytes.Pan_tro_3.0.cds.all.fa.gz
gunzip -c Pan_troglodytes.Pan_tro_3.0.cds.all.fa.gz > "$DAT"/Pan_troglodytes.Pan_tro_3.0.cds.all.fa
fi

if [ ! -s "$CHIMP_GTF" ]; then
  wget -c -O Pan_troglodytes.gtf.gz https://ftp.ensembl.org/pub/release-109/gtf/pan_troglodytes/Pan_troglodytes.Pan_tro_3.0.109.gtf.gz
  gunzip -c Pan_troglodytes.gtf.gz > "$CHIMP_GTF"
fi

# =========[ 2. build neandertal pseudogenome ]=========
PSG=${OUT}/neandertal_pseudogenome.fa
if [ ! -s "$PSG" ]; then
  bcftools index -t -f "$NEAN_VCF" || true
  bcftools consensus -f "$REF" -s "$NEAN_SAMPLE" "$NEAN_VCF" > "$PSG"
fi

# =========[ 3. extract CDS for human & neandertal (same coordinates) ]=========
HUMAN_CDS=${OUT}/human/all_CDS.fa
NEAN_CDS=${OUT}/nean/all_CDS.fa

if [ ! -s "$HUMAN_CDS" ]; then
  gffread -g "$REF" -x "$HUMAN_CDS" "$GTF" 2> "${OUT}/logs/gffread_human.log"
fi
if [ ! -s "$NEAN_CDS" ]; then
  gffread -g "$PSG"  -x "$NEAN_CDS" "$GTF"  2> "${OUT}/logs/gffread_nean.log"
fi

# =========[ 4. helper maps (human tx→gene; chimp protein→tx; chimp tx→gene) ]=========
TX2GENE=${OUT}/tx2gene.tsv
if [ ! -s "$TX2GENE" ]; then
  awk -F'\t' '$3=="transcript" && $0 ~ /transcript_biotype "protein_coding"/ {
    if (match($0,/transcript_id "([^"]+)"/,m1) && match($0,/gene_id "([^"]+)"/,m2))
      print m1[1] "\t" m2[1];
  }' "$GTF" | sort -u > "$TX2GENE"
fi

CHIMP_P2T=${OUT}/chimp_protein2tx.tsv
if [ ! -s "$CHIMP_P2T" ]; then
  awk -F'\t' '$3=="CDS"{
    pid=""; tid="";
    if (match($0,/protein_id "([^"]+)"/,m1)) pid=m1[1];
    if (match($0,/transcript_id "([^"]+)"/,m2)) tid=m2[1];
    if (pid!="" && tid!="") print pid "\t" tid;
  }' "$CHIMP_GTF" | sort -u > "$CHIMP_P2T"
fi

CHIMP_T2G=${OUT}/chimp_tx2gene.tsv
if [ ! -s "$CHIMP_T2G" ]; then
  awk -F'\t' '$3=="transcript"{
    tid=""; gid=""; gname="";
    if (match($0,/transcript_id "([^"]+)"/,m1)) tid=m1[1];
    if (match($0,/gene_id "([^"]+)"/,m2))      gid=m2[1];
    if (match($0,/gene_name "([^"]+)"/,m3))  gname=m3[1];
    if (tid!="") print tid "\t" gid "\t" gname;
  }' "$CHIMP_GTF" | sort -u > "$CHIMP_T2G"
fi

# =========[ 5. per-gene build: CDS → AA MSA → codon alignment ]=========
# BioMart orthmap columns expected:
# 1: Gene stable ID (ENSG...) | 2: Gene name | 3: Transcript stable ID (ENST...)
# 4: Chimp gene ID (ENSPTRG...) | 5: Chimp gene name | 6: Chimp protein OR transcript ID (ENSPTRP* OR ENSPTRT*) | 7: homology type
[ -s "$ORTHMAP" ] || { echo "[FATAL] missing $ORTHMAP"; exit 1; }
[ -s "$GENE_LIST" ] || { echo "[FATAL] missing $GENE_LIST"; exit 1; }

# clean pal2nal logs
: > "${OUT}/logs/pal2nal.log"
: > "${OUT}/logs/pal2nal_fail.list"


# process (header skip)
awk 'NR>1 && NF {print $1}' "$GENE_LIST" | sed 's/\r$//' | while read -r GENE; do
  echo "[*] ${GENE}"

  # (A) list all human transcripts for this gene
  awk -v g="$GENE" '$2==g{print $1}' "$TX2GENE" > "${OUT}/human/${GENE}.tx"
  if [ ! -s "${OUT}/human/${GENE}.tx" ]; then
    echo "  - no human transcripts in tx2gene" | tee -a "${OUT}/logs/missing_human_tx.log"
    continue
  fi

  # (B) pull those CDS records (human) and choose the longest
  seqkit grep -n -f "${OUT}/human/${GENE}.tx" "$HUMAN_CDS" > "${OUT}/human/${GENE}.all.fa" || true
  if [ ! -s "${OUT}/human/${GENE}.all.fa" ]; then
    echo "  - no human CDS records" | tee -a "${OUT}/logs/missing_human_tx.log"; continue
  fi
  seqkit fx2tab -nl "${OUT}/human/${GENE}.all.fa" | sort -k2,2nr | head -n1 | cut -f1 > "${OUT}/human/${GENE}.id"
  TX=$(cat "${OUT}/human/${GENE}.id")
  seqkit grep -n -p "$TX" "${OUT}/human/${GENE}.all.fa" > "${OUT}/human/${GENE}.fa"
  [ -s "${OUT}/human/${GENE}.fa" ] || { echo "  - representative human CDS missing"; continue; }

  # (C) neandertal: find the same transcript id
  seqkit grep -n -p "$TX" "$NEAN_CDS" > "${OUT}/nean/${GENE}.fa" || true
  if [ ! -s "${OUT}/nean/${GENE}.fa" ]; then
    echo "  - no Neandertal match for $TX" | tee -a "${OUT}/logs/missing_nean_tx.log"; continue
  fi

  # (D) chimp: resolve from orthmap; prefer col6 (ENSPTRT* or ENSPTRP*)
  # ORTHMAP columns (tab-separated, with header):
  # 1 Gene stable ID | 2 Gene name | 3 Transcript stable ID | 4 Chimp gene stable ID
  # 5 Chimp gene name | 6 Chimp protein or transcript stable ID | 7 Chimp homology type
  read chimp_geneid chimp_sym chimp_token < <(
    awk -v g="$GENE" 'BEGIN{FS=OFS="\t"}
      NR==1 {next}          # skip header
      $1==g {
        cg = ($4==""?"NA":$4);
        cs = ($5==""?"NA":$5);
        ct = ($6==""?"NA":$6);
        print cg, cs, ct; exit
      }' "$ORTHMAP" | tr -d '\r'
  )

  if [ -z "$chimp_geneid" ] && [ -z "$chimp_sym" ] && [ -z "$chimp_token" ]; then
    echo "  - no row in orthmap" | tee -a "${OUT}/logs/missing_chimp.log"
    continue
  fi

  CHIMP_TX=""
  if [[ "$chimp_token" =~ ^ENSPTRT ]]; then
    CHIMP_TX="$chimp_token"
  elif [[ "$chimp_token" =~ ^ENSPTRP ]]; then
    CHIMP_TX=$(awk -v p="$chimp_token" 'BEGIN{FS=OFS="\t"} $1==p{print $2; exit}' "$CHIMP_P2T")
  fi

  # fallback by chimp gene (id -> tx, or symbol -> tx); pick longest
  if [ -z "$CHIMP_TX" ]; then
    rm -f "${OUT}/chimp/${GENE}.tx"
    if [ -n "$chimp_geneid" ] && [ "$chimp_geneid" != "NA" ]; then
      awk -v gid="$chimp_geneid" 'BEGIN{FS=OFS="\t"} $2==gid{print $1}' "$CHIMP_T2G" > "${OUT}/chimp/${GENE}.tx"
    fi
    if [ ! -s "${OUT}/chimp/${GENE}.tx" ] && [ -n "$chimp_sym" ] && [ "$chimp_sym" != "NA" ]; then
      awk -v s="$chimp_sym" 'BEGIN{FS=OFS="\t"} $3==s{print $1}' "$CHIMP_T2G" > "${OUT}/chimp/${GENE}.tx"
    fi

    if [ -s "${OUT}/chimp/${GENE}.tx" ]; then
      seqkit grep -n -f "${OUT}/chimp/${GENE}.tx" "$CHIMP_CDS" > "${OUT}/chimp/${GENE}.all.fa" || true
      if [ -s "${OUT}/chimp/${GENE}.all.fa" ]; then
        seqkit fx2tab -nl "${OUT}/chimp/${GENE}.all.fa" | sort -k2,2nr | head -n1 | cut -f1 > "${OUT}/chimp/${GENE}.id"
        CHIMP_TX=$(cat "${OUT}/chimp/${GENE}.id}")
      fi
    fi
  fi

  if [ -z "$CHIMP_TX" ]; then
    echo "  - could not resolve chimp transcript (token=${chimp_token} geneid=${chimp_geneid} sym=${chimp_sym})" \
      | tee -a "${OUT}/logs/missing_chimp_seq.log"
    continue
  fi

  # fetch chosen chimp transcript CDS
  seqkit grep -nrp "$CHIMP_TX" "$CHIMP_CDS" > "${OUT}/chimp/${GENE}.fa" || true
  if [ ! -s "${OUT}/chimp/${GENE}.fa" ]; then
    echo "  - chimp CDS not found for tx=${CHIMP_TX}" | tee -a "${OUT}/logs/missing_chimp_seq.log"; continue
  fi

  # (E) standardize headers to human/nean/chimp and enforce frame (len%3==0)
  seqkit seq -n -i "${OUT}/human/${GENE}.fa" >/dev/null || true
  for sp in human nean chimp; do
    src="${OUT}/${sp}/${GENE}.fa"
    tgt="${OUT}/${sp}/${GENE}.ren.fa"
    seqkit replace -p '.*' -r "$sp" "$src" > "$tgt"
  done

  # quick frame check (drop if not divisible by 3)
  bad=$( (seqkit fx2tab -nl "${OUT}/human/${GENE}.ren.fa"; \
          seqkit fx2tab -nl "${OUT}/nean/${GENE}.ren.fa";  \
          seqkit fx2tab -nl "${OUT}/chimp/${GENE}.ren.fa") | awk '$2%3!=0' | wc -l )
  if [ "$bad" -ne 0 ]; then
    echo "  - frame not multiple of 3; skipping" | tee -a "${OUT}/logs/frame_issues.log"
    continue
  fi

  # (F) translate → align → back-translate
  seqkit translate -f 1 "${OUT}/human/${GENE}.ren.fa" -o "${OUT}/aa/${GENE}.human.aa.fa" || continue
  seqkit translate -f 1 "${OUT}/nean/${GENE}.ren.fa"  -o "${OUT}/aa/${GENE}.nean.aa.fa"  || continue
  seqkit translate -f 1 "${OUT}/chimp/${GENE}.ren.fa" -o "${OUT}/aa/${GENE}.chimp.aa.fa" || continue

  cat "${OUT}/aa/${GENE}.human.aa.fa" \
      "${OUT}/aa/${GENE}.nean.aa.fa"  \
      "${OUT}/aa/${GENE}.chimp.aa.fa" > "${OUT}/aa/${GENE}.aa.fa"

  mafft --auto "${OUT}/aa/${GENE}.aa.fa" > "${OUT}/aln/${GENE}.aa.aln.fa" 2>> "${OUT}/logs/mafft.log" || continue

  cat "${OUT}/human/${GENE}.ren.fa" \
      "${OUT}/nean/${GENE}.ren.fa"  \
      "${OUT}/chimp/${GENE}.ren.fa" > "${OUT}/cds/${GENE}.cds.fa"

  tmpcodon=$(mktemp "${OUT}/tmp/${GENE}.XXXXXX")
  if "$PAL2NAL" "${OUT}/aln/${GENE}.aa.aln.fa" "${OUT}/cds/${GENE}.cds.fa" -output paml \
        > "$tmpcodon" 2>> "${OUT}/logs/pal2nal.log"; then
    mv -f "$tmpcodon" "${OUT}/codeml/${GENE}.codon.phy"
  else
    echo "[PAL2NAL FAIL] ${GENE}" | tee -a "${OUT}/logs/pal2nal_fail.list"
    rm -f "$tmpcodon"
    continue
  fi

  # (G) default tree (Neandertal sister to human)
  echo "(chimp,(human,nean));" > "${OUT}/tree/${GENE}.nwk"
done
```
After the loop, count how many codon alignments you actually produced:
```bash
ls -1 "$OUT/codeml"/*.codon.phy 2>/dev/null | wc -l
```

```bash
# create missing trees, if any
shopt -s nullglob
for CODON in "$OUT"/codeml/*.codon.phy; do
  G=$(basename "$CODON" .codon.phy)
  T="$OUT/tree/${G}.nwk"
  [ -s "$T" ] || echo "(chimp,(human,nean));" > "$T"
done
shopt -u nullglob
```

run codeml (robust loop)

```bash
# run codeml per gene
#!/usr/bin/env bash
set -euo pipefail

: "${OUT:?Set OUT to your output directory (e.g., export OUT=/path/to/output)}"

mkdir -p "$OUT/logs" "$OUT/tree" "$OUT/codeml"

# Make the for-loop skip when there are no matches
shopt -s nullglob

# A small helper for messages
log() { printf '%s\n' "$*" >&2; }

for CODON in "$OUT"/codeml/*.codon.phy; do
  G=$(basename "$CODON" .codon.phy)

  TREE_RAW="$OUT/tree/${G}.nwk"
  TREE_FG="$OUT/tree/${G}.fg.nwk"
  CTL="$OUT/codeml/${G}.ctl"
  OUTFILE="$OUT/codeml/${G}.codeml.out"
  LOGF="$OUT/logs/codeml_${G}.log"

  # 0) inputs exist?
  if [[ ! -s "$CODON" ]]; then
    log "[WARN] $G: missing codon file $CODON"
    continue
  fi
  if [[ ! -s "$TREE_RAW" ]]; then
    log "[WARN] $G: missing tree $TREE_RAW"
    continue
  fi

  # 1) minimal sanity checks that don't depend on label placement
  NSEQ=$(awk 'NR==1{print $1}' "$CODON")
  if [[ "$NSEQ" != "3" ]]; then
    log "[SKIP] $G: PHYLIP header reports $NSEQ sequences (expected 3)"
    sed -n '1,6p' "$CODON" | sed 's/^/   CODON: /' >&2
    echo "[SKIP] $G: PHYLIP header $NSEQ (expected 3)" >> "$OUT/logs/bad_phylip_header.log"
    continue
  fi

  if ! grep -q 'human' "$TREE_RAW"; then
    log "[SKIP] $G: 'human' not found in tree"
    log "   tree=$(cat "$TREE_RAW")"
    echo "[SKIP] $G: human not in tree" >> "$OUT/logs/tree_name_mismatch.log"
    continue
  fi

  # 2) mark human as foreground (avoid duplicating #1)
  if grep -q 'human#1' "$TREE_RAW"; then
    cp "$TREE_RAW" "$TREE_FG"
  else
    sed 's/human/human#1/' "$TREE_RAW" > "$TREE_FG"
  fi

  # 3) write ctl using RELATIVE paths (we cd into codeml/)
  # IMPORTANT: The heredoc is CLOSED exactly by EOF on its own line. Do not add anything after EOF.
  cat > "$CTL" <<EOF
seqfile = $(basename "$CODON")
treefile = ../tree/$(basename "$TREE_FG")
outfile  = $(basename "$OUTFILE")

noisy    = 9
verbose  = 1
runmode  = 0

seqtype  = 1       * 1:codons
CodonFreq= 2       * F3x4
clock    = 0
model    = 2       * branch model (allows different omega on FG)
NSsites  = 0
icode    = 0

fix_kappa= 0
kappa    = 2

fix_omega= 0
omega    = 0.2

cleandata= 1
getSE    = 0
Small_Diff = .5e-6
EOF

  # 4) quick sanity check the ctl
  if ! grep -q '^outfile' "$CTL"; then
    log "[ERR] $G: ctl malformed (missing outfile)"; head -n 30 "$CTL" >&2; continue
  fi
  if ! grep -q '^model[[:space:]]*=' "$CTL"; then
    log "[ERR] $G: ctl malformed (missing model)"; head -n 30 "$CTL" >&2; continue
  fi

  # 5) run codeml from within the codeml dir so relative paths resolve
  (
    cd "$OUT/codeml"
    # Remove any stale output so we know if codeml produced a fresh file
    rm -f "$(basename "$OUTFILE")"
    # Run codeml
    if ! codeml "$(basename "$CTL")" > "$LOGF" 2>&1; then
      log "[ERR] $G: codeml exited non-zero (see $LOGF)"
      continue
    fi
  )

  # 6) basic result validation
  if [[ ! -s "$OUTFILE" ]]; then
    log "[ERR] $G: codeml produced no outfile (see $LOGF)"
    continue
  fi

  # Common failure mode: 0 usable sites with cleandata=1
  if grep -qi '0 usable sites' "$LOGF"; then
    log "[WARN] $G: 0 usable sites after cleandata — alignment likely too gappy (see $LOGF)"
  else
    log "[OK]   $G: finished (outfile $(basename "$OUTFILE"))"
  fi

done

# restore default globbing behavior for the shell/session
shopt -u nullglob
```
