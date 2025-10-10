```bash
micromamba activate bioinfo
set -uex -o pipefail

cd ~/work/language_gwas/dn_ds
# A) Files to download
cd ~/work/language_gwas
mkdir -p data
cd data
# Download human reference genome FASTA
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
gunzip -c Homo_sapiens.GRCh37.dna.toplevel.fa.gz > homo_sapiens.GRCh37.dna.fa

# Download human reference genome GTF
wget ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
gunzip -c Homo_sapiens.GRCh37.87.gtf.gz > homo_sapiens.GRCh37.dna.gtf
# Download Neandertal VCF from max planck
base_url="http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19"

for chr in {1..22} X; do
    wget "${base_url}/chr${chr}_mq25_mapab100.vcf.gz"
    wget "${base_url}/chr${chr}_mq25_mapab100.vcf.gz.tbi"
done

bcftools concat -a -Oz \
  chr{1..22}_mq25_mapab100.vcf.gz chrX_mq25_mapab100.vcf.gz \
  -o Vindija33_19.merged.vcf.gz
bcftools index -t Vindija33_19.merged.vcf.gz

# Download Chimp reference CDS FASTA (or cDNA) and ortholog mapping
wget https://ftp.ensembl.org/pub/current_fasta/pan_troglodytes/cdna/Pan_troglodytes.Pan_tro_3.0.cdna.all.fa.gz
gunzip -c Pan_troglodytes.Pan_tro_3.0.cdna.all.fa.gz > Pan_troglodytes.Pan_tro_3.0.cdna.all.fa

wget -O Pan_troglodytes.gtf.gz \
  https://ftp.ensembl.org/pub/release-109/gtf/pan_troglodytes/Pan_troglodytes.Pan_tro_3.0.109.gtf.gz
gunzip -c Pan_troglodytes.gtf.gz > Pan_troglodytes.gtf
# For making the ortholog map table use the BioMart ensemble website (https://www.ensembl.org/biomart/martview/6f33a9c67e1e24558a26ce74e673619f). after that, filter it with the gene ID list you have from the GWAS investigation. then choose the chimp relatives and download the results as .tsv file. change the name of the donwloaded folr to "orth_hum-chimp.tsv". The dataset used for getting this file is from GRCh38, in contrary to the human ref genome used in this project which is GRCh37 to better the Nean VCF be aligned with the reference genome.
# GWAS language-linked candidates (Ensembl Gene IDs or gene symbols)
```

# Tools needed
# you have to install gffread and paml. the micromamba does not contain those two packages.
# Multiple sequence alignment: mafft
# Back-translation for codons: pal2nal.pl (PAL2NAL) or codon aligner (MACSE/PRANK)
# dN/dS: PAML/codeml (or HyPhy, but I’ll show PAML here)


# B) Unix pipeline to build codon alignments per gene
# This creates a Neandertal pseudo-genome by applying the VCF to the human reference, extracts CDS per gene for Human & Neandertal using the same GTF, fetches the Chimp ortholog CDS, aligns protein, then back-translates to codons.



```bash
# 0) path
REF=homo_sapiens.GRCh37.dna.fa
GTF=homo_sapiens.GRCh37.dna.gtf
NEAN_VCF=Vindija33_19.merged.vcf.gz # Neandertal VCF (indexed)
NEAN_SAMPLE=Vindija33.19 # Sample name inside VCF (if multi-sample)
GENE_LIST=GWAS_Suggestive_Genes_ID.txt # One gene ID per line (e.g., ENSG...)
OUT=output        # Output directory
HUMAN_CDS=${OUT}/human/all_CDS.fa
NEAN_CDS=${OUT}/nean/all_CDS.fa
CHIMP_CDS=Pan_troglodytes.Pan_tro_3.0.cdna.all.fa
ORTHMAP=orth_hum_chimp.tsv
GENE_LIST=GWAS_Suggestive_Genes_ID.txt
CHIMP_GTF=Pan_troglodytes.gtf

mkdir -p ${OUT}/{logs,fa,cds,human,nean,chimp,aa,aln,tree,codeml}

# 1) index reference
samtools faidx ${REF}

# 2) make Nenderthal pseudo-genome from VCF (consensus)
# add filters/masks for low quality sites
bcftools index -t -f ${NEAN_VCF}
bcftools consensus -f ${REF} -s ${NEAN_SAMPLE} ${NEAN_VCF} > ${OUT}/neandertal_pseudogenome.fa

# 3) extract CDS fasta per gene for Human & Neandertal using same coordinates
# Requires gffread (https://ccb.jhu.edu/software/stringtie/gff.shtml)

# For HUMAN CDS (concatenate exons in CDS, keep transcript per gene)
gffread -g ${REF} -x ${OUT}/human/all_CDS.fa ${GTF} \
  2> ${OUT}/logs/gffread_human.log

# For NEANDERTAL CDS using the pseudogenome
gffread -g ${OUT}/neandertal_pseudogenome.fa -x ${OUT}/nean/all_CDS.fa ${GTF} \
  2> ${OUT}/logs/gffread_nean.log

#check
grep -m1 '^>' ${OUT}/human/all_CDS.fa
grep -m1 '^>' ${OUT}/nean/all_CDS.fa

# 4) Pull gene-ortholog mapping Human↔Chimp (via BioMart) to a TSV: gene_id, chimp_transcript_id

# Build tx2gene map from your GTF
awk -F'\t' '
  $3=="transcript" && $0 ~ /transcript_biotype "protein_coding"/ {
    if (match($0,/transcript_id "([^"]+)"/,m1) && match($0,/gene_id "([^"]+)"/,m2)) {
      print m1[1] "\t" m2[1]
    }
}' homo_sapiens.GRCh37.dna.gtf | sort -u > ${OUT}/tx2gene.tsv

# quick peek
head ${OUT}/tx2gene.tsv


# Build a chimp transcript→gene-name map (so symbols can be resolved to transcripts)
# Now extract transcript_id ↔ gene_name:
# Make chimp tx -> gene_name map
awk -F'\t' '
  $3=="transcript" {
    if (match($0,/transcript_id "([^"]+)"/,m1) && match($0,/gene_name "([^"]+)"/,m2)) {
      print m1[1] "\t" m2[1]
    }
}' Pan_troglodytes.gtf | sort -u > ${OUT}/chimp_tx2gene.tsv

head ${OUT}/chimp_tx2gene.tsv
# ENSPTRT0000xxxxx    GBP5

# (a) protein_id -> transcript_id (from CDS features)
awk -F'\t' '
  $3=="CDS" {
    pid=""; tid="";
    if (match($0,/protein_id "([^"]+)"/,m1)) pid=m1[1];
    if (match($0,/transcript_id "([^"]+)"/,m2)) tid=m2[1];
    if (pid!="" && tid!="") print pid "\t" tid;
}' "$CHIMP_GTF" | sort -u > ${OUT}/chimp_protein2tx.tsv

# (b) transcript_id -> gene_id and gene_name (from transcript features)
awk -F'\t' '
  $3=="transcript" {
    tid=""; gid=""; gname="";
    if (match($0,/transcript_id "([^"]+)"/,m1)) tid=m1[1];
    if (match($0,/gene_id "([^"]+)"/,m2))      gid=m2[1];
    if (match($0,/gene_name "([^"]+)"/,m3))  gname=m3[1];
    if (tid!="") print tid "\t" gid "\t" gname;
}' "$CHIMP_GTF" | sort -u > ${OUT}/chimp_tx2gene.tsv

# quick sanity checks
head ${OUT}/chimp_protein2tx.tsv
head ${OUT}/chimp_tx2gene.tsv


# 5) For each gene in your gene list, extract representative CDS for Human/Nean/Chimp
# Strategy: choose the longest protein-coding transcript per gene to keep orthology simple.
# Skip a header if present; strip blank lines & DOS endings
# Clean gene list: drop header/blank lines & DOS endings
set -Eeuo pipefail
# ==== Prep: ensure folders & helper maps exist ====
# Human transcript->gene map (used to list human transcripts per ENSG)
if [ ! -s "${OUT}/tx2gene.tsv" ]; then
  awk -F'\t' '$3=="transcript"{
    if (match($0,/transcript_id "([^"]+)"/,m1) &&
        match($0,/gene_id "([^"]+)"/,m2))
      print m1[1] "\t" m2[1];
  }' "${HUMAN_GTF}" | sort -u > "${OUT}/tx2gene.tsv"
fi

# Chimp protein->transcript map (for ENSPTRP* -> ENSPTRT*)
if [ ! -s "${OUT}/chimp_protein2tx.tsv" ]; then
  awk -F'\t' '$3=="CDS"{
    pid=""; tid="";
    if (match($0,/protein_id "([^"]+)"/,m1)) pid=m1[1];
    if (match($0,/transcript_id "([^"]+)"/,m2)) tid=m2[1];
    if (pid!="" && tid!="") print pid "\t" tid;
  }' "${CHIMP_GTF}" | sort -u > "${OUT}/chimp_protein2tx.tsv"
fi

# Chimp transcript->(gene_id, gene_name) map (fallback by gene)
if [ ! -s "${OUT}/chimp_tx2gene.tsv" ]; then
  awk -F'\t' '$3=="transcript"{
    tid=""; gid=""; gname="";
    if (match($0,/transcript_id "([^"]+)"/,m1)) tid=m1[1];
    if (match($0,/gene_id "([^"]+)"/,m2))      gid=m2[1];
    if (match($0,/gene_name "([^"]+)"/,m3))  gname=m3[1];
    if (tid!="") print tid "\t" gid "\t" gname;
  }' "${CHIMP_GTF}" | sort -u > "${OUT}/chimp_tx2gene.tsv"
fi

CHIMP_P2T="${OUT}/chimp_protein2tx.tsv"
CHIMP_T2G="${OUT}/chimp_tx2gene.tsv"

# pal2nal path (adjust if needed)
PAL2NAL="./pal2nal.pl"

# ==== Main loop over ENSG IDs ====
awk 'NF && $1!="gene_id"{print $1}' "${GENE_LIST}" | sed 's/\r$//' | while read -r GENE; do
  echo "[*] Processing ${GENE}"

  # 1) all HUMAN transcripts for this gene (from tx2gene)
  awk -v g="${GENE}" '$2==g{print $1}' "${OUT}/tx2gene.tsv" > "${OUT}/human/${GENE}.tx"
  if [ ! -s "${OUT}/human/${GENE}.tx" ]; then
    echo "  - No human transcripts listed in tx2gene for ${GENE}" | tee -a "${OUT}/logs/missing_human_tx.log"
    continue
  fi

  # 2) pull those transcripts from HUMAN CDS
  seqkit grep -n -f "${OUT}/human/${GENE}.tx" "${HUMAN_CDS}" > "${OUT}/human/${GENE}.all.fa" || true
  if [ ! -s "${OUT}/human/${GENE}.all.fa" ]; then
    echo "  - No human CDS records found for ${GENE}" | tee -a "${OUT}/logs/missing_human_tx.log"
    continue
  fi

  # 3) choose longest HUMAN transcript
  seqkit fx2tab -nl "${OUT}/human/${GENE}.all.fa" | sort -k2,2nr | head -n1 | cut -f1 > "${OUT}/human/${GENE}.id"
  TX=$(cat "${OUT}/human/${GENE}.id")
  seqkit grep -n -p "${TX}" "${OUT}/human/${GENE}.all.fa" > "${OUT}/human/${GENE}.fa"
  if [ ! -s "${OUT}/human/${GENE}.fa" ]; then
    echo "  - Representative human transcript sequence missing for ${GENE}" | tee -a "${OUT}/logs/missing_human_tx.log"
    continue
  fi

  # 4) NEANDERTAL: find the SAME transcript ID in Nean CDS
  seqkit grep -n -p "${TX}" "${NEAN_CDS}" > "${OUT}/nean/${GENE}.fa" || true
  if [ ! -s "${OUT}/nean/${GENE}.fa" ]; then
    echo "  - No Neandertal match for ${TX} (${GENE})" | tee -a "${OUT}/logs/missing_nean_tx.log"
    continue
  fi

  # 5) CHIMP resolution from orthology map (7-col BioMart export)
  #    Prefer col6 (chimp transcript ENSPTRT* or protein ENSPTRP*), else fallback by chimp gene (col4/col5)
  line=$(awk -v g="${GENE}" 'BEGIN{FS=OFS="\t"} NR>1 && $1==g{print; exit}' "${ORTHMAP}")
  if [ -z "$line" ]; then
    echo "  - No row in orth map for ${GENE}" | tee -a "${OUT}/logs/missing_chimp.log"
    continue
  fi

  chimp_token=$(awk 'BEGIN{FS=OFS="\t"} {print $6}' <<< "$line")   # ENSPTRT* or ENSPTRP* or empty
  chimp_geneid=$(awk 'BEGIN{FS=OFS="\t"} {print $4}' <<< "$line")  # ENSPTRG...
  chimp_sym=$(awk 'BEGIN{FS=OFS="\t"} {print $5}' <<< "$line")     # symbol

  CHIMP_TX=""

  if [[ "$chimp_token" =~ ^ENSPTRT ]]; then
    CHIMP_TX="$chimp_token"
  elif [[ "$chimp_token" =~ ^ENSPTRP ]]; then
    CHIMP_TX=$(awk -v p="$chimp_token" '$1==p{print $2; exit}' "${CHIMP_P2T}")
  fi

  # Fallback: by chimp gene (ID then symbol), pick longest transcript
  if [ -z "$CHIMP_TX" ]; then
    rm -f "${OUT}/chimp/${GENE}.tx"
    if [ -n "$chimp_geneid" ]; then
      awk -v gid="$chimp_geneid" '$2==gid{print $1}' "${CHIMP_T2G}" > "${OUT}/chimp/${GENE}.tx"
    fi
    if [ ! -s "${OUT}/chimp/${GENE}.tx" ] && [ -n "$chimp_sym" ]; then
      awk -v s="$chimp_sym" '$3==s{print $1}' "${CHIMP_T2G}" > "${OUT}/chimp/${GENE}.tx"
    fi

    if [ -s "${OUT}/chimp/${GENE}.tx" ]; then
      seqkit grep -n -f "${OUT}/chimp/${GENE}.tx" "${CHIMP_CDS}" > "${OUT}/chimp/${GENE}.all.fa" || true
      if [ -s "${OUT}/chimp/${GENE}.all.fa" ]; then
        seqkit fx2tab -nl "${OUT}/chimp/${GENE}.all.fa" | sort -k2,2nr | head -n1 | cut -f1 > "${OUT}/chimp/${GENE}.id"
        CHIMP_TX=$(cat "${OUT}/chimp/${GENE}.id")
      fi
    fi
  fi

  if [ -z "$CHIMP_TX" ]; then
    echo "  - Could not resolve chimp transcript for ${GENE} (token=${chimp_token}, geneid=${chimp_geneid}, sym=${chimp_sym})" \
      | tee -a "${OUT}/logs/missing_chimp_seq.log"
    continue
  fi

  # Fetch the chosen chimp transcript CDS
  seqkit grep -nrp "${CHIMP_TX}" "${CHIMP_CDS}" > "${OUT}/chimp/${GENE}.fa" || true
  if [ ! -s "${OUT}/chimp/${GENE}.fa" ]; then
    echo "  - Chimp CDS not found for ${GENE} (tx=${CHIMP_TX})" | tee -a "${OUT}/logs/missing_chimp_seq.log"
    continue
  fi

  # 6) Translate AA (frame check)
  seqkit translate -f 1 "${OUT}/human/${GENE}.fa" -o "${OUT}/aa/${GENE}.human.aa.fa" || continue
  seqkit translate -f 1 "${OUT}/nean/${GENE}.fa"  -o "${OUT}/aa/${GENE}.nean.aa.fa"  || continue
  seqkit translate -f 1 "${OUT}/chimp/${GENE}.fa" -o "${OUT}/aa/${GENE}.chimp.aa.fa" || continue

  # 7) Protein align, back-translate to codons
  cat "${OUT}/aa/${GENE}.human.aa.fa" \
      "${OUT}/aa/${GENE}.nean.aa.fa"  \
      "${OUT}/aa/${GENE}.chimp.aa.fa" > "${OUT}/aa/${GENE}.aa.fa"

  mafft --auto "${OUT}/aa/${GENE}.aa.fa" > "${OUT}/aln/${GENE}.aa.aln.fa" 2>> "${OUT}/logs/mafft.log" || continue

  cat "${OUT}/human/${GENE}.fa" \
      "${OUT}/nean/${GENE}.fa"  \
      "${OUT}/chimp/${GENE}.fa" > "${OUT}/cds/${GENE}.cds.fa"

  # after you create ${OUT}/aln/${GENE}.aa.aln.fa and ${OUT}/cds/${GENE}.cds.fa
TMP=$(mktemp)
if "${PAL2NAL}" "${OUT}/aln/${GENE}.aa.aln.fa" "${OUT}/cds/${GENE}.cds.fa" -output paml \
    > "$TMP" 2>> "${OUT}/logs/pal2nal.log"; then
  mv -f "$TMP" "${OUT}/codeml/${GENE}.codon.phy"
else
  echo "[PAL2NAL FAIL] ${GENE}" | tee -a "${OUT}/logs/pal2nal_fail.list"
  rm -f "$TMP"
  continue
fi

  # 8) Tree for codeml (Neanderthal sister to Human)
  echo "(chimp,(human,nean));" > "${OUT}/tree/${GENE}.nwk"

done

```

# C) dN/dS with PAML/codeml
# Below is a codeml control file template (codeml.ctl) to test branch model with the human lineage foreground (or test the human-Neandertal split). You’ll run codeml per gene.

# Make a foreground-labeled tree if you want branch-specific ω (model=2). For example, mark the human branch with #1:
# You can adjust parameters as needed
# seqfile = input codon alignment
# treefile = input tree
# outfile = output file

```bash
# Build a default tree for every .codon.phy that lacks one
shopt -s nullglob
for CODON in "$OUT"/codeml/*.codon.phy; do
  G=$(basename "$CODON" .codon.phy)
  T="$OUT/tree/${G}.nwk"
  [ -s "$T" ] || echo "(chimp,(human,nean));" > "$T"
done
```

# Run codeml for each gene
```bash
for CODON in "$OUT"/codeml/*.codon.phy; do
  G=$(basename "$CODON" .codon.phy)
  TREE="$OUT/tree/${G}.nwk"
  FG_TREE="$OUT/tree/${G}.fg.nwk"
  OUTFILE="$OUT/codeml/${G}.codeml.out"
  CTL="$OUT/codeml/${G}.ctl"
  LOG="$OUT/logs/codeml_${G}.log"

  # skip empty/missing inputs
  if [ ! -s "$CODON" ]; then
    echo "[WARN] Missing/empty codon file for $G"
    continue
  fi
  if [ ! -s "$TREE" ]; then
    echo "[WARN] Missing tree for $G (expected $TREE)"
    continue
  fi

  # mark human foreground
  sed 's/human/human#1/' "$TREE" > "$FG_TREE"

  # write ctl
  cat > "$CTL" <<EOF
seqfile = $CODON
treefile = $FG_TREE
outfile  = $OUTFILE
noisy    = 9
verbose  = 1
runmode  = 0
seqtype  = 1
CodonFreq= 2
clock    = 0
model    = 2
NSsites  = 0
icode    = 0
fix_kappa= 0
kappa    = 2
fix_omega= 0
omega    = 0.2
cleandata= 1
getSE    = 0
EOF

  ( cd "$OUT/codeml" && codeml "$(basename "$CTL")" ) > "$LOG" 2>&1 || {
    echo "[ERR] codeml failed for $G (see $LOG)"
    continue
  }

  echo "[OK] $G"
done
shopt -u nullglob
```








```
      seqfile = your_codon_alignment.phy
     treefile = your_tree.nwk
      outfile = your_codeml_output.txt

        noisy = 9
      verbose = 1
      runmode = 0

      seqtype = 1       * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2       * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

        model = 2       * 0:one w; 1:separate w for branches; 2:user tree
      NSsites = 0       * 0:one w;1:neutral;2:selection;3:discrete;4:freqs;
                        * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                        *10:beta&gamma+1;11:beta&gamma+2;12:0&2normal>1;
                        *13:3normal>1;14:3normal>2;15:cladeC

        icode = 0       * 0:universal code; 1:mammalian mt; 2-10:see below
    fix_kappa = 0       * 1:kappa fixed, 0:kappa to be estimated
        kappa = 2       * initial or fixed kappa

    fix_omega = 0       * 1: omega or omega_1 fixed, 0: estimate
        omega = 1       * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1       * do not estimate alpha
        alpha = 0.       * initial or fixed alpha, 0:infinity (constant rate)

       Malpha = 0       * different alphas for genes

        ncatG = 10      * # of categories in the dG or AdG models of rates

        getSE = 0       * do not want them, too time consuming
 RateAncestor = 0       * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6   * default value for Newton-R
   cleandata = 1       * remove sites with ambiguity data (1:yes, 0:no)?
```
