micromamba activate bioinfo
set -uex -o pipefail

cd ~/work/language_gwas/dn_ds
# A) Files to download
cd ~/work/language_gwas
mkdir -p data
cd data
# Download human reference genome FASTA
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz > homo_sapiens.dna.fa
# Download human reference genome GTF
wget https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
gunzip Homo_sapiens.GRCh38.115.gtf.gz > homo_sapiens.dna.gtf
# Download Neandertal VCF from max planck
base_url="http://ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF"
for chr in {1..22} X; do
    wget "${base_url}/chr${chr}.noRB.vcf.gz.tbi"
done
# Download Chimp reference CDS FASTA (or cDNA) and ortholog mapping
wget https://ftp.ensembl.org/pub/current_fasta/pan_troglodytes/cdna/Pan_troglodytes.Pan_tro_3.0.cdna.all.fa.gz
gunzip Pan_troglodytes.Pan_tro_3.0.cdna.all.fa.gz > pan_troglodytes.cdna.fa
# For making the ortholog map table use the R code in the ortholog_mapping.Rmd and the name of the output file is "human_chimp_orthologs.csv"
# 119 language-linked candidates (Ensembl Gene IDs or gene symbols)

# Tools needed
# bcftools, samtools, gffread (from gffread / Cufflinks), seqkit
# Multiple sequence alignment: mafft
# Back-translation for codons: pal2nal.pl (PAL2NAL) or codon aligner (MACSE/PRANK)
# dN/dS: PAML/codeml (or HyPhy, but I’ll show PAML here)


# B) Unix pipeline to build codon alignments per gene
