############################################
# Visualizing GWAS and exploring data
############################################

# Steps in Terminal:
# ex. SNPs associated with stomach cancer downloaded from the pan-UKBB
# using wget https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/categorical-20001-both_sexes-1018.tsv.bgz
# gunzip –c categorical-20001-both_sexes-1018.tsv.bgz > categorical-20001-both_sexes-1018.tsv

# Steps in RStudio

#Set your working directory
#setwd("~/Downloads")

##### Installing and loading packages
#install.packages("qqman")
#install.packages("data.table")

# These packages help us clean up our data
library(dplyr)
library(tidyr)
library(data.table)

# These packages help us visualize our data 
library(ggplot2)
library(qqman)

##### Reading in the data and QC
# May take a few minutes
#df <- read.table("categorical-20001-both_sexes-1018.tsv", sep = '\t', header = T)

# List all 8 TSV file paths (update with your actual file paths or use list.files)
files <- c("gwas-association-downloaded_2025-06-04-EFO_0003926.tsv", 
           "gwas-association-downloaded_2025-06-04-HP_0025268.tsv",
           "gwas-association-downloaded_2025-06-04-pubmedId_21051773.tsv",
           "gwas-association-downloaded_2025-06-04-pubmedId_23423138.tsv", 
           "gwas-association-downloaded_2025-06-04-pubmedId_30741946.tsv",
           "gwas-association-downloaded_2025-06-04-pubmedId_36266505.tsv", 
           "gwas-association-downloaded_2025-06-04-pubmedId_34861174.tsv", 
           "gwas-association-downloaded_2025-06-04-HP_0100033.tsv")

# Read and combine them
df1 <- do.call(rbind, lapply(files, function(f) read.delim(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)))
colnames(df1)
nrow(df1)
ncol(df1)

df <- read.table("GCST90297560.tsv", sep = '\t', header = T)
colnames(df)

df <- read.table("GCST90558312.tsv", sep = '\t', header = T, fill = TRUE)
colnames(df)

# Looking at the first two rows of our data
head(df, 2)

# Sometimes a file format will have multiple types of data mashed into one column and we need to separate them
# ex. Separate column 1 into 4 individual columns
# Convert to data.table (if not already)
# setDT(df)
# Split the 'variant' column
# df[, variant := as.character(variant)]
# df[, c("Chr", "Pos", "Allele_1", "Allele_2") := tstrsplit(variant, ":", fixed = TRUE)]


# Different data sources will give their data's columns slightly different names. Here we will rename
# each column so that the script recognizes the data properly.

# colnames(df)[?] <- "Chr"
# colnames(df)[?] <- "Pos"
# colnames(df)[?] <- "pval"
# colnames(df)[?] <- "Beta"
# colnames(df)[?] <- "minor_AF"


# View your dataframe (df).
# Some columns are easy to identify and just have slightly different names, but others are more complicated.
# ex.chr, pos, and beta_EUR are easy to identify as Chr, Pos, and Beta, but here the Pval and minor_AF are
# harder to identify.

colnames(df)[1] <- "Chr"
colnames(df)[2] <- "Pos"
colnames(df)[7] <- "Beta"

# This dataset gives the p-value as "neglog10_pval_Eur" so we have to reverse the neglog10 to get the p-value.
df$pval <- 10^(-(df$neglog10_pval_EUR))

# The dataset does not give a minor allele frequency, only the effect (alternate) allele frequency for the control and trait
colnames(df)[5] <- "effect_AF"

# Make sure that AF column is numeric
df$effect_AF <- as.numeric(df$effect_AF)

# Filter to keep only those variants with AF > 0.05
df <- df %>% filter(effect_AF > 0.05)

# Look at your new column names
head(df, 2)

# Filter to keep Chr X in the chromosome column 
length(which(df$Chr == "X")) # Checking how many rows (SNPs) are on Chr X

df <- df %>% filter(Chr != "X")


# Make sure other columns are numeric 
df$pval <- as.numeric(df$pval)
df$Chr <- as.numeric(df$Chr)
df$Pos <- as.numeric(df$Pos)

# Make ID column
df$ID <- paste(df$Chr, df$Pos, sep = ":")

##### Generating Manhattan Plot
# Selecting the 4 columns required for a Manhattan plot
man <- select(df, Chr, Pos, pval, ID) # chromosome, position, pvalue, SNP ID

# Renaming the columns to the headers recognized by the Manhattan function
colnames(man)[1] <- "CHR" 
colnames(man)[2] <- "BP"
colnames(man)[3] <- "P"
colnames(man)[4] <- "SNP"

# Creating Manhattan plot and assigning it to an object called manhattan 1
manhattan1 <- manhattan(man)
manhattan1 # run this to see your plot

# we can also sort the chromosomes by color and edit our plot
manhattan2 <- manhattan(man, main = "Stomach Cancer",
                        cex = 0.8, # cex changes the size of the points
                        col = c("cadetblue2", "darkslategray4"))
manhattan2
# May take upwards of 13 minutes to generate

##### Exploring our data 
# Counting how many SNPs are significantly associated with the phenotype 
length(which(df$pval <= 5e-8)) # genome-wide threshold, note how many SNPs 
length(which(df$pval <= 1e-5)) # suggestive threshold 

sig <- df %>% filter(pval <= 5e-8)

write.table(sig, "significant_SNPs_StomachCancer.txt", sep ='\t', row.names = F, quote = F)

# We can plot the beta (effect size) against p-value
ggplot(sample, aes(x=abs(beta), y=-log10(P-value))) +
  geom_point(color = "#0c4c8a") +
  theme_minimal() +
  labs(title = "TITLE", x = "Effect size", y = "-log10(P-value)") + 
  theme(plot.title = element_text(hjust=0.5, size=15, face = "bold"))

# If you're working with two populations you can compare the effect size of significant SNPs
sig <- sample %>% filter(P_EAS <= 1e-5 & P_EUR <= 1e-5) # pulling the SNPs that are significant in BOTH populations

ggplot(sig, aes(x=Beta_EAS, y=Beta_EUR)) +
  geom_point(color = "darkslategray4") +
  theme_minimal() +
  labs(title = "Two Pop Comparison", x = "EAS effect size", y = "EUR effect size") + 
  theme(plot.title = element_text(hjust=0.5, size=15, face = "bold"))

# If you're working with two phenotypes you can compare the effect size of significant SNPs
merged <- merge(pheno1, pheno2, by=c("ID")) # first you would have to merge the two data sets, either by ID (or chr AND position)

# You would then change the column names to keep track of the information for each phenotype
colnames(sig)[] <- "pheno1_pvalue" # renaming the pvalue column for phenotype 1
colnames(sig)[] <- "pheno2_pvalue" # renaming the pvalue column for phenotype 2
colnames(sig)[] <- "pheno1_beta" # renaming the Beta column for phenotype 1
colnames(sig)[] <- "pheno2_beta" # renaming the Beta column for phenotype 2
# repeat this for all of the columns 

# Plotting for two phenotype comparison
sig <- sample %>% filter(pheno1_pvalue <= 1e-5 & pheno2_pvalue <= 1e-5) # pulling the SNPs that are significant in BOTH populations

ggplot(sig, aes(x=pheno1_beta, y=pheno2_beta)) +
  geom_point(color = "darkslategray4") +
  theme_minimal() +
  labs(title = "Two Pop Comparison", x = "Pheno1 effect size", y = "Pheno2 effect size") + 
  theme(plot.title = element_text(hjust=0.5, size=15, face = "bold"))




