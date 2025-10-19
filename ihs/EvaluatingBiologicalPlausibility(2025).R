############################################
# Evaluating biological plausibility
############################################

##### Loading packages
# These packages help us clean up our data
library(dplyr)

# These packages help us visualize our data 
library(ggplot2)

#Set your working directory
setwd("~/Downloads")

##### Gene ontology #####

##### Read in SNPs that are significantly associated with the phenotype (or at least suggestive)
sig <- read.table("suggestive_SNPs_StomachCancer.txt", sep = '\t', header = T)

##### Reading in genes to investigate at biological plausibility 
genes <- read.table("biomaRtprotcoding_genes_GRCh37.txt", sep = '\t', header = T)

head(genes, 2)

colnames(genes)[1] <- "Chr"

# Expanding genes start and stop positions to encompass regulatory regions 
genes <- mutate(genes, start_10000 = (genes$start) - 10000) 
genes <- mutate(genes, end_10000 = (genes$end) + 10000)

head(genes,2)

# Making sure columns are numeric 
genes$start_10000 <- as.numeric(genes$start_10000)
genes$end_10000 <- as.numeric(genes$end_10000)
genes$Chr <- as.numeric(genes$Chr)

sig$Chr <- as.numeric(sig$Chr)
sig$Pos <- as.numeric(sig$Pos)

###### Investigating which genes overlap phenotype-associated SNPs
for(i in 1:nrow(genes)){
  y <- which((sig$Pos >= genes$start_10000[i])&(sig$Pos <= genes$end_10000[i])&(sig$Chr == genes$Chr[i]))
  if(length(y)>=1){
    genes$SNPs[i] <- print("1")
  } else {genes$SNPs[i] <- print("0")}
}

length(which(genes$SNPs == 1))

# New dataframe with the genes that contain at least one 1 SNP
sig.genes <- genes %>% filter(SNPs == "1")

# Taking only the gene ID and gene name
sig.genesID <- dplyr::select(sig.genes, gene_id)
sig.genesNAME <- dplyr::select(sig.genes, gene_name)

# Making sure gene ID and name is read as characters (not numeric)
sig.genesID$gene_id <- as.character(sig.genesID$gene_id)
sig.genesNAME$gene_name <- as.character(sig.genesNAME$gene_name)

# write table with gene IDs for your phenotype
write.table(sig.genesID, "StomachCancer_genes_ID.txt", sep = '\t', row.names = F, quote = F)
write.table(sig.genesNAME, "StomachCancer_genes_NAME.txt", sep = '\t', row.names = F, quote = F)

# In order to view some of the things your genes are associated with, open the "genes_ID" file you just made.
# It has the ENSEMBL IDs for each of the genes you found that ovelap the significantly associated SNPs.

# Go to https://davidbioinformatics.nih.gov/summary.jsp
# On the panel on the right, copy-paste the gene IDs from your file into the "Gene List" section under Step 1
# Under Step 2, choose ENSEMBL_GENE_ID from the dropdown menu
# For Step 3, choose "gene list," then Submit

# You can then click "Functional Annotation Table" to see which genes are associated with what.
# Does anything catch your eye that you think could be related to your trait of interest?
# You can also look up individual genes by copy-pasting the gene ID into the PAN-GO Human Functionome's searchbar at the top of the page here:
# https://functionome.geneontology.org



##### Additional plots

## We can plot the beta (effect size) against p-value
ggplot(df, aes(x=abs(Beta), y=-log10(pval))) +
  geom_point(color = "#0c4c8a") +
  theme_minimal() +
  labs(title = "Controlling for effect size", x = "Effect size", y = "-log10(P-value)") + 
  theme(plot.title = element_text(hjust=0.5, size=15, face = "bold"))

## If you're working with two populations you can compare the effect size of significant SNPs
# We first have to merge the data sets, either by ID (or chr AND position)
merged <- merge(Pop1, Pop2, by=c("ID"))

# You would then change the column names to keep track of the information for each population
colnames(sig)[] <- "pop1_pvalue" # renaming the pvalue column for population 1
colnames(sig)[] <- "pop2_pvalue" # renaming the pvalue column for population 2
colnames(sig)[] <- "pop1_beta" # renaming the Beta column for population 1
colnames(sig)[] <- "pop2_beta" # renaming the Beta column for population 2
# repeat this for all of the columns 

sig <- sample %>% filter(P_Pop1 <= 1e-5 & P_Pop2 <= 1e-5) # pulling the SNPs that are significant in BOTH populations

ggplot(sig, aes(x=Beta_Pop1, y=Beta_Pop2)) +
  geom_point(color = "darkslategray4") +
  theme_bw() +
  labs(title = "Two Pop Comparison", x = "Pop1 effect size", y = "Pop2 effect size") + 
  theme(plot.title = element_text(hjust=0.5, size=15, face = "bold")) +
  theme(
    panel.border = element_blank(), # Remove the plot border
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    axis.line = element_line(color = "black"), # Keep axis lines
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14), # Set axis title size
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold")
  )