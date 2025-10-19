######################################################################################
# Investigating evolutionary history of trait-associated SNPs
# -- testing whether regions of the genome significantly associated with trait 
# -- overlap candidate selection regions
######################################################################################

# Installing and loading packages
install.packages("ivs")

library(dplyr)
library(ivs)

# Reading in data 
sig <- read.table("Phenotype_sigSNPS.txt", sep = '\t', header = T)
ihs <- read.table("EUR.IHSregions.50kbp.txt", sep = '\t', header = T) # candidate selection regions iHS - past 25,000 years

# Generating regions associated with trait to run the evolutionary analysis 
sig$Chr <- as.numeric(sig$Chr)
sig$Pos <- as.numeric(sig$Pos)

# Create a window around each SNP with 25000 bp to either side
sig$win_start <- sig$Pos - 25000 
sig$win_end <- sig$Pos + 25000

length(which(sig$Chr == 22)) # checking to see which chromosomes are represented in the significant SNPs

# Go through each chromosome to create regions from the overlapping windows
# Instead of doing it manually for each chromosome , we can create a loop to generate the regions

for (i in 1:22) { #1:22 because we're working with 22 chromosomes 
  
  chr <- sig[sig$Chr %in% paste(i, sep = ""),] #pull snps in chr "i" into a dataframe called chr
  
  if (i %in% c(1,4,7,12:22)){ # allows loop to continue if we're missing a chromsome (for ex ch21)
    next
  }
  
  # create range from windows
  chr <- chr %>%
    mutate(Range = iv(win_start, win_end), .keep = "unused")
  
  # look for overlapping windows and creating regions, output is in list format
  regions <- iv_groups(chr$Range, abutting = T)
  
  # turn list output into data frame
  regions <- as.data.frame(regions)
  x <- data.frame(Reduce(rbind, regions$regions))
  
  # remove row names 
  rownames(x) <- NULL 
  
  # unlist columns 1 and 2 
  x[1:2] <- lapply(x, unlist) 
  
  # create a new column with chromosome number 
  x$chr <- paste(i, sep = "")
  
  # change name of data frame "x" to the chromosome "i"
  assign(paste('ch', i, sep = ""), x)
  
} # this will repeat 22 times and generate a data frame with chr, start of region, end of region

colnames(ch9)[1] <- "start" #use this to add column headers to those chr with only 1 region

regions <- rbind(ch2, ch3, ch5, ch6, ch8, ch9, ch10, ch11) # bind the regions for each chromosome into 1 data frame

# write out file with regions 
write.table(regions, "sample_regions.txt", sep = '\t', row.names = F, quote = F)

# clean your global environment
rm(list=ls(pattern = "ch"), x)

# Preparing data frames for permutations
regions$start <- as.numeric(regions$start)
regions$end <- as.numeric(regions$end)
regions$chr <- as.numeric(regions$chr)

ihs$start <- as.numeric(ihs$start)
ihs$end <- as.numeric(ihs$end)
ihs$chr <- as.numeric(ihs$chr)

sample$Pos <- as.numeric(sample$Pos)
sample$Chr <- as.numeric(sample$Chr)

regions$ihs_regions <- ""

for (i in 1:nrow(regions)) {
  y <- which((regions$start[i] <= ihs$end)&(ihs$start <= regions$end[i])&(regions$chr[i] == ihs$chr))
  if (length(y)>=1) {
    regions$ihs_regions[i] <- print("1")
  } else {regions$ihs_regions[i] <- print("0")}
}

obs.overlaps <- length(which(regions$ihs_regions == 1)) # 8

# Creating functions
regions$width <- regions$end - regions$start # creating new column with the width of the region

# this function creates new random regions
new.GWAS.regions <- function(snp) { 
  
  for (k in 1:nrow(snp)) {
    snp$end <- snp$Pos
    snp$start <- snp$Pos - regions$width
    GWAS.regions <<- select(snp, Chr, start, end)
  }
  
  GWAS.regions <<- GWAS.regions
  
}

# this function looks for overlaps betweent the random regions and candidate selection regions
overlaps <- function(GWAS.regions, ihs) {
  
  GWAS.regions$ihs_regions <- rep(NA)
  
  for (i in 1:nrow(GWAS.regions)) {
    y <- which((GWAS.regions$start[i] <= ihs$end)&(ihs$start <= GWAS.regions$end[i])&(GWAS.regions$Chr[i] == ihs$chr))
    if (length(y)>=1) {
      GWAS.regions$ihs_regions[i] <- print("1")
    } else {GWAS.regions$ihs_regions[i] <- print("0")}
  }
  
  GWAS.regions <<- GWAS.regions
}
# make sure Chr and Pos have appropriate name in the functions 

# Permutation analysis to see if our significant GWAS regions have more overlaps than expected by chance
PermutationTable <- rep(NA, 5000)

for (j in 1:5000) { 
  
  snp <- sample_n(sample, nrow(regions)) # sampling from our total GWAS data set 
  new.GWAS.regions(snp) # using the function we created to generate new regions from these random samples SNPs
  
  overlaps(GWAS.regions, ihs) # using the function we created to look for overlaps
  
  table <- GWAS.regions %>% filter(ihs_regions == "1") # pull the regions that overlap into a table
  PermutationTable[j] <- nrow(table) # count the number of overlapping regions and add them to permutation table
  print(j)
  
}

PermutationTable <- as.data.frame(PermutationTable)

# Visualizing permutation
CalcPercent <- PermutationTable %>% filter(PermutationTable>=obs.overlaps)
pval <- nrow(CalcPercent)/nrow(PermutationTable) 

Percent <- quantile(PermutationTable$PermutationTable, 0.95)
Percent2 <- quantile(PermutationTable$PermutationTable, 1-pval)

ggplot(PermutationTable, aes(x=PermutationTable, y=)) +
  theme(text=element_text(size=50)) +
  geom_histogram(binwidth = 1) + 
  geom_vline(xintercept=Percent, color = "olivedrab4", linetype="dashed", linewidth = 1.25) +
  geom_point(aes(x=Percent2, y=0), color="olivedrab4", size=5.5) +
  scale_x_continuous(name = "Number of overlapping regions") + 
  scale_y_continuous(name="Number of permutations") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Permutation analysis iHS - EAS only") +
  theme(plot.title = element_text(hjust = 0.5))