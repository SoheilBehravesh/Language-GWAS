######Finding common and unique genes among multiple populations###

# Finding common genes overlapping phenotype-associated SNPs among two populations (East Asians and Europeans)
# For Ngyut
common.genes <- intersect(sigEAS.genesNAME[[1]], sigEUR.genesNAME[[1]])
# If none, this could have interesting implications that diabetes has different genetic mechanisms for East Asians and Europeans
# Want to focus on genes overlapping phenotype-associated SNPs in East Asians? (Indonesians?)
# Make table with columns with number of genes significant in each population, and in common


# Finding common genes overlapping phenotype-associated SNPs among three populations (Africans, Central-South Asians, and Europeans)
# For Adele (check names of gene name file, I was working from memeory haha)
# Make vectors
afr <- sigAFR.genesNAME[[1]]
csa <- sigCSA.genesNAME[[1]]
eur <- sigEUR.genesNAME[[1]]

# Common to AFR, CSA, and EUR
common_all_three <- Reduce(intersect, list(afr, csa, eur))

# Genes in any of the pops
all <- union(afr, csa, eur)

# Common to any two but not all three
common_afr_csa <- intersect(afr, csa)
common_afr_eur <- intersect(afr, eur)
common_csa_eur <- intersect(csa, eur)
# Combine pairs, remove those common to all three
common_two_only <- unique(c(common_afr_csa, common_afr_eur, common_csa_eur))
common_two_only <- setdiff(common_two_only, common_all_three)

#Unique to each group
unique_afr <- setdiff(afr, union(csa, eur))
unique_csa <- setdiff(csa, union(afr, eur))
unique_eur <- setdiff(eur, union(afr, csa))