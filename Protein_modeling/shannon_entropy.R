rm(list=ls())
install.packages("bio3d", dependencies=TRUE)


library(msa)
library(seqinr)
library(bio3d)

setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/adaptive_sweep")

alpha_in_file <- readAAStringSet("alpha_entire_AA_locus2.fa")
beta_in_file <- readAAStringSet("beta_entire_AA_locus2.fa")

#align
alpha_align <- msa(alpha_in_file, "ClustalOmega", order = "input")
beta_align <- msa(beta_in_file, "ClustalOmega", order = "input")

#convert to matrix
alpha_aln_matrix <- as.matrix(alpha_align)
beta_aln_matrix <- as.matrix(beta_align)

# calculate entropy
alpha_h <- entropy(alpha_aln_matrix)
beta_h <- entropy(beta_aln_matrix)

#entropy table
alpha_entropy_values <- data.frame(alpha_normalized_entropy = (alpha_h$H/max(alpha_h$H)))
beta_entropy_values <- data.frame(beta_normalized_entropy = (beta_h$H/max(beta_h$H)))

# entropy values with names
alpha_con <- consensus(alpha_aln_matrix)
names(alpha_h$H) = alpha_con$seq
alpha_entropy_values_with_residues <- c(alpha_h$H)
write.csv(alpha_entropy_values, "alpha_entropy_values.csv")


beta_con <- consensus(beta_aln_matrix)
names(beta_h$H) = beta_con$seq
beta_entropy_values_with_residues <- c(beta_h$H)
write.csv(beta_entropy_values, "beta_entropy_values.csv")
