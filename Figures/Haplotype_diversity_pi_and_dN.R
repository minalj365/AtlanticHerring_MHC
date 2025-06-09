############# Haplotype diversity ################

######### pi #########
setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/2024/haplotype_diversity")
samples=c("Reference", "CS2_hap1", "CS2_hap2", "CS4_hap1", "CS4_hap2", "CS5_hap1", "CS5_hap2", "CS7_hap1", "CS7_hap2", "CS8_hap1", "CS8_hap2", "CS10_hap1", "CS10_hap2", "BS1_hap1", "BS1_hap2", "BS2_hap1", "BS2_hap2", "BS3_hap1", "BS3_hap2", "BS4_hap1", "BS4_hap2", "BS5_hap1", "BS5_hap2", "BS6_hap1", "BS6_hap2", "NSSH2_hap1", "NSSH2_hap2", "NSSH10_hap1", "NSSH10_hap2")
alpha_samples=c("Reference", "CS2_hap1", "CS2_hap2", "CS4_hap1", "CS4_hap2", "CS5_hap2", "CS7_hap1", "CS7_hap2", "CS8_hap1", "CS8_hap2", "CS10_hap1", "CS10_hap2", "BS1_hap1", "BS1_hap2", "BS2_hap1", "BS2_hap2", "BS3_hap1", "BS3_hap2", "BS4_hap1", "BS4_hap2", "BS5_hap1", "BS5_hap2", "BS6_hap1", "BS6_hap2", "NSSH2_hap1", "NSSH2_hap2", "NSSH10_hap1", "NSSH10_hap2")

## save aligned sequences ##
for(i in alpha_samples){
  # read in the sequence file for current gene and exon
  in_file <- readDNAStringSet(paste0("Sequences/Raw/", i, "_alpha.fa"))
  d <- as.dist(pwalign::stringDist(in_file, ))
  align <- msa(in_file, "ClustalOmega", order = "input")
  Biostrings::writeXStringSet(DNAStringSet(align), paste0("Sequences/Aligned/", i, "_alpha.fa"))
}

for(i in samples){
  # read in the sequence file for current gene and exon
  in_file <- readDNAStringSet(paste0("Sequences/Raw/", i, "_beta.fa"))
  d <- as.dist(pwalign::stringDist(in_file, ))
  align <- msa(in_file, "ClustalOmega", order = "input")
  Biostrings::writeXStringSet(DNAStringSet(align), paste0("Sequences/Aligned/", i, "_beta.fa"))
}


### calculate pi ###
# beta
beta_pi_df <- data.frame()
for(i in samples){
# read in the sequence file for current gene and exon
in_file <- readDNAStringSet(paste0(i, "_beta.fa"))
d <- as.dist(pwalign::stringDist(in_file, ))
align <- msa(in_file, "ClustalOmega", order = "input")
CDS_pi <- round(nuc.div(as.DNAbin(align)), digits = 3)
# assign pi value to appropriate cell in dataframe
 beta_pi_df[i, "beta_pi"] <- CDS_pi
 }

alpha_pi_df <- data.frame()
for(i in alpha_samples){
  # read in the sequence file for current gene and exon
  in_file <- readDNAStringSet(paste0(i, "_alpha.fa"))
  d <- as.dist(pwalign::stringDist(in_file, ))
  align <- msa(in_file, "ClustalOmega", order = "input")
  CDS_pi <- round(nuc.div(as.DNAbin(align)), digits = 3)
  # assign pi value to appropriate cell in dataframe
  alpha_pi_df[i, "alpha_pi"] <- CDS_pi
}

### plot pi for haplotypes ##
data <- read_tsv("hapl_div.tsv", col_names = T)
data_long <- pivot_longer(data, cols = c(alpha_pi, beta_pi), names_to = "group", values_to = "pi") %>%
  drop_na(pi)

haplotype_pi <- ggplot(data_long, aes(x = group, y = pi)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "#990000", fatten = 0) +
  annotate("segment", x = 0.8, xend = 1.2, y = 0.064, yend = 0.064, color = "#0066CC", linewidth = 0.4) +
  annotate("segment", x = 1.8, xend = 2.2, y = 0.083, yend = 0.083, color = "#0066CC", linewidth = 0.4) +
  coord_cartesian(ylim = c(0.06, 0.13)) + 
  theme_classic() +
  scale_x_discrete(labels = c("alpha", "beta")) +
  labs(x = "Locus  2 genes",
       y = expression(paste("Nucleotide diversity (", italic("\u03C0"), ")")),
       tag = "a") 

########## similar plot with dN values ###########
## first, towards getting average dN value  (to be plotted in blue)
# align
  # read in the sequence file for current gene and exon
  in_file <- readDNAStringSet("Sequences/Raw/all_alpha_not_54.fa")
  d <- as.dist(pwalign::stringDist(in_file, ))
  align <- msa(in_file, "ClustalOmega", order = "input")
  Biostrings::writeXStringSet(DNAStringSet(align), "Sequences/Aligned/all_alpha_not_54.fa")

  
  in_file <- readDNAStringSet("Sequences/Raw/all_beta_not_60.fa")
  d <- as.dist(pwalign::stringDist(in_file, ))
  align <- msa(in_file, "ClustalOmega", order = "input")
  Biostrings::writeXStringSet(DNAStringSet(align), "Sequences/Aligned/all_beta_not_60.fa")

  
## With MEGA, calculate dN. Output is:
  # dN for alpha = 0.0624; dN for beta = 0.071
## plot ##
data <- read_tsv("dN.tsv", col_names = T)

haplotype_dN <- ggplot(data, aes(x = group, y = dN)) +
  geom_jitter(width = 0.1, alpha = 0.6) +
  stat_summary(fun = mean, geom = "errorbar",  aes(ymin = after_stat(y), ymax = after_stat(y)), width = 0.4, color = "#990000", ) +
  annotate("segment", x = 0.8, xend = 1.2, y = 0.0624, yend = 0.0624, color = "#0066CC", linewidth = 0.4) +
  annotate("segment", x = 1.8, xend = 2.2, y = 0.071, yend = 0.071, color = "#0066CC", linewidth = 0.4) +
  theme_classic() +
  scale_x_discrete(labels = c("alpha", "beta")) +
  labs(x = "Locus  2 genes",
       y = expression(paste("Non-synonymous substitution rate (", italic("dN"), ")"))) 

ggsave("dN_hapl_div.png", haplotype_dN, width = 4, height = 4, dpi = 300)
ggsave("dn_hapl_div.pdf", haplotype_dN, width = 4, height = 4, dpi = 300)




#haplotype_pi/haplotype_dN
hapl_div <- plot_grid(haplotype_pi, haplotype_dN, nrow=1)
# Fig5 <- alpha_pi_plot / beta_pi_plot    ## not working! Forgot which package

ggsave("hapl_div.png", hapl_div, width = 7.08, height = 4, dpi = 300)
ggsave("hapl_div.pdf", hapl_div, width = 7.08, height = 4, dpi = 300)


######## plot for all sequences (54 and 49 alpha AND 60 and 56 beta) #####


  # read in the sequence file for current gene and exon
  in_file <- readDNAStringSet("Sequences/Raw/all_alpha.fa")
  d <- as.dist(pwalign::stringDist(in_file, ))
  align <- msa(in_file, "ClustalOmega", order = "input")
  #Biostrings::writeXStringSet(DNAStringSet(align), paste0("Sequences/Aligned/", i, "_alpha.fa"))
  CDS_pi <- round(nuc.div(as.DNAbin(align)), digits = 3)
  print(CDS_pi)  # 0.066


# read in the sequence file for current gene and exon
  in_file <- readDNAStringSet("Sequences/Raw/all_beta.fa")
  d <- as.dist(pwalign::stringDist(in_file, ))
  align <- msa(in_file, "ClustalOmega", order = "input")
  CDS_pi <- round(nuc.div(as.DNAbin(align)), digits = 3)
  print(CDS_pi)  #0.083

  
########## not 54 and not 60 ##########
  # read in the sequence file for current gene and exon
  in_file <- readDNAStringSet("Sequences/Raw/all_alpha_not_54.fa")
  d <- as.dist(pwalign::stringDist(in_file, ))
  align <- msa(in_file, "ClustalOmega", order = "input")
  #Biostrings::writeXStringSet(DNAStringSet(align), paste0("Sequences/Aligned/", i, "_alpha.fa"))
  CDS_pi <- round(nuc.div(as.DNAbin(align)), digits = 3)
  print(CDS_pi) ## 0.064
  
  
  # read in the sequence file for current gene and exon
  in_file <- readDNAStringSet("Sequences/Raw/all_beta_not_60.fa")
  d <- as.dist(pwalign::stringDist(in_file, ))
  align <- msa(in_file, "ClustalOmega", order = "input")
  CDS_pi <- round(nuc.div(as.DNAbin(align)), digits = 3)
  print(CDS_pi) ##0.083
