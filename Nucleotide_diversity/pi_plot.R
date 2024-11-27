rm(list=ls())
library(msa)
library(seqinr)
library(ape)
library(pegas)
library(phytools)
library(pwalign)
library(ggtree)
library(ggmsa)
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)
library(RColorBrewer)
library(zoo)
#install.packages("tinytex") # latex - to save alignment in pdf file

setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/Sequence_analysis/pi")

#genes=c("DAA1.1", "DAA1.2", "DAA1.3", "DAA2.1", "DAA2.2", "DAA2.3", "DAA4.1", "DAA4.2", "DAA4.3", "DAA5.1", "DAA7.1", "DAA7.2", "DAA8.1", "DFA3.1", "DFA4.1", "DFA4.2", "DAB1.1", "DAB1.2", "DAB1.3", "DAB2.1", "DAB2.2", "DAB2.3", "DAB4.1", "DAB4.2", "DAB4.3", "DAB5.1a", "DAB5.1b", "DAB7.1", "DAB7.2", "DAB8.1", "DFB2.1", "DFB3.1", "DFB4.1", "DFB4.2", "DFB6.1")
genes=c("DAA1.1", "DAA1.2", "DAA1.3", "DAA2.1", "DAA2.2", "DAA2.3", "DAA4.1", "DAA4.2", "DAA5.1", "DAA7.1" ,"DAA8.1", "DAA8.2", "DAA9.1", 
        "DBA3.1", "DBA4.1", "DBA4.2", "DBA6.1", 
        "DAB1.1", "DAB1.2", "DAB1.3", "DAB2.1", "DAB2.2", "DAB2.3", "DAB4.1", "DAB4.2", "DAB5.1a", "DAB5.1b", "DAB7.1" ,"DAB8.1", "DAB8.2", "DAB9.1", 
        "DBB2.1", "DBB3.1", "DBB4.1", "DBB4.2", "DBB6.1")

exons=c("E1", "E2", "E3", "E4")

pi_df <- data.frame()
# loop through each gene and exon
for(i in genes){
  for(j in exons){
    # read in the sequence file for current gene and exon
    in_file <- readDNAStringSet(paste0(i, "_", j, "_CDS.fa"))
    d <- as.dist(stringDist(in_file, ))
    align <- msa(in_file, "ClustalOmega", order = "input")
    CDS_pi <- round(nuc.div(as.DNAbin(align)), digits = 3)
    # assign pi value to appropriate cell in dataframe
    pi_df[i, j] <- CDS_pi
  }
}

# label rows and columns of dataframe using gene and exon names
rownames(pi_df) <- genes
colnames(pi_df) <- c("Exon 1", "Exon 2", "Exon 3", "Exon 4")


# count the number of sequences per gene and write that in a last column "count"
count <- data.frame()
for(i in genes){
  value <- max(seq(readDNAStringSet(paste0(i, "_", j, "_CDS.fa"))))
  count <- rbind(count, value)
}
rownames(count) <- genes
colnames(count) <- c("No.seq")

# combine tables containing pi values and counts
pi_df <- cbind(pi_df, count)

print(pi_df)

# convert rownames into a proper column "Genes"
pi_df <- cbind(Genes = rownames(pi_df), pi_df)
rownames(pi_df) <- NULL

# split the table into two for alpha aand beta genes
AlphaGenes=c("DAA1.1", "DAA1.2", "DAA1.3", "DAA2.1", "DAA2.2", "DAA2.3", "DAA4.1", "DAA4.2", "DAA5.1", "DAA7.1" ,"DAA8.1", "DAA8.2", "DAA9.1", 
             "DBA3.1", "DBA4.1", "DBA4.2", "DBA6.1")
alpha = pi_df[pi_df$Genes %in% AlphaGenes, ]
beta = pi_df[!(pi_df$Genes %in% AlphaGenes), ]

## save the output
write.csv(alpha, "alpha_pi.csv")
write.csv(beta, "beta_pi.csv")


########## plot pi values ##########
# alpha
alpha_tidy_df <- gather(alpha, Exons, pi, -Genes, -No.seq)

alpha_pi_plot <- ggplot(alpha_tidy_df, aes(x = Genes, y = pi, fill = Exons)) +
  geom_bar(stat = "identity", position="dodge") +
  geom_text(data = alpha_tidy_df[alpha_tidy_df$Exons == "Exon 2", ], aes(label=No.seq), vjust=-0.5, size=3) +
  theme_classic(base_size = 15) +
  theme(axis.text.x=element_text(angle=45, hjust=1 , face="italic")) +
  labs(x = "Alpha genes",
       y = expression(paste("Nucleotide diversity (", italic("\u03C0"), ")")),
       tag = "a") +
  ylim(0, 0.2) +
  scale_fill_manual(values = c("#c1272d", "#eecc16", "#008176", "#b3b3b3"))
# scale_fill_brewer(palette = "Set2")
#theme(strip.background = element_blank(),
#  strip.text.x = element_text(hjust = 0.5, face = "bold", size = 14),
# plot.tag = element_text())

# beta
beta_tidy_df <- gather(beta, Exons, pi, -Genes, -No.seq)

beta_pi_plot <- ggplot(beta_tidy_df, aes(x = Genes, y = pi, fill = Exons)) +
  geom_bar(stat = "identity", position="dodge") +
  geom_text(data = beta_tidy_df[beta_tidy_df$Exons == "Exon 2", ], aes(label=No.seq), vjust=-0.5, size=3) +
  theme_classic(base_size = 15) +
  theme(axis.text.x=element_text(angle=45, hjust=1 , face="italic")) +
  labs(x = "Beta genes",
       y = expression(paste("Nucleotide diversity (", italic("\u03C0"), ")")),
       tag = "b") +
  ylim(0, 0.2) +
  scale_fill_manual(values = c("#c1272d", "#eecc16", "#008176", "#b3b3b3"))

plot <- plot_grid(alpha_pi_plot, beta_pi_plot, nrow=2)
# Fig5 <- alpha_pi_plot / beta_pi_plot    ## not working! Forgot which package

ggsave("Fig5.png", plot, width = 7.08, height = 7, dpi = 300)
ggsave("Fig5.pdf", plot, width = 7.08, height = 7, dpi = 300)
