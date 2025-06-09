rm(list=ls())
library(msa)
library(Biostrings)
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
library(cowplot)
library(scales)
library(Cairo)
library(ggrepel)

#install.packages("tinytex") # latex - to save alignment in pdf file

setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/Figures_final/Fig4_pi_dNdS")

#genes=c("DAA1.1", "DAA1.2", "DAA1.3", "DAA2.1", "DAA2.2", "DAA2.3", "DAA4.1", "DAA4.2", "DAA4.3", "DAA5.1", "DAA7.1", "DAA7.2", "DAA8.1", "DFA3.1", "DFA4.1", "DFA4.2", "DAB1.1", "DAB1.2", "DAB1.3", "DAB2.1", "DAB2.2", "DAB2.3", "DAB4.1", "DAB4.2", "DAB4.3", "DAB5.1a", "DAB5.1b", "DAB7.1", "DAB7.2", "DAB8.1", "DFB2.1", "DFB3.1", "DFB4.1", "DFB4.2", "DFB6.1")
genes=c("DAA1.1", "DAA1.2", "DAA1.3", "DAA2.1", "DAA2.2", "DAA2.3", "DAA4.1", "DAA4.2", "DAA5.1", "DAA7.1" ,"DAA8.1", "DAA8.2", "DAA9.1", 
        "DBA3.1", "DBA4.1", "DBA6.1", 
        "DAB1.1", "DAB1.2", "DAB1.3", "DAB2.1", "DAB2.2", "DAB2.3", "DAB4.1", "DAB4.2", "DAB5.1a", "DAB5.1b", "DAB7.1" ,"DAB8.1", "DAB8.2", "DAB9.1", 
        "DBB2.1", "DBB3.1", "DBB4.1", "DBB4.2", "DBB6.1")

exons=c("E1", "E2", "E3", "E4")

pi_df <- data.frame()
# loop through each gene and exon
for(i in genes){
  for(j in exons){
    # read in the sequence file for current gene and exon
    in_file <- readDNAStringSet(paste0("./Raw_seq/", i, "_", j, "_CDS.fa"))
    d <- as.dist(pwalign::stringDist(in_file, ))
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
  value <- max(seq(readDNAStringSet(paste0("./Raw_seq/", i, "_", j, "_CDS.fa"))))
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

# split the table into two for alpha and beta genes
AlphaGenes=c("DAA1.1", "DAA1.2", "DAA1.3", "DAA2.1", "DAA2.2", "DAA2.3", "DAA4.1", "DAA4.2", "DAA5.1", "DAA7.1" ,"DAA8.1", "DAA8.2", "DAA9.1", 
             "DBA3.1", "DBA4.1", "DBA6.1")
alpha = pi_df[pi_df$Genes %in% AlphaGenes, ]
beta = pi_df[!(pi_df$Genes %in% AlphaGenes), ]

## save the output
write.csv(alpha, "alpha_pi.csv")
write.csv(beta, "beta_pi.csv")


########## plot pi values ##########
# alpha
alpha_tidy_df <- gather(alpha, Exons, pi, -Genes, -No.seq)

Fig4a <- ggplot(alpha_tidy_df, aes(x = Genes, y = pi, fill = Exons)) +
  geom_bar(stat = "identity", position="dodge") +
  geom_text(data = alpha_tidy_df[alpha_tidy_df$Exons == "Exon 2", ], aes(label=No.seq), vjust=-2.2, size=2.5) +
  theme_classic(base_size = 8) +
  theme(axis.text.x=element_text(angle=45, hjust=1 , face="italic")) +
  labs(x = "Alpha genes",
       y = expression(paste("Nucleotide diversity (", italic("\u03C0"), ")")),
       tag = "a") +
  ylim(0, 0.2) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#c1272d", "#eecc16", "#008176", "#b3b3b3"))
# scale_fill_brewer(palette = "Set2")
#theme(strip.background = element_blank(),
#  strip.text.x = element_text(hjust = 0.5, face = "bold", size = 14),
# plot.tag = element_text())

# beta
beta_tidy_df <- gather(beta, Exons, pi, -Genes, -No.seq)

Fig4b <- ggplot(beta_tidy_df, aes(x = Genes, y = pi, fill = Exons)) +
  geom_bar(stat = "identity", position="dodge") +
  geom_text(data = beta_tidy_df[beta_tidy_df$Exons == "Exon 2", ], aes(label=No.seq), vjust=-2.2, size=2.5) +
  theme_classic(base_size = 8) +
  theme(axis.text.x=element_text(angle=45, hjust=1 , face="italic")) +
  labs(x = "Beta genes",
       y = expression(paste("Nucleotide diversity (", italic("\u03C0"), ")")),
       tag = "b") +
  ylim(0, 0.2) +
  theme(legend.position = c(0.9, 0.7)) +
  scale_fill_manual(values = c("#c1272d", "#eecc16", "#008176", "#b3b3b3"))


######################## Figures 4c and 4d: dN/dS  ###################
data <- read_tsv("dNdS_MEGA_for_plotting.tsv", col_names = T) %>% select(Gene, Exons, omega)
Fig4c <- data  %>% 
  filter(grepl("DAA", Gene) | grepl("DBA", Gene)) %>%
  ggplot(aes(x = Gene, y = omega, fill = Exons)) +
  geom_bar(stat = "identity", position="dodge") +
  theme_classic(base_size = 8) +
  theme(axis.text.x=element_text(angle=45, hjust=1 , face="italic")) +
  labs(x = "Alpha genes",
       y = expression(italic("dN/dS")),
       tag = "c") +
  ylim(0,3.5) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1, linetype='dashed', col = 'red', linewidth = 1.5) +
  scale_fill_manual(values = c("#eecc16", "#008176"))

Fig4d <- data  %>% 
  filter(grepl("DAB", Gene) | grepl("DBB", Gene)) %>%
  ggplot(aes(x = Gene, y = omega, fill = Exons)) +
  geom_bar(stat = "identity", position="dodge") +
  theme_classic(base_size = 8) +
  theme(axis.text.x=element_text(angle=45, hjust=1 , face="italic")) +
  labs(x = "Beta genes",
       y = expression(italic("dN/dS")),
       tag = "d") +
  ylim(0,3.5) +
  theme(legend.position = c(0.9, 0.8)) + 
  geom_hline(yintercept = 1, linetype='dashed', col = 'red', linewidth = 1.5) +
  scale_fill_manual(values = c("#eecc16", "#008176"))


######################## Figures 4e and 4f: human-herring dN and dS  ###################
df <- read_delim("human_herring_dNdS.tsv", delim="\t", col_names = T)
alpha_df <- df %>% filter(Chain == "alpha")
beta_df <- df %>% filter(Chain == "beta")

plot_human_herring_dNdS <- function(data){
  p <- ggplot(data, aes(x = dS, y = dN)) +
    geom_point(aes(color = Species), size = 2.2, alpha = 1) +  # dN points, opaque
    geom_point(aes(color = Species), size = 2.2, alpha = 0.5, position = position_nudge(y = -0.01)) +  # dS points, transparent
    geom_text_repel(aes(label = Gene, color = Species), 
                    fontface = "italic", 
                    nudge_y = 0.04,  # Adjust nudge to position labels above points
                    size = 3, 
                    box.padding = unit(0.45, "lines"), 
                    point.padding = unit(0.5, "lines"), 
                    segment.color = 'grey50',  # Adds arrows with grey color
                    segment.size = 0.5,  # Adjust the thickness of the arrow
                    max.overlaps = 10,  # Allows more flexibility in repelling
                    direction = "y",  # Encourages vertical movement over horizontal
                    arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
    scale_shape_manual(name = "Substitution type", 
                       values = c("dN" = 16, "dS" = 17),
                       labels = c(dN = expression(italic("dN")), dS = expression(italic("dS")))) +
    scale_color_manual(values = c("human" = "#0066CC", "herring" = "#990000"), guide = FALSE) +
    theme_classic(base_size = 8) +
    guides(shape = guide_legend(title = "Substitution\ntype")) +
    scale_x_continuous(labels = label_number(accuracy = 0.01)) +  # Format x-axis with two decimals
    scale_y_continuous(labels = label_number(accuracy = 0.01))  # Format y-axis with two decimals
  
  # Return the plot object
  return(p)
}


# alpha
Fig4e <- ggplot(alpha_df, aes(x = dS, y = dN)) +
  geom_point(aes(color = Species), size = 1.8, alpha = 1) +
  geom_text_repel(aes(label = Gene, color = Species), 
            fontface = "italic", 
            max.overlaps = Inf,
            box.padding = 0.4,
            point.padding = 0.3,
            segment.color = "gray40",
            size = 3) +
  labs(x = expression(paste(italic("dS"))),
       y = expression(paste(italic("dN"))),
       tag = "e") +
  scale_shape_manual(name = "Substitution type", 
                     values = c("dN" = 16, "dS" = 17))+
  #                  labels = c(expression(italic(dN)), expression(italic(dS)))) +
  scale_color_manual(values = c("human" = "#0066CC", "herring" = "#990000"), guide = FALSE) +
  theme_classic(base_size = 8) +
  guides(shape = guide_legend(title = "Substitution\ntype")) +
  scale_x_continuous(limits = c(0, 0.11), 
                     expand = c(0, 0),
                     breaks = seq(0.02, 0.11, by = 0.02),  # skips 0
                     labels = label_number(accuracy = 0.01)) +  # Format x-axis with two decimals
  scale_y_continuous(limits = c(0, 0.165),
                     expand = c(0, 0),
                     labels = label_number(accuracy = 0.01))  # Format y-axis with two decimals

##### beta 
Fig4f <- ggplot(beta_df, aes(x = dS, y = dN)) +
  geom_point(aes(color = Species), size = 1.8, alpha = 1) +
  geom_text_repel(aes(label = Gene, color = Species),
                  fontface = "italic",
                  size = 3,
                  max.overlaps = 100,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "gray40",
                  segment.size = 0.3) +
  labs(x = expression(paste(italic("dS"))),
       y = expression(paste(italic("dN"))),
       tag = "f") +
  scale_shape_manual(name = "Substitution type", 
                     values = c("dN" = 16, "dS" = 17))+
  #                  labels = c(expression(italic(dN)), expression(italic(dS)))) +
  scale_color_manual(values = c("human" = "#0066CC", "herring" = "#990000"), guide = FALSE) +
  theme_classic(base_size = 8) +
  guides(shape = guide_legend(title = "Substitution\ntype")) +
  scale_x_continuous(limits = c(0, 0.190), 
                     expand = c(0, 0),
                     breaks = seq(0, 0.190, by = 0.025),
                     labels = label_number(accuracy = 0.01)) +  # Format x-axis with two decimals
  scale_y_continuous(labels = label_number(accuracy = 0.01))  # Format y-axis with two decimals
Fig4f


####### Final plot #########
Fig4 <- plot_grid(Fig4a, Fig4b, Fig4c, Fig4d, Fig4e, Fig4f, nrow=3)
Cairo::CairoPDF(file = "Fig4.pdf", width = 7.08, height = 7.28)
print(Fig4)
dev.off()

ggplot2::ggsave("Fig4.pdf", Fig4, width = 7.08, height = 7.28, device = "pdf")
ggplot2::ggsave("Fig4.png", Fig4, width = 7.08, height = 7.28, device = "png")
