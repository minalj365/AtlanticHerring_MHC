rm(list=ls())
library(gggenomes)
library(tidyverse)
library(msa)
library(seqinr)
library(ape)
library(pegas)
library(phytools)
library(ggtree)
library(purrr)
library(stringr)
library(readr)
library(zoo)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library(ggforce)
library(Cairo)


########## Fig. 1 ##############



######### Fig 1A ###########
setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/gggenome/entire_distribution")
pretended_order <- c("chr7", "chr13", "chr5")

GeneTrack <- read.csv("positions.csv") %>%
  arrange(factor(seq_id, levels=pretended_order)) %>%
  mutate(seq_id = recode(seq_id, "chr7" = "Chromosome 7", "chr13" = "Chromosome 13", "chr5" = "Chromosome 5", "chr8" = "Chromosome 8"))


SeqTrack <- read.csv("lengths.csv") %>%
  arrange(factor(seq_id, levels=pretended_order)) %>%
  mutate(seq_id = recode(seq_id, "chr7" = "Chromosome 7", "chr13" = "Chromosome 13", "chr5" = "Chromosome 5", "chr8" = "Chromosome 8"))

## plot
colorPalette <- c("#88419D", "#78C679")

Fig1A <- gggenomes(seqs=SeqTrack, genes=GeneTrack)  
Fig1A <- Fig1A +
  geom_seq(size=0.7) +         # draw contig/chromosome lines
  geom_bin_label(size = 4.5) +   # label each sequence 
  geom_gene_note(aes(label=size_range), nudge_y = -0.1, size = 4) +
  #geom_gene(aes(fill=locus), size=2) +        # draw genes as arrow
  #geom_gene_label(aes(label=gene, color = gene), fontface = "italic", size = 3, angle = 32, nudge_y = 0.15, nudge_x = -0.2) +
  
  geom_gene(size = 3) +        # draw genes as arrow
  geom_text_repel(aes(x=(x+xend)/2, y=y+.1, label=locus, color = color), data=genes(),
                  size = 4, fontface = "bold", hjust=0, vjust=0, angle=45, direction = "x", 
                  min.segment.length = Inf, max.overlaps = Inf) +
  scale_fill_manual(values= colorPalette) +
  scale_color_manual(values= colorPalette) +
  theme(axis.text.x = element_text(size=12), legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) + 
  labs(tag = "a") +
  scale_x_bp(suffix = "b", sep=" ", limits=c((0-9000000), 34000000))





#ggsave("distribution.png", width = 25, height = 11, dpi = 300, units = "cm")
#ggsave("distribution.pdf", device = cairo_pdf, width = 25, height = 11, dpi = 300, units = "cm")




######### Fig 1B ###########
dir = "C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/Sequence_analysis/Fig1_tree/with_outgroup/"
in_file <- readDNAStringSet(paste0(dir, "Ref_alpha.fa"))
align <- msa(in_file, "ClustalO", order = "aligned")

labels <- data.frame(seq = rownames(align),
                     category = c(rep("DB", 2), "Sasa", rep("DA", 5)))

dist <- dist.dna(as.DNAbin(align), model = "JC69", pairwise.deletion = T) # ape
tree <- bionj(dist)
write.tree(tree, paste0(dir, "DA.nwk")) # open in FigTree to set an outgroup
DA_tree <- read.tree(paste0(dir, "DA_figtree.nwk"))


Fig1B <- ggtree(DA_tree, size=0.7) %<+% labels + geom_tiplab(hjust = -0.05, fontface = "italic", size = 5, aes(color = category)) +
  scale_color_manual(values=c(DA = "#c1272d", DB= "#008176", Sasa="black")) +
  theme(legend.position='none') +
  geom_treescale() + geom_rootedge() + geom_rootpoint(size=1) +
  labs(tag = "b") +
  ggplot2::xlim(0, 0.5)
#geom_text(aes(label=label), hjust=-.1) +
#geom_tiplab(hjust = -0.05) +
#geom_treescale(offset = -0.5) 

library(ggtree)
library(ggplot2)

# Example tree loading and labels preparation
# DA_tree <- read.tree("your_tree_file.newick")
# labels <- data.frame(label=rownames(DA_tree$tip.label), category=sample(c("DA", "DB", "Sasa"), size=length(DA_tree$tip.label), replace=TRUE))



######### Fig 1C ###########
dir = "C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/Sequence_analysis/Fig1_tree/with_outgroup/"

in_file <- readDNAStringSet(paste0(dir, "Ref_beta.fa"))
align <- msa(in_file, "ClustalO", order = "aligned")
dist <- dist.dna(as.DNAbin(align), model = "JC69", pairwise.deletion = T) # ape
tree <- bionj(dist)
write.tree(tree, paste0(dir, "DB.nwk")) # open in FigTree to set an outgroup
DB_tree <- read.tree(paste0(dir, "DB_figtree.nwk"))

labels <- data.frame(seq = rownames(align),category = c(rep("DB", 4), "Sasa", rep("DA", 6)))

Fig1C <- ggtree(DB_tree, size=0.7) %<+% labels + geom_tiplab(hjust = -0.05, fontface = "italic", size = 5, aes(color = category)) +
  scale_color_manual(values=c(DA = "#c1272d",DB= "#008176", Sasa="black")) +
  theme(legend.position='none') +
  geom_treescale() + geom_rootedge() + geom_rootpoint(size=1) +
  labs(tag = "c") +
  ggplot2::xlim(0, 0.5)


#Fig1 <- (Fig1A | Fig1B) / Fig1C
Fig1 <- Fig1A / (Fig1B | Fig1C)



ggsave(paste0("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/Figures/Fig1_250318.png"), Fig1, width = 20, height = 23, dpi = 300, units = "cm")
ggsave(paste0("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/Figures/Fig1_250318.pdf"), Fig1, width = 20, height = 23, dpi = 300, units = "cm")
