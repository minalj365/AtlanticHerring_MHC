##################### April 24, 2024 #################
########### Fig. 2 #############
rm(list=ls())
library(gggenomes)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library(ggforce)
library(Cairo)
library(cowplot)

setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/Figures_final/Fig2_Locus1_gggenome")

desired_order <- c(
  "BS1_hap1_chr7", "BS1_hap2_chr7",
  "BS4_hap1_chr7", "BS4_hap2_chr7",
  "CS8_hap1_chr7", "CS8_hap2_chr7",  
  "CS10_hap1_chr7", "CS10_hap2_chr7",
  "NSSH2_hap1_chr7", "NSSH2_hap2_chr7"
)

GeneTrack <- read.csv("Locus1_positions.csv") %>%
  # chage start and end to start with 0
  group_by(seq_id) %>%
  mutate(value = min(start),
         start = (start - value),
         end = (end - value)) %>%
  # adding a colors column
  mutate(gene = str_remove_all(gene, "G")) %>%
  mutate(lineage = gene,
         number = gene) %>%
  mutate(lineage = str_sub(lineage, 2, 2)) %>%
  mutate(number = str_sub(number, start = 4)) %>%
  mutate(seq_id = factor(seq_id, levels = desired_order)) %>%
  unite('colors', lineage:number, sep="")
# strand 
#group_by(seq_id) %>%
#mutate(value = ifelse(strand == "+", min(start), max(end)),
#start = ifelse(strand == "+", (start - value), abs(start - value)),
# end = ifelse(strand == "+", (end - value), abs(end - value)))


# preparing length track
SeqTrack <- GeneTrack %>%
  group_by(seq_id) %>%
  summarise(length=max(end)) %>%
  mutate(seq_id = factor(seq_id, levels = desired_order))

## plot
gggenome <- gggenomes(seqs=SeqTrack, genes=GeneTrack)  
colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")

Fig2 <- gggenome +
  geom_seq(size=0.4) +         # draw contig/chromosome lines
  geom_bin_label(size = 2.5) +   # label each sequence 
  #geom_gene_label(aes(label=gene, color = gene), fontface = "italic", size = 3, angle = 32, nudge_y = 0.15, nudge_x = -0.2) +
  
  geom_gene(size = 1.8) +        # draw genes as arrow
  geom_text_repel(aes(x=(x+xend)/2, y=y+.1, label=gene, color = colors), data=genes(),
                  size = 2.3, fontface = "bold.italic", hjust=0, vjust=0, angle=30, direction = "x", min.segment.length = Inf, box.padding = 0.1) +
  theme(axis.text.x = element_text(size=8), legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) + 
  scale_color_manual(values = colors) +
  scale_x_continuous(
    name = "Position (kb)",
    breaks = seq(0, 150000, by = 50000),
    limits = c(-25000, 155000),
    labels = function(x) paste0(x/1000, " kb"))
Fig2

ggsave(paste0("Fig2.png"), Fig2, width = 7.5, height = 4.5, dpi = 300)
ggsave(paste0("Fig2.pdf"), Fig2, width = 7.5, height = 4.5, dpi = 300)

