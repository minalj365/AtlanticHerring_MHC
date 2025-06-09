rm(list=ls())
library(gggenomes)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library(ggforce)
library(Cairo)

setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/gggenome/regions/Ragtag/")
file_prefix=c("Locus1", "Locus2", "Locus3", "Locus4", "Locus5", "Locus6", "Locus7", "Locus8", "Locus9")
#pdf(file="try.pdf", width = 8.5, height = 11)
#plots_list <- list()
#pdf("All_loci.pdf", width = 16, height = 18, dpi = 300, units = "cm")
pdf("All_loci.pdf", width = 20, height = 25, dpi = 300, units = "cm")
for (i in file_prefix) {
  
  GeneTrack <- read.csv(paste0(i,"_positions.csv")) %>%
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
    unite('colors', lineage:number, sep="")
  # strand 
  #group_by(seq_id) %>%
  #mutate(value = ifelse(strand == "+", min(start), max(end)),
  #start = ifelse(strand == "+", (start - value), abs(start - value)),
  # end = ifelse(strand == "+", (end - value), abs(end - value)))
  
  
  # preparing length track
  SeqTrack <- GeneTrack %>%
    group_by(seq_id) %>%
    summarise(length=max(end))
  
  ## plot
  gggenome <- gggenomes(seqs=SeqTrack, genes=GeneTrack)  
  plot <- gggenome +
    geom_seq(size=0.5) +         # draw contig/chromosome lines
    geom_bin_label(size = 2) +   # label each sequence 
    #geom_gene(aes(fill=gene), size=2) + 
    #geom_gene(aes(fill=gene), size=2) +        # draw genes as arrow
    #geom_gene_label(aes(label=gene, color = gene), fontface = "italic", size = 3, angle = 32, nudge_y = 0.15, nudge_x = -0.2) +
    
    geom_gene(size = 1.5) +        # draw genes as arrow
    geom_text_repel(aes(x=(x+xend)/2, y=y+.1, label=gene, color = colors), data=genes(),
                    size = 2, fontface = "bold", hjust=0, vjust=0, angle=25, direction = "x", min.segment.length = Inf) +
    theme(axis.text.x = element_text(size=12), legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) + 
    ggtitle(i) +
 scale_x_bp(suffix = "b", sep=" ")
  ggsave(paste0(i,".pdf"), path = "./plots_pdf", width = 20, height = 25, dpi = 300, units = "cm")
  ggsave(paste0(i,".png"), path = "./plots_png", width = 20, height = 25, dpi = 300, units = "cm")
  ### Following not working yet
  # Append the current plot to the list
  #grid.arrange(plot, ncol = 1)
  #plots_list[[i]] <- plot
  #grid.arrange(plot, ncol = 1)

  #print(file_prefix[[i]])
}
#dev.off()
