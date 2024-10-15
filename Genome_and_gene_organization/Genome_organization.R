rm(list=ls())
library(gggenomes)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library(ggforce)
library(Cairo)

setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/gggenome/assemblies/RagTag")
file_prefix=c("Reference")#c("CS2_hap1", "CS2_hap2", "CS4_hap1", "CS4_hap2", "CS5_hap1", "CS5_hap2", "CS7_hap1", "CS7_hap2", "CS8_hap1", "CS8_hap2", "CS10_hap1", "CS10_hap2", "BS1_hap1", "BS1_hap2", "BS2_hap1", "BS2_hap2", "BS3_hap1", "BS3_hap2", "BS4_hap1", "BS4_hap2", "BS5_hap1", "BS5_hap2", "BS6_hap1", "BS6_hap2", "NSSH2_hap1", "NSSH2_hap2", "NSSH10_hap1", "NSSH10_hap2", "Reference")
for (i in file_prefix) {
  pretended_order <- c("chr7", "chr13", "chr5", "chr8")
  
  GeneTrack <- read.csv(paste0(i,"_positions.csv")) %>%
    separate_wider_delim(cols=seq_id, delim="_", names=c("sample", "haplotype", "chr")) %>%
    unite(sample:haplotype, col="assembly", sep="_") %>%
    arrange(factor(chr, levels=pretended_order)) %>%
    # adding a colors column
    mutate(gene = str_sub(gene, end = -2)) %>%
    mutate(lineage = gene,
           number = gene) %>%
    mutate(lineage = str_sub(lineage, 2, 2)) %>%
    mutate(number = str_sub(number, start = 4)) %>%
    unite('colors', lineage:number, sep="") %>%
    
    select(seq_id=chr, start, end, strand, gene, colors) %>%
    mutate(seq_id = recode(seq_id, "chr7" = "Chromosome 7", "chr13" = "Chromosome 13", "chr5" = "Chromosome 5", "chr8" = "Chromosome 8"))
  
  
  SeqTrack <- read.csv(paste0(i,"_lengths.csv")) %>%
    separate_wider_delim(cols=seq_id, delim="_", names=c("sample", "haplotype", "chr")) %>%
    unite(sample:haplotype, col="assembly", sep="_") %>%
    arrange(factor(chr, levels=pretended_order)) %>%
    select(seq_id=chr, length, assembly) %>%
    mutate(seq_id = recode(seq_id, "chr7" = "Chromosome 7", "chr13" = "Chromosome 13", "chr5" = "Chromosome 5", "chr8" = "Chromosome 8"))

  ## plot
  gggenome <- gggenomes(seqs=SeqTrack, genes=GeneTrack)  
  gggenome +
    geom_seq(size=0.5) +         # draw contig/chromosome lines
    geom_bin_label(size = 4) +   # label each sequence 
    #geom_gene(aes(fill=gene), size=2) + 
    #geom_gene(aes(fill=gene), size=2) +        # draw genes as arrow
    #geom_gene_label(aes(label=gene, color = gene), fontface = "italic", size = 3, angle = 32, nudge_y = 0.15, nudge_x = -0.2) +
    
    geom_gene(size = 2) +        # draw genes as arrow
    geom_text_repel(aes(x=(x+xend)/2, y=y+.1, label=gene, color = colors), data=genes(),
                    size = 2.5, fontface = "bold", hjust=0, vjust=0, angle=45, direction = "x", 
                    min.segment.length = Inf, max.overlaps = Inf) +
    theme(axis.text.x = element_text(size=12), legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) + 
    ggtitle(i) +
    scale_x_continuous(breaks = c(0, 10000000, 20000000, 30000000, 40000000),
               labels = c("0 Mb", "10 Mb", "20 Mb", "30 Mb", "40 Mb"),
               limits = c(-7000000, 39000000))
  
  #scale_x_bp(suffix = "b", sep=" ", limits=c((0-7000000), 39000000))
  ggsave(paste0(i,".pdf"), path = "./plots_pdf", width = 25, height = 11, dpi = 300, units = "cm")
  ggsave(paste0(i,".png"), path = "./plots_png", width = 25, height = 11, dpi = 300, units = "cm")
}


####### reference assembly #########
pretended_order <- c("chr7", "chr13", "chr5")

GeneTrack <- read.csv("Reference_positions.csv") %>%
  separate_wider_delim(cols=seq_id, delim="_", names=c("sample", "chr")) %>%
  arrange(factor(chr, levels=pretended_order)) %>%
  # adding a colors column
  mutate(gene = str_sub(gene, end = -2)) %>%
  mutate(lineage = gene,
         number = gene) %>%
  mutate(lineage = str_sub(lineage, 2, 2)) %>%
  mutate(number = str_sub(number, start = 4)) %>%
  unite('colors', lineage:number, sep="") %>%
  
  select(seq_id=chr, start, end, strand, gene, colors) %>%
  mutate(seq_id = recode(seq_id, "chr7" = "Chromosome 7", "chr13" = "Chromosome 13", "chr5" = "Chromosome 5"))


SeqTrack <- read.csv("Reference_lengths.csv") %>%
  separate_wider_delim(cols=seq_id, delim="_", names=c("sample", "chr")) %>%
  arrange(factor(chr, levels=pretended_order)) %>%
  select(seq_id=chr, length) %>%
  mutate(seq_id = recode(seq_id, "chr7" = "Chromosome 7", "chr13" = "Chromosome 13", "chr5" = "Chromosome 5"))

## plot
gggenome <- gggenomes(seqs=SeqTrack, genes=GeneTrack)  
gggenome +
    geom_seq(size=0.5) +         # draw contig/chromosome lines
    geom_bin_label(size = 4) +   # label each sequence 
    #geom_gene(aes(fill=gene), size=2) + 
    #geom_gene(aes(fill=gene), size=2) +        # draw genes as arrow
    #geom_gene_label(aes(label=gene, color = gene), fontface = "italic", size = 3, angle = 32, nudge_y = 0.15, nudge_x = -0.2) +
    
    geom_gene(size = 2) +        # draw genes as arrow
    geom_text_repel(aes(x=(x+xend)/2, y=y+.1, label=gene, color = colors), data=genes(),
                    size = 2.5, fontface = "bold", hjust=0, vjust=0, angle=45, direction = "x", 
                    min.segment.length = Inf, max.overlaps = Inf) +
    theme(axis.text.x = element_text(size=12), legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")) + 
  ggtitle("Reference") +
    scale_x_continuous(breaks = c(0, 10000000, 20000000, 30000000, 40000000),
               labels = c("0 Mb", "10 Mb", "20 Mb", "30 Mb", "40 Mb"),
               limits = c(-7000000, 39000000))
ggsave("Reference.png", path = "./plots_png", width = 25, height = 11, dpi = 300, units = "cm")
