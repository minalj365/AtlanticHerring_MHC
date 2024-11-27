rm(list=ls())
library(gggenomes)
library(ggrepel)
library(tidyverse)
library(gridExtra)
library(ggforce)
library(Cairo)
library(cowplot)

setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/2024/Figures_and_Tables/Fig2/")

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
Fig2A <- gggenome +
  geom_seq(size=0.6) +         # draw contig/chromosome lines
  geom_bin_label(size = 3) +   # label each sequence 
  #geom_gene(aes(fill=gene), size=2) + 
  #geom_gene(aes(fill=gene), size=2) +        # draw genes as arrow
  #geom_gene_label(aes(label=gene, color = gene), fontface = "italic", size = 3, angle = 32, nudge_y = 0.15, nudge_x = -0.2) +
  
  geom_gene(size = 2.2) +        # draw genes as arrow
  geom_text_repel(aes(x=(x+xend)/2, y=y+.1, label=gene, color = colors), data=genes(),
                  size = 2.8, fontface = "bold", hjust=0, vjust=0, angle=30, direction = "x", min.segment.length = Inf) +
  theme(axis.text.x = element_text(size=10), legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")) + 
  scale_x_bp(suffix = "b", sep=" ")
  #labs(tag = "a")

ggsave(paste0("Fig2A.png"), Fig2A, width = 4, height = 4, dpi = 300)
ggsave(paste0("Fig2A.pdf"), Fig2A, width = 4, height = 4, dpi = 300)







############### Fig 2B ####################
# dotplot of haplotypes (NSSH10)

setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/Figures/dotplot_haplotypes/")

my_thm= list( theme_bw(base_size = 12),
              theme(panel.grid.minor = element_blank()),
              theme(panel.grid.major.y = element_line(color = "grey")),
              theme(axis.line = element_line(colour = "black")),
              theme(aspect.ratio =1,legend.position = "none"),
              theme(axis.text=element_text(size=12)))


mydot <- read.table("dot.txt", header = FALSE) #read the coordinates for the dots
myline <- read.table("line.txt", header = FALSE) #read the coordinates for the lines
myxtics <- read.table("xticks.txt", header = TRUE) #read the x axis ticks
myytics <- read.table("yticks.txt", header = TRUE) #read the y-axis ticks
mycolor <- c("red","blue") #assign the colors to your own color scheme
names(mycolor) <- c("F","R") #name the colors with the forward and reverse codes. F and R should match the first and second color, respectively

Fig2B <- ggplot(mydot, aes(x = V1, y=V2, color=V3)) + geom_point() #create the ggplot object with the dots
#now create the plot. Parameters can be modified
Fig2B <- Fig2B + geom_segment(data = myline, aes(x=V1, y=V2, xend=V3, yend=V4, color = V5)) + 
  scale_y_continuous(labels = function(x)round(x/1000, 2), name = "Haplotype 1 (kb)") +
  scale_x_continuous(labels = function(x)round(x/1000, 2), name = "Haplotype 2 (kb)") +
  scale_color_manual(values = mycolor) +
  theme_bw(base_size = 12) +
  #theme(strip.background = element_blank(),
  # strip.text.x = element_text(face = "bold", size = 14))
  theme(strip.background = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        aspect.ratio =1,legend.position = "none") 
  labs(tag = "b")


ggsave(paste0("Fig2B.png"), Fig2B, width = 2.8, height = 3, dpi = 300)
ggsave(paste0("Fig2B.pdf"), Fig2B, width = 2.8, height = 3, dpi = 300)

plot_grid(Fig2A, Fig2B, labels = c('a', 'b'), label_size = 6, nrow=1, rel_widths = c(2, 1))

all_PLOTS<-plot_grid(p2, p3, p1, align = "v", nrow = 3)
