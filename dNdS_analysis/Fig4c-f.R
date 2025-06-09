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

setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/Figures_final/Fig4_pi_dNdS")

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
