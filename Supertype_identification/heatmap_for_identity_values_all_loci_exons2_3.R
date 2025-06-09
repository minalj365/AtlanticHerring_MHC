rm(list=ls())
library(msa)
library(Biostrings)
library(ape)
library(bio3d)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(magick)

## colors for loci are selected from the following palettes
display.brewer.all(colorblindFriendly = TRUE)
col_L1 <- brewer.pal(n = 9,name = 'YlGn')
col_L2 <- brewer.pal(n = 9,name = 'YlOrRd')
col_L4 <- brewer.pal(n = 9,name = 'BuPu')
col_L5 <- brewer.pal(n = 9,name = 'Greys')
col_L7 <- brewer.pal(n = 9,name = 'GnBu')
col_L8 <- brewer.pal(n = 9,name = 'PuRd')
col_L9 <- brewer.pal(n = 9,name = 'OrRd')


## selected colors
colors_all <- list("Loci" = c('L1' = "#238443",
                              'L2' = "#FD8D3C",
                              'L4' = "#88419D",
                              'L5' = "#969696",
                              'L7' = "#4EB3D3",
                              'L8' = "#DF65B0",
                              'L9' = "#FDBB84"))

################## classical ###############
# Fig 4A: alpha exon 2
# Fig 4B: beta exon 2
# Fig 4C: alpha exon 3
# Fig 4D: beta exon 3

dir <- "C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/Sequence_analysis/clustering/classical/"
setwd(dir)

###################### amino acid entire #################
files <- c("alpha_all_entire_AA", "beta_all_entire_AA", "alpha_all_entire_AA", "beta_all_entire_AA")

align_ <- list() # create empty lists

pdf(file = "plots/Fig4_AA_entire_Identity.pdf") # create an empty pdf file to save plots
# fill up the lists
for (i in files){
  in_file <- readAAStringSet(paste0(dir, i, ".fa"))
  
  # creating alignments
  align_[[i]] <- msa(in_file, "ClustalW", order = "aligned") # msa package
  Biostrings::writeXStringSet(AAStringSet(align_[[i]]), paste0(dir, "alignments/", i, ".fa")) # save
  
  # calculating sequence identity to plot in a heatmap
  identity <- seqidentity(read.fasta(paste0(dir, "alignments/", i, ".fa"))) ## bio3d package
  
  # preparing metadata #####
  metadata <- data.frame(Genes = rownames(align_[[i]]))
  metadata <- metadata %>% mutate(Annotation = case_when(
    str_detect(Genes, ".*DA[AB]1.[1-5]") ~ "L1",
    str_detect(Genes, ".*DA[AB]2.[1-5]") ~ "L2",
    str_detect(Genes, ".*DA[AB]4.[1-5]") ~ "L4",
    str_detect(Genes, ".*DA[AB]5.[1-5]") ~ "L5",
    str_detect(Genes, ".*DA[AB]7.[1-5]") ~ "L7",
    str_detect(Genes, ".*DA[AB]8.[1-5]") ~ "L8",
    str_detect(Genes, ".*DA[AB]9.[1-5]") ~ "L9"))
  
  
  
  ### making annotation file for heatmap
  ann <- data.frame(metadata$Annotation)
  colnames(ann) <- "Loci"
  
  colAnn <- HeatmapAnnotation(df = ann,
                              which = 'col',
                              col = colors_all,
                              annotation_width = unit(c(1, 4), 'cm'),
                              gap = unit(1, 'mm'),
                              show_legend = FALSE)
  #annotation_legend_param = list(ann = list(direction = "horizontal")))
  
  RowAnn <- rowAnnotation(df = ann,
                          col = colors_all,
                          annotation_width = unit(c(1, 4), 'cm'),
                          gap = unit(1, 'mm'))
  # annotation_legend_param = list(ann = list(direction = "horizontal")))
  
  # plotting
  plot <- Heatmap(as.matrix(identity), top_annotation = colAnn, left_annotation = RowAnn,
                  column_title = i,
                  row_names_side = "left",
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 10),
                  row_names_gp = gpar(fontsize = 8),
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  use_raster = TRUE,
                  #row_dend_reorder = TRUE,
                  #column_dend_reorder = TRUE,
                  name = "Identity")
  #heatmap_legend_param = list(direction = "horizontal"))
  
  draw(plot) #, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
}
dev.off() # save plots

##############################exon 2 and exon 3 ############################

files <- c("alpha_all_E2", "beta_all_E2", "alpha_all_E3", "beta_all_E3")

align_ <- list() # create empty lists

pdf(file = "plots/Fig4_all_Identity.pdf") # create an empty pdf file to save plots
# fill up the lists
for (i in files){
  in_file <- readDNAStringSet(paste0(dir, i, ".fa"))
  
  # creating alignments
  align_[[i]] <- msa(in_file, "ClustalW", order = "aligned") # msa package
  Biostrings::writeXStringSet(DNAStringSet(align_[[i]]), paste0(dir, "alignments/", i, ".fa")) # save
  
  # calculating sequence identity to plot in a heatmap
  identity <- seqidentity(read.fasta(paste0(dir, "alignments/", i, ".fa"))) ## bio3d package
  
  # preparing metadata #####
  metadata <- data.frame(Genes = rownames(align_[[i]]))
  metadata <- metadata %>% mutate(Annotation = case_when(
    str_detect(Genes, ".*DA[AB]1.[1-5]") ~ "L1",
    str_detect(Genes, ".*DA[AB]2.[1-5]") ~ "L2",
    str_detect(Genes, ".*DA[AB]4.[1-5]") ~ "L4",
    str_detect(Genes, ".*DA[AB]5.[1-5]") ~ "L5",
    str_detect(Genes, ".*DA[AB]7.[1-5]") ~ "L7",
    str_detect(Genes, ".*DA[AB]8.[1-5]") ~ "L8",
    str_detect(Genes, ".*DA[AB]9.[1-5]") ~ "L9"))
  
  
  
  ### making annotation file for heatmap
  ann <- data.frame(metadata$Annotation)
  colnames(ann) <- "Loci"
  
  colAnn <- HeatmapAnnotation(df = ann,
                              which = 'col',
                              col = colors_all,
                              annotation_width = unit(c(1, 4), 'cm'),
                              gap = unit(1, 'mm'),
                              show_legend = FALSE)
                              #annotation_legend_param = list(ann = list(direction = "horizontal")))
                              
  RowAnn <- rowAnnotation(df = ann,
                          col = colors_all,
                          annotation_width = unit(c(1, 4), 'cm'),
                          gap = unit(1, 'mm'))
                         # annotation_legend_param = list(ann = list(direction = "horizontal")))
  
  # plotting
  plot <- Heatmap(as.matrix(identity), top_annotation = colAnn, left_annotation = RowAnn,
                  column_title = i,
                  row_names_side = "left",
                  column_names_side = "top",
                  column_names_gp = gpar(fontsize = 10),
                  row_names_gp = gpar(fontsize = 8),
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  use_raster = TRUE,
                  #row_dend_reorder = TRUE,
                  #column_dend_reorder = TRUE,
                  name = "Identity")
                  #heatmap_legend_param = list(direction = "horizontal"))
  
  draw(plot) #, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
}
dev.off() # save plots
