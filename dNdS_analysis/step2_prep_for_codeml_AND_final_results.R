rm(list=ls())
library(msa)
library(seqinr)
library(ape)
library(readr)
library(stringr)
library(tidyverse)
library(ggmsa)
library(cowplot)
stringsAsFactors=FALSE

setwd("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/dNdS")
# genes with at least 5 alleles/sequences
genes=c("DAA1.1", "DAA1.2", "DAA1.3", "DAA2.1", "DAA2.2", "DAA2.3", "DAA4.1", "DAA4.2", "DAA5.1", "DAA7.1", 
        "DAB1.1", "DAB1.2", "DAB1.3", "DAB2.1", "DAB2.2", "DAB2.3", "DAB4.1", "DAB4.2", "DAB5.1a", "DAB5.1b", "DAB7.1",
        "DBA3.1", "DBA4.1", 
        "DBB2.1", "DBB3.1", "DBB4.1", "DBB6.1")
exons=c("E2", "E3")

############## make alignments ready for pal2nal and trees ready for codeml ############

for(i in genes){
  for(j in exons){
    in_cds <- readDNAStringSet(paste0("Raw_sequences/", i, "_", j, "_CDS.fa"))
    in_aa <- readAAStringSet(paste0("Raw_sequences/", i, "_", j, "_AA.fa"))
    
    align_cds <- msa(in_cds, "ClustalW", order = "aligned")
    align_aa <- msa(in_aa, "ClustalW", order = "aligned")
    
    dist_cds <- dist.aa(as.DNAbin(align_cds)) # ape
    dist_aa <- dist.aa(as.AAbin(align_aa)) # ape
    
    align_cds_clustal<-DNAStringSet(align_cds) 
    Biostrings::writeXStringSet(align_cds_clustal, paste0("alignments/", i, "_", j, "_CDS.fa"))
    
    align_aa_clustal<-AAStringSet(align_aa) 
    Biostrings::writeXStringSet(align_aa_clustal, paste0("alignments/", i, "_", j, "_AA.fa"))
    
    tree_cds <- bionj(dist_cds)
    tree_aa <- bionj(dist_aa)
    
    plot.phylo(tree_cds, type = "unrooted", edge.width = 2, main = "CDS distance")
    write.tree(tree_cds, file = paste0("trees/", i, "_", j, "_CDS.tre"))
    
    plot.phylo(tree_aa, type = "unrooted", edge.width = 2, main = "AA distance")
    write.tree(tree_aa, file = paste0("trees/", i, "_", j, "_AA.tre"))
  }
}

## use these outputs to run codeml ######

########## processing codeml NS01278 output to make a table ##########
genes=c("DAA1.1_E2", "DAA1.2_E2", "DAA1.3_E2", "DAA2.1_E2", "DAA2.2_E2", "DAA2.3_E2", "DAA4.1_E2", "DAA4.2_E2", "DAA7.1_E2",
        "DAB1.1_E2", "DAB1.2_E2", "DAB1.3_E2", "DAB2.1_E2", "DAB2.2_E2", "DAB2.3_E2", "DAB4.1_E2", "DAB4.2_E2", "DAB5.1a_E2", "DAB5.1b_E2", "DAB7.1_E2",
        "DBA3.1_E2", "DBA4.1_E2",
        "DBB2.1_E2", "DBB3.1_E2", "DBB4.1_E2", "DBB6.1_E2",
        "DAA1.1_E3", "DAA1.2_E3", "DAA1.3_E3", "DAA2.1_E3", "DAA2.2_E3", "DAA2.3_E3", "DAA4.1_E3", "DAA4.2_E3", "DAA7.1_E3",
        "DAB1.1_E3", "DAB1.2_E3", "DAB1.3_E3", "DAB2.1_E3", "DAB2.2_E3", "DAB2.3_E3", "DAB4.1_E3", "DAB4.2_E3", "DAB5.1a_E3", "DAB5.1b_E3", "DAB7.1_E3",
        "DBA3.1_E3", "DBA4.1_E3",
        "DBB2.1_E3", "DBB3.1_E3", "DBB4.1_E3", "DBB6.1_E3")
# Create an empty list to store the data frames
dataframe_list <- list()

for(i in genes){
  OutFile <- readLines(paste0("codeml/DAA1.1_E2/DAA1.1_E2.out"))
  OutFile = readLines(paste0("codeml/", i, "/", i, ".out"))
  patternA <- "NSsites"
  pattern0 <- "omega"
  pattern1 <- "w:"
  pattern2 <- "^p:"
  pattern3 <- "p0 ="
  pattern4 <- "p ="
  pattern5 <- "p1 ="
  pattern6 <- "lnL"
  pattern7 <- "tree length for dN:"
  pattern8 <- "tree length for dS:"
  
  ### extract all lines that has output information of interest ###
  values <- grep(paste(patternA, pattern0, pattern1, pattern2, pattern3, pattern4, pattern5, pattern6, pattern7, pattern8, sep = "|"), OutFile, value = TRUE)
  
  ### extract ln values so that I can calculate 2*deltaln
  lnL_M1 <- str_extract(values[5], "(-)\\d+\\.\\d+")
  lnL_M2 <- str_extract(values[8], "(-)\\d+\\.\\d+")
  lnL_M7 <- str_extract(values[12], "(-)\\d+\\.\\d+")
  lnL_M8 <- str_extract(values[16], "(-)\\d+\\.\\d+")
  ### calculate 2*deltaln
  LRT_M1_vs_M2 <- 2*(as.numeric(lnL_M2)-as.numeric(lnL_M1))
  LRT_M7_vs_M8 <- 2*(as.numeric(lnL_M8)-as.numeric(lnL_M7))
  
  # if 2deltalnL is >13.8155, then P-value is <0.001
  
  ### extract dN/dS values
  omega <- as.numeric(str_extract(values[2], "\\d+\\.\\d+"))
  
  ### extract dN and dS values
  dN <- as.numeric(str_extract(values[3], "\\d+\\.\\d+"))
  dS <- as.numeric(str_extract(values[4], "\\d+\\.\\d+"))
  
  ### Extract parameter values  
  M1_par <- paste(values[6], values[7])
  M2_par <- paste(values[9], values[10])
  M7_par <- c(values[13])
  M8_par <- paste(values[17], values[18])
  
  ### process parameter values to show in the table
  # First, extract only numbers separacted by comma
  M1_numbers <- gsub("[^0-9.]+", " ", M1_par) 
  M2_numbers <- gsub("[^0-9.]+", " ", M2_par) 
  M7_numbers <- gsub("[^0-9.]+", " ", M7_par) 
  M8_numbers <- gsub("[^0-9.]+", " ", M8_par)
  M8_numbers <- gsub("\\b(0|1)\\b", "", M8_numbers)
  M8_numbers <- gsub("\\s+", " ", M8_numbers)

  # Split numbers into individual elements
  M1_numbers <- strsplit(trimws(M1_numbers), " ")[[1]]
  M2_numbers <- strsplit(trimws(M2_numbers), " ")[[1]]
  M7_numbers <- strsplit(trimws(M7_numbers), " ")[[1]]
  M8_numbers <- strsplit(trimws(M8_numbers), " ")[[1]]
  
  ### define labels for each parameter
  M1_labels <- c("p0", "p1", "w0", "w1")
  M2_labels <- c("p0", "p1", "p2", "w0", "w1", "w2")
  M7_labels <- c("p", "q")
  M8_labels <- c("p0", "p", "q", "p1", "w")
  
  # Apply character transformations
  M1_numbers_with_labels <- paste(M1_labels[1:length(M1_numbers)], M1_numbers, sep = "=")
  M2_numbers_with_labels <- paste(M2_labels[1:length(M2_numbers)], M2_numbers, sep = "=")
  M7_numbers_with_labels <- paste(M7_labels[1:length(M7_numbers)], M7_numbers, sep = "=")
  M8_numbers_with_labels <- paste(M8_labels[1:length(M8_numbers)], M8_numbers, sep = "=")
  
  ### insert comma and space between parameter values
  M1_numbers_with_labels <- paste(M1_numbers_with_labels, collapse = ", ")
  M2_numbers_with_labels <- paste(M2_numbers_with_labels, collapse = ", ")
  M7_numbers_with_labels <- paste(M7_numbers_with_labels, collapse = ", ")
  M8_numbers_with_labels <- paste(M8_numbers_with_labels, collapse = ", ")
  
  df <- data.frame(Gene = i,
                   omega = omega,
                   dN = dN,
                   dS = dS,
                   LRT_M1_vs_M2 = LRT_M1_vs_M2,
                   LRT_M7_vs_M8 = LRT_M7_vs_M8,
                   M1_parameters = M1_numbers_with_labels,
                   M2_parameters = M2_numbers_with_labels,
                   M7_parameters = M7_numbers_with_labels,
                   M8_parameters = M8_numbers_with_labels)
  # Add the data frame to the list
  dataframe_list[[i]] <- df
}

# Combine the data frames into one by column names
combined_df <- bind_rows(dataframe_list)
write.csv(combined_df, "results_codeml.csv")

#qchisq(0.999, df=2)

### plotting omega ###
plot_data <- separate(combined_df, col = Gene, into = c("Gene", "Exons"), sep = "_") %>%
  select(Gene, Exons, omega) %>%
  mutate(Exons = recode(Exons, "E2" = "Exon 2", "E3" = "Exon 3"))

    #pivot_wider(names_from = Exon, values_from = omega) %>%
    #rename(`Exon 2` = E2,
           #`Exon 3` = E3)


alpha_plot_data <- filter(plot_data, grepl("D[AB]A.*", Gene))
alpha_plot <-  ggplot(alpha_plot_data, aes(x = Gene, y = omega, fill = Exons)) +
  geom_bar(stat = "identity", position="dodge") +
  theme_classic(base_size = 15) +
  theme(axis.text.x=element_text(angle=45, hjust=1 , face="italic")) +
  labs(x = "Alpha genes",
       y = expression(italic("dN/dS")),
       tag = "a") +
  ylim(0,5.5) +
  geom_hline(yintercept = 1, linetype='dashed', col = 'red', linewidth = 2) +
  scale_fill_manual(values = c("#eecc16", "#008176"))




beta_plot_data <- filter(plot_data, grepl("D[AB]B.*", Gene))
beta_plot <-  ggplot(beta_plot_data, aes(x = Gene, y = omega, fill = Exons)) +
  geom_bar(stat = "identity", position="dodge") +
  theme_classic(base_size = 15) +
  theme(axis.text.x=element_text(angle=45, hjust=1 , face="italic")) +
  labs(x = "Beta genes",
       y = expression(italic("dN/dS")),
       tag = "b") +
  ylim(0,5.5) +
  geom_hline(yintercept = 1, linetype='dashed', col = 'red', linewidth = 2) +
  scale_fill_manual(values = c("#eecc16", "#008176"))

FigOmega <- alpha_plot / beta_plot
FigOmega <- plot_grid(alpha_plot, beta_plot, nrow=2)    # plot_grid is a function in cowplot


ggsave("FigOmega.png", FigOmega, width = 20, height = 23, dpi = 300, units = "cm")
ggsave("FigOmega.pdf", FigOmega, width = 20, height = 23, dpi = 300, units = "cm")





############## extracting selected amino acids to use in ggmsa highlight_positions
# These genes did not have any positively selected site in M2 and M8 so removing them from the below gene lost: DBB2.1_E2, DBB3.1_E2, DBB4.1_E2, DBB6.1_E2, DAB5.1a_E3, DAB5.1b_E3, DBA3.1_E3, DBA4.1_E3, DBB3.1_E3
genes=c("DAA1.1_E2", "DAA1.2_E2", "DAA1.3_E2", "DAA2.1_E2", "DAA2.2_E2", "DAA2.3_E2", "DAA4.1_E2", "DAA4.2_E2", "DAA7.1_E2",
        "DAB1.1_E2", "DAB1.2_E2", "DAB1.3_E2", "DAB2.1_E2", "DAB2.2_E2", "DAB2.3_E2", "DAB4.1_E2", "DAB4.2_E2", "DAB5.1a_E2", "DAB5.1b_E2", "DAB7.1_E2",
        "DBA3.1_E2", "DBA4.1_E2",
        "DAA1.1_E3", "DAA1.2_E3", "DAA1.3_E3", "DAA2.1_E3", "DAA2.2_E3", "DAA2.3_E3", "DAA4.1_E3", "DAA4.2_E3", "DAA7.1_E3",
        "DAB1.1_E3", "DAB1.2_E3", "DAB1.3_E3", "DAB2.1_E3", "DAB2.2_E3", "DAB2.3_E3", "DAB4.1_E3", "DAB4.2_E3", "DAB7.1_E3",
        "DBB2.1_E3", "DBB4.1_E3", "DBB6.1_E3")
site_list <- list()

while (!is.null(dev.list())) dev.off()
pdf("msa_plots.pdf")
for(i in genes){
  ## read and align amino acid files
  AA_file <- readAAStringSet(paste0("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/dNdS/Raw_sequences/",
  i, "_AA.fa"))
  align <- msa(AA_file, "ClustalW", order = "aligned")
  class(align) <- "AAMultipleAlignment"
  
  ## read and process codeml M8 NSSites output to extract positions to highlight in the alignment
  M8_codeml_out = readLines(paste0("C:/Users/minal03/OneDrive - Texas A&M University/U-Drive/MHC/dNdS/M8_codeml/",
                                   i, "/", i, ".out"))
  # Define the starting and ending patterns
  start_pattern <- "Bayes Empirical Bayes"
  end_pattern <- "The grid"
  # Use grep to find the line numbers matching the patterns
  start_line_number <- grep(start_pattern, M8_codeml_out)
  end_line_number <- grep(end_pattern, M8_codeml_out)
  # Extract the lines between the matched line numbers
  extracted_lines <- M8_codeml_out[(start_line_number + 1):(end_line_number - 1)]
  # Specify the number of lines to remove from the beginning
  lines_to_remove <- 5
  # Remove the specified number of lines from the beginning
  lines <- extracted_lines[(lines_to_remove + 1):length(extracted_lines)]
  # Prepare the content to write in a file
  heading <- paste0(i, "\n", "\t", "Site", "\t", "Pr(w>1)", "\t", "Posterior mean +- SE for w", "\n")
  lines_for_file <- paste0(lines, collapse = "\n")
  site_list[[i]] <- paste0(heading, lines_for_file)
  
  
  ######## converting character lines into a dataframe ########
  # Split the character elements by spaces
  split_elements <- strsplit(lines, "\\s+")
  # Convert the list of split elements into a data frame
  dataframe <- as.data.frame(do.call(rbind, split_elements), stringsAsFactors = FALSE)[,2:4]
  # Assign column names
  colnames(dataframe) <- c("position", "AA", "Prob")
  position_highlight <- as.vector(dataframe$position)
  
  ## plot alignments and highlight selected residues
  msa_plot <- ggmsa(align, char_width = 0.5, seq_name = TRUE, use_dot = TRUE, 
        consensus_views = TRUE, disagreement = TRUE, ref = rownames(align)[1],
        border = NA, font = "helvetical",
        position_highlight = position_highlight)
  msa_plot <- msa_plot +
    ggtitle(i) +
    theme(axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 6))
  print(msa_plot)
}
# save plots
dev.off()
# write sites into a file and save file
site_list_file <- unlist(site_list)
writeLines(site_list_file, "site_list.txt")
