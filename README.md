# Atlantic Herring MHC Class II Diversity

This repository contains scripts and resources associated with the study:

**Unheralded high MHC Class II polymorphism in the abundant Atlantic herring resolved by long-read sequencing**  now available on bioRxiv with doi - 10.1101/2025.06.08.658498v1

## 📜 Overview

This study characterizes the extraordinary diversity of MHC Class II genes in Atlantic herring (*Clupea harengus*) using long-read PacBio HiFi sequencing data from 14 individuals across three geographic regions (Celtic Sea, Baltic Sea, and Norwegian Sea) along with the reference genome assembly (Ch_v2.0.2v2). By analyzing phased haplotype assemblies, we uncover:

- **Nine MHC Class II loci** distributed across four chromosomes
- Two distinct gene lineages: **DA (polymorphic)** and **DB (non-polymorphic)** and a total of 195 DAA, 220 DAB, 53 DBA, and 107 DBB sequences
- Extensive **copy number variation** and **divergent allelic supertypes**
- Strong **positive selection** in exon 2, encoding the peptide-binding region, as well as **high dN and dS than human counterparts**
- non-random association of alpha-beta genes pairs, indicating **co-evolution of alpha-beta gene pairs**
- **haplotype diversity** maintained by divergent supertypes

This work sets a new benchmark for understanding vertebrate MHC gene diversity.

---

## 🔬 Project Goals

- Identify and annotate MHC Class II genes across haplotype assemblies
- Measure nucleotide and amino acid diversity (π, dN/dS)
- Model MHC Class II protein structure and peptide-binding residues
- Group alleles into supertypes and assess their co-segregation
- Provide reproducible scripts and workflows for the above analyses

---

## 📂 Repository Structure
AtlanticHerring_MHC/
├── Annotation/ # Scripts for gene annotation and curation of MHC class II genes in reference and PacBio assemblies
├── Genome_and_gene_organization/ # MHC II organization in 29 haploid genomes (genome organization) and all 9 loci for each haploid genome (gene organization)
├── Nucleotide_diversity/ # Scripts for extracting sequences and calculate nucleotide diversity parameter π
├── de_novo_assemblies/ # Scripts to build de novo genome assemblies from PacBio data
├── Protein_modeling/ # Script to calculate shannon entropy values to plot on the predicted protein structure
├── Supertype_identification/ # Initial analysis to find major allelic groups in Locus 2 and 4
├── dNdS_analysis/ # Scripts for extracting sequences and perform dN/dS analysis in PAML, and plot the positively selected residues on the amino acid alignment for each gene.
└── README.md
