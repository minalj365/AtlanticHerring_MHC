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
Below is an overview of the directory structure:

| Folder                          | Description                                                                                                                                           |
| ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Annotation/`                   | Scripts for manual annotation and curation of MHC class II genes in the Atlantic herring reference genome and PacBio haplotype assemblies.            |
| `Genome_and_gene_organization/` | Genome-wide organization of MHC class II genes across 29 haploid assemblies, including visualization of all 9 loci and their structural arrangements. |
| `Nucleotide_diversity/`         | Scripts to extract MHC gene sequences and calculate nucleotide diversity (π) across coding regions.                                                   |
| `de_novo_assemblies/`           | Scripts for generating de novo genome assemblies from PacBio HiFi reads using `hifiasm`.                                                              |
| `Protein_modeling/`             | Shannon entropy analysis and AlphaFold-based structural modeling of MHC class II proteins, highlighting peptide-binding diversity.                    |
| `Supertype_identification/`     | Identification and clustering of allelic supertypes at Locus 2 and Locus 4 based on amino acid similarity.                                            |
| `dNdS_analysis/`                | Scripts to calculate dN/dS ratios using PAML, visualize positively selected residues, and compare selection signatures with human MHC genes.          |
| `README.md`                     | This file. Provides a description of the repository and its organization.                                                                             |


README.md
This file — provides project description and directory guide.AtlanticHerring_MHC/
├── Annotation/               # Scripts for gene annotation and curation of MHC class II genes in reference and PacBio assemblies
├── Genome_and_gene_organization/  
│   └── Genome organization in 29 haploid genomes and gene organization at all 9 MHC loci per haplotype
├── Nucleotide_diversity/     # Scripts to extract sequences and calculate nucleotide diversity (π)
├── de_novo_assemblies/       # Scripts to build de novo genome assemblies from PacBio HiFi data
├── Protein_modeling/         # Scripts to calculate Shannon entropy and map it on predicted protein structures
├── Supertype_identification/ # Analysis to identify major allelic groups (supertypes) at Locus 2 and 4
├── dNdS_analysis/            # Scripts to perform dN/dS analysis (PAML) and highlight positively selected residues in alignments
└── README.md                 # Project documentation
