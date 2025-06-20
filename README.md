# Atlantic Herring MHC Class II Diversity

This repository contains scripts and resources associated with the study:

**Unheralded high MHC Class II polymorphism in the abundant Atlantic herring resolved by long-read sequencing**  now available on bioRxiv - [https://www.biorxiv.org/content/10.1101/2025.06.08.658498v1]

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

## 📂 Repository Structure

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
