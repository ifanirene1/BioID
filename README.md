# BioID network propagation identified TWIST1-CHD8 module in cranial neural crest specification and ectomesenchymal potential

## Table of Contents

- [Project title](#BioID-network-propagation-identified-TWIST1-CHD8-module)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Results](#results)
    - [Proteomic data analysis and visualization](#proteomic-data-analysis-and-visualization)
    - [ChIP-seq data analysis](#chIP-seq-data-analysis)
    - [Cell migration assays](#cell-migration-assays)

## Introduction

Using BioID-proximity-labelling approach and network propagation analysis, we identified a chromatin regulatory module in the functional network of TWIST1, a key molecular nexus in cranial neural crest cell (NCC) development. Genomic occupancy and transcriptome data of chromatin regulatory module members CHD7, CHD8 and WHSC1 suggested that they share roles in activating NCC-related neurogenic and cell migratory programs. Compound heterozygous loss of function mutants further revealed their dynamic expression and stage-specific functional cooperativity for NCC specification: In early ectoderm differentiation, the TWIST1-chromatin modifier module potentiates NCC fate in bipotent neural progenitors while repels SOX2-driven cortical neural stem cell fate. Later, TWIST1-CHD8 interaction maintains the ectomesenchymal potential in NC stem cells. Combined mutations of Twist1 and Chd8 in mouse embryos also led to malformation of craniofacial tissues. Altogether, our results have suggested a role of TWIST1 in pioneering the assembly of an epigenetic complex responsible in establishing the promoter/enhancer repertoires for NCC fate progression and ectomesenchymal potential. Our findings also suggests a common mechanism underlying the neural crest-related pathologies associated with TWIST1, CHD7, CHD8 and WHSC1.

## Methods

### Proteomic data analysis and visualization
* Experimental Overview

The comprehensive protein interactome of TWIST1 was characterised using the BioID technique (Roux et al., 2012) in three mouse cell types: the cranial neural crest stem cell line (O9-1) (Ishii et al., 2012), the C3H10T1/2 multipotent mesenchymal cell line (Katagiri et al., 1990; Shea et al., 2003) and the 3T3 fibroblast cell line. 

![BioID-workflow](Figures/BioID-workflow.png)

* Computational analysis Overview

    * [BioID analysis by EdgeR](https://github.com/ifanirene1/BioID/blob/master/Merged_bioID_Reanalysis.R) 

### ChIP-seq data analysis

* [ChIP-seq analysis pipeline](https://github.com/ifanirene1/BioID/blob/master/Grimes_ChIP-seq_IF.sh) 


### Fluidigm high-througput qPCR analysis

* [high-througput qPCR analysis code](https://github.com/ifanirene1/BioID/blob/master/Biomark_TWIST_NPC_IF.R) 

### Cell migration assays

* [Time-lapse imaging analysis in ImageJ: NPC](https://github.com/ifanirene1/BioID/blob/master/NPC_migration_ImageJ_analysis.ijm) 
* [Time-lapse imaging analysis in ImageJ: MSC Scratch Assay](https://github.com/ifanirene1/BioID/blob/master/ScratchAssay_sensitive_ImageJ_analysis.ijm) 

