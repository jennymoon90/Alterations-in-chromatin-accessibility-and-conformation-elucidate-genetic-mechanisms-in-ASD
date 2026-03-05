# Alterations in chromatin accessibility and conformation elucidate genetic mechanisms in ASD

Code, data, and plots accompanying the manuscript

# Pipeline Stages:
1.	Genome-wide chromatin accessibility profiling (ATAC-seq) and differential analysis (Figure 1)
2.	Transcription factor (TF) binding and regulatory network inference (Figure 2 & 3)
3.	Common variant and genetic risk modeling (Figure 4)
4.	Rare variant impact analysis on gene regulation (Figure 4)
5.	3D genome (Hi-C) differential interaction and domain analysis (Figure 5)
6.	Multi-modal integration of 3D genome architecture and DNA methylation (Figure 5)
7.	Quantifying Individual and joint contributions of chromatin accessibility, 3D genome organization, and DNA methylation to gene expression changes (Figure 6)

# Script Number:
Refers to the script within the pipeline stage. 

# File Locations:
**data/** contains the ATAC, Hi-C, RNA-seq, and H3K27ac datasets, alongside all relevant metadata for the analysis.

**scripts/** contains the complete suite of computational scripts used throughout the study.

**results/** and plots/ contain representative results and figures generated from the pipeline.

**utils/** contains a suite of genomic utilities supporting reproducible multi-omic integration. Key assets include the hg19 reference genome, the ENCODE-standard blacklist regions, and the JASPAR 2022 motif library, specifically filtered for brain-expressed transcription factors.

*The comprehensive set of results and data tables will be uploaded in their entirety following the formal publication of the manuscript.
