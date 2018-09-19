## CRISPRanto

A (Snakemake) pipeline for CRISPR offtarget detection

Used in Raitskin et al. (2018?)


### Requirements

* Python3.5+ with latest Snakemake
* jellyfish 2 (e.g. jellyfish-2.2.6)
* bbmap >= 38.06
* bedtools >= 2.26.0
* if on Mac: GNU versions of sed, awk, join

### Data for Raitskin et al.

* TAIR10 gene annotation https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
* TAIR10 chromosomes https://www.arabidopsis.org/download_files/Sequences/TAIR10_blastsets/TAIR10_seq_20101214_updated
* Use these data to extract all sequences of CDS (i.e. coding part of exon) features
* Functional annotation
E.g. from https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_functional_descriptions
* Need 2-column tab-separated file: AGI,Function(with space replaced by "_")
