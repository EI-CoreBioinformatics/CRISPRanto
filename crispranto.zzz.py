import sys
import glob
from os.path import join, basename, dirname, exists
import csv

# CRISPRanto - pipeline for CRISPR offtarget mining
# authors: Christian Schudoma, Oleg Raitskin
# contact: christian.schudoma@earlham.ac.uk
# (C) 2017-2018 Earlham Institute
# License: MIT

def readCasParams(sheet):
    keys = ["pam", "k", "pamsite"]
    return {row[0]: dict(zip(keys, row[1:])) for row in csv.reader(open(sheet), delimiter=",")}

cas_params = readCasParams("cassheet.csv")

GENE_GFF = "TAIR10_GFF3_cds_only.gff" 
GENE_DESC = "AGILIST_all_gene_descriptions.tsv.2"

def get_kmer_file_for_cas(wc):
    return "TAIR10_cds_k" + cas_params[wc.cas]["k"] + ".fasta"

TARGETS = list()
for cas in cas_params:
    TARGETS.append(cas + ".targets_offtargets.genes.tsv")
    TARGETS.append(cas + ".targets_offtargets.genes.dedup.tsv")

rule all: input: TARGETS

# rule count_kmers:
#    input:
# i don't know how to transform this into a snakemake rule currently
# need to start with counting kmers externally
#
#
# source jellyfish-2.2.6
# TAIR10CDS=TAIR10_cds_20101214_updated.txt

# for k in 23 24 27; do
#  jellyfish count -m ${k} -t 8 -o TAIR10_cds_k${k}.jf -s 2G ${TAIR10CDS} && jellyfish dump TAIR10_cds_k${k}.jf > TAIR10_cds_k${k}.fasta";
# done

rule extract_uniqmers:
    input:
        get_kmer_file_for_cas
    output:
        "TAIR10_{cas}.candidates.fq"
    params:
        pam = lambda wildcards: cas_params[wildcards.cas]["pam"],
        pamsite = lambda wildcards: ("--upstream-pam" if cas_params[wildcards.cas]["pamsite"] == "5p" else ""),
        cas = lambda wildcards: wildcards.cas
    shell:
        "python uniqmer2fq.py --cas {params.cas} --pam {params.pam} {params.pamsite} {input} > {output}" 

rule map_uniqmers_to_genome:
    input:
        "TAIR10_{cas}.candidates.fq" 
    output:
        "TAIR10_{cas}.candidates.bam"
    log:
        "TAIR10_{cas}.bbmap_genome.log"
    params:
        tmp_bam = lambda wildcards: "TAIR10_" + wildcards.cas + ".candidates.bam.tmp"
    threads:
        8
    shell:
        # "set +u && source jre-8u144 && source bbmap-38.06 && source samtools-1.7 && " + \
        "bbmap.sh in={input} out={params.tmp_bam} k=8" + \
        " minid=0.75 threads={threads} -Xmx4g ambig=all mappedonly=t secondary=t sssr=0.75 ssao=t saa=f" + \
        " mdtag=t nhtag=t xmtag=t amtag=t nmtag=t xstag=t indelfilter=0 subfilter=4 editfilter=4" + \
        " && samtools sort -@ {threads} -o {output} {params.tmp_bam}" + \
        " && samtools index {output}" + \
        " && rm {params.tmp_bam}" + \
        " 2> {log}"
   
rule filter_gene_uniqmers:
    input:
        "TAIR10_{cas}.candidates.bam"
    output:
        "TAIR10_{cas}.candidates.filtered.sam"
    log:
        "TAIR10_{cas}.genefilter.log"
    params:
        gene_gff = GENE_GFF
    threads:
        4
    shell:
        # "set +u && ml bedtools/2.27.1 && ml samtools/1.7 && " + \
        "bedtools intersect -wa -a {input} -b {params.gene_gff} | samtools sort -@ {threads} -n - | samtools view -h - > {output}"

rule offtarget_scan:
    input:
        "TAIR10_{cas}.candidates.filtered.sam"
    output:
        "{cas}.targets.gff",
        "{cas}.offtargets.gff",
        "{cas}.targets_offtargets.tsv"
    params:
        pam = lambda wildcards: cas_params[wildcards.cas]["pam"],
        pamsite = lambda wildcards: ("--upstream-pam" if cas_params[wildcards.cas]["pamsite"] == "5p" else ""),
        k = lambda wildcards: cas_params[wildcards.cas]["k"],
        cas = lambda wildcards: wildcards.cas
    shell:
        "python otscan.py --pam {params.pam} -k {params.k} --cas {params.cas} {params.pamsite} {input}"
        
#rule filter_gene_targets_old:
#    input:
#        targets_gff = "{cas}.targets.gff",
#        offtargets_gff = "{cas}.offtargets.gff",
#        full_annotation = "{cas}.targets_offtargets.tsv"
#    output:
#        "{cas}.targets_offtargets.genes.tsv"
#    params:
#        gene_gff = GENE_GFF,
#        gene_desc = GENE_DESC,
#        guidelen = lambda wildcards: cas_params[wildcards.cas]["k"]
#    shell:
#        "set +u && ml bedtools/2.27.1 && ml samtools/1.7" + \
#        " && bedtools intersect -wo -a {input.targets_gff} -b {params.gene_gff} | awk -v guidelen={params.guidelen} -v FS=\"\\t\" -v OFS=\"\\t\" '$19 == guidelen'" + \
#        " | cut -f 2- -d = | cut -f 1,8,10 | cut -f 1 -d \; | sed \"s/ID=//\" | sort -u > {input.targets_gff}.agi" + \
#        " && bedtools intersect -wo -a {input.offtargets_gff} -b {params.gene_gff} | awk -v guidelen={params.guidelen} -v FS=\"\\t\" -v OFS=\"\\t\" '$19 == guidelen'" + \
#        " | cut -f 2- -d = | cut -f 1,8,10 | cut -f 1 -d \; | sed \"s/ID=//\" | sort -u > {input.offtargets_gff}.agi" + \
#        " && join -1 3 -2 1 <(sort -k3,3 {input.targets_gff}.agi) <(sort -k1,1 {params.gene_desc}) | tr \" \" \"\\t\" > {input.targets_gff}.agi.desc" + \
#        " && join -1 3 -2 1 <(sort -k3,3 {input.offtargets_gff}.agi) <(sort -k1,1 {params.gene_desc}) | tr \" \" \"\\t\" > {input.offtargets_gff}.agi.desc" + \
#        " && join -1 1 -2 1 <(sort -k1,1 {input.full_annotation}) <(join -1 2 -2 2 <(sort -k2,2 {input.targets_gff}.agi.desc) <(sort -k2,2 {input.offtargets_gff}.agi.desc) | tr \" \" \"\\t\" | sort -k1,1) | tr \" \" \"\\t\" | sort -u > {output}" + \
#        "" #" && rm {input.targets_gff}.agi* {input.offtargets_gff}.agi*"


rule annotate_and_filter_gene_targets:
    input:
        targets_gff = "{cas}.targets.gff",              
        offtargets_gff = "{cas}.offtargets.gff",
        full_annotation = "{cas}.targets_offtargets.tsv"
    output:
        "{cas}.targets_offtargets.genes.tsv"
    params:
        gene_gff = GENE_GFF,
        gene_desc = GENE_DESC,
        guidelen = lambda wildcards: cas_params[wildcards.cas]["k"]
    shell:
        # "set +u && ml bedtools/2.27 && ml samtools/1.7 && " + \
        "awk -v FS=\"\\t\" -v OFS=\"\\t\" '{{ key=$1\":\"$2\":\"$3\":\"$4\":\"$5; print key,$0; }}' {input.full_annotation} | sort -u | sort -k1,1 > {input.full_annotation}.with_tkey && " + \
        "bedtools intersect -wo -a {input.targets_gff} -b {params.gene_gff} | awk -v guidelen={params.guidelen} -v FS=\"\t\" -v OFS=\"\t\" '{{ if ($19 == guidelen) {{ key=gensub(\"ID=\", \"\", \"g\", $9)\":\"$1\":\"$4\":\"$5\":\"$7;" + \
        " $18=substr($18,4,length($18)-4); print key,$0; }} }}' | sort -u | sort -k1,1 > {input.targets_gff}.with_geneid && " + \
        "join -1 1 -2 1 -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,2.19 {input.full_annotation}.with_tkey {input.targets_gff}.with_geneid | tr \" \" \"\\t\" | sort -u > {input.full_annotation}.targets_done" + \
        " && " + \
        "awk -v FS=\"\\t\" -v OFS=\"\\t\" '{{ key=$1\":\"$7\":\"$8\":\"$9\":\"$10; print key,$0; }}' {input.full_annotation}.targets_done | sort -u | sort -k1,1 > {input.full_annotation}.with_okey && " + \
        "bedtools intersect -wo -a {input.offtargets_gff} -b {params.gene_gff} | awk -v guidelen={params.guidelen} -v FS=\"\\t\" -v OFS=\"\\t\" '{{ if ($19 == guidelen) {{ key=gensub(\"ID=\", \"\", \"g\", $9)\":\"$1\":\"$4\":\"$5\":\"$7;" + \
        " $18=substr($18,4,length($18)-4); print key,$0; }} }}' | sort -u | sort -k1,1 > {input.offtargets_gff}.with_geneid && " + \
        "join -1 1 -2 1 -o 1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,2.19 {input.full_annotation}.with_okey {input.offtargets_gff}.with_geneid | tr \" \" \"\\t\" | sort -u > {input.full_annotation}.offtargets_done" + \
        " && " + \
        "join -1 16 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,2.2 <(sort -k16,16 {input.full_annotation}.offtargets_done) {params.gene_desc} | tr \" \" \"\\t\" | sort -k17,17 > {input.full_annotation}.tdesc && " + \
        "join -1 17 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,2.2 {input.full_annotation}.tdesc {params.gene_desc} | tr \" \" \"\\t\" | sort -k1,1 -k2,2 -k3,3 | sed \"s/function=//g\" > {output}" + \
        " && " + \
        "rm {input.full_annotation}.with_tkey {input.targets_gff}.with_geneid {input.full_annotation}.targets_done {input.full_annotation}.with_okey {input.offtargets_gff}.with_geneid {input.full_annotation}.offtargets_done {input.full_annotation}.tdesc"


# the following rule is necessary as the join-commands in filter_gene_targets will generate ambiguous rows in the table if there is no duplicate-filtering in extract_uniqmers 

if True: 
    rule duplicate_removal_patch:
        input:
            full_list = "{cas}.targets_offtargets.genes.tsv"
        output:
            filtered_list = "{cas}.targets_offtargets.genes.dedup.tsv"
        run:
            import sys
            import csv

            dupes = dict()
            with open(input.full_list) as _in, open(output.filtered_list, "wt") as _out:
                for row in csv.reader(_in, delimiter="\t"):
                    if dupes.get(row[5], row[0]) != row[0]:
                        # i.e. iff there is a record with the same target sequence but different target-id then ignore that record
                        # this only affects a minority of target sequences 
                        # if no guide+PAM candidate is found in forward direction then the reverse complement (rc) is checked
                        # in such a case, however, it is not tested whether the rc-sequence occurs in the jellyfish output as well,
                        # which will generate a duplicate target
                        # i have now added a duplicate check to uniqmer2fq.py, but for the Raitskin et al.-study we ran it with this patch-rule.
                        continue
                    print(*row, sep="\t", file=_out)
                    dupes[row[5]] = row[0]










