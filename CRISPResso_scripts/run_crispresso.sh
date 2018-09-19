#!/bin/bash

function runCRISPRessoSingle {
 r1=$(basename $(dirname $1))/$(basename $1);
 amplicon=$2
 guide=$3
 docker run -v ${PWD}/oleg:/oleg -w /oleg -i lucapinello/crispresso CRISPResso -r1 $r1 -a $amplicon -g $guide $4 $5
}

amplicon=GATAAGGAATTTTGCATAGTATTTAGGTTCACAAGTGGGACAATCTTCTTACACTGAAATATCTTTATGTCAGGCTTAATTTACTGCTATCTTGTTCAATAAAATGCCCCAAATTGGACTTGTTTCTGCCGTTAATTTGAGAGTCCAAGGTAATTCAGCTTATCTTTGGAGCTCGAGGTCTTCGTTGGGAACTGAAAGTCAAGATGTTTGCTTGCAAAGGAATTTGTTATGTTTTGGTAGTAGCGACTCCATGGGGCATAAGTTAAGGATTCGTACTCCAAGTGCCACGACCCGAAGATTGACAAAGGACTTTAATCCT
guide=CCGTTAATTTGAGAGTCCA

for r1 in oleg/Sample_*part[12]*/*bbduk.fq; do
 runCRISPRessoSingle $r1 $amplicon $guide
done
exit

# for r1 in oleg/Sample_*_part2_rep*/*R1.fastq.gz
# for r1 in oleg/Sample_ctrl*/*R1.fastq.gz
for r1 in oleg/Sample_*part1*/*R1.fastq.gz
do
  r1=$(basename $(dirname $r1))/$(basename $r1);
  echo $r1
  r2=$(dirname $r1)/$(basename $r1 _R1.fastq.gz)_R2.fastq.gz;
  # r2=$(echo $r1 | gsed "s/
  echo $r2
  docker run -v ${PWD}/oleg:/oleg -w /oleg -i lucapinello/crispresso CRISPResso -r1 $r1 -r2 $r2 -a GATAAGGAATTTTGCATAGTATTTAGGTTCACAAGTGGGACAATCTTCTTACACTGAAATATCTTTATGTCAGGCTTAATTTACTGCTATCTTGTTCAATAAAATGCCCCAAATTGGACTTGTTTCTGCCGTTAATTTGAGAGTCCAAGGTAATTCAGCTTATCTTTGGAGCTCGAGGTCTTCGTTGGGAACTGAAAGTCAAGATGTTTGCTTGCAAAGGAATTTGTTATGTTTTGGTAGTAGCGACTCCATGGGGCATAAGTTAAGGATTCGTACTCCAAGTGCCACGACCCGAAGATTGACAAAGGACTTTAATCCT -g CCGTTAATTTGAGAGTCCA 
  # break
done
exit
Sample_ctrl__1/ctrl__1_AAGAGGCAAAGGAGTA_L001_R1.fastq.gz
Sample_ctrl__1/ctrl__1_AAGAGGCAAAGGAGTA_L001_R2.fastq.gz
Sample_ctrl__2/ctrl__2_AAGAGGCACTAAGCCT_L001_R1.fastq.gz
Sample_ctrl__2/ctrl__2_AAGAGGCACTAAGCCT_L001_R2.fastq.gz
Sample_ctrl__3/ctrl__3_AAGAGGCACGTCTAAT_L001_R1.fastq.gz
Sample_ctrl__3/ctrl__3_AAGAGGCACGTCTAAT_L001_R2.fastq.gz
Sample_ctrl__4/ctrl__4_AAGAGGCATCTCTCCG_L001_R1.fastq.gz
Sample_ctrl__4/ctrl__4_AAGAGGCATCTCTCCG_L001_R2.fastq.gz
Sample_ctrl__5/ctrl__5_GTAGAGGATCTCTCCG_L001_R1.fastq.gz
Sample_ctrl__5/ctrl__5_GTAGAGGATCTCTCCG_L001_R2.fastq.gz
Sample_ctrl_part1__1/ctrl_part1__1_TACGCTGCAAGGAGTA_L001_R1.fastq.gz
Sample_ctrl_part1__1/ctrl_part1__1_TACGCTGCAAGGAGTA_L001_R2.fastq.gz
Sample_ctrl_part1__2/ctrl_part1__2_TACGCTGCCTAAGCCT_L001_R1.fastq.gz
Sample_ctrl_part1__2/ctrl_part1__2_TACGCTGCCTAAGCCT_L001_R2.fastq.gz
Sample_ctrl_part1__3/ctrl_part1__3_TACGCTGCCGTCTAAT_L001_R1.fastq.gz
Sample_ctrl_part1__3/ctrl_part1__3_TACGCTGCCGTCTAAT_L001_R2.fastq.gz
Sample_ctrl_part1__4/ctrl_part1__4_TACGCTGCTCTCTCCG_L001_R1.fastq.gz
Sample_ctrl_part1__4/ctrl_part1__4_TACGCTGCTCTCTCCG_L001_R2.fastq.gz
