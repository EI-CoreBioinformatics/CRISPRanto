#!/bin/bash

function runCRISPResso {
 r1=$(basename $(dirname $1))/$(basename $1);
 echo $r1
 #r2=$(dirname $r1)/$(basename $r1 _R1.fastq.gz)_R2.fastq.gz;
 r2=$(dirname $r1)/$(basename $r1 _R1.bbduk.fq)_R2.bbduk.fq;

 echo $r2
 amplicon=$2
 guide=$3 
 docker run -v ${PWD}/oleg:/oleg -w /oleg -i lucapinello/crispresso CRISPResso -r1 $r1 -r2 $r2 -a $amplicon -g $guide $4 $5 --keep_intermediate               
}

function runCRISPRessoSingle {
 r1=$(basename $(dirname $1))/$(basename $1);
 amplicon=$2
 guide=$3 
 docker run -v ${PWD}/oleg:/oleg -w /oleg -i lucapinello/crispresso CRISPResso -r1 $r1 -a $amplicon -g $guide $4 $5                
}

amplicon_30=TTCTCTCTGTCACACCGATGTTTACTTCTGGGAAGCTAAGGTAGAGTAATCAATTTATTACACTCCAAATTCATAATCAAGTTCTAATTTTTTTAGAATTCTAATTTTTTATCTAAAAAAATTCAACCTTTTTGATTCCACAGGGACAAACACCGTTGTTTCCACGTATCTTCGGCCATGAAGCTGGAGGGTAATAGAAACACTAATCTTCTTTGCTTCGTTTTGGATATTTTTAAGGTTTTAGAGATTCAAGGTCGTTTTTTTTGTTGTTGTGTAGGATTGTTGAGAGTGTTGGAGAAGGAGTGACTGATCTTCAGCCA
# guide_30=TATCTTCGGCCATGAAGCTG
guide_30=ATCTTCGGCCATGAAGCTGG
 
amplicon_31=TTCTCTCTGTCACACCGATGTTTACTTCTGGGAAGCTAAGGTAGAGTAATCAATTTATTACACTCCAAATTCATAATCAAGTTCTAATTTTTTTAGAATTCTAATTTTTTATCTAAAAAAATTCAACCTTTTTGATTCCACAGGGACAAACACCGTTGTTTCCACGTATCTTCGGCCATGAAGCTGAGGGTAATAGAAACACTAATCTTCTTTGCTTCGTTTTGGATATTTTTAAGGTTTTAGAGATTCAAGGTCGTTTTTTTTGTTGTTGTGTAGGATTGTTGAGAGTGTTGGAGAAGGAGTGACTGATCTTCAGCCA
guide_31=TATCTTCGGCCATGAAGCTG
 
amplicon_32=TTCTCTCTGTCACACCGATGTTTACTTCTGGGAAGCTAAGGTAGAGTAATCAATTTATTACACTCCAAATTCATAATCAAGTTCTAATTTTTTTAGAATTCTAATTTTTTATCTAAAAAAATTCAACCTTTTTGATTCCACAGGGACAAACACCGTTGTTTCCACGTATCTTCGGCCATGAAGCTGAGAGTAATAGAAACACTAATCTTCTTTGCTTCGTTTTGGATATTTTTAAGGTTTTAGAGATTCAAGGTCGTTTTTTTTGTTGTTGTGTAGGATTGTTGAGAGTGTTGGAGAAGGAGTGACTGATCTTCAGCCA
guide_32=TATCTTCGGCCATGAAGCTG

amplicon_33=TTCTCTCTGTCACACCGATGTTTACTTCTGGGAAGCTAAGGTAGAGTAATCAATTTATTACACTCCAAATTCATAATCAAGTTCTAATTTTTTTAGAATTCTAATTTTTTATCTAAAAAAATTCAACCTTTTTGATTCCACAGGGACAAACACCGTTGTTTCACGTATCTTCGGCCATGAAGCTGGAGGGTAATAGAAACACTAATCTTCTTTGCTTCGTTTTGGATATTTTTAAGGTTTTAGAGATTCAAGGTCGTTTTTTTTGTTGTTGTGTAGGATTGTTGAGAGTGTTGGAGAAGGAGTGACTGATCTTCAGCCA
guide_33=ACGTATCTTCGGCCATGAAGCTG
 
amplicon_36=TTCTCTCTGTCACACCGATGTTTACTTCTGGGAAGCTAAGGTAGAGTAATCAATTTATTACACTCCAAATTCATAATCAAGTTCTAATTTTTTTAGAATTCTAATTTTTTATCTAAAAAAATTCAACCTTTTTGATTCCACAGGGACAAACACCGTTGTTTCCACGTATCTTCGGCCATGAAGCTGAGCGTAATAGAAACACTAATCTTCTTTGCTTCGTTTTGGATATTTTTAAGGTTTTAGAGATTCAAGGTCGTTTTTTTTGTTGTTGTGTAGGATTGTTGAGAGTGTTGGAGAAGGAGTGACTGATCTTCAGCCA
guide_36=TATCTTCGGCCATGAAGCTG





for r2 in $(ls oleg/*/*and_31*R2.bbduk.fq)
do
  echo $r2
  runCRISPRessoSingle $r2 $amplicon_31 $guide_31
done

for r2 in $(ls oleg/*/*and_32*R2.bbduk.fq)
do
  echo $r2
  runCRISPRessoSingle $r2 $amplicon_32 $guide_32
done

for r2 in $(ls oleg/*/*and_36*R2.bbduk.fq)
do
  echo $r2
  runCRISPRessoSingle $r2 $amplicon_36 $guide_36
done

exit






for r2 in $(ls oleg/Sample_30_level1*/*R2.bbduk.fq)
do
  echo $r2
  runCRISPRessoSingle $r2 $amplicon_30 $guide_30
done

for r2 in $(ls oleg/Sample_31_level1*/*R2.bbduk.fq)
do
  echo $r2
  runCRISPRessoSingle $r2 $amplicon_31 $guide_31
done

for r2 in $(ls oleg/Sample_32_level1*/*R2.bbduk.fq)
do
  echo $r2
  runCRISPRessoSingle $r2 $amplicon_32 $guide_32
done

for r2 in $(ls oleg/Sample_33_level1*/*R2.bbduk.fq)
do
  echo $r2
  runCRISPRessoSingle $r2 $amplicon_33 $guide_33 --cleavage_offset 1
done

for r2 in $(ls oleg/Sample_36_level1*/*R2.bbduk.fq)
do
  echo $r2
  runCRISPRessoSingle $r2 $amplicon_36 $guide_36
done
exit
for r2 in $(ls oleg/Sample_38_and_33*/*R2.bbduk.fq)
do
  echo $r2
  runCRISPRessoSingle $r2 $amplicon_33 $guide_33 --cleavage_offset 1
done

for r2 in $(ls oleg/Sample_39_and_33*/*R2.bbduk.fq)
do
  echo $r2
  runCRISPRessoSingle $r2 $amplicon_33 $guide_33 --cleavage_offset 1
done

for r2 in $(ls oleg/Sample_40_and_30*/*R2.bbduk.fq)
do
  echo $r2
  # runCRISPResso $r1 $amplicon_30 $guide_30
  runCRISPRessoSingle $r2 $amplicon_30 $guide_30
done

for r2 in $(ls oleg/Sample_*and_31*/*R2.fastq.gz)
do
  echo $r2 
  # runCRISPRessoSingle $r2 $amplicon_31 $guide_31
done

for r2 in $(ls oleg/Sample_*_and_32*/*R2.fastq.gz)
do
  echo $r2 
  # runCRISPRessoSingle $r2 $amplicon_32 $guide_32
done

for r2 in $(ls oleg/Sample_*_and_36*/*R2.fastq.gz)
do
  echo $r2 
  # runCRISPRessoSingle $r2 $amplicon_36 $guide_36
done                                      

exit




#for r1 in $(ls oleg/Sample_30_level1*/*R1.fastq.gz)
#do
#  echo $r1
#  runCRISPResso $r1 $amplicon_30 $guide_30
#done
#
#for r1 in $(ls oleg/Sample_31_level1*/*R1.fastq.gz)
#do
#  echo $r1
#  runCRISPResso $r1 $amplicon_31 $guide_31
#done
#for r1 in $(ls oleg/Sample_32_level1*/*R1.fastq.gz)
#do
#  echo $r1
#  runCRISPResso $r1 $amplicon_32 $guide_32
#done
#for r1 in $(ls oleg/Sample_33_level1*/*R1.fastq.gz)
#do
#  echo $r1
#  runCRISPResso $r1 $amplicon_33 $guide_33
#done
#for r1 in $(ls oleg/Sample_36_level1*/*R1.fastq.gz)
#do
#  echo $r1
#  runCRISPResso $r1 $amplicon_36 $guide_36
#done
#
#
#
#exit


for r1 in $(ls oleg/Sample_38_and_33*/*R1.fastq.gz)
do
  echo $r1
  runCRISPResso $r1 $amplicon_33 $guide_33
done

for r1 in $(ls oleg/Sample_39_and_33*/*R1.fastq.gz)
do
  echo $r1
  runCRISPResso $r1 $amplicon_33 $guide_33
done

for r1 in $(ls oleg/Sample_40_and_30*/*R1.fastq.gz)
do
  echo $r1
  runCRISPResso $r1 $amplicon_30 $guide_30
done

exit



for r1 in $(ls oleg/*/*and_31*R1.fastq.gz)
do
  echo $r1 
  runCRISPResso $r1 $amplicon_31 $guide_31
done

for r1 in $(ls oleg/*/*and_32*R1.fastq.gz)
do
  echo $r1 
  runCRISPResso $r1 $amplicon_32 $guide_32
done

for r1 in $(ls oleg/*/*and_36*R1.fastq.gz)
do
  echo $r1 
  runCRISPResso $r1 $amplicon_36 $guide_36
done                                      


exit 


















oleg/Sample_32_and_31_part3_rep1/32_and_31_part3_rep1_GTAGAGGACTCTCTAT_L001_R1.fastq.gz
oleg/Sample_32_and_31_part3_rep1/32_and_31_part3_rep1_GTAGAGGACTCTCTAT_L001_R2.fastq.gz
oleg/Sample_32_and_31_part3_rep2/32_and_31_part3_rep2_GCTCATGATCGACTAG_L001_R1.fastq.gz
oleg/Sample_32_and_31_part3_rep2/32_and_31_part3_rep2_GCTCATGATCGACTAG_L001_R2.fastq.gz
oleg/Sample_32_and_31_part3_rep3/32_and_31_part3_rep3_ACTCGCTACTCTCTAT_L001_R1.fastq.gz
oleg/Sample_32_and_31_part3_rep3/32_and_31_part3_rep3_ACTCGCTACTCTCTAT_L001_R2.fastq.gz
oleg/Sample_32_and_31_part3_rep4/32_and_31_part3_rep4_GTAGCTCCTCGACTAG_L001_R1.fastq.gz
oleg/Sample_32_and_31_part3_rep4/32_and_31_part3_rep4_GTAGCTCCTCGACTAG_L001_R2.fastq.gz
oleg/Sample_32_and_32_part3_rep1/32_and_32_part3_rep1_GTAGAGGATCGACTAG_L001_R1.fastq.gz
oleg/Sample_32_and_32_part3_rep1/32_and_32_part3_rep1_GTAGAGGATCGACTAG_L001_R2.fastq.gz
oleg/Sample_32_and_32_part3_rep2/32_and_32_part3_rep2_ATCTCAGGCTCTCTAT_L001_R1.fastq.gz
oleg/Sample_32_and_32_part3_rep2/32_and_32_part3_rep2_ATCTCAGGCTCTCTAT_L001_R2.fastq.gz
oleg/Sample_32_and_32_part3_rep3/32_and_32_part3_rep3_ACTCGCTATCGACTAG_L001_R1.fastq.gz
oleg/Sample_32_and_32_part3_rep3/32_and_32_part3_rep3_ACTCGCTATCGACTAG_L001_R2.fastq.gz
oleg/Sample_32_and_32_part3_rep4/32_and_32_part3_rep4_GCGTAGTACTCTCTAT_L001_R1.fastq.gz
oleg/Sample_32_and_32_part3_rep4/32_and_32_part3_rep4_GCGTAGTACTCTCTAT_L001_R2.fastq.gz
oleg/Sample_32_and_36_part3_rep1/32_and_36_part3_rep1_GCTCATGACTCTCTAT_L001_R1.fastq.gz
oleg/Sample_32_and_36_part3_rep1/32_and_36_part3_rep1_GCTCATGACTCTCTAT_L001_R2.fastq.gz
oleg/Sample_32_and_36_part3_rep2/32_and_36_part3_rep2_ATCTCAGGTCGACTAG_L001_R1.fastq.gz
oleg/Sample_32_and_36_part3_rep2/32_and_36_part3_rep2_ATCTCAGGTCGACTAG_L001_R2.fastq.gz
oleg/Sample_32_and_36_part3_rep3/32_and_36_part3_rep3_GTAGCTCCCTCTCTAT_L001_R1.fastq.gz
oleg/Sample_32_and_36_part3_rep3/32_and_36_part3_rep3_GTAGCTCCCTCTCTAT_L001_R2.fastq.gz
oleg/Sample_32_and_36_part3_rep4/32_and_36_part3_rep4_GCGTAGTATCGACTAG_L001_R1.fastq.gz
oleg/Sample_32_and_36_part3_rep4/32_and_36_part3_rep4_GCGTAGTATCGACTAG_L001_R2.fastq.gz
oleg/Sample_33_and_31_part3_rep1/33_and_31_part3_rep1_GTAGAGGATATCCTCT_L001_R1.fastq.gz
oleg/Sample_33_and_31_part3_rep1/33_and_31_part3_rep1_GTAGAGGATATCCTCT_L001_R2.fastq.gz
oleg/Sample_33_and_31_part3_rep2/33_and_31_part3_rep2_GCTCATGATTCTAGCT_L001_R1.fastq.gz
oleg/Sample_33_and_31_part3_rep2/33_and_31_part3_rep2_GCTCATGATTCTAGCT_L001_R2.fastq.gz
oleg/Sample_33_and_31_part3_rep3/33_and_31_part3_rep3_ACTCGCTATATCCTCT_L001_R1.fastq.gz
oleg/Sample_33_and_31_part3_rep3/33_and_31_part3_rep3_ACTCGCTATATCCTCT_L001_R2.fastq.gz
oleg/Sample_33_and_31_part3_rep4/33_and_31_part3_rep4_GTAGCTCCTTCTAGCT_L001_R1.fastq.gz
oleg/Sample_33_and_31_part3_rep4/33_and_31_part3_rep4_GTAGCTCCTTCTAGCT_L001_R2.fastq.gz
oleg/Sample_33_and_32_part3_rep1/33_and_32_part3_rep1_GTAGAGGATTCTAGCT_L001_R1.fastq.gz
oleg/Sample_33_and_32_part3_rep1/33_and_32_part3_rep1_GTAGAGGATTCTAGCT_L001_R2.fastq.gz
oleg/Sample_33_and_32_part3_rep2/33_and_32_part3_rep2_ATCTCAGGTATCCTCT_L001_R1.fastq.gz
oleg/Sample_33_and_32_part3_rep2/33_and_32_part3_rep2_ATCTCAGGTATCCTCT_L001_R2.fastq.gz
oleg/Sample_33_and_32_part3_rep3/33_and_32_part3_rep3_ACTCGCTATTCTAGCT_L001_R1.fastq.gz
oleg/Sample_33_and_32_part3_rep3/33_and_32_part3_rep3_ACTCGCTATTCTAGCT_L001_R2.fastq.gz
oleg/Sample_33_and_32_part3_rep4/33_and_32_part3_rep4_GCGTAGTATATCCTCT_L001_R1.fastq.gz
oleg/Sample_33_and_32_part3_rep4/33_and_32_part3_rep4_GCGTAGTATATCCTCT_L001_R2.fastq.gz
oleg/Sample_33_and_36_part3_rep1/33_and_36_part3_rep1_GCTCATGATATCCTCT_L001_R1.fastq.gz
oleg/Sample_33_and_36_part3_rep1/33_and_36_part3_rep1_GCTCATGATATCCTCT_L001_R2.fastq.gz
oleg/Sample_33_and_36_part3_rep2/33_and_36_part3_rep2_ATCTCAGGTTCTAGCT_L001_R1.fastq.gz
oleg/Sample_33_and_36_part3_rep2/33_and_36_part3_rep2_ATCTCAGGTTCTAGCT_L001_R2.fastq.gz
oleg/Sample_33_and_36_part3_rep3/33_and_36_part3_rep3_GTAGCTCCTATCCTCT_L001_R1.fastq.gz
oleg/Sample_33_and_36_part3_rep3/33_and_36_part3_rep3_GTAGCTCCTATCCTCT_L001_R2.fastq.gz
oleg/Sample_33_and_36_part3_rep4/33_and_36_part3_rep4_GCGTAGTATTCTAGCT_L001_R1.fastq.gz
oleg/Sample_33_and_36_part3_rep4/33_and_36_part3_rep4_GCGTAGTATTCTAGCT_L001_R2.fastq.gz
oleg/Sample_34_and_31_part3_rep1/34_and_31_part3_rep1_GTAGAGGAGTAAGGAG_L001_R1.fastq.gz
oleg/Sample_34_and_31_part3_rep1/34_and_31_part3_rep1_GTAGAGGAGTAAGGAG_L001_R2.fastq.gz
oleg/Sample_34_and_31_part3_rep2/34_and_31_part3_rep2_GCTCATGACCTAGAGT_L001_R1.fastq.gz
oleg/Sample_34_and_31_part3_rep2/34_and_31_part3_rep2_GCTCATGACCTAGAGT_L001_R2.fastq.gz
oleg/Sample_34_and_31_part3_rep3/34_and_31_part3_rep3_ACTCGCTAGTAAGGAG_L001_R1.fastq.gz
oleg/Sample_34_and_31_part3_rep3/34_and_31_part3_rep3_ACTCGCTAGTAAGGAG_L001_R2.fastq.gz
oleg/Sample_34_and_31_part3_rep4/34_and_31_part3_rep4_GTAGCTCCCCTAGAGT_L001_R1.fastq.gz
oleg/Sample_34_and_31_part3_rep4/34_and_31_part3_rep4_GTAGCTCCCCTAGAGT_L001_R2.fastq.gz
oleg/Sample_34_and_32_part3_rep1/34_and_32_part3_rep1_GTAGAGGACCTAGAGT_L001_R1.fastq.gz
oleg/Sample_34_and_32_part3_rep1/34_and_32_part3_rep1_GTAGAGGACCTAGAGT_L001_R2.fastq.gz
oleg/Sample_34_and_32_part3_rep2/34_and_32_part3_rep2_ATCTCAGGGTAAGGAG_L001_R1.fastq.gz
oleg/Sample_34_and_32_part3_rep2/34_and_32_part3_rep2_ATCTCAGGGTAAGGAG_L001_R2.fastq.gz
oleg/Sample_34_and_32_part3_rep3/34_and_32_part3_rep3_ACTCGCTACCTAGAGT_L001_R1.fastq.gz
oleg/Sample_34_and_32_part3_rep3/34_and_32_part3_rep3_ACTCGCTACCTAGAGT_L001_R2.fastq.gz
oleg/Sample_34_and_32_part3_rep4/34_and_32_part3_rep4_GCGTAGTAGTAAGGAG_L001_R1.fastq.gz
oleg/Sample_34_and_32_part3_rep4/34_and_32_part3_rep4_GCGTAGTAGTAAGGAG_L001_R2.fastq.gz
oleg/Sample_34_and_36_part3_rep1/34_and_36_part3_rep1_GCTCATGAGTAAGGAG_L001_R1.fastq.gz
oleg/Sample_34_and_36_part3_rep1/34_and_36_part3_rep1_GCTCATGAGTAAGGAG_L001_R2.fastq.gz
oleg/Sample_34_and_36_part3_rep2/34_and_36_part3_rep2_ATCTCAGGCCTAGAGT_L001_R1.fastq.gz
oleg/Sample_34_and_36_part3_rep2/34_and_36_part3_rep2_ATCTCAGGCCTAGAGT_L001_R2.fastq.gz
oleg/Sample_34_and_36_part3_rep3/34_and_36_part3_rep3_GTAGCTCCGTAAGGAG_L001_R1.fastq.gz
oleg/Sample_34_and_36_part3_rep3/34_and_36_part3_rep3_GTAGCTCCGTAAGGAG_L001_R2.fastq.gz
oleg/Sample_34_and_36_part3_rep4/34_and_36_part3_rep4_GCGTAGTACCTAGAGT_L001_R1.fastq.gz
oleg/Sample_34_and_36_part3_rep4/34_and_36_part3_rep4_GCGTAGTACCTAGAGT_L001_R2.fastq.gz
oleg/Sample_35_and_31_part3_rep1/35_and_31_part3_rep1_GTAGAGGAACTGCATA_L001_R1.fastq.gz
oleg/Sample_35_and_31_part3_rep1/35_and_31_part3_rep1_GTAGAGGAACTGCATA_L001_R2.fastq.gz
oleg/Sample_35_and_31_part3_rep2/35_and_31_part3_rep2_GCTCATGAGCGTAAGA_L001_R1.fastq.gz
oleg/Sample_35_and_31_part3_rep2/35_and_31_part3_rep2_GCTCATGAGCGTAAGA_L001_R2.fastq.gz
oleg/Sample_35_and_31_part3_rep3/35_and_31_part3_rep3_ACTCGCTAACTGCATA_L001_R1.fastq.gz
oleg/Sample_35_and_31_part3_rep3/35_and_31_part3_rep3_ACTCGCTAACTGCATA_L001_R2.fastq.gz
oleg/Sample_35_and_31_part3_rep4/35_and_31_part3_rep4_GTAGCTCCGCGTAAGA_L001_R1.fastq.gz
oleg/Sample_35_and_31_part3_rep4/35_and_31_part3_rep4_GTAGCTCCGCGTAAGA_L001_R2.fastq.gz
oleg/Sample_35_and_32_part3_rep1/35_and_32_part3_rep1_GTAGAGGAGCGTAAGA_L001_R1.fastq.gz
oleg/Sample_35_and_32_part3_rep1/35_and_32_part3_rep1_GTAGAGGAGCGTAAGA_L001_R2.fastq.gz
oleg/Sample_35_and_32_part3_rep2/35_and_32_part3_rep2_ATCTCAGGACTGCATA_L001_R1.fastq.gz
oleg/Sample_35_and_32_part3_rep2/35_and_32_part3_rep2_ATCTCAGGACTGCATA_L001_R2.fastq.gz
oleg/Sample_35_and_32_part3_rep3/35_and_32_part3_rep3_ACTCGCTAGCGTAAGA_L001_R1.fastq.gz
oleg/Sample_35_and_32_part3_rep3/35_and_32_part3_rep3_ACTCGCTAGCGTAAGA_L001_R2.fastq.gz
oleg/Sample_35_and_32_part3_rep4/35_and_32_part3_rep4_GCGTAGTAACTGCATA_L001_R1.fastq.gz
oleg/Sample_35_and_32_part3_rep4/35_and_32_part3_rep4_GCGTAGTAACTGCATA_L001_R2.fastq.gz
oleg/Sample_35_and_36_part3_rep1/35_and_36_part3_rep1_GCTCATGAACTGCATA_L001_R1.fastq.gz
oleg/Sample_35_and_36_part3_rep1/35_and_36_part3_rep1_GCTCATGAACTGCATA_L001_R2.fastq.gz
oleg/Sample_35_and_36_part3_rep2/35_and_36_part3_rep2_ATCTCAGGGCGTAAGA_L001_R1.fastq.gz
oleg/Sample_35_and_36_part3_rep2/35_and_36_part3_rep2_ATCTCAGGGCGTAAGA_L001_R2.fastq.gz
oleg/Sample_35_and_36_part3_rep3/35_and_36_part3_rep3_GTAGCTCCACTGCATA_L001_R1.fastq.gz
oleg/Sample_35_and_36_part3_rep3/35_and_36_part3_rep3_GTAGCTCCACTGCATA_L001_R2.fastq.gz
oleg/Sample_35_and_36_part3_rep4/35_and_36_part3_rep4_GCGTAGTAGCGTAAGA_L001_R1.fastq.gz
oleg/Sample_35_and_36_part3_rep4/35_and_36_part3_rep4_GCGTAGTAGCGTAAGA_L001_R2.fastq.gz
oleg/Sample_36_and_31_part3_rep1/36_and_31_part3_rep1_GTAGAGGAAAGGAGTA_L001_R1.fastq.gz
oleg/Sample_36_and_31_part3_rep1/36_and_31_part3_rep1_GTAGAGGAAAGGAGTA_L001_R2.fastq.gz
oleg/Sample_36_and_31_part3_rep2/36_and_31_part3_rep2_GCTCATGACTATTAAG_L001_R1.fastq.gz
oleg/Sample_36_and_31_part3_rep2/36_and_31_part3_rep2_GCTCATGACTATTAAG_L001_R2.fastq.gz
oleg/Sample_36_and_31_part3_rep3/36_and_31_part3_rep3_ACTCGCTAAAGGAGTA_L001_R1.fastq.gz
oleg/Sample_36_and_31_part3_rep3/36_and_31_part3_rep3_ACTCGCTAAAGGAGTA_L001_R2.fastq.gz
oleg/Sample_36_and_31_part3_rep4/36_and_31_part3_rep4_GTAGCTCCCTATTAAG_L001_R1.fastq.gz
oleg/Sample_36_and_31_part3_rep4/36_and_31_part3_rep4_GTAGCTCCCTATTAAG_L001_R2.fastq.gz
oleg/Sample_36_and_32_part3_rep1/36_and_32_part3_rep1_GTAGAGGACTATTAAG_L001_R1.fastq.gz
oleg/Sample_36_and_32_part3_rep1/36_and_32_part3_rep1_GTAGAGGACTATTAAG_L001_R2.fastq.gz
oleg/Sample_36_and_32_part3_rep2/36_and_32_part3_rep2_ATCTCAGGAAGGAGTA_L001_R1.fastq.gz
oleg/Sample_36_and_32_part3_rep2/36_and_32_part3_rep2_ATCTCAGGAAGGAGTA_L001_R2.fastq.gz
oleg/Sample_36_and_32_part3_rep3/36_and_32_part3_rep3_ACTCGCTACTATTAAG_L001_R1.fastq.gz
oleg/Sample_36_and_32_part3_rep3/36_and_32_part3_rep3_ACTCGCTACTATTAAG_L001_R2.fastq.gz
oleg/Sample_36_and_32_part3_rep4/36_and_32_part3_rep4_GCGTAGTAAAGGAGTA_L001_R1.fastq.gz
oleg/Sample_36_and_32_part3_rep4/36_and_32_part3_rep4_GCGTAGTAAAGGAGTA_L001_R2.fastq.gz
oleg/Sample_36_and_36_part3_rep1/36_and_36_part3_rep1_GCTCATGAAAGGAGTA_L001_R1.fastq.gz
oleg/Sample_36_and_36_part3_rep1/36_and_36_part3_rep1_GCTCATGAAAGGAGTA_L001_R2.fastq.gz
oleg/Sample_36_and_36_part3_rep2/36_and_36_part3_rep2_ATCTCAGGCTATTAAG_L001_R1.fastq.gz
oleg/Sample_36_and_36_part3_rep2/36_and_36_part3_rep2_ATCTCAGGCTATTAAG_L001_R2.fastq.gz
oleg/Sample_36_and_36_part3_rep3/36_and_36_part3_rep3_GTAGCTCCAAGGAGTA_L001_R1.fastq.gz
oleg/Sample_36_and_36_part3_rep3/36_and_36_part3_rep3_GTAGCTCCAAGGAGTA_L001_R2.fastq.gz
oleg/Sample_36_and_36_part3_rep4/36_and_36_part3_rep4_GCGTAGTACTATTAAG_L001_R1.fastq.gz
oleg/Sample_36_and_36_part3_rep4/36_and_36_part3_rep4_GCGTAGTACTATTAAG_L001_R2.fastq.gz
oleg/Sample_37_and_31_part3_rep1/37_and_31_part3_rep1_GTAGAGGACTAAGCCT_L001_R1.fastq.gz
oleg/Sample_37_and_31_part3_rep1/37_and_31_part3_rep1_GTAGAGGACTAAGCCT_L001_R2.fastq.gz
oleg/Sample_37_and_31_part3_rep2/37_and_31_part3_rep2_GCTCATGAAAGGCTAT_L001_R1.fastq.gz
oleg/Sample_37_and_31_part3_rep2/37_and_31_part3_rep2_GCTCATGAAAGGCTAT_L001_R2.fastq.gz
oleg/Sample_37_and_31_part3_rep3/37_and_31_part3_rep3_ACTCGCTACTAAGCCT_L001_R1.fastq.gz
oleg/Sample_37_and_31_part3_rep3/37_and_31_part3_rep3_ACTCGCTACTAAGCCT_L001_R2.fastq.gz
oleg/Sample_37_and_31_part3_rep4/37_and_31_part3_rep4_GTAGCTCCAAGGCTAT_L001_R1.fastq.gz
oleg/Sample_37_and_31_part3_rep4/37_and_31_part3_rep4_GTAGCTCCAAGGCTAT_L001_R2.fastq.gz
oleg/Sample_37_and_32_part3_rep1/37_and_32_part3_rep1_GTAGAGGAAAGGCTAT_L001_R1.fastq.gz
oleg/Sample_37_and_32_part3_rep1/37_and_32_part3_rep1_GTAGAGGAAAGGCTAT_L001_R2.fastq.gz
oleg/Sample_37_and_32_part3_rep2/37_and_32_part3_rep2_ATCTCAGGCTAAGCCT_L001_R1.fastq.gz
oleg/Sample_37_and_32_part3_rep2/37_and_32_part3_rep2_ATCTCAGGCTAAGCCT_L001_R2.fastq.gz
oleg/Sample_37_and_32_part3_rep3/37_and_32_part3_rep3_ACTCGCTAAAGGCTAT_L001_R1.fastq.gz
oleg/Sample_37_and_32_part3_rep3/37_and_32_part3_rep3_ACTCGCTAAAGGCTAT_L001_R2.fastq.gz
oleg/Sample_37_and_32_part3_rep4/37_and_32_part3_rep4_GCGTAGTACTAAGCCT_L001_R1.fastq.gz
oleg/Sample_37_and_32_part3_rep4/37_and_32_part3_rep4_GCGTAGTACTAAGCCT_L001_R2.fastq.gz
oleg/Sample_37_and_36_part3_rep1/37_and_36_part3_rep1_GCTCATGACTAAGCCT_L001_R1.fastq.gz
oleg/Sample_37_and_36_part3_rep1/37_and_36_part3_rep1_GCTCATGACTAAGCCT_L001_R2.fastq.gz
oleg/Sample_37_and_36_part3_rep2/37_and_36_part3_rep2_ATCTCAGGAAGGCTAT_L001_R1.fastq.gz
oleg/Sample_37_and_36_part3_rep2/37_and_36_part3_rep2_ATCTCAGGAAGGCTAT_L001_R2.fastq.gz
oleg/Sample_37_and_36_part3_rep3/37_and_36_part3_rep3_GTAGCTCCCTAAGCCT_L001_R1.fastq.gz
oleg/Sample_37_and_36_part3_rep3/37_and_36_part3_rep3_GTAGCTCCCTAAGCCT_L001_R2.fastq.gz
oleg/Sample_37_and_36_part3_rep4/37_and_36_part3_rep4_GCGTAGTAAAGGCTAT_L001_R1.fastq.gz
oleg/Sample_37_and_36_part3_rep4/37_and_36_part3_rep4_GCGTAGTAAAGGCTAT_L001_R2.fastq.gz

















