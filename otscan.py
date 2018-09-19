import sys
import csv
import re
from argparse import ArgumentParser

# OTScan - (offtarget)-scan
# CRISPRanto - pipeline for CRISPR offtarget mining
# authors: Christian Schudoma, Oleg Raitskin
# contact: christian.schudoma@earlham.ac.uk
# (C) 2017-2018 Earlham Institute
# License: MIT

"""
SaCas9_target407_+	16	Chr2	16981270	1	27=	*	0	0	ACCCATCAAAATCTCATAACATGTCGT	IIIIIIIIIIIIIIIIIIIIIIIIIII	XT:A:R	NM:i:0AM:i:1	XM:i:2	MD:Z:27	NH:i:2
SaCas9_target407_+	256	Chr4	17943700	1	6=1X20=	*	0	0	ACGACATGTTATGAGATTTTGATGGGT	IIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:1	AM:i:1XM:i:2	MD:Z:6G20	NH:i:2
"""
IUPAC_CODE = {"A": "A", "C": "C", "G": "G", "T": "T",
              "N": "[ACGT]",
              "R": "[AG]", "Y": "[CT]",
              "W": "[AT]", "S": "[CG]",
              "M": "[AC]", "K": "[GT]",
              "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]"}


def compilePAM(pam, is_upstream=False):
    up = "^" if is_upstream else ""
    down = "$" if not is_upstream else ""
    return re.compile(up + "".join(IUPAC_CODE[base] for base in pam) + down)

from uniqmer2fq import reverseComplement #, compilePAM, checkGuide

def makeOffTargetSequence(seq, MD):
    saveMD=MD
    oseq = ""
    p = 0
    while MD:
        matches = re.match("[0-9]+", MD)
        if not matches:
            raise ValueError("MD contains weird characters: " + MD + "," + saveMD) 
        m = int(matches.group())
        oseq += seq[p:p+m]
        p += m
        pp = int(matches.end())
        while MD[pp:] and not MD[pp].isdigit():
            oseq += MD[pp] 
            pp += 1
            p += 1
           
        MD = MD[pp:]
    return oseq


def processHits(hits, outfile, tgff, ogff, pam, rpam, upstream_pam, tlen, plen, allowed_mismatches):
    def tag2dict(tags):
        return dict((tag.split(":")[0], tag.split(":")[2]) for tag in tags)
    def locateMutations(cigar):
        mpos, p = list(), 1
        while cigar:
            len_op = re.match("[0-9]+[^0-9]", cigar)    
            length, op = int(len_op.group()[:-1]), len_op.group()[-1]
            if op == "X":
                mpos.extend(range(p, p + length))                
            p += length
            cigar = cigar[len(len_op.group()):] 
        return mpos

    # evaluate the tag fields in order to find the target
    # or eliminate the current guide+PAM if there
    # are more than 1 perfect hits
    tags = list(map(tag2dict, (hit[11:] for hit in hits)))
    try:
        has_no_mismatches = list(map(lambda x:int(x.get("NM")) == 0, tags))
    except:
        raise ValueError("Hit missing NM-tag: {}".format(str(tags)))

    perfect_hits = list(filter(lambda x:x, has_no_mismatches))
    if len(perfect_hits) > 1:
        pass # used to ignore those, keeping this check in case it will become relevant again
    if not perfect_hits:
        return None # something went wrong between CDS/exon and genomic sequence! -> ignore

    #Figure out which of the hits is the target - should be first hit
    #since all others are secondary
    # as of 2018-09-11 this has been relaxed to include perfect matches,
    # alignments should be ordered by genomic position, so "first" perfect alignment of identical sequences
    # will be considered as target
    try:
        target_index = has_no_mismatches.index(True)
    except:
        return None 
    target, target_tags = hits[target_index], tags[target_index]
    target_strand = "+"
    if int(target[1]) & 16:
        target[9] = reverseComplement(target[9])
        target_strand = "-"
    target_line = [target[0], target[2], int(target[3]), int(target[3]) + tlen, target_strand, target[9]]

    for i, otdata in enumerate(zip(hits, tags)):
        offtarget, offtarget_tags = otdata
        if i != target_index and int(offtarget_tags.get("NM")) <= allowed_mismatches:
            mutations = locateMutations(offtarget[5])
            offtarget_strand = "+"
            offtarget[9] = makeOffTargetSequence(offtarget[9], offtarget_tags.get("MD"))

            if int(offtarget[1]) & 16:
                mutations = [tlen - m + 1 for m in mutations][::-1]
                offtarget_strand = "-"
                offtarget[9] = reverseComplement(offtarget[9])

            dh3 = 0
            if upstream_pam:
                # pam is at 5'
                pam_mutations = [m for m in mutations if m <= plen]
                check_pam = pam.search(offtarget[9])
                if not pam_mutations or pam.search(offtarget[9]): 
                    dh3 = len([m for m in mutations if m >= tlen - 2])
                else:
                    continue
            else:
                # pam is at 3'
                pam_mutations = [m for m in mutations if m >= tlen - plen + 1]
                check_pam = pam.search(offtarget[9])
                if not pam_mutations or pam.search(offtarget[9]):
                    dh3 = len([m for m in mutations if m <= 3])
                else:
                    continue
            print("MUTATIONS", mutations)
            print("PAM MUTATIONS", pam_mutations)
            dh = len(mutations) - len(pam_mutations)
            if dh > 1:
                continue

            offtarget_line = [offtarget[2], int(offtarget[3]), int(offtarget[3]) + tlen, offtarget_strand, offtarget[9], ",".join(map(str, mutations)) if mutations else "perfect", dh, dh3, ",".join(map(str, pam_mutations)) if pam_mutations else "n/a", str(check_pam.group())]
            print(*target_line, *offtarget_line, sep="\t", file=outfile)
            # Chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
            print(target[2], ".", ".", int(target[3]), int(target[3]) + tlen, ".", "+" if int(target[1]) & 16 == 0 else "-", ".", "ID=" + target[0], sep="\t", file=tgff)
            print(offtarget[2], ".", ".", int(offtarget[3]), int(offtarget[3]) + tlen, ".", "+" if int(offtarget[1]) & 16 == 0 else "-", ".", "ID=" + target[0], sep="\t", file=ogff)

            print(*target_line, *offtarget_line, sep="\t")
  
            
    return target, offtarget

        


if __name__ == "__main__":
    ap = ArgumentParser()
    ap.add_argument("input_sam", type=str)
    ap.add_argument("--pam", type=str)
    ap.add_argument("--upstream-pam", action="store_true")
    ap.add_argument("--target-length", "-k", type=int) 
    ap.add_argument("--cas", type=str)
    args = ap.parse_args()

    pam = compilePAM(args.pam, is_upstream=args.upstream_pam)
    rpam = compilePAM(reverseComplement(args.pam), is_upstream=args.upstream_pam) 

    allowed_mismatches = len(re.sub("[ACGT]", "", args.pam)) + 1

    last = list()
    with open(args.cas + ".targets_offtargets.tsv", "w") as outfile, open(args.cas + ".targets.gff", "w") as gff_t, open(args.cas + ".offtargets.gff", "w") as gff_ot:
        for row in csv.reader(open(args.input_sam), delimiter="\t"):
            if row[0].startswith("@"):
                continue
            if row[2] == "ChrM" or row[2] == "ChrC":
                continue
            if not last or last[0][0] == row[0]:
                pass
            else:
                if len(last) > 1:
                    res = processHits(last, outfile, gff_t, gff_ot, pam, rpam, args.upstream_pam, args.target_length, len(args.pam), allowed_mismatches)
                last = list()
            last.append(row)
    
        if len(last) > 1:
            res = processHits(last, outfile, gff_t, gff_ot, pam, rpam, args.upstream_pam, args.target_length, len(args.pam), allowed_mismatches)

