import sys
import re
from argparse import ArgumentParser

from ktio.ktio import readFasta


# uniqmer2fq - Find CRISPR target candidates in exons
# CRISPRanto - pipeline for CRISPR offtarget mining
# authors: Christian Schudoma, Oleg Raitskin
# contact: christian.schudoma@earlham.ac.uk
# (C) 2017-2018 Earlham Institute
# License: MIT


def reverseComplement(seq, alphabet='ACGT'):
    """
    Returns the reverse complement of nucleic acid sequence input.
    """
    compl= dict(zip('ACGTNRYWSMKBHDV', 'TGCANYRWSKMVDHB'))
    return ''.join([compl[base]
                    for base in seq.upper().replace('U', 'T')])[::-1]

IUPAC_CODE = {"A": "A", "C": "C", "G": "G", "T": "T", 
              "N": "[ACGT]",
              "R": "[AG]", "Y": "[CT]",
              "W": "[AT]", "S": "[CG]",
              "M": "[AC]", "K": "[GT]",
              "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]"}


def compilePAM(pam):
    return re.compile("".join(IUPAC_CODE[base] for base in pam))
def validGuide(guide, plen, pam, rpam, is_upstream_pam):
    if is_upstream_pam:
        return pam.match(guide[:plen]) or rpam.match(guide[-plen:])
    else:
        return pam.match(guide[-plen:]) or rpam.match(guide[:plen])
def checkGuide(seq, plen, pam, rpam, is_upstream_pam):
    """ replaces validGuide, convert sequence to +-strand if pam is found on --strand """
    if is_upstream_pam:
        if pam.match(seq[:plen]):
            yield seq, "+"
        if rpam.match(seq[-plen:]):
            yield reverseComplement(seq), "-"
    else:
        if pam.match(seq[-plen:]):
            yield seq, "+"
        if rpam.match(seq[:plen]):
            yield reverseComplement(seq), "-"
    #yield "", ""
    



if __name__ == "__main__":
    ap = ArgumentParser()
    ap.add_argument("kmer_counts", type=str)
    ap.add_argument("--pam", "-p", type=str, default="NGG")
    ap.add_argument("--upstream-pam", "-u", action="store_true")
    ap.add_argument("--cas", "-c", type=str, default="SpCas9-wt")
    ap.add_argument("--duplicate-check", action="store_true")  # added after Raitskin et al analysis
    args = ap.parse_args()

    seqct = 0
    seen = set()
 
    pam_compiled, rpam_compiled = map(compilePAM, (args.pam, reverseComplement(args.pam)))
    for _id, _seq in readFasta(args.kmer_counts, headless=True):
        if True: # int(_id) == 1:
            for target, strand in checkGuide(_seq, len(args.pam), pam_compiled, rpam_compiled, args.upstream_pam):
            # if validGuide(_seq, len(args.pam), pam_compiled, rpam_compiled, args.upstream_pam): 
                if target:
                    if not args.duplicate_check or target not in seen:
                        seqct += 1
                        print("@{}_target{}_{}\n{}\n+\n{}".format(args.cas, seqct, strand, target, "I" * len(_seq)))
                        if args.duplicate_check: # do not need to allocate memory if duplicate-check is False
                            seen.add(target)
