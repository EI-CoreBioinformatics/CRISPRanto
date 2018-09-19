import sys

from collections import Counter

from ktio.ktio import readFastq

unique_candidates = Counter()

for _id,_seq,_qual in readFastq(sys.argv[1]):    
    unique_candidates[_seq] += 1

print(sys.argv[1], sum(unique_candidates.values()), len(unique_candidates), sep="\t")
