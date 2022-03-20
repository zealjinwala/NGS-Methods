# By Zeal Jinwala
# Date: January 27, 2022
# This program is an algorithm used in bioinformatics to GLOBALY align protein or nucleotide sequences. 
# It is an application of dynamic programming to compare biological sequences.

# imoprting libraries
import numpy as np
from Bio import pairwise2
import argparse

# dynamic proogramming function
def highestScore(subMat,gap):
    rows = len(subMat)
    cols = len(subMat[0])
    dynamicTable = np.zeros((rows,cols),dtype=int)
    top = gap + gap*max(rows,cols)
    diag = gap + gap*max(rows,cols)
    for r in range(rows):
        for c in range(cols):
            if r + c > 0:
                curr = subMat[r,c]
                if c>0 and r>0:
                    diag = curr + dynamicTable[r-1,c-1]
                if r>0:
                    top = dynamicTable[r-1,c] + gap
                if c>0:
                    left = dynamicTable[r,c-1] + gap
                maxScore = max(diag,top,left)
                dynamicTable[r,c] = maxScore
    score = dynamicTable[rows-1, cols-1]
    return score 
    
def seqAlign(seq1, seq2, match, mismatch, gap):
    if not match:
        match = 5
    if not mismatch:
        mismatch = -4
    if not gap:
        gap = -5
    match = int(match)
    mismatch = int(mismatch)
    gap = int(gap)
    lenSeq1 = len(seq1)
    lenSeq2 = len(seq2)
    subMat = np.zeros((lenSeq1+1,lenSeq2+1), dtype=int)
    for i in range(lenSeq1):
        for j in range(lenSeq2):
            if seq1[i] == seq2[j]: 
                subMat[i+1,j+1] = match
            else:
                subMat[i+1,j+1] = mismatch
    score = highestScore(subMat, gap)
    return score

# comparing with Biopyhton's tools
def seqAlignBiopython(seq1, seq2, match, mismatch, gap):
    if not match:
        match = 5
    if not mismatch:
        mismatch = -4
    if not gap:
        gap = -5
    match = int(match)
    mismatch = int(mismatch)
    gap = int(gap)
    alignments = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap, gap)
    return alignments


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sequence1",
        type=str,
        required=True,
        help="First Sequence",
    )

    parser.add_argument(
        "--sequence2",
        type=str,
        required=True,
        help="Second Sequence",
    )

    parser.add_argument(
        "--match",
        type=int,
        required=False,
        help="Match Score",
    )

    parser.add_argument(
        "--mismatch",
        type=int,
        required=False,
        help="Mismatch score",
    )

    parser.add_argument(
        "--gap",
        type=int,
        required=False,
        help="Gap score",
    )

    args = parser.parse_args()
    score = seqAlign(
        args.sequence1, 
        args.sequence2, 
        args.match, 
        args.mismatch, 
        args.gap
        )
    print("Alignment score:", score)