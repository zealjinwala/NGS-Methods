# By Zeal Jinwala
# Date: January 27, 2022
# This program is algorithm used in bioinformatics to GLOBALY align protein or nucleotide sequences. 
# It is an application of dynamic programming to compare biological sequences.

import pandas as pd
import numpy as np
from Bio import SeqIO

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

# def seqAlignBiopython():
#     print("Hello World 3!")

def main():
    seq1 = input("Enter the first sequence: ")
    seq2 = input("Enter the second sequence: ")
    match = input("Enter the match score: ")
    mismatch = input("Enter the mismatch score: ")
    gap = input("Enter the gap penalty: ")
    score = seqAlign(seq1, seq2, match, mismatch, gap)
    print(score)


if __name__ == "__main__":
    main()
