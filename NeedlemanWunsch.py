# By Zeal Jinwala
# Date: January 27, 2022

# This program is algorithm used in bioinformatics to align protein or nucleotide sequences. 
# It is an application of dynamic programming to compare biological sequences.
import pandas as pd
import numpy as np
from Bio import SeqIO

# Ask for seq 1, seq 2, match, mismatch, and gap
def parseArgs():
    seq1 = input("Enter the first sequence: ")
    seq2 = input("Enter the second sequence: ")
    match = input("Enter the match score: ")
    mismatch = input("Enter the mismatch score: ")
    gap = input("Enter the gap penalty: ")
    score = seqAlign(seq1, seq2, match, mismatch, gap)
    return score

def highestScore(subMat,gap):
    rows = len(subMat)
    cols = len(subMat[0])
    weight = subMat[1,1]
    print(rows)
    dynamicTable = np.zeros((rows,cols), dtype=int)
    for r in range(1,rows):
        for c in range(1,cols):
            curr = subMat[r,c]
            if r+c != 2:
                if c>1: 
                    fromLeft = dynamicTable[r,c-1] + gap
                if r>1: 
                    fromTop = dynamicTable[r-1,c] + gap
                if (r>1) and (c>1): 
                    topLeft = dynamicTable[r-1,c-1]
                    fromDiag = curr + topLeft
                # FIX ME: -----------------                
                if fromLeft and fromTop and fromDiag:
                    weight = max([fromLeft, fromTop, fromDiag])
                    del fromTop
                    del fromLeft
                    del fromDiag
                elif fromTop and fromDiag:
                    weight = fromLeft
                    del fromLeft
                elif fromLeft and fromDiag:
                    weight = fromTop
                    del fromTop
                # ----------------- 
            dynamicTable[r,c] = weight
    score = dynamicTable[rows, cols]
    

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
                subMat[i+1,j+1] = 1
            else:
                subMat[i+1,j+1] = 0
    highestScore(subMat, gap)
    

# def seqAlignBiopython():
#     print("Hello World 3!")

def main():
    parseArgs()

if __name__ == "__main__":
    main()