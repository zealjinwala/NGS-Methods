
# By Zeal Jinwala 
# Date: March 19, 2022

# Description:
# This program take an XML file (BLAST search results) as input and answers the following questions:
#   1) Prints the name of the top protein hit.
#   2) Prints a unique list of species names of all the hits.
#   3) Shows a histogram of percent identities of the first HSPs from all the hits.

# Query Protein: <br>
# SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIG
# HSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVS
# FCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKY
# NYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQ

# import argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def parseXML(filepath):
    '''takes in the XML(BLAST result), outputs list in csv format'''
    from Bio.Blast import NCBIXML
    resultHandle = open(filepath)
    blastRecords = NCBIXML.parse(resultHandle)
    df = pd.DataFrame(columns = [  'Query',
                'Name', 'Length', 'Score', 'Expect',
                'Query Seq', 'Subject Sequence'
            ])
    count = 0
    for blastRecord in blastRecords:
        count +=1
        query = blastRecord.query
        for alignment in blastRecord.alignments:
            name = alignment.title
            length = alignment.length
            hsp = alignment.hsps[0]
            score = hsp.score
            expect = hsp.expect
            querySeq = hsp.query
            subSeq = hsp.sbjct
            df = df.append({'Query':query,
                'Name':name , 'Length': length,'Score': score, 'Expect':expect,'Query Seq':querySeq, 'Subject Sequence':subSeq}, ignore_index = True)
    resultHandle.close()
    return df

def getInfo(blastResults):
    '''takes in list of results, outputs results '''
    # get the top hit
    hitName = str(blastResults['Name'].iloc[0])
    sep = ">"
    ids = [pos for pos, char in enumerate(hitName) if char == sep]
    topHit = hitName[ids[0]+1:ids[1]]
    # unique list of all proteins found
    descriptions = blastResults.Name
    uniqueList = getUniqueList(descriptions)
    # get percent idents
    queurySeqs = blastResults['Query Seq'] 
    subSeq = blastResults['Subject Sequence'] 
    percentIdents = getPercentIdents(queurySeqs,subSeq)
    return topHit, uniqueList, percentIdents

def getUniqueList(descriptions):
    '''takes in list of all hit proteins, outputs a unique, cleaned list '''
    speciesList = []
    for i in range(len(descriptions)): 
        description = descriptions[i]
        # find all pdb matches for this protein
        sep = ">"
        ids = [pos for pos, char in enumerate(description) if char == sep]
        print(ids)
        if len(ids)>0: # if more than one pdb matches are found
            if len(ids) == 1: # if only 1 pdb match is found
                substring = description[ids[0]+1:len(description)]
                brackStart = description.find("[")
                brackEnd = description.find("]")
                species = description[brackStart:brackEnd]
                speciesList.append(substring)
            else: # if more than one pdb matches are found
                # get all [ NOT followed by a (
                subString = description[ids[0]+1:len(description)] 
                brackStart = [pos for pos, char in enumerate(description) if char == "["]
                print(brackStart)
                brackEnd = [pos for pos, char in enumerate(description) if char == "]"]
                ### FIX ME ###
        else: # if the protein hit has no pdb matches
            brackStart = description.find("[")
            brackEnd = description.find("]")
            species = description[brackStart:brackEnd]
            speciesList.append(substring)
    speciesArray = np.array(speciesList)
    uniqueList = np.unique(speciesArray)
    uniqueList = list(uniqueList)
    return uniqueList

def getPercentIdents(queurySeqs,subSeq):
    percentIdentities = []
    for i in range(queurySeqs.size): 
        matches = sum (a==b for a, b in zip(queurySeqs[i], subSeq[i]))
        percentIdentity = matches/len(queurySeqs[i]) *100
        percentIdentities.append(percentIdentity)
    return percentIdentities

# histogram of percent identities of the first HSPs from all the hits
if __name__ == "__main__":
    filepath = "/Users/zsj24/GitHub/NGS-Methods/query_blastpdb.xml"
    blastResults = parseXML(filepath)
    [topHit, uniqueList, PercIdentHist] = getInfo(blastResults)
    print(uniqueList)
    plt.hist(PercIdentHist)
    plt.title("Percent identities of the first HSPs from all the hits")
    plt.xlabel("Percent Identity")
    plt.ylabel("Frequency")
    plt.show()

