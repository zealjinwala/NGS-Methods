
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
from matplotlib import pyplot as plt
# from typing import List

def parseXML(filepath):
    """takes in the XML(BLAST result), outputs list in csv format
    Args:
        filepath (str): path to the XML file containing blast results
    Returns:
        blastresults (df): extracted results      
    """
    from Bio.Blast import NCBIXML
    resultHandle = open(filepath)
    blastRecords = NCBIXML.parse(resultHandle)
    blastResults = pd.DataFrame(columns = [  'Query',
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
            blastResults = blastResults.append({'Query':query,
                'Name':name , 'Length': length,'Score': score, 
                'Expect':expect,'Query Seq':querySeq, 
                'Subject Sequence':subSeq}, ignore_index = True)
    resultHandle.close()
    return blastResults

def getInfo(blastResults):
    """ takes in list of results, outputs results 
    Args:
        blastResults (pandas dataframe): blast results dataframe
    Returns:
        topHit (str): top protein hit
        uniqueList (pandas series): Unique list of species matched
        percentIdents (list): percent identities of the first HSPs from all the hits
    """
    hitName = str(blastResults['Name'].iloc[0])
    sep = ">"
    ids = [pos for pos, char in enumerate(hitName) if char == sep]
    topHit = hitName[ids[0]+1:ids[1]]
    descriptions = blastResults.Name
    uniqueList = getUniqueList(descriptions)
    queurySeqs = blastResults['Query Seq'] 
    subSeq = blastResults['Subject Sequence'] 
    percentIdents = getPercentIdents(queurySeqs,subSeq)
    return topHit, uniqueList, percentIdents

def getUniqueList(descriptions):
    """takes in list of all hit proteins, outputs a unique, cleaned list
    Args:
        descriptions (pandas series): Descriptions of all pdb hits
    Returns:
        uniqueList (pandas series): Unique list of species matched
    """
    speciesList = []
    for i in range(len(descriptions)): 
        description = descriptions[i]
        # find all pdb matches for this protein
        sep = ">"
        ids = [pos for pos, char in enumerate(description) if char == sep]
        if len(ids)>0: # if more than one pdb matches are found
            if len(ids) == 1: # if only 1 pdb match is found
                substring = description[ids[0]+1:len(description)]
                brackStart = substring.find("[")
                brackEnd = substring.find("]")
                species = description[brackStart:brackEnd]
                speciesList.append(species)
            else: # if more than one pdb matches are found
                brackStart = [pos for pos, char in enumerate(description) if char == "["]
                brackEnd = [pos for pos, char in enumerate(description) if char == "]"]
                for i in range(len(brackStart)):
                    if description[brackStart[i]+1] == "(":
                        species = description[brackStart[i+1]:brackEnd[i+1]]
                        speciesList.append(species)
        else: # if the protein hit has no pdb matches
            brackStart = description.find("[")
            brackEnd = description.find("]")
            species = description[brackStart:brackEnd]
            speciesList.append(species)
    speciesList = pd.DataFrame(speciesList, columns = ['Species'],dtype = str)
    uniqueList = speciesList['Species'].unique()
    return uniqueList

def getPercentIdents(queurySeqs,subSeq):
    """takes in list of all hit proteins, outputs a unique, cleaned list
    Args:
        queurySeqs (pandas series): query sequnece of the protien hit
        subSeq (pandas series): reference/subject sequnece of the protien hit
    Returns:
        percentIdentities (pandas series): percent identities of the first HSPs from all the hits
    """
    percentIdentities = []
    for i in range(queurySeqs.size): 
        matches = sum (a==b for a, b in zip(queurySeqs[i], subSeq[i]))
        percentIdentity = matches/len(queurySeqs[i]) *100
        percentIdentities.append(percentIdentity)
    return percentIdentities

if __name__ == "__main__":
    filepath = "/Users/zsj24/GitHub/NGS-Methods/query_blastpdb.xml"
    blastResults = parseXML(filepath)
    [topHit, uniqueList, PercIdentHist] = getInfo(blastResults)
    print("Top protein Hit: \n", topHit)
    print("Unique list of species: \n", uniqueList)
    plt.hist(PercIdentHist)
    plt.title("Percent identities of the first HSPs from all the hits")
    plt.xlabel("Percent Identity")
    plt.ylabel("Frequency")
    plt.show()



