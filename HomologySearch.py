
# By Zeal Jinwala 
# Date: March 19, 2022

# Description:
# This program take an XML file (BLAST search results) as input and answers the following questions:
#   1) Prints the name of the top protein hit.
#   2) Prints a unique list of species names of all the hits.
#   3) Finds the top scoring hit with a mouse protein. Prints the local sequence alignment of the query with this mouse protein.
#   4) Shows a histogram of percent identities of the first HSPs from all the hits.

# Query Protein: <br>
# SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIG
# HSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVS
# FCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKY
# NYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQ

# import argparse
import pandas as pd

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
    # top hit
    hitNames = str(blastResults.iloc[0,1])
    sep = ">"
    ids = [pos for pos, char in enumerate(hitNames) if char == sep]
    topHit = hitNames[ids[0]:ids[1]]
    # unique list of proteins
    descriptions = blastResults.Name
    uniqueList = getUniqueList(descriptions)
    # percentIdents = getPercentIdents()
    return topHit, uniqueList

def getUniqueList():

# def getPercentIdents():



    # top scoring hit with a mouse protein
    # histogram of percent identities of the first HSPs from all the hits
    return topHit, uniqueList, PercIdentHist

if __name__ == "__main__":
    filepath = "/Users/zsj24/GitHub/NGS-Methods/query_blastpdb.xml"
    blastResults = parseXML(filepath)
    print(blastResults.head(5))
    print((blastResults.Alignment[0]))
    [topHit, uniqueList, PercIdentHist] = getInfo(blastResults)

    # parser = argparse.ArgumentParser()

    # parser.add_argument(
    # "--xmlFileName",
    # type=str,
    # required=True,
    # help="Name of the xml file (BLAST result)",
    # )

    # args = parser.parse_args()
    # blastresults =  parseXML(
    #     args.xmlFileName
    #     )
