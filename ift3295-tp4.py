'''
Binary tree visualization tool : https://www.cs.usfca.edu/~galles/visualization/BST.html
Newick to tree tool : http://etetoolkit.org/treeview/#
'''

import struct
import sys
import newick_binaryTree


# Reads newick files and filters out newlines. Returns a list
def readNewickFile(filename):
    return list(filter(None, open(filename).read().splitlines()))


# Reads fasta files and returns a dictionnary
def readFastaFile(filename):
    fastaDict = {}
    file = open(filename).read().splitlines()
    names = []
    prots = []
    nextSeq = ""
    for line in range(len(file)):
        if(file[line][0] == ">"):
            names.append(file[line][1:])
            if len(nextSeq) > 0:
                prots.append(nextSeq)
                nextSeq = ""
        else:
            nextSeq = nextSeq + file[line]
    if len(nextSeq) > 0:
                prots.append(nextSeq)

    for name in range(len(names)):
        fastaDict.update({names[name]:prots[name]})

    return fastaDict


def getValues(dict):
    values = []
    for i in dict:
        values.append(dict[i])
    return values



def removeGaps(seqs):
    newSeqs = seqs.copy()
    gapIndexTab = []
    for i in range(len(seqs)):
        indexTab = []
        for j in range(len(seqs[0])):
            if seqs[i][j] == "-":
                indexTab.append(j)
        gapIndexTab = gapIndexTab + indexTab
    
    gapIndexTab = list(dict.fromkeys(gapIndexTab))
    gapIndexTab.sort(reverse=True)

    for i in range(len(gapIndexTab)):
        for seq in range(len(seqs)):
            newSeqs[seq] = newSeqs[seq][:gapIndexTab[i]] + newSeqs[seq][gapIndexTab[i]+1:]

    return newSeqs
        


arbres = readNewickFile("arbres.nw")
proteines = readFastaFile("proteines.fasta")
# print(proteines)
proteinesValues = getValues(proteines)
# print(proteinesValues)
noGapProteins = removeGaps(proteinesValues)
print(noGapProteins)


# mini tests
# tab = ["ab-vd","d-uvg","5-555","-9999","00--0","11111"]
# noGapsTab = removeGaps(tab)
# print(noGapsTab)