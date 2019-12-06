'''
Binary tree visualization tool : https://www.cs.usfca.edu/~galles/visualization/BST.html
Newick to tree tool : http://etetoolkit.org/treeview/#
'''

import struct
import sys
import newick_binaryTree


""" Functions """

# Reads newick files and filters out newlines. Returns a list
def readNewickFile(filename):
    return list(filter(None, open(filename).read().splitlines()))


# Reads fasta files and returns a dictionary
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


# Returns a list of the values of a dictionary
def getValues(dict):
    values = []
    for i in dict:
        values.append(dict[i])
    return values


# Returns a list of the names of a dictionary
def getNames(dict):
    return list(dict.keys())


# Makes a dictionary from a list of keys and a list of respective values
def makeDictionary(keys, values):
    newDict = {}
    for element in range(len(keys)):
        newDict.update({keys[element]: values[element]})
    return newDict


# Removes gaps from a list of multiple sequences
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
        

# Converts a newick string (keys -> values) from dictionary
def proteinKeyToValue_NewickString(nw, dictionary):
    for i in range(len(dictionary)):
        nw = nw.replace(list(dictionary.keys())[i],list(dictionary.values())[i])
    return nw


# Converts a list of newick strings (keys -> values) from dictionary
def proteinKeyToValue_NewickList(nwList, dictionary):
    for i in range(len(nwList)):
        nwList[i] = proteinKeyToValue_NewickString(nwList[i], dictionary)
    return nwList



""" Gestion des proteines: Creation de dictionaire sans gaps """
proteines = readFastaFile("proteines.fasta")
# print(proteines)
proteinesValues = getValues(proteines)
# print(proteinesValues)
proteinesNames = getNames(proteines)
# print(proteinesNames, proteinesValues)
noGapProteinesValues = removeGaps(proteinesValues)
# print(noGapProteinesValues)
proteinesDictionary = makeDictionary(proteinesNames, noGapProteinesValues)
# print(proteinesDictionary)



""" Gestion des arbres: Remplacement (key -> value) dans les arbres newick """
arbres = readNewickFile("arbres.nw")
# print(arbres)
arbresProteines = proteinKeyToValue_NewickList(arbres,proteinesDictionary)
# print(arbresProteines)



""" Tests """
# tab = ["ab-vd","d-uvg","5-555","-9999","00--0","11111"]
# noGapsTab = removeGaps(tab)
# print(noGapsTab)