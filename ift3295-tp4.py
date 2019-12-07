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


""" CODE QUE NICO A AJOUTÉ """

seqs = ["INDNPPVFRGREQIIFIPESSRFPIEDADIGANALLTYTLELSKSLWLELRKYLDREETPELHLLLTATDGGKPELQGTVELLITVLDVNDNAPLFDQAVYRVHLL"
        "ETLVTTLNASDADEGVNGEVVFSFDSGISRDFKVDSSSGEIRLIDKLKVLDVNDNAPELAVTSLYLPIREDAPLSTVIALITVSRDSGANGQVTCSLMPHVPFKLVS",
        "MNDDGKVNASSEGYFILVGFSNWPHLEVVIFVVVLIFYLMTLIGNLFIIILSYLDSHLHTPMYFFLSYTTSSIPQYAGCMIVLLVVMSYDRYAAVCRPLHYTVLMH"
        "PRFCHLLAVASWVSGFTNSSSFTFWVPLCGHRQVDHFFCEVPALLRLSCVDTHVNELTLMITSSIFVIPLILILTAIRAVLRMQSTTGLFGTCGAHLMAVSLFFIVT",
        "INDNPPIFKGSEQRIFIPENSRFPLEDADIGANSLLTYTLELSKSLSLELRKSLDREETPELQLLLTATDGGKPELEGAVRLQITVLDVNDNAPVFDQAVYRAQLT"
        "ESLVTTLNATDADEGVNGEVVFSFGNDVSPDFKVDSISGEIRVIGDLKVLDVNDNVPELMITSLSLPIKEDAPLNTVVALIKVSIDSGVNGQVTCSLSPHLPFKLVS",
        "INDNPPVFKGAEQRIFIPENSRFPLEDADIGANSLLTYTLELSKSLSLELRKSLDREETPELQLLLTATDGGKPELEGTVRLQITVLDVNDNAPLFDQAIYRAQLV"
        "ESLVTTLNATDADEGVNGEVVFSFGNDVSLDFNVDSLSGEIRVIGDLKVLDINDNAPELSITSLSLPIKEDTPLNTIIALIKVSIDSGVNGQVTCSLTPHVPFKLVS",
        "INDNPPVFRGREQIIFIPESSRFPIEDADIGANALLTYTLELSKSLWLELRKSLDREETPELHLLLTATDGGKPELQGTVELLITVLDVNDNAPLFDQAVYRVLLL"
        "ETLMTTLNASDADEGVNGEVVFSFDSGISRDFKVDSSSGEIRLIDKLKVLDVNDNAPELAVTSLYLPIREDAPLSTVIALITVSRDSGANGQVTCSLMPHVPFKLVS"]

weight = "BLOSUM62.txt"


# Calculates distance matrix of sequences
def getDistanceMatrix(seqs, weight):
    weightDict = {}
    file = open(weight, "r")
    line = file.readline()
    weightHead = line.split()[1:]
    while line:
        line = file.readline()
        lineArr = line.split()
        if len(lineArr) > 0:
            weightDict[lineArr[0]] = lineArr[1:]
    distMat = [[float('inf') for i in range(len(seqs))] for j in range(len(seqs))]
    for i in range(len(distMat)-1):
        for j in range(i+1, len(distMat[0])):
            distMat[i][j] = calculateD(seqs[i], seqs[j], weightHead, weightDict)
            distMat[j][i] = distMat[i][j]
    return distMat


# calculates d(i,j)
def calculateD(seqi, seqj, header, dict):
    p, qi, qj = 0, 0, 0
    for k in range(len(seqi)):
        ai = seqi[k]
        aj = seqj[k]
        indexi = header.index(ai)
        indexj = header.index(aj)
        p += float(dict[ai][indexj])
        qi += float(dict[ai][indexi])
        qj += float(dict[aj][indexj])
    return 1 - (p / max(qi, qj))


# Neighbor-Joining algorithm
def neighborJoining(matrix, seqs, newNJNode):
    nicePrint(matrix, "\nItération : " + str(newNJNode))
    if len(seqs) == 2:
        #node = struct.Node(newNJNode)
        print("Créer noeud ayant pour enfant gauche : " + str(0))
        print(min(seqs[0], seqs[1]))
        #node.left = min(seqs[0], seqs[1])
        print("et pour enfant droit : " + str(1))
        print(max(seqs[0], seqs[1]))
        #node.right = max(seqs[0], seqs[1])
        print("\n --------  Fin de Neighbor-Joining  -------- ")
        return tree
    #totalDistances = getTotalDistances(matrix)
    mini, minj = findMinimum(matrix)
    #limbLengthi = 0.5 * (matrix[mini][minj] + (totalDistances[mini] - totalDistances[minj]) / (len(seqs)-2))
    #limbLengthj = matrix[mini][minj] - limbLengthi
    #node = struct.Node(newNJNode)
    print("Créer noeud ayant pour enfant gauche : " + str(min(mini, minj)))
    print(seqs[min(mini, minj)])
    #node.left = seqs[min(mini, minj)]
    #leaves[mini].data = limbLengthi
    print("et pour enfant droit : " + str(max(mini, minj)))
    print(seqs[max(mini, minj)])
    #node.right = (seqs[max(mini, minj)])
    #leaves[minj].data = limbLengthj
    #tree.append(node)
    newMatrix, newSeqs, newNJNode = createNewNJ(matrix, seqs, mini, minj, newNJNode)
    return neighborJoining(newMatrix, newSeqs, newNJNode)


# Returns the sum of each line of an array
def getTotalDistances(matrix):
    total = []
    for s in matrix:
        total.append(sum(s))
    return total


# Creates the matrix used for Neighbor-Joining algorithm
def getNJMatrix(matrix, totDist):
    njMat = [[0 for i in range(len(matrix))] for j in range(len(matrix[0]))]
    for i in range(len(njMat)-1):
        for j in range(i+1, len(njMat[0])):
            njMat[i][j] = (len(matrix)-2) * matrix[i][j] - totDist[i] - totDist[j]
            njMat[j][i] = njMat[i][j]
    return njMat


# Returns the minimal value of matrix and its indexes
def findMinimum(matrix):
    mins = []
    for l in matrix:
        mins.append(min(l))
    minimum = min(mins)
    mini = mins.index(minimum)
    minj = matrix[mini].index(minimum)
    return mini, minj


def createNewNJ(matrix, seqs, indexi, indexj, newNJNode):
    newMatrix = []
    lastLine = []
    for i in range(len(matrix)):
        if i != indexi and i != indexj:
            line = []
            for j in range(len(matrix[0])):
                if j != indexi and j != indexj:
                    line.append(matrix[i][j])
            newCell = (matrix[indexi][i] + matrix[indexj][i] - matrix[indexi][indexj]) * 0.5
            line.append(newCell)
            lastLine.append(newCell)
            newMatrix.append(line)
    lastLine.append(float('inf'))
    newMatrix.append(lastLine)
    seqs.pop(max(indexi, indexj))
    seqs.pop(min(indexi, indexj))
    seqs.append("NewNJNode_" + str(newNJNode))
    newNJNode += 1
    return newMatrix, seqs, newNJNode


# Creates an array containing all the original leaves
def createLeaves(array):
    leaves = []
    for i in range(len(array)):
        leaves.append(struct.Node(0))
    return leaves


def nicePrint(matrix, title):
    print(title)
    newMatrix = []
    for i in range(len(matrix)):
        newLine = []
        for j in range(len(matrix[0])):
            newLine.append(float("{0:.3f}".format(matrix[i][j])))
        newMatrix.append(newLine)
    for line in newMatrix:
        print(line)
    print(" ")
    return


# q2.1
distMat = getDistanceMatrix(seqs, weight)
nicePrint(distMat, "Matrice de distance")
# q2.2
leaves = createLeaves(seqs)
tree = leaves
print("\n --------  Début de Neighbor-Joining  -------- ")
neighborJoining(distMat, seqs, 1)
# Tests

""" FIN DU CODE QUE NICO A AJOUTÉ """

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