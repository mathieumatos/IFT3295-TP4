import struct

# Returns a binary tree from a newick string
def treeBuilder(nw):
    if(nw[len(nw)-1]==";"):
        nw = nw[:-1]

    if nw[0] == "(":
        node = struct.Node("1")
        nextLeft = ""
        nextRight = ""
        nbParenthese = 0
        rightIndex = 0
        for s in range(1,len(nw)):
            if nw[s] == "(":
                nbParenthese = nbParenthese + 1
            elif nw[s] == ")":
                nbParenthese = nbParenthese - 1
            
            if (nw[s] == ",") and (nbParenthese == 0):
                if nextLeft[0] == "(":
                    node.left = treeBuilder(nextLeft)
                else:
                    node.left = struct.Node(nextLeft)
                rightIndex = s+2
                break
 
            nextLeft = nextLeft + nw[s]

        nextRight = nw[rightIndex:-1]
        if nextRight[0] == "(":
            node.right = treeBuilder(nextRight)
        else:
            node.right = struct.Node(nextRight)

    return node


# Returns all the leafs of a tree
def getLeafs(tree):
    leafs = []
    for node in tree.postOrderTraversal(tree):
        if node.data != "1":
            leafs.append(node.data)
        else:
            continue
    return set(leafs)


# Returns if node is leaf or not
def isLeaf(node):
    return node.left == None and node.right == None


# Returns bipartitions from node in format [[all leafs from node],[all leafs from node - leafs from node i]]
def getBipartitions(node):

    allLeafs = getLeafs(node)
    # print("ALL LEAFS : "+str(allLeafs))
    bipartitions = []

    for nd in node.postOrderTraversal(node):
        # print("STARTING FROM NODE "+str(nd.data))

        if(isLeaf(nd)):
            # print("IS LEAF, NEXT NODE")
            continue
  
        nodeLeafs = getLeafs(nd)
        # print("LEAFS FROM NODE ND : "+str(nodeLeafs))
        notLeafs = allLeafs-nodeLeafs
        # print("NOT LEAFS : "+str(notLeafs))

        if(not len(notLeafs) <= 1):
            bipartitions.append([sorted(allLeafs),sorted(notLeafs)])
            # print("ADDED BIPARTITION "+str([sorted(allLeafs),sorted(notLeafs)]))
    # print("ALL BIPARTITIONS : "+str(bipartitions))
    # print("BIPARTITIONS LEN : "+str(len(bipartitions)))
    return bipartitions


# Returns RF score between two arrays of bipartitions
def robinson_foulds(bipartition1,bipartition2):
    rf = len(bipartition1) + len(bipartition2)
    for partition1 in bipartition1:
        for partition2 in bipartition2:

            if(partition1 == partition2):
                rf = rf - 2

    return rf


# Tests
# str1 = "(((PCDHA1_Humain, PCDHA1_Rat), (PCDHA1_Souris, PCDHA1_Bonobo)), OR2J3_Humain);"
# str2 = "(((OR2J3_Humain, PCDHA1_Rat), (PCDHA1_Souris, PCDHA1_Bonobo)), PCDHA1_Humain);"

# tree = treeBuilder(str1)
# tree2 = treeBuilder(str2)

# print("RF : "+str(robinson_foulds(getBipartitions(tree),getBipartitions(tree2))))

