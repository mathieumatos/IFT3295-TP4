import struct

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


# TESTS
str = "((PCDHA1_Humain, OR2J3_Humain), ((PCDHA1_Rat, PCDHA1_Souris), PCDHA1_Bonobo));"

tree = treeBuilder(str)
tree.printStruct()
print(tree.postOrderTraversal(tree))
print("TREE MANUAL SEARCH")
print("tree.data = " +(tree.data))
print("tree.left.data = " +(tree.left.data))
print("tree.right.data = " +(tree.right.data))
print("tree.left.left.data = " +(tree.left.left.data))
print("tree.left.right.data = " +(tree.left.right.data))
print("tree.right.left.data = " +(tree.right.left.data))
print("tree.right.right.data = " +(tree.right.right.data))
print("tree.right.left.left.data = "+(tree.right.left.left.data))
print("tree.right.left.right.data = "+(tree.right.left.right.data))