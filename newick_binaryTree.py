import struct

str = "((PCDHA1_Humain, OR2J3_Humain), ((PCDHA1_Rat, PCDHA1_Souris), PCDHA1_Bonobo));"

# root = struct.Node(1)

def treeBuilder(nw):
    if nw[0] == "(":
        node = struct.Node(1)
        nextLeft = ""
        nextRight = ""
        parenthCounter = 0
        rightIndex = 0
        for s in range(1,len(nw)):
            if nw[s] == "(":
                parenthCounter = parenthCounter + 1
            elif nw[s] == ")":
                parenthCounter = parenthCounter - 1
            
            if (nw[s] == ",") and (parenthCounter == 0):
                node.left = treeBuilder(nextLeft)
                rightIndex = s+2
                break

            nextLeft = nextLeft + nw[s]

        nextRight = nw[rightIndex:-1]
        print("NEXT LEFT = "+nextLeft)
        print("NEXT RIGHT = "+nextRight)
            



treeBuilder(str)

# root.printStruct()