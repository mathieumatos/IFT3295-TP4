'''
Binary tree visualization tool : https://www.cs.usfca.edu/~galles/visualization/BST.html
'''

import struct
import sys

root = struct.Node(100)
root.insert(13)
root.insert(33)
root.insert(102)
root.printStruct()
print(root.postOrderTraversal(root))