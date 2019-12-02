'''
Class code taken from https://www.tutorialspoint.com/python_data_structure/python_tree_traversal_algorithms.htm
Used for post-order traversal, each child (left/right) is visited before its respective parent.
'''

class Node:

    def __init__(self, data):
        self.left = None
        self.right = None
        self.data = data

    def insert(self, data):
        if self.data:
            if data < self.data:
                if self.left is None:
                    self.left = Node(data)
                else:
                    self.left.insert(data)
            elif data > self.data:
                if self.right is None:
                    self.right = Node(data)
                else:
                    self.right.insert(data)
        else:
            self.data = data

    # prints in ascending order
    def printStruct(self):
        if self.left:
            self.left.printStruct()
        print(self.data)
        if self.right:
            self.right.printStruct()

    # returns array corresponding to the post-order traversal to the root
    def postOrderTraversal(self, root):
        res = []
        if root:
            res = self.postOrderTraversal(root.left)
            res = res + self.postOrderTraversal(root.right)
            res.append(root.data)
        return res

# Tests
# root = Node(1)
# root.left = Node(1)
# root.right = Node(1)
# root.left.left = Node("Left Left")
# root.left.right = Node("Left Right")


# root.printStruct()
# print(root.postOrderTraversal(root))