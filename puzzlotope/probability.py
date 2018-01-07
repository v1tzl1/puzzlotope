import copy

class TreeNode:
	parent=None
	children=[]
	label=''
 	
	def __init__(self, parent=None, label='', prob=1.0):
		self.parent = parent
		self.children = []
		self.label = label
		self.prob = prob if self.isRoot() else parent.prob*prob
	
	def isRoot(self):
		return (self.parent is None)
	
	def getRoot(self):
		return self if self.isRoot() else self.parent.getRoot()
	
	def getNodes(self):
		stack = [self.getRoot()]
		while stack:
			node = stack.pop()
			yield node
			for child in node.children:
				stack.append(child)
	
	def getLeafs(self):
		return [n for n in list(self.getNodes()) if n.isLeaf()]
 		 		 		
	def isLeaf(self):
		return len(self.children)==0
 	
	""" return array with labels along the tree. root label is omitted """
	def getBranchString(self):
		ret = []
		if self.isRoot():
			pass
		elif self.parent.isRoot():
			ret.append(self.label)
		else:
			ret += self.parent.getBranchString()
			ret.append(self.label)
		return ret
 	
	def delete(self):
		for i in range(len(self.parent.children)):
			if self.parent.children[i].label == self.label:
				self.parent.children.pop(i)
				break
		if len(self.parent.children)==0:
			self.parent.delete()
	
	def absorb(self, other):
		self.prob += other.prob
		other.delete()
 	
	""" add new branches to current node """
	def addBranch(self, probs, labels):
		for i in range(len(probs)):
			self.children.append(TreeNode(self, labels[i], probs[i]))
 	
	""" add new branches to all leaves """
	def addLevel(self, probs, labels):
		#print('add %s to tree' % str(labels))
		root=self.getRoot()
		#root.printTree()
		leafs=root.getLeafs()
		for l in leafs:
			#print('  addng branch to leaf %s' % l)
			l.addBranch(probs, labels)
		#print('new tree')
		#root.printTree()
		root.simplify()
		#print('simplified')
		#root.printTree()
	
	def cutoff(self, cutoff):
		for l in self.getLeafs():
			if l.prob<=cutoff:
				l.delete()
	
	def simplify(self):
		leafs=copy.copy(self.getLeafs())
		while leafs:
			l = leafs.pop()
			bs=sorted(l.getBranchString())
			flag=True
			while flag:
				flag=False
				for i in range(len(leafs)):
					if sorted(leafs[i].getBranchString())==bs:
						l.absorb(leafs[i])
						leafs.pop(i)
						flag=True
						break
						
	def printTree(self):
		print('Tree with leafs: ' + str([l for l in self.getLeafs()]))
		
	def __repr__(self):
		#print('bs: %s' % str(self.getBranchString()))
		if self.isRoot():
			return 'Root node'
		else:
			#print('repr for node. label=%s bs=%s' % (self.label, self.getBranchString()))
			return 'Node %s (%f%%)' % ('->'.join(self.getBranchString()), self.prob*100.0)

if __name__=='__main__':
	tree = TreeNode()
	tree.printTree()
	tree.addLevel([0.5, 0.5], ['a1','a2'])
	tree.printTree()
	tree.addLevel([0.5, 0.5], ['a1','a2'])
	tree.printTree()
	tree.addLevel([0.9, 0.1],['b1','b2'])
	tree.printTree()