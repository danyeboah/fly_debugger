# Daniel Yeboah-Kordieh
# This program just reads in a newick file and displays it

# import necessary modules

# To run...use python display_tree.py <name of newick file>
from Bio import Phylo
from Bio import motifs

# import sys
import sys

# read in name of input newick file and motif file
inputnewick_file = sys.argv[1]


# read in tree
tree = Phylo.read(inputnewick_file, "newick")

# convert tree to phyloxml format...this allows manipulation of the tree in terms of visualization
tree_xml = tree.as_phyloxml()
print (tree.total_branch_length())

	
Phylo.draw(tree_xml, branch_labels=lambda c: c.branch_length)

