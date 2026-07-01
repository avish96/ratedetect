## Python script to remove in-line nodes from phylopomp transmission trees
## This results in a phylogenetic tree
##

####################################################################

## Load modules
import os
import sys
from ete3 import Tree # Install ete3 if not in system

import importlib.metadata
assert importlib.metadata.version("ete3") >= "3.1.2"

## Function to recursively remove in-line nodes in a binary tree
def trim_inline(node):
    # Only consider non-leaf nodes
    if not node.is_leaf():
        # Node should be removed if it does not start with 'm_'
        if not node.name.startswith(('m_')):
            # Add branch lengths of node and its child
            if len(node.children) == 1:  # For binary tree, 1 or 2 children
                child = node.children[0]
                child.dist += node.dist
                # Remove node by detaching it from its parent and connecting it directly to its child
                parent = node.up
                node.detach()
                if parent: # Add child of detached node to upstream parent node
                    parent.add_child(child)
                else:  # If the current node is root
                    t.set_outgroup(child)

        # Recursively call all child nodes
        for child in node.children:
            trim_inline(child)


## Read files in directory
file_list = []
for x in os.listdir(sys.argv[1]):
    if x.endswith(".nwk"):
        file_list.append(x)


## Remove in-line nodes across files
for ind in range(0,len(file_list)):
    file = open(file_list[ind])
    newick_string = file.read()
    file.close()

    print(file_list[ind])

    # Load the tree
    t = Tree(newick_string, format=1)

    # Start the collapsing process from the root
    tlen_change = 1
    while tlen_change > 0:
        tlen = len(t.write(format=1))
        trim_inline(t)
        tlen1 = len(t.write(format=1))
        tlen_change = tlen - tlen1


    # Obtain the modified Newick string
    newick_string_modified = t.write(format=1)

    # Save the modified phylogeny as nwk file 
    file_list_split = str.split(file_list[ind],'.nwk')
    f = open(file_list_split[0] + '_skip.nwk',"a")
    f.write(newick_string_modified)
    f.close()

