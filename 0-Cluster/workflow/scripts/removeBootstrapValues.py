""""
RAxML output tree includes bootstrap values in square brackets, which cannot be read
in Newick format. Remove square brackets and all values contained within them, and write
tree to new file.
"""

import sys
import re

ORIGINAL_TREE_F = sys.argv[1]
OUTPUT_TREE_F = sys.argv[2]

with open(ORIGINAL_TREE_F, 'r') as f:
    tree = f.read()

    edited_tree = re.sub("[\[].*?[\]]","",tree)

with open(OUTPUT_TREE_F, 'w') as f:
    f.write(edited_tree)
