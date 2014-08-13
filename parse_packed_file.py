import numpy as np
import sys

W = 64 # Block Width
E = 5  # Element Width

# Blocks have r-bit fixed width elements, packed from 0th bit
def get_element_from_block(block, i, r=E):
    return (int(block) & (31 << (r * ((W/r)-i-1)))) >> (r * ((W/r)-i-1))

# Find which block the given index is in, then extract the element
def get_element(blocks, i, r=E):
    block_idx = i/(W/r)
    return get_element_from_block(blocks[block_idx], i%(W/r))

get_sym   = lambda x: "$acgt"[x >> 2]
get_flags = lambda x: (x & 2 >> 1, x & 1)

if len(sys.argv) != 2:
    print "Please provide an input file!"
    exit()
filename = sys.argv[1]

with open(filename) as f:
    data = np.fromfile(f, dtype=np.uint64)
# Last 6 uint64s are the 5 ($,a,c,g, and t) cumulative counts for the 2nd last column, and k respectively
counts = data[-6:-1]
blocks = data[:-6]
k = data[-1]
total = counts[-1]

last = 1
for i in xrange(total):
    edge = get_element(blocks, i)
    sym = get_sym(edge)
    first_start_node, first_end_node = get_flags(edge)
    print first_start_node, sym, first_end_node

print " ".join(map(str,counts))
print k
