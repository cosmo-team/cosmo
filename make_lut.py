import numpy as np

src = """
#include "lut.h"

const unsigned char reverse_8_lut[256] = {
  %s
};

const unsigned char revcomp_8_lut[256] = {
  %s
};
"""

# Generate tables
reverse_8 = lambda x: ((x & 0xC0) >> 6) | (x << 6) | ((x & 0x30) >> 2) | ((x & 0x0C) << 2)
complement = lambda x: ~x

all_bytes = np.array(range(256), np.uint8)
rev_bytes = reverse_8(all_bytes)
rev_comps = complement(rev_bytes)

print src % (",\n ".join(map(hex, rev_bytes.tolist())),
             ",\n  ".join(map(hex,rev_comps.tolist())))
