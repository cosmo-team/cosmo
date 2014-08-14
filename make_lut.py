import numpy as np

src = """#pragma once
#ifndef LUT_HPP
#define LUT_HPP

#define reverse_8(x) (reverse_8_lut[(x)])
#define revcomp_8(x) (revcomp_8_lut[(x)])

static const uint8_t reverse_8_lut[256] = {
  %s
};

static const uint8_t revcomp_8_lut[256] = {
  %s
};

#endif
"""

# Generate tables
reverse_8 = lambda x: ((x & 0xC0) >> 6) | (x << 6) | ((x & 0x30) >> 2) | ((x & 0x0C) << 2)
complement = lambda x: ~x

all_bytes = np.array(range(256), np.uint8)
rev_bytes = reverse_8(all_bytes)
rev_comps = complement(rev_bytes)

print src % (",\n ".join(map(hex, rev_bytes.tolist())),
             ",\n  ".join(map(hex,rev_comps.tolist())))
