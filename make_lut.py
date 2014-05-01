import numpy as np

doc = """
#pragma once
#ifndef LUT_H
#define LUT_H

#define revcomp_8(x) (revcomp_8_lut[(x)])

unsigned char revcomp_8_lut[256] = {
  %s
};

#endif"""

reverse_8 = lambda x: ((x & 0xC0) >> 6) | (x << 6) | ((x & 0x30) >> 2) | ((x & 0x0C) << 2)
complement = lambda x: ~x

bytes     = np.array(range(256), np.uint8)
rev_comps = complement(reverse_8(bytes))

print doc % (",\n  ".join(map(hex,rev_comps.tolist())))
