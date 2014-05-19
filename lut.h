#ifndef LUT_H
#define LUT_H

#define reverse_8(x) (reverse_8_lut[(x)])
#define revcomp_8(x) (revcomp_8_lut[(x)])

const unsigned char reverse_8_lut[256];

const unsigned char revcomp_8_lut[256];

#endif
