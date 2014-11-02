import numpy as np
from bisect import bisect_right
import argparse, os

# TODO: try ideas from http://scikit-learn.org/stable/auto_examples/cluster/plot_lena_compress.html

parser = argparse.ArgumentParser(description='Quantize LCS vector to compress better.')
parser.add_argument('filename')
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument('-q', dest='levels', metavar='N', nargs='*',
                    help='a quantization level.', type=int)
mode.add_argument('-t', dest='threshold', metavar='<threshold>', help='remove all below this minimum value.', type=int)

args = parser.parse_args()
filename, levels, threshold = args.filename, args.levels, args.threshold

# Add check for file
if not os.path.exists(filename):
  print "Error: file %s does not exist" % (filename)

lcs = np.fromfile(filename, dtype=np.uint8)
#max_val = lcs.max()

if threshold:
  print "Removing values below %s." % (threshold)
  quantize = np.vectorize(lambda x: x if x >= threshold else 0, otypes=[np.uint8])
if levels:
  if 0 not in levels: levels = [0]+sorted(levels)
  print "Quantizing into these buckets: %s" % (levels)
  quantize = np.vectorize(lambda x: levels[bisect_right(levels, x)-1], otypes=[np.uint8])

qzd = quantize(lcs)
qzd.tofile(filename+".qzd")
