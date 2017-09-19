
import os
import re
import sys
import optparse
import tempfile
import subprocess
import shutil
p = optparse.OptionParser()

p.add_option("--input1", action="store", dest="input1", help="sequence data in fasta format")
# p.set_defaults(query_order="2")

opts,args = p.parse_args()
#query = opts.query



# install directory
instdir = os.path.dirname(os.path.realpath(sys.argv[0]))

tempdir = tempfile.mkdtemp()
print("# Storing temporary files in directory", tempdir)

cmd = instdir + "/3rd_party_src/bin/kmc -ci0 -k10 -fm " + opts.input1 + " " + opts.intput1 + ".kmc " + tempdir


shutil.rmtree(tempdir)
