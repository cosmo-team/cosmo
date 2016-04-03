import sys
import Bio.SeqIO

with open(sys.argv[1], "w") as out:

        for seq in Bio.SeqIO.parse(sys.stdin, "fasta"):
            #Bio.SeqIO.write(seq, out, "fasta")
            out.write(">" + seq.id + "\n")
            out.write(str(seq.seq) + "\n")
        
