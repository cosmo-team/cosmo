# Download the zip files
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
/usr/bin/python3 fetch_ecoli.py

# Decompress them and join the lines of each sequence since KMC2 assumes reads, not contigs, and expects them to be one line only
# Don't forget to change KMC/kmer_counter/splitter.h line to allow longer lines; by default it will only accept sequences up to 2^14 bp
# I used this for e.coli:
# template <bool QUAKE_MODE> uint32 CSplitter<QUAKE_MODE>::MAX_LINE_SIZE = 1 << 23; //14
ls -1 --color=no *.fna.gz |xargs -l -i echo "gzip -dc {} |/usr/bin/python3 join_fasta_lines_stream.py {}.flat" >flatten.sh
source flatten.sh

# KMC2 needs a temp dir
mkdir kmc_temp

# count kmers
ls -1 --color=no *.flat |xargs -l -i echo "~/git/KMC/bin/kmc -ci0 -fa -k31 -cs300 {} {}.kmc kmc_temp" >kmercount.sh
source kmercount.sh

# sort them
ls -1 --color=no *.flat |xargs -l -i echo "~/git/KMC/bin/kmc_tools sort {}.kmc {}.kmc.sorted " >kmercountsort.sh
source kmercountsort.sh

# make a list of these files which will be the cosmo input
ls -1 --color=no *.flat |xargs -l -i echo "{}.kmc.sorted" >kmc2_list


# build the BOSS
numactl --interleave=all /bin/time -v ~/git/cosmo/cosmo-pack -k kmc2_list

# build the RRR
numactl --interleave=all /bin/time -v /s/chopin/l/grad/muggli/git/cosmo/pack-color kmc2_list.colors 3765

# And finally, find bubbles!
/bin/time -v  ~/git/cosmo/cosmo-color kmc2_list.packed kmc2_list.colors.rrr >kmc2_list.bubbles


