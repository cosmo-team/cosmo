ls -1 --color=no *.fna |xargs -l -i echo "cat {} |python3 ~/git/cosmo/experiments/join_fasta_lines_stream.py {}.flat" >flatten.sh
source flatten.sh 
ls -1 --color=no *.flat |xargs -l -i echo "~/git/KMC/bin/kmc -ci0 -fa -k31 -cs300 {} {}.kmc kmc_temp" >kmercount.sh
source kmercount.sh
ls -1 --color=no *.flat |xargs -l -i echo "~/git/KMC/bin/kmc_tools sort {}.kmc {}.kmc.sorted " >kmercountsort.sh
source kmercountsort.sh

numactl --interleave=all /bin/time -v ~/git/cosmo/cosmo-pack -k ecoli6_kmc2_list
numactl --interleave=all /bin/time -v ~/git/cosmo/pack-color ecoli6_kmc2_list.colors 6
/bin/time -v  ~/git/cosmo/cosmo-color ecoli6_kmc2_list.packed ecoli6_kmc2_list.colors.rrr
