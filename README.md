
                                     o8o  
                                     `"'  
     oooo    ooo  .oooo.   oooo d8b oooo  
      `88.  .8'  `P  )88b  `888""8P `888  
       `88..8'    .oP"888   888      888  
        `888'    d8(  888   888      888  
         `8'     `Y888""8o d888b    o888o
                                ver 0.7.0


# VARI

[**Version**][semver]: 0.7.0

Cosmo is a fast, low-memory DNA assembler that uses a [succinct de Bruijn graph][succ].

**VARI** is an extension to Cosmo and supports offline construction of succinct colored de Bruijn graphs.

If you use VARI in scholarly work, please cite:

_Muggli, M. D., Bowe, A., Noyes, N. R., Morley, P., Belk, K., Raymond, R., Gagie, T., Puglisi, S. J., and Boucher, C. (2017). Succinct colored de bruijn graphs. Bioinformatics. doi: 10.1093/bioinformatics/btx067_

```
@article{muggli2017vari,
  title={Succinct Colored de Bruijn Graphs},
  author={Muggli, Martin D and Bowe, Alexander and Noyes, Noelle R and Morley, Paul and Belk, Keith and Raymond, Robert and  Gagie, Travis and Puglisi, Simon J and Boucher, Christina},
  journal={Bioinformatics},
  year={2017},
  publish={Oxford Univ Press}
}
```

## Building notes

Five third party packages are required for VARI. All should be cloned within the 3rd_party_src directory.  Any 3rd party software may change in incompatible ways.  The revisions known to work for the published results are included.


1. KMC2 --  'git clone https://github.com/refresh-bio/KMC' (commit f090276855a3f7c0b14e9f3abc8c99d3213247b3)
2. sdsl-lite -- 'git clone https://github.com/cosmo-team/sdsl-lite.git' (commit 9fa981958a9d2ddade12d083548f2b09939514fb)
3. stxxl -- 'git clone https://github.com/stxxl/stxxl' (commit 5b9663e6b769748f3b3d3a9a779b4b89e24d7a27)
4. tclap -- 'git clone https://github.com/eile/tclap' (commit f41dcb5ce3d063c9fe95623193bba693338f3edb)
5. Boost 1.54* -- 'wget http://sourceforge.net/projects/boost/files/boost/1.54.0/boost_1_54_0.tar.bz2'

* VARI fails to compile with later versions of Boost.  For the time being, it is necessary to download and compile Boost 1.54 and update BOOST_PATH in the Makefile to reflect the installed directory.  See Issue [#7](/../../issues/7).

They should be configured and built following their own instructions and set to install their files in a 3rd_party_inst subdirectory which is a sibling of 3rd_party_src.  The following sequence of commands should build the required parts.  Compilation errors may or may not affect the functionality of VARI, as VARI doesn't use all functionality of 3rd party sources.   Please email me if you run into trouble. I'm intermitently working on streamlining the process. -MDM May 17, 2017

**Note**: Change "/home/martin_muggli/git/test/cosmo" to wherever your cosmo working tree ends up.

    # Fetch software and setup directories
    git clone https://github.com/cosmo-team/cosmo/
    cd cosmo/
    git checkout VARI
    mkdir 3rd_party_src
    mkdir -p 3rd_party_inst/boost
    cd 3rd_party_src
    git clone https://github.com/refresh-bio/KMC
    git clone https://github.com/cosmo-team/sdsl-lite.git
    git clone https://github.com/stxxl/stxxl
    git clone https://github.com/eile/tclap
    wget http://sourceforge.net/projects/boost/files/boost/1.54.0/boost_1_54_0.tar.bz2
    tar -xjf boost_1_54_0.tar.bz2

    # Build the five dependencies
    cd boost_1_54_0
    ./bootstrap.sh --prefix=../../3rd_party_inst/boost
    ./b2 install
    cd ..
    
    cd sdsl-lite/
    /usr/bin/time sh install.sh /home/martin_muggli/git/test/cosmo/3rd_party_inst
    cd ..

    cd stxxl
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/martin_muggli/git/test/cosmo/3rd_party_inst -DBUILD_STATIC_LIBS=ON
    make
    make install
    cd ../..

    cd KMC
    make
    cd ..

    cd tclap/
    autoreconf -fvi
    ./configure --prefix=/home/martin_muggli/git/test/cosmo/3rd_party_inst
    make
    make install
    cd ..
    
    # Build VARI
    make


## Usage

Below is an example of using the succinct colored de Bruijn graph for bubble calling.  

### External vs Internal memory
cosmo-build can use external memory for various arrays during construction.  See the STXXL website and documentation for configuring whether these arrays are in internal or external memory and where the external data should be placed on various drives/directories.  

### Color Matrix Compression
The color matrix can be compressed either with RRR or Elias-Fano encoding.  The current version uses Elias-Fano and will stream the uncompressed color matrix (.colors file emitted by cosmo-build) from disk and build the succinct version online.  The pack-color program in earlier commits uses RRR and reads the same .colors file and only takes the number of colors as an argument.  To use this flow, change the sdsl::sd_vector types to sdsl::rrr_vector.  

### Input files
cosmo-build can, in addition to the streaming KMC2 flow, accept a multi-colored cortex de Bruijn graph binary file (.ctx).  This flow was used during development before the implementation of the KMC2 flow. It requires that the host machine have sufficient RAM to store the non-succinct colored de Bruijn graph in memory.  

### k-mer size
k-mer size is determined by the -k parameter for KMC2.  The revision of KMC2 used has a bug where large-ish k values (in the neighborhood of 64) may result in corrupted k-mer values.  This will manifest as an abnormally large number of dummy nodes and likely running out of memory or disk space.  There is a workaround in io.hpp which currently must be manually uncommented.   See Issue [#6](/../../issues/6).

### Colored de Bruijn graph example:
```sh
# Grab 6 E. coli assemblies:
git clone https://github.com/cosmo-team/e_coli6
cd e_coli6

# Use KMC2 to k-mer count the FASTA (*.fna) files
$ mkdir kmc_temp
$ ls -1 --color=no *.fna |xargs -l -i  ~/kmc -ci0 -fm -k32 -cs300 {} {}_kmc kmc_temp
$ ls -1 --color=no *.fna |xargs -l -i  ~/kmc_tools sort {}_kmc {}_kmc_sorted_kmc.kmc
$ ls -1 --color=no *.fna |xargs -l -i echo "{}_kmc_sorted_kmc.kmc" >ecoli6_kmc2_list

# Build the succinct de Bruijn graph and permute uncompresed color matrix accordingly
# cosmo-build -d <KMC2_count_names> # KMC2_count_names list base names for k-mer counts produced by KMC2 (i.e. no .kmc_pre/.kmc_suf)
$ cosmo-build -d ecoli6_kmc2_list

# Make succinct color matrix
# pack-color -input <filename.colors>  <num colors> <total bits> <set bits> 
#     The sdsl-lite Elias Fano encoder must know ahead of time the size of the vector and number of 1s.
#     pack-color will fail if these are wrong, but it will tell you the actual number it found during loading, so you can double
#     check.  cosmo-build reports total bits and set bits in its output
$ pack-color ecoli6_kmc2_list.colors 6  55539132 54489174

# Run bubble caller to load and traverse the succinct colored de Bruijn graph
# cosmo-color  [-b <color_mask2>] [-a <color_mask1>] [-o <output_prefix>] [--] [--version] [-h] <input_file> <color_file> # BubbleCaller
$ cosmo-color -a 1 -b 2 ecoli6_kmc2_list.dbg ecoli6_kmc2_list.sd_vector >bubbles
```
# Legacy Information
Note: Information below this point was in sync with previous states of this software but may be out of date now.  

practical example using the cortex front end:
```sh
$ cd /s/oak/b/nobackup/muggli/src/CORTEX_release_v1.0.5.21/demo/example4_using_reference_genome_to_exclude_paralogs
$ ../../bin/cortex_var_31_c2 --kmer_size 17 --colour_list colours  --dump_binary both.ctx
$ cd ~/git/cosmo
$ ./cosmo-pack -c /s/oak/b/nobackup/muggli/src/CORTEX_release_v1.0.5.21/demo/example4_using_reference_genome_to_exclude_paralogs/both.ctx
$ ./pack-color both.ctx.colors 2
$ ./cosmo-color both.ctx.packed both.ctx.colors.rrr
```

A practical example of the KMC2 front end can be found in cosmo/experiments/ecoli6.sh.  Note that the programs are those in github.com/mmuggli/cosmo.  Some of the commands have changed and will be documented further in the next day or two - MDM, Sept 19. 2016


## Caveats

Here are some things that you don't want to let surprise you:

### DSK Only

Currently Cosmo only supports [DSK][dsk] files with k <= 64 (so, 128 bit or less blocks).
Support is planned for [DSK][dsk] files with larger k, and possibly output from other k-mer
counters.

### Definition of "k-mer"

Note that since our graph is edge-based, k defines the length of our edges, hence our nodes are only k-1 symbols long.
If you want to construct a [Succinct de Bruijn Graph][succ] where the nodes are k-mers, you will need to run [DSK][dsk]
with k set to k+1. E.g. using output from `$ dsk <input_file> 27` will actually build a 26-dimension de Bruijn graph.

*Note: Both even and odd k values should work with this assembler due to our loop-immune traversal.*

Furthermore, most de Bruijn graph based assemblers add edges between *all* nodes that overlap. Instead, we are taking the
k-mers as our edges (of two k-1-length nodes), so we only have edges that were *directly represented in the read set*
(this makes more sense to us, though, as it reduces unnecessary branching). I may add support for the standard way in the
future if anyone wants it (it would be similar to the dummy edge adding code).


### Graph Traversal

We currently only output the unitigs (paths between branching nodes).


## Compilation

There is an included Makefile - just type `make` to build it (assuming you have the dependencies listed below).
To build with "Variable order mode", use the `varord=1` flag.

### Dependencies  
- A compiler that supports C++11,
- [Boost][boost] - ranges and range algorithms, zip iterator, tuple comparison, lots of good stuff,
- [SDSL-lite][sdsl-lite] - low level succinct data structures (For now you will have to use my branch if you want to use variable order
graphs: clone [this](https://github.com/alexbowe/sdsl-lite) and checkout the `develop` branch before compiling),
- [TClap][tclap] - command line parsing,
- [DSK][dsk] - k-mer counting (we need this for input),
- Optionally (for developers): [Python][python] and [NumPy][numpy] - rebuilding the lookup tables,
- [STXXL][stxxl] - external merging (not actually required yet though)
- [KMC2][kmc2] - k-mer counting with sorted output, used to generate the union for colored de Bruijn graph

Many of these are all installable with a package manager (e.g. `(apt-get | yum | brew ) install boost libstxxl tclap`).
However, you will have to download and build these manually: [DSK][dsk] and [SDSL-lite][sdsl-lite].


## Authors

Implemented by [Alex Bowe][abowe]. Original concept and prototype by [Kunihiko Sadakane][ksadakane].
Colored extension prototyped by Robert Raymond and substantially extended by Martin D. Muggli, with help from Alex Bowe.

These people also proved *incredibly* helpful: [Rayan Chikhi][rchikhi], [Simon Puglisi][spuglisi],
[Travis Gagie][tgagie], [Christina Boucher][cboucher], [Simon Gog][sgog], [Dominik Kempa][dkempa].


## Contributing

Your help is more than welcome! Please fork and send a pull request, or contact me directly :)


## Why "Cosmo"?

**Cosmos** */ˈkɑz.moʊs/ (n)* : ["An ordered, harmonious whole."](http://en.wiktionary.org/wiki/cosmos).

If that doesn't suit an assembly program then I don't know what does. The last s was dropped because it was nicer to say.
Furthermore, it is a reference to the Seinfeld character Cosmo Kramer (whose last name I'm often reminded of while working on
this stuff).


## License

This software is copyright (c) Alex Bowe 2014 (bowe dot alexander at gmail dot com).
It is released under the GNU General Public License (GPL) version 3.


[dsk]: http://minia.genouest.org/dsk/
[minia]: http://minia.genouest.org/
[abyss]: https://github.com/bcgsc/abyss
[succ]: http://alexbowe.com/succinct-debruijn-graphs
[debby]: http://github.com/alexbowe/debby

[boost]: http://www.boost.org
[bgl]: http://www.boost.org/doc/libs/1_56_0/libs/graph/doc/
[sdsl-lite]: https://github.com/simongog/sdsl-lite
[networkx]: https://networkx.github.io/
[stxxl]: http://stxxl.sourceforge.net/
[python]: https://www.python.org/
[numpy]: http://www.numpy.org/
[tclap]: http://tclap.sourceforge.net/

[semver]: http://semver.org/
[nucleotides]: http://nucleotid.es/
[tci]: https://travis-ci.org

[abowe]: https://github.com/alexbowe
[cboucher]: http://christinaboucher.com/
[tgagie]: http://www.cs.helsinki.fi/u/gagie/
[ksadakane]: http://researchmap.jp/sada/
[spuglisi]: http://www.cs.helsinki.fi/u/puglisi/
[dkempa]: http://www.cs.helsinki.fi/u/dkempa/
[rchikhi]: https://github.com/rchikhi
[sgog]: https://github.com/simongog/
