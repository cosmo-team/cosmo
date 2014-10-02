

                    .ooooo.   .ooooo.   .oooo.o ooo. .oo.  .oo.    .ooooo.  
                   d88' `"Y8 d88' `88b d88(  "8 `888P"Y88bP"Y88b  d88' `88b 
                   888       888   888 `"Y88b.   888   888   888  888   888 
                   888   .o8 888   888 o.  )88b  888   888   888  888   888 
                   `Y8bod8P' `Y8bod8P' 8""888P' o888o o888o o888o `Y8bod8P' 
                                                                  ver 0.4.4


# Cosmo

[**Version**][semver]: 0.4.4

Cosmo is a fast, low-memory DNA assembler that uses a [succinct de Bruijn graph][succ].


## Usage

After [compiling](#compilation), you can run Cosmo like so:

```sh
$ pack-edges <input_file> # this adds reverse complements and dummy edges, and packs them
$ cosmo-build <input_file>.packed # compresses and builds indices
$ cosmo-assemble <input_file>.packed.dbg # output: <input_file>.packed.dbg.fasta
```

Where `input_file` is the binary output of a [DSK][dsk] run. Each program has a `--help` option for a more
detailed description of how to use them.


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

*Note: it has only been tested on Mac OS X. Changes to work on any *NIX should be minor.*

### Dependencies  
- A compiler that supports C++11,
- [Boost][boost] - ranges and range algorithms, zip iterator, and tuple comparison),
- [STXXL][stxxl] - external merging,
- [SDSL-lite][sdsl-lite] - low level succinct data structures,
- [TClap][tclap] - command line parsing,
- [DSK][dsk] - k-mer counting (we need this for input),
- Optionally (for developers): [Python][python] and [NumPy][numpy] - rebuilding the lookup tables.

Many of these are all installable with a package manager (e.g. `(apt-get | yum | brew ) install boost libstxxl tclap`).
However, you will have to download and build these manually: [DSK][dsk] and [SDSL-lite][sdsl-lite].


## Overview and Performance

Here is a general overview of each program (details in the upcoming paper):

### pack-edges  
- Ignoring memory requirements, main operations are map - O(m), reduce - O(m), radix sort - O(mk), set difference - O(m), and merging - O(m), for m edges.
- Generating all incoming dummies (e.g. {$ACG} -> {$ACG, $$AC, $$$A}, so we don't lose any node label data from storing only the last two symbols) - O(dk) for d dummies.
- Adding each dummy shift means that incoming dummies have to be sorted again: O(dk * k) = O(dk^2).
- Since in the worst case d = m, total is O(mk^2), but since usually d << m, O(mk) in practice.
- If this was all implemented in memory, the space requirement would be m * k * 2 * 2 (we add reverse complements and use a copy-based radix sort) + d * k * 2 = 4mk + 2dk
nucleotides, so 8mk + 4dk bits (which might sound like a lot, but...),
- Using the copy-based radix sort is actually a speed optimisation, since it lets us save the second last iteration which we need for the set difference calculations (how we find the required dummies).
- The sort, merge, set difference, map and reduce design of this means it is easy to distribute or make external. Besides, [DSK][dsk] reduces the memory requirement drastically as it is.
- In the output `.packed` file, each edge is represented as five bits (edge symbol + flags) in 64-bit blocks (with four bits wasted per block).

### cosmo-build  
- Constructs the de Bruijn graph *in memory* using succinct data structures that each have linear time construction algorithms - O(m).

### cosmo-assemble  
- Iterates over every edge to build a compressed bit-vector that marks nodes that branch in or out - O(m), since indegree and outdegree are O(1);
- Selects to each branching edge, and follows subsequent edges until it reaches another branch node - O(m).


## .plan

### Upcoming Release

- [ ] Handle the first k symbols of incoming tips (backtrack until $ or k-1 - we want the first edge to make a kmer)
- [ ] Handle palindromic edges (detect palindrome in last k edges followed and exit)
- [ ] Write contig postprocessing system (store only one of the pairs of contigs - compare start to twin(end), store min)
- [ ] Rewrite edge vector so it is faster (currently a wavelet tree - fine for general case)
  - [ ] Vector with four bits for each node, with rank/select only on the non-minus flagged edges? (potential problem with sampling)
  - [ ] Canonical Huffman coding (2^4 ints for frequencies, map the prefix code)
- [ ] Add support for external sorting (for large data sets)

### Future Releases

- Set up Docker image for [nucleotid.es][nucleotides],
- Add [Boost Graph Library][bgl] style API (to get merged into a heavyweight assembler),
- Prepare documentation for said API,
- Improve memory use for dummy edge generation (many of the shifted dummies are repeated... could build a trie instead),
- Add support for indirect sorting (to let people attach k-mer counts or colours or whatever people want...) accessible like node/edge properties in [Boost Graph Library][bgl],
- Improve assembly and add error correction (iterative construction),
- Implement dynamic version (necessary for online construction and dynamic error correction),
- Remove alphabet limitation (currently only supports DNA),
- Write unit tests and refactor (I have some IPython notebooks that have tests in them, so this wasn't completely duct-taped together),
- Set up continuous integration for [Travis CI][tci],
- Add Python wrapper (for learning purposes and Python pipelines) with [NetworkX][networkx] style API.


## Authors

Implemented by [Alex Bowe][abowe]. Original concept and prototype by [Kunihiko Sadakane][ksadakane].

These people also proved incredibly helpful:

- [Rayan Chikhi][rchikhi] - endless advice regarding de Bruijn graphs and assembly in general,
- [Simon Puglisi][spuglisi] - fruitful discussions regarding optimisation,
- [Simon Gog][sgog] - help with [SDSL-lite][sdsl-lite],
- [Dominik Kempa][dkempa] - help with [STXXL].


## Contributing

Your help is more than welcome! Please fork and send a pull request, or contact me directly :)


## Why "Cosmo"?

**Cosmos** */ˈkɑz.moʊs/ (n)* : ["An ordered, harmonious whole."](http://en.wiktionary.org/wiki/cosmos).

If that doesn't suit an assembly program then I don't know what does. The last s was dropped because it was nicer to say.
Furthermore, it is a reference to the Seinfeld character Cosmo Kramer (whose last name I'm often reminded of while working on
this stuff). Finally, it's a nod to the [ABySS][abyss] assembler's name, which also makes me think of The Universe.


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
[ksadakane]: http://researchmap.jp/sada/
[spuglisi]: http://www.cs.helsinki.fi/u/puglisi/
[dkempa]: http://www.cs.helsinki.fi/u/dkempa/
[rchikhi]: https://github.com/rchikhi
[sgog]: https://github.com/simongog/
