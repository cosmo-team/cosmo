
                    .ooooo.   .ooooo.   .oooo.o ooo. .oo.  .oo.    .ooooo.  
                   d88' `"Y8 d88' `88b d88(  "8 `888P"Y88bP"Y88b  d88' `88b 
                   888       888   888 `"Y88b.   888   888   888  888   888 
                   888   .o8 888   888 o.  )88b  888   888   888  888   888 
                   `Y8bod8P' `Y8bod8P' 8""888P' o888o o888o o888o `Y8bod8P' 
                                                                  ver 0.4.3


# Cosmo

[**Version**][semver]: 0.4.3

Cosmo is a fast, low-memory DNA assembler that uses a [succinct de Bruijn graph][succ].


## Usage

After compiling, you can run Cosmo like so:

```sh
$ pack-edges <input_file> # this adds reverse complements and dummy edges, and packs them
$ cosmo-build <input_file>.packed # compresses and builds indices
$ cosmo-assemble <input_file>.packed.dbg # output: <input_file>.packed.dbg.fasta
```

Where `input_file` is the binary output of a [DSK][dsk] run. Each program has a `--help` option for a more
detailed description of how to use them.


## Things to be aware of

Here are some things that you don't want to let surprise you:

### DSK Only

Currently Cosmo only supports [DSK][dsk] files with k <= 64 (so, 128 bit or less blocks).
Support is planned for [DSK][dsk] files with larger k, and possibly output from other k-mer
counters.

### Definition of "k-mer"

Note that since our graph is edge-based, k defines the length of our edges, hence our nodes are only k-1 symbols long.
If you want to construct a [Succinct de Bruijn Graph][succ] where the nodes are k-mers, you will need to run [DSK][dsk]
with k set to k+1. E.g. using output from `$ dsk <input_file> 27` will actually build a 26-dimension de Bruijn graph.

Furthermore, most de Bruijn graph based assemblers add edges between *all* nodes that overlap. We are taking k-mer counter
output for our edges, so we only have edges that were directly represented in the read set (this makes more sense to us, though).
I may add support for the standard way in the future.

### Graph Traversal

The traversal strategy is currently fairly primitive. We only output the unitigs (paths between branches).
Unlike [Minia][minia] (for example), we don't treat each node as equal to its reverse complement (each way
has their pros and cons). In fact, this is actually wrong :/ we need to fix the way it handles reverse complements.


## Overview and Performance

Here is a general overview of each program (details in the upcoming paper):

### pack-edges  
- Ignoring memory requirements, main operations are map - `O(m)`, reduce - `O(m)`, radix sort - `O(mk)`, set difference `O(m)` and merging `O(m)`, for `m` edges;
- Generating all incoming dummies (e.g. `{$ACG} -> {$ACG, $$AC, $$$A}`, so we don't lose any node label data from storing only the last two symbols): `O(dk)` for `d` dummies;
- Adding each dummy shift means that incoming dummies have to be sorted again: `O(dk * k) = O(dk^2)`;
- Since in the worst case `d = m`, total is `O(mk^2)`, but since usually `d << m`, `O(mk)` in practice;
- If this was all implemented in memory, the space requirement would be `m * k * 2 * 2` (we add reverse complements and use a copy-based radix sort)
`+ d * k * 2 = 4mk + 2dk` nucleotides, so `8mk + 4dk` bits (wait...);
- Using the copy-based radix sort is actually a speed optimisation, since it lets us save the second last iteration which we need for the set difference calculations (how we find the required dummies);
- The sort, merge, set difference, map and reduce design of this means it is easy to distribute or make external. Besides, [DSK][dsk] reduces the memory requirement drastically as it is;
- In the output `.packed` file, each edge is represented as five bits (edge symbol + flags) in 64-bit blocks (with four bits wasted per block).

### cosmo-build  
- Constructs the de Bruijn graph *in memory* using succinct data structures that each have linear time construction algorithms - `O(m)`.

### cosmo-assemble  
- Iterates over every edge to build a compressed bit-vector that marks nodes that branch in or out - `O(m)`, since indegree and outdegree are `O(1)`;
- Selects to each branching edge, and follows subsequent edges until it reaches another branch node - `O(m)`.

## Compilation

There is an included Makefile - just type `make` to build it (assuming you have the dependencies listed below).

*Note: it has only been tested on Mac OS X. Changes to work on any *NIX should be minor.*

### Dependencies  
- a compiler that supports C++11,
- [Boost][boost] - ranges and range algorithms, zip iterator, and tuple comparison),
- [STXXL][stxxl] - external merging,
- [SDSL-lite][sdsl-lite] - low level succinct data structures,
- [TClap][tclap] - (command line parsing) libraries installed,
- [DSK][dsk] - we use its output...
- Optionally (for developers): [Python][python] and [NumPy][nympy] - rebuilding the lookup tables.

Many of these are all installable with a package manager (e.g. `(apt-get | yum | brew) boost libstxxl tclap`).
However, you will have to download and build these manually: [DSK][dsk] and [SDSL-lite][sdsl-lite].


## Plan

In no particular order, these are features that I'd like to add:

- Set up Docker image for [nucleotid.es][nucleotides],
- Add [Boost Graph Library][bgl] style API,
- Add support for indirect sorting (to let people attach k-mer counts or colours or whatever people want...) accessible like node/edge properties in [Boost Graph Library][bgl],
- Improve assembly and add error correction (iterative construction),
- Implement dynamic version (necessary for online construction and dynamic error correction),
- Remove alphabet limitation (currently only supports DNA),
- Write unit tests (I have some IPython notebooks that have tests in them, so wasn't completely duct-taped together),
- Set up continuous integration for [Travis CI][tci],
- Add Python wrapper (for learning purposes and Python pipelines) with [NetworkX][networkx] style API.


## Authors

Implemented by [Alex Bowe][abowe]. Original concept and prototype by [Kunihiko Sadakane][ksadakane].

These people also proved incredibly helpful:

- [Rayan Chikhi][rchikhi] - endless advice regarding de Bruijn graphs and assembly in general,
- [Simon Puglisi][spuglisi] - fruitful discussions regarding optimisation,
- [Simon Gog][sgog] - help with [SDSL-lite][sdsl-lite],
- [Dominik Kempa] - help with [STXXL].


## Contributing

Your help is more than welcome! Please fork and send a pull request, or contact me directly :)


## Why "Cosmo"?

I called an earlier version of this Kramer because of its focus on k-mers, and
because I was a fan of the Kramer character on Seinfeld. Kramer's first name happens to be
Cosmo, which sounded cooler to me. It is also a nod to Cosmos (the show that makes science/maths accessible to
regular peeps), and ABySS (which to me sounds spacey as well) - I feel that people exploring genomes are like
small-scale cosmonauts.

You have permission to scoff ;)


## License

This software is copyright (c) Alex Bowe 2014 (bowe dot alexander at gmail dot com).
It is released under the GNU General Public License (GPL) version 3.


[dsk]: http://minia.genouest.org/dsk/
[minia]: http://minia.genouest.org/
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
