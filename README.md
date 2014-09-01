
                    .ooooo.   .ooooo.   .oooo.o ooo. .oo.  .oo.    .ooooo.  
                   d88' `"Y8 d88' `88b d88(  "8 `888P"Y88bP"Y88b  d88' `88b 
                   888       888   888 `"Y88b.   888   888   888  888   888 
                   888   .o8 888   888 o.  )88b  888   888   888  888   888 
                   `Y8bod8P' `Y8bod8P' 8""888P' o888o o888o o888o `Y8bod8P' 


# Cosmo

version 1.0


## Description

Cosmo is a fast, low-memory DNA assembler using a [Succinct de Bruijn Graph][succ].


## Usage

After compiling, you can run Cosmo as simply as:

    $ pack-edges <input_file> # this adds reverse complements and dummy edges, and packs them
    $ cosmo-build <input_file>.packed # compresses and builds indices
    $ cosmo-assemble <input_file>.packed.dbg # output: <input_file>.packed.dbg.fasta

Where `input_file` is the binary output of a [DSK][dsk] run. Each program has a `--help` option for a more
detailed description.


## Things to be aware of

Here are some things that you don't want to let surprise you:

### DSK Only

Currently Cosmo only supports [DSK][dsk] files with k <= 64 (so, 128 bit or less blocks).
Support is planned for [DSK][dsk] files with larger k, and possibly output from other kmer
counters.

### Definition of "K-mer"

Note that since our graph is edge-based, `k` defines the length of our edges, hence our nodes are only `k-1` symbols long.
If you want to construct a [Succinct de Bruijn Graph][succ] where the nodes are `k`-mers, you will need to run [DSK][dsk]
with k set to k+1. E.g. using output from `$ dsk <input_file> 27` will actually build a `26`-dimension de Bruijn graph.

At the time of writing, [DSK][dsk] doesn't support even `k` values though... so I had better hurry up and support different
k-mer counters.

Furthermore, most de Bruijn graph based assemblers add edges between *all* nodes that overlap. We are taking k-mer counter
output for our edges, so we only have edges that were directly represented in the read set (this makes more sense to us, though -
why add information that wasn't there, in a non-controlled way? Yes, we are *slightly* more dependent on read quality, but adding all
possible edges isn't a sensible way to error correct).

### Graph Traversal

The traversal strategy is currently fairly primitive. We only output the unitigs (paths between branches).
Unlike [Minia][minia] (for example), we don't treat each node as equal to its reverse complement (each way
has their pros and cons).


## Overview and Performance

Here is a general overview of each program (details in the upcoming paper):

### pack-edges  
- Ignoring memory requirements, main operations are map - `O(m)`, reduce - `O(m)`, radix sort - `O(mk)`, set difference `O(m)` and merging `O(m)`, for `m` edges;
- Generating all incoming dummies (e.g. `{$ACG} -> {$ACG, $$AC, $$$A}`, so we don't lose any node label data from storing only the last two symbols): `O(dk)` for `d` dummies;
- Adding each dummy shift means that incoming dummies have to be sorted again: `O(dk * k) = O(dk^2)`;
- Since in the worst case `d = m`, total is `O(mk^2)`, but since usually `d << m`, `O(mk)` in practice;
- If this was all implemented in memory, the space requirement would be `m * k * 2 * 2` (we add reverse complements and use a copy-based radix sort)
`+ d * k * 2 = 4mk + 2dk` nucleotides, so `8mk + 4dk` bits (wait...);
- Using the copy-based radix sort is actually a speed optimization, since it lets us save the second last iteration which we need for the set difference calculations (how we find the required dummies);
- The sort, merge, set difference, map and reduce design of this means it is easy to distribute or make external. Besides, [DSK][dsk] reduces the memory requirement drastically as it is;
- In the output `.packed` file, each edge is represented as `5` bits (edge symbol + flags) in `64`-bit blocks (with `4` bits wasted per block).

### cosmo-build  
- Constructs the de Bruijn graph *in memory* using succinct data structures that each have linear time construction algorithms - `O(m)`.

### cosmo-assemble  
- Iterates over every edge to build a compressed bit-vector that marks nodes that branch in or out - `O(m)`, since indegree and outdegree are `O(1)`;
- Selects to each branching edge, and follows subsequent edges until it reaches another branch node - `O(m)`.

## Compilation

There is an included Makefile - just type `make` to build it.

You will need a compiler that supports C++11, the `Boost` (ranges and range algorithms, zip iterator, and tuple comparison) `libstxxl` (external merging), 
`sdsl-lite` (low level succinct data structures), and `TClap` (command line parsing) libraries installed,
and optionally `Python` and `numpy` (to rebuild the lookup tables).

These are all installable with any good package manager (e.g. `apt-get`, `yum` or `brew`), except for `sdsl-lite`, which you should [download][sdsl] and build manually.


## Authors

Implemented by Alex Bowe. Original concept and prototype by Kunihiko Sadakane.

These people also proved incredibly helpful:

- Simon Puglisi - Fruitful discussions regarding optimimisation
- Dominik Kempa - Help with STXXL
- Rayan Chikhi - Endless advice regarding de Bruijn graphs and assembly in general
- Simon Gog - support with SDSL


## Contributing

Your help is more than welcome! Please fork and send a pull request, or contact me directly :)


## Why "Cosmo"?

I called an earlier version of this Kramer because of its focus on k-mers, and
because I was a fan of the Kramer character on Seinfeld. Kramer's first name happens to be
Cosmo, which sounded cooler to me. It is also a nod to Cosmos (the show that makes science/maths accessible to
regular peeps), and ABySS (which to me sounds spacey as well) - I feel that people exploring genomes are contemporary
cosmonauts.

You have permission to scoff ;)


## License

This software is copyright (c) Alex Bowe 2014 (bowe dot alexander at gmail dot com).
It is released under the GNU General Public License (GPL) version 3.


[dsk]: http://minia.genouest.org/dsk/
[minia]: http://minia.genouest.org/
[succ]: http://alexbowe.com/succinct-debruijn-graphs
[debby]: http://github.com/alexbowe/debby
[sdsl]: https://github.com/simongog/sdsl-lite

