
                   .ooooo.   .ooooo.   .oooo.o ooo. .oo.  .oo.    .ooooo.  
                  d88' `"Y8 d88' `88b d88(  "8 `888P"Y88bP"Y88b  d88' `88b 
                  888       888   888 `"Y88b.   888   888   888  888   888 
                  888   .o8 888   888 o.  )88b  888   888   888  888   888 
                  `Y8bod8P' `Y8bod8P' 8""888P' o888o o888o o888o `Y8bod8P' 


# Cosmo

version 1.0


## Description

Cosmo is a tool to convert kmer counter outputs into the packed form suitable for constructing a
[Succinct de Bruijn Graph][succ]. This includes sorting in the required order, and finding and
adding necessary dummy edges.

Currently Cosmo only supports [DSK][dsk] files with k <= 64 (so, 128 bit or less blocks).
Support is planned for [DSK][dsk] files with larger k, and possibly output from other kmer
counters.


## Usage

After building, you can run Cosmo as simply as:

    $ pack_edges <input_file>
    $ cosmo-build <input_file>.packed
    $ cosmo-assemble <input_file>.packed.dbg

Where `input_file` is the binary output of a [DSK][dsk] run (more kmer counters at a later date).
Kramer will detect the k value, and outputs to `<input_file>.packed` by default (which can be altered
with the `-o` option).

Note that if you want to construct a [Succinct de Bruijn Graph][succ] where the nodes are k-mers, you
will need to set DSK's k to k+1. *I might abuse the word "k-mer" in this document... Hopefully it's
clear from context.*

**Example:**

    $ dsk readfile.fq 35 -t 3
    $ kramer readfile.solid_kmers_binary

Will result in a file `readfile.solid_kmers_binary.packed` in the working directory.


## Format

For each edge (a k+1-mer), we need to output the outgoing edge label (one symbol in "$acgt"), a bit flag
which indicate whether this edge's start node (the first k bases) is the first occurence as a start node,
and a bit flag which indicates whether this edge's end node (the last k bases) is the first occurence as an
end node. In total, five bits per edge, packed into 64 bit blocks (12 edges each, with four bits wasted).

The four wasted bits will be at left half of the MSB, and the 12th edge will be in the LSB. The ith edge for a 
block will always be in the same position, even if we dont fill the block up completely.

The edges are followed by five 64 bit integers to store the cumulative counts (which serve as start/end pointers to
the *sorted* second-last symbol in each kmer), and a 64 bit integer representing k (64 bits was chosen to make it
extra easy to parse the file).

A Python script (`parse_kramer.py`) has been provided to demonstrate how to use the output. This can be loaded into
[Debby][debby] (to simulate a [Succinct de Bruijn graph][succ] and recover your full-length kmers if desired).
When we release our fast de Bruijn graph library you will be able to use it with that too.


## Details

Currently Kramer allocates enough space to store all kmers in memory (an external version will be completed later).
All reverse complements are added (as is required for [our graph][succ]'s representation), and then radix sorted into the
correct order - O(mk^2) for m edges - keeping the last and second last iteration - the kmers sorted on start node (say, A),
and the kmers sorted on end node (B) respectively. This may seem wasteful but it let's us very quickly discover the required dummy edges.

Using the two resulting tables, we can detect required incoming dummy edges from the set difference A - B (nodes which appear as start nodes
of edges, but *not* end nodes). Likewise, we can detect required outgoing dummy edges from B - A (nodes which appear as end nodes, but not start nodes).
Both set differences are completed in O(m) time (since the two tables are sorted, and we store kmers in such a way that they can be compared
as integers). We then add all right-shifts of the incoming dummy edges (e.g. $acgt, $$acg, $$$ac, $$$$a) - O(dk) - and sort them - O(dk^2) for d dummies.

In the worst case, there are as many dummy edges as there are edges, which makes our complexity O(mk + mk^2) = O(mk^2), but in reality
there are usually *very few* dummies (especially after outputting only the "solid" kmers, having frequency >= 3, for example).

Finally we 3-way-merge the dummies with the edges sorted on start node, adding the bit flags as we output - O(m).

Since usually d << m, the total run time is O(mk), in 4*m*2k = O(mk) space.

Due to the linear nature of this approach, we can easily break it into mergable chunks or multiple passes, which will allow us to
support files larger than memory at a later date. This also means it should be easy to distribute across multiple computers, or utilise
a GPU, Intel TBB, etc...


## Compilation

There is an included Makefile, so just type `make` to build it.

You will need a compiler that supports C++11, the `Boost` and `TClap` libraries installed, and optionally `Python` (if you want to rebuild the lookup tables)
and `numpy`. These are all installable with any good package manager (e.g. `apt-get`, `yum` or `brew`).


## Contributing

Your help is more than welcome! Please fork and send a pull request, or contact me directly :)


## License

This software is copyright (c) Alex Bowe 2014, bowe.alexander at gmail dot com.
It is released under the GNU General Public License (GPL) version 3.


[dsk]: http://minia.genouest.org/dsk/
[succ]: http://alexbowe.com/succinct-debruijn-graphs
[debby]: http://github.com/alexbowe/debby

