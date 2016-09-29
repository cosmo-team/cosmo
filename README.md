                                                                  88                                 o8o  
                                                                 .8'                                 `"'  
     .ooooo.   .ooooo.   .oooo.o ooo. .oo.  .oo.    .ooooo.     .8'  oooo    ooo  .oooo.   oooo d8b oooo  
    d88' `"Y8 d88' `88b d88(  "8 `888P"Y88bP"Y88b  d88' `88b   .8'    `88.  .8'  `P  )88b  `888""8P `888  
    888       888   888 `"Y88b.   888   888   888  888   888  .8'      `88..8'    .oP"888   888      888  
    888   .o8 888   888 o.  )88b  888   888   888  888   888 .8'        `888'    d8(  888   888      888  
    `Y8bod8P' `Y8bod8P' 8""888P' o888o o888o o888o `Y8bod8P' 88          `8'     `Y888""8o d888b    o888o 
                                                                                                      
                                                                                                ver 0.5.1


# Cosmo

[**Version**][semver]: 0.5.1

Cosmo is a fast, low-memory DNA assembler that uses a [succinct de Bruijn graph][succ].

**VARI**, a succinct colored de Bruijn graph, can be found in the [VARI](https://github.com/cosmo-team/cosmo/tree/VARI) branch.


## Usage

After [compiling](#compilation), you can run Cosmo like so:

```sh
$ pack-edges <input_file> # this adds reverse complements and dummy edges, and packs them
$ cosmo-build <input_file>.packed # compresses and builds indices
$ cosmo-assemble <input_file>.packed.dbg # output: <input_file>.packed.dbg.fasta # NOT IMPLEMENTED YET
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

Many of these are all installable with a package manager (e.g. `(apt-get | yum | brew ) install boost libstxxl tclap`).
However, you will have to download and build these manually: [DSK][dsk] and [SDSL-lite][sdsl-lite].


## Authors

Implemented by [Alex Bowe][abowe]. Original concept and prototype by [Kunihiko Sadakane][ksadakane].

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
