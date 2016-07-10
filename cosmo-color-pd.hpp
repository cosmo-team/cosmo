#ifndef COSMO_COLOR_H
#define COSMO_COLOR_H

struct parameters_t {
    std::string input_filename = "";
    std::string color_filename = "";
    std::string output_prefix = "";
    std::string ref_color = "";
    std::string sample_mask = "";
    std::string ref_fasta = "";
    std::string output_matrix = "";
};

int getMilliCount();
int getMilliSpan(int nTimeStart);
void parse_arguments(int argc, char **argv, parameters_t & params);
void test_symmetry(debruijn_graph<> dbg);
void dump_nodes(debruijn_graph<> dbg, uint64_t * colors);
void dump_edges(debruijn_graph<> dbg, uint64_t * colors);
void find_bubbles(debruijn_graph<> dbg, rrr_vector<63> &colors, uint64_t color_mask1, uint64_t color_mask2);
#endif
