#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debruijn_graph_shifted.hpp"
#include "algorithm.hpp"
#include "vari-merge.hpp"

//using namespace std;
//using namespace sdsl;

#include <sys/timeb.h>


int getMilliCount(){
  timeb tb;
  ftime(&tb);
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}


int getMilliSpan(int nTimeStart){
  int nSpan = getMilliCount() - nTimeStart;
  if(nSpan < 0)
    nSpan += 0x100000 * 1000;
  return nSpan;
}

void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd(BANNER, ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input", ".dbg file.", true, "", "graph_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg2("input2", ".dbg file.", true, "", "graph_file", cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.input_filename2  = input_filename_arg2.getValue();

}


void dump_edges(const debruijn_graph_shifted<> &dbg) {
  for (size_t i = 0; i < dbg.size(); i++) {
      cout <<  /* << "e:" <<*/ dbg.edge_label(i) << std::endl;
  }
}




// Returns the number of 0's in a run beginning at set_start
uint64_t length(const uint64_t set_start, const std::vector<bool> &sets)
{
    uint64_t i = 0;
    for ( i= set_start; sets[i] != true; ++i)
        ;
    return i - set_start;
}

// Returns the next set_start given some set_start
uint64_t advance(const uint64_t set_start, const std::vector<bool> &sets)
{
    return length(set_start, sets) + 1;
}

// Returns true if the given set_start is in the last set
bool last(const uint64_t set_start, const std::vector<bool> &sets)
{

    return set_start + length(set_start, sets)  == sets.size() - 1;
}

// Returns the number of sets in positions greater or equal to start
uint64_t num_sets_ge(const uint64_t set_start, const std::vector<bool> &sets)
{
    assert (set_start <= sets.size() - 1);
    uint64_t count = 0;
    for (uint64_t i = set_start; !last(i, sets); i += advance(i, sets)) {
        count++;
    }
    return count + 1;

}

template<class E>
int dump_range(const uint64_t start, const uint64_t end, const std::vector<E> &col)
{
    std::cout << "[";
    for (uint64_t i = start; i < end; ++i)
        std::cout << col[i] << ", ";
    std::cout << "]";
    return 0;
}


uint64_t num_sets(const std::vector<bool> &sets)
{
    return num_sets_ge(0, sets);
}


uint64_t validate(const std::vector<bool> &sets)
{
    uint64_t sum = 0;
    uint64_t num = num_sets(sets);
    uint64_t i = 0;
    uint64_t pos = 0;
        
    for (; i < num - 1; ++i) {
        sum += length(pos, sets);
        pos += advance(pos, sets);
    }
    assert(last(pos, sets));
    return sum + length(pos, sets);
}
//FIXME: this is pretty inefficient; we should be able to scan through the columns only once
int subdivide(const std::vector<char> &g1_col, const uint64_t g1_ptr, const uint64_t g1_num,
              const std::vector<char> &g2_col, const uint64_t g2_ptr, const uint64_t g2_num,
              /*const int colno,*/ std::vector<bool> &g1_out_set, std::vector<bool> &g2_out_set, int &active_alpha_size)
{    
    std::set<char> g1_chars;
    for (uint64_t i = 0; i <  g1_num; ++i) {
        g1_chars.insert(g1_col[g1_ptr + i]);
    }
    std::set<char> g2_chars;    
    for (uint64_t i = 0; i <  g2_num; ++i) {
        g2_chars.insert(g2_col[g2_ptr +i]);
    }
    std::vector<char> g1_g2_union;
    std::set_union(g1_chars.begin(), g1_chars.end(),
                   g2_chars.begin(), g2_chars.end(),                  
                   std::back_inserter(g1_g2_union));
    
    // # for each $ in column, if the corresponding edge label exists in the other set, we can delete it
            
    // for letter in active_alphabet:
     for (auto letter: g1_g2_union) {
         g1_out_set.insert(g1_out_set.end(), 
                           count_if (g1_col.begin()+g1_ptr, g1_col.begin()+g1_ptr+g1_num/*+1*/,  [letter](char i){return i  == letter;}),
                           false);
         g1_out_set.push_back(true);
         g2_out_set.insert(g2_out_set.end(), 
                           count_if (g2_col.begin()+g2_ptr, g2_col.begin()+g2_ptr+g2_num/*+1*/,  [letter](char i){return i  == letter;}),
                           false);
         g2_out_set.push_back(true);
     }
     active_alpha_size = g1_g2_union.size();
     return 0;
}


int refine_sets(const std::vector<char> &g1_col, const std::vector<char> &g2_col,
                const std::vector<bool>& g1_sets, const std::vector<bool> &g2_sets,
                const int colno,
                std::vector<bool> &g1_out_set, std::vector<bool> &g2_out_set, std::vector<bool> &Lval)
{

    uint64_t g1_ptr = 0;
    uint64_t g2_ptr = 0;    
    
    uint64_t g1_set_start = 0;
    uint64_t g2_set_start = 0;
    std::cout << "subdividing " << num_sets(g1_sets) << " and " << num_sets(g2_sets) << " sets" << std::endl;
    do {
        uint64_t g1_num = length(g1_set_start, g1_sets);
        uint64_t g2_num = length(g2_set_start, g2_sets);
        int active_alpha_size = 0;
        std::cout << "subdividing ranges (col " << colno <<" ) " << std::endl << "\t" << g1_ptr << ":+" << g1_num << " = " ;
        dump_range(g1_ptr, g1_ptr+ g1_num, g1_col);
        std::cout << std::endl << "\t" << g2_ptr << ":+" << g2_num;
        dump_range(g2_ptr, g2_ptr+ g2_num, g2_col);
        std::cout <<std::endl;
        uint64_t g1_out_set_initsize = g1_out_set.size();
        uint64_t g2_out_set_initsize = g2_out_set.size();        
        subdivide(g1_col, g1_ptr, g1_num,
                  g2_col, g2_ptr, g2_num,
//                  colno,
                  g1_out_set,
                  g2_out_set,
                  active_alpha_size);
        std::cout << "\talpha size: " << active_alpha_size << std::endl;
        std::cout << "\tg1_out_set size: " << g1_out_set.size() << ", g2_out_set size: " << g2_out_set.size() << std::endl;
        std::cout << "validations  " << validate(g1_out_set) << " " << validate(g2_out_set) << std::endl;
        std::cout << "\tg1_out_set additional sets: " << ((g1_out_set_initsize > g1_out_set.size()) ? (num_sets_ge(g1_out_set_initsize + 1, g1_out_set)) : 0)
                  << ", g2_out_set additional sets: " << ((g2_out_set_initsize > g2_out_set.size()) ? num_sets_ge(g2_out_set_initsize + 1, g2_out_set) : 0) << std::endl;        
                  
        if (colno == 0) {
            Lval.push_back(true);
            Lval.insert(Lval.end(), active_alpha_size - 1, false);
        }
                

        g1_ptr += g1_num;
        g2_ptr += g2_num;
        assert(g1_ptr <= g1_col.size());
        assert(g2_ptr <= g2_col.size());
            
        if (last(g1_set_start, g1_sets) || last(g2_set_start, g2_sets)) break;
        
        g1_set_start += advance(g1_set_start, g1_sets);
        g2_set_start += advance(g2_set_start, g2_sets);
    } while (true);

    return 0;
}

int get_column(const debruijn_graph_shifted<> &g, const int col_num, std::vector<char> &g_col)
{
    assert (g_col.size() == 0);
    for (uint64_t i = 0; i < g.num_edges(); ++i) {
        if (col_num == 0 ) {
            g_col.push_back(g._map_symbol(g._strip_edge_flag(g.m_edges[i])));
        } else {
            const char *lab = g.edge_label(i).c_str();
            if (g.k - col_num - 1>= strlen(lab)) {
                g_col.push_back('$');
            } else {
                char c = lab[g.k - col_num - 1];
                g_col.push_back(c ? c : '$');
            }
        }
    }
    return 0;
}

class Flags {

public:
    Flags(const std::vector<bool> &g1_flagsets, const std::vector<bool> &g2_flagsets);
    const std::vector<bool> &g1_flagsets;
    const std::vector<bool> &g2_flagsets;
    bool seen(char nt);
    void add(char nt);
    void adv_g1();
    void adv_g2();
private:
    std::set<char> newflags;
    void check_flagsets();
    uint64_t g1_ptr = 0;
    uint64_t g2_ptr = 0;
    uint64_t g1_base = 0;
    uint64_t g2_base = 0;
};

Flags::Flags(const std::vector<bool> &g1_flagsets_a, const std::vector<bool> &g2_flagsets_a) : g1_flagsets(g1_flagsets_a), g2_flagsets(g2_flagsets_a)
{

}

bool Flags::seen(char c)
{
    return newflags.count(c);
}


void Flags::add(char c)
{
    newflags.insert(c);
}


void Flags::check_flagsets()
{
    if (g1_flagsets[g1_base + g1_ptr] && g2_flagsets[g2_base + g2_ptr]) {
        g1_base += g1_ptr + 1;
        g2_base += g2_ptr + 1;
        g1_ptr = 0;
        g2_ptr = 0;
        newflags.clear();
    }
}
    
void    Flags::adv_g1()
{
    g1_ptr += 1;
    check_flagsets();
}

void    Flags::adv_g2()
{
    g2_ptr += 1;
    check_flagsets();
}



void fill_Lcol(const debruijn_graph_shifted<> &g, std::vector<bool> &Lcol)
{
    size_t p = 1;
    size_t Lindex = 0;
    do {
        Lindex = g.m_node_select(p);
        Lcol[Lindex] = true;
        p += 1;
    } while (Lindex < Lcol.size() - 16 /*FIXME: this is assumed to be realted to other truncation*/);
}

int mainmerge(const debruijn_graph_shifted<> &g1, const debruijn_graph_shifted<> &g2)
{

    std::vector<bool> g1_sets(g1.num_edges() + 1, false);
    g1_sets[g1_sets.size() - 1] = true;
    std::cout << "created g1_sets with " << num_sets(g1_sets) << " sets of size " << length(0, g1_sets) << std::endl;

    std::vector<bool> g2_sets(g2.num_edges() + 1, false);
    g2_sets[g2_sets.size() - 1] = true;
    std::cout << "created g2_sets with " << num_sets(g2_sets) << " sets of size " << length(0, g2_sets) << std::endl;
    std::cout << std::endl;

    std::vector<bool> L;

    assert(g1.k == 10); // FIXME: remove after prototyping, just want to make sure behavior matches Debby's
    

    std::vector<int> cols;
    for (int i = 1; i < g1.k ; ++i) {
        cols.push_back(i);
    }
    cols.push_back(0);

    std::vector<bool> g1_flagsets;//(/*g1_sets*/);
    std::vector<bool> g2_flagsets;//(/*g2_sets*/);

    std::vector<bool> Lcol(g1.num_edges(), false);
    fill_Lcol(g1, Lcol);
    for (auto col: cols) {


        std::vector<char> g1_col; // FIXME: change 'char' type to something less static
        get_column(g1, col, g1_col);
        std::vector<char> g2_col;
        get_column(g2, col, g2_col);


        if (col == 1) {
            std::vector<char> g1_col1(g1_col);
            std::vector<char> g2_col1(g2_col);
        }
        
        std::vector<bool> g1_new_sets;
        std::vector<bool> g2_new_sets;
        std::cout << "going to refine sets in bool vectors of sizes " << g1_sets.size() << ", " << g2_sets.size() << ";   " << std::endl << "calling refine_sets()\n";
        refine_sets(g1_col, g2_col, g1_sets, g2_sets, col, g1_new_sets, g2_new_sets, L);
        std::cout <<"    got back new vectors of bools of sizes " << g1_new_sets.size() << ", " << g2_new_sets.size() << ";   " << std::endl << std::endl;
        assert(g1_col.size() == validate(g1_new_sets));
        assert(g2_col.size() == validate(g2_new_sets));
        std::cout << "validation passed g1: " << validate(g1_new_sets)
                  << " g2: " << validate(g2_new_sets) << std::endl;
        //FIXME: avoid copying values in this next part
        g1_sets.assign(g1_new_sets.begin(), g1_new_sets.end()); 
        g2_sets.assign(g2_new_sets.begin(), g2_new_sets.end()); 

        if (col == g1.k - 2) {
            g1_flagsets.insert(g1_flagsets.end(), g1_sets.begin(), g1_sets.end());
            g2_flagsets.insert(g2_flagsets.end(), g2_sets.begin(), g2_sets.end());
        }
    }        
    //FIXME: find some other name than pointer, since really not
    uint64_t g1_ptr = 0;
    uint64_t g2_ptr = 0;
    
    Flags flags(g1_flagsets, g2_flagsets);
        std::cout << "flags: ";
        dump_range(0, g1_flagsets.size(), g1_flagsets);
        std::cout << std::endl;
        dump_range(0, g2_flagsets.size(), g2_flagsets);
        std::cout << std::endl;
        

        std::map<char, uint64_t> ntcounts;
        ntcounts['$'] = 0;
        ntcounts['A'] = 0;
        ntcounts['C'] = 0;
        ntcounts['G'] = 0;
        ntcounts['T'] = 0;



        uint64_t g1_set_ptr = 0;
        uint64_t g2_set_ptr = 0;
        uint64_t out_ptr = 0;
        std::cout << "debug Lcol.size = " << Lcol.size() << " L.size = " << L.size() << std::endl;
        do {
            while (!g1_sets[g1_set_ptr]) ++g1_set_ptr;
            while (!g2_sets[g2_set_ptr]) ++g2_set_ptr;

            if (g1_set_ptr > 0 && !g1_sets[g1_set_ptr - 1]) {
                
                std::cout << (int)Lcol[g1_ptr] << " :: ";
                std::cout << L[out_ptr] << " " << g1._map_symbol(g1._strip_edge_flag(g1.m_edges[g1_ptr])) << " "; // FIXME: assuming no flags in m_edges
                ntcounts[g1._map_symbol(g1._strip_edge_flag(g1.m_edges[g1_ptr]))] += 1;
                if (flags.seen(g1._map_symbol(g1._strip_edge_flag(g1.m_edges[g1_ptr])))) {
                    std::cout << 1 << " == " << g1.edge_label(g1_ptr) << " " << g1_ptr /*<< " : " << (int)(g1.m_edges[g1_ptr] & 0x1) */ << std::endl;
                } else {
                    std::cout << 0 << " == " << g1.edge_label(g1_ptr) << " " << g1_ptr/*<< " : " << (g1.m_edges[g1_ptr] & 0x1) */ << std::endl;
                }
                flags.add(g1._map_symbol(g1._strip_edge_flag(g1.m_edges[g1_ptr])));
                g1_ptr += 1;
                flags.adv_g1();
                if (g2_set_ptr > 0 && !g2_sets[g2_set_ptr - 1]) {
                    g2_ptr += 1;
                    flags.adv_g2();
                }

            } else {
                std::cout << (int)Lcol[g2_ptr] << " :: ";
                std::cout << L[out_ptr] << g2._map_symbol(g2._strip_edge_flag(g2.m_edges[g2_ptr])) << " ";
                ntcounts[g2._map_symbol(g2._strip_edge_flag(g2.m_edges[g2_ptr]))] += 1;
                if (flags.seen(g2._map_symbol(g2._strip_edge_flag(g2.m_edges[g2_ptr])))) {
                    std::cout << 1 /*<< (int)(g1.m_edges[g1_ptr] & 0x1)*/ << std::endl;
                } else {
                    std::cout << 0 /*<< (int)(g1.m_edges[g1_ptr] & 0x1)*/ << std::endl;
                }
                flags.add(g2._map_symbol(g2._strip_edge_flag(g2.m_edges[g2_ptr])));
                g2_ptr += 1;
                flags.adv_g2();
            }
            ++out_ptr;
            ++g1_set_ptr;
            ++g2_set_ptr;
        } while (g1_set_ptr != g1_sets.size()  && g2_set_ptr != g2_sets.size() );
        std::cout << "set pointers " << g1_set_ptr << "/" <<g1_sets.size() << " " << g2_set_ptr << "/" << g2_sets.size() << std::endl;
        std::cout << "elements in family: " << validate(g1_sets) << " " << validate(g2_sets) << std::endl;
        assert (g1_set_ptr == g1_sets.size() );
        assert ( g2_set_ptr == g2_sets.size() );
            

        std::cout << 0 << " "
                  << ntcounts['$'] << " "
                  << ntcounts['$'] + ntcounts['A'] << " "
                  << ntcounts['$'] + ntcounts['A'] + ntcounts['C']   << " "
                  << ntcounts['$'] + ntcounts['A'] + ntcounts['C'] + ntcounts['G']  << " "
                  << std::endl;
        std::cout << "m_symbol_ends for g1: ";
        std::cout  << g1.k << std::endl;
            
    return 0;
}


void dumpcolumns( const debruijn_graph_shifted<> &dbg, const  parameters_t &p)
{
    for (int i = 0; i < dbg.k; ++i) {
        std::stringstream fname;
        fname <<  p.input_filename;
        fname << i;
        std::string s = fname.str();
        ofstream f(s.c_str());
        std::vector<char> g_col;
        get_column(dbg, i, g_col);
        for (auto c: g_col) {
            f << c << std::endl;
        }
        f.close();
    }

}
int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);
  cerr << "pack-color compiled with supported colors=" << NUM_COLS << std::endl;
  //ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  // Can add this to save a couple seconds off traversal - not really worth it.
  cerr << "loading dbg "<< p.input_filename  << std::endl;
  debruijn_graph_shifted<> dbg;
  load_from_file(dbg, p.input_filename);
  //input.close();

  cerr << "loading dbg2 " << p.input_filename2 << std::endl;
  debruijn_graph_shifted<> dbg2;
  load_from_file(dbg2, p.input_filename2);

  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
  dumpcolumns(dbg, p);
  dump_edges(dbg);
  mainmerge(dbg, dbg2);

}
