#include "kmer.hpp"
#include "utility.hpp"
#include "debug.hpp"
#include "multi_bit_vector.hpp"
#include "bgl_sdb_adapter.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <boost/filesystem.hpp>
#include <boost/concept/assert.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/breadth_first_search.hpp>
/* the following #define tells the boost::test framework to generate a main() for this module. only one cpp file in the test suite should have it */
#define BOOST_TEST_MODULE BGL SDB Adpater Test
#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>

/*
  Test graph without dummy shifts.

  input = "ATCG", "AACGACGTCGACT", k = 4

       i v d n lbl W
       0 0 0 1 CGA C // A T
       1 1 1 1 AAC G  << add an extra edge to this (AACT)
       2 2 0 1 GAC g
       3   x 0 GAC T
       4 3 1 1 ATC G
       5 4 0 1 GTC g
       6 5 0 1 ACG A
       7   x 0 ACG T
       8 6 0 1 TCG a
       9 7 0 1 ACT $
      10 8 0 1 CGT C
*/

const int sigma    = 4;
string edges       = "CGgTGgATa$C";
string node_flags  = "11101110111";
string dummy_flags = "01011001000"; // mark extra edges with 1s -> can get dummy_idx by subtracting from node count + select works better
vector<string> dummies = vector<string>({ "AAC", "ATC" });
std::array<size_t, 1+sigma> symbol_ends{0, 1, 6, 9, 11};

const size_t g_size = 11;
const size_t num_edges = 10; // outgoing dummy edge is not a real edge
const size_t num_nodes = 9;
string expected_edge_labels[g_size] = {
  "CGAC", "AACG", "GACG", "GACT", "ATCG", "GTCG", "ACGA", "ACGT", "TCGA", "ACT$", "CGTC"
};
string expected_node_labels[num_nodes] = {
//    0     1     2     3     4     5     6     7     8
  "CGA","AAC","GAC","ATC","GTC","ACG","TCG","ACT","CGT"
};

                                         // 0  1  2  3  4  5  6  7  8
const int expected_outdegree[num_nodes] = { 1, 1, 2, 1, 1, 2, 1, 0, 1 };
const int  expected_indegree[num_nodes] = { 2, 0, 1, 0, 1, 2, 2, 1, 1 };

const ssize_t expected_forward[g_size] = {
     // i v d n lbl W
  2, // 0 0 0 1 CGA C
  6, // 1 1 1 1 AAC G
  6, // 2 2 0 1 GAC g
  9, // 3   x 0 GAC T
  8, // 4 3 1 1 ATC G
  8, // 5 4 0 1 GTC g
  0, // 6 5 0 1 ACG A
 10, // 7   x 0 ACG T
  0, // 8 6 0 1 TCG a
 -1, // 9 7 0 1 ACT $
  5, //10 8 0 1 CGT C
};

const ssize_t expected_backward[g_size] = {
     // i v d n lbl W
  6, // 0 0 0 1 CGA C
 -1, // 1 1 1 1 AAC G
  0, // 2 2 0 1 GAC g
  0, // 3   x 0 GAC T
 -1, // 4 3 1 1 ATC G
 10, // 5 4 0 1 GTC g
  1, // 6 5 0 1 ACG A
  1, // 7   x 0 ACG T
  4, // 8 6 0 1 TCG a
  3, // 9 7 0 1 ACT $
  7, //10 8 0 1 CGT C
};

/* test the graph by directly calling the debruijn_graph methods */
template <class DBG>
void test_graph_directly(const DBG & db)
{
  // Edge tests
  for (size_t i = 0;i < g_size;i ++) {
    auto label = db.edge_label(i);
    BOOST_CHECK_MESSAGE(label == expected_edge_labels[i],"edge_label(" << i << ") : " << label << " != " << expected_edge_labels[i]);
    // TODO : also check node values returned by incoming and outgoing
    BOOST_CHECK_MESSAGE(db._backward(i) == expected_backward[i],"backward(" << i << ") : " << db._backward(i) << " != " << expected_backward[i]);
    BOOST_CHECK_MESSAGE(db._forward(i) == expected_forward[i],"forward(" << i << ") : " << db._forward(i) << " != " << expected_forward[i]);
  }

  // Node tests
  for (size_t v = 0; v < num_nodes; v ++) {
    auto label = db.node_label(v);
    BOOST_CHECK_MESSAGE(label == expected_node_labels[v],"node_label(" << v << ") : " << label << " != " << expected_node_labels[v]);

    // for outdegree we will also count the outgoing paths we find via the outgoing function
    int outdegree = db.outdegree(v);
    int indegree = db.indegree(v);
    BOOST_CHECK_MESSAGE(outdegree == expected_outdegree[v],"outdegree(" << v << ") : " << outdegree << " != " << expected_outdegree[v]);
    BOOST_CHECK_MESSAGE(indegree == expected_indegree[v],"indegree(" << v << ") : " << indegree << " != " << expected_indegree[v]);

    // and count the actual outgoing paths
    for (int x = 1;x <= sigma; x++) {
      if (db.outgoing(v,x) >= 0) outdegree--;
      if (db.incoming(v,x) >= 0) {
        //COSMO_LOG(info) << ">>>incoming: " << v << " " << x;
        indegree--;
      }
    }

    BOOST_CHECK_MESSAGE(outdegree == 0, "The outgoing call returns different number of paths to outdegree(" << v << "): " << outdegree);
    // TODO: fix this
    //BOOST_CHECK_MESSAGE(indegree == 0, "The incoming call returns different number of paths to indegree(" << v << "): " << indegree);
  }
}

/* this BFS visitor collects the names of each edge as it is added to the tree */
/*
template<typename EdgeNameMap>
class bfs_edge_collector : public boost::default_bfs_visitor
{
public:
  bfs_edge_collector(EdgeNameMap n_map, std::stringstream & ss) : m_name_map(n_map),ss(ss) {}
  template<typename Edge, typename Graph>
  void tree_edge(Edge e, const Graph&) const {
    ss << get(m_name_map, e);
  }
  std::string str() const { return ss.str(); }
private:
  EdgeNameMap m_name_map;
  std::stringstream & ss;
};
*/

/* test the graph by calling it via boost graph library methods */
using namespace boost;
template<typename Graph>
void test_graph_via_boost(Graph & db)
{
  BOOST_CONCEPT_ASSERT(( GraphConcept<Graph> ));
  BOOST_CONCEPT_ASSERT(( IncidenceGraphConcept<Graph> ));
  BOOST_CONCEPT_ASSERT(( VertexListGraphConcept<Graph> ));
  BOOST_CONCEPT_ASSERT(( ReadablePropertyMapConcept<typename property_map<Graph, vertex_name_t>::type,size_t> ));
  BOOST_CONCEPT_ASSERT(( ReadablePropertyMapConcept<typename property_map<Graph, edge_name_t>::type,size_t> ));
  
  for (size_t i = 0;i < num_nodes; ++i)
  {
    // the easy way to read the vertex_name_t property
    BOOST_CHECK_MESSAGE(get(vertex_name_t(),db,i) == expected_node_labels[i],"node_label(" << i << ") failed");

    int outdegree = out_degree(i, db);
    //int indegree  = in_degree(i, db);
    BOOST_CHECK_MESSAGE(outdegree == expected_outdegree[i],"outdegree(" << i << ") failed");
    //BOOST_CHECK_MESSAGE(in_degree == expected_indegree[i],"indegree(" << i << ") failed");
  }

  //auto edge_name_map = get(edge_name_t(), db);
  /*
  for (size_t i = 0; i < g_size; ++i) {
    std::cerr << get(edge_name_map, i);
  }
  */

  // NOTE: BFS commented out as it is hard to know which edges will be chosen first in a bigger test
  // prepare the visitor for a BFS
  //std::stringstream ss;
  //typename property_map<Graph, edge_name_t>::type edge_name_map = get(edge_name_t(), db);
  //bfs_edge_collector<typename property_map<Graph , boost::edge_name_t>::type > edge_visitor(edge_name_map,ss);

  // do the BFS
  //breadth_first_search(db, 0, boost::visitor(edge_visitor));
  //BOOST_CHECK_EQUAL(ss.str(), "AACGATCCTG");
}

/* builds the graph and calls the test functions */
BOOST_AUTO_TEST_CASE ( create_graph_using_constructor )
{
  size_t k = 4;
  std::string alphabet = "$ACGT";

  bit_vector node_bv(node_flags.length(), 1);
  bit_vector dummy_bv(dummy_flags.length(), 0);
  int_vector<8> edges_v(edges.length(), 0);

  for (unsigned int i = 0; i < edges.length(); ++i) {
    if (node_flags[i] == '1') node_bv[i] = 0;
    if (dummy_flags[i] == '1') dummy_bv[i] = 1;
    edges_v[i] = flagged_nt_to_int(edges[i]);
  }
  sd_vector<> node_flags(node_bv);
  sd_vector<> dummy_flags(dummy_bv);
  typedef wt_huff<rrr_vector<63>> edge_idx_t;
  //typedef cosmo::multi_bit_vector<> edge_idx_t;
  typedef debruijn_graph<edge_idx_t> dbg_t;
  edge_idx_t edge_idx;
  construct_im(edge_idx, edges_v);

  vector<kmer_t> dummy_v;
  dummy_v.reserve(dummies.size());
  for (auto dummy : dummies) {
    auto x = string_to_kmer<kmer_t>(dummy);
    dummy_v.push_back(x);
  }

  // construct the graph!
  dbg_t db(k, node_flags, edge_idx, symbol_ends, alphabet, dummy_bv, dummy_v);

  test_graph_directly(db);

  // Serialization
  std::string temp_file = boost::filesystem::unique_path().native();
  dbg_t db2;
  sdsl::store_to_file(db, temp_file);
  sdsl::load_from_file(db2, temp_file);
  boost::filesystem::remove(temp_file);
  test_graph_directly(db2);

  test_graph_via_boost(db);
}
