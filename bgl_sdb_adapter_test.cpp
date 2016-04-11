/*

  Thoroughly tests the debruijn graph for the graph in the paper
  http://alexbowe.com/succinct-debruijn-graphs/
 
  Does it both directly via calls to the debruijn_graph methods, and also via the boost graph library interfaces
 
*/

#include "kmer.hpp"
#include "utility.hpp"
#include "debug.hpp"
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
#include <boost/test/included/unit_test.hpp>

/*
  Test graph without dummy shifts.

  input = "AACGACGTCGACT", k = 4

       i v d n lbl W
       0 0 0 1 CGA C
       1 1 1 1 AAC G
       2 2 0 1 GAC g
       3   0 0 GAC T
       4 3 0 1 GTC G
       5 4 0 1 ACG A
       6   0 0 ACG T
       7 5 0 1 TCG a
       8 6 0 1 ACT $
       9 7 0 1 CGT C
*/

const int sigma    = 4;
string edges       = "CGgTGATa$C";
string node_flags  = "1110110111";
string dummy_flags = "0100000000";
vector<string> dummies = vector<string>({ "AAC" });
std::array<size_t, 1+sigma> symbol_ends{0, 1, 5, 8, 10};

const size_t g_size = 10;
const size_t num_edges = 9; // outgoing dummy edge is not a real edge
const size_t num_nodes = 8;
string expected_node_labels[num_nodes] = {
  "CGA","AAC","GAC","GTC","ACG","TCG","ACT","CGT"
};

const int expected_outdegree[num_nodes] = { 1, 1, 2, 1, 2, 1, 0, 1};
const int  expected_indegree[num_nodes] = { 2, 0, 1, 1, 2, 1, 1, 1};

const ssize_t expected_forward[g_size] = {
     // i v d n lbl W
  2, // 0 0 0 1 CGA C
  5, // 1 1 1 1 AAC G
  5, // 2 2 0 1 GAC g
  8, // 3   0 0 GAC T
  7, // 4 3 0 1 GTC G
  0, // 5 4 0 1 ACG A
  9, // 6   0 0 ACG T
  0, // 7 5 0 1 TCG a
 -1, // 8 6 0 1 ACT $
  7, // 9 7 0 1 CGT C
};

const ssize_t expected_backward[g_size] = {
     // i v d n lbl W
  5, // 0 0 0 1 CGA C
 -1, // 1 1 1 1 AAC G
  0, // 2 2 0 1 GAC g
  0, // 3   0 0 GAC T
  9, // 4 3 0 1 GTC G
  1, // 5 4 0 1 ACG A
  1, // 6   0 0 ACG T
  4, // 7 5 0 1 TCG a
  3, // 8 6 0 1 ACT $
  6, // 9 7 0 1 CGT C
};

/* test the graph by directly calling the debruijn_graph methods */
template <class DBG>
void test_graph_directly(const DBG & db)
{
  for (size_t i = 0;i < num_nodes;i ++) {
    BOOST_CHECK_MESSAGE(db.node_label(i) == expected_node_labels[i],"node_label(" << i << ") failed");
    // TODO : also check node values returned by incoming and outgoing
    BOOST_CHECK_MESSAGE(db._backward(i) == expected_backward[i],"backward(" << i << ") failed");
    BOOST_CHECK_MESSAGE(db._forward(i) == expected_forward[i],"forward(" << i << ") : " << db._forward(i) << " != " << expected_forward[i]);

    // for outdegree we will also count the outgoing paths we find via the outgoing function
    int outdegree = db.outdegree(i);
    int indegree = db.indegree(i);
    BOOST_CHECK_MESSAGE(outdegree == expected_outdegree[i],"outdegree(" << i << ") : " << outdegree << " != " << expected_outdegree[i]);
    BOOST_CHECK_MESSAGE(indegree == expected_indegree[i],"indegree(" << i << ") failed");

    // and count the actual outgoing paths
    for (int x = 0;x <= sigma; x++) {
      if (db.outgoing(i,x) >= 0) outdegree--;
      if (db.incoming(i,x) >= 0) indegree--;
    }

    BOOST_CHECK_MESSAGE(outdegree == 0, "The outgoing call returns different number of paths to outdegree(" << i << ")");
    BOOST_CHECK_MESSAGE(indegree == 0, "The incoming call returns different number of paths to indegree(" << i << ")");
  }
}

/* this BFS visitor collects the names of each edge as it is added to the tree */
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
  
  for (size_t i = 0;i < num_nodes;i ++)
  {
    // the easy way to read the vertex_name_t property
    BOOST_CHECK_MESSAGE(get(vertex_name_t(),db,i) == expected_node_labels[i],"node_label(" << i << ") failed");
    
    int outdegree = out_degree(i, db);
    BOOST_CHECK_MESSAGE(outdegree == expected_outdegree[i],"outdegree(" << i << ") failed");
  }
  
  // prepare the visitor for a BFS
  std::stringstream ss;
  typename property_map<Graph, edge_name_t>::type edge_name_map = get(edge_name_t(), db);
  bfs_edge_collector<typename property_map<Graph , boost::edge_name_t>::type > edge_visitor(edge_name_map,ss);

  // do the BFS
  breadth_first_search(db, 0, boost::visitor(edge_visitor));
  
  BOOST_CHECK_EQUAL(ss.str(), "AACGATCCTG");
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
  wt_huff<rrr_vector<63>> wt;
  construct_im(wt, edges_v);

  vector<kmer_t> dummy_v;
  dummy_v.reserve(dummies.size());
  for (auto dummy : dummies) {
    auto x = string_to_kmer<kmer_t>(dummy);
    dummy_v.push_back(x);
  }

  // construct the graph!
  debruijn_graph<> db(k, node_flags, wt, symbol_ends, alphabet, dummy_bv, dummy_v);

  test_graph_directly(db);

  // Serialization
  std::string temp_file = boost::filesystem::unique_path().native();
  debruijn_graph<> db2;
  sdsl::store_to_file(db, temp_file);
  sdsl::load_from_file(db2, temp_file);
  boost::filesystem::remove(temp_file);
  test_graph_directly(db2);

  test_graph_via_boost(db);
}
