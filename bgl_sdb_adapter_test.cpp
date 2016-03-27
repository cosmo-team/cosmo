/*

  Thoroughly tests the debruijn graph for the graph in the paper
  http://alexbowe.com/succinct-debruijn-graphs/
 
  Does it both directly via calls to the debruijn_graph methods, and also via the boost graph library interfaces
 
*/

// this define tweaks the debruijn_graph _forward(i) function to work when W[i] is a flagged edge
#define FIX_FORWARD

// this define enables another change I had to make to stop outgoing() from crashing in certain test cases
#define FIX_OUTGOING

#include "bgl_sdb_adapter.hpp"
#include <boost/concept/assert.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/breadth_first_search.hpp>
/* the following #define tells the boost::test framework to generate a main() for this module. only one cpp file in the test suite should have it */
#define BOOST_TEST_MODULE BGL SDB Adpater Test
#include <boost/test/included/unit_test.hpp>

/*

  The data for the graph is summarised in the matrix below
 
       i i+1 L L' v label  W  W'   L*
       0  1  1 1  0  $$$   T  T+   0
       1  2  1 1  1  CGA   C  C+   1
       2  3  1 1  2  $TA   C  C+   0
       3  4  0 1  3  GAC   G  G   (3)
       4  5  1 0  3  GAC   T  T+   2
       5  6  1 1  4  TAC   G- G+   1
       6  7  1 1  5  GTC   G  G+   0
       7  8  0 1  6  ACG   A  A   (3)
       8  9  1 0  6  ACG   T  T+   2
       9  10 1 1  7  TCG   A- A+   0
       10 11 1 1  8  $$T   A  A+   1
       11 12 1 1  9  ACT   $  $+   1
       12 13 1 1  10 CGT   C  C+
 
*/
const char * symbols =    "TCCGTGGATAA$C";
#ifdef FIX_FORWARD
const char * end_flag =   "1111101110111"; // W flagging indicates this is the first edge that has this symbol & target node (i.e. not -)
#else
// my experiment to try and find an alternative end_flag interpretation that gets the graph to work without the FIX_FORWARD...
// this would make sense from the way that _forward assumes the non-flagged edge comes _last_, but it breaks lots of other things.
const char * end_flag =   "1110111011111"; // W'+ flagging indicates this is the last edge that has this symbol & target node
#endif
const char * start_flag = "1111011101111"; // L' = is this the first edge of the node (unlike L in the paper which is last edge)
const int sigma = 4;

/* The test data that provides the various attributes of the graph that we test for */
const size_t num_nodes = 11;
const char * expected_node_labels[num_nodes] = {
  "$$$","CGA","$TA","GAC","TAC","GTC","ACG","TCG","$$T","ACT","CGT"
};
const ssize_t expected_forward[13] =
{
  10, // 0 $$$T -> $$TA
  3,  // 1 CGAC -> GACT
  5,  // 2 $TAC -> TACG-
  7,  // 3 GACG -> ACGT
  11, // 4 GACT -> ACT$
  
  7,  // 5 TACG- -> ACGT
  9,  // 6 GTCG -> TCGA-
  1,  // 7 ACGA -> CGAC
  12, // 8 ACGT -> CGTC
  1,  // 9 TCGA -> CGAC
  2,  //10 $$TA -> $TAC
  (ssize_t)-1,//11 ACT$ -> -1
  6,  //12 CGTC -> GTCG
};
const size_t expected_backward[13] =
{
  0,  // 0 $$$T -> (invalid)
  7,  // 1 CGAC -> ACG(A)
  10, // 2 $TAC -> $$T
  1,  // 3 GACG -> CGA
  1,  // 4 GACT -> CGA
  
  2,  // 5 TACG -> $TA
  12, // 6 GTCG -> CGT
  3,  // 7 ACGA -> GAC
  3,  // 8 ACGT -> GAC
  6,  // 9 TCGA -> GTC
  0,  //10 $$TA -> $$$
  3,  //11 ACT$ -> GAC
  7,  //12 CGTC -> ACG
};

// indegree
// look at the graph in http://alexbowe.com/succinct-debruijn-graphs/
// everything except v=0 $$$ has 1 incoming edge, except
// CGA (1) has 2
// ACG (6) has 2
const int expected_indegree[11] = { 0, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1 };

// outdegree
// GAC (3) has 2
// ACG (6) has 2
// ACT (9) has 0
// everything else 1
const int expected_outdegree[11] = { 1, 1, 1, 2, 1, 1, 2, 1, 1, 0, 1 };

/* test the graph by directly calling the debruijn_graph methods */
void test_graph_directly(debruijn_graph<> & db)
{
  for (size_t i = 0;i < num_nodes;i ++)
  {
    BOOST_CHECK_MESSAGE(db.node_label(i) == expected_node_labels[i],"node_label(" << i << ") failed");
    BOOST_CHECK_MESSAGE(db._backward(i) == expected_backward[i],"backward(" << i << ") failed");
    BOOST_CHECK_MESSAGE(db._forward(i) == expected_forward[i],"forward(" << i << ") failed");
    
    // for outdegree we will also count the outgoing paths we find via the outgoing function
    int outdegree = db.outdegree(i);
    int indegree = db.indegree(i);
    BOOST_CHECK_MESSAGE(outdegree == expected_outdegree[i],"outdegree(" << i << ") failed");
    BOOST_CHECK_MESSAGE(indegree == expected_indegree[i],"indegree(" << i << ") failed");
    
    // and count the actual outgoing paths
    for (int x = 0;x <= sigma; x++)
    {
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
    ss << boost::get(m_name_map, e);
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
  
  BOOST_CHECK_EQUAL(ss.str(), "TACGATCCTG");
}

/* builds the graph and calls the test functions */
BOOST_AUTO_TEST_CASE ( create_graph_from_packed_edges )
{
  std::string alphabet = "$ACGT";
  std::fstream packed_edges("test.packed",ios_base::in | ios_base::out | ios_base::binary | ios_base::trunc);;
  
  // write the packed edges
  uint64_t block = 0;
  size_t len   = 0;
  for (int i = 0;i < 13;i ++)
  {
    // io.hpp pack_edge defines the end_flag as meaning 'does not have a -' and start_flag === L ...
    append_packed_edge(block, pack_edge(alphabet.find(symbols[i]),start_flag[i]=='1',end_flag[i]=='1'));
    if (++len == PACKED_CAPACITY)
    {
      block <<= ((PACKED_CAPACITY - len) * PACKED_WIDTH);
      packed_edges.write((char*)&block, sizeof(uint64_t));
      block = 0;
      len = 0;
    }
  }
  if (len)
  {
    block <<= ((PACKED_CAPACITY - len) * PACKED_WIDTH);
    packed_edges.write((char*)&block, sizeof(uint64_t));
  }
  
  // write footer
  std::array<size_t,1+sigma> counts {{1,3,7,10,13}};
  uint64_t k = 4;
  packed_edges.write((char*)&counts[0], (sigma+1) * sizeof(uint64_t));
  packed_edges.write((char*)&k, sizeof(uint64_t));
  packed_edges.flush();
  
  // create graph
  debruijn_graph<> db = debruijn_graph<>::load_from_packed_edges(packed_edges,alphabet);
 
  test_graph_directly(db);
  test_graph_via_boost(db);
}