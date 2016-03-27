/*

  Thoroughly tests the debruijn graph for the graph in the paper
  http://alexbowe.com/succinct-debruijn-graphs/
 
  Does it both directly via calls to the debruijn_graph methods, and also via the boost graph library interfaces
 
*/
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
 
                                $=0,A=2,C=4,G=6,T=8, flag=1
       i i+1 L v label W    L*  encoded edge
       0  1  1 0  $$$   T    0   8
       1  2  1 1  CGA   C    1   4
       2  3  1 2  $TA   C    0   4
       3  4  0 3  GAC   G   (3)  6
       4  5  1 3  GAC   T    2   8
       5  6  1 4  TAC   G-   1   7
       6  7  1 5  GTC   G    0   6
       7  8  0 6  ACG   A   (3)  2
       8  9  1 6  ACG   T    2   8
       9  10 1 7  TCG   A-   0   3
       10 11 1 8  $$T   A    1   2
       11 12 1 9  ACT   $    1   0
       12 13 1 10 CGT   C        4
 
*/
const char * symbols =    "TCCGTGGATAA$C";
const char * end_flag =   "1111101110111";
const char * start_flag = "1111011101111";
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
  
  // prepare the property maps for extracting node labels (strings) and edge characters (chars)
  typename property_map<Graph, edge_name_t>::type edge_name_map = get(edge_name_t(), db);
  typename property_map<Graph, vertex_name_t>::type vertex_name_map = get(vertex_name_t(), db);
  
  for (size_t i = 0;i < num_nodes;i ++)
  {
    BOOST_CHECK_MESSAGE(get(vertex_name_map,i) == expected_node_labels[i],"node_label(" << i << ") failed");
    
    int outdegree = out_degree(i, db);
    BOOST_CHECK_MESSAGE(outdegree == expected_outdegree[i],"outdegree(" << i << ") failed");
  }
  
  // and do a BFS
  std::stringstream ss;
  bfs_edge_collector< boost::property_map<debruijn_graph<> , boost::edge_name_t>::type > edge_visitor(edge_name_map,ss);
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