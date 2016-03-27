/*

 Boost Graph Library - Succinct de Bruijn Adapter
 
 Allows you to use BGL algorithms on debruijn_graph<> objects
 
 The graph supports these concepts (see http://www.boost.org/doc/libs/1_60_0/libs/graph/doc/graph_concepts.html)
   IncidenceGraph
   PropertyGraph
   VertexListGraph

 And implements these properties
   vertex_name_t (std::string) - returns the label for a given node in the graph
   vertex_index_t (size_t) - returns the index of the vertex (used by BFS)
   edge_name_t (char) - returns the character associated with the given edge in the graph
 
*/
#ifndef BGL_SDB_ADAPTER_HPP
#define BGL_SDB_ADAPTER_HPP

#include <cassert>
#include <utility>
#include <cstddef>
#include <boost/iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/irange.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>
#include <algorithm>
#include "debruijn_graph.hpp"


/* Because the identifying type of the graph is quite complex, the following def's make things much more concise */
#define SUCCINCT_TEMPLATE template <size_t t_sigma,\
class  t_bit_vector_type,\
class  t_bv_rank_type,\
class  t_bv_select_type,\
class  t_edge_vector_type,\
class  t_symbol_type,\
class  t_label_type>
#define SUCCINCT_TEMPLATE_ARGUMENTS t_sigma,t_bit_vector_type,t_bv_rank_type,t_bv_select_type,t_edge_vector_type,t_symbol_type,t_label_type
#define SUCCINCT_TYPE debruijn_graph< SUCCINCT_TEMPLATE_ARGUMENTS >

/* Install the edge name map in the boost namespace */
namespace boost {

  /*
   
    Define the map that will extract the label (char) associated with an edge in the graph 
   
    Has to be in the same namespace as the get function otherwise the call to boost::get won't find
    the templated map type argument

  */
  SUCCINCT_TEMPLATE
  class sdb_edge_name_map
  {
  public:
    typedef size_t key_type;
    typedef typename t_label_type::value_type value_type;
    typedef typename t_label_type::value_type reference;
    typedef boost::readable_property_map_tag category;
    
    inline sdb_edge_name_map(const SUCCINCT_TYPE & graph) : graph(graph) {}
    value_type get (const key_type i) const
    {
      return graph.edge_symbol(i);
    }
    
  protected:
    const SUCCINCT_TYPE & graph;
  };


  /* The get function to return the label from the map */
  SUCCINCT_TEMPLATE
  typename t_label_type::value_type get(const sdb_edge_name_map<SUCCINCT_TEMPLATE_ARGUMENTS> & map, const size_t i)
  {
    return map.get(i);
  }
  
  /* The get function to return the map for a graph */
  SUCCINCT_TEMPLATE
  sdb_edge_name_map<SUCCINCT_TEMPLATE_ARGUMENTS>
  get(edge_name_t, SUCCINCT_TYPE& graph)
  {
    return sdb_edge_name_map<SUCCINCT_TEMPLATE_ARGUMENTS>(graph);
  }
  
  /* Install the property map type by specialisaing property_map */
  SUCCINCT_TEMPLATE
  struct property_map<SUCCINCT_TYPE, boost::edge_name_t>
  {
    typedef sdb_edge_name_map<SUCCINCT_TEMPLATE_ARGUMENTS> type;
    typedef sdb_edge_name_map<SUCCINCT_TEMPLATE_ARGUMENTS> const_type;
  };
}

/* Install the vertex name map in the boost namespace */
namespace boost {

  /* 
   
    Define the map that will extract the label (string) associated with a vertex in the graph 
    
    Has to be in the same namespace as the get function otherwise the call to boost::get won't find
    the templated map type argument
  */
  SUCCINCT_TEMPLATE
  struct sdb_vertex_name_map
  {
  public:
    typedef size_t key_type;
    typedef t_label_type value_type;
    typedef t_label_type reference;
    typedef boost::readable_property_map_tag category;
    
    inline sdb_vertex_name_map(const SUCCINCT_TYPE & graph) : graph(graph) {}
    value_type get (const key_type i) const
    {
      return graph.node_label(i);
    }
    
  protected:
    const SUCCINCT_TYPE & graph;
  };

  /* PropertyMap concept for vertex name */
  SUCCINCT_TEMPLATE
  t_label_type get(const sdb_vertex_name_map<SUCCINCT_TEMPLATE_ARGUMENTS> & map, const size_t i)
  {
    return map.get(i);
  }

  /* ReadablePropertyGraph concepts */
  
  /* The get function to return the map for a graph */
  SUCCINCT_TEMPLATE
  sdb_vertex_name_map<SUCCINCT_TEMPLATE_ARGUMENTS>
  get(vertex_name_t, const SUCCINCT_TYPE& graph)
  {
    return sdb_vertex_name_map<SUCCINCT_TEMPLATE_ARGUMENTS>(graph);
  }
  
  SUCCINCT_TEMPLATE
  t_label_type
  get(vertex_name_t, const SUCCINCT_TYPE& graph, const size_t i)
  {
    return graph.node_label(i);
  }
  
  /* Install the property map type by specialisaing property_map */
  SUCCINCT_TEMPLATE
  struct property_map<SUCCINCT_TYPE, boost::vertex_name_t>
  {
    typedef sdb_vertex_name_map<SUCCINCT_TEMPLATE_ARGUMENTS> type;
    typedef sdb_vertex_name_map<SUCCINCT_TEMPLATE_ARGUMENTS> const_type;
  };
  SUCCINCT_TEMPLATE
  struct property_map<const SUCCINCT_TYPE, boost::vertex_name_t>
  {
    typedef sdb_vertex_name_map<SUCCINCT_TEMPLATE_ARGUMENTS> type;
    typedef sdb_vertex_name_map<SUCCINCT_TEMPLATE_ARGUMENTS> const_type;
  };}


/* 
  
   Define the types for the out edge iterator 
 
*/
SUCCINCT_TEMPLATE
struct succinct_out_edge_iterator
{
  typedef boost::counting_iterator<size_t> type;
};

SUCCINCT_TEMPLATE
struct succinct_out_edge_range
{
  typedef typename succinct_out_edge_iterator<SUCCINCT_TEMPLATE_ARGUMENTS>::type Iter;
  typedef std::pair<Iter,Iter> type;
};

/* 
 
   Install the graph in the boost namespace 
 
*/
namespace boost {
  
  /* The graph traits... */
  SUCCINCT_TEMPLATE
  struct graph_traits< SUCCINCT_TYPE >
  {
    typedef SUCCINCT_TYPE array_type;
    
    typedef size_t vertex_descriptor;
    typedef size_t edge_descriptor;
    typedef directed_tag directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category;
    
    // Indicate which 'models' are supported by the graph
    struct traversal_category : public virtual incidence_graph_tag, public virtual vertex_list_graph_tag {};
    
    // IncidenceGraph model
    typedef typename succinct_out_edge_iterator<SUCCINCT_TEMPLATE_ARGUMENTS>::type out_edge_iterator;
    typedef size_t degree_size_type;

    // VertexListGraph model (seems to be assumed by breadth_first_search even though it's not checked)
    typedef size_t vertices_size_type;
    typedef boost::counting_iterator<size_t> vertex_iterator;

    static size_t null_vertex() {return size_t(-1);}
  };
  
}

/* 
 
   IncidenceGraph concept methods

*/
namespace boost {
  
  /* iterate the out edges of a node */
  SUCCINCT_TEMPLATE
  typename succinct_out_edge_range<SUCCINCT_TEMPLATE_ARGUMENTS>::type //graph_traits< SUCCINCT_TYPE >::out_edge_iterator //detail::val_out_edge_ret<EdgeList>::type
  out_edges(typename graph_traits< SUCCINCT_TYPE >::vertex_descriptor v, //aka size_t
            const SUCCINCT_TYPE & g)
  {
    typedef typename succinct_out_edge_iterator<SUCCINCT_TEMPLATE_ARGUMENTS>::type Iter;
    typedef typename succinct_out_edge_range<SUCCINCT_TEMPLATE_ARGUMENTS>::type return_type;
    std::pair<size_t, size_t> rng = g._node_range(v);
    
    // Skip terminator edge
    // The iterator needs to skip all terminators, not just the first one if it happens to be the terminator.
    // It depends on the way the graph is built.  If we never have a mixture of terminators and non-terminators on a node,
    // this shouldn't be an issue.
    if (g.edge_symbol(rng.first) == g.m_alphabet[0]) rng.first++;
    
    return return_type(Iter(rng.first),Iter(rng.second+1));
  }
  
  /* return the out degree of a node */
  SUCCINCT_TEMPLATE
  size_t
  out_degree(size_t v, const SUCCINCT_TYPE & g)
  {
    return g.outdegree(v);
  }
  
  /* return the source vertex for a given edge */
  SUCCINCT_TEMPLATE
  size_t source(size_t e, const SUCCINCT_TYPE & g)
  {
    return g._edge_to_node(e);
  }
  
  /* return the target vertex for a given edge */
  SUCCINCT_TEMPLATE
  size_t target(size_t e, const SUCCINCT_TYPE & g)
  {
    size_t next = g._forward(e);
    //if (next == SUCCINCT_TYPE::npos) next = 0; //
    return g._edge_to_node(next);
  }
  

/* 
 
  VertexListGraph concept methods
 
*/
  
  /* Return the list of vertices */
  SUCCINCT_TEMPLATE
  std::pair<boost::counting_iterator<size_t>,
    boost::counting_iterator<size_t> >
  vertices(const SUCCINCT_TYPE & g)
  {
    typedef boost::counting_iterator<size_t> Iter;
    return std::make_pair(Iter(0), Iter(g.num_nodes()));
  }
  
  /* Count the vertices */
  SUCCINCT_TEMPLATE
  size_t
  num_vertices(const SUCCINCT_TYPE & g)
  {
    return g.num_nodes();
  }

  /* Define the vertex index property map to just return identity */
  SUCCINCT_TEMPLATE
  struct property_map<SUCCINCT_TYPE, vertex_index_t>
  {
    typedef identity_property_map type;
    typedef type const_type;
  };
  
/*
 
  vertex index property
 
*/

  /* BFS also assumes a vertex_index_t property exists for vertices */
  SUCCINCT_TEMPLATE
  identity_property_map
  get(vertex_index_t, const SUCCINCT_TYPE&)
  {
    /* Fortunately the vertex descriptor is the index so we can just use the identity map */
    return identity_property_map();
  }

} // namespace boost


#endif