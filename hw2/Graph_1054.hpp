#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <functional>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E = char>
class Graph {

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;
  
  using node_idx_type = size_type;
  using edge_idx_type = size_type;
  
  using node_value_type = V;
  using edge_value_type = E;

 private:
 
  using node_id_type = size_type;
  using edge_id_type = size_type;
 
  /** Invalid node and edge id value */
  static node_id_type constexpr invalid_node_id = node_id_type(-1);
  static edge_id_type constexpr invalid_edge_id = edge_id_type(-1);

 public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
  {
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
   public:
    
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node()
    : gp(nullptr),
      nid(invalid_node_id)
    {
    }
    
    /** Return this node's index, a number in the range [0, graph_size). */
    node_idx_type index() const
    {
      return gp->all_nodes_ind[nid];
    }

    /** Return const reference to this node's position. */
    const Point & position() const
    {
      assert(is_valid());
      return gp->npos[index()];
    }
    
    /** Return reference to this node's position. */
    Point & position()
    {
      assert(is_valid());
      return gp->npos[index()];
    }
    
    /** Return this node's value */
    node_value_type & value()
    {
      assert(is_valid());
      return gp->all_nodes[index()];
    }
    
    /** Return this node's value as const */
    const node_value_type & value() const
    {
      assert(is_valid());
      return gp->all_nodes[index()];
    }
    
    /* Return the number of incident edges of this node*/
    size_type degree() const
    {
      assert(is_valid());
      return gp->nodes_incd[index()].size();
    }
    
    
    /* Returns an iterator to the beginning of incident edges*/
    incident_iterator edge_begin() const
    {
      assert(is_valid());
      return incident_iterator(gp, nid, gp->nodes_incd[index()].begin());
    }
    
    /* Returns an iterator to the end of incident edges*/
    incident_iterator edge_end() const
    {
      assert(is_valid());
      return incident_iterator(gp, nid, gp->nodes_incd[index()].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node & n) const {
      // HW0: YOUR CODE HERE
      //(void) n;          // Quiet compiler warning
      if(gp == n.gp && nid == n.nid)
      {
        return true; 
      }
      return false; 
    }
    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node & n) const {
      // HW0: YOUR CODE HERE
      //(void) n;           // Quiet compiler warning
      if(gp < n.gp)
      {
        return true; 
      }
      if(gp == n.gp && nid < n.nid)
      {
        return true; 
      }
      return false; 
    }

   private:
    
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    /** Check node for validity.
     * @return @p true if node is valid (i.e. represents an actual graph node)
     * or @p false otherwise.
     */ 
    bool is_valid() const
    {
      return gp != nullptr &&
             nid != invalid_node_id &&
             nid < gp->all_nodes_ind.size() &&
             index() < gp->num_nodes() &&
             gp->nodes_ids[index()] == nid;
    }
    
    /** Convenience accessor for id **/
    node_id_type id() const
    {
      return nid;
    }
    
    //Node Constructor
    Node(const Graph * curr_gp, const node_id_type curr_nind){
      gp = const_cast<Graph*> (curr_gp);
      nid = curr_nind;
    }

    Graph * gp;
    node_id_type nid;
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return num_nodes();
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return nodes_ids.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point & position, node_value_type value = node_value_type{})
  {
    const node_id_type id = all_nodes_ind.size();
    const node_idx_type idx = nodes_ids.size();
    
    all_nodes_ind.push_back(idx);
    nodes_ids.push_back(id);
    npos.push_back(position);
    all_nodes.emplace_back(std::move(value)); 
    nodes_incd.push_back({});
    
    return Node(this, id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //(void) n;            // Quiet compiler warning
    if(n.gp ==this && n.is_valid())
    {
      return true; 
    }
    return false;   
  }


  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(node_idx_type i) const
  {
    return Node(const_cast<Graph *>(this), nodes_ids[i]);
  }
  
  /** Remove a node from the graph and all incident edges **/
  size_type remove_node(const Node & n)
  {    
    // remove edges
    for(auto it = n.edge_begin(); it != n.edge_end();)
    {
      remove_edge(*it);
    }
    
    const node_idx_type nidx = n.index();
    
    // pop removed node data
    using std::swap; 
    swap(all_nodes[nidx], all_nodes.back()); all_nodes.pop_back();
    std::swap(npos[nidx], npos.back()); npos.pop_back();
    std::swap(nodes_ids[nidx], nodes_ids.back()); nodes_ids.pop_back();
    std::swap(nodes_incd[nidx], nodes_incd.back()); nodes_incd.pop_back();
    all_nodes_ind[nodes_ids[nidx]] = nidx;

    return 1;
  }
  
  /** Remove a node from the graph and all incident edges*/
  node_iterator remove_node(node_iterator n_it)
  {
    remove_node(*n_it);
    return n_it; 
  }


  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */    
   Edge(){
      gp = nullptr;
      nid = edge_id_type(-1); 
      nid1 = node_id_type(-1); 
      nid2 = node_id_type(-1); 
    }
    /** Return index of the edge **/
    edge_idx_type index() const
    {
      return gp->all_edge_ind[nid];
    }

    /** Return the first node of the Edge */
    Node node1() const
    {
      return Node(gp, nid1);
    }

    /** Return the second node of the Edge */
    Node node2() const
    {
      return Node(gp, nid2);
    }
    
    /** Returns length of a valid edge **/
    double length() const
    {
      assert(is_valid());
      return norm_2(node1().position() - node2().position());
    }
    
    /** Return const reference to edge's value **/
    const edge_value_type & value() const
    {
      assert(is_valid());
      return gp->all_edges[index()];
    }
    
    /** Return reference to edge's value **/
    edge_value_type & value()
    {
      assert(is_valid());
      return gp->all_edges[index()];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge & e) const
    {
      return (gp == e.gp && nid == e.nid);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge & e) const
    {
      if(gp < e.gp)
      {
        return true;
      }
      if(gp == e.gp)
      {
        if(nid < e.nid)
        {
          return true;
        }
      }
      return false;
    }

   private:
    
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    /** Check if the edge is valid.
     * @return @p true if edge is valid (i.e. represents an actual graph edge)
     * or @p false otherwise.
     */ 
    bool is_valid() const
    {
      return gp != nullptr &&
             nid != invalid_edge_id &&
             nid < gp->all_edge_ind.size() &&
             index() < gp->num_edges() &&
             gp->edges_ids[index()] == nid &&
             nid1 != invalid_node_id && 
             nid2 != invalid_node_id;
    }
    
    //Constructor for Edge
    Edge(const Graph * curr_gp, const edge_id_type curr_nid,
         const node_id_type curr_nid1, 
         const node_id_type curr_nid2)
    {
      gp = const_cast<graph_type*>(curr_gp);
      nid = curr_nid; 
      nid1 = curr_nid1;
      nid2 = curr_nid2;
    }

    Edge(Graph * const graph,
         const edge_id_type id)
    : Edge(graph, id, 
           graph->all_edges_nodes[graph->all_edge_ind[id]].first,
           graph->all_edges_nodes[graph->all_edge_ind[id]].second)
    {}
    
    Graph * gp;
    edge_id_type nid;
    node_id_type nid1;
    node_id_type nid2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return edges_ids.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    return Edge(const_cast<Graph *>(this), edges_ids[i], 
                all_edges_nodes[i].first, all_edges_nodes[i].second);
  }


  /** Find an edge between nodes a and b. */
  std::pair<bool, std::vector<edge_id_type>::iterator> 
  edge_lookup(const node_id_type a, const node_id_type b)
  {
    auto & edges = nodes_incd[all_nodes_ind[a]]; 
    auto it = std::find_if(edges.begin(), edges.end(),[&](edge_id_type const & euid)
       { 
         auto const & nodes = all_edges_nodes[all_edge_ind[euid]];
         return ((nodes.first == a && nodes.second == b) || (nodes.first == b && nodes.second == a));
       });
    return std::make_pair(it != edges.end(), it);
  }
  
  std::pair<bool, std::vector<edge_id_type>::const_iterator> 
  edge_lookup(const node_id_type a, const node_id_type b) const
  {
    auto const lookup = const_cast<Graph *>(this)->edge_lookup(a, b);
    return std::pair<bool, std::vector<edge_id_type>::const_iterator>(lookup.first, lookup.second);
  }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node & a, const Node & b) const
  {
    return edge_lookup(a.id(), b.id()).first;
  }
  
  /** Test whether a given edge exists in the graph.
   * This is pretty much the same as the edge is being valid.
   * @param[in] e the edge to test
   */
  bool has_edge(const Edge & e) const
  {
    return (e.is_valid() && e.gp == this);
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node & a, const Node & b, edge_value_type value = edge_value_type{})
  {    
    const auto lookup = edge_lookup(a.id(), b.id());
    const edge_id_type euid = lookup.first ? *lookup.second : all_edge_ind.size();
    
    if (!lookup.first)
    {
      all_edge_ind.push_back(edges_ids.size());
      edges_ids.push_back(euid);
      all_edges_nodes.emplace_back(a.id(), b.id());
      all_edges.push_back(value);
      nodes_incd[a.id()].push_back(euid);
      nodes_incd[b.id()].push_back(euid);
    }
    
    return Edge(this, euid, a.id(), b.id());
  }
  
  /** Remove the edge connecting two nodes, if one exists.   */
  size_type remove_edge(const Node & a, const Node & b)
  {
    assert(a.gp == this && a.is_valid());
    assert(b.gp == this && b.is_valid());
    
    if (!has_edge(a, b))
    {
      return 0;
    }
    
    // check if edge exists
    const auto it1 = edge_lookup(a.id(), b.id()).second;
    const auto it2 = edge_lookup(b.id(), a.id()).second;
    assert(*it1 == *it2); // sanity check
    
    const edge_id_type euid = *it1;
    const edge_idx_type eidx = all_edge_ind[euid];
    
    std::swap(*it1, nodes_incd[a.index()].back()); nodes_incd[a.index()].pop_back();
    std::swap(*it2, nodes_incd[b.index()].back()); nodes_incd[b.index()].pop_back();
    
    // pop removed edge data
    using std::swap; // ADL lookup in case value type has a custom swap function
    swap(all_edges[eidx], all_edges.back()); all_edges.pop_back();
    std::swap(all_edges_nodes[eidx], all_edges_nodes.back()); all_edges_nodes.pop_back();
    std::swap(edges_ids[eidx], edges_ids.back()); edges_ids.pop_back();
    all_edge_ind[edges_ids[eidx]] = eidx;
    
    return 1;
  }
  
  /** Remove an edge from the graph */
  size_type remove_edge(const Edge & e)
  {
    assert(has_edge(e));
    return remove_edge(e.node1(), e.node2());
  }
  
  /** Remove an edge from the graph.
   * 
   * @param[in] e_it iterator to the edge to remove
   * @return iterator pointing to the next edge after removed one,
   *         or past-the-end if removed edge was last (by index)
   *
   * @pre @a e_it is a valid edge iterator for this graph.
   * @post same as remove_edge(*e_it)
   * 
   * Complexity and behavior same as remove_edge(*e_it).
   */
  edge_iterator remove_edge(edge_iterator e_it)
  {
    remove_edge(*e_it);
    return e_it; // pointing to next edge now
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    all_nodes_ind.clear();
    nodes_ids.clear();
    npos.clear();
    all_nodes.clear();
    nodes_incd.clear();
    all_edge_ind.clear();
    edges_ids.clear();
    all_edges_nodes.clear();
    all_edges.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
    {
    }
    
    /** Dereference operator. 
     * @return a proxy object representing pointed-to node
     * @pre iterator is valid (not default-constructed or past-the-end)
     **/
    value_type operator*() const
    {
      assert(is_valid());
      return Node(gp, *niter);
    }
    
    /** Prefix increment operator.
     * @return reference to this iterator
     * @pre iterator is valid (not default-constructed or past-the-end)
     * @post iterator is pointing to the next node in graph's node collection,
     *       or past-the-end if previously pointed-to node was last in graph
     **/
    NodeIterator & operator++()
    {
      assert(is_valid());
      ++niter;
      return *this;
    }
    
    /** Equality comparison operator. 
     * @return @p true if iterators are identical, i.e. pointing to
     *         the same graph object and same node in that graph,
     *         or both are default-constructed, or @p false otherwise
     **/
    bool operator==(const NodeIterator & other) const
    {
      return gp == other.gp && niter == other.niter;
    }

   private:
   
    using iterator = std::vector<node_id_type>::const_iterator;
   
    friend class Graph;
    
    /** Check whether the iterator is valid (not default-constructed) **/
    bool is_valid() const
    {
      return gp != nullptr;
    }
    
    /** Construct a valid NodeIterator. Used privately by Graph. */
    NodeIterator(Graph * const curr_gp, const iterator curr_iter)
    {
      gp = const_cast<graph_type*>(curr_gp);
      niter = curr_iter; 
    }
    
    Graph * gp;
    iterator niter;
    
  };

  /** Node begin iterator */
  node_iterator node_begin() const
  {
    return node_iterator(const_cast<Graph *>(this), nodes_ids.begin());
  }
  
  /** Node end iterator */
  node_iterator node_end() const
  {
    return node_iterator(const_cast<Graph *>(this), nodes_ids.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator()
    {
    }

    /** Dereference operator */
    Edge operator*() const
    {
      const auto & nodes = gp->all_edges_nodes[gp->all_edge_ind[*niter]];
      const node_id_type dst_node = (c_node == nodes.first) ? nodes.second : nodes.first;
      return Edge(gp, *niter, c_node, dst_node);
    }
    
    /** Prefix increment operator */
    IncidentIterator & operator++()
    {
      ++niter;
      return *this;
    }
    
    /** Equality comparison operator **/
    bool operator==(const IncidentIterator & other) const
    {
      return ( niter == other.niter && gp == other.gp && c_node == other.c_node);
    }

   private:
   
    friend class Graph;
    
    using iterator = std::vector<edge_id_type>::const_iterator;
    
    /** Check whether the iterator is valid**/
    bool is_valid() const
    {
      return gp != nullptr && c_node != invalid_node_id;
    }
    
    IncidentIterator(Graph * const curr_gp, const node_id_type curr_node,
      const iterator curr_iter)
    {
      gp = const_cast<graph_type*>(curr_gp);
      niter = curr_iter; 
      c_node = curr_node; 
    }
    
    Graph * gp;
    iterator niter;
    node_id_type c_node;  
    
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() 
    {
    }
    
    /** Dereference operator
     * @return a proxy Edge object representing the pointed-to graph edge
     * @pre the iterator is valid (not past-the-end or default-constructed)
     */
    Edge operator*() const
    {
      return Edge(gp, *niter);
    }
    
    /** Prefix increment operator
     * @return reference to this iterator
     * @pre the iterator is valid (not past-the-end or default-constructed)
     * @post iterator is pointing to the next edge in the edge collection,
     *       or past-the-end if currently pointed-to edge is the last one
     */
    EdgeIterator & operator++()
    {
      ++niter;
      return *this;
    }
    
    /** Equality comparison operator
     * @return @p true if iterators are identical, i.e. pointing to
     *         the same graph object and same edge in that graph,
     *         or both are default-constructed, or @p false otherwise
     */
    bool operator==(const EdgeIterator & other) const
    {
      return (niter == other.niter && gp == other.gp);
    }

   private:
   
    friend class Graph;
    
    using iterator = std::vector<edge_id_type>::const_iterator;
    
    bool is_valid() const
    {
      return gp != nullptr;
    }
    
    
    EdgeIterator(Graph * const curr_gp, const iterator curr_iter)
    {
      gp = const_cast<graph_type*>(curr_gp);
      niter = curr_iter; 
    }

    Graph * gp;
    iterator niter;
  };
  
  /** Edge begin iterator **/
  edge_iterator edge_begin() const
  {
    return edge_iterator(const_cast<Graph *>(this), edges_ids.begin());
  }
  
  /** Edge end iterator **/
  edge_iterator edge_end() const
  {
    return edge_iterator(const_cast<Graph *>(this), edges_ids.end());
  }

 private:
           
  std::vector<Point> npos;
  std::vector<node_value_type> all_nodes;
  std::vector<node_id_type> nodes_ids;
  std::vector<node_idx_type> all_nodes_ind;
  std::vector<std::vector<edge_id_type>> nodes_incd;
  std::vector<edge_id_type> edges_ids;
  std::vector<edge_idx_type> all_edge_ind;
  std::vector<std::pair<node_id_type, node_id_type>> all_edges_nodes;
  std::vector<edge_value_type> all_edges;

};

#endif // CME212_GRAPH_HPP
