#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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

  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
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
    Node() {
      // HW0: YOUR CODE HERE
      gp = nullptr;
      nind = size_type(-1);
    }

    /** Return this node's position. */
    const Point & position() const {
      // HW0: YOUR CODE HERE
      return gp->all_nodes[nind];//Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return nind; //size_type(-1);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();

    node_value_type & value()
    {
      return gp->nval[nind];
    }

    // const node_value_type& value() const;
    const node_value_type & value() const
    {
      return gp->nval[nind];
    }

    // size_type degree() const;
    size_type degree() const
    {
      return gp->all_nodes_list[nind].size();
    }


    // incident_iterator edge_begin() const;
    incident_iterator edge_begin() const
    {
      return IncidentIterator(gp, 0, nind);
    }


    // incident_iterator edge_end() const;
    incident_iterator edge_end() const
    {
      return IncidentIterator(gp, this->degree(), nind);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //(void) n;          // Quiet compiler warning
      if(gp == n.gp && nind == n.nind)
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
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE
      //(void) n;           // Quiet compiler warning
      if(gp < n.gp)
      {
        return true; 
      }
      if(gp == n.gp && nind < n.nind)
      {
        return true; 
      }
      return false; 
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    //Node Constructor
    Node(const graph_type* curr_gp, size_type curr_nind){
      gp = const_cast<graph_type*> (curr_gp);
      nind = curr_nind;
    }

    graph_type* gp;
    size_type nind;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return all_nodes.size();;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    //(void) position;      // Quiet compiler warning
    Node new_node;
    new_node.gp = this;
    new_node.nind = num_all_nodes;
    all_nodes.push_back(position);
    num_all_nodes++; 
    return new_node; 
    //return Node();        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    if(n.gp ==this && n.index() < num_all_nodes)
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
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    return Node(this, i);        // Invalid node
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      gp = nullptr;
      nind1 = size_type(-1);
      nind2 = size_type(-1);
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(gp, nind1);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(gp, nind2);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return gp == e.gp && nind1 == e.nind1 && nind2 == e.nind2;
      //return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      if(gp < e.gp)
      {
        return true;
      }
      if(gp == e.gp)
      {
        if(nind1 < e.nind1)
        {
          return true;
        }
        if( (nind1 == e.nind1) and (nind2 < e.nind2))
        {
          return true;
        }
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type*  gp;
    size_type nind1;
    size_type nind2; 

    //Constructor for Edge

    Edge(const graph_type* curr_gp, size_type n1, size_type n2)
    {
      gp = const_cast<graph_type*>(curr_gp);
      nind1 = n1;
      nind2 = n2;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return all_edges.size(); 
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning

    return Edge(this, all_edges[i].first, all_edges[i].second);  
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //(void) a; (void) b;   // Quiet compiler warning
    //curr_n = all_edges.size(); 
    for (unsigned i=0; i < num_all_edges; i++)
    {
      if( all_edges[i].first == a.index() && all_edges[i].second == b.index())
      {
        return true;
      }
      if(all_edges[i].first== b.index() && all_edges[i].second == a.index())
      {
        return true;
      }
            
    }
    return false;
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    //(void) a, (void) b;   // Quiet compiler warning
    //Checking if the edge already exists
    //curr_n = all_edges.size(); 
    const auto new_node = get_node_indices(a, b);
    for (unsigned i=0; i < num_all_edges; i++)
    {
      if( all_edges[i].first == a.index() && all_edges[i].second == b.index())
      {
        return Edge(this, new_node.first, new_node.second);      
      }
      if(all_edges[i].first== b.index() && all_edges[i].second == a.index())
      {
        return Edge(this,  new_node.second, new_node.first);      
      }
            
    }
    // If the edge does not exist already
    
    all_edges.push_back(new_node);
    return Edge(this, new_node.first, new_node.second);      
    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    all_edges.clear();
    all_nodes.clear();
    all_nodes_list.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      gp = nullptr;
      nind = size_type(-1);
    }

    NodeIterator(const graph_type* c_graph, const size_type c_ind):
     gp(const_cast<graph_type*>(c_graph)), nind(c_ind) {};
    // :gp(graph),
    //   nind(c_ind)
    // {
    // }
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    Node operator*() const {
      return Node(gp, nind); 
    }

    // NodeIterator& operator++()

    NodeIterator& operator++() {
      nind ++; 
      return *this;
    }

    // bool operator==(const NodeIterator&) const
    bool operator==(const NodeIterator& niter) const {
      return (nind == niter.nind && gp == niter.gp);
    }


   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type * gp;
    size_type nind;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  // node_iterator node_end() const
  node_iterator node_end() const {
      return node_iterator(this, num_all_nodes);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    IncidentIterator()
    : gp(nullptr),
      inciter{}
    {
      //check this
    }

    IncidentIterator(const graph_type* c_gp, size_type c_ed, size_type c_nd):
        gp(const_cast<graph_type*>(c_gp)), inciter(c_ed, node_id(c_nd){
    }



    Edge operator*() const
    {
      const auto & nodes = gp->all_edges[*inciter];
      const size_type cur_node = (list_node == nodes.first) ? nodes.second : nodes.first;
      return Edge(gp, list_node, cur_node);
    }
    

    // // HW1 #3: YOUR CODE HERE
    // // Supply definitions AND SPECIFICATIONS for:
    // IncidentIterator(const graph_type* curr_gp, size_type curr_edge, 
    //     size_type curr_n):
    //     graph_ptr(const_cast<graph_type*>(curr_gp)), inciter(curr_edge), node_id(curr_n) {
    // }   

    // // Edge operator*() const
    // Edge operator*() const {
    //     return Edge(gp, nind, gp->adj_list[nind][inciter].first,
    //                 gp->adj_list[nind][inciter].second);
    // }
    
    // IncidentIterator& operator++()
    IncidentIterator& operator++() {
        inciter++;
        return *this;
    }

    // bool operator==(const IncidentIterator&) const
    bool operator==(const IncidentIterator& iter) const {
        return (gp == iter.gp && inciter == iter.inciter
                && nind == iter.nind);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    using edge_list_iterator = std::vector<size_type>::const_iterator;
    IncidentIterator(graph_type* const graph, 
                     size_type lnode,
                     size_type iter)
    : gp(graph),
      list_node(lnode),
      inciter(iter)
    {
    }

    graph_type* gp;
    size_type nind; 
    size_type list_node; 
    edge_list_iterator inciter;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      gp = nullptr;
      nind = size_type(-1);
    }

    EdgeIterator(const graph_type* c_graph, size_type c_edge):
                 gp(const_cast<graph_type*>(c_graph)), nind(c_edge) {};

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    Edge operator*() const {
      //return gp->Edge(nind); 
      return Edge(gp, gp->all_edges[nind].first, 
        gp->all_edges[nind].second);
    }

    // EdgeIterator& operator++()
    EdgeIterator& operator++(){
      if(nind < gp->num_all_edges){
        nind++; 
      }
      return *this;
    } 

    // bool operator==(const EdgeIterator&) const
    bool operator==(const EdgeIterator & eiter) const
    {
      return (nind == eiter.nind && gp == eiter.gp);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* gp;
    size_type nind;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  edge_iterator edge_begin() const {
      return EdgeIterator(this, 0);
  }


  // edge_iterator edge_end() const
  edge_iterator edge_end() const {
      return EdgeIterator(this, num_all_edges);
  }

   private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> all_nodes;
  using edge_st = std::pair<size_type, size_type>;
  edge_st get_node_indices(const Node & a, const Node & b) const
  {
    if(a.index() > b.index())
    {
      return {b.index(), a.index()};
    }
    return {a.index(), b.index()};
  }

  std::vector<edge_st> all_edges;
  std::map<edge_st, size_type> all_edges_list;
  size_type num_all_nodes;
  size_type num_all_edges;

  std::vector<std::vector<edge_st>> all_nodes_list;
  std::vector<node_value_type> nval;

};

#endif // CME212_GRAPH_HPP
