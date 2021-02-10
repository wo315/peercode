#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <list> // linked list from STL 

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

  /** Predeclaration of Internal Node type. */
  struct internal_node_;
  /** Synonym for InternalNode_ (following STL conventions) **/
  using internal_node_type = internal_node_;

  /** Predeclaration of Internal edge type */
  struct internal_edge_;
  /** Synonym for internal_edge_ (following STL conventions) **/
  using internal_edge_type = internal_edge_;

public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Synonym for value of a node */ 
  using node_value_type = V;

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

  /** Alias for the type of iterator used to traverse the std::vector container for nodes */ 
  using stl_node_iterator = typename std::vector<internal_node_type>::const_iterator; 

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Alias for the type of iterator used to traverse the std::vector container for edges */
  using stl_edge_iterator = typename std::vector<internal_edge_type>::const_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Alias for the type of iterator used to traverse the std::vector container for incident edges */
  using stl_incident_iterator = typename std::vector<unsigned>::const_iterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    n_nodes_ = 0;
    n_edges_ = 0; 
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
   * This class is lightweight and thus sizeof(Graph::Node) <= 16 bytes. 
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
      my_graph_ = nullptr;
      idx_ = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().pos_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx_;
    }

    /** Return this node's degree. */ 
    size_type degree() const{
      return fetch().incident_edges_.size(); 
    }

    /** Return this node's value. */
    node_value_type& value() {
      return fetch().value_; 
    };

    /** Return this node's value (overloaded const version of the above method). */
    const node_value_type& value() const {
      return fetch().value_;
    };

    /** return an iterator pointing at the start of the incident edge list. */
    incident_iterator edge_begin() const
    {
      return IncidentIterator(fetch().incident_edges_.begin(), fetch().my_graph_, idx_); 
    }

    /** return an iterator pointing at the end (aka one past the last edge) of the incident edge list. */
    incident_iterator edge_end() const
    {
      return IncidentIterator(fetch().incident_edges_.end(), fetch().my_graph_, idx_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((this->my_graph_ == n.my_graph_) & (this->idx_ == n.idx_));
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
      // if from same graph, order by node index. Otherwise, order by 
      // address of each graph. Inspired by HW0 peercode #276.
      return (my_graph_ != n.my_graph_) ? (my_graph_ < n.my_graph_) : (this->index() < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the parent graph. 8 bytes of memory. 
    Graph* my_graph_;

    // node's ID in the global list of nodes for the graph. 4 bytes of memory. 
    size_type idx_;

    /** Private Constructor */
    Node(const Graph *graph, size_type idx)
        : my_graph_(const_cast<Graph*>(graph)), idx_(idx)
    {
    }

    // Helper method to return the corresponding internal node. 
    // Inspired by peercode for HW0
    internal_node_type& fetch() const
    {
      return my_graph_->nodes_[idx_];
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return n_nodes_; 
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
  Node add_node(const Point& position, const node_value_type& value=node_value_type()) {

    // add node to list of nodes
    nodes_.push_back(internal_node_(this, n_nodes_, position, value));

    // increment size of graph
    n_nodes_ += 1;

    // Returns a proxy to the new node
    return Node(this, n_nodes_ - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return ((this == n.my_graph_) & (n.idx_ < n_nodes_)); 
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < n_nodes_); 
    return Node(this, i); // return a proxy object for element i, as in proxy_example.cpp 
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   * This class is lightweight and thus sizeof(Graph::Edge) <= 32 bytes. 
   */
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      my_graph_ = nullptr;
      idx_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // return a proxy node. use reverse_order_ flag to make sure we return the correct node 
      return reverse_order_ ? Node(my_graph_, fetch().node_b) : Node(my_graph_, fetch().node_a); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // return a proxy node. use reverse_order_ flag to make sure we return the correct node
      return reverse_order_ ? Node(my_graph_, fetch().node_a) : Node(my_graph_, fetch().node_b);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((this->my_graph_ == e.my_graph_) & (this->idx_ == e.idx_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // if from same graph, order by edge index. Otherwise, order by
      // address of each graph. Inspired by HW0 peercode #276.
      return (my_graph_ != e.my_graph_) ? (my_graph_ < e.my_graph_) : (this->idx_ < e.idx_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer back to the parent graph
    Graph *my_graph_;

    // edge's ID in the global list of edge for the graph
    size_type idx_;

    // whether to return nodes in ascending or descending order of their index
    // we need this so that we can satisfy e.node1() == @a a and e.node2() == @a b
    // regardless of the relative indices of the input nodes in add_edge() 
    bool reverse_order_; 

    /** Private Constructor */
    Edge(const Graph *graph, size_type idx, bool reverse_order = false)
        : my_graph_(const_cast<Graph*>(graph)), idx_(idx), reverse_order_(reverse_order)
    {
    }

    // Helper method to return the corresponding internal edge. 
    // Inspired by peercode for HW0 
    internal_edge_type fetch() const
    {
      return my_graph_->edges_[idx_];
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return n_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i); // return a proxy object for element i, as in proxy_example.cpp
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    // figure out which of a and b has the smaller and larger index
    size_type small_idx = std::min(a.index(), b.index());
    size_type large_idx = std::max(a.index(), b.index());

    // loop through the edges of node a checking to see if the candidate edge
    // connects node a and node b. 
    size_type vec_a_size =  a.fetch().incident_edges_.size(); 
    for (size_type i = 0; i < vec_a_size; i++)
    {
      size_type current_edge_idx = a.fetch().incident_edges_[i];
      if (edges_[current_edge_idx].node_a == small_idx && edges_[current_edge_idx].node_b == large_idx)
      {
        return true;
      }
    }

    // If we reach this point, then an appropriate edge was not found 
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
  // Note: revised after HW0 using ideas from peercode to meet runtime complexity 
  Edge add_edge(const Node& a, const Node& b) {

    // check that a and b are present in the graph and that they are not the same node
    assert((this->has_node(a) && this->has_node(b)) && !(a == b));

    // take note of what the new edges index would be
    size_type new_edge_idx = n_edges_; 

    // set flag that keeps track of whether or not node a has the smaller index
    bool reverse_flag = (a.index() < b.index()) ? false : true;

    // loop through the edges of node a checking to see if the candidate edge
    // connects node a and node b. If so, return proxy edge. 
    if (has_edge(a, b)){
      return Edge(this, new_edge_idx, reverse_flag);
    }

    // If we reach this point, then the edge doesn't exist. 
    // In this case, we should push_back the new edge
    // and update the incident_edges of node a and node b
    // we also need to update the number of edges
    edges_.push_back(internal_edge_(this, new_edge_idx, a.index(), b.index()));
    a.fetch().incident_edges_.push_back(new_edge_idx);
    b.fetch().incident_edges_.push_back(new_edge_idx);
    n_edges_ +=  1; 

    // Return a proxy to the new edge
    return Edge(this, new_edge_idx, reverse_flag);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

    // Zero out node and edge count 
    n_nodes_ = 0;
    n_edges_ = 0;

    // Empty out containers
    nodes_.empty(); 
    edges_.empty(); 

    // Call clear methods
    nodes_.clear(); 
    edges_.clear(); 
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Defines equality between two node iterators as having the same 
     *    underlying STL iterator. */ 
    bool operator==(const node_iterator& node_iter) const
    {
      return node_iter.it_ == it_; 
    }

    /** Defines inequality between two node iterators as having different
     *    underlying STL iterators. */
    bool operator!=(const node_iterator& node_iter) const 
    {
      return node_iter.it_ != it_;
    }

    /* Defines the dereference operator for the node iterator. 
    *   Returns a proxy node having the same graph and index. */ 
    node_type operator*() const 
    {
      return Node(it_->my_graph_, it_->idx_); // return proxy node
    }

    /** Defines the increment operator for the node iterator.
     *    Increments the underlying STL iterator. */
    node_iterator& operator++()
    {
      ++it_; 
      return *this; 
    }

   private:
    friend class Graph;

    /** STL iterator used for the node container */
    stl_node_iterator it_; 

    /** Private constructor that can be accessed by the Graph class. */ 
    NodeIterator(stl_node_iterator it) : it_(it) {}
  };

  // Returns a node iterator pointing at the start of the node list.
  node_iterator node_begin() const
  {
    return node_iterator(nodes_.begin());
  }

  // Returns a node iterator pointing at the end (aka one past the last node) of the node list.
  node_iterator node_end() const
  {
    return node_iterator(nodes_.end()); 
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Defines equality between two incident iterators as having the same 
     *    underlying STL iterator. */
    bool operator==(const incident_iterator& incident_iter) const
    {
      return incident_iter.it_ == it_;
    }

    /** Defines inequality between two incident iterators as having different
     *    underlying STL iterators. */
    bool operator!=(const incident_iterator& incident_iter) const
    {
      return incident_iter.it_ != it_;
    }

    /* Defines the dereference operator for the incident iterator. 
    *   Returns a proxy edge having the same graph and edge index. 
    *   Care is taken to ensure that the edge returned is such that
    *   Edge::node1() will always return the node used to create 
    *   this iterator. */
    edge_type operator*() const
    {
      // To note:
      // - By default node1() returns the node with the smaller index of the edge
      // - This can be flipped by passing the flag "reversed" 
      // Thus, to ensure Edge::node1() returns the spawning node index in this case of 
      // IncidentIterator(), we need to turn on the "reversed" flag if it's the larger index
      size_type current_edge_idx = *it_;
      size_type small_idx = graph_->edges_[current_edge_idx].node_a; // guranteed to be the smaller node index
      size_type large_idx = graph_->edges_[current_edge_idx].node_b; // guranteed to be the larger node index 

      // if the spawning node's index is the smaller index, then nothing needs to be done.
      // otherwise we need to return an edge with the reverse flag turned on. 
      bool reversed = (spawning_node_idx_ == small_idx) ? false : true;

      // return proxy edge
      return Edge(graph_, *it_, reversed); 
    }

    /** Defines the increment operator for the incident iterator.
     *    Increments the underlying STL iterator. */
    incident_iterator operator++()
    {
      ++it_;
      return *this; 
    }

   private:
    friend class Graph;

    // STL iterator for the container for incident edges
    stl_incident_iterator it_;

    // Pointer to the parent graph. Necessary for dereferencing since the incident edge container only stores indices. 
    const graph_type* graph_;

    // We need to keep track of the index of the spawning node so we can make sure to always return it via Edge::node1()
    size_type spawning_node_idx_;

    // Private constructor that can be accessed by the Graph class.
    IncidentIterator(stl_incident_iterator it, const graph_type *graph, size_type node_idx) 
      : it_(it), graph_(graph), spawning_node_idx_(node_idx) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Defines equality between two edge iterators as having the same 
     *    underlying STL iterator. */
    bool operator==(const edge_iterator& edge_iter) const
    {
      return edge_iter.it_ == it_;
    }

    /** Defines inequality between two edge iterators as having different
     *    underlying STL iterators. */
    bool operator!=(const edge_iterator& edge_iter) const
    {
      return edge_iter.it_ != it_;
    }

    /* Defines the dereference operator for the edge iterator. 
    *   Returns a proxy edge having the same graph and edge index. */
    edge_type operator*() const
    {
      return Edge(it_->my_graph_, it_->idx_); // return proxy edge
    }

    /** Defines the increment operator for the edge iterator.
     *    Increments the underlying STL iterator. */
    edge_iterator& operator++()
    {
      ++it_;
      return *this; 
    }

   private:
    friend class Graph;

    // STL iterator for the edge container
    stl_edge_iterator it_;

    // Private constructor that can be accessed by the Graph class.
    EdgeIterator(stl_edge_iterator it) : it_(it) {}
  };

  // Returns an iterator pointing at the start of the edge list.
  edge_iterator edge_begin() const
  {
    return EdgeIterator(edges_.begin());
  }

  // Returns an iterator pointing at the end (aka one past the last edge) of the edge list.
  edge_iterator edge_end() const
  {
    return EdgeIterator(edges_.end());
  }

 private:

  size_type n_nodes_;
  size_type n_edges_;

  // Internal type for nodes 
  // - Inspired by peercode from HW 0 
  struct internal_node_
  {
    // Pointer back to the parent graph
    Graph *my_graph_;

    // edge's ID in the global list of edge for the graph
    size_type idx_;

    // position of the node
    Point pos_;

    // value of this node
    node_value_type value_;

    // keep track of incident edges for this node
    std::vector<size_type> incident_edges_; 

    // constructor for internal_node_
    internal_node_(const Graph *graph, const size_type idx, const Point pos, const node_value_type value = node_value_type())
        : my_graph_(const_cast<Graph *>(graph)), idx_(idx), pos_(pos), value_(value), incident_edges_() {}
  };

  // Internal type for edges
  struct internal_edge_
  {
    // Pointer back to the parent graph
    Graph *my_graph_;

    // edge's ID in the global list of edge for the graph
    size_type idx_;

    size_type node_a; // id of the first node of the edge. Node with the smaller index. 
    size_type node_b; // id of the second node of the edge. Node with the larger index.

    // constructor for internal_edge_. Inspired by peercode for HW0. 
    internal_edge_(const Graph *graph, const size_type idx, const size_type node_a, const size_type node_b)
        : my_graph_(const_cast<Graph *>(graph)), idx_(idx) 
        {
          // this ensures internally node_a always refers to the node with the smaller index and node_b always
          // refers to the larger index 
          this->node_a = std::min(node_a, node_b);
          this->node_b = std::max(node_a, node_b);
        }
  };

  // Containers for nodes (their positions) and edges
  std::vector<internal_node_type> nodes_; 
  std::vector<internal_edge_type> edges_;     

};

#endif // CME212_GRAPH_HPP
