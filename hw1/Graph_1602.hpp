#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
#include <tuple>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int>
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
  /** Synonym for template value type */
  using node_value_type = V;
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->point_vec_.at(index_);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    /** Return this node's value. */
    node_value_type& value(){
      return graph_->node_values_.at(index_);
    }
    const node_value_type& value() const{
      return graph_->node_values_.at(index_);
    }
    
    /** Return the degree of this node. */
    size_type degree() const{
      //Check if valid index
      return graph_->incident_vec_.at(index_).size();
    }
    
    /** IncidentIterator begin(). */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, index_, 0);
    }

    /** IncidentIterator end(). */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, index_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (this->graph_ == n.graph_) & (this->index_ == n.index_);
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     * Ordering: graph
     */
    bool operator<(const Node& n) const {
      // Can assume this->graph_ == n.graph_ (https://piazza.com/class/kjx6t0a2rv4589?cid=75)
      return this->index_ < n.index_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    uint index_;
    // Private Constructor
    Node(const Graph* graph, int index) : graph_(const_cast<Graph*>(graph)), index_(index) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return point_vec_.size();
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

  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()){
    point_vec_.push_back(position);
    // Add new incident_vec_ entry for the new node
    incident_vec_.push_back(std::vector<uint>{});
    node_values_.push_back(node_value);
    return Node(this, (point_vec_.size()-1));
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //Check that it's a member of this graph and graph hasn't been cleared
    return (this == n.graph_) & (point_vec_.size() > n.index_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
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
  class Edge : private totally_ordered<Edge>  {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node_a_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node_b_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     * Checks wether the edges are from the same graph, then whether they are the same edge.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ != e.graph_){
        return false;
      }
      return ((this->node_a_ == e.node_a_ and this ->node_b_ == e.node_b_) or (this->node_a_ == e.node_b_ and this ->node_b_ == e.node_a_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //Assume this->graph_ == e.graph_ (https://piazza.com/class/kjx6t0a2rv4589?cid=75)
      //Find smaller node in each edge
      uint a_min_ = std::min(this->node_a_, this->node_b_);
      uint b_min_ = std::min(e.node_a_, e.node_b_);
      //If they have the same first node, check the second node (by sum)
      if (a_min_ == b_min_){
        return (this->node_a_ + this-> node_b_ < e.node_a_ + e.node_b_);
        }
      //Else, return the edge with the smaller first node
      return (a_min_ < b_min_);
      }
      

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    uint node_a_;
    uint node_b_;
    // Private constructor
    Edge(const Graph* graph, uint node_a_, uint node_b_) : graph_(const_cast<Graph*>(graph)), node_a_(node_a_), node_b_(node_b_){};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_vec_.size();
    }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    uint node_a_ = std::get<0>(edges_vec_.at(i));
    uint node_b_ = std::get<1>(edges_vec_.at(i));
    return Edge(this, node_a_, node_b_);  
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Actual complexity: O(a.degree())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Verify valid edge for first call
    return (std::find(incident_vec_.at(a.index_).begin(), incident_vec_[a.index_].end(), b.index_) != incident_vec_[a.index_].end());
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
    if (has_edge(a, b)){
      return Edge(this, a.index_, b.index_);
    }
    // Here, we don't use at because these indices are guaranteed by has_edge.
    incident_vec_[a.index_].push_back(b.index_);
    incident_vec_[b.index_].push_back(a.index_);
    edges_vec_.push_back(std::make_tuple(a.index_, b.index_));
    return Edge(this, a.index_,b.index_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0s
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    point_vec_.clear();
    incident_vec_.clear();
    edges_vec_.clear();
    node_values_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>  {
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

    /** Dereference NodeIterator. */
    Node operator*() const{
      return Node(graph_, it_);
    }

    /** Increment NodeIterator. */
    NodeIterator& operator++(){
      ++it_;
      return *this;
    }

    /** Check NodeIterator equivalence.
     * Checks that it points to the same graph and the same iterator.
     */
    bool operator==(const NodeIterator& node_it_b) const{
      return ((graph_ == node_it_b.graph_) & (it_ == node_it_b.it_));
    }

   private:
    friend class Graph;
    Graph* graph_;
    /**  We only call nodes by index (not by a reference to an object), so iterate over indices.
     *   The index refers to both the point (point_vec_) and the value (node_values_).
     *   It can also be used to call incident edges to a given node.
     */
    uint it_;
    // Private constructor
    NodeIterator (const Graph* graph, uint it = 0) : graph_(const_cast<Graph*>(graph)), it_(it) {}

  };

  /** Return a NodeIterator pointing to the first node. */
  node_iterator node_begin() const{
    return(NodeIterator(this, 0));
  }

  /** Return a NodeIterator pointing past the last node. */
  node_iterator node_end() const{
    return(NodeIterator(this, point_vec_.size()));
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
    IncidentIterator() {
    }

    /** Derefernce IncidentIterator. */
    Edge operator*() const{
      return Edge(graph_, node_1_index_, graph_->incident_vec_[node_1_index_][node_2_index_]);
    }

    /** Increment IncidentIterator. */
    IncidentIterator& operator++(){
      node_2_index_++;
      return *this;
    }

    /** Check for IncidentIterator Equality.
     *  Checks that it points to the same graph and the same indices.
     */
    bool operator==(const IncidentIterator& incident_2) const{
      return ((graph_ == incident_2.graph_) & (node_1_index_ == incident_2.node_1_index_)&(node_2_index_ == incident_2.node_2_index_));
    }
   private:
    friend class Graph;
    Graph* graph_;
    // The first node is constant and is also the index of node_1.
    const uint node_1_index_;
    /** Here, we iterate over the indices of <b>a vector of nodes incident to node_1 </b>
     * We must do it this way because we only store incident nodes and therefore the actual index
     * is not ordered.
     * Ex: incident_vec_[node_1_index_][node_2_index_] = index of (node_2_index_)th node incident to node_1_.
     */
    uint node_2_index_;
    // Private constructor
    IncidentIterator(const Graph* graph, const uint node_1_index = 0, uint node_2_index = 0) : graph_(const_cast<Graph*>(graph)), node_1_index_(node_1_index),node_2_index_(node_2_index) {}
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
    EdgeIterator() {
    }
    /** Dereference EdgeIterator. */
    Edge operator*() const{
      std::tuple<uint, uint> tmp_edge = graph_->edges_vec_[edge_index_];
      return Edge(graph_, std::get<0>(tmp_edge), std::get<1>(tmp_edge));
    }

    /** Increment EdgeIterator. */
    EdgeIterator& operator++(){
      edge_index_++;
      return *this;
    }

    /** Check for EdgeIterator Equality.
     *  Checks that it points to the same graph and the same edges.
     */
    bool operator==(const EdgeIterator& edge_iterator_2) const{
      return ((graph_ == edge_iterator_2.graph_) & (edge_index_ == edge_iterator_2.edge_index_));
    }

   private:
    friend class Graph;
    Graph* graph_;
    /**  We can call edges by index (rather than by its nodes), so iterate over indices.
     *   If we iterated over incident edges, we would have to check to not double count.
     *   By calling by index, we guarantee to count each edge once and only once.
     */
    uint edge_index_;
    // Private onstructor
    EdgeIterator(const Graph* graph, uint edge_index = 0) : graph_(const_cast<Graph*>(graph)), edge_index_(edge_index) {}
  };

  /** Return EdgeIterator pointing to first edge */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Return EdgeIterator pointing past last edge */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_vec_.size());
  }

 private:
  // Vector containing node's point objects
  std::vector<Point> point_vec_;
  // Vector containing node's values
  std::vector<node_value_type> node_values_;
  // Vector of vectors for incident nodes.
  // incident_vec_[a] contains the index of all nodes incident to a.
  std::vector<std::vector<uint>> incident_vec_;
  // Vector of edges, as tuple. Used for calling edges by index.
  std::vector<std::tuple<uint, uint>> edges_vec_;
};

#endif // CME212_GRAPH_HPP
