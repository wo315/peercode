#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

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

  /** Type of node value  */
  using node_value_type = V;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  private:


    // node_vec_ is a vector of pointers to Point variables,
    // representing the nodes of the graph
    std::vector< const Point *> node_vec_;

    // node_val_ is a vector containing the values of the nodes
    std::vector<node_value_type> node_val_;

    // adj_ is a vector such that edge_vec_[i] contains the indexes of the
    // neighbors of the node of index i. It is basically an adjacency list
    std::vector<std::vector<size_type>> adj_;

    // edge_vec_ is a vector of pairs of indexes such that, at index i,
    // there is the pair of indexes of the nodes that form edge of index i
    std::vector<std::pair<size_type, size_type>> edge_vec_;

  public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(): node_vec_(std::vector<const Point*>()),
           node_val_(std::vector<node_value_type>()),
           adj_(std::vector<std::vector<size_type>>()),
           edge_vec_(std::vector<std::pair<size_type, size_type>>())
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
  class Node: private totally_ordered<Node> {
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
      graph_ = nullptr;
      index_ = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // Look into node_vec_ with the appropriate index
      return *((*graph_).node_vec_[index_]);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    /** Return this node's value */
    node_value_type& value(){
      return (*graph_).node_val_[index_];
    }

    /** Return this node's value */
    const node_value_type& value() const{
      return (*graph_).node_val_[index_];
    }

    /** Return this node's degree (number of incident edges)
     * @post degree() = size(adj_[index_])
     *
     * Complexity: O(1)
    */
    size_type degree() const{
      return (*graph_).adj_[index_].size();
    }

    /** Return an incident iterator at the first incident edge at this node */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, index_, 0);
    }

    /** Return an incident iterator after the last incident edge at this node */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, index_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // As stated in @75 in Piazza, we can assume the nodes to
      return (graph_ == n.graph_ and index_ == n.index_);
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
      // As stated in @75 in Piazza, we can assume the nodes to be from
      // the same graph
      return index_ < n.index_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Define a node with a pointer to its graph and its index
    graph_type* graph_;
    size_type index_;

    /** Valid node constructor */
    Node(const graph_type* graph, size_type index) :
         graph_(const_cast<graph_type*>(graph)), index_(index) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_vec_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_value The new node's value (optional)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() == node_value
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node (const Point & position,
                 const node_value_type & node_value = node_value_type()){
    // Add the Point pointer to node_vec_
    Point* p = new Point;
    *p = position;
    node_vec_.push_back(p);

    // Add the value to the vector of values
    node_val_.push_back(node_value);

    // Add an empty vector to the adjacency list
    adj_.push_back(std::vector<size_type>());
    return Node(this, size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_ and n.index_ < size()) ;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i);
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      graph_ = nullptr;
      index1_ = 0;
      index2_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, index1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, index2_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Equality <=> same graph and same nodes IN ANY ORDER
      return graph_ == e.graph_ and
             ((index1_ == e.index1_ and index2_ == e.index2_)
             or (index1_ == e.index2_ and index2_ == e.index1_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // As stated in @75 in Piazza, we assume the edges come from same graph
      // Lexicographic order between edges
      // represented as (smallest index, largest index)
      size_type this_min = std::min(index1_, index2_);
      size_type this_max = std::max(index1_, index2_);
      size_type e_min = std::min(e.index1_, e.index2_);
      size_type e_max = std::max(e.index1_, e.index2_);

      return (this_min < e_min or (this_min == e_min and this_max < e_max));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // An Edge is represented with a pointer to its graph and the two indexes
    // of the nodes it connects

    graph_type* graph_;
    size_type index1_;
    size_type index2_;


    // Edge constructor
    Edge(const graph_type* graph, size_type index1, size_type index2) :
         graph_(const_cast<graph_type*>(graph)),
         index1_(index1), index2_(index2){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_vec_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // It can easily be found in the vector edge_vec_
    return Edge(this, edge_vec_[i].first, edge_vec_[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // It is sufficent to look into the neighbors of a.
    // https://stackoverflow.com/questions/571394/

    return std::find(adj_[a.index_].begin(),adj_[a.index_].end(), b.index_)
           != adj_[a.index_].end() ;
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
    else {
      // We add a to the neighbors of b and b to the neighbors of a
      adj_[a.index_].push_back(b.index_);
      adj_[b.index_].push_back(a.index_);

      // Add the index pair to edge_vec_
      edge_vec_.push_back(std::make_pair(a.index_, b.index_));

      // e.node1() == a  and e.node2() == b is verified
      return Edge(this, a.index_, b.index_);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // We use the clear method on vectors

    adj_.clear();
    edge_vec_.clear();
    node_val_.clear();

    // Avoid memory leak
    for (auto p: node_vec_){
        delete p;
    }
    node_vec_.clear();
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
      graph_ = nullptr;
      index_ = 0;
    }

    /** Dereference operator */
    Node operator*() const {
      return Node(graph_, index_);
    }

    /** Increments to the next Node in the graph */
    NodeIterator& operator++() {
      index_++;
      return *this;
    }

    /** Defines equality between iterators
     *  @param[in] ni Node iterator to which we compare this node iterator
    */
    bool operator==(const NodeIterator& ni) const{
      return (graph_ == ni.graph_ and index_ == ni.index_);
    }

    /** Defines inequality between iterators
     *  @param[in] ni Node iterator to which we compare this node iterator
    */
    bool operator!=(const NodeIterator& ni) const{
      return (!(*this == ni));
    }


   private:
    friend class Graph;

    // Pointer to the graph on which we iterate
    graph_type* graph_;

    // Index of current node
    size_type index_;

    /** Valid NodeIterator constructor */
    NodeIterator(const graph_type* graph, size_type index) :
         graph_(const_cast<graph_type*>(graph)), index_(index) {
    }

  };

  /** Return a node iterator at the first node of this graph */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** Return a node iterator after the last node of this graph */
  node_iterator node_end() const{
    return NodeIterator(this, this->node_vec_.size());
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
      graph_ = nullptr;
      main_index_ = 0;
      neighbor_index_  = 0;
    }

    /** Dereference operator */
    Edge operator*() const {
      return Edge(graph_, main_index_,
                  (*graph_).adj_[main_index_][neighbor_index_]);
    }

    /** Increments to the next edge incident to this node */
    IncidentIterator& operator++() {
      neighbor_index_++;
      return *this;
    }

    /** Defines equality between iterators
     *  @param[in] ii Incident iterator to which we compare this node iterator
    */
    bool operator==(const IncidentIterator& ii) const {
      return (graph_ == ii.graph_
              and main_index_ == ii.main_index_
              and neighbor_index_ == ii.neighbor_index_);
    }

    /** Defines inequality between iterators
     *  @param[in] ii Incident iterator to which we compare this node iterator
    */
    bool operator!=(const IncidentIterator& ii) const {
      return (!(*this == ii));
    }


   private:
    friend class Graph;

    // Pointer to the graph
    graph_type* graph_;

    // Index of the main node (the one of which we consider the neighbors)
    size_type main_index_;

    // Index in the adjacency list (!= node.index()) of the current neighbor
    // of node(main_index_)
    size_type neighbor_index_;

    /** Valid IncidentIterator constructor */
    IncidentIterator(const graph_type* graph, size_type i1, size_type i2) :
         graph_(const_cast<graph_type*>(graph)),
         main_index_(i1), neighbor_index_(i2){
    }

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
      graph_ = nullptr;
      index_ = 0;
    }


    /** Dereference operator */
    Edge operator*() const {
      return Edge(graph_,
                  (*graph_).edge_vec_[index_].first,
                  (*graph_).edge_vec_[index_].second);
    }

    /** Increments to the next edge in the graph */
    EdgeIterator& operator++() {
      index_++;
      return *this;
    }

    /** Defines equality between iterators
     *  @param[in] ei Edge iterator to which we compare this node iterator
    */
    bool operator==(const EdgeIterator& ei) const{
      return (graph_ == ei.graph_ and index_ == ei.index_);
    }

    /** Defines inequality between iterators
     *  @param[in] ei Edge iterator to which we compare this node iterator
    */
    bool operator!=(const EdgeIterator& ei) const{
      return (!(*this == ei));
    }

   private:
    friend class Graph;
    // Note that we could have combined NodeIterator and IncidentIterator to
    // construct EdgeIterator. Instead, we take advantage of edge_vec_ which
    // allows a simple implementation of EdgeIterator.

    // Pointer to the graph
    graph_type* graph_;

    // Index of current edge in edge_vec_
    size_type index_;

    /** Valid EdgeIterator constructor */
    EdgeIterator(const graph_type* graph, size_type index) :
         graph_(const_cast<graph_type*>(graph)), index_(index) {
    }

  };

  /** Return an edge iterator at the first edge of this graph */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Return an edge iterator after the last edge of this graph */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

};

#endif // CME212_GRAPH_HPP
