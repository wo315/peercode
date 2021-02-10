#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <algorithm>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 * @pre   The value type _V_ must have a default initializer.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int>
class Graph {
 private:

  /** Predeclaration of Internal Node type. */
  struct InternalNode_;
  /** Synonym for InternalNode_ (following STL conventions) **/
  using internal_node_type = InternalNode_;

  /** Predeclaration of Internal Edge type. **/
  struct InternalEdge_;
  /** Synonym for InternalEdge_ (following STL conventions) **/
  using internal_edge_type = InternalEdge_;

  /** Predeclaration of type of container for InteranlNode_, i.e.,
   * Graph::nodes_.
   */
  using node_container_type = std::vector<internal_node_type>;

  /** Predeclaration of type of container for InternalEdge_, i.e.,
   * Graph::edges_.
   *
   * NOTE: Edges are stored in a vector, rather than, e.g., a map with two node
   * indices as keys, to facilitate fast lookup by uid_ in Edge::fetch(), which
   * is required for this implementation of the proxy design pattern.
   */
  using edge_container_type = std::vector<internal_edge_type>;

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

  /** Synonym for V. */
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

 private:

  /** Predeclaration of type of container for incident edges to node, i.e.,
   * Graph::InternalNode_::incident_edges.
   */
  using incident_container_type = std::vector<size_type>;

 public:

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

    /** Return this node's position in 3-space. */
    const Point& position() const {
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size)
     * recording the order in which it was added to the graph. */
    size_type index() const {
      return uid_;
    }
    
    /** @brief Fetch a reference to the value of a node. */
    node_value_type& value() {
      return fetch().value;
    }

    /** @brief Fetch a constant reference to the value of a node. */
    const node_value_type& value() const {
      return fetch().value;
    }

    /** @brief Return the degree of a node. (I.e., the number of incident
     *         edges.)
     */
    size_type degree() const {
      return fetch().incident_edges.size();
    }

    /** @brief  yields an iterator pointing to the first edge incident to a
     *          node.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, fetch().incident_edges.begin());
    }

    /** @brief  yields an iterator pointing to the last edge incident to a
     *          node.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, fetch().incident_edges.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && uid_ == n.uid_);
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
      // Dictionary order in terms of graph location and then index.
      return graph_ != n.graph_ ? graph_ < n.graph_ : uid_ < n.uid_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    size_type uid_;             // This element's unique identifier.
    graph_type * const graph_;  // Pointer to parent graph.

    /** Private Constructor */
    Node(size_type uid, const graph_type* graph)
      : uid_(uid), graph_(const_cast<graph_type*>(graph)) {
    }

    /** Helper method to return the corresponding InternalNode_ */
    internal_node_type& fetch() const {
      return graph_->nodes_[uid_];
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
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
  node_type add_node(const Point& position, const node_value_type& value =
                node_value_type {}) {
    size_type i = num_nodes();
    // Data for newest node goes to back of vector.
    nodes_.push_back(InternalNode_(i, position, value)); 
    return Node(i, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // A node is only in the graph if and only if its parent graph is this graph
    // and its index is in the range [0, num_nodes())
    return (this == n.graph_ && n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(const size_type i) const {
    // We can only return a node if the index is valid, i.e., in the range [0,
    // num_nodes()).
    assert(i < num_nodes());
    return Node(i, this);
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
    Edge() {}

    /** Return a node of this Edge. By default, returns the node with the
     * smaller index, unless `reverse_` is set. */
    node_type node1() const {
      // Use parent graph `node()` method to get node by index.
      return reverse_
        ? graph_->node(fetch().node2)
        : graph_->node(fetch().node1);
    }

    /** Return the other node of this Edge. By default, returns the node with
     * the larger index, unless `reverse_` is set. */
    node_type node2() const {
      // Use parent graph `node()` method to get node by index. If reversed, get
      // smaller-indexed node instead of larger.
      return reverse_
        ? graph_->node(fetch().node1)
        : graph_->node(fetch().node2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Delegate to operator==(const InternalEdge_& e).
      return (fetch() == e.fetch());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Delegate to operator<(const InternalEdge_& e).
      return (fetch() < e.fetch());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    // NOTE: The "orientation" of an edge---i.e., which node is returned by
    // node1() and which by node2()---is not a property of the underlying data,
    // which is symmetric in the two nodes, but rather of this wrapper. This
    // orientation is stored in `reverse_`.

    size_type uid_;             // This element's unique identifier.
    graph_type * const graph_;  // Pointer to parent graph.
    bool reverse_;              // Reverse order nodes are output if set.

    /** Private Constructor */
    Edge(size_type uid, const graph_type* graph, const bool reverse = false)
      : uid_(uid), graph_(const_cast<graph_type*>(graph)), reverse_(reverse) {
    }

    /** Helper method to return the corresponding InternalEdge_ */
    internal_edge_type& fetch() const {
      return graph_->edges_[uid_];
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  edge_type edge(const size_type i) const {
    return Edge(i, this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Loop through all of the edges. If one of the edges joins the correct
    // nodes, return true. Otherwise, return false.
    
    // If one of the nodes is not in the graph, return false.
    if (a.graph_ != this || b.graph_ != this) return false;

    size_type min_index = std::min(a.index(), b.index());
    size_type max_index = std::max(a.index(), b.index());
    for(InternalEdge_ e : edges_) {
      if (e.node1 == min_index && e.node2 == max_index) {
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
  edge_type add_edge(const Node& a, const Node& b) {
    // `a` and `b` must be distinct nodes of this graph.
    assert(a.graph_ == this && b.graph_ == this && a != b);

    // Since Graph is undirected, InternalEdge_ has canonical representation in
    // terms of pair of node indices where smaller index always is first.
    size_type min_index = std::min(a.index(), b.index());
    size_type max_index = std::max(a.index(), b.index());

    // Next index is current number of edges.
    size_type i = num_edges();

    // If the index of `a` is greater than the index of `b`, reverse to ensure
    // that a is return.node1() and b is return.node2()
    bool reverse = a.index() > b.index();

    // Check that edge does not already exist; if it does, return early. Only
    // have to loop through edges in one of the nodes, since both nodes have to
    // be incident for the edge to join them.
    for (size_type j : a.fetch().incident_edges) {
      if (edges_[j].node1 == min_index && edges_[j].node2 == max_index) {
        return Edge(j, this, reverse);
      }
    }

    // If edge does not exist, push a new edge to the end of the edge container,
    // and add that it is incident to both nodes.
    edges_.push_back(InternalEdge_(i, a.index(), b.index()));
    a.fetch().incident_edges.push_back(i);
    b.fetch().incident_edges.push_back(i);

    return Edge(i, this, reverse);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
  }

  /** @brief  Sets a default value for all nodes.
   *  @param  v Default value for nodes.
   *  @post   n.value() == v for every node n.
   */
  void default_value(const V& v) {
    for (internal_node_type& n : nodes_) n.value = v;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: public equality_comparable<node_iterator> {
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

    /** Dereference a NodeIterator to get the underlying Node. */
    node_type operator*() const {
      return graph_->node(it_->uid);
    }

    /** Advance the NodeIterator to the next node. */
    node_iterator& operator++() {
      ++it_;
      return *this;
    }

    /** Compare two iterators to see if they're equal. */
    bool operator==(const node_iterator& it) const {
      return graph_ == it.graph_ && it_ == it.it_;
    }

   private:
    friend class Graph;

    // Pointer to parent graph.
    const graph_type * graph_;
    // Iterator pointing to current node.
    typename node_container_type::const_iterator it_;

    /** Private constructor from iterator over Graph::node_ */
    NodeIterator(const graph_type * graph, typename
                 node_container_type::const_iterator it) :
      graph_(graph), it_(it) {
    }
  };

  /** @brief  Yields an iterator pointing to the first node in the graph. */
  node_iterator node_begin() const {
    return NodeIterator(this, nodes_.begin());
  }

  /** @brief  Yields an iterator pointing to the last node in the graph. */
  node_iterator node_end() const {
    return NodeIterator(this, nodes_.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: public equality_comparable<incident_iterator> {
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

    /** Dereference an IncidentIterator to get the underlying Node. */
    edge_type operator*() {
      edge_type edge = graph_->edge(*it_);
      // Returned edge must have root node as first node.
      edge.reverse_ = graph_->node(node_uid_) > edge.node1();
      return edge;
    }

    /** Advance the NodeIterator to the next node. */
    incident_iterator& operator++() {
      ++it_;
      return *this;
    }

    /** Compare two iterators to see if they're equal. */
    bool operator==(const incident_iterator& it) const {
      return (graph_ == it.graph_
              && node_uid_ == it.node_uid_
              && it_ == it.it_);
    }

   private:
    friend class Graph;

    // Pointer to parent graph.
    graph_type * const graph_;
    // Unique identifier of relevant node.
    size_type node_uid_;
    // Iterator pointing to current edge.
    typename incident_container_type::iterator it_;

    /** Private constructor from iterator over
     * Graph::InternalNode_::incident_edges */
    IncidentIterator(graph_type* graph, size_type node_uid, typename
                     incident_container_type::iterator it) :
      graph_(graph), node_uid_(node_uid), it_(it) {
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: public equality_comparable<edge_iterator> {
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

    /** Dereference a EdgeIterator to get the underlying edge. */
    edge_type operator*() const {
      return graph_->edge(it_->uid);
    }

    /** Advance the EdgeIterator to the next edge. */
    edge_iterator& operator++() {
      ++it_;
      return *this;
    }

    /** Compare two iterators to see if they're equal. */
    bool operator==(const edge_iterator& it) const {
      return graph_ == it.graph_ && it_ == it.it_;
    }

   private:
    friend class Graph;

    // Pointer to parent graph.
    const graph_type * graph_;
    // Iterator pointing to current node.
    typename edge_container_type::const_iterator it_;

    /** Private constructor from iterator over Graph::node_ */
    EdgeIterator(const graph_type * graph, typename
        edge_container_type::const_iterator it) :
      graph_(graph), it_(it) {
    }
  };

  /** @brief  Yields an iterator pointing to the first edge in the graph. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edges_.begin());
  }

  /** @brief  Yields an iterator pointing to the last edge in the graph. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.end());
  }

 private:

  //
  // Internal Node
  //

  /** @class Graph::InternalNode_
   * @brief Class containing node internals. Proxied by Graph::Node.
   *
   * NOTE: Since InternalNode_ does not store a reference to containing graph,
   * there is no notion of equality: two nodes with identical internals should
   * not be equal if they belong to different graphs. Therefore, equality is
   * only implemented at the level of the Node wrapper.
   */
  struct InternalNode_ {
    const size_type uid;                          // The unique identifier for a node.
    const Point position;                         // The position of the node.
    node_value_type value;                        // The value contained by the node.
    std::vector<size_type> incident_edges;        // Indices of incident edges.

    /** Constructor for InternalNode_ */
    InternalNode_(const size_type& uid, const Point& position, node_value_type
                  value)
      : uid(uid), position(position), value(value), incident_edges() {
    }
  };

  //
  // Internal Edge
  //

  /** @class Graph::InternalEdge_
   * @brief Class containng edge internals. Proxied by Graph::Edge */
  struct InternalEdge_ : private totally_ordered<internal_edge_type> {

    // NOTE: InternalEdges_ represent edges uniquely by always putting the node
    // with the smaller index in `node1`. The "orientation" of an edge is
    // properly a property of the Edge proxy, rather than the underlying data,
    // represented by the Edge::reverse_ member.

    // Unique identifier of edge.
    const size_type uid;
    // The uid of the first node in the edge. This is the node with the smaller index.
    const size_type node1;
    // The uid of the second node in the edge. This is the node with the larger index.
    const size_type node2;

    /** Constructor for InternalNode_ */
    InternalEdge_(const size_type uid, const size_type node1, const size_type node2) : 
      uid(uid), node1(std::min(node1, node2)), node2(std::max(node1, node2)) {
    }

    /** @brief  Test whether this edge is less than @a e in a graph-specific order.
     *  @pre    This edge and _e_ must be in the same graph.
     */
    bool operator<(const InternalEdge_& e) const {
      // Dictionary order in terms of the smaller node index and then the larger
      // node index. Since edges have a unique representation as @a
      // InternalEdge_, this makes the definition independent of the order the
      // edges are given in.
      return node1 != e.node1 ? node1 < e.node1 : node2 < e.node2;
    }

    /** @brief  Test whether this edge is the same as @a e in a global order. To
     *          be equal, edges must join the same nodes in the same graph.
     *  @pre  This node and _e_ must be in the same graph.
     */
    bool operator==(const InternalEdge_& e) const {
      return node1 == e.node1 && node2 == e.node2;
    }
  };

  node_container_type nodes_;  // Collection of nodes in graph.
  edge_container_type edges_;  // Collection of edges in graph.
};

#endif // CME212_GRAPH_HPP
