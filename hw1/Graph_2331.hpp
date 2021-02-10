#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <utility>
#include <vector>

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
  /** Synonym for value type of Node. */
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
  // Internal type for a graph node (pair of location and value)
  using InternalNode = std::pair<Point, node_value_type>;
  // Internal type for a graph edge (pair of node indexes)
  using InternalEdge = std::pair<size_type, size_type>;
  // Internal type for storage of edges incident to a single node
  using NodeEdges = std::vector<InternalEdge const*>;

 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_(), edges_(), edge_idx_refs_(), node_edges_() {
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
  class Node : totally_ordered<Node> {
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
      : graph_(nullptr), idx_(0) {
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[idx_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx_;
    }

    /** Return the node's value
     * @pre This node is valid and a member of a Graph
     * @return reference to this node's value
     *
     * Complexity: O(1)
     */
    node_value_type& value() {
      return graph_->nodes_[idx_].second;
    }

    /** Return the node's value
     * @pre This node is valid and a member of a Graph
     * @return const reference to this node's value
     *
     * Complexity: O(1)
     */
    const node_value_type& value() const {
      return graph_->nodes_[idx_].second;
    }

    /** Return the number of edges incident to this node
     * @pre This node is valid and a member of a Graph
     *
     * Complexity: O(1)
     */
    size_type degree() const {
      return graph_->node_edges_[idx_].size();
    }

    /** Returns iterator pointing to the first edge from this node
     * @pre This node is valid and a member of a Graph
     *
     * Complexity: O(1)
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, idx_, graph_->node_edges_[idx_].begin());
    }

    /** Returns iterator pointing one past the last edge from this node
     * @pre This node is valid and a member of a Graph
     *
     * Complexity: O(1)
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, idx_, graph_->node_edges_[idx_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (index() == n.index());
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
      return (graph_ < n.graph_)
        || ((graph_ == n.graph_) && (index() < n.index()));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // This nodes's index
    size_type idx_;
    // Private constructor
    Node(const Graph* graph, size_type idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx) {
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
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type()) {
    nodes_.push_back(make_internal_node(position, value));
    node_edges_.emplace_back();

    return node(num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_;
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
  class Edge : totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge()
      : graph_(nullptr), idx1_(0), idx2_(0) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node(idx1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node(idx2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((node1() == e.node1()) && (node2() == e.node2()))
        || ((node1() == e.node2()) && (node2() == e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (node1() < e.node1())
        || ((node1() == e.node1()) && (node2() < e.node2()));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Pointer to the Graph container
    Graph* graph_;
    // Index of nodes defining edge
    size_type idx1_, idx2_;
    // Private constructor
    Edge(const Graph* graph, size_type idx1, size_type idx2)
      : graph_(const_cast<Graph*>(graph)), idx1_(idx1), idx2_(idx2) {
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
  Edge edge(size_type i) const {
    InternalEdge const* edge_ref = edge_idx_refs_[i];

    return Edge(this, edge_ref->first, edge_ref->second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    InternalEdge edge = make_internal_edge(a, b);

    return (edges_.count(edge) > 0);
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
      InternalEdge edge = make_internal_edge(a, b);

    // Insert edge if it doesn't exist
    // Store reference to inserted edge for quick index lookup
    // Associate both nodes with the edge
    if (edges_.count(edge) == 0) {
      auto insert_result = edges_.insert(edge);
      InternalEdge const& inserted_edge = *(insert_result.first);

      edge_idx_refs_.push_back(&inserted_edge);
      node_edges_[a.index()].push_back(&inserted_edge);
      node_edges_[b.index()].push_back(&inserted_edge);
    }

    return Edge(this, a.index(), b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    edge_idx_refs_.clear();
    node_edges_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
      : graph_(nullptr), idx_(0) {
    }

    /** Return the node at the current iterator position
     * @pre The NodeIterator is valid and points to a valid Node
     * @return a node of the graph
     *
     * Complexity: O(1)
     */
    Node operator*() const {
      return graph_->node(idx_);
    }

    /** Increment the iterator position by one
     * @pre The NodeIterator is valid
     * @return a reference to the modified NodeIterator
     *
     * Complexity: O(1)
     */
    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }

    /** Test whether this iterator and @a nit are equal.
     *
     * Equal iterators have the same graph and the same index.
     */
    bool operator==(const NodeIterator& nit) const {
      return (graph_ == nit.graph_) && (idx_ == nit.idx_);
    }

   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // Iterator to underlying graph's node storage
    int idx_;
    // Private constructor
    NodeIterator(const Graph* graph, int idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx) {
    }
  };

  /** Returns iterator pointing to the first node of the graph
   *
   * Complexity: O(1)
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Returns iterator pointing one past the last node of the graph
   *
   * Complexity: O(1)
   */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : node_idx_(0), edges_it_() {
    }

    /** Return the edge at the current iterator position
     * @pre The IncidentIterator is valid and points to a valid Edge
     * @return an Edge object e with e.node1() equal to the node that created the IncidentIterator
     *
     * Complexity: O(1)
     */
    Edge operator*() const {
      InternalEdge const* edge_ref = *edges_it_;

      size_type second_node = (edge_ref->first == node_idx_)
        ? edge_ref->second
        : edge_ref->first;

        return Edge(graph_, node_idx_, second_node);
    }

    /** Increment the iterator position by one
     * @pre The IncidentIterator is valid
     * @return a reference to the modified IncidentIterator
     *
     * Complexity: O(1)
     */
    IncidentIterator& operator++() {
      ++edges_it_;
      return *this;
    }

    /** Test whether this iterator and @a iit are equal.
     *
     * Equal iterators refer to the same edge in the same set of incident edges.
     */
    bool operator==(const IncidentIterator& iit) const {
      return (graph_ == iit.graph_)
        && (node_idx_ == iit.node_idx_)
        && (edges_it_ == iit.edges_it_);
    }

   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // Index of node to which edges are incident
    size_type node_idx_;
    // Internal type for iterator over NodeEdges
    using edges_iterator_ = NodeEdges::const_iterator;
    // Iterator to underlying graph's edge storage
    edges_iterator_ edges_it_;
    // Private constructor
    IncidentIterator(const Graph* graph, int node_idx, edges_iterator_ edges_it)
      : graph_(const_cast<Graph*>(graph)), node_idx_(node_idx),
        edges_it_(edges_it) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(nullptr), idx_() {
    }

    /** Return the edge at the current iterator position
     * @pre The EdgeIterator is valid and points to a valid Edge
     * @return an edge of the graph
     *
     * Complexity: O(1) due to implementation of Graph::edge()
     */
    Edge operator*() const {
      return graph_->edge(idx_);
    }

    /** Increment the iterator position by one
     * @return a reference to the modified EdgeIterator
     *
     * Complexity: O(1)
     */
    EdgeIterator& operator++() {
      ++idx_;
      return *this;
    }

    /** Test whether this iterator and @a eit are equal.
     *
     * Equal iterators refer to the same edge in the same graph.
     */
    bool operator==(const EdgeIterator& eit) const {
      return (graph_ == eit.graph_) && (idx_ == eit.idx_);
    }

   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // Iterator to underlying graph's edge storage
    int idx_;
    // Private constructor
    EdgeIterator(const Graph* graph, int edge_idx)
      : graph_(const_cast<Graph*>(graph)), idx_(edge_idx) {
    }
  };

    /** Returns iterator pointing to the first edge of the graph
   *
   * Complexity: O(1)
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Returns iterator pointing one past the last edge of the graph
   *
   * Complexity: O(1)
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:
  /** Helper function to create internal node type
   * @param[in] p The node's position
   * @param[in] val The node's value
   * @return an InternalNode object
   * Complexity: O(1)
   */
  InternalNode make_internal_node(const Point& p, const node_value_type& val) {
    return std::make_pair(p, val);
  }

  /** Helper function to create internal edge type with ordered node indexes
   * @pre @a a and @a b are valid nodes
   * @return an InternalEdge object p with p.first < p.second
   * Complexity: O(1)
   */
  InternalEdge make_internal_edge(const Node& a, const Node& b) const {
    size_type a_idx = a.index();
    size_type b_idx = b.index();

    return a_idx < b_idx ? std::make_pair(a_idx, b_idx)
      : std::make_pair(b_idx, a_idx);
  }
  // Hash function for internal edge type (for use with unordered_set)
  // Adapted from https://www.geeksforgeeks.org/how-to-create-an-unordered_map-of-pairs-in-c/
  struct hash_internaledge {
    size_t operator()(const InternalEdge& p) const {
        size_t hash1 = std::hash<size_type>{}(p.first);
        size_t hash2 = std::hash<size_type>{}(p.second);
        return hash1 ^ hash2;
    }
  };

  // Vector of graph's nodes
  std::vector<InternalNode> nodes_;
  // Unordered set of graph's undirected edges (stored as pairs of node indexes)
  std::unordered_set<InternalEdge, hash_internaledge> edges_;
  // Vector of references to graph's edges, ordered by insertion order
  std::vector<InternalEdge const*> edge_idx_refs_;
  // Vector of NodeEdges to holds graph's edges per node
  std::vector<NodeEdges> node_edges_;
};

#endif // CME212_GRAPH_HPP
