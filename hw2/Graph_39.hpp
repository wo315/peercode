#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <unordered_map>
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
template <typename V, typename E>
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
  /** Synonym for value type of Edge. */
  using edge_value_type = E;

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
  // Internal type for a graph node
  struct InternalNode {
    Point pos_;  // Node position
    node_value_type val_;  // Node value
    size_type idx_;  // Node index
    bool is_removed_;  // Whether node has been removed
  };
  // Internal type for edge uid (pair of node uids)
  using EdgeUid = std::pair<size_type, size_type>;
  // Hash function for internal edge uid type (for use with unordered_set)
  // Adapted from https://www.geeksforgeeks.org/how-to-create-an-unordered_map-of-pairs-in-c/
  struct hash_edgeuid {
    size_t operator()(const EdgeUid& p) const {
        size_t hash1 = std::hash<size_type>{}(p.first);
        size_t hash2 = std::hash<size_type>{}(p.second);
        return hash1 ^ hash2;
    }
  };
  // Internal type for a graph edge
  struct InternalEdge {
    edge_value_type val_;  // Edge value
    size_type idx_;  // Edge index
  };
  // Internal type for storage of edges incident to a single node
  using NodeEdges = std::unordered_set<EdgeUid, hash_edgeuid>;

 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_(), idx2uid_(), edges_(), edgeidx2uid_(), node_edges_() {
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //
  struct InvalidNodeError {};

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
      : graph_(nullptr), uid_(0) {
    }

    /** Return the node's position
     * @pre This node is valid and a member of a Graph
     * @return reference to this node's position
     *
     * Complexity: O(1)
     */
    Point& position() {
      return graph_->nodes_[uid_].pos_;
    }

    /** Return this node's position.
     * @pre This node is valid and a member of a Graph
     * @return reference to this node's position
     *
     * Complexity: O(1)
     */
    const Point& position() const {
      return graph_->nodes_[uid_].pos_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      if (graph_ == nullptr || graph_->nodes_[uid_].is_removed_)
        throw InvalidNodeError();

      return graph_->nodes_[uid_].idx_;
    }

    /** Return the node's value
     * @pre This node is valid and a member of a Graph
     * @return reference to this node's value
     *
     * Complexity: O(1)
     */
    node_value_type& value() {
      return graph_->nodes_[uid_].val_;
    }

    /** Return the node's value
     * @pre This node is valid and a member of a Graph
     * @return const reference to this node's value
     *
     * Complexity: O(1)
     */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].val_;
    }

    /** Return the number of edges incident to this node
     * @pre This node is valid and a member of a Graph
     *
     * Complexity: O(1)
     */
    size_type degree() const {
      return graph_->node_edges_[uid_].size();
    }

    /** Returns iterator pointing to the first edge from this node
     * @pre This node is valid and a member of a Graph
     *
     * Complexity: O(1)
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, graph_->node_edges_[uid_].begin());
    }

    /** Returns iterator pointing one past the last edge from this node
     * @pre This node is valid and a member of a Graph
     *
     * Complexity: O(1)
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, graph_->node_edges_[uid_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (uid_ == n.uid_);
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
        || ((graph_ == n.graph_) && (uid_ < n.uid_));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // This nodes's uid
    size_type uid_;
    // Private constructor
    Node(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), uid_(uid) {
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return idx2uid_.size();
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
    size_type next_idx = num_nodes();
    size_type next_uid = nodes_.size();

    nodes_.push_back({ position, value, next_idx, false });
    idx2uid_.push_back(next_uid);
    node_edges_.emplace_back();

    return Node(this, next_uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Test idx <-> uid match to determine if uid has been deleted
    return this == n.graph_ && !nodes_[n.uid_].is_removed_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, idx2uid_[i]);
  }

  /**
   * @brief Removes the node given by @a n if it exists
   * @param n Node to remove
   * @return size_type The number of nodes removed
   * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1
   *       Else,                  new num_nodes() == old num_nodes()
   * @post has_node(@a n) == false
   * @post has_edge(@a n, <any node>) == false
   * @post has_edge(<any node>, @a n) == false
   *
   * If the node doesn't exist, there is no effect.
   * Otherwise,
   *   Node @a n is invalidated.
   *   All nodes not equal to @a n remain valid.
   *   Can invalidate node indexes -- in other words, old node(@a i) might not
   *     equal new node(@a i).
   *   Existing NodeIterators are invalidated.
   *
   *   Edges involving @a n are invalidated.
   *   All edges not involving @a n remain valid edges.
   *   Can invalidate edge indexes -- in other words, old edge(@a i) might not
   *     equal new edge(@a i).
   *   Existing IncidentIterators created from @a n or a neighbor of @a n are invaliated.
   *   All other IncidentIterators remain valid.
   *   EdgeIterators are invalidated if @a n.degree() > 0.
   *
   * Complexity: Average O(@a n.degree())
   *             Worst-case O(num_edges()) if rebucketing occurs
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n))
      return 0;

    // Remove edges from node
    auto n_it = n.edge_begin();
    while (n_it != n.edge_end()) {
      remove_edge(*n_it);

      n_it = n.edge_begin();
    }

    // Swap last node (by index) into idx being vacated
    size_type idx = n.index();
    size_type last_node_uid = idx2uid_[num_nodes() - 1];
    nodes_[last_node_uid].idx_ = idx;
    idx2uid_[idx] = last_node_uid;
    idx2uid_.pop_back();

    // Set node has removed
    nodes_[n.uid_].is_removed_ = true;

    return 1;
  }

  /**
   * @brief Removes the node given by iterator @a n_it
   * @pre Iterator @a n_it points to a valid node of the graph
   * @param n_it Iterator pointing to node to be removed
   * @return node_iterator Iterator pointing to next node
   *
   * Delegated to remove_node(const Node& n) for @a n=(*n_it).
   * See that method for more documentation.
   */
  node_iterator remove_node(node_iterator n_it) {
    size_type node_uid = *(n_it.nodes_it_);
    size_type removed_idx = nodes_[node_uid].idx_;

    remove_node(*n_it);

    return NodeIterator(this, idx2uid_.begin() + removed_idx);
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
      : graph_(nullptr), uid1_(0), uid2_(0) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, uid2_);
    }

    /** Return the edge's value
     * @pre This edge is valid and a member of a Graph
     * @return reference to this edge's value
     *
     * Complexity: O(1)
     */
    edge_value_type& value() {
      EdgeUid edge_uid = graph_->make_edge_uid(node1(), node2());

      return graph_->edges_[edge_uid].val_;
    }

    /** Return the edge's value
     * @pre This edge is valid and a member of a Graph
     * @return const reference to this edge's value
     *
     * Complexity: O(1)
     */
    const edge_value_type& value() const {
      EdgeUid edge_uid = graph_->make_edge_uid(node1(), node2());

      return graph_->edges_[edge_uid].val_;
    }

    /** Return the resting length of the edge
     * @pre This edge is valid and a member of a Graph
     *
     * Complexity: O(1)
     */
    double length() const {
      return norm(node1().position() - node2().position());
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
    // Uids of nodes defining edge
    size_type uid1_, uid2_;
    // Private constructor
    Edge(const Graph* graph, size_type uid1, size_type uid2)
      : graph_(const_cast<Graph*>(graph)), uid1_(uid1), uid2_(uid2) {
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return edgeidx2uid_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    const EdgeUid edge_uid = edgeidx2uid_[i];

    return Edge(this, edge_uid.first, edge_uid.second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: Average O(1). Worst-case O(num_edges()).
   */
  bool has_edge(const Node& a, const Node& b) const {
    EdgeUid edge_uid = make_edge_uid(a, b);

    return (edges_.count(edge_uid) > 0);
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
   * Complexity: Average O(1). Worst-case O(num_edges())
   */
  Edge add_edge(const Node& a, const Node& b) {
      EdgeUid edge_uid = make_edge_uid(a, b);

    // Insert edge if it doesn't exist
    // Store reference to inserted edge for quick index lookup
    // Associate both nodes with the edge
    if (edges_.count(edge_uid) == 0) {
      size_type edge_idx = num_edges();

      InternalEdge edge = {edge_value_type(), edge_idx};

      edges_.insert({edge_uid, edge});
      edgeidx2uid_.push_back(edge_uid);
      node_edges_[a.uid_].insert(edge_uid);
      node_edges_[b.uid_].insert(edge_uid);
    }

    return Edge(this, a.uid_, b.uid_);
  }

  /**
   * @brief Removes the edge between @a a and @a b from the graph if it exists
   * @param a First node of the edge
   * @param b Second node of the edge
   * @return size_type The number of edges removed
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1
   *       Else,                        new num_edges() == old num_edges()
   * @post has_edge(@a a, @a b) == false
   *
   * If the edge doesn't exist, there is no effect.
   * Otherwise,
   *   All existing nodes are unaffected.
   *   Existing NodeIterators remain valid.
   *   All edges not between @a a and @a b remain valid edges.
   *   Edges between @a a and @a b are invalidated.
   *   Can invalidate edge indexes -- in other words, old edge(@a i) might not
   *     equal new edge(@a i).
   *   Existing IncidentIterators created from @a a and @a b are invaliated.
   *   All other IncidentIterators remain valid.
   *   All EdgeIterators are invalidated.
   * 
   * Complexity: Average O(1)
   *             Worst-case O(num_edges()) if rebucketing occurs
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (!has_edge(a, b))
      return 0;

    // Get edge's uid and idx
    EdgeUid edge_uid = make_edge_uid(a, b);
    size_type edge_idx = edges_[edge_uid].idx_;

    // Swap last edge (by index) into idx being vacated
    EdgeUid last_edge_uid = edgeidx2uid_[num_edges() - 1];
    edges_[last_edge_uid].idx_ = edge_idx;
    edgeidx2uid_[edge_idx] = last_edge_uid;
    edgeidx2uid_.pop_back();

    // Remove edge from adjacency list of both nodes
    node_edges_[edge_uid.first].erase(edge_uid);
    node_edges_[edge_uid.second].erase(edge_uid);

    // Remove edge from edges unordered map
    edges_.erase(edge_uid);

    return 1;
  }

  /**
   * @brief Removes the edge @a e from the graph if it exists
   * @param e Edge to be removed
   *
   * Delegated to remove_edge(const Node& a, const Node& b)
   *  for @a a=e.node1() and @a b=e.node2().
   * See that method for more documentation.
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /**
   * @brief Removes the edge given by iterator @a e_it
   * @pre Iterator @a e_it points to a valid edge of the graph
   * @param e_it Iterator pointing to edge to be removed
   * @return edge_iterator Iterator pointing to next edge
   *
   * Delegated to remove_edge(const Edge& e) for @a e=(*e_it).
   * See that method for more documentation.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    EdgeUid edge_uid = *(e_it.edges_it_);
    size_type removed_idx = edges_[edge_uid].idx_;

    remove_edge(*e_it);

    return EdgeIterator(this, edgeidx2uid_.begin() + removed_idx);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   *
   * Complexity: O(num_nodes() + num_edges())
   */
  void clear() {
    nodes_.clear();
    idx2uid_.clear();

    edges_.clear();
    edgeidx2uid_.clear();
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
      : graph_(nullptr), nodes_it_() {
    }

    /** Return the node at the current iterator position
     * @pre The NodeIterator is valid and points to a valid Node
     * @return a node of the graph
     *
     * Complexity: O(1)
     */
    Node operator*() const {
      return Node(graph_, *nodes_it_);
    }

    /** Increment the iterator position by one
     * @pre The NodeIterator is valid
     * @return a reference to the modified NodeIterator
     *
     * Complexity: O(1)
     */
    NodeIterator& operator++() {
      ++nodes_it_;
      return *this;
    }

    /** Test whether this iterator and @a nit are equal.
     *
     * Equal iterators have the same graph and the same index.
     */
    bool operator==(const NodeIterator& nit) const {
      return (graph_ == nit.graph_) && (nodes_it_ == nit.nodes_it_);
    }

   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // Internal type for iterator over nodes
    using nodes_iterator_ = std::vector<size_type>::const_iterator;
    // Iterator to underlying graph's node storage
    nodes_iterator_ nodes_it_;
    // Private constructor
    NodeIterator(const Graph* graph, nodes_iterator_ nodes_it)
      : graph_(const_cast<Graph*>(graph)), nodes_it_(nodes_it) {
    }
  };

  /** Returns iterator pointing to the first node of the graph
   *
   * Complexity: O(1)
   */
  node_iterator node_begin() const {
    return NodeIterator(this, idx2uid_.begin());
  }

  /** Returns iterator pointing one past the last node of the graph
   *
   * Complexity: O(1)
   */
  node_iterator node_end() const {
    return NodeIterator(this, idx2uid_.end());
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
    IncidentIterator() : graph_(nullptr), node_uid_(0), edges_it_() {
    }

    /** Return the edge at the current iterator position
     * @pre The IncidentIterator is valid and points to a valid Edge
     * @return an Edge object e with e.node1() equal to the node that created the IncidentIterator
     *
     * Complexity: O(1)
     */
    Edge operator*() const {
      EdgeUid edge_uid = *edges_it_;

      size_type second_node_uid = (edge_uid.first == node_uid_)
        ? edge_uid.second
        : edge_uid.first;

        return Edge(graph_, node_uid_, second_node_uid);
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
        && (node_uid_ == iit.node_uid_)
        && (edges_it_ == iit.edges_it_);
    }

   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // Uid of node to which edges are incident
    size_type node_uid_;
    // Internal type for iterator over NodeEdges
    using edges_iterator_ = typename NodeEdges::const_iterator;
    // Iterator to underlying graph's edge storage
    edges_iterator_ edges_it_;
    // Private constructor
    IncidentIterator(const Graph* graph, int node_uid, edges_iterator_ edges_it)
      : graph_(const_cast<Graph*>(graph)), node_uid_(node_uid),
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
    EdgeIterator() : graph_(nullptr), edges_it_() {
    }

    /** Return the edge at the current iterator position
     * @pre The EdgeIterator is valid and points to a valid Edge
     * @return an edge of the graph
     *
     * Complexity: O(1) due to implementation of Graph::edge()
     */
    Edge operator*() const {
      EdgeUid edge_uid = *edges_it_;

      return Edge(graph_, edge_uid.first, edge_uid.second);
    }

    /** Increment the iterator position by one
     * @return a reference to the modified EdgeIterator
     *
     * Complexity: O(1)
     */
    EdgeIterator& operator++() {
      ++edges_it_;
      return *this;
    }

    /** Test whether this iterator and @a eit are equal.
     *
     * Equal iterators refer to the same edge in the same graph.
     */
    bool operator==(const EdgeIterator& eit) const {
      return (graph_ == eit.graph_) && (edges_it_ == eit.edges_it_);
    }

   private:
    friend class Graph;
    // Pointer back to the Graph container
    Graph* graph_;
    // Internal type for iterator over edges
    using edges_iterator_ = std::vector<EdgeUid>::const_iterator;
    // Iterator to underlying graph's edge storage
    edges_iterator_ edges_it_;
    // Private constructor
    EdgeIterator(const Graph* graph, edges_iterator_ edges_it)
      : graph_(const_cast<Graph*>(graph)), edges_it_(edges_it) {
    }
  };

    /** Returns iterator pointing to the first edge of the graph
   *
   * Complexity: O(1)
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edgeidx2uid_.begin());
  }

  /** Returns iterator pointing one past the last edge of the graph
   *
   * Complexity: O(1)
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edgeidx2uid_.end());
  }

 private:
  /** Helper function to create internal edge uid with ordered node uids
   * @pre @a a and @a b are valid nodes
   * @return an EdgeUid object p with p.first < p.second
   * Complexity: O(1)
   */
  EdgeUid make_edge_uid(const Node& a, const Node& b) const {
    return a.uid_ < b.uid_ ? std::make_pair(a.uid_, b.uid_)
      : std::make_pair(b.uid_, a.uid_);
  }

  // Vector of graph's nodes, indexed by uid
  std::vector<InternalNode> nodes_;
  // Map the node index to uid, indexed by idx
  std::vector<size_type> idx2uid_;

  // Unordered map of graph's undirected edges, keyed by uid
  std::unordered_map<EdgeUid, InternalEdge, hash_edgeuid> edges_;
  // Map the edge index to graph edge, indexed by idx
  std::vector<EdgeUid> edgeidx2uid_;
  // Vector of NodeEdges to holds graph's edges per node,
  // indexed by node uid of first node in edge
  std::vector<NodeEdges> node_edges_;
};

#endif // CME212_GRAPH_HPP
