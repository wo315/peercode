#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
 public:
  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  using node_type = Node;  // Synonym for Node (following STL conventions).
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  using edge_type = Edge;  // Synonym for Edge (following STL conventions).
  using edge_value_type = E;

  class NodeIterator;                  // Iterator over all graph nodes.
  using node_iterator = NodeIterator;  // Synonym for NodeIterator

  class EdgeIterator;                  // Iterator over all graph edges.
  using edge_iterator = EdgeIterator;  // Synonym for EdgeIterator

  class IncidentIterator;  // Iterator over node incident edges.
  using incident_iterator = IncidentIterator;  // Synonym for IncidentIterator

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  /** Construct an empty graph. [HW0] */
  Graph() {}

  /** Default destructor */
  ~Graph() = default;

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node. [HW0]
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
    Node() {}

    /** Return this node's position. [HW2] */
    Point& position() { return g_->nodes_[uid_].position; }

    /** Return this node's position (read-only). [HW0] */
    const Point& position() const { return g_->nodes_[uid_].position; }

    /** Return this node's index, a number in the range [0, graph_size).
     * [HW0, HW2]
     */
    size_type index() const { return g_->nodes_[uid_].index; }

    /** Return a reference to this node's value. [HW1] */
    node_value_type& value() { return g_->nodes_[uid_].value; }

    /** Return a reference to this node's value (read-only). [HW1] */
    const node_value_type& value() const { return g_->nodes_[uid_].value; }

    /** Return the number of incident edges. [HW1] */
    size_type degree() const { return g_->nodes_[uid_].inc_edges.size(); }

    /** Return an iterator pointing to the first incident edge. [HW1] */
    incident_iterator edge_begin() const {
      return IncidentIterator(g_, uid_, 0);
    }

    /** Return an iterator pointing past the last incident edge. [HW1] */
    incident_iterator edge_end() const {
      return IncidentIterator(g_, uid_, g_->nodes_[uid_].inc_edges.size());
    }

    /** Test whether this node and @a n are equal. [HW0]
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return g_ == n.g_ && index() == n.index();
    }

    /** Test whether this node is less than @a n in a global order. [HW0]
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     *
     */
    bool operator<(const Node& n) const {
      return (uid_ != n.uid_) ? (uid_ < n.uid_) : (g_ < n.g_); 
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* g_;
    size_type uid_;  // Unique internal node ID.

    /** Private constructor accessed by Graph. [HW0] */
    Node(const Graph* g, size_type uid)
        : g_(const_cast<Graph*>(g)), uid_(uid) {}

    bool valid() const {
      return uid_ < g_->nodes_.size()        // uid in range.
             && index() < g_->ni2u_.size()   // idx in range.
             && g_->ni2u_[index()] == uid_;  // uid in sync.
    }
  };

  /** Return the number of nodes in the graph. [HW0]
   *
   * Complexity: O(1).
   */
  size_type size() const { return ni2u_.size(); }

  /** Synonym for size(). */
  size_type num_nodes() const { return size(); }

  /** Return the node with index @a i. [HW0]
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const { return Node(this, ni2u_[i]); }

  /** Determine if a Node belongs to this Graph. [HW0]
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return Node(this, ni2u_[n.index()]) == n;
  }

  /** Add a node to the graph, returning the added node. [Hw0]
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1).
   */
  Node add_node(const Point& position) {
    size_type next_uid = nodes_.size();
    size_type next_index = size();
    ni2u_.push_back(next_uid);
    nodes_.emplace_back(next_index, position);
    assert(Node(this, next_uid).valid());
    return Node(this, next_uid);
  }

  /** Add a node to the graph, returning the added node. [HW1]
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1).
   */
  Node add_node(const Point& position, const node_value_type& value) {
    size_type next_uid = nodes_.size();
    size_type next_index = num_nodes();
    ni2u_.push_back(next_uid);
    nodes_.emplace_back(next_index, position, value);
    assert(Node(this, next_uid).valid());
    return Node(this, next_uid);
  }

  /** Remove a node and its incident edges from the graph.
   * @param [in] n The node to be removed.
   * @return 1 if an node was removed, 0 otherwise.
   * @post If old has_node(@a a, @a b), new num_nodes() == old num_nodes() - 1.
   *       Else,                        new num_nodes() == old num_nodes().
   * @post graph.node(i).index() == i for i in [0, new num_nodes())
   *
   * Can invalidate node indexes, node iterators, edge iterators, and incident
   * edge iterators.
   *
   * Complexity: O(n.degree()).
   */
  size_type remove_node(const Node& n) {
    if (!n.valid()) return 0;  // Node doesn't exist.

    // Remove all incident edges.
    for (auto ei = n.edge_begin(); ei != n.edge_end();) {
      remove_edge(*ei);
    }

    nodes_[ni2u_.back()].index = n.index();  // Update last node's index.
    remove(ni2u_, n.index());  // Remove node and replace with last node.

    assert(!n.valid());
    return 1;
  }

  /** Remove a node and its incident edges from the graph.
   * @param [in] n_it An iterator pointing to the edge to be removed.
   * @returns A valid node iterator or the end iterator if the removed node was
   * the last node.
   */
  node_iterator remove_node(node_iterator n_it) {
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
    /** Construct an invalid Edge. [HW0] */
    Edge() {}

    /** Return a node of this Edge [HW0] */
    Node node1() const { return g_->node_uid(g_->edges_[uid_].nid1); }

    /** Return the other node of this Edge [HW0] */
    Node node2() const { return g_->node_uid(g_->edges_[uid_].nid2); }

    /** Return this node's index, a number in the range [0, num_edges). */
    size_type index() const { return g_->edges_[uid_].index; }

    /** Return a reference to this edge's value. [HW2] */
    edge_value_type& value() { return g_->edges_[uid_].value; }

    /** Return a reference to this edge's value (read-only). [HW2] */
    const edge_value_type& value() const { return g_->edges_[uid_].value; }

    /** Compute the length of the edge. [HW2] */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Test whether this edge and @a e are equal. [HW0]
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1() == e.node1() && node2() == e.node2()) ||
             (node1() == e.node2() && node2() == e.node1());
    }

    /** Test whether this edge is less than @a e in a global order. [HW0]
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     *
     */
    bool operator<(const Edge& e) const { 
      return (uid_ != e.uid_) ? (uid_ < e.uid_) : (g_ < e.g_); 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* g_;
    size_type uid_;

    /** Private constructor accessed by Graph. [HW0] */
    Edge(const Graph* g, size_type uid)
        : g_(const_cast<Graph*>(g)), uid_(uid) {}

    bool valid() const {
      return uid_ < g_->edges_.size()       // uid in range.
             && index() < g_->ei2u_.size()  // idx in range.
             && g_->ei2u_[index()] == uid_  // uid in sync.
             && node1().valid() && node2().valid();
    }
  };

  /** Return the total number of edges in the graph. [HW0]
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { return ei2u_.size(); }

  /** Return the edge with index @a i. [HW0]
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const { return Edge(this, ei2u_[i]); }

  /** Test whether two nodes are connected by an edge. [HW0]
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const { return index(a, b) > -1; }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * [HW0]
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: At most O(num_nodes() + num_edges()).
   */
  Edge add_edge(const Node& a, const Node& b) {
    int edge_index = index(a, b);
    if (edge_index > -1)
      return edge(edge_index, a.uid_, b.uid_);  // Edge already exists

    size_type next_uid = edges_.size();
    size_type next_index = num_edges();
    ei2u_.push_back(next_uid);
    edges_.emplace_back(next_index, a.uid_, b.uid_);

    edge_map_[pair(a, b)] = next_index;

    nodes_[a.uid_].inc_edges.push_back(next_uid);
    nodes_[b.uid_].inc_edges.push_back(next_uid);

    return edge(next_index);
  }

  /** Remove an edge from the graph between two nodes.
   * @param [in] a The first node of the edge to be removed.
   * @param [in] b The second node of the edge to be removed.
   * @return 1 if an edge was removed, 0 otherwise.
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate edge indexes, edge iterators, and incident edge iterators.
   *
   * Complexity: O(num_edges()).
   */
  size_type remove_edge(const Node& a, const Node& b) {
    int edge_index = index(a, b);
    if (edge_index == -1) return 0;  // Edge doesn't exist.

    // Remove from node incidence lists.
    size_type inc_index = 0;
    for (auto ei = a.edge_begin(); ei != a.edge_end(); ++ei) {
      if ((*ei).node2() == b) {
        remove(nodes_[a.uid_].inc_edges, inc_index);
        break;
      }
      inc_index++;
    }
    inc_index = 0;
    for (auto ei = b.edge_begin(); ei != b.edge_end(); ++ei) {
      if ((*ei).node2() == a) {
        remove(nodes_[b.uid_].inc_edges, inc_index);
        break;
      }
      inc_index++;
    }

    // Update last edge's index.
    edges_[ei2u_.back()].index = edge_index;

    // Remove edge and replace with last edge.
    remove(ei2u_, edge_index);
    Edge e = edge(edge_index);
    edge_map_[pair(e.node1(), e.node2())] = edge_index;
    edge_map_.erase(pair(a, b));

    return 1;
  }

  /** Remove an edge from the graph.
   * @param [in] e An edge to be removed.
   * @return 1 if an edge was removed, 0 otherwise.
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from the graph.
   * @param [in] e_it An iterator pointing to the edge to be removed.
   * @returns A valid edge iterator or the end iterator if the removed edge was
   * the last edge.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return e_it;
  }

  /** Remove all nodes and edges from this graph. [HW0]
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    ni2u_.clear();
    ei2u_.clear();
    edge_map_.clear();
  }

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                            // Element type
    using pointer = Node*;                              // Pointers to elements
    using reference = Node&;                            // Reference to elements
    using difference_type = std::ptrdiff_t;             // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** Dereference operator. [HW1] */
    Node operator*() const { return g_->node(index_); }

    /** Increments (prefix) to the next node in the graph. [Hw1] */
    NodeIterator& operator++() {
      ++index_;
      return *this;
    }

    /** Defines equality between two iterators. [HW1] */
    bool operator==(const NodeIterator& iter) const {
      return (g_ == iter.g_) && (index_ == iter.index_);
    }

   private:
    // Allow Graph to access NodeIterator's private member data and functions.
    friend class Graph;

    Graph* g_;         // Pointer to graph.
    size_type index_;  // Node index.

    /** Private constructor accessed by Graph. [HW1] */
    NodeIterator(const Graph* g, size_type index)
        : g_(const_cast<Graph*>(g)), index_(index) {}
  };

  /** Return an iterator pointing to the first node in the graph. [Hw1] */
  node_iterator node_begin() const { return NodeIterator(this, 0); }

  /** Return an iterator pointing past the last node in the graph. [Hw1] */
  node_iterator node_end() const { return NodeIterator(this, num_nodes()); }

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                            // Element type
    using pointer = Edge*;                              // Pointers to elements
    using reference = Edge&;                            // Reference to elements
    using difference_type = std::ptrdiff_t;             // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /** Dereference operator. [HW1]
     *
     * @return An edge incident to node1().
     */
    Edge operator*() const {
      size_type uid = g_->ni2u_[nodei_];
      size_type edge_uid = g_->nodes_[uid].inc_edges[index_];
      size_type nid1 = g_->edges_[edge_uid].nid1;
      size_type nid2 = g_->edges_[edge_uid].nid2;

      if (uid == nid1) return g_->edge(edge_uid, nid1, nid2);

      return g_->edge(edge_uid, nid2, nid1);
    }

    /** Increments (prefix) to the next incident edge. [HW1] */
    IncidentIterator& operator++() {
      ++index_;
      return *this;
    }

    /** Defines equality between two iterators. [HW1] */
    bool operator==(const IncidentIterator& iter) const {
      return (g_ == iter.g_) && (nodei_ == iter.nodei_) &&
             (index_ == iter.index_);
    }

   private:
    // Allow Graph to access IncidentIterator's private member data and
    // functions.
    friend class Graph;

    Graph* g_;         // Pointer to graph.
    size_type nodei_;  // Node index [0, num_nodes).
    size_type index_;  // Incident edge index.

    /** Private constructor accessed by Graph. [HW1] */
    IncidentIterator(const Graph* g, size_type nodei, size_type index)
        : g_(const_cast<Graph*>(g)), nodei_(nodei), index_(index) {}
  };

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                            // Element type
    using pointer = Edge*;                              // Pointers to elements
    using reference = Edge&;                            // Reference to elements
    using difference_type = std::ptrdiff_t;             // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /** Dereference operator. [HW1] */
    Edge operator*() const { return g_->edge(index_); }

    /** Increments (prefix) to the next edge in the Graph. [HW1] */
    EdgeIterator& operator++() {
      ++index_;
      return *this;
    }

    /** Defines equality between two iterators. [HW1] */
    bool operator==(const EdgeIterator& iter) const {
      return (g_ == iter.g_) && (index_ == iter.index_);
    }

   private:
    // Allow Graph to access EdgeIterator's private member data and functions.
    friend class Graph;
    Graph* g_;         // Pointer to graph.
    size_type index_;  // Edge index.

    /** Private constructor accessed by Graph. [HW1] */
    EdgeIterator(const Graph* g, size_type index)
        : g_(const_cast<Graph*>(g)), index_(index) {}
  };

  /** Return an iterator pointing to the first edge in the graph. [HW1] */
  edge_iterator edge_begin() const { return EdgeIterator(this, 0); }

  /** Return an iterator pointing past the last edge in the graph. [HW1] */
  edge_iterator edge_end() const { return EdgeIterator(this, num_edges()); }

 private:
  /* Graph internals */
  struct NodeData {
    size_type index;  // User-facing index.
    Point position;
    node_value_type value;
    std::vector<size_type> inc_edges;  // Incident edge uids.
    NodeData(size_type i, Point p) : index(i), position(p), value() {}
    NodeData(size_type i, Point p, node_value_type v)
        : index(i), position(p), value(v) {}
  };

  struct EdgeData {
    size_type index;  // User-facing index.
    size_type nid1;   // First node uid.
    size_type nid2;   // Second node uid.
    edge_value_type value;
    EdgeData(size_type i, size_type n1, size_type n2)
        : index(i), nid1(n1), nid2(n2), value() {}
  };

  /** Stores data for any node added. Indexed by node uid. */
  std::vector<NodeData> nodes_;

  /** Stores data for any edge added. Indexed by edge uid. */
  std::vector<EdgeData> edges_;

  /** Maps node indices to uid. */
  std::vector<size_type> ni2u_;

  /** Maps edge indices to uid. */
  std::vector<size_type> ei2u_;

  /** Maps pairs of node uids to edge indices. */
  std::map<std::pair<size_type, size_type>, size_type> edge_map_;

  /** Return the node with the given uid. */
  Node node_uid(size_type uid) const { return Node(this, uid); }

  /** Return the index of the edge connecting two nodes. [HW0]
   * @pre @a a and @a b are valid nodes of this graph
   * @return An edge index in the range [0, num_edges) if it exists, else -1.
   */
  int index(const Node& a, const Node& b) const {
    auto search = edge_map_.find(pair(a, b));
    if (search != edge_map_.end()) return search->second;
    return -1;
  }

  /** Return the edge with index @a i and nodes @a a and @a b. [HW0]
   * @pre 0 <= @a i < num_edges()
   * @post result_edge.node1().uid_ == a
   * @post result_edge.node2().uid_ == b
   *
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   */
  Edge edge(size_type index, size_type nid1, size_type nid2) {
    size_type uid = ei2u_[index];
    edges_[uid].nid1 = nid1;
    edges_[uid].nid2 = nid2;
    return Edge(this, uid);
  }

  /** Removes an element from a vector, without preserving order.
   * @param [in] vec Reference to the vector to be modifed.
   * @param [in] i Index of the element to be removed.
   *
   * Complexity: O(1).
   */
  template <typename T>
  void remove(std::vector<T>& vec, size_type i) {
    vec[i] = vec.back();
    vec.pop_back();
  }

  /** Creates a unique pair of uids for two nodes.
   * @return A std::pair (a.uid_, b.uid_).
   * @post std::get<0>(result) < std::get<1>(result)
   */
  std::pair<size_type, size_type> pair(const Node& a, const Node& b) const {
    return (a < b) ? std::make_pair(a.uid_, b.uid_)
                   : std::make_pair(b.uid_, a.uid_);
  }
};

#endif  // CME212_GRAPH_HPP
