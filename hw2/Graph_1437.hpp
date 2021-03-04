#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = double>
class Graph {
 private:

  // Use this space for declarations of important internal types we need
  // later in the Graph's definition.

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

  /** Synonym for Node value, */
  using node_value_type = V;

  /** Predeclaration of Node information. */
  struct NodeInfo;
  /** Synonym for Node information (following STL conventions). */
  using node_info = NodeInfo;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Synonym for Edge value, */
  using edge_value_type = E;

  /** Predeclaration of Edge information. */
  struct EdgeInfo;
  /** Synonym for Edge information (following STL conventions). */
  using edge_info = EdgeInfo;

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
  Graph() : num_nodes_(0), node_i2u(), nodes_(), nodes_adj(),
            edges_(), num_edges_(0) {}

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** Custom structure of information to store with Nodes */
  struct NodeInfo {
    Point            position;        //< Node's position
    size_type        index;           //< Node's index
    node_value_type  value;           //< Node's value

    NodeInfo() : position(Point()), index(0), value(node_value_type()) {}

    /** Construct a NodeInfo by the given _node_pos_ , _node_idx_ , and
     *  _node_value.
     * @param[in] node_pos    Node position.
     * @param[in] node_idx    Node index.
     * @param[in] node_value  Node value.
     */
    NodeInfo(const Point& node_pos, size_type node_idx,
             const node_value_type& node_value) :
      position(node_pos), index(node_idx), value(node_value) {}
  };

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
      graph_ = nullptr;
      n_uid_ = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      if (!valid()) {
        throw std::runtime_error("The node is not valid.");
      }
      return graph_->nodes_.at(n_uid_).position;
    }

    /**
     * @brief Return the position of the node.
     * @return A reference to the position of _node_.
     *
     * @pre The node is a valid node.
     * @post The returned position can be modified.
     *
     * Complexity: O(1).
     */
    Point& position() {
      if (!valid()) {
        throw std::runtime_error("The node is not valid.");
      }
      return graph_->nodes_.at(n_uid_).position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodes_.at(n_uid_).index;
    }

    /**
     * @brief Return the value of the node.
     * @return A reference to the value of _node_.
     *
     * @pre The node is a valid node.
     * @post The returned value can be modified.
     *
     * Complexity: O(1).
     */
    node_value_type& value() {
      if (!valid()) {
        throw std::runtime_error("The node is not valid.");
      }
      return graph_->nodes_.at(n_uid_).value;
    }

    /**
     * @brief Return the value of the node.
     * @return A reference to the value of _node_.
     *
     * @pre The node is a valid node.
     * @post The returned value cannot be modified.
     *
     * Complexity: O(1).
     */
    const node_value_type& value() const {
      if (!valid()) {
        throw std::runtime_error("The node is not valid.");
      }
      return graph_->nodes_.at(n_uid_).value;
    }

    /**
     * @brief Return the degree of the node.
     * @return The degree of _node_.
     *
     * @pre The node is a valid node.
     */
    size_type degree() const {
      if (!valid()) {
        throw std::runtime_error("The node is not valid.");
      }
      return graph_->nodes_adj.at(n_uid_).size();
    }

    /**
     * @brief Return beginning of incident iterator to the current node.
     * @return The beginning of incident iterator to _node_.
     *
     * @pre The node is a valid node.
     */
    incident_iterator edge_begin() const {
      if (!valid()) {
        throw std::runtime_error("The node for incident iterator is"
                                 "not valid.");
      }
      return IncidentIterator(graph_, n_uid_,
                              graph_->nodes_adj.at(n_uid_).begin());
    }

    /**
     * @brief Return end of incident iterator to the current node.
     * @return The end of incident iterator to _node_.
     *
     * @pre The node is a valid node.
     */
    incident_iterator edge_end() const {
      if (!valid()) {
        throw std::runtime_error("The node for incident iterator is"
                                 "not valid.");
      }
      return IncidentIterator(graph_, n_uid_,
                              graph_->nodes_adj.at(n_uid_).end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same uid.
     */
    bool operator==(const Node& n) const {
      return graph_ == n.graph_ && n_uid_ == n.n_uid_;
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
      return !(*this == n) && (n_uid_ < n.n_uid_ || graph_ < n.graph_);
    }

    /** Set operator "=" so that we can assign a valid node
     * to an invalid node later.
     */
    Node& operator=(const Node& n) = default;

    /** Default destructor */
    ~Node() = default;

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to Graph
    graph_type* graph_;
    // The unique index for the node
    size_type n_uid_;

    /** Construct a valid node by the given _graph_ and _uid_.
     * @param[in] graph A pointer to the given Graph
     * @param[in] uid   A unique entry index for the node
     */
    Node(const graph_type* graph, size_type uid) {
      graph_ = const_cast<graph_type*>(graph);
      n_uid_ = uid;
    }

    /** Check whether the current node is valid.
     * There are three constraints to indicate whether a node is valid.
     * 1. uid in range: 0 <= node.uid < num_nodes_,
     * 2. idx in range: 0 <= node.idx < node_i2u.size(),
     * 3. uid in sync:  node_i2u.at(node.idx) = node.uid,
     *
     * Complexity: O(1).
     */
     bool valid() const {
      return n_uid_ >= 0 && n_uid_ < graph_->num_nodes_ &&
             graph_->nodes_.at(n_uid_).index < graph_->node_i2u.size() &&
             graph_->node_i2u.at(graph_->nodes_.at(n_uid_).index) == n_uid_;
    }

  };

  /** Return the number of valid nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_i2u.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_val The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& node_val = node_value_type()) {
    node_type new_node = Node(this, num_nodes_);
    nodes_.push_back(NodeInfo(position, node_i2u.size(), node_val));
    node_i2u.push_back(num_nodes_);
    nodes_adj.push_back(std::vector<edge_info>());
    num_nodes_++;
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (!n.valid()) {
      return false;
    }
    return this == n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < 0 || i >= num_nodes()) {
      throw std::invalid_argument("Invalid index to get the edge.");
    }
    return Node(this, node_i2u.at(i));
  }

  //
  // EDGES
  //

  /** Custom structure of information to store with Edges */
  struct EdgeInfo {
    size_type        node2_uid;       //< Node 2's uid
    edge_value_type  value;           //< Edge's value

    EdgeInfo() : node2_uid(size_type()), value(edge_value_type()) {}

    /** Construct a EdgeInfo by the given _node2_uidx_ and _edge_value.
     * @param[in] node2_uidx  Node2's uid.
     * @param[in] edge_value  Edge value.
     */
    EdgeInfo(size_type node2_uidx, const edge_value_type& edge_value) :
      node2_uid(node2_uidx), value(edge_value) {}
  };

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
      graph_loc = nullptr;
      node1_uid_ = 0;
      node2_uid_ = 0;
    }

    /** Return a node of this Edge. */
    Node node1() const {
      return Node(graph_loc, node1_uid_);
    }

    /** Return the other node of this Edge. */
    Node node2() const {
      return Node(graph_loc, node2_uid_);
    }

    /** Return the length of this Edge. */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /**
     * @brief Return the value of the Edge.
     * @return A reference to the value of _edge_.
     *
     * @pre The edge is a valid edge, i.e., its two nodes are valid.
     * @post The returned value can be modified (but only with edge (a, b),
     *       the edge value of (b, a) remains the same).
     *
     * Complexity: O(d) where d is the largest degree of a node in the graph.
     */
    edge_value_type& value() {
      if (!node1().valid() || !node2().valid()) {
        throw std::invalid_argument("The edge is not a valid edge.");
      }
      for (size_type i = 0; i < graph_loc->nodes_adj.at(node1_uid_).size();
           i++) {
        if (graph_loc->nodes_adj.at(node1_uid_).at(i).node2_uid == node2_uid_) {
          return graph_loc->nodes_adj.at(node1_uid_).at(i).value;
        }
      }
      std::string err = "Value of this edge (" + std::to_string(node1_uid_) +
                        ", " + std::to_string(node2_uid_) +
                        ") in uid is not found.";
      throw std::runtime_error(err);
    }

    /**
     * @brief Return the value of the edge.
     * @return A reference to the value of _edge_.
     *
     * @pre The edge is a valid edge, i.e., its two nodes are valid.
     * @post The returned value cannot be modified.
     *
     * Complexity: O(d) where d is the largest degree of a node in the graph.
     */
    const edge_value_type& value() const {
      if (!node1().valid() || !node2().valid()) {
        throw std::invalid_argument("The edge is not a valid edge.");
      }
      for (size_type i = 0; i < graph_loc->nodes_adj.at(node1_uid_).size();
           i++) {
        if (graph_loc->nodes_adj.at(node1_uid_).at(i).node2_uid == node2_uid_) {
          return graph_loc->nodes_adj.at(node1_uid_).at(i).value;
        }
      }
      std::string err = "Value of this edge (" + std::to_string(node1_uid_) +
                        ", " + std::to_string(node2_uid_) +
                        ") in uid is not found.";
      throw std::runtime_error(err);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph_loc == e.graph_loc &&
             ((node1_uid_ == e.node1_uid_ && node2_uid_ == e.node2_uid_) ||
              (node1_uid_ == e.node2_uid_ && node2_uid_ == e.node1_uid_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return !(*this == e) && ((node1_uid_ < e.node1_uid_) ||
              (node1_uid_ == e.node1_uid_ && node2_uid_ < e.node2_uid_) ||
              graph_loc < e.graph_loc);
    }

    /** Set operator "=" so that we can assign a valid edge
     * to an invalid edge later.
     */
    edge_type& operator=(const edge_type& edge) = default;

    /** default destructor */
    ~Edge() = default;

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // A pointer to the graph
    graph_type* graph_loc;
    // The first node uid of Edge
    size_type node1_uid_;
    // The second node uid of Edge
    size_type node2_uid_;

    /** Construct a valid edge by the given @a node1 and @a node2.
     * @param[in] node1  The first node
     * @param[in] node2  The second node
     */
    Edge(const node_type& node1, const node_type& node2) {
      graph_loc = node1.graph_;
      node1_uid_ = graph_loc->node_i2u.at(node1.index());
      node2_uid_ = graph_loc->node_i2u.at(node2.index());
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i < 0 || i >= num_edges_) {
      throw std::invalid_argument("Invalid index to get the edge.");
    }
    return edges_.at(i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (!a.valid() || !b.valid()) {
      throw std::invalid_argument("Node a or Node b is not valid node.");
    }
    for (auto a_adj_it = a.edge_begin(); a_adj_it != a.edge_end(); ++a_adj_it) {
      if ((*a_adj_it).node2() == b) {
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
  Edge add_edge(const Node& a, const Node& b,
                const edge_value_type& edge_val = edge_value_type()) {
    if (!a.valid() || !b.valid()) {
      throw std::invalid_argument("Node a or Node b is not a valid node.");
    } else if (a == b) {
      throw std::invalid_argument("Node a and Node b are not distinct nodes.");
    } else if (a.graph_ != b.graph_) {
      throw std::invalid_argument("Node a and Node b are not in the"
                                  "same graph");
    }

    edge_type edge = Edge(a, b);
    if (!has_edge(a, b)) {
      // Add a new edge.
      edges_.push_back(edge);
      num_edges_++;
      // Increase degrees for the two corresponding nodes.
      nodes_adj[a.index()].push_back(EdgeInfo(b.index(), edge_val));
      nodes_adj[b.index()].push_back(EdgeInfo(a.index(), edge_val));
    }
    return edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    num_edges_ = 0;
    num_nodes_ = 0;
    nodes_.clear();
    node_i2u.clear();
    nodes_adj.clear();
    edges_.clear();
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
    NodeIterator() {
      graph_loc = nullptr;
      index = 0;
    }

    /**
     * @brief Return the node which the node iterator points to.
     * @return _node_ which the node iterator points to.
     *
     * Complexity: O(1).
     */
    Node operator*() const {
      return Node(graph_loc, graph_loc->node_i2u.at(index));
    }

    /**
     * @brief Increment the node iterator pointing to the next node.
     * @return a reference to the current node iterator.
     *
     * Complexity: O(1).
     */
    NodeIterator& operator++() {
      index++;
      return *this;
    }

    /**
     * @brief Define the equality for two node iterators.
     * @return A boolean value which represents whether the two node iterators
     *         are equal. Return true if they are equal; false otherwise.
     * @param[in] nodeIter The other node iterator in comparison.
     *
     * Complexity: O(1).
     *
     * Check whether the two node iterators pointing at the same node
     * of the same graph.
     */
    bool operator==(const NodeIterator& nodeIter) const {
      return graph_loc == nodeIter.graph_loc && index == nodeIter.index;
    }

   private:
    friend class Graph;

    graph_type* graph_loc;
    size_type index;

    /** Construct a node iterator by the given _graph_ and _idx_.
     * @param[in] graph A pointer to the given Graph
     * @param[in] idx   A unique index for the node
     */
    NodeIterator(const graph_type* graph, size_type idx) {
      graph_loc = const_cast<graph_type*>(graph);
      index = idx;
    }

  };

  /**
   * @brief Return beginning of node iterator to the current graph.
   * @return The beginning of node iterator to _graph_.
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /**
   * @brief Return end of node iterator to the current graph.
   * @return The end of node iterator to _graph_.
   */
  node_iterator node_end() const {
    return NodeIterator(this, node_i2u.size());
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

    // An alias of const_iterator in vector.
    using inci_iter = typename std::vector<edge_info>::const_iterator;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      graph_loc = nullptr;
      node1_uid = 0;
      node2_ptr = inci_iter();
    }

    /**
     * @brief Return the edge which the incident iterator points to.
     * @return _edge_ which the incident iterator points to.
     *
     * @post The node that spawns the incident iterator is returned by node1()
     *       of _edge_ and the adjacent node was returned by node2(). The
     *       returned _edge_ is valid.
     *
     * Complexity: O(1).
     */
    Edge operator*() const {
      return Edge(Node(graph_loc, node1_uid),
                  Node(graph_loc, (*node2_ptr).node2_uid) );
    }

    /**
     * @brief Increment the incident iterator pointing to the next valid edge.
     * @return a reference to the current incident iterator.
     *
     * Complexity: O(1).
     */
    IncidentIterator& operator++() {
      ++node2_ptr;
      // Keep incrementing until finding the next valid adjacent node.
      while (node2_ptr != graph_loc->nodes_adj.at(node1_uid).end() &&
             !(Node(graph_loc, (*node2_ptr).node2_uid)).valid()) {
        ++node2_ptr;
      }
      return *this;
    }

    /**
     * @brief Define the equality for two incident iterators.
     * @return A boolean value which represents whether the two incident
     *         iterators are equal. Return true if they are equal;
     *         false otherwise.
     * @param[in] inciIter The other incident iterator in comparison.
     *
     * Complexity: O(1).
     *
     * Check whether the two incident iterators are spawned by the same node
     * and are pointing at the same adjacent node.
     */
    bool operator==(const IncidentIterator& inciIter) const {
      return graph_loc == inciIter.graph_loc && node1_uid == inciIter.node1_uid
             && node2_ptr == inciIter.node2_ptr;
    }

   private:
    friend class Graph;

    graph_type* graph_loc;
    size_type node1_uid;
    inci_iter node2_ptr;

    /** Construct an incident iterator by the given _node_.
     * @param[in] graph          A graph where the nodes are in.
     * @param[in] node1_uid      The uid of a node which spawns the incident
     *                           iterator.
     * @param[in] incident_ptr   A constant vector iterator pointing to
     *                           the adjacent nodes of _node_.
     */
    IncidentIterator(const graph_type* graph, size_type node_uid,
                     const inci_iter incident_ptr) {
      graph_loc = const_cast<graph_type*>(graph);
      node1_uid = node_uid;
      node2_ptr = incident_ptr;
      // Keep incrementing until finding the next valid adjacent node.
      while (node2_ptr != graph_loc->nodes_adj.at(node1_uid).end() &&
             !(Node(graph_loc, (*node2_ptr).node2_uid)).valid()) {
        ++node2_ptr;
      }
    }
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

    // An alias of const_iterator in vector.
    using edge_iter = typename std::vector<edge_type>::const_iterator;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      graph_loc = nullptr;
      edge_ptr = edge_iter();
    }

    /**
     * @brief Return the edge which the edge iterator points to.
     * @return _edge_ which the edge iterator points to.
     *
     * Complexity: O(1).
     */
    Edge operator*() const {
      return *edge_ptr;
    }

    /**
     * @brief Increment the edge iterator pointing to the next edge.
     * @return a reference to the current edge iterator.
     *
     * @post The _edge_ pointed by edge iterator is a valid edge, which
     *       means the two nodes as its endpoints are valid.
     *
     * Complexity: O(1).
     */
    EdgeIterator& operator++() {
      ++edge_ptr;
      // Keep incrementing until finding the next valid edge.
      while (edge_ptr != graph_loc->edges_.end() &&
             !graph_loc->has_edge((*edge_ptr).node1(), (*edge_ptr).node2())) {
        ++edge_ptr;
      }

      return *this;
    }

    /**
     * @brief Define the equality for two edge iterators.
     * @return A boolean value which represents whether the two edge
     *         iterators are equal. Return true if they are equal;
     *         false otherwise.
     * @param[in] edgeIter The other edge iterator in comparison.
     *
     * Complexity: O(1).
     *
     * Check whether the two edge iterators are spawned by the same node
     * and are pointing at the same adjacent node.
     */
    bool operator==(const EdgeIterator& edgeIter) const {
      return graph_loc == edgeIter.graph_loc && edge_ptr == edgeIter.edge_ptr;
    }

   private:
    friend class Graph;

    graph_type* graph_loc;
    edge_iter edge_ptr;

    /** Construct an edge iterator by the given root _node_.
     * @param[in] graph          A graph where the nodes are in.
     * @param[in] edge_iter1     An edge iterator which visits all nodes in
     *                           the graph.
     */
    EdgeIterator(const graph_type* graph, edge_iter edge_iter1) {
      graph_loc = const_cast<graph_type*>(graph);
      edge_ptr = edge_iter1;
      // Move the iterator pointing to the next valid edge or the "end" edge.
      while (edge_ptr != graph_loc->edges_.end() &&
             !graph_loc->has_edge((*edge_ptr).node1(), (*edge_ptr).node2())) {
        ++edge_ptr;
      }
    }
  };

  /**
   * @brief Return beginning of edge iterator to the current graph.
   * @return The beginning of edge iterator to _graph_.
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edges_.begin());
  }

  /**
   * @brief Return end of edge iterator to the current graph.
   * @return The end of edge iterator to _graph_.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.end());
  }

  /** Remove a Node _node_ from the graph.
   * @return An unsigned value (0 or 1) where 1 implies the node is removed
   *         and 0 implies the node is not removed.
   * @param[in] node        A node to be removed.
   * @post  If the graph contains this (valid) node, Node _node_ is
   *        invalidated in the graph with no neighbors. Although
   *        its information is kept in node_, an index-uid vector for nodes
   *        denoted as node_i2u will not assign an index to _node_.
   *        All edges incident to _node_ will be invalidated
   *        in the graph (removed from corresponding data structure that
   *        stored the edges).
   * @post  node_i2u is still consecutive in index. For all nodes n in node_i2u,
   *        they are valid, i.e.,
   *          1. uid in range: 0 <= n.uid < num_nodes_,
   *          2. idx in range: 0 <= n.idx < node_i2u.size(),
   *          3. uid in sync:  node_i2u.at(n.idx) = n.uid,
   *        where num_nodes_ denotes the number of nodes ever added to the
   *        graph.
   *
   * Complexity: O(d * (m+n)) where d is the largest degree of a node in
   *             the graph (O(n * (n + m)) as d < n) where
   *             G = (V, E) with |V| = n and |E| = m.
   */
  size_type remove_node(const Node& node) {
    if (has_node(node)) {
      // Remove edges incident to the node.
      for (unsigned i = 0; i < nodes_adj.at(node.n_uid_).size(); ) {
        auto front_uid = nodes_adj.at(node.n_uid_).at(0).node2_uid;
        remove_edge(node, Node(node.graph_, front_uid));
      }
      // Update node_i2u by swapping the last element's uid and _node_'s uid
      // and delete the latter in node_i2u.
      if (node.index() < node_i2u.size() - 1) {
        auto currIndex = node.index();
        auto tempUid = std::move(node_i2u.at(currIndex));
        node_i2u.at(currIndex) = std::move(node_i2u.back());
        node_i2u.back() = std::move(tempUid);

        // Update the swapped last element's new index in nodes_'s index.
        nodes_.at(node_i2u.at(currIndex)).index = currIndex;
      }
      node_i2u.pop_back();
      return 1;
    }
    return 0;
  }

  /** Remove a Node _node_ from the graph by a Node iterator _n_it_.
   * @return A Node iterator pointing to the closest Node. If _node_ is
   *         removed, we have three conditions: if _node_ is the
   *         last node, then pointing to node_end(); if _node_ is the final
   *         node of node set V (Graph G = (V, E) and |V| > 1), pointing to
   *         its previous element; if _node_ is a normal element, pointing
   *         to the newly swapped node at the same index. Otherwise,
   *         increment _n_it_.
   * @param[in] n_it      A node iterator pointing to the node to be removed.
   *
   * @post  If the graph contains this (valid) node, Node _node_ is
   *        invalidated in the graph with no neighbors. Although
   *        its information is kept in node_, an index-uid vector for nodes
   *        denoted as node_i2u will not assign an index to _node_.
   *        All edges incident to _node_ will be invalidated
   *        in the graph (removed from corresponding data structure that
   *        stored the edges).
   * @post  node_i2u is still consecutive in index. For all nodes n in node_i2u,
   *        they are valid, i.e.,
   *          1. uid in range: 0 <= n.uid < num_nodes_,
   *          2. idx in range: 0 <= n.idx < node_i2u.size(),
   *          3. uid in sync:  node_i2u.at(n.idx) = n.uid,
   *        where num_nodes_ denotes the number of nodes ever added to the
   *        graph.
   *
   * Complexity: O(d * (m+n)) where d is the largest degree of a node in
   *             the graph (O(n * (n + m)) as d < n) where
   *             G = (V, E) with |V| = n and |E| = m.
   */
  node_iterator remove_node(node_iterator n_it) {
    auto diff_index = n_it.index - 0;
    size_type isRemoved = remove_node(*n_it);
    if (!isRemoved) {
      // No node is removed.
      return ++n_it;
    } else {
      if (node_i2u.size() == 0) {
        // The removed node is the last node (|V| = 1) in the graph.
        return this->node_end();
      } else if (diff_index == node_i2u.size()) {
        // The removed node is the final element of node set V (|V| > 1)
        // in the graph (G = (V, E)).
        return NodeIterator(this, diff_index - 1);
      } else {
        // The removed node is not the final element of node set.
        return NodeIterator(this, diff_index);
      }
    }
  }

  /** Remove an Edge _edge_ from the graph by its two nodes _node1_ and _node2_.
   * @return An unsigned value (0 or 1) where 1 implies the edge is removed
   *         and 0 implies the edge is not removed.
   * @param[in] node1        Node1 of the edge to be removed.
   * @param[in] node2        Node2 of the edge to be removed.
   * @post  If the graph contains this (valid) edge, Edge _edge_ is
   *        invalidated in the graph, i.e., removed from the adjacency
   *        list nodes_adj and index list edges_.
   * @post  edges_ is still consecutive in index. For all edges in edges_,
   *        they are valid (i.e., the two nodes of an edge are both valid).
   *        same property holds for nodes_adj.
   *
   * Complexity: O(m + n) where G = (V, E) with |V| = n and |E| = m.
   */
  size_type remove_edge(const Node& node1, const Node& node2) {
    // Delete the corresponding edge information in nodes_adj.
    if (has_edge(node1, node2)) {
      remove_edge_helper(node1.n_uid_, node2.n_uid_);
      remove_edge_helper(node2.n_uid_, node1.n_uid_);
      // Delete the corresponding edge information in edges_.
      remove_edge_in_edges_(node1.n_uid_, node2.n_uid_);
      num_edges_--;
      return 1;
    }
    return 0;
  }

  /** Remove an Edge _edge_ from the graph.
   * @return An unsigned value (0 or 1) where 1 implies the edge is removed
   *         and 0 implies the edge is not removed.
   * @param[in] edge        The Edge to be removed.
   * @post  If the graph contains this (valid) edge, Edge _edge_ is
   *        invalidated in the graph, i.e., removed from the adjacency
   *        list nodes_adj and index list edges_.
   * @post  edges_ is still consecutive in index. For all edges in edges_,
   *        they are valid (i.e., the two nodes of an edge are both valid).
   *        same property holds for nodes_adj.
   *
   * Complexity: O(m + n) where G = (V, E) with |V| = n and |E| = m.
   */
  size_type remove_edge(const Edge& edge) {
    return remove_edge(edge.node1(), edge.node2());
  }

  /** Remove an Edge _edge_ from the graph by an edge iterator _e_it_.
   * @return An Edge iterator pointing to the closest Edge. If _edge_ is
   *         removed, we have three conditions: if _edge_ is the
   *         last edge, then pointing to edge_end(); if _edge_ is the final
   *         edge of edge set E (Graph G = (V, E) and |E| > 1), pointing to
   *         its previous element; if _edge_ is a normal element, pointing
   *         to the newly swapped edge at the same index. Otherwise,
   *         increment _e_it_.
   * @param[in] e_it_     An edge iterator pointing to the edge to be removed.
   * @post  If the graph contains this (valid) edge, Edge _edge_ is
   *        invalidated in the graph, i.e., removed from the adjacency
   *        list nodes_adj and index list edges_.
   * @post  edges_ is still consecutive in index. For all edges in edges_,
   *        they are valid (i.e., the two nodes of an edge are both valid).
   *        same property holds for nodes_adj.
   *
   * Complexity: O(m + n) where G = (V, E) with |V| = n and |E| = m.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    auto diff = e_it.edge_ptr - edges_.begin();
    size_type isRemoved = remove_edge(*e_it);
    if (!isRemoved) {
      return ++e_it;
    } else {
      if (edges_.size() == 0) {
        // The removed edge is the last edge (|E| = 1) in the graph.
        return this->edge_end();
      } else if (diff == edges_.end() - edges_.begin()) {
        // The removed edge is the final element of edge set E (|E| > 1)
        // in the graph (G = (V, E)).
        return EdgeIterator(this, edges_.begin() + diff - 1);
      } else {
        // The removed edge is not the final element of edge set.
        return EdgeIterator(this, edges_.begin() + diff);
      }
    }
  }


 private:
  // The number of nodes ever added into the graph.
  size_type num_nodes_;
  // A vector that keeps track of node uid by its node index.
  std::vector<size_type> node_i2u;
  // A vector that keeps track of node value and position.
  std::vector<node_info> nodes_;
  // A vector that keeps track of a node's incident edge's information.
  std::vector<std::vector<edge_info>> nodes_adj;
  // A vector that keeps track of the edges.
  std::vector<edge_type> edges_;
  // The number of edges in the graph.
  size_type num_edges_;

  /** A helper function to remove an Edge _edge_ from the adjacency list
   * nodes_adj by the indices of its two nodes.
   * @param[in] n1_uid    Node1 index of the edge.
   * @param[in] n2_uid    Node2 index of the edge.
   * @post  For all edges in nodes_adj, they are valid (i.e., the two nodes
   *        of an edge are both valid).
   *
   * Complexity: O(d) where d is the largest degree of a node in the graph.
   */
  void remove_edge_helper(size_type n1_uid, size_type n2_uid) {
    for (unsigned i = 0; i < nodes_adj.at(n1_uid).size(); i++) {
      if (nodes_adj.at(n1_uid).at(i).node2_uid == n2_uid) {
        if (i < nodes_adj.at(n1_uid).size() - 1) {
          auto temp = std::move(nodes_adj.at(n1_uid).at(i));
          nodes_adj.at(n1_uid).at(i) = std::move(nodes_adj.at(n1_uid).back());
          nodes_adj.at(n1_uid).back() = std::move(temp);
        }
        nodes_adj.at(n1_uid).pop_back();
        break;
      }
    }
  }

  /** A helper function to remove an Edge _edge_ from the index list of edges
   * edges_ by the indices of its two nodes.
   * @param[in] n1_uid    Node1 index of the edge.
   * @param[in] n2_uid    Node2 index of the edge.
   * @post  The indices of edges_ is still consecutive. For all edges
   *        in edges_, they are valid (i.e., the two nodes of an edge are
   *        both valid).
   *
   * Complexity: O(m) where G = (V, E) and |E| = m.
   */
  void remove_edge_in_edges_(size_type n1_uid, size_type n2_uid) {
    for (unsigned i = 0; i < edges_.size(); i++) {
      auto edge = edges_.at(i);
      if ((edge.node1().n_uid_ == n1_uid && edge.node2().n_uid_ == n2_uid) ||
          (edge.node2().n_uid_ == n1_uid && edge.node1().n_uid_ == n2_uid)) {
        if (i < edges_.size() - 1) {
          auto temp_edge = std::move(edges_.at(i));
          edges_.at(i) = std::move(edges_.back());
          edges_.back() = std::move(temp_edge);
        }
        edges_.pop_back();
        break;
      }
    }
  }

};

#endif // CME212_GRAPH_HPP
