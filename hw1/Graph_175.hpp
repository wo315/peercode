#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
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
template <typename V = int>
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
  Graph() : nodes_pt(), num_nodes_(0), nodes_val(), nodes_adj(),
            edges_(), edges_val(), num_edges_(0) {}

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
      graph_ = nullptr;
      n_idx_ = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_pt.at(n_idx_);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return n_idx_;
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
      return graph_->nodes_val.at(n_idx_);
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
      return graph_->nodes_val.at(n_idx_);
    }

    /**
     * @brief Return the degree of the node.
     * @return The degree of _node_.
     *
     * @pre The node is a valid node.
     *
     * Complexity: O(1).
     */
    size_type degree() const {
      return graph_->nodes_adj[*this].size();
    }

    /**
     * @brief Return beginning of incident iterator to the current node.
     * @return The beginning of incident iterator to _node_.
     *
     * @pre The node is a valid node.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(*this, graph_->nodes_adj[*this].begin());
    }

    /**
     * @brief Return end of incident iterator to the current node.
     * @return The end of incident iterator to _node_.
     *
     * @pre The node is a valid node.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(*this, graph_->nodes_adj[*this].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return graph_ == n.graph_ && n_idx_ == n.n_idx_;
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
      return n_idx_ < n.n_idx_;
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
    size_type n_idx_;

    /** Construct a valid node by the given _graph_ and _idx_.
     * @param[in] graph A pointer to the given Graph
     * @param[in] idx   A unique index for the node
     */
    Node(const graph_type* graph, size_type idx) {
      graph_ = const_cast<graph_type*>(graph);
      n_idx_ = idx;
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_nodes_;
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
    nodes_pt.push_back(position);
    nodes_val.push_back(node_val);
    nodes_adj[new_node] = std::vector<node_type>();
    num_nodes_++;
    return new_node;
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      node1_ = Node();
      node2_ = Node();
    }

    /** Return a node of this Edge */
    Node node1() const {
      return node1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1_ == e.node1_ && node2_ == e.node2_) ||
             (node1_ == e.node2_ && node2_ == e.node1_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (node1_ < e.node1_) || (node1_ == e.node1_ && node2_ < e.node2_);
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

    // The first node of Edge
    node_type node1_;
    // The second node of Edge
    node_type node2_;

    /** Construct a valid edge by the given @a node1 and @a node2.
     * @param[in] node1 The first node
     *            node2 The second node
     */
    Edge(const node_type& node1, const node_type& node2) {
      node1_ = const_cast<node_type&>(node1);
      node2_ = const_cast<node_type&>(node2);
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
    return edges_val.at(i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return edges_.find(Edge(a, b)) != edges_.end() ||
           edges_.find(Edge(b, a)) != edges_.end();
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
    edge_type edge = Edge(a, b);
    if (!has_edge(a, b)) {
      // Add a new edge.
      edges_[edge] = num_edges_;
      edges_val.push_back(edge);
      num_edges_++;
      // Increase degrees for the two corresponding nodes.
      nodes_adj[a].push_back(b);
      nodes_adj[b].push_back(a);

      return edge;
    } else {
      // Find an old edge.
      auto tempEdgePair = edges_.find(edge);
      if (tempEdgePair != edges_.end()) {
        return tempEdgePair->first;
      } else {
        return edges_.find(Edge(b, a))->first;
      }
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges_.clear();
    num_edges_ = 0;
    nodes_pt.clear();
    num_nodes_ = 0;
    nodes_val.clear();
    nodes_adj.clear();
    edges_val.clear();
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
      return Node(graph_loc, index);
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
    return NodeIterator(this, num_nodes_);
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
    using inci_iter = typename std::vector<node_type>::const_iterator;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      node1 = Node();
      node2_ptr = inci_iter();
    }

    /**
     * @brief Return the edge which the incident iterator points to.
     * @return _edge_ which the incident iterator points to.
     *
     * @post The node that spawns the incident iterator is returned by node1()
     *       of _edge_ and the adjacent node was returned by node2().
     *
     * Complexity: O(1).
     */
    Edge operator*() const {
      return Edge(node1, *node2_ptr);
    }

    /**
     * @brief Increment the incident iterator pointing to the next edge.
     * @return a reference to the current incident iterator.
     *
     * Complexity: O(1).
     */
    IncidentIterator& operator++() {
      node2_ptr++;
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
      return node1 == inciIter.node1 && node2_ptr == inciIter.node2_ptr;
    }

   private:
    friend class Graph;

    node_type node1;
    inci_iter node2_ptr;

    /** Construct an incident iterator by the given _node_.
     * @param[in] node           A node which spawns the incident iterator.
     * @param[in] incident_ptr   A constant vector iterator pointing to
     *                           the adjacent nodes of _node_.
     */
    IncidentIterator(const node_type& node, const inci_iter incident_ptr) {
      node1 = node;
      node2_ptr = incident_ptr;
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

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      node_ptr = NodeIterator();
      inci_ptr = IncidentIterator();
    }

    /**
     * @brief Return the edge which the edge iterator points to.
     * @return _edge_ which the edge iterator points to.
     *
     * Complexity: O(1).
     */
    Edge operator*() const {
      return *inci_ptr;
    }

    /**
     * @brief Increment the edge iterator pointing to the next edge.
     * @return a reference to the current edge iterator.
     *
     * @post The _edge_ pointed by edge iterator is either an "end" edge
     *       which is defined to be the end of the last node's adjacency
     *       node list or a valid edge in the graph.
     * @post The _edge_ pointed by edge iterator has node1() < node2()
     *
     * Complexity: O(1).
     *
     * We only count the edge which has node1() < node2() as a candidate
     * for next edge so that each edge in the undirected graph is visited
     * once.
     */
    EdgeIterator& operator++() {
      incrementInPlaceHelper(inci_ptr);
      // Before reaching the last node of the adjacency list and its last
      // adjacent node, if node1() of the edge is larger than node2(),
      // we will continue pointing to the next edge (i.e. {node, adjacent node}
      // pair) so that each edge is visited once.
      while (!((inci_ptr) == (*node_ptr).edge_end() &&
               (*node_ptr).index() == node_num - 1) &&
             (*inci_ptr).node1() > (*inci_ptr).node2()) {
        incrementInPlaceHelper(inci_ptr);
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
      return inci_ptr == edgeIter.inci_ptr;
    }

   private:
    friend class Graph;

    NodeIterator node_ptr;
    IncidentIterator inci_ptr;
    size_type node_num;

    /** Construct an edge iterator by the given root _node_.
     * @param[in] node_iter   A node iterator which visits all nodes in
     *                        the graph.
     * @param[in] inci_iter   A incident iterator which visits all incident
     *                        edges to a node.
     * @param[in] num_node    The number of nodes in the graph.
     */
    EdgeIterator(const NodeIterator node_iter, const IncidentIterator inci_iter,
                 size_type num_node) {
      node_ptr = node_iter;
      inci_ptr = inci_iter;
      node_num = num_node;
      // Move the iterator pointing to a valid edge or the "end" edge.
      while ((inci_ptr) == (*node_ptr).edge_end() &&
             (*node_ptr).index() < node_num - 1) {
        // We do not check whether the edge is valid in the constructor
        // because the first edge is {u, v} with u < v. Consider that
        // node_ptr visits by index++, node comparison uses node index.
        // If the first valid edge is {v, u} with v > u, then
        // it leads to a contradiction because we check u by node_ptr
        // before v, so we must find {u, v} before finding {v, u}.
        ++node_ptr;
        inci_ptr = (*node_ptr).edge_begin();
      }
    }

    /**
     * @brief Increment the edge iterator to the next edge.
     * @param[in,out] incident_ptr   A incident iterator which visits all
     *                               incident edges to a node.
     *
     * @post The _edge_ pointed by _incident_ptr_ is either an "end" edge
     *       which is defined to be the end of the last node's adjacency
     *       node list or a valid edge in the graph.
     *
     * The next edge is not necessarily a "valid" edge for edge iterator
     * which satisfies the property that edge has node1() < node2()).
     */
    void incrementInPlaceHelper(IncidentIterator& incident_ptr) {
      IncidentIterator inci_ptr_temp = incident_ptr;
      if ((++inci_ptr_temp) == (*node_ptr).edge_end() &&
          (*node_ptr).index() < node_num - 1) {
        // The iterator reaches the end of adjacency node list of
        // a not-last node.
        ++node_ptr;
        incident_ptr = (*node_ptr).edge_begin();
      } else if ((*node_ptr).index() <= node_num - 1) {
        // The iterator reaches an existed edge in the adjacency node list.
        ++incident_ptr;
      }
    }
  };

  /**
   * @brief Return beginning of edge iterator to the current graph.
   * @return The beginning of edge iterator to _graph_.
   */
  edge_iterator edge_begin() const {
    node_iterator n0_iter = node_begin();
    return EdgeIterator(n0_iter, (*n0_iter).edge_begin(), num_nodes_);
  }

  /**
   * @brief Return end of edge iterator to the current graph.
   * @return The end of edge iterator to _graph_.
   */
  edge_iterator edge_end() const {
    node_iterator nf_iter = NodeIterator(this, num_nodes_ - 1);
    return EdgeIterator(nf_iter, (*nf_iter).edge_end(), num_nodes_);
  }


 private:

  // A vector that keeps track of Point value of
  // the node
  std::vector<Point> nodes_pt;
  // The number of nodes in the graph
  size_type num_nodes_;
  // A vector that keeps track of node value
  std::vector<node_value_type> nodes_val;
  // A ordered map with keys that keeps track of the node and
  // values in a vector that keeps track of its adjacent nodes
  std::map<node_type, std::vector<node_type>> nodes_adj;
  // An ordered map with keys that keep track of the edge
  // and values that keep track of the edge index
  std::map<edge_type, size_type> edges_;
  // A vector that keeps track of the edges
  std::vector<edge_type> edges_val;
  // The number of edges in the graph
  size_type num_edges_;

};

#endif // CME212_GRAPH_HPP
