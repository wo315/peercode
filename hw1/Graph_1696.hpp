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
template  <typename V>
class Graph {
//  private:

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  using node_type = Node; // Synonym for Node (following STL conventions).
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  using edge_type = Edge; // Synonym for Edge (following STL conventions).

  class NodeIterator;                 // Iterator over all graph nodes.
  using node_iterator = NodeIterator; // Synonym for NodeIterator

  class EdgeIterator;                 // Iterator over all graph edges.
  using edge_iterator = EdgeIterator; // Synonym for EdgeIterator

  class IncidentIterator;             // Iterator over node incident edges.
  using incident_iterator = IncidentIterator; // Synonym for IncidentIterator

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;


  /* CONSTRUCTORS AND DESTRUCTOR */

  /** Construct an empty graph. [HW0] */
  Graph() {}

  /** Default destructor */
  ~Graph() = default;

  
  /* NODES */

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

    /** Return this node's position. [HW0] */
    const Point& position() const {
      return graph_->nodes_[index()].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). [HW0] */
    size_type index() const {
      return index_;
    }

    /* Incident edge iterator [HW1] */

    /** Return a reference to this node's value. [HW1] */
    node_value_type& value() {
      return graph_->nodes_[index()].value;
    }

    /** Return a reference to this node's value (read-only). [HW1] */
    const node_value_type& value() const {
      return graph_->nodes_[index()].value;
    }

    /** Return the number of incident edges. [HW1] */
    size_type degree() const {
      return graph_->nodes_[index()].inc_edges.size();
    }

    /** Return an iterator pointing to the first incident edge. [HW1] */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, index(), 0);
    }

    /** Return an iterator pointing past the last incident edge. [HW1] */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, index(), graph_->nodes_[index()].inc_edges.size());
    }

    /** Test whether this node and @a n are equal. [HW0]
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return graph_ == n.graph_ && index() == n.index();
    }

    /** Test whether this node is less than @a n in a global order. [HW0]
     * 
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     * 
     * Defined for nodes of the same graph (Piazza @75).
     */
    bool operator<(const Node& n) const {
      return index() < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type index_;

    /** Private constructor accessed by Graph. [HW0] */
    Node(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)),
        index_(index) {}
  };

  /** Return the number of nodes in the graph. [HW0]
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

  /** Add a node to the graph, returning the added node. [Hw0]
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations. [HW0]
   */
  Node add_node(const Point& position) {
    size_type next_index = size();
    nodes_.emplace_back(position);
    return Node(this, next_index);
  }

  /** Add a node to the graph, returning the added node. [HW1]
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * 
   */
  Node add_node(const Point& position, const node_value_type& value) {
    size_type next_index = size();
    nodes_.emplace_back(position, value);
    return Node(this, next_index);
  }

  /** Determine if a Node belongs to this Graph. [HW0]
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return Node(this, n.index()) == n;
  }

  /** Return the node with index @a i. [HW0]
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, i);
  }


  /* EDGES */

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
    Node node1() const {
      return graph_->edges_[index_].node1;
    }

    /** Return the other node of this Edge [HW0] */
    Node node2() const {
      return graph_->edges_[index_].node2;
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
     * Defined for edges of the same graph (Piazza @75).
     */
    bool operator<(const Edge& e) const {
      return index_ < e.index_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type index_;

    /** Private constructor accessed by Graph. [HW0] */
    Edge(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)),
        index_(index) {}
  };

  /** Return the total number of edges in the graph. [HW0]
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i. [HW0]
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge. [HW0]
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return index(a, b) > -1;
  }

  /** Add an edge to the graph, or return the current edge if it already exists. [HW0]
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
    int edge_index = index(a, b);
    if (edge_index > -1) return edge(edge_index, a, b);   // Edge already exists

    size_type next_index = num_edges();
    edges_.emplace_back(a, b);
    edge_map_[std::make_pair(a.index(), b.index())] = next_index;

    nodes_[a.index()].inc_edges.push_back(next_index);
    nodes_[b.index()].inc_edges.push_back(next_index);

    return edge(next_index);
  }

  /** Remove all nodes and edges from this graph. [HW0]
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    edge_map_.clear();
  }


  /* Node iterator */

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
    NodeIterator() {}

    /** Dereference operator. [HW1] */
    Node operator*() const {
      return graph_->node(index_);
    }

    /** Increments (prefix) to the next node in the graph. [Hw1] */
    NodeIterator& operator++() {
      ++index_;
      return *this;
    }

    /** Defines equality between two iterators. [HW1] */
    bool operator==(const NodeIterator& iter) const {
      return (graph_ == iter.graph_) && (index_ == iter.index_);
    }

   private:
    // Allow Graph to access NodeIterator's private member data and functions.
    friend class Graph;

    Graph* graph_;      // Pointer to graph.
    size_type index_;   // Node index.

    /** Private constructor accessed by Graph. [HW1] */
    NodeIterator(const Graph* graph, size_type index) 
        : graph_(const_cast<Graph*>(graph)), 
          index_(index) {}
  };

  /** Return an iterator pointing to the first node in the graph. [Hw1] */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return an iterator pointing past the last node in the graph. [Hw1] */
  node_iterator node_end() const {
    return NodeIterator(this, nodes_.size());
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
    IncidentIterator() {}

    /** Dereference operator. [HW1] 
     * 
     * @return An edge incident to node1().
     */
    Edge operator*() const {
      int edge_index = graph_->nodes_[nodei_].inc_edges[index_];
      Node n1 = graph_->edges_[edge_index].node1;
      Node n2 = graph_->edges_[edge_index].node2;

      if (nodei_ == n1.index()) {
        return graph_->edge(edge_index, n1, n2);
      }

      return graph_->edge(edge_index, n2, n1);
    }

    /** Increments (prefix) to the next incident edge. [HW1] */
    IncidentIterator& operator++() {
      ++index_;
      return *this;
    }

    /** Defines equality between two iterators. [HW1] */
    bool operator==(const IncidentIterator& iter) const {
      return (graph_ == iter.graph_) &&
          (nodei_ == iter.nodei_) &&
          (index_ == iter.index_);
    }

   private:
    // Allow Graph to access IncidentIterator's private member data and functions.
    friend class Graph;

    Graph* graph_;      // Pointer to graph.
    size_type nodei_;   // Node index.
    size_type index_;   // Incident edge index.

    /** Private constructor accessed by Graph. [HW1] */
    IncidentIterator(const Graph* graph, size_type nodei, size_type index) 
        : graph_(const_cast<Graph*>(graph)), 
          nodei_(nodei),
          index_(index) {}
  };


  /* Edge iterator */

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
    EdgeIterator() {}

    /** Dereference operator. [HW1] */
    Edge operator*() const {
      return graph_->edge(index_);
    }

    /** Increments (prefix) to the next edge in the Graph. [HW1] */
    EdgeIterator& operator++() {
      ++index_;
      return *this;
    }

    /** Defines equality between two iterators. [HW1] */
    bool operator==(const EdgeIterator& iter) const {
      return (graph_ == iter.graph_) && (index_ == iter.index_);
    }

   private:
    // Allow Graph to access EdgeIterator's private member data and functions.
    friend class Graph;
    Graph* graph_;      // Pointer to graph.
    size_type index_;   // Edge index.

    /** Private constructor accessed by Graph. [HW1] */
    EdgeIterator(const Graph* graph, size_type index) 
        : graph_(const_cast<Graph*>(graph)), 
          index_(index) {}
  };

  /** Return an iterator pointing to the first edge in the graph. [HW1] */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Return an iterator pointing past the last edge in the graph. [HW1] */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.size());
  }

 private:
  /* Graph internals */
  struct NodeData {
    Point position;
    std::vector<int> inc_edges;
    node_value_type value;
    NodeData(Point p) : position(p) {}
    NodeData(Point p, node_value_type v) : position(p), value(v) {}
  };

  struct EdgeData {
    Node node1;
    Node node2;
    EdgeData(Node n1, Node n2) : node1(n1), node2(n2) {}
  };

  std::vector<NodeData> nodes_;
  std::vector<EdgeData> edges_;
  std::map<std::pair<size_type, size_type>, size_type> edge_map_;

  /** Return the index of the edge connecting two nodes. [HW0]
   * @pre @a a and @a b are valid nodes of this graph
   * @return An edge index in the range [0, num_edges) if it exists, else -1.
   */
  int index(const Node& a, const Node& b) const {
    // search for (a, b)
    auto search1 = edge_map_.find(std::make_pair(a.index(), b.index()));
    if (search1 != edge_map_.end()) return search1->second; 

    // search for (b, a)
    auto search2 = edge_map_.find(std::make_pair(b.index(), a.index()));
    if (search2 != edge_map_.end()) return search2->second;

    // not found
    return -1;
  }


  /** Return the edge with index @a i and nodes @a a and @a b. [HW0]
   * @pre 0 <= @a i < num_edges()
   * @post result_edge.node1() == a
   * @post result_edge.node2() == b
   * 
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   */
  Edge edge(size_type i, const Node& a, const Node& b) {
    edges_[i].node1 = a;
    edges_[i].node2 = b;
    return Edge(this, i);
  }

};

#endif // CME212_GRAPH_HPP
