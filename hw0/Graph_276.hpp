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
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
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
  class Node {
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

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && uid_ == n.uid_);
    }

    /** Test whether this node and @a n are unequal.
     */
    bool operator!=(const Node& n) const {
      return ! (*this == n);
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
      // Dictionary order in terms of graph location and then index
      return graph_ != n.graph_
          ? graph_ < n.graph_ 
          : uid_ < n.uid_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    size_type uid_; // This element's unique identifier.
    Graph* graph_;  // Pointer to parent graph.

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
  Node add_node(const Point& position) {
    size_type i = num_nodes();
    // Data for newest node goes to back of vector.
    nodes_.push_back(InternalNode_(i, this, position)); 
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
  Node node(size_type i) const {
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge. By default, returns the node with the
     * smaller index, unless `reverse_` is set. */
    Node node1() const {
      // Use parent graph `node()` method to get node by index.
      return reverse_
        ? graph_->node(fetch().node2)
        : graph_->node(fetch().node1);
    }

    /** Return the other node of this Edge. By default, returns the node with
     * the larger index, unless `reverse_` is set. */
    Node node2() const {
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

    size_type uid_;     // This element's unique identifier.
    Graph* graph_;      // Pointer to parent graph.
    bool reverse_;     // Reverse order nodes are output if set.

    /** Private Constructor */
    Edge(size_type uid, const graph_type* graph, const bool reverse = false)
      : uid_(uid), graph_(const_cast<graph_type*>(graph)), reverse_(reverse) {
    }

    /** Helper method to return the corresponding InternalEdge_ */
    internal_edge_type fetch() const {
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
  Edge edge(size_type i) const {
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
    for(std::vector<InternalEdge_>::const_iterator it = edges_.begin(); it !=
        edges_.end(); ++it) {
      if (it->node1 == min_index && it->node2 == max_index) {
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
  Edge add_edge(const Node& a, const Node& b) {
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
    std::vector<size_type>::const_iterator j;
    for (j = a.fetch().incident_edges.begin(); j !=
        a.fetch().incident_edges.end(); ++j) {
      if (edges_[*j].node1 == min_index && edges_[*j].node2 == max_index) {
        return Edge(*j, this, reverse);
      }
    }

    // If edge does not exist, push a new edge to the end of the edge container,
    // and add that it is incident to both nodes.
    edges_.push_back(InternalEdge_(i, this, a.index(), b.index()));
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:

  //
  // Internal Node
  //

  /** @class Graph::InternalNode_
   * @brief Class containing node internals. Proxied by Graph::Node */
  struct InternalNode_ {
    size_type uid;                          // The unique identifier for a node.
    Graph* graph_;                          // Pointer to the containing graph.
    Point position;                         // The position of the node.
    std::vector<size_type> incident_edges;  // Indices of incident edges.

    /** Constructor for InternalNode_ */
    InternalNode_(const size_type uid, const Graph* graph, const Point position)
      : uid(uid), graph_(const_cast<Graph*>(graph)), position(position),
      incident_edges() {
    }
  };

  //
  // Internal Edge
  //

  /** @class Graph::InternalEdge_
   * @brief Class containng edge internals. Proxied by Graph::Edge */
  struct InternalEdge_ {

    // NOTE: InternalEdges_ represent edges uniquely by always putting the node
    // with the smaller index in `node1`. The "orientation" of an edge is
    // properly a property of the Edge proxy, rather than the underlying data,
    // represented by the Edge::reverse_ member.

    // Unique identifier of edge.
    size_type uid;
    // Pointer to the parent graph.
    Graph* graph;
    // The uid of the first node in the edge. This is the node with the smaller index.
    size_type node1;
    // The uid of the second node in the edge. This is the node with the larger index.
    size_type node2;

    /** Constructor for InternalNode_ */
    InternalEdge_(const size_type uid, const Graph* graph, const size_type
        node1, const size_type node2) : 
      uid(uid), graph(const_cast<Graph*>(graph)) {
        this->node1 = std::min(node1, node2);
        this->node2 = std::max(node1, node2);
    }

    /* Test whether this edge is less than @a e in a global order.
     *
     * This order function allows us to use std::set<> to enable fast lookup of
     * whether an edge exists based on the nodes it joins.
     */
    bool operator<(const InternalEdge_& e) const {
      // Dictionary order in terms of the parent graph, the smaller node index, and then the
      // larger node index. Since edges have a unique representation as @a
      // InternalEdge_, this makes the definition independent of the order the
      // edges are given in.
      if (graph != e.graph) {
        return graph < e.graph;
      } else if (node1 != e.node1) {
        return node1 < e.node1;
      } else {
        return node2 < e.node2;
      }
    }

    /* Test whether this edge is the same as @a e in a global order. To be
     * equal, edges must join the same nodes in the same graph.
     */
    bool operator==(const InternalEdge_& e) const {
      return graph == e.graph && node1 == e.node1 && node2 == e.node2;
    }
  };

  std::vector<InternalNode_> nodes_;  // Collection of nodes in graph.
  // NOTE: Edges are stored in a vector, rather than, e.g., a map with two node
  // indices as keys, to facilitate fast lookup by uid_ in Edge::fetch(), which
  // is required for this implementation of the proxy design pattern.
  std::vector<InternalEdge_> edges_;  // Collection of edges in graph.
};

#endif // CME212_GRAPH_HPP
