#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
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
  using node_value_type = V;

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
    // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph ? graph->points[index] : Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Returns the type that is stored in the node.*/
    node_value_type & value () {
      return node_value;
    }

    /** Returns the type that is stored in the node.*/
    const node_value_type & value () const {
      return node_value;
    }

    /** Returns the number of neighboring nodes. */
    size_type degree() const {
      return graph->node_to_nodes[index].count();
    }

    /** Returns the beginning of the iterator that iterates through all the neighboring edges. */
    incident_iterator edge_begin() const {
      IncidentIterator iter = IncidentIterator();
      iter.graph = graph;
      iter.node_index = index;
      iter.index = 0;
      return iter;
    }

    /** Returns the beginning of the iterator that iterates through all the neighboring edges. */
    incident_iterator edge_end() const {
      IncidentIterator iter = IncidentIterator();
      iter.graph = graph;
      iter.node_index = index;
      iter.index = degree();
      return iter;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return graph == n.graph && index == n.index;
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
      // HW0: YOUR CODE HERE
      return (index < n.index);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph = nullptr;
    size_type index = -1;
    node_value_type node_value;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_nodes;
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
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    Node node = Node();

    node.node_value = node_value;
    node.index = num_nodes;
    node.graph = this;

    num_nodes += 1;
    nodes.emplace_back(node);
    points.emplace_back(position);

    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return nodes[i];
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph ? graph->nodes[index1] : Node();
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph ? graph->nodes[index2] : Node();
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return graph == e.graph && index1 == e.index1 && index2 == e.index2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (index1 == e.index1)
        return index2 < e.index2;
      return index1 < e.index1;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph = nullptr;
    size_type index1 = -1;
    size_type index2 = -1;
    size_type index = -1;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edges[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return node_to_nodes[a.index].contains(b.index);
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
    Edge e = Edge();

    if (has_edge(a, b)) {
      for (size_type i : node_to_edges[a.index]) {
        e = edges[i];
        if (e.node1() == b || e.node2() == b)
          return e;
      }
    }

    if (a.index < b.index) {
      e.index1 = a.index;
      e.index2 = b.index;
    } else {
      e.index1 = b.index;
      e.index2 = a.index;
    }

    e.index = num_edges;
    num_edges += 1;
    edges.emplace_back(e);

    // Update the maps.
    if (node_to_nodes.count(a.index))
      node_to_nodes.insert(a.index, std::set<size_type>());
    node_to_nodes[a.index].insert(b.index);

    if (node_to_nodes.count(b.index))
      node_to_nodes.insert(b.index, std::set<size_type>());
    node_to_nodes[b.index].insert(a.index);

    if (node_to_edges.count(a.index))
      node_to_edges.insert(a.index, std::vector<size_type>());
    node_to_edges[a.index].emplace_back(e.index);

    if (node_to_edges.count(b.index))
      node_to_edges.insert(b.index, std::vector<size_type>());
    node_to_edges[b.index].emplace_back(e.index);

    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear();
    node_to_nodes.clear();
    node_to_edges.clear();
    num_nodes = 0;
    num_edges = 0;
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
    Node operator*() const {
      return &graph->node(index);
    }

    NodeIterator& operator++() {
      index += 1;
    }
    bool operator==(const NodeIterator& n) const {
      return graph == n.graph && index == n.index;
    }

   private:
    friend class Graph;
    Graph* graph = nullptr;
    size_type index = -1;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Returns the beginning of the node iterator of this graph. */
  node_iterator node_begin() const {
    NodeIterator iter = NodeIterator();
    iter.index = 0;
    iter.graph = this;
    return iter;
  }

  /** Returns the beginning of the node iterator of this graph. */
  node_iterator node_end() const {
    NodeIterator iter = NodeIterator();
    iter.index = num_nodes;
    iter.graph = this;
    return iter;
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
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Returns the current edge of the iterator. */
    Edge operator*() const {
      *curr_edge = graph->edges[graph->node_to_edges[node_index][index]];
      if (curr_edge.index1 == node_index)
        return curr_edge;
      size_type temp_index = curr_edge.index1;
      curr_edge.index1 = curr_edge.index2;
      curr_edge.index2 = temp_index;
      return curr_edge;
    }

    /** Updates the iterator to move to the next edge. */
    IncidentIterator& operator++() {
      index += 1;
      return *this;
    }

    /** Returns true only if the other incident iterator is at the same edge. */
    bool operator==(const IncidentIterator& i) const {
      return graph == i.graph && node_index = i.node_index && index == i.index;
    }

   private:
    friend class Graph;
    Graph* graph = nullptr;
    Edge* curr_edge = Edge();
    size_type node_index = -1;
    size_type index = -1;
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

    /** Returns the current edge of the iterator. */
    Edge operator*() const {
      return &graph->edge(index);
    }

    /** Updates the iterator to the next edge. */
    EdgeIterator& operator++() {
      index += 1;
    }

    /** Returns true only if |e| is at the same edge of the same graph. */
    bool operator==(const EdgeIterator& e) const {
      return graph == e.graph && index == e.index;
    }

   private:
    friend class Graph;
    Graph* graph = nullptr;
    size_type index = -1;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    EdgeIterator iter = EdgeIterator();
    iter.index = 0;
    iter.graph = this;
    return iter;
  }
  edge_iterator edge_end() const {
    EdgeIterator iter = EdgeIterator();
    iter.index = num_edges;
    iter.graph = this;
    return iter;
  }

 private:
  std::vector<Node> nodes = {};
  std::vector<Point> points = {};
  std::vector<Edge> edges = {};
  std::map<size_type, std::set<size_type>> node_to_nodes;
  std::map<size_type, std::vector<size_type>> node_to_edges≈;
  size_type num_nodes = 0;
  size_type num_edges = 0;
};

#endif // CME212_GRAPH_HPP
