#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include <iostream>


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
  int length;
  std::vector<Point> points;
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  using node_value_type = V;
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
  Graph(): length(0), points() {
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
  class Node: private totally_ordered<Node>  {
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
    node_point = nullptr;
    }

    /** Return this node's value*/
    node_value_type& value() {
      return node_point->opt_values.at(index_loc);
    }

    const node_value_type& value() const {
      return node_point->opt_values.at(index_loc);
    }


    /** Return this node's position. */
    const Point& position() const {
      return find();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_loc;
    }


    /** Return this node's degree, a number in the range [0, edge_size)*/
    size_type degree() {
      node_incid = node_point->incidents.at(index_loc);
      return node_incid.size();
    }

    /** Return an Incident Iterator starting at the beginning of this node's incident edges*/
    incident_iterator edge_begin() {
      node_incid = node_point->incidents.at(index_loc);
      return IncidentIterator(const_cast<Edge*>(&node_incid[0]));
    }

    /** Return an Incident Iterator starting at the end of this node's incident edges*/
    incident_iterator edge_end() {
      node_incid = node_point->incidents.at(index_loc);
      return IncidentIterator(const_cast<Edge*>(&node_incid[node_incid.size()]));
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
      return find().data() == n.find().data();
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
      if(index_loc < n.index()){
        return true;
      } else {
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    Graph* node_point;
    size_type index_loc;

    /** Node constructor for the Graph class, sets a Graph pointer and index value
     * @param[in] Graph object pointer
     * @param[in] index of the Node
     */
    Node(const Graph* node_point, size_type index): node_point(const_cast<Graph*>(node_point)), index_loc(index) {
    }

    /** Function to return the node's associated Point object*/
    Point &find() const {
      return node_point->points[index_loc];
      assert(false);
    }

    std::vector<edge_type> node_incid;
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return (size_type)length;
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
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    points.push_back(position);
    ++length;
    std::vector<edge_type> incident_edges;
    incidents.push_back(incident_edges);
    opt_values.push_back(node_value_type());
    Node new_node = Node(this, length - 1);
    nodes.push_back(new_node);
    return new_node;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if(n.index() < nodes.size()){
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return nodes.at(i);        // Invalid node
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
  class Edge: private totally_ordered<Edge>  {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    graph = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph->node(a_point);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph->node(b_point);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(e.graph == graph && e.a_point == a_point && e.b_point == b_point)
        return true;
      if(e.graph == graph && e.a_point == b_point && e.b_point == a_point)
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(get_index() < e.get_index()){
        return true;
      } else {
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    Graph* graph;
    size_type a_point;
    size_type b_point;

    /** Private constructor for the Graph class, sets a Graph pointer and endpoint nodes
     * @param[in] gr Graph pointer object
     * @param[in] node1 index of node 1
     * @param[in] node2 index of node 2
     *
     * Edge is undirected
     */
    Edge(const Graph* gr, size_type node1, size_type node2):
        graph(const_cast<Graph*>(gr)), a_point(node1), b_point(node2) {
    }

    /** function to find the index of the edge*/
    size_type get_index() {
      std::pair<size_type, size_type> endpoints(a_point, b_point);
      if(b_point < a_point){
        std::pair<size_type, size_type> new_pair(b_point, a_point);
        endpoints = new_pair;
      }
      auto it = std::find(node_pairs.begin(), node_pairs.end(), endpoints);
      return (size_type)(it - node_pairs.begin());
    }

    friend class Graph;
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
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edges.at(i);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    std::pair<size_type, size_type> edge_nodes(a.index(), b.index());

    if(a.index() > b.index()){
      std::pair<size_type, size_type> new_pair(b.index(), a.index());
      edge_nodes = new_pair;
    }
    auto it = std::find(node_pairs.begin(), node_pairs.end(), edge_nodes);
    return it != node_pairs.end();
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
    if(has_edge(a,b))
      return Edge(this, a.index(),b.index());

    Edge new_edge = Edge(this, a.index(), b.index());

    std::pair<size_type, size_type> nodep;

    if(a.index() > b.index()) {
      std::pair<size_type, size_type> new_p(b.index(), a.index());
      nodep = new_p;
    }
    node_pairs.push_back(nodep);
    edges.push_back(new_edge);

    incidents.at(a.index()).push_back(new_edge);
    incidents.at(b.index()).push_back(Edge(this, b.index(), a.index()));
    return new_edge;        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    points.clear();
    nodes.clear();
    edges.clear();
    node_pairs.clear();
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
    NodeIterator(pointer ptr): m_ptr(ptr) {
    }

    /** Dereference operator, returns Node object*/
    Node operator*() const{
      return *m_ptr;
    }

    /** Addition operator, moves the iterator forward by 1*/
    NodeIterator& operator++() {
      m_ptr++;
      return *this;
    }

    /** Equal comparison operator, checks if two iterators are equal, returns bool
     * @param[in] x, NodeIterator object
     */
    bool operator==(const NodeIterator& x) const {
      return m_ptr == x.m_ptr;
    }

    /** Unequal comparison operator, checks if two iterators are not equal, returns bool
     * @param[in] x, NodeIterator object
     */
    bool operator!=(const NodeIterator& x) const {
      return m_ptr != x.m_ptr;
    }
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    pointer m_ptr;
    // HW1 #2: YOUR CODE HERE
  };

  /** Creates a Node Iterator located at the start of the nodes vector*/
  node_iterator node_begin() const {
    return NodeIterator(const_cast<Node*>(&nodes[0]));
  }

  /** Creates a Node Iterator located at the end of the nodes vector*/
  node_iterator node_end() const {
    return NodeIterator(const_cast<Node*>(&nodes[nodes.size()]));
  }
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
    IncidentIterator(pointer pt): i_ptr(pt) {
    }

    /** Dereference operator, returns Edge object */
    Edge operator*() const {
      return *i_ptr;
    }

    /** Addition operator, moves the iterator forward by 1 */
    IncidentIterator& operator++() {
      i_ptr++;
      return *this;
    }

    /** Equal comparison operator, checks if two iterators are equal, returns bool
     * @param[in] y, IncidentIterator object
     */
    bool operator==(const IncidentIterator& y) const {
      return i_ptr == y.i_ptr;
    }

    /** Unequal comparison operator, checks if two iterators are not equal, returns bool
     * @param[in] y, IncidentIterator object
     */
    bool operator!=(const IncidentIterator& y) const {
      return i_ptr != y.i_ptr;
    }


    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;
    pointer i_ptr;
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
    EdgeIterator(pointer p): e_ptr(p) {
    }

    /** Dereference operator, returns Edge object*/
    Edge operator*() const {
      return *e_ptr;
    }

    /** Addition operator, moves the iterator forward by 1*/
    EdgeIterator& operator++() {
      e_ptr++;
      return *this;
    }

    /** Equal comparison operator, checks if two iterators are equal, returns bool
     * @param[in] z, EdgeIterator object
     */
    bool operator==(const EdgeIterator& z) const {
      return e_ptr == z.e_ptr;
    }

    /** Unequal comparison operator, checks if two iterators are not equal, returns bool
     * @param[in] z, EdgeIterator object
     */
    bool operator!=(const EdgeIterator& z) const {
      return e_ptr != z.e_ptr;
    }
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    pointer e_ptr;
    // HW1 #5: YOUR CODE HERE
  };

  /**Creates an EdgeIterator object located at the start of the edges vector */
  edge_iterator edge_begin() const {
    return EdgeIterator(const_cast<Edge*>(&edges[0]));
  }

  /**Creates an EdgeIterator object located at the end of the edges vector */
  edge_iterator edge_end() const {
    return EdgeIterator(const_cast<Edge*>(&edges[edges.size()]));
  }
  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:
  size_type edge_count;

  std::vector<node_type> nodes;
  std::vector<edge_type> edges;
  std::vector<std::pair<size_type, size_type>> node_pairs;

  std::vector<std::vector<edge_type>> incidents;

  std::vector<node_value_type> opt_values;
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
