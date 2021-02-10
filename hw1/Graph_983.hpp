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
 *
 * @paramt V This template shows the type used for node values.
 */
template <typename V> class Graph {
public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of node value **/
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
  // A constructor to put number of nodes and edges in this
  // graph to zero also initialize the data for nodes and edges
  Graph()
      : nodes_num(0), edges_num(0), all_nodes(), all_locs(), all_values(),
        all_edges() {}

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   *
   * This class privately inherits from totally_ordered struct, which
   * with the help of comparisons defined here, would define all
   * possible comparisons between two Node objects.
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
      // We here put the node id to value of 0.
      // Although in reality this may not happen but since
      // we know that we would never use a Node object that is
      // defined invalid (for direct calculations)
      // in our graph we can do this.
      nodei_id = 0;
      // We also set the attribute pointing to the
      // graph as nullptr for now
      nodei_graph = nullptr;
    }

    /** Return this node's position.
     * @pre node lives on the graph.
     * @return The position of the node.
     */
    const Point &position() const {
      // We make sure that user does not use an invalid node
      // and try to take its position
      assert(nodei_graph != nullptr);
      // Since order of the all_locs is the same as index
      // we can easily access the location of node with index id
      return nodei_graph->all_locs[nodei_id];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // (For futrue HWs this might be a reference)
      // We make sure that user does not use an invalid node
      // and try to take its index
      assert(nodei_graph != nullptr);
      // Returning the attribute.
      return nodei_id;
    }

    // Return the value of this node
    node_value_type &value() {
      assert(nodei_graph != nullptr);
      // Since order of the all_values is the same as index
      // we can easily access the location of node with index id
      return nodei_graph->all_values[nodei_id];
    }

    // Return the value of this node as a const object
    const node_value_type &value() const {
      assert(nodei_graph != nullptr);
      // Since order of the all_values is the same as index
      // we can easily access the location of node with index id
      return nodei_graph->all_values[nodei_id];
    }

    // Return number of edges connecting to this node
    size_type degree() const {
      return (size_type)nodei_graph->edge_incident[nodei_id].size();
    };

    // Starting iterator of edges connected to this node
    incident_iterator edge_begin() const { return IncidentIterator(this, 0); }
    // Ending iterator that goes over edges connected to this node
    incident_iterator edge_end() const {
      return IncidentIterator(this, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {
      return (nodei_graph == n.nodei_graph and nodei_id == n.nodei_id);
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node &n) const {
      // First we make sure both nodes live on
      // the same graph
      assert(nodei_graph == n.nodei_graph);
      return (nodei_id < n.nodei_id);
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // A number showing the index of this node
    size_type nodei_id;
    // A pointer showing the Graph to which this particular node
    // belongs, we can use this to see if two nodes are the same
    // (as for this, one of the requirements is that they both belong
    // to the same graph)
    graph_type *nodei_graph;
    /** Ordinary private constructor of Node class used inside the Graph class
     * while making valid nodes
     *
     * @param[in] id The id of this particular node (we assign this number)
     * @param graph A pointer to an object of class Graph, showing
     * the graph to which this node belongs.
     */
    Node(const size_type id, graph_type *graph)
        : nodei_id(id), nodei_graph(graph) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const { return nodes_num; }

  /** Synonym for size(). */
  size_type num_nodes() const { return size(); }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position,
                const node_value_type &nvalue = node_value_type()) {
    // A node that would be made by this process
    // Since each added node is the last in the list
    // we can use total number of nodes in the graph
    // as index
    // Setting the graph to which this node belongs
    // which is the graph that we are currently on
    Node n{nodes_num, this};
    // Add this node to the vector
    // of all nodes of this graph
    all_nodes.push_back(n);
    // Add the location of the node to vector of all
    // positions
    all_locs.push_back(position);
    // Add the value of the node to vector of all
    // values
    all_values.push_back(nvalue);
    // Now after this one we have one more node in the
    // graph
    ++nodes_num;
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    // HW0: YOUR CODE HERE
    return (n.nodei_graph == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // First check the bounds
    assert(i >= 0);
    assert(i < nodes_num);
    // We know that in the all_nodes we have entered
    // the node already in order, so node with index i
    // is already at position i of the all_nodes list
    return all_nodes[i];
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   *
   * This class privately inherits from totally_ordered struct, which
   * with the help of comparisons defined here, would define all
   * possible comparisons between two Edge objects.
   */
  class Edge : private totally_ordered<Edge> {
  public:
    /** Construct an invalid Edge.
     * See below for definitions of the attributes
     */
    Edge()
        : edgei_id(0), node1_index(0), node2_index(0), edgei_graph(nullptr) {}

    /** Return a node of this Edge */
    Node node1() const {
      // Make sure that the edge is valid
      assert(edgei_graph != nullptr);
      // Return the first node with the help of edge pointer
      // and index of its node.
      return edgei_graph->all_nodes[node1_index];
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Make sure that the edge is valid
      assert(edgei_graph != nullptr);
      // See above for explanation
      return edgei_graph->all_nodes[node2_index];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      // First we make sure that both edges belong to the same
      // graph
      if (edgei_graph != e.edgei_graph)
        return false;
      // Here we make sure to check the possibility where
      // an edge got nodes [a,b] while the other one
      // has the nodes [b, a] we know that these two are
      // equal.
      // Here we just need to check the indices of nodes (as we know that in
      // first place when we add an edge between two nodes, we make sure that
      // those nodes belong to the graph)
      return ((node1_index == e.node1_index and node2_index == e.node2_index) or
              (node1_index == e.node2_index and node2_index == e.node1_index));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {
      // Make sure that both edges belong to the
      // same graph
      assert(edgei_graph == e.edgei_graph);
      return (edgei_id < e.edgei_id);
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // A number showing the index of this edge
    size_type edgei_id;
    // Two variables showing the index of two nodes of this edge
    size_type node1_index;
    size_type node2_index;
    // A pointer to the graph this edge belongs to
    Graph *edgei_graph;
    /** Ordinary private constructor of Edge class used inside the Graph class
     * to make valid edge objects
     *
     * @param[in] id The id of this particular edge (given by us)
     * @param[in] index1 Index of the first node
     * @param[in] index2 Index of the second node
     * @param graph Pointer to the graph to which this edge belongs
     */
    Edge(const size_type id, const size_type index1, const size_type index2,
         Graph *graph)
        : edgei_id(id), node1_index(index1), node2_index(index2),
          edgei_graph(graph) {}
    // We define IncidentIterator class friend of Edge class, since
    // we know that when we dereference an IncidentIterator of node
    // _a_ we need the first node of the returned edge to be _a_.
    // So to do this in dereferencing an IncidentIterator we use the
    // private edge constructor to make the appropriate edge and
    // return it to the user.
    friend class IncidentIterator;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_num;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i >= 0);
    assert(i < edges_num);
    return all_edges[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {
    // First we check and make sure that node a and b
    // are (distinct) valid nodes
    if (!this->has_node(a) or !this->has_node(b) or a == b)
      return false;
    std::vector<size_type> nodes;
    // First we order the nodes (see the deceleration of node_edge
    // for more information)
    // For short: node_edge is a mapping between nodes and edge
    // connecting those nodes
    if (a < b)
      nodes = {a.index(), b.index()};
    else
      nodes = {b.index(), a.index()};
    // Search and see if these two nodes are connected by an edge
    // in this graph
    if (node_edge.find(nodes) != node_edge.end()) {
      // Found an edge, so they are connected
      return true;
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
  Edge add_edge(const Node &a, const Node &b) {
    // @note has-edge method would make sure that a and b are valid
    // nodes
    // First we make sure that a and b are not the same node
    assert(a != b);
    // First we see if there is an edge between these two nodes
    if (!has_edge(a, b)) {
      // If there is no edge already
      std::vector<size_type> nodes;
      // First we order the nodes for use in our mapping
      if (a < b)
        nodes = {a.index(), b.index()};
      else
        nodes = {b.index(), a.index()};
      // Add this nodes and edge to our mapping
      // (See definition of node_edge for more information)
      node_edge[nodes] = edges_num;
      Edge e(edges_num, a.index(), b.index(), this);
      // Add this edge for edge_incident mapping
      // (See definition of edge_incident for more information)
      edge_incident[a.index()].push_back(edges_num);
      edge_incident[b.index()].push_back(edges_num);
      // Add the edge to vector of all edges
      all_edges.push_back(e);
      // Add one to all the edge number
      ++edges_num;
      return e;
    }
    // This means that we already have an edge that
    // connects these two given nodes.
    // So we simply make an edge with starting node _a_
    // and ending node _b_.
    // @note Here we give the edge id (lest say i here) to be the edge that
    // connected _a_ and _b_. At this point however if the user ask for edge(i)
    // then the first node of that edge might not be _a_ and it might be _b_.
    // @note extra_edge is given as the id, which is the id of the edge already
    // in the graph with same nodes found in get_edgeindex method.
    auto extra_edge = get_edgeindex(a, b);
    return Edge(extra_edge, a.index(), b.index(), this);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_num = 0;
    edges_num = 0;
    all_nodes.clear();
    all_locs.clear();
    all_edges.clear();
    all_values.clear();
    node_edge.clear();
    edge_incident.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator.
   *
   * @note Inherit the != opertator from equality_comparable
   */
  class NodeIterator : private equality_comparable<NodeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : p(nullptr), id(0) {}

    Node operator*() const { return (*p)[id]; };
    NodeIterator &operator++() {
      ++id;
      return *this;
    }
    bool operator==(const NodeIterator &iter) const {
      return (p == iter.p and id == iter.id);
    }

  private:
    friend class Graph;
    // A pointer to the vector of nodes in which this iterator
    // belongs to
    const std::vector<node_type> *p = nullptr;
    // A number (id) showing the position of this ptr in the vector
    size_type id = 0;
    /* A private constructor that would be called by Graph class
     *
     * @param[in] ptr A pointer to the vector of all nodes in the Graph
     * @param given_id The id to be assigned to the iterator.
     */
    NodeIterator(const std::vector<node_type> *ptr, size_type given_id)
        : p{ptr}, id{given_id} {}
  };

  // node_iterator methods to be used from Graph class for iterating over
  // all the nodes.
  node_iterator node_begin() const { return NodeIterator(&all_nodes, 0); }
  node_iterator node_end() const {
    return NodeIterator(&all_nodes, all_nodes.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   *
   * @note Inherit the != opertator from equality_comparable
   */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : p(nullptr), id(0) {}

    Edge operator*() const {
      // The id of the node that this iterator is working on
      size_type id_node = (*p).nodei_id;
      // We know that each node is pointing to a graph in which
      // it belongs, we can find the pointer to that graph
      graph_type *graph_ptr = (*p).nodei_graph;
      // This is the id (number) of the edge connected to this node
      // that we are currently on
      size_type id_edge = graph_ptr->edge_incident[id_node][id];
      // Read the data from vector of all edges saved in Graph
      // class
      Edge e = graph_ptr->all_edges[id_edge];
      size_type first_id = e.node1().index();
      // We know that the returned edge, should have its first edge
      // as this node that we are currently on, so we check and see
      // if that is not the case we make a new edge which holds
      // that requirement.
      if (first_id != id_node) {
        return Edge(id_edge, id_node, first_id, graph_ptr);
      }
      return e;
    }
    IncidentIterator &operator++() {
      ++id;
      return *this;
    }
    bool operator==(const IncidentIterator &iter) const {
      return (p == iter.p and id == iter.id);
    }

  private:
    friend class Graph;
    // A pointer to the edge_incident mapping that we made, we use
    // this to see what edges are connected to each node
    const node_type *p = nullptr;
    // An id (number) showing the position of iterator in the
    // vector of all edges connected to a node
    size_type id = 0;
    /* Ordinary private constructor
     *
     * @param[in] ptr A pointer to the node that this iterator belongs to.
     * @param given_id The id assigned to the iterator.
     */
    IncidentIterator(const node_type *ptr, size_type given_id)
        : p{ptr}, id{given_id} {}
    // We friend this class with Node class so we can use the private
    // constructor to make IncidentIterators
    friend class Node;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator.
   *
   * @note Inherit the != opertator from equality_comparable
   */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : p(nullptr), id(0) {}

    Edge operator*() const { return (*p)[id]; }
    EdgeIterator &operator++() {
      ++id;
      return *this;
    }
    bool operator==(const EdgeIterator &iter) const {
      return (p == iter.p and id == iter.id);
    }

  private:
    friend class Graph;
    // A pointer to the vector of edges in which this iterator
    // belongs to
    const std::vector<edge_type> *p = nullptr;
    // A number (id) showing the position of this ptr in the vector
    size_type id = 0;
    /* A private constructor that would be called by Graph class
     *
     * @param[in] ptr A pointer to the vector of edges in the Graph.
     * @param given_id The id assigned to the iterator.
     */
    EdgeIterator(const std::vector<edge_type> *ptr, size_type given_id)
        : p{ptr}, id{given_id} {}
  };

  // edge_iterator methods to be used from Graph class for iterating over
  // all the edges.
  edge_iterator edge_begin() const { return EdgeIterator(&all_edges, 0); }
  edge_iterator edge_end() const {
    return EdgeIterator(&all_edges, all_edges.size());
  }

private:
  // ######## data attributes ######## //
  // A number showing the number of all nodes in this graph
  size_type nodes_num;
  // A number showing the number of all edges in this graph
  size_type edges_num;
  // A vector of all the nodes in
  // this particular graph
  std::vector<node_type> all_nodes;
  // A vector of location of all nodes in
  // the same order as all_nodes list shown above
  std::vector<Point> all_locs;
  // A vector of values of all the nodes in the
  // same order as all_nodes list
  std::vector<node_value_type> all_values;
  // A vector of all the edges
  std::vector<edge_type> all_edges;
  // This is a mapping that shows what two nodes are connected to what
  // edge of this graph. Note that here we use just the indices of nodes
  // and edges as that's all we need. E.g. node_edge[{a,b}] = i, means
  // that nodes with index a and b are connected with an edge with index i
  // Also to make the finding in the map easier when we populate the map
  // we make sure that a < b (so the key is always in order)
  std::map<std::vector<size_type>, size_type> node_edge;
  // A mapping showing what edges are connected to each node.
  // For example edge_incident[_a_] = {_i_, _j_, _k_} means that edges
  // _i_, _j_, and _k_ are connected to node _a_.
  std::map<size_type, std::vector<size_type>> edge_incident;

  // ######## methods ######## //
  /* A funtion to give back the id of the duplicate
   * edge.
   *
   * @param[in] a One the two nodes of the edge
   * @param[in] b The other node
   * @returns The index of the edge connecting these two nodes.
   * @note: Based on the place that we run this method, we already
   * know that _a_ and _b_ are valid nodes, and they are definitly
   * connected.
   */
  size_type get_edgeindex(const node_type &a, const node_type &b) {
    std::vector<size_type> nodes;
    if (a < b)
      nodes = {a.index(), b.index()};
    else
      nodes = {b.index(), a.index()};
    return node_edge[nodes];
  }
};

#endif // CME212_GRAPH_HPP
