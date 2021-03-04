#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

//--functionality_0
//--great job!
//--END

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
 * @tparam V This template shows the type used for node values.
 * @tparam E Type used for edge values.
 */
template <typename V, typename E> class Graph {
public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of node value **/
  using node_value_type = V;

  /** Type of edge value **/
  using edge_value_type = E;

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

  /** Construct an empty graph.
   * A constructor to put number of nodes and edges in this
   * graph to zero also initialize the data for nodes and edges
   */
  Graph()
      : nodes_num(0), unique_node(0), node_i2u(), edges_num(0), unique_edge(0),
        edge_i2u(), all_nodes(), all_locs(), all_values(), all_edges(),
        edge_values() {}

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
     * @pre node lives on the graph and is valid.
     * @return The position of the node.
     */
    const Point &position() const {
      // We make sure that user does not use an invalid node
      // and try to take its position
      assert(this->is_valid());
      // Since order of the all_locs is the same as index
      // we can easily access the location of node with
      // its unique id
      size_type uid = nodei_graph->node_i2u[nodei_id];
      return nodei_graph->all_locs[uid];
    }

    /** Return this node's position as a reference.
     * @pre node is valid.
     * @return The position of the node.
     *
     * @note This is just like the previous method defined.
     */
    Point &position() {
      assert(this->is_valid());
      size_type uid = nodei_graph->node_i2u[nodei_id];
      return nodei_graph->all_locs[uid];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // We make sure that user does not use an invalid node
      // and try to take its index
      assert(this->is_valid());
      // Returning the attribute.
      return nodei_id;
    }

    /** Return the value of this node */
    node_value_type &value() {
      assert(this->is_valid());
      size_type uid = nodei_graph->node_i2u[nodei_id];
      return nodei_graph->all_values[uid];
    }

    /** Return the value of this node as a const object */
    const node_value_type &value() const {
      assert(this->is_valid());
      size_type uid = nodei_graph->node_i2u[nodei_id];
      return nodei_graph->all_values[uid];
    }

    /** Return number of edges connecting to this node */
    size_type degree() const {
      assert(this->is_valid());
      size_type uid = nodei_graph->node_i2u[nodei_id];
      return (size_type)nodei_graph->edge_incident[uid].size();
    };

    /** Starting iterator of edges connected to this node */
    incident_iterator edge_begin() const {
      assert(this->is_valid());
      return IncidentIterator(this, 0);
    }

    /** Ending iterator that goes over edges connected to this node */
    incident_iterator edge_end() const {
      assert(this->is_valid());
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
      if (nodei_id != n.nodei_id)
        return (nodei_id < n.nodei_id);
      else
        // If two nodes got same id we check their pointer to
        // graph
        return nodei_graph < n.nodei_graph;
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
    /** A method to check if this node is a valid one
     * Here we make sure that id of the node is in the
     * valid range and it had not been removed before.
     *
     * @note By id we mean both user id and unique id
     * of the node and check if their relation holds
     */
    bool is_valid() const {
      // Unique id of thid node
      size_type uid = nodei_graph->node_i2u[nodei_id];
      // Make sure node belong to a graph
      if (nodei_graph == nullptr or
          nodei_graph->all_nodes[uid].nodei_graph == nullptr)
        return 0;
      // Check the size and relations
      return nodei_id >= 0 and nodei_id <= nodei_graph->node_i2u.size() and
             uid >= 0 and uid <= nodei_graph->all_nodes.size() and
             nodei_graph->node_i2u[nodei_graph->all_nodes[uid].nodei_id] == uid;
    }
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
    // Add the unique id of this node to vector of active nodes
    node_i2u.push_back(unique_node);
    // Now after this one we have one more node in the
    // graph
    ++nodes_num;
    ++unique_node;
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    return n.nodei_graph == this and n.is_valid();
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
    // the node already in order, so we simply map
    // the user id to the unique one and find the
    // value
    return all_nodes[node_i2u[i]];
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
      assert(this->is_valid());
      // Return the first node with the help of edge pointer
      // and index of its node.
      // @note This node1_index is the unique node id
      // so we do not need to map that
      return edgei_graph->all_nodes[node1_index];
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Make sure that the edge is valid
      assert(this->is_valid());
      // See above for explanation
      return edgei_graph->all_nodes[node2_index];
    }

    /** Return the value of this node */
    edge_value_type &value() {
      // First we make sure that the edge is valid
      assert(this->is_valid());
      size_type uid = edgei_graph->edge_i2u[edgei_id];
      return edgei_graph->edge_values[uid];
    }

    /** Return the value of this node */
    const edge_value_type &value() const {
      // First we make sure that the edge is valid
      assert(this->is_valid());
      size_type uid = edgei_graph->edge_i2u[edgei_id];
      return edgei_graph->edge_values[uid];
    }

    /** A method to return the length of this edge based on
     * the euclidean distance of its ending nodes
     */
    double length() const {
      assert(this->is_valid());
      auto delta_x = this->node1().position() - this->node2().position();
      return norm(delta_x);
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
      if (edgei_id != e.edgei_id)
        return (edgei_id < e.edgei_id);
      else
        // If they got same id we check their pointers to graph
        return edgei_graph < e.edgei_graph;
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
    /** A method to check if this edge is a valid one
     * Here we make sure that id of the edge is in the
     * valid range and it had not been removed before.
     *
     * @note By id we mean both user id and unique id
     * of the node and check if their relation holds
     * (see node's is_valid for more info)
     */
    bool is_valid() const {
      // Unique id of thid edge
      size_type uid = edgei_graph->edge_i2u[edgei_id];
      if (edgei_graph == nullptr or
          edgei_graph->all_edges[uid].edgei_graph == nullptr)
        return 0;
      return edgei_id >= 0 and edgei_id <= edgei_graph->edge_i2u.size() and
             uid >= 0 and uid <= edgei_graph->all_edges.size() and
             edgei_graph->edge_i2u[edgei_graph->all_edges[uid].edgei_id] == uid;
    }
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
  size_type num_edges() const { return edges_num; }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i >= 0);
    assert(i < edges_num);
    return all_edges[edge_i2u[i]];
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
    if (node_i2u[a.index()] < node_i2u[b.index()])
      nodes = {node_i2u[a.index()], node_i2u[b.index()]};
    else
      nodes = {node_i2u[b.index()], node_i2u[a.index()]};
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
      if (node_i2u[a.index()] < node_i2u[b.index()])
        nodes = {node_i2u[a.index()], node_i2u[b.index()]};
      else
        nodes = {node_i2u[b.index()], node_i2u[a.index()]};
      // Add this nodes and edge to our mapping
      // (See definition of node_edge for more information)
      node_edge[nodes] = unique_edge;
      Edge e(edges_num, node_i2u[a.index()], node_i2u[b.index()], this);
      // Add this edge for edge_incident mapping
      // (See definition of edge_incident for more information)
      edge_incident[node_i2u[a.index()]].push_back(unique_edge);
      edge_incident[node_i2u[b.index()]].push_back(unique_edge);
      // Add the edge to vector of all edges
      all_edges.push_back(e);
      // Ådd the appropriate value of this edge
      // For now we just add the default value
      edge_values.push_back(edge_value_type());
      // Add unique id of this edge to vector of all active edges
      edge_i2u.push_back(unique_edge);
      // Add one to all the edge number
      ++edges_num;
      ++unique_edge;
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
    // Find the user id of that edge
    auto extra_edge_id = all_edges[extra_edge].edgei_id;
    return Edge(extra_edge_id, node_i2u[a.index()], node_i2u[b.index()], this);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_num = 0;
    unique_node = 0;
    node_i2u.clear();
    edges_num = 0;
    unique_edge = 0;
    edge_i2u.clear();
    all_nodes.clear();
    all_locs.clear();
    all_edges.clear();
    edge_values.clear();
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
    Node operator*() const {
      // First we find the unique node id
      size_type uid = p->node_i2u[id];
      return p->all_nodes[uid];
    };
    NodeIterator &operator++() {
      ++id;
      return *this;
    }
    bool operator==(const NodeIterator &iter) const {
      return (p == iter.p and id == iter.id);
    }

  private:
    friend class Graph;
    // A pointer to graph in which this iterator
    // belongs to
    const graph_type *p = nullptr;
    // A number (id) showing the position of this ptr in the vector
    // @note This id would just go over active nodes
    size_type id = 0;
    /* A private constructor that would be called by Graph class
     *
     * @param[in] ptr A pointer to the Graph
     * @param given_id The id to be assigned to the iterator.
     */
    NodeIterator(const graph_type *ptr, size_type given_id)
        : p{ptr}, id{given_id} {}
  };

  // node_iterator methods to be used from Graph class for iterating over
  // all the nodes.
  node_iterator node_begin() const { return NodeIterator(this, 0); }
  node_iterator node_end() const { return NodeIterator(this, node_i2u.size()); }

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
      // Unique id of this node
      size_type uid = graph_ptr->node_i2u[id_node];
      // This is the unique id (number) of the edge connected to this node
      // that we are currently on
      size_type id_edge = graph_ptr->edge_incident[uid][id];
      // Read the data from vector of all edges saved in Graph
      // class
      Edge e = graph_ptr->all_edges[id_edge];
      size_type first_id = e.node1().index();
      // We know that the returned edge, should have its first edge
      // as this node that we are currently on, so we check and see
      // if that is not the case we make a new edge which holds
      // that requirement.
      if (first_id != id_node) {
        return Edge(e.edgei_id, uid, graph_ptr->node_i2u[first_id], graph_ptr);
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
    // A pointer to the node, we use
    // this to see what edges are connected to each node
    const node_type *p = nullptr;
    // An id (number) showing the position of iterator in the
    // vector of all (active) edges connected to a node
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
    Edge operator*() const {
      size_type uid = p->edge_i2u[id];
      return p->all_edges[uid];
    }
    EdgeIterator &operator++() {
      ++id;
      return *this;
    }
    bool operator==(const EdgeIterator &iter) const {
      return (p == iter.p and id == iter.id);
    }

  private:
    friend class Graph;
    // A pointer to graph in which this iterator
    // belongs to
    const graph_type *p = nullptr;
    // A number (id) showing the position of this ptr in the vector
    // @note This id goes over active edges
    size_type id = 0;
    /* A private constructor that would be called by Graph class
     *
     * @param[in] ptr A pointer to Graph.
     * @param given_id The id assigned to the iterator.
     */
    EdgeIterator(const graph_type *ptr, size_type given_id)
        : p{ptr}, id{given_id} {}
  };

  // edge_iterator methods to be used from Graph class for iterating over
  // all the edges.
  edge_iterator edge_begin() const { return EdgeIterator(this, 0); }
  edge_iterator edge_end() const { return EdgeIterator(this, edge_i2u.size()); }

  // #############
  // Removal methods
  // ############

  /** A method to remove a node in the graph
   *
   * @param[in] n Given node by the user to remove
   * @returns 1 if method was able to remove such node
   *          0 otherwise
   *
   * @pre Given node is a valid node
   * @post This method invalidate all previously instantiated nodes
   * and edges (and all their copies, so user should no use them)
   * @post This method invalidate node and edge and incident iterators
   * instantiated before the call to remove_node()
   *
   * @note Time complexity: O(removing connected edges to the node)
   * or
   * O(_n_log(_d_)) where _n_ is number of connected edges to this node (way
   * smaller than total number of edges _n_ << _d_) and _d_ is total number of
   * edges in the graph so we can say O(log(_d_))
   */
  size_type remove_node(const Node &n) {
    // First make sure this node belongs to the graph
    // @note n.is_valid() is checked inside of has_node()
    if (this->has_node(n)) {
      // This is the id that users care about
      // @note there is no gap in this id and it starts from 0
      // all the way to one less than the current number of edges in
      // the graph.
      size_type idx = n.index();
      // The unique id of the node
      size_type uid = node_i2u[idx];

      // If this node was removed before return 0
      // (just to be safe)
      if (all_nodes[uid].nodei_graph == nullptr)
        return 0;

      // We go over all edges connected to this node and remove them
      // Time complexity: O(num of connected edges to this node)
      std::vector<size_type> connected_edges = edge_incident[uid];
      for (auto &e_remove : connected_edges)
        remove_edge(all_edges[e_remove]);

      // Remove this node from list of active nodes
      // Time complexity: O(1)
      size_type last_elem_uid = node_i2u.back();
      std::swap(node_i2u[idx], node_i2u.back());
      node_i2u.pop_back();

      // Now we set the id of nodes to what uesrs care about
      // Time complexity: O(1)
      all_nodes[last_elem_uid].nodei_id = idx;

      // Invalide this edge
      all_nodes[uid].nodei_graph = nullptr;

      --nodes_num;
      return 1;
    }
    return 0;
  }

  /** A method to remove a node
   *
   * @param n_it node iterator to remove
   * @returns A valid node_iterator
   *
   * @pre n_it is a valid iterators
   *
   * @note See above for more info and post conditions
   */
  node_iterator remove_node(node_iterator n_it) {
    this->remove_node(*n_it);
    return n_it;
  }

  /** A method to remove an edge in the graph
   *
   * @param[in] node1 Given node by the user
   * @param[in] node2 Given node by the user
   * @returns 1 if method was able to remove an edge
   * 					connecting _node1_ and _node2_
   *         	0 otherwise
   *
   * @pre Given nodes are valid
   * @post This method invalidate all previously instantiated edges
   * (and all their copies, so user should no use them)
   * @post This method invalidate edge and incident iterators.
   *
   * @note Time complexity: O(log(_d_))
   * where _d_ is number of all edges and _l_ number of
   * edges connected to _node1_ and _node2_ which we assume is l << d so Time
   * complexity: O(log(_d_) + _l_) = O(log(_d_))
   *
   */
  size_type remove_edge(const Node &node1, const Node &node2) {
    // Make sure an edge exists in the first place
    if (this->has_edge(node1, node2)) {
      // Note that we need to change id of the nodes
      // back to unique ids
      size_type a = node_i2u[node1.index()];
      size_type b = node_i2u[node2.index()];

      std::vector<size_type> e_nodes;
      if (a < b)
        e_nodes = {a, b};
      else
        e_nodes = {b, a};

      // Find the unique id of this edge
      size_type uid = node_edge[e_nodes];

      // If this edge was removed before return 0
      // (just to be safe)
      if (all_edges[uid].edgei_graph == nullptr)
        return 0;

      // This is the id that users care about
      // @note there is no gap in this id and it starts from 0
      // all the way to one less than the current number of edges in
      // the graph.
      size_type idx = all_edges[uid].edgei_id;

      // Now we delete this edge from list of
      // our connections, meaning that it's node
      // are not connected anymore
      // Time complexity: Logarithmic in the size of the container
      auto ite = node_edge.find(e_nodes);
      // Time complexity: Amortized constant
      node_edge.erase(ite);

      // Now we remove the edge from list of its node connections
      // Time complexity: At most last - first of the containter
      auto e1 =
          std::find(edge_incident[a].begin(), edge_incident[a].end(), uid);
      // Time complexity: O(1)
      std::swap(*e1, edge_incident[a].back());
      edge_incident[a].pop_back();

      auto e2 =
          std::find(edge_incident[b].begin(), edge_incident[b].end(), uid);
      // Time complexity: O(1)
      std::swap(*e2, edge_incident[b].back());
      edge_incident[b].pop_back();

      // Remove this edge from list of active edges
      // Time complexity: O(1)
      size_type last_elem_uid = edge_i2u.back();
      std::swap(edge_i2u[idx], edge_i2u.back());
      edge_i2u.pop_back();

      // Now we set the id to edges to what uesrs care about
      // Time complexity: O(1)
      all_edges[last_elem_uid].edgei_id = idx;

      // Invalide this edge
      all_edges[uid].edgei_graph = nullptr;

      --edges_num;
      return 1;
    }
    return 0;
  }

  /** A method to remove an edge
   *
   * @param[in] e Given edge by the user.
   * @returns 1 If method was able to remove the edge.
   * 					0 otherwise
   *
   * @pre _e_ is a valid edge
   *
   * @note See above for more info and post conditions
   */
  size_type remove_edge(const Edge &e) {
    size_type uid = edge_i2u[e.edgei_id];
    // First we make sure that this edge is a valid one
    if (all_edges[uid].edgei_graph == nullptr or e.edgei_graph == nullptr)
      return 0;
    return remove_edge(e.node1(), e.node2());
  }

  /** A method to remove an edge
   *
   * @param e_it edge iterator to be removed
   * @returns A valid edge_iterator
   *
   * @pre e_it is a valid iterator
   *
   * @note See above for more info and post conditions
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    this->remove_edge(*e_it);
    return e_it;
  }

private:
  // ######## data attributes ######## //
  // A number showing the number of all nodes in this graph
  size_type nodes_num;
  // A unique number to id all the nodes added to the graph (even if
  // later we remove them)
  size_type unique_node;
  // A vector of all the active nodes (this vector saves unique id of
  // each node)
  std::vector<size_type> node_i2u;
  // A number showing the number of all edges in this graph
  size_type edges_num;
  // A unique number to id all the edges ever added to the graph
  size_type unique_edge;
  // A vector of all active edges
  std::vector<size_type> edge_i2u;

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
  // A vector to hold all the values of edges in the same order
  // as all_edges vector
  std::vector<edge_value_type> edge_values;

  // This is a mapping that shows what two nodes are connected to what
  // edge of this graph. Note that here we use just the indices of nodes
  // and edges as that's all we need. E.g. node_edge[{a,b}] = i, means
  // that nodes with index a and b are connected with an edge with index i
  // Also to make the finding in the map easier when we populate the map
  // we make sure that a < b (so the key is always in order)
  // @note We use unique node and edge ids in this mapping
  // this is important, if later we want to read data from this
  // container we need to first map any user id to unique one
  // and then use this map
  std::map<std::vector<size_type>, size_type> node_edge;
  // A mapping showing what edges are connected to each node.
  // For example edge_incident[_a_] = {_i_, _j_, _k_} means that edges
  // _i_, _j_, and _k_ are connected to node _a_.
  // @note We use unique node and edge ids in this mapping
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
    if (node_i2u[a.index()] < node_i2u[b.index()])
      nodes = {node_i2u[a.index()], node_i2u[b.index()]};
    else
      nodes = {node_i2u[b.index()], node_i2u[a.index()]};
    return node_edge[nodes];
  }
};

#endif // CME212_GRAPH_HPP
