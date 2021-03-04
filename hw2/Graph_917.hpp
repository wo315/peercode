#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
//--functionality_0
//--missing header file <map>
//--START
#include <map>
//--END
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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // A vector holding the locations of the nodes
  // Should be indexed using raw_ind
  std::vector<Point> location;

  // A vector holding the value of the nodes
  // Should be indexed using raw_ind
  std::vector<V> node_values;

  // If a node is the nth one that's ever added to the graph,
  // and it is currently the ith node in the graph,
  // Then ind_to_raw_ind[i] = n
  std::vector<unsigned> ind_to_raw_ind;

  // The reverse mapping of ind_to_raw_ind
  std::vector<unsigned> raw_ind_to_ind;

  // Current number of nodes
  unsigned node_count;

  // Number of nodes that are ever added to this graph
  unsigned ever_added_node_count;

  // Current number of edges
  unsigned edge_count;

  // Number of edges that's ever added to this graph
  unsigned ever_added_edge_count;

  // If (i, j) is an edge in this graph with i < j
  // Then edge_data[i] has key j, all are raw_index
  std::vector<std::set<unsigned>> edge_data;

  // If (i, j) is an edge in this graph with i > j
  // Then edge_data[i] has key j, all are raw_index
  // This is useful to get the degree and traverse the
  // edges adjacent to a node
  std::vector<std::set<unsigned>> edge_data_transpose;

  // If (i, j) is an edge in this graph with i < j
  // Then edge_data[i] has key j and the value of (i, j)
  std::vector<std::map<unsigned, E>> edge_values;

  // @return the location of the node at index @a ind
  const Point& get_location(unsigned ind) const {
    assert(ind < location.size());
    unsigned raw_index = ind_to_raw_ind[ind];
    return location[raw_index];
  }

  const Point& get_location_from_raw(unsigned raw_ind) const {
    return location[raw_ind];
  }

  // @return the raw index of the node at index @a ind
  unsigned get_raw_from_ind(unsigned ind) const {
    assert(ind < ind_to_raw_ind.size());
    return ind_to_raw_ind[ind];
  }

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

  using set_iterator = std::set<size_type>::const_iterator;
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    node_count = 0;
    ever_added_node_count = 0;
    edge_count = 0;
    ever_added_edge_count = 0;
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //
  using node_value_type = V;
  using edge_value_type = E;
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
    Node() : graph_(nullptr) {
      // HW0: YOUR CODE HERE
      raw_index = static_cast<size_type>(-1);
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(valid_node());
      return graph_->get_location_from_raw(raw_index);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->raw_ind_to_ind[raw_index];
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /**
     * @brief Return a reference to the value of this node
     * @pre this is a valid node.
     * @return node_value_type&
     */
    node_value_type& value() {
      return const_cast<Graph*>(graph_)->node_values[raw_index];
    }

    /**
     * @brief Return a const reference to the value of this node
     * @pre this is a valid node.
     * @return node_value_type&
     */
    const node_value_type& value() const { return graph_->node_values[raw_index]; }

    /**
     * @brief Return the number of edges adjacent to this node
     * @pre this is a valid node.
     * @return size_type, the degree of this node.
     */
    size_type degree() const {
      return (graph_->edge_data_transpose[raw_index].size()) +
             (graph_->edge_data[raw_index].size());
    }

    /**
     * @brief Get the incident iterator corresponding to the first edge adjacent to this node
     * @pre This node must be a valid node
     * @return incident_iterator corresponding to the first edge adjacent to this node
     */
    incident_iterator edge_begin() const {
      size_type ind = index();
      return incident_iterator(graph_, ind, graph_->edge_data[raw_index].begin(),
                               graph_->edge_data_transpose[raw_index].begin());
    }

    /**
     * @brief Get the incident iterator corresponding to the end of all edges adjacent to this node
     * @pre This node must be a valid node
     * @return incident_iterator corresponding to the end of all edges adjacent to this node
     */
    incident_iterator edge_end() const {
      size_type ind = index();
      return incident_iterator(graph_, ind, graph_->edge_data[raw_index].end(),
                               graph_->edge_data_transpose[raw_index].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.get_graph_ptr()) && (raw_index == n.get_raw_index());
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
      // Return the order by which the node is added
      // if they are on the same graph
      if(graph_ == n.get_graph_ptr()){
        return raw_index < n.get_raw_index();
      }
      return (graph_ <  n.get_graph_ptr());
    }

    /**
     * @brief Get a reference to the position of this node, so we are able to
     * modify it.
     *
     * @return Point&, a reference to the position of this node
     */
    Point& position() { return const_cast<Graph*>(graph_)->location[raw_index]; }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    const Graph * graph_;

    // raw_index is the number of nodes that's ever added to
    // graph when this node is added, it is unique!
    size_type raw_index;

    Node(const Graph* g, size_type raw) : graph_(g), raw_index(raw) {}

    const Graph* get_graph_ptr() const { return this->graph_; }
    size_type get_raw_index() const { return this->raw_index; }

    //@return True if this Node is a valid node in a graph
    bool valid_node() const {
      return (graph_ != nullptr) &&
             (graph_->raw_ind_to_ind[raw_index] < graph_->num_nodes());
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_count;
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
    // HW0: YOUR CODE HERE
    location.emplace_back(position.x, position.y, position.z);
    edge_data.emplace_back();
    edge_data_transpose.emplace_back();
    edge_values.emplace_back();
    node_values.emplace_back(value);
    ind_to_raw_ind.push_back(ever_added_node_count);
    raw_ind_to_ind.push_back(node_count);
    node_count++;
    ever_added_node_count++;
    return Node(this, ever_added_node_count - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.get_raw_index() >= ever_added_node_count) {
      return false;
    }
    return (this == n.get_graph_ptr()) && (raw_ind_to_ind[n.get_raw_index()] < node_count);
  }
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < node_count);
    return Node(this, ind_to_raw_ind[i]);
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
    Edge() : graph_(nullptr) {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, raw_ind1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, raw_ind2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Quiet compiler warning
      // HW0: YOUR CODE HERE
      if (graph_ != e.get_graph_ptr()) {
        return false;
      }
      size_type raw1, raw2;
      if (raw_ind1 < raw_ind2) {
        raw1 = raw_ind1;
        raw2 = raw_ind2;
      } else {
        raw1 = raw_ind2;
        raw2 = raw_ind1;
      }
      size_type other1, other2;
      if (e.raw_ind1 < e.raw_ind2) {
        other1 = e.raw_ind1;
        other2 = e.raw_ind2;
      } else {
        other1 = e.raw_ind2;
        other2 = e.raw_ind1;
      }
      return (raw1 == other1) && (raw2 == other2);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Quiet compiler warning
      // HW0: YOUR CODE HERE
      if (graph_ == e.get_graph_ptr()) {
        size_type raw1, raw2;
        if (raw_ind1 < raw_ind2) {
          raw1 = raw_ind1;
          raw2 = raw_ind2;
        } else {
          raw1 = raw_ind2;
          raw2 = raw_ind1;
        }
        size_type other1, other2;
        if (e.raw_ind1 < e.raw_ind2) {
          other1 = e.raw_ind1;
          other2 = e.raw_ind2;
        } else {
          other1 = e.raw_ind2;
          other2 = e.raw_ind1;
        }
        if (raw1 == other1) {
          return (raw2 < other2);
        }
        return raw1 < other1;
      }
      return (graph_ < e.get_graph_ptr());
    }

    /**
     * @brief Get the Euclidean distance of this edge
     *
     * @return double, the Euclidean distance of this edge
     */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /**
     * @brief Get a reference to the value of this edge.
     *
     * @return edge_value_type& reference to the value of this edge.
     */
    edge_value_type& value() {
      if (raw_ind1 < raw_ind2) {
        return const_cast<Graph*>(graph_)->edge_values[raw_ind1].at(raw_ind2);
      }
      return const_cast<Graph*>(graph_)->edge_values[raw_ind2].at(raw_ind1);
    }

    /**
     * @brief Get a const reference to the value of this edge.
     *
     * @return const edge_value_type&, the const reference to this edge's value.
     */
    const edge_value_type& value() const {
      if (raw_ind1 < raw_ind2) {
        return graph_->edge_values[raw_ind1].at(raw_ind2);
      }
      return graph_->edge_values[raw_ind2].at(raw_ind1);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects


    const Graph * graph_;

    // A unique index for an edge even after removal
    size_type raw_ind1;
    size_type raw_ind2;


    Edge(const Graph* g, size_type input_node1, size_type input_node2)
        : graph_(g), raw_ind1(input_node1), raw_ind2(input_node2) {}

    const Graph* get_graph_ptr() const { return graph_; }
    size_type get_ind1() const { return raw_ind1; }
    size_type get_ind2() const { return raw_ind2; }
    };


  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_count;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type cumu_edge_prev = 0;
    size_type cumu_edge = 0;
    for(size_type j = 0; j < node_count; ++j){
      unsigned raw_ind = ind_to_raw_ind[j];
      cumu_edge += edge_data[raw_ind].size();
      if(cumu_edge > i) {
        auto it = edge_data[raw_ind].begin();
        std::advance(it, i - cumu_edge_prev);
        return Edge(this, j, *it);
      }
      cumu_edge_prev = cumu_edge;

    }
    assert(false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if(!(has_node(a) && has_node(b))){
      return false;
    }
    size_type ind1 = a.index();
    size_type ind2 = b.index();
    if(ind2 == ind1){
      return false;
    }

    if((ind1 >= node_count) || ind2 >= node_count){
      return false;
    }
    size_type raw1 = a.get_raw_index();
    size_type raw2 = b.get_raw_index();
    if(raw1 < raw2){
      return (edge_data[raw1].count(raw2) > 0);
    }
    return (edge_data[raw2].count(raw1) > 0);
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
    assert((has_node(a) && has_node(b)));
    // HW0: YOUR CODE HERE
    size_type ind1 = a.index();
    size_type ind2 = b.index();
    assert(ind1 < node_count);
    assert(ind2 < node_count);
    assert(ind1 != ind2);

    if(has_edge(a, b)){
      return Edge(this, a.get_raw_index(), b.get_raw_index());
    }
    size_type raw1 = a.get_raw_index();
    size_type raw2 = b.get_raw_index();
    if(raw1 < raw2){
      edge_data[raw1].emplace(raw2);
      edge_data_transpose[raw2].emplace(raw1);
      edge_values[raw1].emplace(raw2, edge_value_type{});
    } else {
      edge_data[raw2].emplace(raw1);
      edge_data_transpose[raw1].emplace(raw2);
      edge_values[raw2].emplace(raw1, edge_value_type{});
    }
    edge_count++;
    ever_added_edge_count++;
    return Edge(this, raw1, raw2);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    edge_count = 0;
    ever_added_edge_count = 0;
    node_count = 0;
    ever_added_node_count = 0;
    ind_to_raw_ind.resize(0);
    raw_ind_to_ind.resize(0);
    edge_data.resize(0);
    edge_data_transpose.resize(0);
    edge_values.resize(0);
    location.resize(0);
    node_values.resize(0);
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
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

    /**
     * @brief Dereference this iterator to get the node it poitns to
     * @pre This iterator points to a valid node. In particular,
     * this->graph_ is not a nullptr and node_index < graph_->node_count
     * @return Node, the node the iterator currently points to.
     */
    Node operator*() const {
      return Node(graph_, graph_->get_raw_from_ind(node_index));
    }

    /**
     * @brief Increment this iterator to go to the next node.
     * @return NodeIterator&, a reference to this iterator after update.
     * @post This iterator moves to the next node.
     */
    NodeIterator& operator++() {
      node_index++;
      return *this;
    }

    /**
     * @brief compare if two NodeIterators refer to the same node.
     * @param other Another NodeIterator to compare to.
     * @return true if this iterator and other points to the same node.
     * @return false if this iterator and other do not point to the same node.
     */
    bool operator==(const NodeIterator& other) const {
      return (graph_ == other.graph_) && (node_index == other.node_index);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph * graph_;
    size_type node_index;
    NodeIterator(const Graph* graph, size_type ind)
        : graph_(graph), node_index(ind) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
   * @brief Get the iterator to the first node of this graph
   *
   * @return node_iterator at the first node
   */
  node_iterator node_begin() const { return NodeIterator(this, 0); }

  /**
   * @brief Get the iterator to the end of the nodes
   *
   * @return node_iterator that represents the end of all nodes.
   */
  node_iterator node_end() const { return NodeIterator(this, node_count); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>  {
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

    /**
     * @brief Dereference this IncidentIterator to get the corresponding edge
     * @pre This IncidentIterator must correspond to a valid edge
     * @return Edge that is referred by this IncidentIterator. The output e satisfies
     * that e.node1() is the node that spawns this iterator.
     */
    Edge operator*() const {
      size_type raw_ind2;
      size_type raw_ind = graph_->ind_to_raw_ind[node_index];
      if (m_it == graph_->edge_data[raw_ind].end()) {
        // if this happens ind2 < node_index
        raw_ind2 = *m_trans_it;
        return Edge(graph_, raw_ind, raw_ind2);
      }
      raw_ind2 = *m_it;
      return Edge(graph_, raw_ind, raw_ind2);
    }

    /**
     * @brief Increment this iterator to point to the next incident edge.
     * @pre This iterator is not at the end of the container.
     * @post The iterator is incremented to the next adjacent edge.
     * @return IncidentIterator&, a reference to this iterator after incrementation
     */
    IncidentIterator& operator++() {
      size_type raw_ind = graph_->ind_to_raw_ind[node_index];
      if (m_it == graph_->edge_data[raw_ind].end()) {
        m_trans_it++;
        return *this;
      }
      m_it++;
      return *this;
    }

    /**
     * @brief compare if two IncidentIterator refer to the same node.
     * @param other Another IncidentIterator to compare to.
     * @return true If both iterator refer to the same edge
     * @return false If the two iterators do not refer to the same edge.
     */
    bool operator==(const IncidentIterator& other) const {
      return (graph_ == other.graph_) && (node_index == other.node_index) &&
             (m_it == other.m_it) && (m_trans_it == other.m_trans_it);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph * graph_;
    size_type node_index;
    set_iterator m_it;
    set_iterator m_trans_it;
    IncidentIterator(const Graph* graph, size_type ind,
                     const set_iterator& m_it, const set_iterator& m_trans_it)
        : graph_(graph), node_index(ind), m_it(m_it), m_trans_it(m_trans_it) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
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

    /**
     * @brief Deference an EdgeIterator to get the edge this iterator refers to
     * @pre This must be an iterator to a valid edge
     * @return The edge this iterator refers to.
     */
    Edge operator*() const {
      size_type raw1 = graph_->ind_to_raw_ind[node_index];
      size_type raw2 = *m_it;
      return Edge(graph_, raw1, raw2);
    }

    /**
     * @brief Increment the iterator to the next edge.
     * @pre The iterator must not be at the end of all edges.
     * @post this iterator is incremented to the next edge.
     * @return The updated EdgeIterator.
     */
    EdgeIterator& operator++() {
      size_type raw = graph_->ind_to_raw_ind[node_index];
      if (m_it != graph_->edge_data[raw].end()) {
        m_it++;
      }
      while (m_it == graph_->edge_data[raw].end()) {
        if (node_index == (graph_->num_nodes() - 1)) {
          break;
        }
        node_index++;
        raw = graph_->ind_to_raw_ind[node_index];
        m_it = graph_->edge_data[raw].begin();
      }
      return *this;
    }

    /**
     * @brief Check if two EdgeIterators are the same.
     *
     * @param other Another EdgeIterator to compare to.
     * @return true If the two EdgeIterator refer to the same edge.
     * @return false If the two EdgeIterator do not refer to the same edge.
     */
    bool operator==(const EdgeIterator& other) const {
      return (graph_ == other.graph_) && (node_index == other.node_index) &&
             (m_it == other.m_it);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph * graph_;
    size_type node_index;
    set_iterator m_it;
    EdgeIterator(const Graph* graph, size_type ind, const set_iterator& m_it)
        : graph_(graph), node_index(ind), m_it(m_it) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
   * @brief Get an interator to the first edge in this graph.
   * @pre This graph must have at least one node.
   * @return edge_iterator to the first edge in the graph.
   */
  edge_iterator edge_begin() const {
    for (size_type i = 0; i < node_count; ++i) {
      size_type raw = ind_to_raw_ind[i];
      if (edge_data[raw].size() > 0) {
        return EdgeIterator(this, i, edge_data[raw].begin());
      }
    }
    size_type raw = ind_to_raw_ind[node_count - 1];
    return EdgeIterator(this, node_count - 1, edge_data[raw].end());
  }

  /**
   * @brief Get an interator to the last edge in this graph.
   * @pre This graph must have at least one node.
   * @return edge_iterator to the last edge in the graph.
   */
  edge_iterator edge_end() const {
    size_type raw = ind_to_raw_ind[node_count - 1];
    return EdgeIterator(this, node_count-1, edge_data[raw].end());
  }

  // HW2 code here

  /**
   * @brief Function that removes the edge (n1, n2) in the graph, if exists
   *
   * @param n1 One endpoint of the edge to be removed
   * @param n2 The other endpoint of the edge to be removed
   * @post If the function returns 1, the the ege (n1, n2) is removed. All nodes
   * and node itreators are still valid. All other edges are still valid. Edge
   * iterators whose endpoints is one of (n1, n2) are invalidated. Incident
   * iterators of nodes other than n1 and n2 are still valid. If the function
   * returns 0, no edge is removed and all objects are valid.
   *
   * Complexity: O(log(n1.degree()) + log(n2.degree()))
   * @return size_type, 1 if an edge is removed, 0 if the edge does not exist in
   * the graph.
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (has_edge(n1, n2)) {
      size_type ind1 = n1.get_raw_index();
      size_type ind2 = n2.get_raw_index();
      size_type tmp = ind1;
      // make sure ind1 < ind2
      if (ind2 < ind1) {
        ind1 = ind2;
        ind2 = tmp;
      }
      edge_data[ind1].erase(ind2);
      edge_values[ind1].erase(ind2);
      edge_data_transpose[ind2].erase(ind1);
      edge_count--;
      return 1;
    }
    return 0;
  }

  /**
   * @brief Function that removes the edge e in the graph, if exists
   *
   * @param e The edge to be removed, if it exists
   * @post If the function returns 1, the the ege (n1, n2) is removed. All nodes
   * and node itreators are still valid. All other edges are still valid. Edge
   * iterators whose endpoints is one of (n1, n2) are invalidated. Incident
   * iterators of nodes other than n1 and n2 are still valid. If the function
   * returns 0, no edge is removed and all objects are valid.
   *
   * Complexity: O(log(n1.degree()) + log(n2.degree())), where e = (n1, n2)
   * @return size_type. 1 if the edge is removed, 0 if the edge does not exist
   * in the graph.
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /**
   * @brief Function that removes the edge from an edge_iterator e_it, if the
   * edge exists.
   *
   * @param e_it Iterator to the edge to be removed.
   * @post  If the function returns 1, the the edge (n1, n2) is removed, where
   * *e_it = (n1, n2) All nodes and node itreators are still valid. All other
   * edges are still valid. Edge iterators whose endpoints is one of (n1, n2)
   * are invalidated. Incident iterators of nodes other than n1 and n2 are still
   * valid. If the function returns 0, no edge is removed and all objects are
   * valid.
   * @post If an edge is removed the input iterator e_it is no longer valid, but
   * the returned iterator is always valid (possiblt the end).
   *
   * Complexity: O(log(n1.degree()) + log(n2.dergee())), where *e_it = (n1, n2),
   * assuming the number of isolated nodes is O(1).
   * @return edge_iterator An edge iterator to a valid edge after e_it.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    if(e_it == this->edge_end()){
      return e_it;
    }
    size_type node_index = e_it.node_index;
    size_type ind1 = ind_to_raw_ind[node_index];
    size_type ind2 = (*e_it).node2().get_raw_index();
    set_iterator new_it = edge_data[ind1].erase(e_it.m_it);
    edge_values[ind1].erase(ind2);
    edge_data_transpose[ind2].erase(ind1);
    while (new_it == edge_data[ind1].end()) {
      if (ind1 == (node_count - 1)) {
        break;
      }
      node_index++;
      ind1 = ind_to_raw_ind[node_index];
      new_it = edge_data[ind1].begin();
    }
    edge_count--;
    return EdgeIterator(this, node_index, new_it);
  }

  /**
   * @brief Function that removes a node and all its adjacent edges from the
   * graph, if the node exists.
   *
   * @param n The node to be removed
   * @post If the function returns 1, then the node @a n and all of its adjacent
   * edges are removed from the graph. All nodes other than @a n are still
   * valid, but the node iterators to these node might be reordered. All edges
   * that are not adjacent to @a n are still valid, so are their corresponding
   * edge_iterators. edge_iterators adjacent to @a n are invalidated.
   * Incident iterators for nodes that are not adjacent to @a n are still valid.
   *
   * Complexity: Let M = max_{n is a node in this graph} n.degree(). Then the
   * complexity is O(Mlog(M)). If the degree of nodes are uniformly bounded,
   * then the complexity is O(1).
   * @return size_type 1 if the node is removed. 0 if the node does not exists.
   */
  size_type remove_node(const Node& n) {
    if (has_node(n)) {
      size_type raw_ind = n.get_raw_index();
      // First remove edges attached to this node
      for(auto it = edge_data[raw_ind].begin(); it != edge_data[raw_ind].end(); ++it){
        size_type ind2 = *it;
        edge_data_transpose[ind2].erase(raw_ind);
        edge_count--;
      }

      for(auto it = edge_data_transpose[raw_ind].begin(); it != edge_data_transpose[raw_ind].end(); ++it){
        size_type ind2 = *it;
        edge_data[ind2].erase(raw_ind);
        edge_count--;
      }

      // Update ind <-> raw_ind mapping
      size_type ind = n.index();
      raw_ind_to_ind[raw_ind] = static_cast<size_type>(-1);
      ind_to_raw_ind[ind] = ind_to_raw_ind[node_count - 1];
      raw_ind_to_ind[ind_to_raw_ind[ind]] = ind;
      ind_to_raw_ind.pop_back();
      node_count--;
      return 1;
    }
    return 0;
  }

  /**
   * @brief Function that removes a node and all its adjacent edges from the
   * graph, if the node exists. The input is an iterator to the node.
   *
   * @param n_it Iterator to the node to be removed.
   * @post If the function returns 1, then the node @a n_it and all of its
   * adjacent edges are removed from the graph. All nodes other than @a n_it are
   * still valid, but the node iterators to these node might be reordered. All
   * edges that are not adjacent to @a n_it are still valid, so are their
   * corresponding edge_iterators. edge_iterators adjacent to @a n_it are
   * invalidated. Incident iterators for nodes that are not adjacent to @a n_it
   * are still valid.
   *
   * Complexity: Let M = max_{n is a node in this graph} n.degree(). Then the
   * complexity is O(Mlog(M)). If the degree of nodes are uniformly bounded,
   * then the complexity is O(1).
   * @return node_iterator An iterator to a valid node after removal, possibly
   * end.
   */
  node_iterator remove_node(node_iterator n_it){
    if(n_it == node_end()){
      return n_it;
    }
    size_type node_index = (*n_it).index();
    remove_node(*n_it);
    return NodeIterator(this, node_index);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
