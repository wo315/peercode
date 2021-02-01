#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>

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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
  // A vector holding the locations of the nodes
  std::vector<Point> location;

  // If a node is the nth one that's ever added to the graph,
  // and it is currently the ith node in the graph,
  // Then ind_to_raw_ind[i] = n
  std::vector<unsigned> ind_to_raw_ind;

  // Current number of nodes
  unsigned node_count;

  // Number of nodes that are ever added to this graph
  unsigned ever_added_node_count;

  // Current number of edges
  unsigned edge_count;

  // Number of edges that's ever added to this graph
  unsigned ever_added_edge_count;

  // If (i, j) is an edge in this graph with i < j
  // Then edge_data[i] has key j and value 
  // equal to the number of edge ever added when
  // (i, j) is added (its raw_index)
  std::vector<std::map<unsigned,unsigned>> edge_data;

  // @return the location of the node at index @a ind
  const Point& get_location(unsigned ind) const {
    assert(ind < location.size());
    return location[ind];
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
    Node() : graph_(nullptr) {
      // HW0: YOUR CODE HERE
      raw_index = static_cast<size_type>(-1);
      ind = static_cast<size_type>(-1);
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(valid_node());
      return graph_->get_location(ind);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(valid_node());
      return ind;
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
    // Potential index of this node in the graph
    size_type ind;

    Node(const Graph* g, size_type raw, size_type ind)
        : graph_(g), raw_index(raw), ind(ind) {}

    const Graph* get_graph_ptr() const { return this->graph_; }
    size_type get_raw_index() const { return this->raw_index; }

    //@return True if this Node is a valid node in a graph
    bool valid_node() const {
      return (graph_ != nullptr) && (ind < graph_->num_nodes()) &&
             (graph_->get_raw_from_ind(ind) == raw_index);
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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    location.emplace_back(position.x, position.y, position.z);
    edge_data.emplace_back();
    ind_to_raw_ind.push_back(ever_added_node_count);
    node_count++;
    ever_added_node_count++;
    return Node(this, ever_added_node_count - 1, node_count-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.index() >= node_count) {
      return false;
    }
    if (ind_to_raw_ind[n.index()] != n.get_raw_index()) {
      return false;
    }
    return (this == n.get_graph_ptr());
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
    return Node(this, ind_to_raw_ind[i], i);
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
    Edge() : graph_(nullptr) {
      // HW0: YOUR CODE HERE
      raw_index = static_cast<size_type>(-1);
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      if(raw_index == static_cast<size_type>(-1)){
        return Node();
      }
      return Node(graph_,graph_->get_raw_from_ind(ind1) ,ind1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if(raw_index == static_cast<size_type>(-1)){
        return Node();
      }
      return Node(graph_,graph_->get_raw_from_ind(ind2) ,ind2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
                 // Quiet compiler warning
      //HW0: YOUR CODE HERE
      if(graph_ != e.get_graph_ptr()){
        return false;
      }
      return (raw_index == e.get_raw_index());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
                // Quiet compiler warning
      //HW0: YOUR CODE HERE
      if(graph_ == e.get_graph_ptr()){
        return raw_index < e.get_raw_index();
      }
      return (graph_ < e.get_graph_ptr());
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
    size_type raw_index;
    size_type ind1;
    size_type ind2;


    Edge(const Graph* g, size_type input_node1, size_type input_node2,
         size_type raw)
        : graph_(g), raw_index(raw), ind1(input_node1), ind2(input_node2) {}

    const Graph* get_graph_ptr() const { return graph_; }
    size_type get_ind1() const { return ind1; }
    size_type get_ind2() const { return ind2; }
    size_type get_raw_index() const { return raw_index; }
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
      cumu_edge += edge_data[j].size();
      if(cumu_edge > i) {
        auto it = edge_data[j].begin();
        std::advance(it, i - cumu_edge_prev);
        return Edge(this, j, it->first, it->second);
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
    if(ind1 < ind2){
      return (edge_data[ind1].count(ind2) > 0);
    }
    return (edge_data[ind2].count(ind1) > 0);
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
      size_type raw_index;
      if(ind1 < ind2){
        raw_index = edge_data[ind1].at(ind2);
      } else {
        raw_index = edge_data[ind2].at(ind1);
      }
      return Edge(this, ind1, ind2, raw_index);
    }
    if(ind1 < ind2){
      edge_data[ind1].emplace(ind2, ever_added_edge_count);
    } else {
      edge_data[ind2].emplace(ind1, ever_added_edge_count);
    }
    edge_count++;
    ever_added_edge_count++;
    return Edge(this, ind1, ind2, ever_added_edge_count-1);
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
    edge_data.resize(0);
    location.resize(0);
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

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
