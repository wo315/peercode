#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

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

  // stores node positions
  std::vector<Point> nodes_;
  // stores number of edges
  unsigned int num_edges_;
  // edge map
  // nested map of node id -> connected node ids -> edge id
  std::unordered_map<unsigned int, 
                     std::unordered_map<unsigned int, 
                                        unsigned int>> nodeid_map;
  // map edge id to node pair
  std::unordered_map<unsigned int, 
                     std::pair<unsigned int, 
                               unsigned int>> edgeid_map;

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
    num_edges_ = 0;
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
      // HW0: YOUR CODE HERE
      graph_ = nullptr;
      nodeid_ = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // if node is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Node is invalid." << std::endl;
        exit(0);
      }

      return this->graph_->nodes_[nodeid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // if node is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Node is invalid." << std::endl;
        exit(0);
      }
      // if node id is invalid, throw error
      if (this->graph_->nodes_.size() < 0 or
          this->nodeid_ >= this->graph_->nodes_.size()) {
        std::cerr << "Error: Node index " << this->nodeid_ << " is invalid.\n"
          << "Graph size is " << this->graph_->nodes_.size()
          << " nodes." <<std::endl;
        exit(0);
      }

      return this->nodeid_;
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
      (void) n; // Quiet compiler warning
      if (this->graph_ == n.graph_ and this->nodeid_ == n.nodeid_) {
        return true;
      }
      return false;
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
      (void) n; // Quiet compiler warning
      if (this->nodeid_ < n.nodeid_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Pointer to Graph 
    graph_type* graph_;
    // ID of node in Graph
    size_type nodeid_;
    
    /** Private Constructor */
    Node(const graph_type* graph, size_type nodeid)
        : graph_(const_cast<graph_type*>(graph)), nodeid_(nodeid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return this->nodes_.size();
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
    (void) position; // Quiet compiler warning
    this->nodes_.push_back(position);
    return Node(this, this->nodes_.size() - 1); // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n; // Quiet compiler warning
    if (this == n.graph_ and this->nodes_.size() > n.nodeid_) {
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
    // HW0: YOUR CODE HERE
    (void) i; // Quiet compiler warning
    // if index > graph size, throw error
    if (i >= this->nodes_.size()) {
      std::cerr << "Error: Invalid index.\n" 
        << "Index: " << i << "\n"
        << "Graph Size: " << this->nodes_.size() << std::endl;
      exit(0);
    }

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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      graph_ = nullptr;
      edgeid_ = 0;
      nodeids_.first = 0;
      nodeids_.second = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // if edge is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Edge is invalid." << std::endl;
        exit(0);
      }

      return Node(this->graph_, this->nodeids_.first);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // if edge is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Edge is invalid." << std::endl;
        exit(0);
      }

      return Node(this->graph_, this->nodeids_.second);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e; // Quiet compiler warning
      //HW0: YOUR CODE HERE
      if (this->graph_ != e.graph_) {
        return false;
      }

      if (this->nodeids_.first == e.nodeids_.first and 
          this->nodeids_.second == e.nodeids_.second) {
        return true;
      }
      else if (this->nodeids_.first == e.nodeids_.second and 
               this->nodeids_.second == e.nodeids_.first) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e; // Quiet compiler warning
      //HW0: YOUR CODE HERE
      if (this->edgeid_ < e.edgeid_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // pointer to graph
    graph_type* graph_;
    // edge id
    size_type edgeid_;
    // pair of IDs for nodes in edge
    std::pair<size_type, size_type> nodeids_;

    /** Private Constructor */
    Edge(const graph_type* graph, 
         size_type edgeid, 
         std::pair<size_type, size_type> nodeids)
        : graph_(const_cast<graph_type*>(graph)), 
          edgeid_(edgeid),
          nodeids_(nodeids) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i; // Quiet compiler warning
    // if index is not less than number of edges, throw error
    if (i < 0 or i >= this->num_edges_) {
      std::cerr << "Error: Invalid index.\n"
      << "Index: " << i << "\n"
      << "Number of edges in graph: " << this->num_edges_ << std::endl;
      exit(0);
    }

    return Edge(this, i, this->edgeid_map.at(i));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // if node a or b is invalid, throw an error
    (void) a; (void) b; // Quiet compiler warning
    if (a.graph_ == nullptr) {
      std::cerr << "Error: Node a is invalid." << std::endl;
      exit(0);
    }
    if (b.graph_ == nullptr) {
      std::cerr << "Error: Node b is invalid." << std::endl;
      exit(0);
    }

    if (this->nodeid_map.count(a.nodeid_) == 1) {
      if (this->nodeid_map.at(a.nodeid_).count(b.nodeid_) == 1) {
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
    // HW0: YOUR CODE HERE
    (void) a, (void) b; // Quiet compiler warning
    // if node a or b is invalid, throw an error
    if (a.graph_ == nullptr) {
      std::cerr << "Error: Node a is invalid." << std::endl;
      exit(0);
    }
    if (b.graph_ == nullptr) {
      std::cerr << "Error: Node b is invalid." << std::endl;
      exit(0);
    }

    // check if edge exists in graph; if so, return that edge
    if (this->has_edge(a, b)) {
      size_type edgeid = this->nodeid_map[a.nodeid_][b.nodeid_];
      std::pair<size_type, size_type> nodeids = this->edgeid_map[edgeid];
      return Edge(this, edgeid, nodeids);
    }

    // add edge to node id map of graph
    this->nodeid_map[a.nodeid_].insert({b.nodeid_, this->num_edges_});
    this->nodeid_map[b.nodeid_].insert({a.nodeid_, this->num_edges_});

    // add edge to edge id map of graph
    size_type edgeid = this->num_edges_;
    std::pair<size_type, size_type> nodes;
    nodes.first = a.nodeid_;
    nodes.second = b.nodeid_;
    this->edgeid_map.insert({edgeid, nodes});

    // increment num_edges_
    this->num_edges_ += 1;
    
    // return edge object
    return Edge(this, edgeid, nodes);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
  // HW0: YOUR CODE HERE
  this->nodes_.clear();
  this->num_edges_ = 0;
  this->nodeid_map.clear();
  this->edgeid_map.clear();
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
