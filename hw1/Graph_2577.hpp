#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <unordered_map>
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
  /** Synonym for V. */
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
  Graph() {
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
  class Node: private totally_ordered <Node> {
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
      node_id_ = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      // if node is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Node is invalid." << std::endl;
        exit(0);
      }

      return this->graph_->nodes_[node_id_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // if node is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Node is invalid." << std::endl;
        exit(0);
      }
      // if node id is invalid, throw error
      if (this->node_id_ >= this->graph_->nodes_.size()) {
        std::cerr << "Error: Node index " << this->node_id_ << " is invalid.\n"
                  << "Graph size is " << this->graph_->nodes_.size()
                  << " nodes." <<std::endl;
        exit(0);
      }

      return this->node_id_;
    }

    /** Return this node's value, of type V. */
    node_value_type& value() {
      return this->graph_->nodes_[this->node_id_].second;
    }
    /** Return this node's value, of type V; const version */
    const node_value_type& value() const {
      return this->graph_->nodes_[this->node_id_].second;
    }
    
    /** Return number of edges connected to this node */
    size_type degree() const {
      if (this->graph_ == nullptr || this->graph_->node_id_map_.size() == 0) {
        return 0;
      }

      return this->graph_->node_id_map_[this->node_id_].size();
    }
    
    /** Returns iterator at beginning of map of adjacent nodes */
    incident_iterator edge_begin() const {
      return IncidentIterator(this->graph_, this->node_id_, 
                              this->graph_->node_id_map_[this->node_id_]
                                .begin());
    }
    /** Returns iterator at end of map of adjacent nodes */
    incident_iterator edge_end() const {
      return IncidentIterator(this->graph_, this->node_id_, 
                              this->graph_->node_id_map_[this->node_id_]
                                .end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (this->graph_ == n.graph_ && this->node_id_ == n.node_id_) {
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
      if (this->graph_ < n.graph_) {
        return true;
      }
      if (this->graph_ == n.graph_) {
        if (this->node_id_ < n.node_id_) {
          return true;
        }
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Pointer to Graph 
    graph_type* graph_;
    // ID of node in Graph
    size_type node_id_;
    
    /** Private Constructor */
    Node(const graph_type* graph, size_type nodeid)
        : graph_(const_cast<graph_type*>(graph)), node_id_(nodeid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
  Node add_node(const Point& position, 
                const node_value_type& v = node_value_type()) {
    std::pair<Point, node_value_type> node_pos_val;
    node_pos_val.first = position;
    node_pos_val.second = v;
    this->nodes_.push_back(node_pos_val);

    return Node(this, this->nodes_.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (this == n.graph_ && this->nodes_.size() > n.node_id_) {
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
  class Edge: private totally_ordered <Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      graph_ = nullptr;
      edge_id_ = 0;
      node_ids_.first = 0;
      node_ids_.second = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // if edge is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Edge is invalid." << std::endl;
        exit(0);
      }
      return Node(this->graph_, this->node_ids_.first);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // if edge is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Edge is invalid." << std::endl;
        exit(0);
      }
      return Node(this->graph_, this->node_ids_.second);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ == e.graph_ && this->edge_id_ == e.edge_id_) {
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
      if (this->graph_ < e.graph_) {
        return true;
      }
      if (this->graph_ == e.graph_) {
        if (this->edge_id_ < e.edge_id_) {
          return true;
        }
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // pointer to graph
    graph_type* graph_;
    // edge id
    size_type edge_id_;
    // pair of IDs for nodes in edge
    std::pair<size_type, size_type> node_ids_;

    /** Private Constructor */
    Edge(const graph_type* graph, 
         size_type edgeid, 
         std::pair<size_type, size_type> node_ids)
        : graph_(const_cast<graph_type*>(graph)), 
          edge_id_(edgeid),
          node_ids_(node_ids) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return this->num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    // if index is not less than number of edges, throw error
    if (i >= this->num_edges_) {
      std::cerr << "Error: Invalid index.\n"
                << "Index: " << i << "\n"
                << "Number of edges in graph: "
                << this->num_edges_ << std::endl;
      exit(0);
    }
    return Edge(this, i, this->edge_id_map_.at(i));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // if node a or b is invalid, throw an error
    if (a.graph_ == nullptr) {
      std::cerr << "Error: Node a is invalid." << std::endl;
      exit(0);
    }
    if (b.graph_ == nullptr) {
      std::cerr << "Error: Node b is invalid." << std::endl;
      exit(0);
    }

    if (this->node_id_map_.count(a.node_id_) == 1) {
      if (this->node_id_map_.at(a.node_id_).count(b.node_id_) == 1) {
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
      size_type edgeid = this->node_id_map_[a.node_id_][b.node_id_];
      std::pair<size_type, size_type> nodes;
      nodes.first = a.node_id_;
      nodes.second = b.node_id_;

      return Edge(this, edgeid, nodes);
    }

    // add edge to node id map of graph
    this->node_id_map_[a.node_id_].insert({b.node_id_, this->num_edges_});
    this->node_id_map_[b.node_id_].insert({a.node_id_, this->num_edges_});

    // add edge to edge id map of graph
    size_type edgeid = this->num_edges_;
    std::pair<size_type, size_type> nodes;
    nodes.first = a.node_id_;
    nodes.second = b.node_id_;
    this->edge_id_map_.insert({edgeid, nodes});

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
  this->nodes_.clear();
  this->num_edges_ = 0;
  this->node_id_map_.clear();
  this->edge_id_map_.clear();
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
      graph_ = nullptr;
      node_id_ = 0;
      graph_len_ = 0;
    }

    /** Dereference operator for NodeIterator
     * Returns a Node
     */
    value_type operator*() const {
      if (this->graph_ == nullptr || this->node_id_ >= this->graph_len_) {
        std::cerr << "Error: Attempting to dereference to invalid node." 
                  << std::endl;
        exit(0);   
      }

      return Node(this->graph_, this->node_id_);
    }

    /** Increment operator for NodeIterator
     * Returns a NodeIterator reference
     */
    NodeIterator& operator++() {
      if (this->node_id_ == this->graph_len_) {
        std::cerr << "Error: Attempting to increment at end." << std::endl;
        exit(0);
      }

      this->node_id_ += 1;
      return (*this);
    }

    /** Is Equal operator for NodeIterator
     * Returns a boolean
     */
    bool operator==(const NodeIterator& node_iter) const {
      if (this->graph_ == node_iter.graph_
          && this->node_id_ == node_iter.node_id_) {
            return true;
          }
      return false;
    }

    /** Not Equal operator for NodeIterator
     * Returns a boolean
     */
    bool operator!=(const NodeIterator& node_iter) const {
      return !((*this) == node_iter);
    }

   private:
    friend class Graph;
    // Graph
    graph_type* graph_;
    // Node id
    size_type node_id_;
    // Length of node id vector
    size_type graph_len_;

    // private constructor
    NodeIterator(const graph_type* graph, size_type node_id) 
      : graph_(const_cast<graph_type*>(graph)), node_id_(node_id) {
      graph_len_ = this->graph_->nodes_.size();
    }
  };

  /** Returns iterator at beginning of node vector */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  /** Returns iterator at end of node vector */
  node_iterator node_end() const {
    return NodeIterator(this, this->num_nodes());
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
      this->graph_ = nullptr;
      this->cent_node_ = 0;
    }

    /** Dereference operator for NodeIterator
     * Returns an Edge
     */
    Edge operator*() const {
      if (this->graph_ == nullptr || 
          this->cent_node_ >= this->graph_->num_nodes()) {
        std::cerr << "Error: Attempting to dereference to invalid edge."
                  << std::endl;
        exit(0);
      }

      std::pair<size_type, size_type> nodes;
      nodes.first = this->cent_node_;

      // get nodes associated with current edge
      std::pair<size_type, size_type> edge_node_pair = 
        this->graph_->edge_id_map_[(*(this->map_iterator_)).second];

      // find other node from mapping in graph
      if (edge_node_pair.first == this->cent_node_) {
        nodes.second = edge_node_pair.second;
      }
      else {
        nodes.second = edge_node_pair.first;
      }

      return Edge(this->graph_, (*(this->map_iterator_)).second, nodes);
    }

    /** Increment operator for NodeIterator
     * Returns an IncidentIterator reference
     */
    IncidentIterator& operator++() {
      ++(this->map_iterator_);
      return (*this);
    }

    /** Is Equal operator for NodeIterator
     * Returns a boolean
     */
    bool operator==(const IncidentIterator& inc_iter) const {
      if (this->graph_ == inc_iter.graph_
          && this->cent_node_ == inc_iter.cent_node_
          && this->map_iterator_ == inc_iter.map_iterator_) {
            return true;
      }
      return false;
    }

    /** Not Equal operator for NodeIterator
     * Returns a boolean
     */
    bool operator!=(const IncidentIterator& inc_iter) const {
      return !((*this) == inc_iter);
    }

   private:
    friend class Graph;
    // graph
    graph_type* graph_;
    // central node
    size_type cent_node_;
    // const_iterator for unordered map
    std::unordered_map<size_type, size_type>::const_iterator map_iterator_;
    // private constructor
    IncidentIterator(const graph_type* graph, 
                     size_type cent_node, 
                     std::unordered_map<size_type, size_type>::const_iterator
                      map_iterator)
                     : graph_(const_cast<graph_type*>(graph)),
                       cent_node_(cent_node), 
                       map_iterator_(map_iterator) {
    }
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
      graph_ = nullptr;
      edge_id_ = 0;
      graph_len_ = 0;
    }

    /** Dereference operator for EdgeIterator
     * Returns an Edge
     */
    Edge operator*() const {
      if (this->graph_ == nullptr || 
          this->edge_id_ >= this->graph_len_) {
        std::cerr << "Error: Attempting to dereference to invalid edge."
                  << std::endl;
        exit(0);
      }
      std::pair<size_type, size_type> nodes;
      nodes = this->graph_->edge_id_map_[this->edge_id_];
      
      return Edge(this->graph_, this->edge_id_, nodes);
    }

    /** Increment operator for EdgeIterator
     * Returns an EdgeIterator reference
     */
    EdgeIterator& operator++() {
      if (this->edge_id_ == this->graph_len_) {
        std::cerr << "Error: Attempting to increment at end." << std::endl;
        exit(0);
      }
      
      edge_id_ += 1;
      return (*this);
    }

    /** Is Equal operator for EdgeIterator
     * Returns a boolean
     */
    bool operator==(const EdgeIterator& edge_iter) const {
      if (this->graph_ == edge_iter.graph_
          && this->edge_id_ == edge_iter.edge_id_) {
            return true;
      }
      return false;
    }

    /** Not Equal operator for EdgeIterator
     * Returns a boolean
     */
    bool operator!=(const EdgeIterator& edge_iter) const {
      return !((*this) == edge_iter);
    }

   private:
    friend class Graph;
    // graph
    graph_type* graph_;
    // edge id
    size_type edge_id_;
    // number of edges in graph
    size_type graph_len_;

    // private constructor
    EdgeIterator(const graph_type* graph, size_type edge_id) 
      : graph_(const_cast<graph_type*>(graph)), edge_id_(edge_id) {
      graph_len_ = graph_->edge_id_map_.size();
    }    
  };

  /** Returns iterator at beginning of edge map */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  /** Returns iterator at end of edge map */
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->edge_id_map_.size());
  }

 private:
  // stores node positions and values
  std::vector<std::pair<Point, node_value_type>> nodes_;

  // stores number of edges
  size_type num_edges_;

  // edge map
  // nested map of node id -> connected node ids -> edge id
  std::unordered_map<size_type, 
                     std::unordered_map<size_type, 
                                        size_type>> node_id_map_;
  // map edge id to node pair
  std::unordered_map<size_type, 
                     std::pair<size_type, 
                               size_type>> edge_id_map_;
};

#endif // CME212_GRAPH_HPP
