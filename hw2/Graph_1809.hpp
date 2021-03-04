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
template <typename V = int, typename E = int>
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
  typedef V node_value_type;
  /** Synonym for E. */
  typedef E edge_value_type;

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
    num_nodes_ = 0;
    num_node_deletes_ = 0;
    num_edge_deletes_ = 0;
    num_node_adds_ = 0;
    num_edge_adds_ = 0;

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
      node_idx_ = 0;
    }

    /** Return this node's position; const. */
    const Point& position() const {
      // if node is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Node is invalid." << std::endl;
        exit(0);
      }

      return this->graph_->nodes_.at(this->get_uid()).position_;
    }

    /** Return this node's position; non-const. */
    Point& position() {
      // if node is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Node is invalid." << std::endl;
        exit(0);
      }

      return this->graph_->nodes_[this->get_uid()].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // if node is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Node is invalid." << std::endl;
        exit(0);
      }
      // if node id is invalid, throw error
      if (this->node_idx_ >= this->graph_->size()) {
        std::cerr << "Error: Node index " << this->node_idx_ << " is invalid.\n"
                  << "Graph size is " << this->graph_->size()
                  << " nodes." <<std::endl;
        exit(0);
      }

      return this->node_idx_;
    }

    /** Return this node's value, of type V. */
    node_value_type& value() {
      return this->graph_->nodes_[this->get_uid()].value_;
    }
    /** Return this node's value, of type V; const version */
    const node_value_type& value() const {
      return this->graph_->nodes_[this->get_uid()].value_;
    }
    
    /** Return number of edges connected to this node */
    size_type degree() const {
      if (this->graph_ == nullptr || this->graph_->node_id_map_.size() == 0) {
        return 0;
      }

      return this->graph_->nodes_[this->get_uid()].degree_;
    }

    /** Return this node's unique id */
    size_type get_uid() const {
      return this->graph_->node_idx_uid_map_[this->node_idx_];
    }
    
    /** Returns iterator at beginning of map of adjacent nodes */
    incident_iterator edge_begin() const {
      return IncidentIterator(this->graph_, this->node_idx_, 
                              this->graph_->node_id_map_[this->get_uid()]
                                .begin());
    }
    /** Returns iterator at end of map of adjacent nodes */
    incident_iterator edge_end() const {
      return IncidentIterator(this->graph_, this->node_idx_, 
                              this->graph_->node_id_map_[this->get_uid()]
                                .end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (this->graph_ == n.graph_ && this->node_idx_ == n.index()) {
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
        if (this->node_idx_ < n.node_idx_) {
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
    size_type node_idx_;
    
    /** Private Constructor */
    Node(const graph_type* graph, size_type node_idx)
        : graph_(const_cast<graph_type*>(graph)), node_idx_(node_idx) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->num_nodes_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return this->size();
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

    size_type curr_uid = num_node_adds_;
    const size_type node_idx = this->num_nodes_;

    // add node to nodes_
    this->nodes_.insert({curr_uid, 
                         NodeData(position, v, node_idx, curr_uid, true, 0)});
    
    // add node to node_id_map
    this->node_id_map_.insert({curr_uid, {}});

    // add node to node_idx_uid_map_
    this->node_idx_uid_map_.insert({node_idx, curr_uid});

    // increment num_nodes_
    this->num_nodes_ += 1;
    this->num_node_adds_ += 1;

    return Node(this, node_idx);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (this == n.graph_ && this->num_nodes_ > n.node_idx_) {
      return this->nodes_.at(n.get_uid()).valid_;
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
    if (i >= this->num_nodes_) {
      std::cerr << "Error: Invalid index.\n" 
                << "Index: " << i << "\n"
                << "Graph Size: " << this->num_nodes_ << std::endl;
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
      edge_idx_ = 0;
      node1_idx_ = 0;
      node2_idx_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // if edge is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Edge is invalid." << std::endl;
        exit(0);
      }
      return Node(this->graph_, this->node1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // if edge is invalid, throw error
      if (this->graph_ == nullptr) {
        std::cerr << "Error: Edge is invalid." << std::endl;
        exit(0);
      }
      return Node(this->graph_, this->node2_idx_);
    }

    /** Return the Euclidean distance between the edge's nodes */
    double length() {
      Point position1 = this->node1().position();
      Point position2 = this->node2().position();

      return norm(position1 - position2);
    }

    /** Return this edge's value, of type E. */
    edge_value_type& value() {
      return this->graph_->edges_[this->get_uid()].value_;
    }
    /** Return this edge's value, of type E; const version */
    const edge_value_type& value() const {
      return this->graph_->edges_[this->get_uid()].value_;
    }

    /** Return this edge's unique id */
    size_type get_uid() const {
      return this->graph_->edge_idx_uid_map_[this->edge_idx_];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ == e.graph_ && this->edge_idx_ == e.edge_idx_) {
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
        if (this->edge_idx_ < e.edge_idx_) {
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
    size_type edge_idx_;
    // node1 id
    size_type node1_idx_;
    // node2 id
    size_type node2_idx_;

    /** Private Constructor */
    Edge(const graph_type* graph, 
         size_type edge_idx,
         size_type node1_idx,
         size_type node2_idx)
        : graph_(const_cast<graph_type*>(graph)), 
          edge_idx_(edge_idx),
          node1_idx_(node1_idx),
          node2_idx_(node2_idx) {
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
    size_type node1_idx = this->edges_.at(
                            this->edge_idx_uid_map_.at(i)).node1_idx_;
    size_type node2_idx = this->edges_.at(
                            this->edge_idx_uid_map_.at(i)).node2_idx_;
    return Edge(this, i, node1_idx, node2_idx);
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

    if (!(this->has_node(a)) || !(this->has_node(b))) {
      return false;
    }

    size_type a_uid = a.get_uid();
    size_type b_uid = b.get_uid();
    if (this->node_id_map_.count(a_uid) == 1) {
      if (this->node_id_map_.at(a_uid).count(b_uid) == 1) {
        size_type edge_uid = this->node_id_map_.at(a_uid).at(b_uid);
        if (this->edges_.at(edge_uid).valid_) {
          return true;
        }
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
   * Complexity: O(1)
   */
  Edge add_edge(const Node& a, const Node& b, 
                const edge_value_type& e = edge_value_type()) {
    // if node a or b is invalid, throw an error
    if (a.graph_ == nullptr) {
      std::cerr << "Error: Node a is invalid." << std::endl;
      exit(0);
    }
    if (b.graph_ == nullptr) {
      std::cerr << "Error: Node b is invalid." << std::endl;
      exit(0);
    }

    size_type a_index = a.index();
    size_type b_index = b.index();

    // check if edge exists in graph; if so, return that edge
    if (this->has_edge(a, b)) {
      size_type edge_uid = this->node_id_map_.at(a.get_uid()).at(b.get_uid());
      size_type edge_idx = this->edges_[edge_uid].index_;

      return Edge(this, edge_idx, a_index, b_index);
    }

    size_type curr_uid = num_edge_adds_;

    size_type edge_idx = this->num_edges_;

    size_type a_uid = a.get_uid();
    size_type b_uid = b.get_uid();

    // add edge to node uid map of graph
    if (this->node_id_map_[a_uid].count(b_uid) == 1) {
      this->node_id_map_[a_uid][b_uid] = curr_uid;
      this->node_id_map_[b_uid][a_uid] = curr_uid;
    }
    else {
      this->node_id_map_[a_uid].insert({b_uid, curr_uid});
      this->node_id_map_[b_uid].insert({a_uid, curr_uid});
    }

    // add edge to edges_
    this->edges_.insert({curr_uid, EdgeData(a_index, b_index, 
                                            e, edge_idx, 
                                            curr_uid, true)});

    // add edge to edge_idx_uid_map_
    this->edge_idx_uid_map_.insert({edge_idx, curr_uid});

    // increment degree of nodes
    this->nodes_[a.get_uid()].degree_ += 1;
    this->nodes_[b.get_uid()].degree_ += 1;

    // increment num_edges_
    this->num_edges_ += 1;
    this->num_edge_adds_ += 1;

    // return edge object
    return Edge(this, edge_idx, a_index, b_index);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    this->nodes_.clear();
    this->num_nodes_ = 0;
    this->node_id_map_.clear();
    this->edges_.clear();
    this->num_edges_ = 0;
    this->node_idx_uid_map_.clear();
    this->edge_idx_uid_map_.clear();
    this->num_node_deletes_ = 0;
    this->num_edge_deletes_ = 0;
    this->num_node_adds_ = 0;
    this->num_edge_adds_ = 0;
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
      node_idx_ = 0;
      graph_len_ = 0;
    }

    /** Dereference operator for NodeIterator
     * Returns a Node
     */
    value_type operator*() const {
      const size_type node_uid = this->graph_
                                     ->node_idx_uid_map_[this->node_idx_];

      if (this->graph_ == nullptr || 
          this->node_idx_ >= this->graph_len_ ||
          !(this->graph_->nodes_[node_uid].valid_)) {
        std::cerr << "Error: Attempting to dereference to invalid node." 
                  << std::endl;
        exit(0);   
      }

      return Node(this->graph_, this->node_idx_);
    }

    /** Increment operator for NodeIterator
     * Returns a NodeIterator reference
     */
    NodeIterator& operator++() {
      if (this->node_idx_ == this->graph_len_) {
        std::cerr << "Error: Attempting to increment at end." << std::endl;
        exit(0);
      }

      this->node_idx_ += 1;
      // make sure node is valid
      while ((this->node_idx_ != this->graph_->num_nodes()) && 
              !(this->graph_->nodes_[
                            this->graph_->node_idx_uid_map_[this->node_idx_]]
                            .valid_)) {
        this->node_idx_ += 1;
      }      

      return (*this);
    }

    /** Is Equal operator for NodeIterator
     * Returns a boolean
     */
    bool operator==(const NodeIterator& node_iter) const {
      if (this->graph_ == node_iter.graph_
          && this->node_idx_ == node_iter.node_idx_) {
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
    size_type node_idx_;
    // Length of node id vector
    size_type graph_len_;

    // private constructor
    NodeIterator(const graph_type* graph, size_type node_idx) 
      : graph_(const_cast<graph_type*>(graph)), node_idx_(node_idx) {
      graph_len_ = this->graph_->size();

      // make sure node is valid
      while ((this->node_idx_ < this->graph_->num_nodes()) && 
              !(this->graph_->nodes_[
                this->graph_->node_idx_uid_map_[this->node_idx_]]
                .valid_)) {

        this->node_idx_ += 1;
      }
    }
  };

  /** Returns iterator at beginning of node vector */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  /** Returns iterator at end of node vector */
  node_iterator node_end() const {
    return NodeIterator(this, this->num_nodes_);
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

    /** Dereference operator for IncidentIterator
     * Returns an Edge
     */
    Edge operator*() const {
      if (this->graph_ == nullptr || 
          this->cent_node_ >= this->graph_->num_nodes() ||
          !(this->graph_->edges_[(*(this->map_iterator_)).second].valid_)) {
        std::cerr << "Error: Attempting to dereference to invalid edge."
                  << std::endl;
        exit(0);
      }
      const size_type edge_uid = (*(this->map_iterator_)).second;
      const size_type edge_idx = this->graph_->edges_[edge_uid].index_;
      const size_type node_uid = (*(this->map_iterator_)).first;
      const size_type node_idx = this->graph_->nodes_[node_uid].index_;

      return Edge(this->graph_, edge_idx, this->cent_node_, node_idx);
    }

    /** Increment operator for IncidentIterator
     * Returns an IncidentIterator reference
     */
    IncidentIterator& operator++() {

        ++(this->map_iterator_);
      
        // make sure node is valid
        const size_type node_uid = this->graph_->
                                   node_idx_uid_map_.at(this->cent_node_);
          while (this->map_iterator_ != 
                 this->graph_->node_id_map_[node_uid].end() &&
                 !(this->graph_->edges_
                   [(*(this->map_iterator_)).second].valid_)) {
          ++(this->map_iterator_);
          }
      
      return (*this);
    }

    /** Is Equal operator for IncidentIterator
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

    /** Not Equal operator for IncidentIterator
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
      
      // make sure nodes is valid
      const size_type node_uid = this->graph_->
                                 node_idx_uid_map_.at(this->cent_node_);
      while (this->map_iterator_ != 
             this->graph_->node_id_map_[node_uid].end() &&
             !(this->graph_->edges_[
              (*(this->map_iterator_)).second].valid_)) {     
        ++(this->map_iterator_);
      }
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
      edge_idx_ = 0;
      graph_len_ = 0;
    }

    /** Dereference operator for EdgeIterator
     * Returns an Edge
     */
    Edge operator*() const {
      const size_type edge_uid = this->graph_
                                     ->edge_idx_uid_map_[this->edge_idx_];
      if (this->graph_ == nullptr || 
          this->edge_idx_ >= this->graph_len_ ||
          !(this->graph_->edges_[edge_uid].valid_)) {
        std::cerr << "Error: Attempting to dereference to invalid edge."
                  << std::endl;
        exit(0);
      }

      const size_type node1_idx = this->graph_->edges_[edge_uid].node1_idx_;
      const size_type node2_idx = this->graph_->edges_[edge_uid].node2_idx_;

      return Edge(this->graph_, this->edge_idx_, node1_idx, node2_idx);
    }

    /** Increment operator for EdgeIterator
     * Returns an EdgeIterator reference
     */
    EdgeIterator& operator++() {
      if (this->edge_idx_ == this->graph_len_) {
        std::cerr << "Error: Attempting to increment at end." << std::endl;
        exit(0);
      }
      
      edge_idx_ += 1;
      // make sure node is valid
      while ((this->edge_idx_ != this->graph_->num_edges()) && 
              !(this->graph_->edges_[
                this->graph_->edge_idx_uid_map_[this->edge_idx_]]
                .valid_)) {

        this->edge_idx_ += 1;
      }

      return (*this);
    }

    /** Is Equal operator for EdgeIterator
     * Returns a boolean
     */
    bool operator==(const EdgeIterator& edge_iter) const {
      if (this->graph_ == edge_iter.graph_
          && this->edge_idx_ == edge_iter.edge_idx_) {
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
    size_type edge_idx_;
    // number of edges in graph
    size_type graph_len_;

    // private constructor
    EdgeIterator(const graph_type* graph, size_type edge_idx) 
      : graph_(const_cast<graph_type*>(graph)), edge_idx_(edge_idx) {
      graph_len_ = graph_->num_edges();

      // make sure node is valid
      while ((this->edge_idx_ != this->graph_->num_edges()) && 
              !(this->graph_->edges_[
                this->graph_->edge_idx_uid_map_[this->edge_idx_]]
                .valid_)) {

        this->edge_idx_ += 1;
      }
    }    
  };

  /** Returns iterator at beginning of edge map */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  /** Returns iterator at end of edge map */
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->num_edges_);
  }

  /** Remove an edge from the graph, if it exists in the graph.
   * @pre @a a and @a b are node objects
   * @return 1 if an edge is removed; 0 otherwise
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() - 1.
   * 
   * Data for this edge will still be present in graph containers, but will be
   * marked invalid and inaccessible.
   * 
   * This operation reindexes existing edges, so old edge(@a i) might not
   * equal new edge(@a i).
   * 
   * edge_iterators and incident_iterators will skip edges that 
   * have been removed.
   *
   * Complexity: O(num_edges()), due to reindexing existing edges
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // check if nodes exist
    if (!(this->has_node(a)) || !(this->has_node(b))) {
      return 0;
    }

    // check if edge exists
    if (!(this->has_edge(a, b))) {
      return 0;
    }

    // invalidate in edges_
    const size_type a_uid = a.get_uid();
    const size_type b_uid = b.get_uid();
    const size_type erase_uid = this->node_id_map_[a_uid][b_uid];

    this->edges_[erase_uid].valid_ = false;

    // update indices
    size_type erase_idx = this->edges_[erase_uid].index_;
    for (size_type i = erase_idx; i < this->num_edges() - 1; i++) {
      // update indices
      size_type new_uid = this->edge_idx_uid_map_[i + 1];
      this->edge_idx_uid_map_[i] = new_uid;
      this->edges_[this->edge_idx_uid_map_[i]].index_ = i;
    }

    // decrement degree of node
    this->nodes_[a_uid].degree_ -= 1;
    this->nodes_[b_uid].degree_ -= 1;

    // remove last entry in edge_idx_uid_map_
    this->edge_idx_uid_map_.erase(this->num_edges() - 1);

    // decrement num_edges_
    this->num_edges_ -= 1;

    // increment num_edge_deletes
    this->num_edge_deletes_ += 1;

    return 1;
  }

  /** Remove an edge from the graph, if it exists in the graph.
   * @pre @a e is an edge object
   * @return 1 if an edge is removed; 0 otherwise
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() - 1.
   * 
   * Data for this edge will still be present in graph containers, but will be
   * marked invalid and inaccessible.
   * 
   * This operation reindexes existing edges, so old edge(@a i) might not
   * equal new edge(@a i).
   * 
   * edge_iterators and incident_iterators will skip edges that 
   * have been removed.
   *
   * Complexity: O(num_edges()), due to reindexing existing edges
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from the graph, if it exists in the graph.
   * @pre @a e_it is an edge_iterator object
   * @return 1 if an edge is removed; 0 otherwise
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() - 1.
   * @post @a e_it is incremented to the next valid edge
   * 
   * Data for this edge will still be present in graph containers, but will be
   * marked invalid and inaccessible.
   * 
   * This operation reindexes existing edges, so old edge(@a i) might not
   * equal new edge(@a i).
   * 
   * edge_iterators and incident_iterators will skip edges that 
   * have been removed.
   *
   * Complexity: O(num_edges()), due to reindexing existing edges
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    ++e_it;
    return e_it;
  }

  /** Remove a node and its indident edges from the graph, if it exists 
   * in the graph.
   * @pre @a n is a node object
   * @return 1 if a node is removed; 0 otherwise
   * @post has_node(@a a) == false
   * @post If old has_node(@a n), new num_nodes() == old num_nodes().
   *                              new size()      == old size().
   *       Else,                  new num_edges() == old num_edges() - 1.
   *                              new size ()     == old size() - 1.
   * @post For all edges e_i incident to @a n, has_edge(e_i) == false
   * @post new num_edges() = old num_edges() - num edges incident to @a n
   * 
   * Data for this node and its incident edges will still be present in 
   * graph containers, but will be marked invalid and inaccessible.
   * 
   * This operation reindexes existing nodes, so old node(@a i) might not
   * equal new node(@a i).
   * 
   * This operation reindexes existing edges, so old edge(@a i) might not
   * equal new edge(@a i).
   * 
   * node_iterators will skip nodes that have been removed.
   * 
   * edge_iterators and incident_iterators will skip edges that 
   * have been removed.
   *
   * Complexity: O(num_nodes()), due to reindexing existing of nodes. This
   * assumes sparsity of the graph, i.e. the maximum degree of a node is 
   * considered constant relative to the number of nodes in the graph.
   */
  size_type remove_node(const Node& n) {
    // check if node exists in graph
    if (!(this->has_node(n))) {
      return 0;
    }

    // remove adjacent edges
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      remove_edge(*it);
    }

    // invalidate in nodes_
    this->nodes_[n.get_uid()].valid_ = false;

    // update indices
    size_type erase_idx = n.index();
    for (size_type i = erase_idx; i < this->num_nodes_ - 1; i++) {
      size_type new_uid = this->node_idx_uid_map_[i + 1];
      this->node_idx_uid_map_[i] = new_uid;
      this->nodes_[this->node_idx_uid_map_[i]].index_ = i;

      // update indices in adjacent edges' EdgeData
      for (auto it = this->node(i).edge_begin(); 
                it != this->node(i).edge_end(); ++it) {
        size_type adj_uid = this->edge_idx_uid_map_[(*it).edge_idx_];
        if (this->edges_[adj_uid].node1_idx_ == i + 1) {
          this->edges_[adj_uid].node1_idx_ = i;
        }
        if (this->edges_[adj_uid].node2_idx_ == i + 1) {
          this->edges_[adj_uid].node2_idx_ = i;
        }
      }
    }

    // remove last entry in edge_idx_uid_map_
    this->node_idx_uid_map_.erase(this->num_nodes_ - 1);

    // decrement num_nodes_
    this->num_nodes_ -= 1;

    // increment num_node_deletes_
    this->num_node_deletes_ += 1;

    return 1;
  }

  /** Remove a node and its indident edges from the graph, if it exists 
   * in the graph.
   * @pre @a n_it is a node_iterator object
   * @return 1 if a node is removed; 0 otherwise
   * @post has_node(@a a) == false
   * @post If old has_node(@a n), new num_nodes() == old num_nodes().
   *                              new size()      == old size().
   *       Else,                  new num_edges() == old num_edges() - 1.
   *                              new size ()     == old size() - 1.
   * @post For all edges e_i incident to @a n, has_edge(e_i) == false
   * @post new num_edges() = old num_edges() - num edges incident to @a n
   * @post @a n_it is incremented to the next valid node
   * 
   * Data for this node and its incident edges will still be present in 
   * graph containers, but will be marked invalid and inaccessible.
   * 
   * This operation reindexes existing nodes, so old node(@a i) might not
   * equal new node(@a i).
   * 
   * This operation reindexes existing edges, so old edge(@a i) might not
   * equal new edge(@a i).
   * 
   * node_iterators will skip nodes that have been removed.
   * 
   * edge_iterators and incident_iterators will skip edges that 
   * have been removed.
   *
   * Complexity: O(num_nodes()), due to reindexing existing of nodes. This
   * assumes sparsity of the graph, i.e. the maximum degree of a node is 
   * considered constant relative to the number of nodes in the graph.
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    ++n_it;
    return n_it;
  }

 private:
  // structure to store node data in graph containers
  struct NodeData {
    Point position_;
    node_value_type value_;
    size_type index_;
    size_type uid_;
    bool valid_;
    size_type degree_;

    NodeData(Point position, node_value_type value, size_type index,
             size_type uid, bool valid, size_type degree)
      : position_(position), value_(value), index_(index),
        uid_(uid), valid_(valid), degree_(degree) {
    }

    NodeData()
      : position_(Point()), value_(node_value_type()), index_(0),
        uid_(0), valid_(false), degree_(0) {
    }
  };

  // structure to store edge data in graph containers
  struct EdgeData {
    size_type node1_idx_;
    size_type node2_idx_;
    edge_value_type value_;
    size_type index_;
    size_type uid_;
    bool valid_;

    EdgeData(size_type node1_idx, size_type node2_idx, edge_value_type value,
             size_type index, size_type uid, bool valid)
      : node1_idx_(node1_idx), node2_idx_(node2_idx), value_(value),
        index_(index), uid_(uid), valid_(valid) {
    }

    EdgeData()
      : node1_idx_(0), node2_idx_(0), value_(edge_value_type()),
        index_(0), uid_(0), valid_(false) {
    }
  };

  // stores number of nodes
  size_type num_nodes_;

  // stores number of edges
  size_type num_edges_;

  // number of edge deletions
  size_type num_node_deletes_;

  // number of edge deletions
  size_type num_edge_deletes_;

  // number of node additions
  size_type num_node_adds_;

  // number of edge additions
  size_type num_edge_adds_;

  // map node uid to NodeData struct
  std::unordered_map<size_type, NodeData> nodes_;
  
  // map edge uid to EdgeData struct
  std::unordered_map<size_type, EdgeData> edges_;

  // edge map
  // nested map of node uid -> connected node uids -> edge uid
  std::unordered_map<size_type, 
                     std::unordered_map<size_type, 
                                        size_type>> node_id_map_;

  // map node index to unique id
  std::unordered_map<size_type, size_type> node_idx_uid_map_;

  // map edge index to unique id
  std::unordered_map<size_type, size_type> edge_idx_uid_map_;
};

#endif // CME212_GRAPH_HPP
