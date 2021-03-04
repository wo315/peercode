#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
//--functionality_0
//--missing header <map>
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

  struct internal_node;
  struct internal_edge;

 private:

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
  using node_value_type = V;
  using edge_value_type = E;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
	: nodes_(), node_i2u_(), num_nodes_(0), edges_(), edge_i2u_(), num_edges_(0), adjacency(){
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph(){
      nodes_.clear();
      edges_.clear();
      node_i2u_.clear();
      edge_i2u_.clear();
      adjacency.clear();
  }

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered<Node>{
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
      return graph_->nodes_[graph_->node_i2u_[idx_]].p;
    }

    /** Return this node's position, which is not a const variable. */
    Point& position(){
        return graph_->nodes_[graph_->node_i2u_[idx_]].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->idx_;
    }

    /**
     * @brief Get the reference to the value of the node.
     *
     * @return The reference to the value of the node.
     *
     * This function can also be used as both setter
     * and getter of the value since reference is
     * returned.
     */
    node_value_type& value(){
        return graph_->nodes_[graph_->node_i2u_[idx_]].val;
    }

    /**
     * @brief Get the reference to the value of the node.
     *
     * @return The reference to the value of the node.
     *
     * Since the value here is const, this function can only
     * be used as the getter of the value, not the setter of
     * the value.
     */
    const node_value_type& value() const{
        return graph_->nodes_[graph_->node_i2u_[idx_]].val;
    }

    /**
     * @brief Get the degree of the node.
     *
     * @return The degree of the node.
     *
     * This function uses the adjacency map, which only contains
     * nodes that are not isolated as keys. Hence will directly
     * return 0 if a node's index is not key in the map.
     */
    size_type degree() const{
        auto curr = graph_->adjacency.find(graph_->node_i2u_[idx_]);
        if(curr == graph_->adjacency.end()){
            return 0;
        }else{
            return curr -> second.size();
        }
    }

    /**
     * @brief Returns an Incident Iterator at first edge in the edge list of the node.
     *
     * @return incident_iterator pointing at the first edge.
     */
    incident_iterator edge_begin(){
        return IncidentIterator(graph_, this, 0);
    }

    /**
     * @brief Returns an Incident Iterator at last edge in the edge list of the node.
     *
     * @return incident_iterator pointing at the last edge.
     */
    incident_iterator edge_end(){
        return IncidentIterator(graph_, this, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if(this->graph_ != n.graph_){
          return false;
      }
      if (this->graph_ == n.graph_ && this->index() == n.index()){
          return true;
      }else{
          return false;
      }
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
      if (this->graph_ != n.graph_){
          return true;
      }
      if(this->index() < n.index()){
            return true;
        }else{
            return false;
        }
    }

   private:
     // Pointer back to the Graph container
     Graph* graph_;
     // This node's index
     size_type idx_;
     /** Private Constructor */
     Node(const Graph* graph, size_type idx)
         : graph_(const_cast<Graph*>(graph)), idx_(idx) {
     }
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
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
    // HW0: YOUR CODE HERE
    return this->num_nodes_;
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
  Node add_node(const Point& position, const node_value_type& nodeval = node_value_type()) {
    // HW0: YOUR CODE HERE
    // HW1: modified hw0 code
    ++num_nodes_;
    internal_node i_node;
    i_node.p = position;
    i_node.uid = this->nodes_.size();
    i_node.val = nodeval;
    i_node.idx = this->node_i2u_.size();
    Node node = Node(this, i_node.idx);
    nodes_.push_back(i_node);
    node_i2u_.push_back(i_node.uid);
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.index() < size()){
          return true;
      }else{
          return false;
      }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_nodes_);
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
  class Edge: private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge(){}

    /** Return a node of this Edge */
    Node node1() const{
      // HW0: YOUR CODE HERE
      return graph_->edges_[graph_->edge_i2u_[idx_]].node1;
    }

    /** Return the other node of this Edge */
    Node node2() const{
      // HW0: YOUR CODE HERE
      return graph_->edges_[graph_->edge_i2u_[idx_]].node2;
    }

    size_type index() const{
      // HW0: YOUR CODE HERE
      return idx_;
    }

    /**
     * @brief Get the reference to the value of the edge.
     *
     * @return The reference to the value of the edge.
     *
     * This function can also be used as both setter
     * and getter of the value since reference is
     * returned.
     */
    edge_value_type& value(){
          return graph_->edges_[graph_->edge_i2u_[idx_]].val;
    }

    /**
     * @brief Get the reference to the value of the edge.
     *
     * @return The reference to the value of the edge.
     *
     * Since the value here is const, this function can only
     * be used as the getter of the value, not the setter of
     * the value.
     */
    const edge_value_type& value() const{
        return graph_->edges_[graph_->edge_i2u_[idx_]].val;
    }

    /**
     * @brief Get the length of the edge.
     *
     * @return The L2 norm of the edge.
     */
    double length() const{
        Point p1 = node1().position();
        Point p2 = node2().position();
        return norm(p1 - p2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const{
        if(this->graph_ != e.graph_){
            return false;
        }
        if(idx_ == e.index()){
            return true;
        }else{
	        return false;
        }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const{
        if(this->graph_ != e.graph_){
            return true;
        }
        if(idx_ < e.index()){
	        return true;
        }else{
            return false;
        }
    }

   private:
     // Pointer back to the Graph container
     Graph* graph_;
     // This edge's index
     size_type idx_;
     /** Private Constructor */
     Edge(const Graph* graph, size_type idx)
         : graph_(const_cast<Graph*>(graph)), idx_(idx){
     }
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
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
    assert(i < num_edges_);
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const{
    // HW0: YOUR CODE HERE
    size_type a_uid = this->node_i2u_[a.index()];
    size_type b_uid = this->node_i2u_[b.index()];
    auto node_a = this->adjacency.find(a_uid);
    if(node_a == this->adjacency.end()){
        return false;
    }
    std::map<size_type, size_type> a_adj = node_a -> second; // map of adjacent node id : edge id
    auto node_b = a_adj.find(b_uid);
    if (node_b == a_adj.end()){
        return false;
    }
    return true;
//    for (size_type i = 0; i < edge_i2u_.size(); ++i){
//        size_type edge_uid = edge_i2u_[i];
//        internal_edge e = edges_[edge_uid];
//        if ((e.node1 == a && e.node2 == b) || (e.node2 == a && e.node1 == b)){
//            return true;
//        }
//    }
//    return false;
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
  Edge add_edge(const Node& a, const Node& b){
    // HW0: YOUR CODE HERE
    if(has_edge(a,b)){
        size_type i = adjacency[node_i2u_[a.index()]][node_i2u_[b.index()]];
        edges_[i].node1 = a;
        edges_[i].node2 = b;
        size_type edge_idx = edges_[i].idx;
        return Edge(this, edge_idx);
    }
    internal_edge i_edge;
    i_edge.node1 = a;
    i_edge.node2 = b;
    i_edge.uid = this->edges_.size();
    i_edge.idx = this->edge_i2u_.size();
    this->edges_.push_back(i_edge);
    this->edge_i2u_.push_back(i_edge.uid);
    this->adjacency[this->node_i2u_[a.index()]][this->node_i2u_[b.index()]] = i_edge.uid;
    this->adjacency[this->node_i2u_[b.index()]][this->node_i2u_[a.index()]] = i_edge.uid;
    this->num_edges_ = edge_i2u_.size();
    Edge edge = Edge(this, i_edge.idx);
    return edge;
  }

  /** Remove edge connected by n1 and n2 if exists
   *
   * @param n1[in] First node of the edge
   * @param n2[in] Second node of the edge
   * @return 1 if an edge is removed, 0 if no edge is removed
   */
  size_type remove_edge(const Node& n1, const Node& n2){
      if(has_edge(n1, n2)){
          size_type n1_uid = this->node_i2u_[n1.index()];
          size_type n2_uid = this->node_i2u_[n2.index()];
          size_type e_uid = this->adjacency[n1_uid][n2_uid];
          size_type e_idx = this->edges_[e_uid].idx;
          // remove edge from edge_i2u
          auto last_edge_uid = edge_i2u_[edge_i2u_.size() - 1];
          edge_i2u_[e_idx] = last_edge_uid;
          edge_i2u_.pop_back();
          edges_[last_edge_uid].idx = e_idx;
          // update adjacency
          this->adjacency.clear();
          for (size_type i = 0; i < this->edge_i2u_.size(); ++i){
              size_type edge_uid = edge_i2u_[i];
              Node node1 = edges_[edge_uid].node1;
              Node node2 = edges_[edge_uid].node2;
              size_type node1_uid = this->node_i2u_[node1.index()];
              size_type node2_uid = this->node_i2u_[node2.index()];
              this->adjacency[node1_uid][node2_uid] = edge_uid;
              this->adjacency[node2_uid][node1_uid] = edge_uid;
          }
          this->num_edges_ = this->edge_i2u_.size();
          return 1;
      }
      return 0;
  }

  /** Remove the input edge from graph if exists
   *
   * @param e[in] Reference to the target edge to be removed
   * @return 1 if an edge is removed, 0 if no edge is removed
   */
  size_type remove_edge(const Edge& e) {
      Node n1 = e.node1();
      Node n2 = e.node2();
      return remove_edge(n1, n2);
  }

  /** Remove the edge pointed by the edge iterator
   *
   * @param e_it edge_iterator pointing at the edge to be removed
   * @return edge_iterator pointing another edge
   */
  edge_iterator remove_edge(edge_iterator e_it){
      auto e = *e_it;
      this->remove_edge(e);
      return this->edge_begin();
  }

  /** Remove node n if exists
   *
   * @param n Reference to the node
   * @return 1 if a node is removed else 0
   */
  size_type remove_node(const Node& n){
      if (has_node(n)){
          // remove edges of the node
          size_type n_uid = this->node_i2u_[n.index()];
          while(this->adjacency.find(n_uid) != this->adjacency.end()){
              auto n_a = this->adjacency.find(n_uid);
              std::map<size_type, size_type> a_adj = n_a -> second;
              size_type edge_uid = a_adj.begin()->second;
              size_type edge_idx = edges_[edge_uid].idx;
              this->remove_edge(edge(edge_idx));
          }
          // decrement the number of nodes
          --this->num_nodes_;
          // remove node from node_i2u_
          auto last_node_uid = node_i2u_[node_i2u_.size() - 1];
          node_i2u_[n.index()] = last_node_uid;
          node_i2u_.pop_back();
          nodes_[last_node_uid].idx = n.index();
          return 1;
      }
      return 0;
  }

  /** Remove node pointed by the node_iterator
    *
    * @param n_it node_iterator point at the node to be removed
    * @return node_iterator pointing at another node
    */
  //--functionality_1
  //--node_begin may already been removed
  //--START
  node_iterator remove_node(node_iterator n_it){
      Node n = *n_it;
      remove_node(n);
      return this->node_begin();
  }
  //--END

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear(){
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
    num_nodes_ = 0;
    num_edges_ = 0;
    adjacency.clear();
    node_i2u_.clear();
    edge_i2u_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable<NodeIterator>{
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
     * @brief Dereference the current node.
     *
     * @return A Node object at the position of the NodeIterator.
     */
    Node operator*() const{
        return Node(graph_, node_idx_);
    }

    /**
     * @brief The NodeIterator goes to the next node.
     *
     * @return Reference to the NodeIterator at the next position of the current node.
     */
    NodeIterator& operator++(){
        node_idx_ ++;
        return *this;
    }

    /**
     * @brief Compares if two NodeIterators are the same.
     *
     * @param n_i[in] Reference to the other NodeIterator.
     * @return A boolean indicating if this and the other NodeIterator are the same one.
     */
    bool operator==(const NodeIterator& n_i) const{
        return (n_i.graph_ == graph_) && (n_i.node_idx_ == node_idx_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const graph_type* graph_;
    size_type node_idx_;

    /**
     * @brief Constructor of the NodeIterator.
     *
     * @param graph[in] A pointer to the const graph_type.
     * @param node_id[in] Id of the node to be pointed at.
     *
     * Initialize graph_ to graph, node_id_ to node_id.
     */
    NodeIterator(const graph_type* graph, size_type node_idx){
        graph_ = graph;
        node_idx_ = node_idx;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
   * @brief A node iterator to the first node in the graph.
   *
   * @return node_iterator pointing to the first node in the graph.
   */
  node_iterator node_begin() const{
      return NodeIterator(this, 0);
  }

  /**
   * @brief A node iterator to the last node in the graph.
   *
   * @return node_iterator pointing to the last node in the graph.
   */
  node_iterator node_end() const{
      return NodeIterator(this, this->num_nodes_);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private equality_comparable<IncidentIterator> {
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
     * @brief Dereference the current edge.
     *
     * @return The edge object that's pointed by the iterator.
     */
    Edge operator*() const{
        size_type node_idx_ = this->node_->index();
        std::map<size_type, size_type> node_adj = graph_->adjacency.find(graph_->node_i2u_[node_idx_])->second;
        auto it = node_adj.begin();
        std::advance(it, incident_id_);
        size_type edge_uid_ = it -> second;
        size_type neighbor_node_uid = it -> first;
        size_type neighbor_node_idx = graph_->nodes_[neighbor_node_uid].idx;
        this->graph_->edges_[edge_uid_].node1 = this->graph_->node(node_idx_);
        this->graph_->edges_[edge_uid_].node2 = this->graph_->node(neighbor_node_idx);
        return Edge(graph_, graph_->edges_[edge_uid_].idx);
    }

    /**
     * @brief The iterator goes to the next edge.
     *
     * @return Reference to the IncidentIterator pointing at the next edge.
     */
    IncidentIterator& operator++(){
        incident_id_++;
        return *this;
    }

    /**
     * @brief Check if the two iterators are the same.
     *
     * @param i_i[in] Reference to the other const IncidentIterator.
     * @return An boolean indicating if this and the other IncidentIterator are pointing at the same edge.
     */
    bool operator==(const IncidentIterator& i_i) const{
        return (i_i.graph_ == graph_) && (i_i.node_ == node_) && (i_i.incident_id_ == incident_id_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    node_type* node_;
    size_type incident_id_;

    /**
     * @brief Constructor of the IncidentIterator.
     * @param graph[in] A pointer to the graph.
     * @param node[in] A pointer to the node.
     * @param incident_id_[in] The number of incident edge.
     */
    IncidentIterator(graph_type* graph, node_type* node, size_type incident_id) {
        graph_ = graph;
        node_ = node;
        incident_id_ = incident_id;
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private equality_comparable<EdgeIterator> {
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
     * @brief Dereference the edge.
     *
     * @return The Edge pointed by the iterator.
     */
    Edge operator*() const{
        return Edge(graph_, edge_idx_);
    }

    /**
     * @brief Iterator goes to the next edge.
     *
     * @return Reference to the iterator pointing at the next edge.
     */
    EdgeIterator& operator++(){
        edge_idx_ ++;
        return *this;
    }

    /**
     * @brief Compare if two EdgeIterators are the same.
     *
     * @param e_i[in] Reference to the other const EdgeIterator
     * @return A boolean indicating if this and the other EdgeIterator are the same.
     */
    bool operator==(const EdgeIterator& e_i) const{
        return (e_i.graph_ == graph_) && (e_i.edge_idx_ == edge_idx_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const graph_type* graph_;
    size_type edge_idx_;

    /**
     * @brief Constructor of the EdgeIterator.
     *
     * @param graph[in] Pointer to the const graph.
     * @param edge_id[in] Edge id.
     */
    EdgeIterator(const graph_type* graph, size_type edge_id){
          graph_ = graph;
          edge_idx_ = edge_id;
      }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
   * @brief Get the edge_iterator pointing to the first edge.
   *
   * @return An edge iterator pointing to the first edge.
   */
  edge_iterator edge_begin() const{
      return EdgeIterator(this, 0);
  }

  /**
   * @brief Get the edge_iterator pointing to the last edge.
   *
   * @return An edge iterator pointing to the last edge.
   */
  edge_iterator edge_end() const{
      return EdgeIterator(this, this->num_edges_);
  }


 private:
   struct internal_node {
        Point p;
        size_type uid; // unique id
        node_value_type val;
        size_type idx; // only nodes not removed have this index
      };
    std::vector<internal_node> nodes_;
    std::vector<size_type> node_i2u_; // idx to uid
    size_type num_nodes_;

   struct internal_edge {
        Node node1;
        Node node2;
	    size_type uid;
	    edge_value_type val;
	    size_type idx;
      };
    std::vector<internal_edge> edges_;
    std::vector<size_type> edge_i2u_;
    size_type num_edges_;
    std::map<size_type, std::map<size_type, size_type>> adjacency; // key: node uid, value: map <node uid, edge id>
};

#endif // CME212_GRAPH_HPP
