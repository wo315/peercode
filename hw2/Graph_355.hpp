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
template <typename V, typename E>
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

  /** Type of pairs of indexes. */
  using pair = std::pair<size_type, size_type>;


  /** Types of node values and edge values  */
  typedef V node_value_type;
  typedef E edge_value_type;

  private:

    // Used to store information about nodes and edges
    struct nodeinfo;
    struct edgeinfo;

    // Hash function for pairs of indexes
    struct hash_pair;

    // nodes_ is a vector representing the nodes of the graph
    std::vector<nodeinfo> nodes_;

    // edge_vec_ is a vector representing the edges of the graph
    std::vector<edgeinfo> edge_vec_;

    // adj_ is a vector such that adj_[i] contains the indexes of the
    // neighbors of the node of index i. It is basically an adjacency list
    std::vector<std::vector<size_type>> adj_;

    // edge_reverse_ maps ordered pairs (i1,i2) (i1<i2) of indexes that compose
    // an edge to the corresponding edge internal index
    std::unordered_map<pair, size_type, hash_pair> edge_reverse_;

    // node_i2u_ is a vector that maps node user_facing indexes (idx_) to
    // internal indexes.
    std::vector<size_type> node_i2u_;

    // edge_i2u_ is a vector that maps edge user_facing indexes (idx_) to
    // internal indexes.
    std::vector<size_type> edge_i2u_;


  public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(): nodes_(std::vector<nodeinfo>()),
           edge_vec_(std::vector<edgeinfo>()),
           adj_(std::vector<std::vector<size_type>>()),
           edge_reverse_(std::unordered_map<pair, size_type, hash_pair>()),
           node_i2u_(std::vector<size_type>()),
           edge_i2u_(std::vector<size_type>())
         {
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
  class Node: private totally_ordered<Node> {
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
      index_ = 0;
    }

    /** Return this node's position. non-const to make node modifiable */
    Point& position(){
      // Look into node_vec_ with the appropriate index
      return (*graph_).nodes_[index_].position;
    }

    /** Return this node's position.*/
    const Point& position() const {
      return (*graph_).nodes_[index_].position;
    }

    /** Return this node's (user facing) index, a number in the range
     * [0, graph_size). */
    size_type index() const {
      return (*graph_).nodes_[index_].idx_;
    }

    /** Return this node's value */
    node_value_type& value(){
      return (*graph_).nodes_[index_].value;
    }

    /** Return this node's value */
    const node_value_type& value() const{
      return (*graph_).nodes_[index_].value;
    }

    /** Return this node's degree (number of incident edges)
     * @post degree() = size(adj_[index_])
     *
     * Complexity: O(1)
    */
    size_type degree() const{
      return (*graph_).adj_[index_].size();
    }

    /** Return an incident iterator at the first incident edge at this node */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, index_, 0);
    }

    /** Return an incident iterator after the last incident edge at this node */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, index_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ and index_ == n.index_);
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
      // As stated in @75 in Piazza, we can assume the nodes to be from
      // the same graph
      return index_ < n.index_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Define a node with a pointer to its graph and its index
    graph_type* graph_;
    size_type index_;

    /** Valid node constructor */
    Node(const graph_type* graph, size_type index) :
         graph_(const_cast<graph_type*>(graph)), index_(index) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_value The new node's value (optional)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.value() == node_value
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node (const Point & position,
                 const node_value_type & node_value = node_value_type()){



    // Add the nodeinfo to the nodes_ vector (user facing idx is the size)
    nodes_.push_back(nodeinfo(position, node_value, size()));

    // Add element to the mapping idx -> uid
    node_i2u_.push_back(nodes_.size() - 1);

    // Add an empty vector to the adjacency list
    adj_.push_back(std::vector<size_type>());
    return Node(this, size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_ and n.index() < size()) ;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,node_i2u_[i]);
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      graph_ = nullptr;
      index1_ = 0;
      index2_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, index1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, index2_);
    }

    /** Return this node's value */
    edge_value_type& value(){

      // Here, edge_reverse_ contains the key (min(indexes), max(indexes))
      size_type id_min = std::min(index1_, index2_);
      size_type id_max = std::max(index1_, index2_);

      size_type id = (*graph_).edge_reverse_.at(std::make_pair(id_min,
                                                               id_max));

      // Retrieve value using edge id
      return (*graph_).edge_vec_[id].edge_value;
    }

    /** Return this edge's value */
    const edge_value_type& value() const{

      // Here, edge_reverse_ contains the key (min(indexes), max(indexes))
      size_type id_min = std::min(index1_, index2_);
      size_type id_max = std::max(index1_, index2_);

      size_type id = (*graph_).edge_reverse_.at(std::make_pair(id_min,
                                                               id_max));

      // Retrieve value using edge id
      return (*graph_).edge_vec_[id].edge_value;

    }

    /** Return this edge's length (Euclidean distance between nodes) */
    double length() const {
      return norm(node1().position() - node2().position());
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Equality <=> same nodes IN ANY ORDER
      return (graph_ == e.graph_ and
             ((index1_ == e.index1_ and index2_ == e.index2_)
             or (index1_ == e.index2_ and index2_ == e.index1_)));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Lexicographic order between edges
      // represented as (graph, smallest index, largest index)
      size_type this_min = std::min(index1_, index2_);
      size_type this_max = std::max(index1_, index2_);
      size_type e_min = std::min(e.index1_, e.index2_);
      size_type e_max = std::max(e.index1_, e.index2_);

      return ((graph_ < e.graph_) or
             ((graph_ == e.graph_) and (this_min < e_min)) or
             ((graph_ == e.graph_) and (this_min == e_min)
                                   and (this_max < e_max)));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // An Edge is represented with a pointer to its graph and the two indexes
    // of the nodes it connects

    graph_type* graph_;
    size_type index1_;
    size_type index2_;



    // Edge constructor
    Edge(const graph_type* graph, size_type index1, size_type index2):
         graph_(const_cast<graph_type*>(graph)),
         index1_(index1), index2_(index2){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return edge_i2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    // It can easily be found in the vector edge_vec_

    //First map user facing id to internal id
    size_type internal_id = edge_i2u_[i];
    return Edge(this, edge_vec_[internal_id].index1,
                      edge_vec_[internal_id].index2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(a.degree())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // It is sufficent to look into the neighbors of a

    for (auto it = adj_[a.index_].begin(); it != adj_[a.index_].end(); ++it){
      if ((*it) == b.index_) return true;
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

    if (!has_edge(a, b)){
      // Add index1/index2/default value/idx  to edge_vec_
      edge_vec_.push_back(edgeinfo(a.index_, b.index_,
                                   edge_value_type(), num_edges()));

      // Map ordered pair of node indexes  to internal edge index
      size_type id_min = std::min(a.index_, b.index_);
      size_type id_max = std::max(a.index_, b.index_);
      edge_reverse_[std::make_pair(id_min,id_max)] = edge_vec_.size()- 1;

      // Map user facing edge idx to internal edge index
      edge_i2u_.push_back(edge_vec_.size()- 1);

      // We add a to the neighbors of b and b to the neighbors of a
      adj_[a.index_].push_back(b.index_);
      adj_[b.index_].push_back(a.index_);
    }
    // e.node1() == a and e.node2() == b is verified
    return Edge(this, a.index_, b.index_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // We use the clear method on all vectors
    nodes_.clear();
    adj_.clear();
    edge_vec_.clear();
    edge_reverse_.clear();
    node_i2u_.clear();
    edge_i2u_.clear();
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
      idx_ = 0;
    }

    /** Dereference operator */
    Node operator*() const {
      // Map user facing idx to internal index
      return Node(graph_, (*graph_).node_i2u_[idx_]);
    }

    /** Increments to the next (valid) Node in the graph */
    NodeIterator& operator++() {
      idx_++;
      return *this;
    }

    /** Defines equality between iterators
     *  @param[in] ni Node iterator to which we compare this node iterator
    */
    bool operator==(const NodeIterator& ni) const{
      return (graph_ == ni.graph_ and idx_ == ni.idx_);
    }

    /** Defines inequality between iterators
     *  @param[in] ni Node iterator to which we compare this node iterator
    */
    bool operator!=(const NodeIterator& ni) const{
      return (!(*this == ni));
    }


   private:
    friend class Graph;

    // Pointer to the graph on which we iterate
    graph_type* graph_;

    // User facing idx of current node
    size_type idx_;

    /** Valid NodeIterator constructor */
    NodeIterator(const graph_type* graph, size_type idx) :
         graph_(const_cast<graph_type*>(graph)), idx_(idx) {
    }

  };

  /** Return a node iterator at the first node of this graph */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** Return a node iterator after the last node of this graph */
  node_iterator node_end() const{
    return NodeIterator(this, this->node_i2u_.size());
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
      graph_ = nullptr;
      main_index_ = 0;
      neighbor_index_  = 0;
    }

    /** Dereference operator */
    Edge operator*() const {
      return Edge(graph_, main_index_,
                  (*graph_).adj_[main_index_][neighbor_index_]);
    }

    /** Increments to the next edge incident to this node */
    IncidentIterator& operator++() {
      neighbor_index_++;
      return *this;
    }

    /** Defines equality between iterators
     *  @param[in] ii Incident iterator to which we compare this node iterator
    */
    bool operator==(const IncidentIterator& ii) const {
      return (graph_ == ii.graph_
              and main_index_ == ii.main_index_
              and neighbor_index_ == ii.neighbor_index_);
    }

    /** Defines inequality between iterators
     *  @param[in] ii Incident iterator to which we compare this node iterator
    */
    bool operator!=(const IncidentIterator& ii) const {
      return (!(*this == ii));
    }


   private:
    friend class Graph;

    // Pointer to the graph
    graph_type* graph_;

    // Index of the main node (the one of which we consider the neighbors)
    size_type main_index_;

    // Index in the adjacency list (!= node.index_) of the current neighbor
    // of the node of internal index main_index_
    size_type neighbor_index_;

    /** Valid IncidentIterator constructor */
    IncidentIterator(const graph_type* graph, size_type i1, size_type i2) :
         graph_(const_cast<graph_type*>(graph)),
         main_index_(i1), neighbor_index_(i2){
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
      idx_ = 0;
    }


    /** Dereference operator */
    Edge operator*() const {
      // Map user facing idx to internal index using edge_i2u_
      return Edge(graph_,
                  (*graph_).edge_vec_[(*graph_).edge_i2u_[idx_]].index1,
                  (*graph_).edge_vec_[(*graph_).edge_i2u_[idx_]].index2);
    }

    /** Increments to the next edge in the graph */
    EdgeIterator& operator++() {
      idx_++;
      return *this;
    }

    /** Defines equality between iterators
     *  @param[in] ei Edge iterator to which we compare this node iterator
    */
    bool operator==(const EdgeIterator& ei) const{
      return (graph_ == ei.graph_ and idx_ == ei.idx_);
    }

    /** Defines inequality between iterators
     *  @param[in] ei Edge iterator to which we compare this node iterator
    */
    bool operator!=(const EdgeIterator& ei) const{
      return (!(*this == ei));
    }

   private:
    friend class Graph;
    // Note that we could have combined NodeIterator and IncidentIterator to
    // construct EdgeIterator. Instead, we take advantage of edge_vec_ which
    // allows a simple implementation of EdgeIterator.

    // Pointer to the graph
    graph_type* graph_;

    // User facing index of current edge
    size_type idx_;

    /** Valid EdgeIterator constructor */
    EdgeIterator(const graph_type* graph, size_type idx) :
         graph_(const_cast<graph_type*>(graph)), idx_(idx) {
    }

  };

  /** Return an edge iterator at the first edge of this graph */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Return an edge iterator after the last (valid) edge of this graph */
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

  /** Remove edge between node a and b
   * @param[in] a A valid node of this graph
   * @param[in] b A valid node of this graph
   *
   * @post if old has_edge(a,b) then new num_edges() = old num_edges() - 1
   * @post if !old has_edge(a,b) then new num_edges() = old num_edges()
   * @post new has_edge(a,b) == false
   *
   * @post Invalidates:  Edge between a and b (if old has_edge(a,b))
   * @post EdgeIterator or IncidentIterator pointing at this edge will point at
   *       next Edge
   *
   * @return 1 if old has_edge(a,b), 0 else
   *
   * Complexity: O(a.degree() + b.degree()) amortized operations.
   */
  size_type remove_edge (const Node& a , const Node& b) {
    if (!has_edge(a,b)){
      return 0;
    }
    else {

      // https://stackoverflow.com/questions/3385229/
      // Remove b from the neighbors of a
      adj_[a.index_].erase(std::remove(adj_[a.index_].begin(),
                                       adj_[a.index_].end(),
                                       b.index_),
                           adj_[a.index_].end());

      // Remove a from the neighbors of b
      adj_[b.index_].erase(std::remove(adj_[b.index_].begin(),
                                       adj_[b.index_].end(),
                                       a.index_),
                           adj_[b.index_].end());


      // Find internal index of current edge

      size_type id_min = std::min(a.index_, b.index_);
      size_type id_max = std::max(a.index_, b.index_);
      size_type id = edge_reverse_.at(std::make_pair(id_min,id_max));

      // User facing index
      size_type idx = edge_vec_[id].idx_;

      // Overwrite uid at idx position with last uid
      edge_i2u_[idx] = edge_i2u_.back();

      // Pop last element
      edge_i2u_.pop_back();

      // Update idx of edge that was copied in place of removed edge
      edge_vec_[edge_i2u_[idx]].idx_ = idx;

      return 1;
    }
  }

  /** Remove edge e
   * @param[in] e An edge of this graph
   *
   * @post if e is valid then new num_edges() = old num_edges() - 1
   * @post if e is not valid then new num_edges() = old num_edges()
   *
   * @post Invalidates:  Edge e (if was valid before)
   * @post EdgeIterator or IncidentIterator pointing at this edge will point at
   *       next Edge
   *
   * @return 1 if e was valid, 0 else
   *
   * Complexity: O(e.node1().degree() + e.node2().degree()) amortized operations
   */
  size_type remove_edge (const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove edge e such that e = *e_it
   * @param[in] e_it An edge iterator
   *
   * @post if *e_it is valid then new num_edges() = old num_edges() - 1
   * @post if *e_it is not valid then new num_edges() = old num_edges()
   *
   * @post Invalidates:  Edge *e_it (if was valid before)
   * @post e_it, if valid before, will stay valid
   *       (and may be equal to graph.edge_end() if e was the last edge)
   *
   * @return 1 if e was valid, 0 else
   *
   * Complexity: O(e.node1().degree() + e.node2().degree()) amortized operations
   */
  edge_iterator remove_edge (edge_iterator e_it){
    auto e = *e_it;
    remove_edge(e);
    return e_it;
  }


  /** Remove node n
   * @param[in] n A valid node of this graph
   *
   * @post if old has_node(n), new num_nodes() = old num_nodes() - 1
   *                           new num_edges() = old num_edges() - n.degree()
   * @post if !old has_node(n) then new num_nodes() = old num_nodes()
   *                                new num_edges() = old num_edges()
   * @post new has_node(n) == false
   *
   * @post Invalidates:  Node n
   *                     Edges incident to n
   *                     IncidentIterators pointing at an edge incident to n
   * @post EdgeIterator or IncidentIterator not incident to n pointing at an
   *       invalidated Edge will point at next Edge.
   *       NodeIterator pointing at n will point at next node
   *
   * @return 1 if old has_node(n), 0 else
   *
   * Complexity: O(n.degree() * max degree) amortized operations.
   */
  size_type remove_node (const Node & n){
    if (!has_node(n)){
      return 0;
    }
    else {
      while (n.degree() > 0) {
        remove_edge(*n.edge_begin());
      }
      size_type idx = n.index();
      node_i2u_[idx] = node_i2u_.back();
      node_i2u_.pop_back();
      nodes_[node_i2u_[idx]].idx_ = idx;
      return 1;
    }
  }

  /** Remove node n such that *n_it = n
   * @param[in] n_it A node iterator
   *
   * @post if n_it is valid, new num_nodes() = old num_nodes() - 1
   *                         new num_edges() = old num_edges() - n.degree()
   * @post ifn_it is not valid then new num_nodes() = old num_nodes()
   *                                new num_edges() = old num_edges()
   * @post Invalidates:  Node n
   *                     Edges incident to n
   *                     IncidentIterators pointing at an edge incident to n
   * @post EdgeIterator or IncidentIterator not incident to n pointing at an
   *       invalidated Edge will point at next Edge.
   *       n_it will point at next Node (possible graph.node_end())
   *
   * @return 1 if n was valid, 0 else
   *
   * Complexity: O(n.degree() * max degree) amortized operations.
  */
  node_iterator remove_node (node_iterator n_it) {
    auto n = *n_it;
    remove_node(n);
    return n_it;
  }

 private:

  /** Structure representing a node */
  struct nodeinfo{

    Point position;
    node_value_type value;
    size_type idx_; // User-facing idx

    nodeinfo(Point p, node_value_type v, size_type idx):
    position(p), value(v), idx_(idx) {}
  };

  /** Structure representing an edge */
  struct edgeinfo{
    size_type index1; // Index of a node of the edge
    size_type index2; // Index of other node
    edge_value_type edge_value;
    size_type idx_; // User-facing idx of edge

    edgeinfo(size_type id1, size_type id2, edge_value_type v, size_type idx):
              index1(id1), index2(id2), edge_value(v), idx_(idx) {}
  };


  //Hash function for pair of ints
  // https://www.geeksforgeeks.org/how-to-create-an-unordered_map-of-pairs-in-c/
  struct hash_pair {
    static_assert(sizeof(int) * 2 == sizeof(size_t), "Hash issue");

    size_t operator()(pair p) const noexcept {
        return size_t(p.first) << 32 | p.second;
    }
  };
};

#endif // CME212_GRAPH_HPP
