#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

//--functionality_0
//--great job!
//--END


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
  typedef V node_value_type;
  typedef E edge_value_type;

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct internal_node;
  struct internal_edge;

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
  Graph(): next_id_(0), num_edges_(0), nodes_(), i2u_(), adj_() {}

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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */

    Point& position() {
      return fetch().position_;
    }

    const Point& position() const {
      return fetch().position_;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return fetch().idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    node_value_type& value(){
      return fetch().value_;
    }

    const node_value_type& value() const {
      return fetch().value_;
    }

    size_type degree() const {
      return graph_->adj_[uid_].size();
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, index(), 0);
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, index(), degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_==n.graph_ && uid_==(n.uid_));
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
      return (graph_<n.graph_) || (graph_==n.graph_ && uid_<(n.uid_));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;

    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {}

    // Return the internal node associated to Node
    internal_node& fetch() const {
      return graph_->nodes_[uid_];
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return i2u_.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // Define the uid_ of the new node using next_id_ (which is unused)
    // and incrementing this value
    size_type uid = next_id_;
    next_id_++;

    // Define the idx_ of the new node which is the size of the graph
    size_type idx = size();

    // Create a new internal node
    internal_node internal_node_(uid, idx, position, val);

    // The internal node is added to the vector: easy retrieval from the uid_
    nodes_.push_back(internal_node_);

    // Put in i2u_
    i2u_.push_back(uid);

    // Intialize the adjacency vector of the new node
    std::map<size_type, size_type> adjacent_edges;
    adj_.push_back(adjacent_edges);

    return Node(this, uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_ && n.index()<size());
  }

  /** Removes a node from the graph
   * @param[in] n the node to remove
   * @return 1 if the node was in the graph, 0 otherwise.
   * @post new num_nodes() == old num_nodes()-1 if n was in the graph
   * @post new num_edges() == old num_edges()-old n.degree() if n was in the graph
   * Invalidates node n
   * Complexity: O(old num_nodes())
   */

  size_type remove_node(const Node& n) {
   if(!has_node(n)){
     return 0;
   }
   for(auto it = adj_[n.uid_].begin(); it != adj_[n.uid_].end(); it++) {
     adj_[it->first].erase(n.uid_);
     num_edges_ -= 1;
   }
   adj_[n.uid_].clear();
   size_type node_idx = n.index();
   size_type last_node_uid = i2u_[i2u_.size()-1];
   nodes_[last_node_uid].idx_ = node_idx;
   i2u_[node_idx] = last_node_uid;
   i2u_.pop_back();
   return 1;
  }

  /** Removes a node from the graph
   * @param[in] n_it the iterator to the node to remove
   * @return a new valid iterator
   * @post new num_nodes() == old num_nodes()-1
   * @post new num_edges() == old num_edges()-old (*n_it).degree()
   * Invalidates node *n_it
   * Complexity: O(old num_nodes())
   */
   node_iterator remove_node(node_iterator n_it) {
     remove_node(*n_it);
     return n_it;
   }

   /** Return the node with index @a i.
    * @pre 0 <= @a i < num_nodes()
    * @post result_node.index() == i
    *
    * Complexity: O(1).
    */
  Node node(size_type i) const {
    return Node(this, i2u_.at(i));
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
    Edge() {
      // HW0: YOUR CODE HERE

    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_uid_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(graph_!=e.graph_){
        return false;
      }
      size_type e1n1_uid = fetch().node1_uid_;
      size_type e1n2_uid = fetch().node2_uid_;
      size_type e2n1_uid = fetch().node1_uid_;
      size_type e2n2_uid = fetch().node2_uid_;
      return (e1n1_uid==e2n1_uid&&e1n2_uid==e2n2_uid)||(e1n1_uid==e2n2_uid&&e1n2_uid==e2n1_uid);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return node1()<e.node1() || (node1()==e.node1() && node2()<e.node2());
    }

    double length() const {
      return norm(node1().position()-node2().position());
    }

    edge_value_type& value() {
      return fetch().value_;
    }

    const edge_value_type& value() const {
      return fetch().value_;
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type node1_uid_;
    size_type node2_uid_;

    Edge(const Graph* graph, size_type node1_uid, size_type node2_uid)
        : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid) {}

    // Returns the internal edge associated to the edge
    internal_edge& fetch() const {
      size_type edge_idx = graph_->adj_.at(node1_uid_).find(node2_uid_)->second;
      return graph_->edges_.at(edge_idx);
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, edges_.at(i).node1_uid_, edges_.at(i).node2_uid_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Test in the adjacency lists
    // Firstly test if a is a node of the graph
    if(has_node(a)&&has_node(b)){
      // Then tests if b is in the adjacency list of a
      return adj_.at(a.uid_).find(b.uid_) != adj_.at(a.uid_).end();
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    if(!has_edge(a, b)){ // Only add an edge if it is not in the graph
      // Create the associated internal edge
      internal_edge internal_edge_(a.uid_, b.uid_, val);
      // Add the edge to the adjacency lists of a and b
      adj_.at(a.uid_).emplace(b.uid_, edges_.size());
      adj_.at(b.uid_).emplace(a.uid_, edges_.size());

      // Add the edge to the list of edges
      edges_.push_back(internal_edge_);
      num_edges_ += 1;
    }
    return Edge(this, a.uid_, b.uid_);
  }

  /** Removes an edge from the graph
   * @param[in] a is a node from the edge to remove
   * @param[in] b is the other node from the edge to remove
   * @return 1 if the edge was in the graph, 0 otherwise.
   * @post new num_edges() == old num_edges()-1 if the edge was in the graph
   * Invalidates edge (a,b)
   * Complexity: O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if(!has_edge(a, b)){
      return 0;
    }
    adj_[a.uid_].erase(b.uid_);
    adj_[b.uid_].erase(a.uid_);
    num_edges_ -= 1;
    return 1;
  }

  /** Removes an edge from the graph
   * @param[in] e the edge to remove
   * @return 1 if the edge was in the graph, 0 otherwise.
   * @post new num_edges() == old num_edges()-1 if the edge was in the graph
   * Invalidates edge e
   * Complexity: O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Removes an edge from the graph
   * @param[in] e_it an iterator to the edge to remove
   * @return a valid edge iterator
   * @post new num_edges() == old num_edges()-1
   * Invalidates edge *e_it
   * Complexity: O(num_nodes() + num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    size_type edge_index = e_it.edge_index_;
    remove_edge(*e_it);
    return EdgeIterator(this, edge_index);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    i2u_.clear();
    adj_.clear();
    next_id_ = 0;
    num_edges_ = 0;
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

    Node operator*() const {
      return Node(graph_, graph_->i2u_[index_]);
    }

    NodeIterator& operator++(){
      index_++;
      return *this;
    }

    bool operator==(const NodeIterator& it) const{
      return (graph_==it.graph_)&&(index_==it.index_);
    }

   private:
    friend class Graph;
    // A NodeIterator is a graph and the index in the graph of the current node
    Graph* graph_;
    size_type index_;
    NodeIterator(const Graph* graph, size_type index)
    : graph_(const_cast<Graph*>(graph)), index_(index) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  node_iterator node_end() const {
    return NodeIterator(this, size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
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

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    Edge operator*() const {
      return Edge(graph_, graph_->i2u_[node_idx_], internal_it_->first);
    }
    IncidentIterator& operator++() {
      index_++;
      ++internal_it_;
      return *this;
    }
    bool operator==(const IncidentIterator& it) const {
      return(graph_==it.graph_)&&(node_idx_==it.node_idx_)&&(index_==it.index_);
    }

   private:
    friend class Graph;
    // An incident IncidentIterator is a graph, the uid_ of the indicent node
    // the index of the current edge in the adjacency list and an internal iterator to the current edge
    Graph* graph_;
    size_type node_idx_;
    size_type index_; // from 0 to degree of the node
    typename std::map<size_type, size_type>::iterator internal_it_;

    IncidentIterator(const Graph* graph, size_type node_idx, size_type index)
    : graph_(const_cast<Graph*>(graph)), node_idx_(node_idx), index_(index), internal_it_() {
      internal_it_ = graph_->adj_[graph_->i2u_[node_idx]].begin();
      for(size_type i=0; i<index; i++) {
        internal_it_++;
      }
    }
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

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    Edge operator*() const {
      return Edge(graph_, internal_it_->node1_uid_, internal_it_->node2_uid_);
    }

    EdgeIterator& operator++(){
      if(edge_index_>=graph_->num_edges()-1){
        edge_index_++;
      }
      else{
        do{
          ++internal_it_;
        } while(!isValid(internal_it_));
        edge_index_++;
      }
      return *this;
    }

    bool operator==(const EdgeIterator& it) const{
      return graph_==it.graph_ && edge_index_==it.edge_index_;
    }

   private:
    friend class Graph;
    // An EdgeIterator is a graph and the index in the graph of the current edge
    Graph* graph_;
    size_type edge_index_;
    typename std::vector<internal_edge>::iterator internal_it_;
    EdgeIterator(const Graph* graph, size_type edge_index)
    : graph_(const_cast<Graph*>(graph)), edge_index_(edge_index), internal_it_() {
      if(edge_index_==graph->num_edges()) return;
      internal_it_ = graph_->edges_.begin();
      while(!isValid(internal_it_)){
        ++internal_it_;
      }
      for(size_type i=0; i<edge_index_; i++){
        do{
          ++internal_it_;
        } while(!isValid(internal_it_));
      }
    }

    bool isValid(typename std::vector<internal_edge>::iterator it){
      return graph_->has_edge(Node(graph_, it->node1_uid_), Node(graph_, it->node2_uid_));
    }

  };

  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  edge_iterator edge_begin(){
    return EdgeIterator(this, 0);
  }

  edge_iterator edge_end(){
    return EdgeIterator(this, num_edges());
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node {
    // Members
    size_type uid_;
    size_type idx_;
    Point position_;
    node_value_type value_;
    // Constructors
    internal_node() {}
    internal_node(size_type uid, size_type idx, Point position, node_value_type value):
    uid_(uid), idx_(idx), position_(position), value_(value) {}
  };

  struct internal_edge {
    // Members
    size_type node1_uid_;
    size_type node2_uid_;
    edge_value_type value_;
    // Constructor
    internal_edge(size_type node1_uid, size_type node2_uid, edge_value_type value):
    node1_uid_(node1_uid), node2_uid_(node2_uid), value_(value) {}
  };

  size_type next_id_; // Store and increment to ensure each node as a unique uid_
  // Nodes contain a uid_ (never modified) and an idx_.
  // If 0 <= idx_ < graph.size(): the node is in the graph
  // If idx_ >= graph.size(): the node was removed

  size_type num_edges_;

  std::vector<internal_node> nodes_;  // Indexed by node uid_

  // i2u_ allows to find quickly the node idexed by i
  std::vector<size_type> i2u_ ; // Maps idx_ to uid_

  // Edges are stored in a vector
  std::vector<internal_edge> edges_;  // Indexed by edge uid_
  // We use adjacency lists for incident iteration
  // Vector of vectors of pair
  //adj_[node_uid1_]: vector of pairs of type (node_uid2_, edge uid_)
  std::vector<std::map<size_type, size_type>> adj_;


};

#endif // CME212_GRAPH_HPP
