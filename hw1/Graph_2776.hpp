#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <list>
#include <vector>


#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V>
class Graph {
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
  using node_value_type = V;
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
  Graph() : 
  node_map(), edge_map(),
  node_idx2uid(), edge_idx2uid(),
  next_node_uid_(0), next_edge_uid_(0),
  num_nodes_(0), num_edges_(0){
    // HW0: YOUR CODE HERE
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
  class Node : private totally_ordered<Node>{
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
    Node():
    graph_(NULL), uid_(-1){
    }

    /** Return this node's position. */
    const Point& position() const {
      return this->graph_->node_map[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return this->graph_->node_map[uid_].idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    node_value_type& value(){
      node_value_type& val = this->graph_->node_map.at(uid_).value_;
      return val;
    }

    const node_value_type& value() const{
      const node_value_type& val = this->graph_->node_map.at(uid_).value_;
      return val;
    }

    size_type degree() const{
      size_type deg = this->graph_->node_map[uid_].incident_edges.size();
      return deg;
    }

    incident_iterator edge_begin() const{
      return IncidentIterator(*this, 0);
    }

    incident_iterator edge_end() const{
      return IncidentIterator(Node(), -1);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      bool graph_eq = (this->graph_ == n.graph_);
      bool idx_eq = (this->uid_ == n.uid_);

      return (graph_eq && idx_eq);
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
      assert(n.graph_ != NULL);
      return (this->uid_ < n.uid_);
    }

   // private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;

    Node(const Graph* graph, size_type uid)
      :graph_(const_cast<Graph*>(graph)), uid_(uid){
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return num_nodes_;
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type(-1)) {
    // HW0: YOUR CODE HERE

    //Get index and uid
    size_type uid = next_node_uid_;
    size_type idx = num_nodes_;

    //Add to node_map internal element
    node_map[uid] = internal_node();
    node_map[uid].idx = idx;
    node_map[uid].position = position;
    node_map[uid].value_ = value;

    //Add to idx2uid 
    node_idx2uid[idx] = uid;

    //Increment counters
    ++next_node_uid_;
    ++num_nodes_;
    return Node(this, uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // if (n.graph_ == this){
    //   return true;
    // } else {
    //   return false;
    // }
    return (n.graph_ == this);
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
    return Node(this, node_idx2uid.at(i));
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge(): 
    graph_(NULL), uid_(), reverse_(){
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      if (this->reverse_ == true){
        return this->graph_->edge_map.at(uid_).node2_;
      } else {
        return this->graph_->edge_map.at(uid_).node1_;
      }
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if (this->reverse_ == true){
        return this->graph_->edge_map.at(uid_).node1_;
      } else {
      return this->graph_->edge_map.at(uid_).node2_;
      }
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      bool graph_eq = (this->graph_ == e.graph_);
      bool node_eq;
      if (this->reverse_ ^ e.reverse_){
        node_eq = ((this->node1() == e.node2()) && (this->node2() == e.node1()));
      } else{
        node_eq = ((this->node1() == e.node1()) && (this->node2() == e.node2()));
      }

      // if ((this->node1() == e.node1()) && (this->node2() == e.node2())){
      //   return true;
      // } else if ((this->node1() == e.node2()) && (this->node2() == e.node1())){
      //   return true;
      // } else {
      //   return false;
      // }

      return (graph_eq && node_eq);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      assert(this->graph_ == e.graph_);
      return (this->uid_ < e.uid_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type uid_;
    bool reverse_;

    Edge(const Graph* graph, size_type uid, bool reverse)
      :graph_(const_cast<Graph*>(graph)), uid_(uid), reverse_(reverse){
      }
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, edge_idx2uid.at(i), false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (size_type i = 0; i < num_edges_; i++){
      const Node& node1 = this->edge_map.at(edge_idx2uid.at(i)).node1_;
      const Node& node2 = this->edge_map.at(edge_idx2uid.at(i)).node2_;

      if ((node1.uid_ == a.uid_) && (node2.uid_ == b.uid_)){
      return true;
      } if ((node1.uid_ == b.uid_) && (node2.uid_ == a.uid_)){
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

    bool edge_found = false;
    for (size_type i = 0; i < a.degree(); i++){
      edge_type adj_edge = a.graph_->node_map[a.uid_].incident_edges[i];
      node_type node2 = adj_edge.node2();
      if (b == node2){
        edge_found = true;
        return adj_edge;
      }
    }

    if (!edge_found) {
    size_type edge_uid = next_edge_uid_;
    size_type idx = num_edges_;

    //Add to node_map internal element
    edge_map[edge_uid] = internal_edge();
    if (a < b){
      edge_map[edge_uid].node1_ = a;
      edge_map[edge_uid].node2_ = b;
    } else {
      edge_map[edge_uid].node1_ = b;
      edge_map[edge_uid].node2_ = a;
    }

    edge_map[edge_uid].idx_ = idx;

    //Add to idx2uid 
    edge_idx2uid[idx] = edge_uid;

    //Add to incident lists
    this->node_map[a.uid_].incident_edges.push_back(Edge(this, edge_uid, !(a < b)));
    this->node_map[b.uid_].incident_edges.push_back(Edge(this, edge_uid, (a < b)));
    //Increment counters
    ++next_edge_uid_;
    ++num_edges_;
    return Edge(this, edge_uid, !(a < b));
    }
    return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_map.clear();
    edge_map.clear();
    node_idx2uid.clear();
    edge_idx2uid.clear();
    next_node_uid_ = 0;
    next_edge_uid_ = 0;
    num_nodes_ = 0;
    num_edges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() :
    graph_(nullptr), idx_(){
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const{
      return (*(this->graph_)).node(this->idx_);
    }

    NodeIterator& operator++(){
      if (this->idx_ == this->graph_->num_nodes_ - 1){
        this->graph_ = nullptr;
        this->idx_ = -1;
      } else{
        this->idx_++;

      }
      return *this;
    }

    bool operator==(const NodeIterator& node_iter) const{
      bool graph_eq = (this->graph_ == node_iter.graph_);
      bool idx_eq = (this->idx_ == node_iter.idx_);
      return graph_eq && idx_eq;
    }

   private:
    friend class Graph;

    const Graph* graph_;
    size_type idx_;

    NodeIterator(const Graph* graph, size_type idx)
    : graph_(graph), idx_(idx){}
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const{
    return NodeIterator(nullptr, -1);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator>{
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
    Edge operator*() const{
      return this->node_.graph_->node_map.at(this->node_.uid_).incident_edges[this->idx_];
    }

    IncidentIterator& operator++(){
      if (this->idx_ == this->node_.degree() - 1){
        this->node_ = Node();
        this->idx_ = -1;
      } else{
        (this->idx_)++;
      }
      return *this;
    }

    bool operator==(const IncidentIterator& incident_iter) const{
      if ((int(this->idx_ ) == -1) && (int(incident_iter.idx_) == -1)){
        return true;
      } else if ((int(this->idx_) == -1) ^ (int(incident_iter.idx_) == -1)){
        return false;
      } else{
        bool node_eq = (this->node_ == incident_iter.node_);
        bool iter_eq = (this->idx_ == incident_iter.idx_);
        return (node_eq && iter_eq);
      }

    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Node node_;
    size_type idx_;

    // Construct a valid incident iterator
    IncidentIterator(const Node node, size_type idx)
    : node_(node), idx_(idx){
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() :
    graph_(nullptr), idx_(){
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
      return (*(this->graph_)).edge(this->idx_);
    }

    EdgeIterator& operator++(){
      if (this->idx_ == this->graph_->num_edges_ - 1){
        this->graph_ = nullptr;
        this->idx_ = -1;
      } else{
        this->idx_++;
      }
      return *this;
    }

    bool operator==(const EdgeIterator& edge_iter) const{
      bool graph_eq = (this->graph_ == edge_iter.graph_);
      bool idx_eq = (this->idx_ == edge_iter.idx_);
      return graph_eq && idx_eq;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_;
    size_type idx_;

    EdgeIterator(const Graph* graph, size_type idx)
    :graph_(graph), idx_(idx){
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const{
    return EdgeIterator(nullptr, -1);
  }


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node{
    Point position;
    size_type idx;
    node_value_type value_;
    std::vector<Edge> incident_edges;
  };

  struct internal_edge{
    Node node1_;
    Node node2_;
    size_type idx_;
  };

  std::map<size_type, internal_node> node_map;
  std::map<size_type, internal_edge> edge_map;

  std::map<size_type, size_type> node_idx2uid;
  std::map<size_type, size_type> edge_idx2uid;

  size_type next_node_uid_;
  size_type next_edge_uid_;

  size_type num_nodes_;
  size_type num_edges_;

};

#endif // CME212_GRAPH_HPP
