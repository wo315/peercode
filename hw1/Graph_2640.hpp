#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>
#include <unordered_map>

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

  /* Nodes bank, STL container containing all the nodes of the graph*/
  struct node_wrapper;
  std::vector<node_wrapper> nodes_;
  /* Edges bank, we are using a simple wrapper to store the indices of both
  nodes as well as the index of the edge*/
  struct edge_wrapper;
  std::vector<edge_wrapper> edges_;
  /* map index of nodes to index of neighbors, nested list */
  std::unordered_map<unsigned, std::vector<unsigned int>> nodes_map_;
  /* map index of edges to index of nodes */
  std::unordered_map<unsigned, std::pair<unsigned int, unsigned int>>edges_map_; 
  

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
    }


    /** Return this node's position. */
    const Point& position() const {
      // return Point();
      return (graph_ -> nodes_[uid_].point);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    node_value_type& value(){
      return graph_->nodes_[uid_].value;
    }

    const node_value_type& value() const{
      return graph_->nodes_[uid_].value;
    }

    size_type degree() const{
      return graph_->nodes_map_[graph_->nodes_[uid_].node_index].size();
    }

    incident_iterator edge_begin() const{
      return IncidentIterator(this->graph_, uid_, 0);
    }

    incident_iterator edge_end() const{
      std::vector<unsigned int> nei;
      nei = (this -> graph_ -> nodes_map_)[this -> uid_];
      return IncidentIterator(this->graph_, uid_, nei.size());
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if(graph_ == n.graph_ and uid_ == n.uid_){
        return true;
      }
      else{
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
      assert(n.graph_ != NULL);
      if(graph_ == n.graph_ and uid_ < n.uid_){
        return true;
      }
      else if(graph_ < n.graph_){ // breaking ties
        return true;
      }
      else{
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_; //SimpleSet
    size_type uid_; //unique id of the node in the graph
    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
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
  Node add_node(const Point& position, const node_value_type& nvt = node_value_type()) {
    Node node(this, size());
    nodes_.push_back(node_wrapper(size(), position, nvt));
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Needs to be updated when we ll add a remove function
    if (!(this == n.graph_)){
      return false;
    }
    else if (this -> size() < n.uid_){
      return false;
    }
    else{
      return true;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      graph_ = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, uid1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, uid2_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((node1() == e.node1()) && (node2() == e.node2())){
        return true;
      }
      else if ((node2() == e.node1()) && (node1() == e.node2())){
        return true;
      }
      else{
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order. 
     * Uses a minimum to ensure that (a, b) and (b, a) have the same 
     * behaviour
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if(graph_<e.graph_){
        return true;
      }
      unsigned int min = std::min(e.node1().uid_, e.node2().uid_);
      if (uid1_ < min || uid2_ < min){
        return true;
      }
      else{
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type uid1_;
    size_type uid2_;

    /* Private constructor */
    Edge(const graph_type* graph, size_type uid1, size_type uid2)
		: graph_(const_cast<graph_type*>(graph)), uid1_(uid1), uid2_(uid2) {
	}
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i < num_edges()){
      return Edge(this, this -> edges_[i].node_id1, this -> edges_[i].node_id2);
    }
    else{
      return Edge();        // Invalid Edge
    }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (nodes_map_.find(a.uid_) == nodes_map_.end()){
      return false;
    }
    else{
      std::vector<size_type> neighbors = nodes_map_.at(a.uid_);
      for (unsigned i = 0; i < neighbors.size(); ++i){
        if (neighbors[i] == b.uid_){
          return true;
        }
      }
      return false;
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
    int number_of_edges = edges_.size();
    if (has_edge(a, b)){
      return Edge(this, a.uid_, b.uid_);
    }
    std::pair<size_type, size_type> pair_nodes;
    pair_nodes.first = a.uid_;
    pair_nodes.second = b.uid_;
    edges_map_[number_of_edges] = pair_nodes;
    if (nodes_map_.find(a.uid_) == nodes_map_.end()){ //node wasn't connected
      std::vector<size_type> neighbors_a;
      neighbors_a.push_back(b.uid_);
      nodes_map_[a.uid_] = neighbors_a;
    }
    else{
      nodes_map_[a.uid_].push_back(b.uid_);
    }
    if (nodes_map_.find(b.uid_) == nodes_map_.end()){
      std::vector<size_type> neighbors_b;
      neighbors_b.push_back(a.uid_);
      nodes_map_[b.uid_] = neighbors_b;
    }
    else{
      nodes_map_[b.uid_].push_back(a.uid_);
    }
    edge_wrapper edge_to_be_added(number_of_edges, a.uid_, b.uid_);
    this -> edges_.push_back(edge_to_be_added);
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    nodes_map_.clear();
    edges_map_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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

    /** Dereference operator for NodeIterator
     * @return an Node object n
     * Complexity: Constant
     */
    Node operator*() const{
      return Node(graph_, node_iterator_index_);
    }

    /** Incrementation operator for NodeIterator
     * @return a reference to a NodeIterator object
     * Complexity: Constant
     */
    NodeIterator& operator++(){
      ++node_iterator_index_;
      return *this;
    }

    /** Equality operator for NodeIterator
     * @param nodeiterator NodeIterator to be compared with this
     * @return a bool that is true if and only if this and @nodeiterator
     * has same index and same graph than this
     * Complexity: Constant
     */
    bool operator==(const NodeIterator& nodeiterator) const{
      return (nodeiterator.graph_ == this->graph_ &&
      this -> node_iterator_index_ == nodeiterator.node_iterator_index_);
    }

   private:
    friend class Graph;
    size_type node_iterator_index_;
    graph_type* graph_;
    NodeIterator(size_type index, const graph_type* graph) : 
    node_iterator_index_(index), graph_(const_cast<graph_type*>(graph)) {
    }
  };

  /**
   * @brief Start NoteIterator with index 0
   * 
   * @return node_iterator with graph_ == this and index 0
   */
  node_iterator node_begin() const{
    return NodeIterator(0, this);
  }
  /**
   * @brief End NoteIterator with index this->num_nodes()
   * 
   * @return node_iterator with graph_ == this and index 0
   */
  node_iterator node_end() const{
    return NodeIterator(this->num_nodes(), this);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator>  {
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
    /** Dereference operator for IncidentIterator
     * @return a Edge object e that is gives e.node1().index() == node_id
     * Complexity: Constant
     */
    Edge operator*() const{
      std::vector<unsigned> neighbors;
      neighbors = this->graph_->nodes_map_[node_id_];
      return Edge(graph_, node_id_, neighbors[neighbor_index_]);
    }

    /** Incrementation operator for IncidentIterator
     * @return a reference to a NodeIterator object where the incident is the 
     * next edge
     * Complexity: Constant
     */
    IncidentIterator& operator++(){
      ++neighbor_index_;
      return *this;
    }
    
    /**
     * @brief test if incidentiterator and this are the same
     * 
     * @param incidentiterator IncidentIterator object to be compared
     * @return true if graph, original node and destination node 
     * are the same
     * @return false otherwise
     */
    bool operator==(const IncidentIterator& incidentiterator) const{
      return this->graph_ == incidentiterator.graph_ && this->node_id_ ==
      incidentiterator.node_id_ && this->neighbor_index_ == 
      incidentiterator.neighbor_index_;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type node_id_; // node_id of node of reference
    size_type neighbor_index_;
    IncidentIterator(const graph_type* graph, size_type node_id, size_type neighbor_index) : 
    graph_(const_cast<graph_type*>(graph)), node_id_(node_id), neighbor_index_(neighbor_index){
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
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

    /**
     * @brief Dereference operator for EdgeIterator
     * 
     * @return Edge object e with e.node1() and e.node2() are the nodes
     * corresponding to the index edge_iterator_index_
     */
    Edge operator*() const{
      std::pair<size_type, size_type> edge;
      edge = this->graph_->edges_map_[this->edge_iterator_index_];
      return Edge(this->graph_, edge.first, edge.second);
    }
    /**
     * @brief Incrementation operator
     * @pre index is smaller than this->num_edges()
     * 
     * @return EdgeIterator& with index incremented by one
     */
    EdgeIterator& operator++(){
      ++edge_iterator_index_;
      return *this;
    }
    /**
     * @brief Equality operator for EdgeIterator
     * 
     * @param edgeiterator EdgeIterator object to be compared with this
     * @return true if graph and index are identical
     * @return false otherwise
     */
    bool operator==(const EdgeIterator& edgeiterator) const{
      return (edgeiterator.graph_ == this->graph_ && 
      edgeiterator.edge_iterator_index_ == this -> edge_iterator_index_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type edge_iterator_index_;
    EdgeIterator(const Graph* graph, size_type id)
    : graph_(const_cast<Graph*>(graph)), edge_iterator_index_(id) {
    }
  };

  /**
   * @brief Start EdgeIterator, the starting index is 0 and the graph is 
   * given by this.
   * 
   * 
   * @return edge_iterator starting at the first edge
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /**
   * @brief End EdgeIterator, the ending index is graph.num_edges()
   * 
   * @return edge_iterator that finishes the iteration
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this, this->num_edges());
  }


 private:

  /* wrapper for the node information */
  struct node_wrapper {
    size_type node_index;
    Point point;
    node_value_type value;
    
    /* public constructor for the node wrapper */
    node_wrapper(size_type i, Point P, node_value_type val) 
			 :node_index(i),point(P), value(val) {}
  };

  /* wrapper for the edge information */
  struct edge_wrapper {
    size_type edge_index;
		size_type node_id1;
		size_type node_id2;

    
    // public constructor of edge wrapper
		edge_wrapper(size_type ei, size_type id_1, size_type id_2) 
			:edge_index(ei), node_id1(id_1), node_id2(id_2){}
	 };
};

#endif // CME212_GRAPH_HPP
