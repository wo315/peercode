#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_set>
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

  // HW0: YOUR CODE HERE
  struct internal_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  // HW1: Modifiable Node Value and Class Templates
  using node_value_type = V;

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
  Graph() 
    : next_uid_(0), nodes_(), edges_(), edges_set_(), adjacency_() {
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
      std::cout << "Invalid node: you should use Graph methods" << std::endl;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_index_;
    }

    // HW1: YOUR CODE HERE
    /** Return this node's value. */
    node_value_type& value(){
      return (*graph_).nodes_.at(node_index_).node_value_;
    }

    /** Return this node's value.
     * We guaranty that the node is not modified and the returned value is const. */
    const node_value_type& value() const{
      return (*graph_).nodes_.at(node_index_).node_value_;
    }

    /** Return this node's degree, the number of edges connected to that node. */
    size_type degree() const{
      return graph_->adjacency_.at(node_index_).size();
    }

    /** Return an IncidentIterator to iterate over the edges connected to that node. 
     * The node that spwans the incident_iterator is returned by node1() of each incident Edge
     * and the adjacent Node is returned by node2().
    */

    /** Return an IncidentIterator pointing to the node's first adjacent edge: 
     * an edge between the node and the first node in graph_->adjacency_.at(node_index_).
     * If dereferenced, node1() returns the node and node2() returns the adjacent node. */   
    IncidentIterator edge_begin() const{
      return IncidentIterator(graph_, node_index_, 0);
    }

    /** Return an IncidentIterator pointing to the node's first adjacent edge: 
     * an edge between the node and the last node in graph_->adjacency_.at(node_index_).
     * If dereferenced, node1() returns the node and node2() returns the adjacent node. */ 
    IncidentIterator edge_end() const{
      return IncidentIterator(this->graph_, node_index_, graph_->adjacency_.at(node_index_).size());
    }    

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_ && this->index() == n.index());
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
      return ( (graph_ < n.graph_) || (graph_ == n.graph_ && this->index() < n.index()) );
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_;  // Pointer back to the Graph
    size_type node_index_;

    // Private constructor to create a valid Node
    Node(const Graph* graph, size_type node_index)
      : graph_(const_cast<Graph*>(graph)), node_index_(node_index){
    }

    /** Return the node at index node_index. */
    internal_node& fetch() const {
      return graph_->nodes_.at(node_index_);
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    return add_node(position, node_value_type());

  }

  /** Add a node to the graph with value node_value, returning the added node. */
  Node add_node(const Point& position, const node_value_type& node_value) {
    // HW1: Your code
    internal_node *node = new internal_node;
    node->position_ = position;
    node->uid_ = next_uid_;
    node->node_value_ = node_value;
    next_uid_++;
    nodes_.push_back(*node);

    std::vector<Node> vec;
    adjacency_.push_back(vec);

    return Node(this, this->size()-1);
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ( this == n.graph_ && n.index() < size() );
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      std::cout << "Invalid edge: you should use Graph methods" << std::endl;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return *node1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return *node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (graph_ == e.graph_) \
        && ( (node1().index() == e.node1().index() && node2().index() == e.node2().index()) \
        || (node1().index() == e.node2().index() && node2().index() == e.node1().index()) );

    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if(graph_ < e.graph_) return true;
      else if(graph_ == e.graph_){
        if(this->node1().index() < e.node1().index()) return true;
        else if(this->node1().index() == e.node1().index()){
          if(this->node2().index() < e.node2().index()) return true;
        }
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

    Graph* graph_;
    const Node* node1_;
    const Node* node2_;

    Edge(const Graph* graph, const Node* node1, const Node* node2) 
      : graph_(const_cast<Graph*>(graph)), node1_(node1), node2_(node2) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return edges_.at(i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    Edge edge(this, &a, &b);

    if(edges_set_.find(edge) != edges_set_.end()){
      return true;
    }
    else{
      Edge edge(this, &b, &a);
      return (edges_set_.find(edge) != edges_set_.end());
    }

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
    Edge* edge = new Edge(this, &a, &b);

    if(has_edge(a, b) == false){
      edges_.push_back(*edge);
      edges_set_.emplace(*edge);

      adjacency_.at(a.index()).push_back(b);
      adjacency_.at(b.index()).push_back(a);  
    }
    return *edge;

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
    edges_set_.clear();
    adjacency_.clear();
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
    /** Return the node to which the node iterator is pointing to. */   
    Node operator*() const{
      // HW1: Your code
      return Node(graph_, node_index_);
    }

    /** Increment the node iterator and dereference it. */   
    NodeIterator& operator++(){
      // HW1: Your code
      ++node_index_;
      return *this;
    }

    /** Check equality between node iterators. */
    bool operator==(const NodeIterator& node_iter) const {
      // HW1: Your code
      return (graph_ == node_iter.graph_) && (node_index_ == node_iter.node_index_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_index_;

    NodeIterator(const Graph* graph, size_type node_index) :
      graph_(const_cast<Graph*>(graph)), node_index_(node_index) {}

  };

  // HW1 #2: YOUR CODE HERE
  /** Return a node iterator for the graph's nodes that points to the first node in nodes_. */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Return a node iterator for the graph's nodes that points to the last node in nodes_. */
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

    // HW1 #3: YOUR CODE HERE
    /** Return the edge to which the incident iterator is pointing to. */
    Edge operator*() const{
      // HW1: Your code
      std::vector<Node> adjacent_nodes = graph_->adjacency_.at(node_index_);
      Node node(graph_, node_index_);
      Node adjacent_node = adjacent_nodes.at(adjacent_node_index_);
      Edge edge(this->graph_, &node, &adjacent_node);
      return edge;
    }

    /** Increment the incident iterator and derefence it. */
    IncidentIterator& operator++(){
      // HW1: Your code
      ++adjacent_node_index_;
      return *this;
    }

    /** Check equality between two incident iterators. */
    bool operator==(const IncidentIterator& incident_iter) const{
      // HW1: Your code
      return (graph_ == incident_iter.graph_) && (node_index_ == incident_iter.node_index_) \
        && (adjacent_node_index_ == incident_iter.adjacent_node_index_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node_index_;
    size_type adjacent_node_index_;

    IncidentIterator(const Graph* graph, size_type node_index, size_type adjacent_node_index) :
      graph_(const_cast<Graph*>(graph)), node_index_(node_index), adjacent_node_index_(adjacent_node_index) {}


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
    /** Return the edge to which the edge iterator is pointing to. */
    Edge operator*() const{
      // HW1: Your code
      return graph_->edges_.at(edge_index_);
    }

    /** Increment the edge iterator and derefence it. */
    EdgeIterator& operator++(){
      // HW1: Your code
      ++edge_index_;
      return *this;
    }

    /** Check equality between two edge iterators. */
    bool operator==(const EdgeIterator& edge_iter) const{
      // HW1: Your code
      return (graph_ == edge_iter.graph_) && (edge_index_ == edge_iter.edge_index_);
    }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type edge_index_;

    EdgeIterator(const Graph* graph, size_type edge_index) :
      graph_(const_cast<Graph*>(graph)), edge_index_(edge_index) {}
    
  };

  // HW1 #5: YOUR CODE HERE
  /** Return an edge iterator for the graph's edges that points to the first edge in edges_. */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /** Return an edge iterator for the graph's edges that points to the last edge in edges_. */
  edge_iterator edge_end() const{
    return EdgeIterator(this, this->edges_.size());
  }
  

  class EdgeHash{
    public:
      size_t operator()(const Edge& t) const{
        return t.node1().node_index_ + 21269*t.node2().node_index_;
      }
  };

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node{
    Point position_;
    size_type uid_; // Unique id
    node_value_type node_value_;
  };

  size_type next_uid_;
  std::vector<internal_node> nodes_;
  std::vector<Edge> edges_;
  std::unordered_set<Edge, EdgeHash> edges_set_;
  std::vector<std::vector<Node>> adjacency_;

};

#endif // CME212_GRAPH_HPP
