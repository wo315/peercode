#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

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

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  using node_value_type = V;

  /** Type of this graph. */
  using graph_type = Graph<V>;

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
  Graph(): points_(), nodes_(), edges_(), values_(), n_size_(0), e_size_(0), adjacency_map(){
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // graph_ points to memory location, to access its members use ->
      return graph_->points_[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }


    /** function returns a reference to the value of a node, 
     *  so you can edit it
     * @return reference of the value of node
     *
     * Complexity: O(1)
    */
    node_value_type& value(){
      return graph_->values_[uid_];
    }


    /** function returns the value of a node
     * @return value of node
     *
     * Complexity: O(1)
    */
    const node_value_type& value() const{
      return graph_->values_[uid_];
    }

    /** function returns the degree of a node
     * @return degree of node
     *
     * Complexity: O(1)
    */
    size_type degree() const{
      return graph_->adjacency_map.at(uid_).size();

    }


    /** function returns iterator to incident nodes
     * @return beginning of IncidentIterator
     *
     * Complexity: O(1)
    */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, graph_->adjacency_map.at(this->uid_).begin(), uid_);
    }


    /** function returns iterator to incident nodes
     * @return end of IncidentIterator
     *
     * Complexity: O(1)
    */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, graph_->adjacency_map.at(this->uid_).end(), uid_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ( n.index() == this->index() and n.graph_ == this->graph_);
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
      return this->index() < n.index();
    }

   private:
      // Use this space to declare private data members and methods for Node
      // that will not be visible to users, but may be useful within Graph.
      // i.e. Graph needs a way to construct valid Node objects

      // Pointer back to the Graph container
      Graph* graph_;
      // This Node's unique identification number
      size_type uid_;

      //Node(const Graph* graph, size_type uid)
      //    : graph_(const_cast<Graph*>(graph)), uid_(uid), value_(0) {
      //}

      Node(const Graph* graph, size_type uid, node_value_type value)
          : graph_(const_cast<Graph*>(graph)), uid_(uid){
            graph_->values_.push_back(value);
      }

      // Allow Graph to access Node's private member data and functions.
      friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return n_size_;
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    points_.push_back(position);
    ++n_size_;
    nodes_.push_back(Node(this, n_size_ - 1, value));
    // initialize entry in map for new node
    adjacency_map[n_size_-1] = std::map<size_type, size_type>();
    return nodes_[n_size_-1];
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    unsigned idx = n.index(); 
    return (idx  < n_size_ and nodes_[idx] == n );
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());  
    return nodes_[i] ;     
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return this->graph_->nodes_[uid_n1_];      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return this->graph_->nodes_[uid_n2_];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      size_type e1_n1 = this->node1().index();
      size_type e1_n2 = this->node2().index();
      size_type e2_n1 = e.node1().index();
      size_type e2_n2 = e.node2().index();
      return ((e2_n1 == e1_n1 and e2_n2 == e1_n2)  or (e2_n1 == e1_n2 and e2_n2 == e1_n1));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (this->uid_e_ < e.uid_e_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

      Graph* graph_;
      // This element's unique identification number
      size_type uid_n1_, uid_n2_, uid_e_;

      Edge(const Graph* graph, size_type uid_n1, size_type uid_n2, size_type uid_e)
          : graph_(const_cast<Graph*>(graph)), uid_n1_(uid_n1), uid_n2_(uid_n2), uid_e_(uid_e) {
      }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return e_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return this->edges_[i];        // valid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type a1 = a.index();
    size_type b1 = b.index();

    if (adjacency_map.count(a1) == 1){

      if (adjacency_map.at(a1).count(b1) == 1){
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
    size_type idx = 0;
    Edge e;
    // see if the edge exists
    size_type a1 = a.index();
    size_type b1 = b.index();

    if  (this->has_edge(a,b) == true){
      

      idx = adjacency_map.at(a1).at(b1);
    }
    else{

      ++e_size_;
      edges_.push_back(Edge(this, a1, b1, e_size_ - 1));
      idx = e_size_ - 1;
      // update for edge (a,b) and (b,a)
      adjacency_map[a1][b1] = idx;
      adjacency_map[b1][a1] = idx;
    }
    return edges_[idx];
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

    std::vector<Point>().swap(points_);
    std::vector<Node>().swap(nodes_);
    std::vector<Edge>().swap(edges_);
    std::vector<node_value_type>().swap(values_);
    adjacency_map = std::map<size_type, std::map<size_type, size_type>>();

    n_size_ = 0;
    e_size_ = 0;
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

  
    /**
      * Dereference iterator
      * @return Node for the current Nodeiterator
      *
      * Complexity: O(1)
    */
    Node operator*() const{
      node_type node = *node_iterator_;
      return node;
    }
    

    /**
      * Return reference to Node iterator+1
      * @return reference for the next Nodeiterator
      *
      * Complexity: O(1)
    */
    NodeIterator& operator++(){
      node_iterator_++;
      return *this;
    }
    
    /**
      * check if two iterators are the same
      *
      * Complexity: O(1)
    */
    bool operator==(const NodeIterator& node_iter2) const{
      return (this->node_iterator_ == node_iter2.node_iterator_);
    }

    /**
      * check if two iterators are the different
      *
      * Complexity: O(1)
    */
    bool operator!=(const NodeIterator& node_iter2) const{
      return (this->node_iterator_ != node_iter2.node_iterator_);
    }

   private:
    friend class Graph;
    
    typename std::vector<Node>::const_iterator node_iterator_;

    Graph* graph_;

    NodeIterator(const Graph* graph, typename std::vector<Node>::const_iterator node_iterator)
          : node_iterator_(node_iterator), graph_(const_cast<Graph*>(graph)) {
    }

  };

  /**
    * @return iterator for the first node
    *
    * Complexity: O(1)
  */
  node_iterator node_begin() const{
    return NodeIterator(this, this->nodes_.begin());
  }


  /**
    * @return iterator for the last node
    *
    * Complexity: O(1)
  */
  node_iterator node_end() const{
    return NodeIterator(this, this->nodes_.end());
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
    }

    
    /**
      * Dereference iterator
      * @return Edge for the current Incident edge iterator
      *
      * Complexity: O(1)
    */
    Edge operator*() const{
      // (this, a1, b1, e_size_ - 1)
      return Edge(graph_,  node_uid_, incid_iterator_->first, incid_iterator_->second);
    }


    /**
      * Return reference to Incident iterator+1
      * @return reference for the next incident iterator
      *
      * Complexity: O(1)
    */
    IncidentIterator& operator++(){
      incid_iterator_++;
      return *this;
    }


    /**
      * check if two iterators are the same
      *
      * Complexity: O(1)
    */
    bool operator==(const IncidentIterator& incid_iter2) const{
      return (this->incid_iterator_ == incid_iter2.incid_iterator_);
    }

    /**
      * check if two iterators are the different
      *
      * Complexity: O(1)
    */
    bool operator!=(const IncidentIterator& incid_iter2) const{
      return (this->incid_iterator_ != incid_iter2.incid_iterator_);
    }

   private:
    friend class Graph;
    
    // iterator for the map
    typename std::map<size_type, size_type>::const_iterator incid_iterator_;

    size_type node_uid_;

    Graph* graph_;

    IncidentIterator(const Graph* graph, typename std::map<size_type, size_type>::const_iterator incid_iterator, size_type  node_uid)
          : incid_iterator_(incid_iterator), node_uid_(node_uid), graph_(const_cast<Graph*>(graph)) {
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
    }

    
    /**
      * Dereference iterator
      * @return Edge for the current Edgeiterator
      *
      * Complexity: O(1)
    */
    Edge operator*() const{
      return *edge_iterator_;
    }


    /**
      * Return reference to Edge iterator+1
      * @return reference for the next Edgeiterator
      *
      * Complexity: O(1)
    */
    EdgeIterator& operator++(){
      edge_iterator_++;
      return *this;
    }


    /**
      * check if two iterators are the same
      *
      * Complexity: O(1)
    */
    bool operator==(const EdgeIterator& edge_iter2) const{
      return (this->edge_iterator_ == edge_iter2.edge_iterator_);
    }

    /**
      * check if two iterators are the different
      *
      * Complexity: O(1)
    */
    bool operator!=(const EdgeIterator& edge_iter2) const{
      return (this->edge_iterator_ != edge_iter2.edge_iterator_);
    }

   private:
    friend class Graph;
    // iterator for the map
    typename std::vector<Edge>::const_iterator edge_iterator_;

    size_type node_uid_;

    Graph* graph_;

    EdgeIterator(const Graph* graph, typename std::vector<Edge>::const_iterator edge_iterator)
          : edge_iterator_(edge_iterator),graph_(const_cast<Graph*>(graph)) {
    }
  };


  /**
    * @return iterator for the first edge
    *
    * Complexity: O(1)
  */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, this->edges_.begin());
  }


  /**
    * @return iterator for the last edge
    *
    * Complexity: O(1)
  */
  edge_iterator edge_end() const{
    return EdgeIterator(this, this->edges_.end());
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::vector<Point> points_;
  std::vector<Node> nodes_;
  std::vector<Edge> edges_;
  std::vector<node_value_type> values_;

  size_type n_size_, e_size_;

  // map: {node_a_idx : {node_b_idx: edge_ab_idx, ...}, ...}
  std::map<size_type, std::map<size_type, size_type>> adjacency_map;

};

#endif // CME212_GRAPH_HPP
