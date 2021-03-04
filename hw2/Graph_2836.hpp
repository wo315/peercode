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
template <typename V, typename E>
class Graph {
 private:

  class NodeFeature;
  class EdgeFeature;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  using node_value_type = V;
  using edge_value_type = E;

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

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
  Graph(): nodes_(), edges_(), n_idx2uid_(), e_idx2uid_(), adjacency_map(){
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
      Point& position() {
        // graph_ points to memory location, to access its members use ->
        return this -> graph_ -> nodes_[uid_].pos;
      }

      /** Return this node's position. */
      const Point& position() const {
        // graph_ points to memory location, to access its members use ->
        return this -> graph_ -> nodes_[uid_].pos;
      }

      /** Return this node's index, a number in the range [0, graph_size). */
      size_type index() const {
        return this -> graph_ -> nodes_[uid_].n_idx;
      }


      /** function returns a reference to the value of a node, 
       *  so you can edit it
       * @return reference of the value of node
       *
       * Complexity: O(1)
      */
      node_value_type& value(){
        return this -> graph_ -> nodes_[uid_].value;
      }


      /** function returns the value of a node
       * @return value of node
       *
       * Complexity: O(1)
      */
      const node_value_type& value() const{
        return this -> graph_ -> nodes_[uid_].value;
      }

      /** function returns the degree of a node
       * @return degree of node
       *
       * Complexity: O(1)
      */
      size_type degree() const{
        return this -> graph_ -> adjacency_map.at(uid_).size();

      }


      /** function returns iterator to incident nodes
       * @return beginning of IncidentIterator
       *
       * Complexity: O(1)
      */
      incident_iterator edge_begin() const{
        return IncidentIterator(graph_, graph_ -> adjacency_map.at(this -> uid_).begin(), this -> uid_);
      }


      /** function returns iterator to incident nodes
       * @return end of IncidentIterator
       *
       * Complexity: O(1)
      */
      incident_iterator edge_end() const{
        return IncidentIterator(graph_, graph_ -> adjacency_map.at(this -> uid_).end(), this -> uid_);
      }

      /** Test whether this node and @a n are equal.
       *
       * Equal nodes have the same graph and the same index.
       */
      bool operator==(const Node& n) const {
        return ( n.graph_ == this->graph_ and n.uid_ == this->uid_ );
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
        return this-> uid_ < n.uid_;
      }

    private:
      // Use this space to declare private data members and methods for Node
      // that will not be visible to users, but may be useful within Graph.
      // i.e. Graph needs a way to construct valid Node objects

      // Allow Graph to access Node's private member data and functions.
      friend class Graph;

      // Pointer back to the Graph container
      Graph* graph_;
      // This Node's unique identification number
      size_type uid_;

      Node(const Graph* graph, size_type uid)
          : graph_(const_cast<Graph*>(graph)), uid_(uid){
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return n_idx2uid_.size();
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

    size_type new_uid = nodes_.size();
    size_type new_idx = n_idx2uid_.size();

    Node node = Node(this, new_uid);

    NodeFeature feat_node = NodeFeature(this, new_uid, new_idx, position, value, true);

    nodes_.push_back(feat_node);
    n_idx2uid_.push_back(new_uid);

    // initialize entry in map for new node
    adjacency_map[new_uid] = std::map<size_type, size_type>();
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return ( this == n.graph_ );
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert( i < num_nodes());
    return Node(this, n_idx2uid_[i]);     
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

    /** Return a length of an Edge */
    double length() const {
      return norm(this->node1().position() - this->node2().position());
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, uid_n1_);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, uid_n2_);
    }

    /** Return this edge's index, a number in the range [0, num_edges()). */
    size_type index() const {
      return this -> graph_ -> edges_[uid_e_].e_idx;
    }

    /** function returns a reference to the value of an edge, 
     *  so you can edit it
     * @return reference of the value of an edge
     *
     * Complexity: O(1)
    */
    edge_value_type& value(){
      return this->graph_->edges_[uid_e_].value;
    }


    /** function returns the value of an edge
     * @return value of an edge
     *
     * Complexity: O(1)
    */
    const edge_value_type& value() const{
      return this->graph_->edges_[uid_e_].value;
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((node1() == e.node1() and node2() == e.node2())  
        or (node2() == e.node1() and node1() == e.node2()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if  (graph_ != e.graph_){
        return graph_ < e.graph_;
      }
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
      size_type uid_e_, uid_n1_, uid_n2_;

      Edge(const Graph* graph, size_type uid_e, size_type uid_n1, size_type uid_n2)
          : graph_(const_cast<Graph*>(graph)), uid_e_(uid_e), uid_n1_(uid_n1), uid_n2_(uid_n2) {
      }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return e_idx2uid_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    
    size_type uid_e = e_idx2uid_[i];
    EdgeFeature feat_e = edges_[uid_e];

    return Edge(this, uid_e, feat_e.n1_uid, feat_e.n2_uid); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type a_uid = a.uid_;
    size_type b_uid = b.uid_;

    if (adjacency_map.count(a_uid) == 1){

      if (adjacency_map.at(a_uid).count(b_uid) == 1){
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value = edge_value_type()) {
    // see if the edge exists
    size_type a_uid = a.uid_;
    size_type b_uid = b.uid_;

    assert (has_node(a) and has_node(b) and a != b);

    if  (this->has_edge(a,b) == true){

      size_type e_uid = adjacency_map.at(a_uid).at(b_uid);
      return Edge(this, e_uid, a_uid, b_uid);
    }
    else{

      size_type e_uid = edges_.size();
      size_type e_idx = e_idx2uid_.size();

      EdgeFeature feat_e = EdgeFeature(this, a_uid, b_uid, e_uid, e_idx, value, true);
      
      edges_.push_back(feat_e);
      e_idx2uid_.push_back(e_uid);

      // update for edge (a,b) and (b,a)
      adjacency_map[a_uid][b_uid] = e_uid;
      adjacency_map[b_uid][a_uid] = e_uid;

      return Edge(this, e_uid, a_uid, b_uid);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

    std::vector<EdgeFeature>().swap(edges_);
    std::vector<NodeFeature>().swap(nodes_);
    std::vector<size_type>().swap(n_idx2uid_);
    std::vector<size_type>().swap(e_idx2uid_);
    adjacency_map = std::map<size_type, std::map<size_type, size_type>>();

  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {//: private totally_ordered<NodeIterator>{
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

      return Node(this->graph_, *node_iterator_);
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
      return not(this->node_iterator_ == node_iter2.node_iterator_);
    }

   private:
    friend class Graph;
    
    typename std::vector<size_type>::const_iterator node_iterator_;

    Graph* graph_;

    NodeIterator(const Graph* graph, typename std::vector<size_type>::const_iterator node_iterator)
          : node_iterator_(node_iterator), graph_(const_cast<Graph*>(graph)) {
    }

  };

  /**
    * @return iterator for the first node
    *
    * Complexity: O(1)
  */
  node_iterator node_begin() const{
    return NodeIterator(this, this->n_idx2uid_.begin());
  }


  /**
    * @return iterator for the last node
    *
    * Complexity: O(1)
  */
  node_iterator node_end() const{
    return NodeIterator(this, this->n_idx2uid_.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {//: private totally_ordered<IncidentIterator>  {
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
      // map: {node_a_uid : {node_b_uid: edge_ab_uid, ...}, ...}
      return Edge(this->graph_, incid_iterator_->second, node_uid_, incid_iterator_->first);
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
      return not(this->incid_iterator_ == incid_iter2.incid_iterator_);
    }

   private:
    friend class Graph;
    
    // iterator for the map
    typename std::map<size_type, size_type>::const_iterator incid_iterator_;

    size_type node_uid_;

    Graph* graph_;

    IncidentIterator(const Graph* graph, typename std::map<size_type, 
      size_type>::const_iterator incid_iterator, size_type  node_uid)
          : incid_iterator_(incid_iterator), node_uid_(node_uid), graph_(const_cast<Graph*>(graph)) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {//: private totally_ordered<EdgeIterator>{
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
      size_type e_uid = *edge_iterator_;
      return Edge(this->graph_, e_uid, graph_->edges_[e_uid].n1_uid, graph_->edges_[e_uid].n2_uid);
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
      return not(this->edge_iterator_ == edge_iter2.edge_iterator_);
    }

   private:
    friend class Graph;
    // iterator for the map
    typename std::vector<size_type>::const_iterator edge_iterator_;

    Graph* graph_;

    EdgeIterator(const Graph* graph, typename std::vector<size_type>::const_iterator edge_iterator)
          : edge_iterator_(edge_iterator),graph_(const_cast<Graph*>(graph)) {
    }
  };


  /**
    * @return iterator for the first edge
    *
    * Complexity: O(1)
  */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, this->e_idx2uid_.begin());
  }


  /**
    * @return iterator for the last edge
    *
    * Complexity: O(1)
  */
  edge_iterator edge_end() const{
    return EdgeIterator(this, this->e_idx2uid_.end());
  }

  
  /**
    * Remove node from the graph 
    * @return 0 if node does not exist, 1 otw
    * @post num_nodes()-1 
    *
    * Complexity: O(num_nodes() + num_edges())
  */
  size_type remove_node(const Node& n){

    Edge e; 
    if (not has_node(n)){
        return 0;
    }
    
    // remove the edges that are incident to node n
    while (n.edge_begin() != n.edge_end()){
        e = * n.edge_begin();
        remove_edge(e); 
    }

    // n is invalid node
    nodes_[n.uid_].valid = false;
    // swap n with the last node
    n_idx2uid_[n.index()] = n_idx2uid_[n_idx2uid_.size() -  1];
    // update index
    nodes_[n_idx2uid_[n.index()]].n_idx = n.index();
    n_idx2uid_.pop_back();
        
    return 1;
  }

  /**
    * Remove node from the graph 
    * @return node iterator pointing to first valid node
    * @post num_nodes()-1 
    *
    * Complexity: O(num_nodes() + num_edges())
  */
  node_iterator remove_node(node_iterator it){

    remove_node(*it);
    return this->node_begin();
  } 

  /**
    * Remove edge from the graph 
    * @return  0 if edge does not exist, 1 otw
    * @post num_edges()-1 
    *
    * Complexity: O(num_nodes() + num_edges())
  */
  size_type remove_edge(const Node& n1, const Node& n2){

    if (not has_edge(n1, n2)){
        return 0;
    }
    size_type e_uid = (adjacency_map.at(n1.uid_)).at(n2.uid_);
    Edge e = Edge(this, e_uid, n1.uid_, n2.uid_);
    
    // edge e is invalid
    edges_[e_uid].valid = false;
    // delete edge in i2u
    e_idx2uid_[e.index()] = e_idx2uid_[e_idx2uid_.size()-1];
    edges_[e_idx2uid_[e.index()]].e_idx = e.index();
    e_idx2uid_.pop_back();
    
    //update adjacency map
    adjacency_map.at(n1.uid_).erase(n2.uid_);
    adjacency_map.at(n2.uid_).erase(n1.uid_);
    
    return 1;
  } 

  /**
    * Remove edge from the graph 
    * @return  0 if edge does not exist, 1 otw
    * @post num_edges()-1 
    *
    * Complexity: O(num_nodes() + num_edges())
  */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  } 

  /**
    * Remove edge from the graph 
    * @return edge iterator  pointing to start
    * @post num_edges()-1 
    *
    * Complexity: O(num_nodes() + num_edges())
  */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return edge_begin();
  }

 private:

  /**
    Class to store Node information, such as whether is deleted or no, 
    its point, uid, idx, value
  */
  class NodeFeature
  {   
      Graph* graph_;
      size_type n_uid;
      size_type n_idx;
      Point pos;
      node_value_type value;
      bool valid;
      /** Construct an invalid NodeFeature. */
      NodeFeature(){ }

    private:
      friend class Graph;

      NodeFeature(const Graph* graph, size_type n_uid, size_type n_idx, Point pos, 
        node_value_type value, bool valid): 
      graph_(const_cast<Graph*>(graph)), n_uid(n_uid), n_idx(n_idx), 
      pos(pos), value(value), valid(valid) {
      }
  };

  /**
    Class to store Edge information, such as whether is deleted or no, 
    its point, uid, idx, value
  */
  class EdgeFeature
  {   
      Graph* graph_;
      size_type n1_uid;
      size_type n2_uid;
      size_type e_uid;
      size_type e_idx;
      edge_value_type value;
      bool valid;
      /** Construct an invalid EdgeFeature. */
      EdgeFeature(){ }

    private:
      friend class Graph;

      EdgeFeature(const Graph* graph, size_type n1_uid, size_type n2_uid, size_type e_uid, 
        size_type e_idx, edge_value_type value, bool valid): 
      graph_(const_cast<Graph*>(graph)), n1_uid(n1_uid), n2_uid(n2_uid), e_uid(e_uid), 
       e_idx(e_idx), value(value), valid(valid) {
      }
  };


  std::vector<NodeFeature> nodes_;
  std::vector<EdgeFeature> edges_;

  // vectors mapping indices to  uid's of Node / Edge
  // uid is a global identification number, idx is an index among valid nodes  / edges 
  std::vector<size_type> n_idx2uid_;
  std::vector<size_type> e_idx2uid_;

  // map: {node_a_uid : {node_b_uid: edge_ab_uid, ...}, ...}
  std::map<size_type, std::map<size_type, size_type>> adjacency_map;

};

#endif // CME212_GRAPH_HPP
