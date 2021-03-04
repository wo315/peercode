#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template<typename V, typename E> 
class Graph {
 private:

  unsigned num_edges_;
  unsigned size_; 

  struct incident_edges;
  struct internal_node;
  struct internal_edge;
  
  std::vector<internal_node> nodes_; // indexed by node uid 
  std::vector<unsigned int> nodes_i2u_; // indexed by node id 

  std::vector<internal_edge> edges_; // indexed by edge uid
  std::vector<unsigned int> edges_i2u_; // indexed by edge id 

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of node value. */
  using node_value_type = V;

  /** Type of edge value. */
  using edge_value_type = E;

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
  {
    num_edges_ = 0;
    size_ = 0;
    nodes_ = std::vector<internal_node>();
    edges_ = std::vector<internal_edge>();
  }

  /** Ensure constructed graph is cleared. */
  ~Graph(){
      clear();
  }

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

    /** Return reference to this node's position so we can modify it. */
    Point& position() {
      return fetch().position_;
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        // ensure graph is not empty and node id is not bigger than graph size
        return fetch().n_id_;
    }

    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
        // ensure graph is not empty and node id is not bigger than graph size
        assert((this->graph_ != NULL) && (index() < graph_->size()));
        return fetch().n_val_;
    }
    
    const node_value_type& value() const{
         // ensure graph is not empty and node id is not bigger than graph size
        assert((this->graph_ != NULL) && (index() < graph_->size()));
        return fetch().n_val_;
    }

    size_type degree() const{
        return fetch().adjacency_list_.size();
    }

    incident_iterator edge_begin() const{
        return IncidentIterator(graph_, n_uid_, 0);
    }

    incident_iterator edge_end() const{
        return IncidentIterator(graph_, n_uid_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((graph_ == n.graph_) && (n_uid_ == n.n_uid_));
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
      return ((graph_ < n.graph_) || (graph_ == n.graph_ && n_uid_ < n.n_uid_));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_;  
    size_type n_uid_;

    internal_node& fetch() const {
      return graph_->nodes_.at(n_uid_);
    }

    Node(const graph_type* graph, size_type n_uid):graph_(const_cast<graph_type*>(graph)), n_uid_(n_uid){
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
      return size_;
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
  Node add_node(const Point& position, const node_value_type& nodeval =node_value_type()) {
    std::vector<incident_edges> adj_list{};
    nodes_.emplace_back(size_, position, nodeval, adj_list);
    nodes_i2u_.push_back(nodes_.size()-1);
    size_++;
    return Node(this, nodes_.size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      return ((this == n.graph_) && (n.index() < size()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // ensure node id is less than number of nodes in the graph (we know i >= 0 because it's an unsigned int)
    assert(i < num_nodes());
    return Node(this, nodes_i2u_[i]);
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
        return Node(graph_, n1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
        return Node(graph_, n2_uid_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        return ((graph_ == e.graph_) && 
      ((node1().index() == e.node1().index() && node2().index() == e.node2().index())
      || (node1().index() == e.node2().index() && node2().index() == e.node1().index())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        return ((graph_ < e.graph_) || ((graph_ == e.graph_) && ((node1().index() < e.node1().index()) 
        || ((node1().index() == e.node1().index()) && (node2().index() < e.node2().index())))) );
    }

    double length() const{
      Point n1_pos = node1().position();
      Point n2_pos = node2().position();
      return norm(n1_pos - n2_pos);
    }

    edge_value_type& value(){
      return graph_->edges_.at(e_uid_).e_val_;
    }

    const edge_value_type& value() const{
      return graph_->edges_.at(e_uid_).e_val_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type e_uid_;
    size_type n1_uid_;
    size_type n2_uid_;

    Edge(const graph_type* graph, size_type e_uid, size_type n1_uid, size_type n2_uid)
    :graph_(const_cast<graph_type*>(graph)), e_uid_(e_uid), n1_uid_(n1_uid), n2_uid_(n2_uid) {
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
    // ensure edge id is less than number of edges in the graph (we know i >= 0 because it's an unsigned int)
    assert(i < num_edges());
    size_type uid = edges_i2u_[i];
    return Edge(this, uid, edges_.at(uid).n1_id_, edges_.at(uid).n2_id_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // ensure nodes a and b are valid nodes of the graph  
    assert((this == a.graph_) && (this == b.graph_));
    // ensure both nodes have index less than the number of nodes in the graph
    assert((a.index() < size()) && (b.index() < size()));

    for (auto & search : nodes_.at(a.n_uid_).adjacency_list_){  
        if (search.adj_n_id_ == b.n_uid_){
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& e_val = edge_value_type()) {
    // ensure nodes a and b are valid nodes of the graph  
    assert(has_node(a) && has_node(b));
    // ensure node ids aren't equal; otherwise we create a self-loop
    assert(a.n_uid_ != b.n_uid_);
    
    size_type aid = a.n_uid_;
    size_type bid = b.n_uid_;
    for (auto & search : nodes_.at(aid).adjacency_list_){
        // if edge{a,b} already exists, return that edge object
        if (search.adj_n_id_ == bid){
            size_type e_id = search.e_id_;
            return Edge(this, e_id, edges_.at(e_id).n1_id_, edges_.at(e_id).n2_id_);
        }
    }
    // edges vector preserves the order of the root and child nodes order
    edges_.emplace_back(num_edges_, aid, bid, e_val);
    edges_i2u_.push_back(edges_.size()-1);
    // adjacency_list_ stores node b as adjacent node to node a and vice-versa (graph is undirected)
    nodes_.at(aid).adjacency_list_.emplace_back(edges_.size()-1, bid);
    nodes_.at(bid).adjacency_list_.emplace_back(edges_.size()-1, aid);
    num_edges_++;
    return Edge(this, edges_.size()-1, aid, bid);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    num_edges_ = 0;
    size_ = 0;
    nodes_.clear();
    nodes_i2u_.clear();
    edges_.clear();
    edges_i2u_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable <NodeIterator> {
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

    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const{
        return Node(graph_, n_id_);
    }

    NodeIterator& operator++(){
        n_id_++;
        return *this;
    }

    bool operator==(const NodeIterator& node_iter) const{
        return ((graph_==node_iter.graph_) && (n_id_==node_iter.n_id_)); 
    }


   private:
    friend class Graph;
    graph_type* graph_;
    size_type n_id_;
    NodeIterator(const graph_type* graph, size_type n_id): 
    graph_(const_cast<graph_type*>(graph)), n_id_(n_id){  
    }
    
  };

  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const{
      return NodeIterator(this, 0);
  }

  node_iterator node_end() const{
    return NodeIterator(this, size());
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

    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
        size_type e_id = graph_->nodes_.at(n_id_).adjacency_list_.at(inc_id_).e_id_;
        size_type adj_n_id = graph_->nodes_.at(n_id_).adjacency_list_.at(inc_id_).adj_n_id_;
        return Edge(graph_, e_id, n_id_, adj_n_id);
    }

    IncidentIterator& operator++(){
        inc_id_++;
        return *this;
    }

    bool operator==(const IncidentIterator& ii) const{
        return ((graph_ == ii.graph_) && (n_id_ == ii.n_id_) && (inc_id_ == ii.inc_id_));
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type n_id_;
    size_type inc_id_;
    
    IncidentIterator(const graph_type* graph, size_type n_id, size_type inc_id): 
    graph_(const_cast<graph_type*>(graph)), n_id_(n_id), inc_id_(inc_id){
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

    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
        return graph_->edge(e_id_);
    }

    EdgeIterator& operator++(){
        e_id_++;
        return *this;
    }

    bool operator==(const EdgeIterator& ei) const{
        return ((graph_ == ei.graph_) && (e_id_==ei.e_id_));
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type e_id_;

    EdgeIterator(const graph_type* graph, size_type e_id): 
    graph_(const_cast<graph_type*>(graph)), e_id_(e_id) {
    }
    
  };

  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const{
      return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const{
      return EdgeIterator(this, this->num_edges());
  }

  /** Remove an edge from the graph, if edge exists, using the edge endpoints @a a and @a b. 
   * Return 1 if one of the edges is removed, otherwise return 0 if no edge is removed.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return 1 if @pre has_edge(@a a, @a b)== true, 0 otherwise
   * @post new has_edge(@a a, @a b) == false
   * @post If new has_edge(@a a, @a b) == false, new num_edges() == old num_edges() - 1.
   *       Else,                                 new num_edges() == old num_edges().
   *
   * Invalidates: edge(graph, eid, aid, bid)
   *              a.adjacency_list(eid, bid) 
   *              b.adjacency_list(eid, aid)
   *              edge_iterators representing the now invalid edge (a,b)
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Node&a, const Node&b){

    if (has_edge(a,b)){
      size_type aid = a.n_uid_;
      size_type bid = b.n_uid_;
      for (size_type i= 0; i<nodes_.at(aid).adjacency_list_.size(); ++i){
        if (nodes_.at(aid).adjacency_list_.at(i).adj_n_id_ == bid){
          size_type e_uid = nodes_.at(aid).adjacency_list_.at(i).e_id_;
          size_type e_id = edges_.at(e_uid).e_id_;
          nodes_.at(aid).adjacency_list_.erase(nodes_.at(aid).adjacency_list_.begin()+i);

          for (size_type j=0; j < nodes_.at(bid).adjacency_list_.size(); ++j){
            if (nodes_.at(bid).adjacency_list_.at(j).adj_n_id_ == aid){
              nodes_.at(bid).adjacency_list_.erase(nodes_.at(bid).adjacency_list_.begin()+j);
            }
          }
          edges_i2u_.erase(edges_i2u_.begin() + e_id);
          for (size_type id = e_id; id < edges_i2u_.size(); ++id){
            edges_.at(edges_i2u_.at(id)).e_id_ = id;
          }
          num_edges_ -= 1;
          return 1;
        }
      }
    }
  return 0;
  }

  /** Remove an edge from the graph if edge exists. 
   * Return 1 if one of the edges is removed, otherwise return 0 if no edge is removed.
   * @pre @a e is a valid edge of this graph
   * @return 1 if @pre has_edge(@a e.node1(), @a e.node2())== true, 0 otherwise
   * @post new has_edge(@a e.node1(), @a e.node2()) == false
   * @post If new has_edge(@a a, @a b) == false, new num_edges() == old num_edges() - 1.
   *       Else,                                 new num_edges() == old num_edges().
   *
   * Invalidates: edge(graph, eid, aid, bid)
   *              a.adjacency_list(eid, bid) 
   *              b.adjacency_list(eid, aid)
   *              edge_iterators representing the now invalid edge (a,b)
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Edge& e){
    Node n1 = e.node1();
    Node n2 = e.node2();
    return remove_edge(n1, n2);
  }

  /** Remove an edge from the graph using an edge iterator. 
   * @pre @a e_it is a valid edge iterator of this graph.
   * @return the next edge iterator for the edge following the removed one
   * @post new has_edge( @a e.node1(), @a e.node2()) == false
   * @post If new has_edge( @a a, @a b) == false, new num_edges() == old num_edges() - 1.
   *       Else,                                 new num_edges() == old num_edges().
   *
   * Invalidates: edge(graph, eid, aid, bid)
   *              a.adjacency_list(eid, bid) 
   *              b.adjacency_list(eid, aid)
   *              edge_iterators representing the now invalid edge (a,b)
   * Complexity: No more than O(num_nodes() + num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it){
    auto e = *e_it;
    remove_edge(e);
    edge_iterator next_eit = e_it;
    return e_it;
  }

  /** Remove a node from the graph.
  * Return 1 if one of the nodes is removed, otherwise return 0 if no node is removed.
  * @pre @a n is a node 
  * @return @a 1 if n is a valid node of this graph, otherwise 0 if not a valid node of this graph.
  * @post new has_node( @a n) == false
  * @post If @a n is removed, new num_nodes() == old num_nodes()-1.
  *       Else,               new num_nodes() == old num_nodess().
  *
  * Invalidates: node(g, n.index())
   *             adjacency list of node @a n, nodes_.at(n.index()).adjacency_list_
   *             also adjacency list of @a n's adjacent nodes
   *             edges with n == e.node1() or n == e.node2()
   *             node_iterators representing the now invalid node @a n
   *             edge_iterators representing the now invalid edge between @a n and adjacent nodes
  * Complexity: No more than O(num_nodes()).
  */
  size_type remove_node(const Node & n){
    if (has_node(n)){
      nodes_i2u_.erase(nodes_i2u_.begin() + n.index());
      for (size_type i = n.index(); i < nodes_i2u_.size(); ++i){
        nodes_.at(nodes_i2u_[i]).n_id_ = i;
      }
		  while (n.degree() > 0) {
			  auto e = *n.edge_begin();
			  remove_edge(e);
      }
      size_ -= 1;
      return 1;
    }
    return 0;
  }


  /** Remove a node from the graph using a node iterator.
  * @pre @a n_it is a node iterator
  * @return the next node iterator for the node following the removed one.
  * @post @a n_it is not a valid node iterator
  * @post If @a n is removed, new num_nodes() == old num_nodes() - 1.
  *       Else,               new num_nodes() == old num_nodess().
  * Invalidates: node(g, n.index())
   *             adjacency list of node @a n, nodes_.at(n.index()).adjacency_list_
   *             also adjacency list of @a n's adjacent nodes
   *             edges with n == e.node1() or n == e.node2()
   *             node_iterators representing the now invalid node @a n
   *             edge_iterators representing the now invalid edge between @a n and adjacent nodes
  * Complexity: No more than O(num_nodes()).
  */
  node_iterator remove_node(node_iterator n_it){
    auto n = *n_it;
    remove_node(n);
    return n_it;
  }

  private:

    // Use this space for your Graph class's internals:
    // helper functions, data members, and so forth.
  
    struct incident_edges{
      size_type e_id_; 
      size_type adj_n_id_;
      incident_edges(size_type e_id, size_type adj_n_id): 
      e_id_(e_id), adj_n_id_(adj_n_id){
      }
    };
    
    struct internal_edge{
      size_type e_id_;
      size_type n1_id_;
      size_type n2_id_;
      edge_value_type e_val_;
      internal_edge(size_type e_id, size_type n1_id, size_type n2_id, edge_value_type e_val)
      : e_id_(e_id), n1_id_(n1_id), n2_id_(n2_id), e_val_(e_val){
      }
    };

    struct internal_node{
        size_type n_id_;
        Point position_;    
        node_value_type n_val_;
        std::vector<incident_edges> adjacency_list_;
        internal_node(size_type n_id, Point position, node_value_type n_val, std::vector<incident_edges> adjacency_list): 
        n_id_(n_id), position_(position), n_val_(n_val), adjacency_list_(adjacency_list){
        }
    };

};

#endif // CME212_GRAPH_HPP