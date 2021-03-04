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
template <typename V, typename E>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  struct neighbors;
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  // HW1: Modifiable Node Value and Class Templates
  using node_value_type = V;
  using edge_value_type = E;

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
  Graph() : nodes_(), i2u_(), edges_map_() {}

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
      std::cout << "Invalid node: you should use Graph methods" << std::endl;
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().position_;
    }

    Point& position() {
      return fetch().position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_index_;
    }

    /** Return this node's value. */
    node_value_type& value(){
      return fetch().node_value_;
    }

    const node_value_type& value() const{
      return fetch().node_value_;
    }

    /** Return this node's degree, the number of edges connected to that node. */
    size_type degree() const{
      return fetch().adjacency_.size();
    }
 
    IncidentIterator edge_begin() const{
      return IncidentIterator(graph_, node_index_, 0);
    }

    IncidentIterator edge_end() const{
      return IncidentIterator(graph_, node_index_, degree());
    }    

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && index() == n.index());
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
      return ((graph_ < n.graph_) || (graph_ == n.graph_ && index() < n.index()));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    Graph* graph_;  // Pointer back to the Graph
    size_type node_index_;

    Node(const Graph* graph, size_type node_index)
      : graph_(const_cast<Graph*>(graph)), node_index_(node_index){
    }

    internal_node& fetch() const {
      return graph_->nodes_.at(graph_->i2u_.at(node_index_));
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
  Node add_node(const Point& position) {
    return add_node(position, node_value_type());

  }

  /** Add a node to the graph with value node_value, returning the added node. */
  Node add_node(const Point& position, const node_value_type& node_value) {
    // create node
    size_type node_index = num_nodes();
    std::vector<neighbors> adjacency;
    nodes_.emplace_back(node_index, position, node_value, adjacency);
    i2u_.emplace_back(nodes_.size()-1);
    return Node(this, node_index); 
  }

  /**
   * @brief Remove a Node in the graph if it is in the graph
   * 
   * @param n   A Node in the graph
   * @return    1 if the node was removed, 0 otherwise
   * 
   * @pre     _n_ is a valid Node in the graph 
   * @pre     has_node(_n_) == true
   * @post    old _n_ is not a valid Node in the graph
   * @post    has_node(old _n_) == false
   * @post    if return 1: new num_nodes() == old num_nodes() - 1
   *          else: new num_nodes() == old num_nodes()
   *  
   * Total complexity: 
   *    O(maximum degree) for remove_adjacent_edges(const Node&)
   *    O(maximum degree) for replace in adjacency_ for all adjacent nodes
   *    (fetch_adjacency_index(idx1, idx2) in O(1))
   *    O(maximum degree) for replace in edges_map_
   *    Total: O(maximum degree) < O(num_nodes()) at most
   */ 
    size_type remove_node(const Node& n){
    if(!has_node(n)){
      return false;
    }
    else{
      remove_adjacent_edges(n);
      edges_map_.erase(n.index());        

      // replace idx=num_nodes()-1 by idx=n.index() before doing swap/pop
      if(n.index() < num_nodes()-1){
        // adjacency_
        auto adjacency = nodes_.at(i2u_.at(num_nodes()-1)).adjacency_;
        for(size_type i=0; i<adjacency.size(); i++){
          size_type adj_adj_idx_ = fetch_adjacency_index(adjacency.at(i).idx_, num_nodes()-1);
          nodes_.at(i2u_.at(adjacency.at(i).idx_)).adjacency_.at(adj_adj_idx_).idx_ = n.index();
        }        

        // edges_map_
        if(edges_map_.find(num_nodes()-1) != edges_map_.end()){
          auto set = edges_map_.find(num_nodes()-1)->second;
          for(auto elem : set){
            edges_map_.find(elem)->second.erase(num_nodes()-1);
            edges_map_.find(elem)->second.emplace(n.index());
          }
          edges_map_.emplace(n.index(), set);
          edges_map_.erase(num_nodes()-1);   
        }
      }

      // swap/pop
      std::iter_swap(i2u_.begin() + n.index(), i2u_.end()-1);
      i2u_.pop_back();

      return !has_node(n);
    }
  }

  /**
   * @brief Remove the Node in the graph pointer to by a node_iterator
   *        if it is in the graph
   * 
   * @param n   A Node in the graph
   * @return    1 if the node was removed, 0 otherwise
   * 
   * @pre     _n_ is a valid Node in the graph 
   * @pre     has_node(_n_) == true
   * @post    old _n_ is not a valid Node in the graph
   * @post    has_node(old _n_) == false
   * @post    if return 1: new num_nodes() == old num_nodes() - 1
   *          else: new num_nodes() == old num_nodes()
   *  
   * Total complexity: O(num_nodes()) at most
   */ 
  node_iterator remove_node(node_iterator n_it){
    remove(*n_it);
    return n_it;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_ && n.index() < size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
      std::cout << "Invalid edge: you should use Graph methods" << std::endl;
    }

    /** Return this edge's value. 
     * Accessing an edge’s value should be as efficient as possible, no more thanO(d) time 
     * where d is the largest degree of anode but hopefully less.
    */
    edge_value_type& value(){
      size_type adj_index = graph_->fetch_adjacency_index(idx1_, idx2_);
      return graph_->nodes_.at(graph_->i2u_.at(idx1_)).adjacency_.at(adj_index).edge_value_;
    }

    const edge_value_type& value() const{
      size_type adj_index = graph_->fetch_adjacency_index(idx1_, idx2_);
      return graph_->nodes_.at(graph_->i2u_.at(idx1_)).adjacency_.at(adj_index).edge_value_; 
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, idx1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, idx2_);
    }

    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_) \
        && ((idx1_ == e.idx1_ && idx2_ == e.idx2_) || (idx1_ == e.idx2_ && idx2_ == e.idx1_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {     
      if(graph_ < e.graph_) 
        return true;
      else if(graph_ == e.graph_){  
        size_type idx_min = std::min(idx1_, idx2_);
        size_type idx_max = std::max(idx1_, idx2_);
        size_type e_idx_min = std::min(e.idx1_, idx2_);
        size_type e_idx_max = std::max(e.idx1_, idx2_);              
        if (idx_min < e_idx_min || (idx_min == e_idx_min && idx_max < e_idx_max)) 
          return true;
      }      
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Graph* graph_;
    size_type idx1_;
    size_type idx2_;

    Edge(const Graph* graph, size_type idx1, size_type idx2) 
      : graph_(const_cast<Graph*>(graph)), idx1_(idx1), idx2_(idx2) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type num_edges = 0;
    for(auto i : nodes_){
      num_edges += i.adjacency_.size();
    }
    return num_edges/2;    
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return *std::next(edge_begin(), i); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    auto search = edges_map_.find(a.index());
    if(search != edges_map_.end()){
      if(search->second.find(b.index()) != search->second.end()){
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
    return add_edge(a, b, edge_value_type());
  }

  Edge add_edge(const Node& a, const Node& b, const edge_value_type& edge_value) {
    Edge edge(this, a.index(), b.index());

    if(has_edge(a, b) == false){
      edges_map_emplace(a, b);
      edges_map_emplace(b, a);

      nodes_.at(i2u_.at(a.index())).adjacency_.emplace_back(b.index(), edge_value);
      nodes_.at(i2u_.at(b.index())).adjacency_.emplace_back(a.index(), edge_value);
    }
    return edge;

  }

  /**
   * @brief Remove an Edge in the graph if it is in the graph
   * 
   * @param a   A Node in the graph
   * @param b   A Node in the graph
   * @return    1 if the edge was removed, 0 otherwise
   * 
   * @pre     _a_ and _b_ are valid Node in the graph 
   * @post    _a_ and _b_ are valid Node in the graph 
   *          and (new _a_).index() == (old _a_).index(), 
   *              (new _b_).index() == (old _b_).index()
   *          and (new _a_).degree() == (old _a_).degree() - 1, 
   *              (new _b_).degree() == (old _b_).degree() - 1
   * @post    has_edge(Edge(this, new _a_.index(), new _b_.index())) == false
   * @post    if return 1: new num_edges() == old num_edges() - 1
   *          else: new num_edges() == old num_edges()
   *   
   * Total complexity: 
   *    O(num_nodes() + num_edges()) for remove_edge(const Edge)
   *    Total: O(num_nodes() + num_edges())
   */  
  size_type remove_edge(const Node& a, const Node& b){
    return remove_edge(Edge(this, a.index(), b.index()));
  }

  /**
   * @brief Remove an Edge in the graph if it is in the graph
   * 
   * @param edge  An Edge in the graph
   * @return      1 if the edge was removed, 0 otherwise
   * 
   * @pre     _edge_ is a valid Edge in the graph
   * @post    old _edge_ is not a valid Edge in the graph
   * @post    has_edge(old _edge_) == false
   * @post    if return 1: new num_edges() == old num_edges() - 1
   *          else: new num_edges() == old num_edges() 
   * @post    (old _edge_).node1() and (old _edge_).node2() 
   *          are valid Node in the graph
   *          of same index (n.index()) in the graph
   *          but with a degree (n.degree()) reduced by one
   *  
   * Total complexity: 
   *    O(num_edges()) for finding the edge (for loop)
   *    O(num_nodes() + num_edges()) for remove_edge(edge_iterator)
   *    Total: O(num_nodes() + num_edges())
   */
    size_type remove_edge(const Edge& edge){
    if(!has_edge(edge.node1(), edge.node2())){
      return false;
    }
    else{
      for(auto e_it=edge_begin(); e_it!=edge_end(); ++e_it){
        *e_it;
        if(*e_it == edge){
          remove_edge(e_it);
          if(!has_edge(edge.node1(), edge.node2())) return true;
          std::cout << "edge was not deleted" << std::endl;
          throw;
        }
      }
      std::cout << "edge not found in edge_iterator" << std::endl;
      throw;
    }
  }

  /**
   * @brief Remove the Edge in the graph pointed to by an edge_iterator
   *        if it is in the graph
   * 
   * @param e_it  An edge_iterator pointing to an edge in the graph
   * @return      An edge_iterator pointing to an edge in the graph
   * 
   * @pre     _e_it_ is a valid edge_iterator in the graph
   * @post    new _e_it_ is a valid edge_iterator in the graph, 
   *          pointing to another Edge
   * @post    has_edge(*(old _e_it_)) == false
   * @post    if return 1: new num_edges() == old num_edges() - 1
   *          else: new num_edges() == old num_edges()
   * @post    *(old _e_it_).node1() and *(old _e_it_).node2() 
   *          are valid Node in the graph
   *          of same index (n.index()) in the graph
   *          but with a degree (n.degree()) reduced by one 
   * 
   * Total complexity: 
   *    O(num_nodes()) for num_edges()
   *    (Note for later: could be implemented more efficiently)
   *    O(num_edges()) at most for remove_from_adjacency(idx1, idx2)
   *    (O(maximum degree of a node), much less than O(num_edges()) if the graph is sparse)
   *    Total: O(num_nodes() + num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it){
    size_type initial_num_edges = num_edges();

    size_type node1_idx_ = (*e_it).node1().index();
    size_type node2_idx_ = (*e_it).node2().index();

    // remove edge in adjacency of node idx1_ and node idx2_ 
    remove_from_adjacency(node1_idx_, node2_idx_);
    remove_from_adjacency(node2_idx_, node1_idx_);

    // remove edge from edges_map_
    edges_map_.find(node1_idx_)->second.erase(node2_idx_);
    edges_map_.find(node2_idx_)->second.erase(node1_idx_);

    if(num_edges()==initial_num_edges-1){
      return e_it;
    }
    std::cout << "edge was not deleted" << std::endl;
    throw;
  }

  size_type fetch_adjacency_index(size_type idx1, size_type idx2){ 
    for(size_type index=0; index<nodes_.at(i2u_.at(idx1)).adjacency_.size(); index++){
      if(nodes_.at(i2u_.at(idx1)).adjacency_.at(index).idx_==idx2){
        return index;
      }
    }
    std::cout << "idx2 not found in adjacency";
    throw;
  }  

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    i2u_.clear();
    edges_map_.clear();
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
      std::cout << "Invalid NodeIterator" << std::endl;
    }

    Node operator*() const{
      return Node(graph_, idx_);
    }

    NodeIterator& operator++(){
      ++idx_;
      return *this;
    }

    bool operator==(const NodeIterator& node_iter) const {
      return (graph_ == node_iter.graph_) && (idx_ == node_iter.idx_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type idx_;

    NodeIterator(const Graph* graph, size_type idx) :
      graph_(const_cast<Graph*>(graph)), idx_(idx) {}

  };

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
      std::cout << "Invalid IncidentIterator" << std::endl;
    }

    Edge operator*() const{
      size_type adj_node_index = graph_->nodes_.at(graph_->i2u_.at(idx_)).adjacency_.at(adj_index_).idx_;
      return Edge(graph_, idx_, adj_node_index);
    }

    IncidentIterator& operator++(){
      ++adj_index_;
      return *this;
    }

    bool operator==(const IncidentIterator& i_iter) const{
      return (graph_ == i_iter.graph_) && (idx_ == i_iter.idx_)
        && (adj_index_ == i_iter.adj_index_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type idx_;
    size_type adj_index_;

    IncidentIterator(const Graph* graph, size_type idx, size_type adj_index) :
      graph_(const_cast<Graph*>(graph)), idx_(idx), adj_index_(adj_index) {}


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

    EdgeIterator() {
      std::cout << "Invalid EdgeIterator" << std::endl;
    }

    Edge operator*() const{
      return Edge(graph_, idx_, graph_->nodes_.at(graph_->i2u_.at(idx_)).adjacency_.at(adjacency_idx_).idx_);
    }

    EdgeIterator& operator++(){
      if(adjacency_idx_ < adjacency_size_ - 1){
        adjacency_idx_+=1;
        return *this;
      }
      if(adjacency_idx_ == adjacency_size_ - 1){
        idx_ += 1;
        adjacency_idx_ = 0;
        skip_empty_adjacency();
      }
      return *this;
    }

    bool operator==(const EdgeIterator& e_iter) const{
      return (graph_ == e_iter.graph_) && (idx_ == e_iter.idx_)
        && (adjacency_idx_ == e_iter.adjacency_idx_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type idx_;
    size_type adjacency_idx_;
    size_type num_nodes_;
    size_type adjacency_size_;

    /** Skips idx_ for which corresponding adjacency_ size is 0 */
    void skip_empty_adjacency(){
      if(idx_ < num_nodes_){
        adjacency_size_ = graph_->nodes_.at(graph_->i2u_.at(idx_)).adjacency_.size(); 
        while(adjacency_size_==0 && idx_ < num_nodes_ - 1){
          idx_+=1;
          adjacency_size_ = graph_->nodes_.at(graph_->i2u_.at(idx_)).adjacency_.size();
        }
        if(adjacency_size_ == 0){
          idx_ += 1;
          adjacency_size_ = 0;
        }                   
      }
      else{
        adjacency_size_ = 0;
      }
    };

    EdgeIterator(const Graph* graph, size_type idx, size_type adjacency_idx) :
      graph_(const_cast<Graph*>(graph)), idx_(idx), adjacency_idx_(adjacency_idx), num_nodes_(graph->num_nodes()), adjacency_size_(0) {
        skip_empty_adjacency();
      }
    
  };

  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0, 0);
  }

  edge_iterator edge_end() const{
    return EdgeIterator(this, num_nodes(), 0);
  }

  /** Remove idx2_ from adjacency_ of idx1_ */
  void remove_from_adjacency(size_type idx1, size_type idx2){
    size_type adj_index = fetch_adjacency_index(idx1, idx2);
    std::iter_swap(nodes_.at(i2u_.at(idx1)).adjacency_.begin() + adj_index, nodes_.at(i2u_.at(idx1)).adjacency_.end()-1);
    nodes_.at(i2u_.at(idx1)).adjacency_.pop_back();
  }

  void remove_adjacent_edges(const Node& n){
    while(nodes_.at(i2u_.at(n.index())).adjacency_.size() > 0){
      remove_edge(n, Node(this, nodes_.at(i2u_.at(n.index())).adjacency_.at(0).idx_));
    }    
  }

  /** Add edge (a, b) to edges_map_ */
  void edges_map_emplace(const Node& a, const Node& b){
    if(edges_map_.find(a.index()) != edges_map_.end()){
      edges_map_.find(a.index())->second.emplace(b.index());
    }
    else{
      edges_map_.emplace(a.index(), std::initializer_list<size_type>{b.index()});
    }
  }

 private:

  struct internal_node{
    size_type idx_;
    Point position_;
    node_value_type node_value_;
    std::vector<neighbors> adjacency_;
    internal_node(size_type idx, Point position, node_value_type node_value, std::vector<neighbors> adjacency)
      : idx_(idx), position_(position), node_value_(node_value), adjacency_(adjacency) {}
  };

  struct neighbors{
    size_type idx_;
    edge_value_type edge_value_;
    neighbors(size_type idx, edge_value_type edge_value) 
      : idx_(idx), edge_value_(edge_value) {}
  };

  std::vector<internal_node> nodes_;  // stores nodes, node values, adjacencies
  std::vector<size_type> i2u_;  // active nodes
  std::unordered_map<size_type, std::unordered_set<size_type>> edges_map_;  // active edges

};

#endif // CME212_GRAPH_HPP