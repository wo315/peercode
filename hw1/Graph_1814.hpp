#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <set>

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

  std::vector <std::pair<Point,V>> nodes_; //vector of points/nodes
  std::unordered_map <unsigned int, std::pair<unsigned int, unsigned int>> edges_; //maps of edges with indexes as keys
  std::unordered_map <unsigned int, std::unordered_map<unsigned int, unsigned int>> n2e_; //maps of nodes to connected nodes and edge index

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
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

    //Construct an invalid node
    Node(): graph_(nullptr), nid_(0) {
    }

    /** Return this node's position. */
    const Point& position() const {
      return this->graph_->nodes_[nid_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(this->graph_!=nullptr);
      return nid_;
    }

    /** Retuen the value of this node */
    node_value_type& value(){
      return this->graph_->nodes_[nid_].second;
    };
    const node_value_type& value() const{
      return this->graph_->nodes_[nid_].second;
    };

    /* Return the number incidents/edges of a node*/
    size_type degree() const{
      return this->graph_->n2e_[this->nid_].size();
    };

    /* begin and end for inciedent iterator */
    incident_iterator edge_begin() const{
      std::unordered_map<size_type, size_type>::const_iterator cit=this->graph_->n2e_[this->nid_].begin();
      return  IncidentIterator(this->graph_, *this, cit);
    };

    incident_iterator edge_end() const{
      std::unordered_map<size_type, size_type>::const_iterator cit=this->graph_->n2e_[this->nid_].end();
      return  IncidentIterator(this->graph_, *this, cit);
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      assert(n.graph_!=nullptr and this->graph_!=nullptr);
      return (this->graph_== n.graph_ and this->nid_==n.nid_);
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
      assert(n.graph_!=nullptr and this->graph_!=nullptr);
      if (graph_ < n.graph_){
        return true;
      }else if (graph_==n.graph_ and nid_<n.nid_){
        return true;
      }else{  
      return false;}
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_; //pointer to graph the node belongs to
    size_type nid_; //node index in the vector
    Node(const graph_type* graph, size_type nid)
        : graph_(const_cast<graph_type*>(graph)), nid_(nid) {
    }
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
    return nodes_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value= node_value_type()) {
    nodes_.push_back(std::pair<Point,node_value_type>(position,value));
    return Node(this, nodes_.size()-1);      
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (this==n.graph_ and n.nid_<size()){
      return true;
    } else {  
      return false;}
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < nodes_.size()) {
      return Node(this,i);
    }else{         
      return Node();}       // Invalid node
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
    Edge(): graph_(nullptr), n1id_(0), n2id_(0), eid_(0){
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(this->graph_,n1id_);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(this->graph_,n2id_);    
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      assert(e.graph_!=nullptr and this->graph_!=nullptr);
      if (graph_== e.graph_ and n1id_==e.n1id_ and n2id_==e.n2id_ ){
        return true;
      }else if(graph_== e.graph_ and n1id_==e.n2id_ and n2id_==e.n1id_){
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
    bool operator<(const Edge& e) const {
      if (graph_ < e.graph_){
        return true;
      }else if (graph_==e.graph_ and eid_<e.eid_){
        return true;
      }else{
        return false;}
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_; //graph where it belongs to
    size_type n1id_; //first  node
    size_type n2id_; //second node
    size_type eid_; //ID of the edge
    Edge(const graph_type* graph, size_type n1id, size_type n2id, size_type eid)
        : graph_(const_cast<graph_type*>(graph)), n1id_(n1id), n2id_(n2id), eid_(eid) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return this->edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i<num_edges()){
      return Edge(this, edges_.at(i).first, edges_.at(i).second,i);
    }else{
      return Edge();
    }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type aid=a.nid_;
    size_type bid=b.nid_;
    auto it = n2e_.find(aid);
    if (it != n2e_.end()){
      return (n2e_.at(aid).count(bid)!=0);
    }else{
      return false;}
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
    unsigned int ida=a.nid_;
    unsigned int idb=b.nid_;
    if (has_edge(a, b)){
      unsigned int eidx=n2e_.at(ida).at(idb); //find the index of the edge in the nested map
      return Edge(this, ida, idb, eidx); //return the current edge if already exists
    }else{
      unsigned int eidx=edges_.size();
      edges_.insert({eidx,std::pair<unsigned int,unsigned int>(ida,idb)});//add new edge to the edges map
      if(n2e_.count(ida) == 0) {
        std::unordered_map<unsigned int, unsigned int> nidmap1 = {{idb, eidx}};
        n2e_.insert({ida,nidmap1});
      }
      if(n2e_.count(idb) == 0){
        std::unordered_map<unsigned int, unsigned int> nidmap2 = {{ida, eidx}};
        n2e_.insert({idb,nidmap2}); 
      }
      n2e_.at(ida).insert({idb, eidx});
      n2e_.at(idb).insert({ida, eidx});
      return Edge(this, ida, idb, eidx); 
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    n2e_.clear();
  }

  //
  // Node Iterator
  //
 
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:public std::iterator<std::forward_iterator_tag,Node> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      graph_=nullptr;
      currentIdx_=0;
    }

    NodeIterator& operator++(){
      currentIdx_++;
      return (*this);
    }

    bool operator==(const NodeIterator& node_iter)const{
      return(graph_==node_iter.graph_ and currentIdx_==node_iter.currentIdx_);
    }

    bool operator!=(const NodeIterator& node_iter)const{
      return !(node_iter == (*this));
    }

    Node operator*() const {
      return Node(graph_,currentIdx_);
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type currentIdx_;

    NodeIterator(const graph_type* graph, size_type currentIdx){
      this->graph_=graph;
      this->currentIdx_=currentIdx;
    }
  };

  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  NodeIterator node_end() const {
    return NodeIterator(this, nodes_.size());
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
      graph_=nullptr;
      node_=NULL;
      cit_=NULL;
    }

    Edge operator*() const{
      return Edge(this->graph_,node_.nid_,(*cit_).first,(*cit_).second);
    }
    
    IncidentIterator& operator++(){
      cit_++;
      return (*this);
    }
    
    bool operator==(const IncidentIterator& inc_iter) const{
      return (this->graph_==inc_iter.graph_ and this->cit_==inc_iter.cit_);
    }

    bool operator!=(const IncidentIterator& inc_iter) const{
      return !(this->cit_==inc_iter.cit_);
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    node_type node_;
    std::unordered_map<size_type, size_type>::const_iterator cit_;

    IncidentIterator(const graph_type* graph, node_type node,std::unordered_map<size_type, size_type>::const_iterator cit){
      this->graph_=graph;
      this->node_=node;
      this->cit_=cit;
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
    EdgeIterator(){
      graph_=nullptr;
      currentIdx_=0;
    }

    Edge operator*() const{
      size_type n1id=graph_->edges_.at(currentIdx_).first;
      size_type n2id=graph_->edges_.at(currentIdx_).second;
      return Edge(this->graph_,n1id,n2id,currentIdx_);
    }
    EdgeIterator& operator++(){
      currentIdx_++;
      return (*this);
    }
    bool operator==(const EdgeIterator& edge_iter) const{
      return(edge_iter.graph_==this->graph_ and edge_iter.currentIdx_==this->currentIdx_);
    }

    bool operator!=(const EdgeIterator& edge_iter) const{
      return !(edge_iter == (*this));
    }

   private:
    friend class Graph;
    const graph_type* graph_;
    size_type currentIdx_;

    EdgeIterator(const graph_type* graph, size_type currentIdx){
      this->graph_=graph;
      this->currentIdx_=currentIdx;
    }
  };

  edge_iterator edge_begin() const{
    return EdgeIterator(this,0);
  }
  edge_iterator edge_end() const{
    return EdgeIterator(this,edges_.size());
  }

 private:


};

#endif // CME212_GRAPH_HPP
