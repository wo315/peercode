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

template<typename V> 
class Graph {
 private:

  unsigned num_edges_;
  unsigned size_; 

  struct incident_edges;
  struct internal_node;
  struct internal_edge;
  
  std::vector<internal_node> nodes_;
  std::vector<internal_edge> edges_;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of node value. */
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

    /** Return this node's position. */
    const Point& position() const {
        return fetch().position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        // ensure graph is not empty and node id is not bigger than graph size
        //todo
        // assert((n_id_ < graph_->size()) && (this->graph_ != NULL));
        return n_id_;
        // return fetch().n_id_;
    }

    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
        // ensure graph is not empty and node id is not bigger than graph size
        assert((this->graph_ != NULL) && (n_id_ < graph_->size()));
        return fetch().n_val_;
    }
    
    const node_value_type& value() const{
         // ensure graph is not empty and node id is not bigger than graph size
        assert((this->graph_ != NULL) && (n_id_ < graph_->size()));
        return fetch().n_val_;
    }

    size_type degree() const{
        return fetch().adjacency_list_.size();
    }

    incident_iterator edge_begin() const{
        return IncidentIterator(graph_, n_id_, 0);
    }

    incident_iterator edge_end() const{
        return IncidentIterator(graph_, n_id_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((graph_ == n.graph_) && (this->index() == n.index()));
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
      return ((graph_ < n.graph_) || (graph_ == n.graph_ && this->index() < n.index()));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_;  
    size_type n_id_;

    internal_node& fetch() const {
      return graph_->nodes_.at(n_id_);
      assert(false);
    }

    Node(const graph_type* graph, size_type n_id):graph_(const_cast<graph_type*>(graph)), n_id_(n_id){
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
    }

    /** Return a node of this Edge */
    Node node1() const {
        return Node(graph_, n1_id_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
        return Node(graph_, n2_id_);
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

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type e_id_;
    size_type n1_id_;
    size_type n2_id_;

    Edge(const graph_type* graph, size_type e_id, size_type n1_id, size_type n2_id)
    :graph_(const_cast<graph_type*>(graph)), e_id_(e_id), n1_id_(n1_id), n2_id_(n2_id) {
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
    return Edge(this, edges_.at(i).e_id_, edges_.at(i).n1_id_, edges_.at(i).n2_id_);
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

    for (auto & search : nodes_.at(a.index()).adjacency_list_){  
        if (search.adj_n_id_ == b.index()){
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
    // ensure nodes a and b are valid nodes of the graph  
    assert(has_node(a) && has_node(b));
    // ensure node ids aren't equal; otherwise we create a self-loop
    assert(a.index() != b.index());
    
    for (auto & search : nodes_.at(a.index()).adjacency_list_){
        // if edge{a,b} already exists, return that edge object
        if (search.adj_n_id_ == b.index()){
            size_type e_id = search.e_id_;
            return Edge(this, e_id, edges_.at(e_id).n1_id_, edges_.at(e_id).n2_id_);
        }
    }

    size_type e_id = num_edges_;
    // edges vector preserves the order of the root and child nodes order
    edges_.emplace_back(e_id, a.index(), b.index());
    // adjacency_list_ stores node b as adjacent node to node a and vice-versa (graph is undirected)
    nodes_.at(a.index()).adjacency_list_.emplace_back(e_id, b.index());
    nodes_.at(b.index()).adjacency_list_.emplace_back(e_id, a.index());
    num_edges_++;
    return Edge(this, e_id, a.index(), b.index());

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
    edges_.clear();
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
  
  private:

    // Use this space for your Graph class's internals:
    // helper functions, data members, and so forth.
  
    struct incident_edges{
      size_type e_id_; 
      size_type adj_n_id_;
      incident_edges(size_type e_id, size_type adj_n_id): e_id_(e_id), adj_n_id_(adj_n_id){
      }
    };
    
    struct internal_edge{
      size_type e_id_;
      size_type n1_id_;
      size_type n2_id_;
      internal_edge(size_type e_id, size_type n1_id, size_type n2_id)
      : e_id_(e_id), n1_id_(n1_id), n2_id_(n2_id){
      }
    };

    struct internal_node{
        size_type n_id_;
        Point position_;    
        node_value_type n_val_;
        std::vector<incident_edges> adjacency_list_;
        internal_node(size_type n_id, Point position, node_value_type n_val, std::vector<incident_edges> adjacency_list): n_id_(n_id), position_(position), n_val_(n_val), adjacency_list_(adjacency_list){}

    };

};

#endif // CME212_GRAPH_HPP