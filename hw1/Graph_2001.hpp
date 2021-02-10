#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */\
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct internal_node;


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

  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes_(), edges_(), edge_hash_(), size_(0), num_edges_(0) {
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
      // HW0: YOUR CODE HERE

    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[idx_].position;
      //return Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_;
      //return size_type(-1);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /**
     * @brief Return this node's value 
     * @return _value_ 
     * @pre the value in the node is defined
     */
    node_value_type& value(){
      return graph_->nodes_[idx_].value;
    }

    /**
     * @brief Return this node's value as a constant
     * @return _value_ 
     * @pre the value in the node is defined
     */
    const node_value_type& value() const{
      return graph_->nodes_[idx_].value;
    }

    /**
     * @brief Returns the degree of this node
     * @return _degree_ 
     */
    size_type degree() const{
      size_type deg = (graph_->adjacency_map_)[idx_].size();
      return deg;
    }

    /**
     * @brief Returns beginning of incident iterator
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_,idx_,0);
    }

    /**
     * @brief Returns end of incident iterator
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, idx_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_) && (idx_ == n.idx_);
      //(void) n;          // Quiet compiler warning
      //return false;
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
      return (idx_ < n.idx_);
      //(void) n;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type idx_;
    Node(const Graph* graph, size_type idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx){
      }
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_;
    //return 0;
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
  Node add_node(const Point& position, const node_value_type& value_ = node_value_type()) {
    // HW0: YOUR CODE HERE
    nodes_.push_back(internal_node());
    nodes_[size_].position = position;
    nodes_[size_].idx = size_;
    nodes_[size_].value = value_;
    adjacency_map_[size_] = std::vector<size_type>();
    size_ += 1;
    return Node(this, size_ - 1);        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.idx_ < size_);
    //return false;
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_,idx1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_,idx2_);     
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if ((idx1_ == e.idx1_) && (idx2_ == e.idx2_)){
        return true;
      }
      if ((idx1_ == e.idx2_) && (idx2_ == e.idx1_)){
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (edge_idx_ < e.edge_idx_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type idx1_;
    size_type idx2_;
    size_type edge_idx_;
    Edge(const Graph* graph, size_type idx1, size_type idx2, 
    size_type edge_idx)
      : graph_(const_cast<Graph*>(graph)), idx1_(idx1), idx2_(idx2), 
      edge_idx_(edge_idx){
      }
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
    size_type idx1 = std::get<0>(edges_[i]);
    size_type idx2 = std::get<1>(edges_[i]);
    return Edge(this,idx1,idx2,i);       
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type idx1 = a.idx_;
    size_type idx2 = b.idx_;
    // std::pair<size_type,size_type> t = std::make_pair(idx1,idx2);
    // std::pair<size_type,size_type> t2 = std::make_pair(idx2,idx1);
    // if (edge_hash_.count(t)>0){
    //   return true;
    // }
    // if (edge_hash_.count(t2)>0){
    //   return true;
    // }
    // return false;
    // Using new adjacency map
    std::vector<size_type> adj_list = adjacency_map_.at(idx1);
    if (std::find(adj_list.begin(),adj_list.end(),idx2) != adj_list.end()){
      return true;
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
    //(void) a, (void) b;   // Quiet compiler warning
    //Check if current edge exists
    if (has_edge(a,b)){
      std::pair<size_type,size_type> t = std::make_pair(a.idx_,b.idx_);
      std::pair<size_type,size_type> t2 = std::make_pair(b.idx_,a.idx_);
      size_type edge_idx;
      if (edge_hash_.count(t)>0){
        edge_idx = edge_hash_[t];
      }
      else {
        edge_idx = edge_hash_[t2];
      }
      return Edge(this,a.idx_,b.idx_,edge_idx);
      
    }
    //Add edge to vector
    std::pair<size_type,size_type> t = std::make_pair(a.idx_,b.idx_);
    edges_.push_back(t);
    //Add edge to hash map
    edge_hash_[t] = num_edges_;
    num_edges_ += 1;
    //Add edge to adjacency map
    adjacency_map_[a.idx_].push_back(b.idx_);
    adjacency_map_[b.idx_].push_back(a.idx_);
    return Edge(this,a.idx_,b.idx_,num_edges_-1);      
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
    edge_hash_.clear();
    adjacency_map_.clear();
    size_ = 0;
    num_edges_ = 0;
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /**
     * @brief Dereferences the node iterator
     * @pre node index is less than size of graph
     * @return _Node_ object at node index
     */
    Node operator*() const{
      return (*graph).node(node_idx); 
    }

    /**
     * @brief Increments node iterator 
     * 
     * @post node index increases by 1
     */
    NodeIterator& operator++(){
      node_idx += 1;
      return *this;
    }

    /**
     * @brief Equality operator for node iterators
     */
    bool operator==(const NodeIterator& node_it) const {
      return ((graph == node_it.graph) && (node_idx == node_it.node_idx));
    }

    /**
     * @brief Inequality operator for node iterators
     */
    bool operator!=(const NodeIterator& node_it) const {
      return ((graph != node_it.graph) || (node_idx != node_it.node_idx));
    }

  
  
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    NodeIterator(const Graph* graph_, size_type node_idx_){
      this->graph = const_cast<Graph*>(graph_);
      this->node_idx = node_idx_;
    }
    Graph* graph; 
    size_type node_idx;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
   * @brief Beginning of node iterator
   * @return NodeIterator object 
   */
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }

  /**
   * @brief End of node iterator
   * @return NodeIterator object 
   */
  node_iterator node_end() const{
    return NodeIterator(this, this->size_);
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /**
     * @brief Deferences incident iterator
     * @pre the edge index is less than number of edges in graph
     * @return Edge object 
     * 
     * The edge object returned will return an orientation pointing
     * outward from the node that the incident iterator belongs to.
     */
    Edge operator*() const{
      size_type n2 = (*graph).adjacency_map_[n1][edge_idx];
      std::pair<size_type,size_type> t = std::make_pair(n1,n2);
      std::pair<size_type,size_type> t2 = std::make_pair(n2,n1);
      size_type e_idx;
      if ((graph->edge_hash_).find(t) == (graph->edge_hash_).end()){
        e_idx = (graph->edge_hash_)[t2];
      }
      else {
        e_idx = (graph->edge_hash_)[t];
      }
      return Edge(graph,n1,n2,e_idx);
    }

    /**
     * @brief Increments incident operator
     * 
     * @post edge index increases by 1
     */
    IncidentIterator& operator++() {
      edge_idx += 1;
      return *this;
    }

    /**
     * @brief Equality operator 
     */
    bool operator==(const IncidentIterator& inc_it) const {
      return ((graph == inc_it.graph) && (n1 == inc_it.n1)
      && (edge_idx == inc_it.edge_idx));
    }

    /**
     * @brief Inequality operator
     */
    bool operator!=(const IncidentIterator& inc_it) const {
      return ((graph != inc_it.graph) || (n1 != inc_it.n1)
      || (edge_idx != inc_it.edge_idx));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    IncidentIterator(const Graph* graph_, size_type n1_, size_type edge_idx_){
      this->graph = const_cast<Graph*>(graph_);
      this->edge_idx = edge_idx_;
      this->n1 = n1_;
    }
    Graph* graph;
    size_type edge_idx;
    size_type n1;

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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** 
     * @brief Dereferences the edge iterator
     * 
     * Calls the dereference operator for the incident iterator
     * used to implement the edge iterator.
     */
    Edge operator*() const{
      return *incident_iter;
    }

    /**
     * @brief Increments edge iterator
     * 
     * @post the iterator points to the next edge that has not
     * already been seen
     * 
     * Iterates through edges until reaching a new one. Map of visited nodes
     * is used to keep track of which edges have already been touched.
     */
    EdgeIterator& operator++(){
      bool is_visited = true;
      do {
      ++incident_iter;
      if (incident_iter == (*node_iter).edge_end()){
        visited_nodes[(*node_iter).index()] = 1;
        ++node_iter;
        if (node_iter != (*graph).node_end()){
          incident_iter = (*node_iter).edge_begin();
        }
      }
      size_type n2 = (*incident_iter).node2().index();
      if (visited_nodes.find(n2) == visited_nodes.end()){
        is_visited = false;
      }
      }
      while ((node_iter != (*graph).node_end()) && (is_visited));
      return *this;
    }

    /**
     * @brief Equality operator for edge iterator
     */
    bool operator==(const EdgeIterator& edge_it) const {
      if ((graph != edge_it.graph) || (node_iter != edge_it.node_iter)){
        return false;
      }
      if (node_iter == (*graph).node_end()){
        return true;
      }
      else {
        return (incident_iter == edge_it.incident_iter);
      }
    }
    /**
     * @brief Inequality operator for edge iterator
     */
    bool operator!=(const EdgeIterator& edge_it) const {
      return !(*this == edge_it);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    EdgeIterator(const Graph* graph_, node_iterator node_iter_){
      this->graph = const_cast<Graph*>(graph_);
      this->node_iter = node_iter_;
      this->visited_nodes = std::unordered_map<size_type,size_type>();
      if (node_iter != (*graph).node_end()){
        this->incident_iter = (*node_iter_).edge_begin();
      }
    }
    Graph* graph;
    std::unordered_map<size_type,size_type> visited_nodes; 
    node_iterator node_iter;
    incident_iterator incident_iter; 

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
   * @brief Beginning of edge iterator
   * @return EdgeIterator object
   * 
   * Sets the node iterator in the beginning edge iterator to 
   * the beginning node iterator of the _graph_
   */
  edge_iterator edge_begin() const{
    node_iterator node_it_begin = node_begin();
    return EdgeIterator(this, node_it_begin);
  }

  /**
   * @brief End of edge iterator
   * @return EdgeIterator object
   * 
   * Sets the node iterator in the end edge iterator to 
   * the end node iterator of the _graph_
   */
  edge_iterator edge_end() const{
    node_iterator node_it_end = node_end();
    return EdgeIterator(this, node_it_end);
  }
  

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct internal_node {
    Point position;
    size_type idx;
    node_value_type value;
  };

  using pair = std::pair<size_type,size_type>;
  std::vector<internal_node> nodes_;
  std::vector<pair> edges_;
  std::unordered_map<pair,size_type,boost::hash<pair>> edge_hash_;
  std::unordered_map<size_type, std::vector<size_type>> adjacency_map_;
  size_type size_; //Size of graph
  size_type num_edges_;

};

#endif // CME212_GRAPH_HPP
