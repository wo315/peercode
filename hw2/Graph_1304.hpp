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

/**
 * My feedback from HW0 indicated that my implementation of a proxy
 * design pattern was inefficient at best and did not meet the overall
 * intent of the assignmnet.
 * So after reviewing some of the peercode, I adopted a different design
 * choice; one that would ensure my Node class would act as an
 * interface with the content of the Graph. In particular I would like to
 * reference Graph_1495.hpp and Graph_163.hpp. These files gave me considerable
 * inspiration for my Graph class.
 */

template <typename V, typename E>
class Graph {
private:
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct EdgeInfo;

public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of the node value */
  typedef V node_value_type;

  /** Type of the edge value*/
  typedef E edge_value_type;

  /** Type of the node values */
  //using node_value_type = V;
  
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

  
private:
  // Use this space for your Graph class's internals:
  // helper functions, data members, and so forth.

 /**
  * Vector that stores pairs, consisting of
  * position and node values. Here node _i_
  * is found at index _i_ within the vector _nodes_
  */
  std::vector<std::pair<Point, V>> nodes;
  
  /**
   * Vector that stores pairs of node indices for each edge.
   * Here the pair of nodes indices at index _i_ forms edge _i_
   */
  std::vector<std::pair<size_type, size_type>> edges;

  /**
   * This Vector of vectors makes up our adjacency list.
   * Within each inner vector we are storing pairs which consist of
   * a node_idx and an edge_idx. This is done such that node _i_
   * can be found when the follwing code is executed: adjacency[node_i_].
   * Here, the corresponding inner vector result gives us all the
   * nodes that are adjacent to node_i_ and the edge index that
   * connects node_i_ with the adjacent node.
   */
  std::vector<std::vector<std::pair<size_type, size_type>>> adjacency;

  //Decided to use Struct just in case we need to add more information
  //in future assignments.
  struct EdgeInfo {
    edge_value_type edge_val;

    EdgeInfo(edge_value_type edge_val_) :
      edge_val(edge_val_) {}
  };

  std::vector<EdgeInfo> edge_vec;
  
public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  // Graphs's private attributes initialized to 0
  Graph() : nodes(std::vector<std::pair<Point, V>>(0)),
            edges(std::vector<std::pair<size_type, size_type>>(0)),
            adjacency(std::vector<std::vector<std::pair<size_type, size_type>>>(0)) {
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
      // if parent_graph is a nullptr then node is invalid
      this->parent_graph = nullptr;
    }

    /** Return a reference to this node's position which can be modified*/
    Point& position(){
      return std::get<0>(parent_graph->nodes[node_idx]);
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(this->parent_graph != NULL);
      return std::get<0>(parent_graph->nodes[node_idx]);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(parent_graph->size() != 0 && parent_graph->size() > node_idx);
      return this->node_idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /**
     * @return The node_value_type of this node.
     */
    node_value_type& value(){
      return std::get<1>(parent_graph->nodes[node_idx]);
    }

    /**
     * @return The const node_value_type of this node
     */
    const node_value_type& value() const{
      return std::get<1>(parent_graph->nodes[node_idx]);
    }

    /**
     * @return The number of nodes adjacent to this node
     */
    size_type degree() const{
      return parent_graph->adjacency[node_idx].size();
    }

    /**
     * @return An IncidentIterator for the current node.
     *         The IncidentIterator points to the first
     *         edge adjacent to this node.
     *
     * @post The inicident_iterator must be initialized
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(parent_graph, this->node_idx, 0);
    }

    /**
      * @return An IncidentIterator for the current node.
      *         The IncidentIterator points to the last
      *         edge adjacent to this node.
      */
    incident_iterator edge_end() const{
      return IncidentIterator(parent_graph, this->node_idx, this->degree());
    }

    /**
     * @brief Test whether this node and @a n are equal.
     *
     * @param[in] _n_  node pf this graph
     * @pre _n_ is a valid node of this graph
     * @return true if _n_'s index is smaller than current node's index
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (parent_graph==n.parent_graph && node_idx == n.node_idx)
        return true;
      else{
        return false;
      }
    }

    /**
     * @brief Test whether this node is less than @a n in a global order.
     *
     * @param[in] _n_ node of this graph
     * @pre _n_ is a valid noe of this graph
     * @return true if _n_ has smaller index than current node
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      assert(n.parent_graph != NULL);
      if(parent_graph < n.parent_graph)
        return true;
      if(parent_graph==n.parent_graph && node_idx < n.node_idx){
        return true;
      }
      else {
        return false;
      }
    }

  private:
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    Graph* parent_graph;
    size_type node_idx;

    /**
     * @brief Private constructor for node object
     *
     * @param[in] _parent_graph_ptr_  A pointer to the current Graph object
     * @param[in] _idx_  The index of a node in the Graph
     * @return A valid Node from the parent Graph class.
     *
     * @pre _parent_graph_ptr_ must not point to a nullptr
     * @pre _idX_ must me a valid node_idx
     */
    Node(const Graph* parent_graph_ptr, size_type idx) :
      parent_graph(const_cast<Graph*>(parent_graph_ptr)),
      node_idx(idx){
  }
      
};

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] next_val The next value for the new node
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    std::pair<Point, node_value_type> new_node(position, value);
    nodes.emplace_back(new_node);
    std::vector<std::pair<size_type, size_type>> new_adj;
    adjacency.emplace_back(new_adj);
    
    return Node(this, nodes.size() - 1);
   }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.parent_graph && n.node_idx < nodes.size());
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(parent_graph, node1_idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(parent_graph, node2_idx);
    }

    /**
     * @return The length of an edge as the Euclidean distance between
     *         the 2 points defining the edge
     */
    double length() const{
      Point p1 = node1().position();
      Point p2 = node2().position();
      return norm(p1 - p2);
    }

    edge_value_type& value() {
      return parent_graph->edge_vec[edge_idx].edge_val;
    }

    const edge_value_type& value() const {
      return parent_graph->edge_vec[edge_idx].edge_val;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return((parent_graph == e.parent_graph) && (node1() == e.node1() && node2() == e.node2()) || (node1() == e.node2() && node2() == e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return((parent_graph < e.parent_graph) || (parent_graph == e.parent_graph && (edge_idx < e.edge_idx)));
    }

   private:
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // Allow Graph to access Edge's private member data and functions.
    
    friend class Graph;
    
    Graph* parent_graph;
    size_type edge_idx;
    size_type node1_idx;
    size_type node2_idx;

     /**
     * @brief Private constructor for edge object
     *
     * @param[in] _parent_graph_ptr_  A pointer to the current Graph object
     * @param[in] _edge_idx_  A node index for the current graph object
     * @param[in] _node1_  One of the node objects required to form an edge
     * @param[in] _node2_  THe other node object required to form an edge
     * @return A valid NodeEdge object from the parent Graph class.
     *
     * @pre _parent_graph_ptr_ must not point to a nullptr
     * @pre _edge_idx_ must me a valid edge index
     * @pre _node1_ must me a valid node object
     * @pre _node2_ must me a valid node object
     */
    Edge(const Graph* parent_graph_ptr, size_type edge_idx_, size_type node1, size_type node2) :
      parent_graph(const_cast<Graph*>(parent_graph_ptr)),
      edge_idx(edge_idx_),
      node1_idx(node1),
      node2_idx(node2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if(i < num_edges())
      return Edge(this, i, std::get<0>(edges[i]), std::get<1>(edges[i]));
    else
      return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if(this == a.parent_graph && this == b.parent_graph) {
      size_type a_idx = a.node_idx;
      size_type b_idx = b.node_idx;

      for(unsigned int i = 0; i < a.degree(); ++i)
        if(b_idx == std::get<0>(adjacency[a_idx][i]))
	  return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type eval = edge_value_type()) {
    assert((this == a.parent_graph && this == b.parent_graph) && (!(a == b)));
    size_type a_idx = a.node_idx;
    size_type b_idx = b.node_idx;

    for(unsigned int i = 0; i < adjacency[a_idx].size(); ++i)
      if(b_idx == std::get<0>(adjacency[a_idx][i]))
	return Edge(this, std::get<1>(adjacency[a_idx][i]), a_idx, b_idx);

    //Build a new edge otherwise
    std::pair<size_type, size_type> node_idx_pair;

    if(a < b){
      std::get<0>(node_idx_pair) = a_idx;
      std::get<1>(node_idx_pair) = b_idx;
      edges.emplace_back(node_idx_pair);
      edge_vec.emplace_back(eval);
    }
    else {
      std::get<0>(node_idx_pair) = b_idx;
      std::get<1>(node_idx_pair) = a_idx;
      edges.emplace_back(node_idx_pair);
      edge_vec.emplace_back(eval);
    }

    // Update adjacency list
    // New edge_idx will be edges.size()-1
    std::pair<size_type, size_type> b_pair(b_idx, edges.size() - 1);
    std::pair<size_type, size_type> a_pair(a_idx, edges.size() - 1);

    adjacency[a_idx].emplace_back(b_pair);
    adjacency[b_idx].emplace_back(a_pair);

    return Edge(this, edges.size() - 1, a_idx, b_idx);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    nodes.clear();
    edges.clear();
    adjacency.clear();
  }

  /** 
   * @brief Remove an edge formed by the two given nodes
   * 
   * @param[in] n1 A valid node at one end of the edge to remove
   * @param[in] n2 The other valid node at one end of the edge to remove
   * @return 1 if the edge was in the graph and has been removed. Return 0 if the edge was
   *         not in the graph.
   * 
   * @pre @a n1 and @a n2 are valid nodes
   * @post old num_edges() = new num_edges() + 1.
   * @post has_edge(@a n1, @a n2) == false;
   * @post All other Edge objects but the deleted one remain valid.
   * @post All outstanding iterators are invalidated
   * Complexity: O(d), d is the maxium degree of nodes in the graph.
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if(!has_edge(n1, n2)) {
      return 0;
    }

    size_type n1_idx = n1.node_idx;
    size_type n2_idx = n2.node_idx;
    
    for(size_type i = 0; i < adjacency[n1.node_idx].size(); ++i) {
      if(std::get<0>(adjacency[n1_idx][i]) == n2_idx) {
	adjacency[n1_idx][i] = adjacency[n1_idx].back();
	adjacency[n1_idx].pop_back();
      }
    }

    for(size_type i = 0; i < adjacency[n2_idx].size(); ++i) {
      if(std::get<0>(adjacency[n2_idx][i]) == n1_idx) {
	adjacency[n2_idx][i] = adjacency[n2_idx].back();
	adjacency[n2_idx].pop_back();
      }
    }

    for(size_type i = 0; i < edges.size(); ++i) {
      if((std::get<0>(edges[i]) == n1_idx && std::get<1>(edges[i]) == n2_idx) ||
	 (std::get<0>(edges[i]) == n2_idx && std::get<1>(edges[i]) == n1_idx)) {
	edges[i] = edges.back();
	edges.pop_back();
      }
    }

    return 1;
  }

  /** 
   * @brief Remove a given edge from the graph. 
   *        This method calls remove_edge(n1, n2).
   *
   * @param[in] e The valid edge to remove
   * @return 1 if the edge was in the graph and has been removed. Return 0 if the edge was
   *         not in the graph.
   * 
   * @pre @a e is valid edge
   * @post new num_edges() == old num_edges() - 1
   * @post has_edge(@a n1, @a n2) == false
   * @post Outstanding EdgeIterators are invalidated.
   * @post Outstanding IncidentIterators associated with @a e.node1() and
   *       @a e.node2() are invalidated.
   * Complexity: O(d), d is the maximum degree of a node in the graph
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** 
   * @brief Remove the edge that the edge iterator points to.
   *
   * @param[in] e_it The valid EdgeIterator to remove
   * @return The EdgeIterator pointing to the edge after the removed edge
   *
   * @pre The @a e_it is valid
   * @post old num_edges() = new num_edges() + 1.
   * @post has_edge(old *e_it) == false;
   * @post All other Edge objects but the deleted one remain valid.
   * @post All outstanding iterator are invalidated
   * Complexity: O(d), d is the maximum degree of a node in the graph
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    node_iterator node_i = e_it.node_it;
    incident_iterator inci_i = e_it.inci_it;
    remove_edge(*e_it);
    if(inci_i == (*node_i).edge_end()) {
      ++e_it;
    }
    return e_it;
  }

  /** 
   * @brief Remove the given node and its associated edges from the graph
   * 
   * @param[in] node The valid node to remove
   * @return If graph contains node return 1 else return 0
   *
   * @pre @a node is valid node
   * @pre graph.node(i).index() == i for all 0 <= i <= num_nodes()
   * @post old num_nodes() = new num_nodes() + 1.
   * @post old num_edges() = new num_edges() + old @a node.degree()
   * @post has_node(@a node) == false;
   * @post All other node objects but the deleted one remain valid.
   * @post All other edge objects but those linked to the node @a node remian valid.
   * @post All outstanding iterator of the graph are invalidated
   * Complexity: O(num_nodes()) given a sparse graph
   */
  size_type remove_node(const Node& node) {
    if(!has_node(node)) {
      return 0;
    }

    while(node.edge_begin() != node.edge_end()) {
      remove_edge(*node.edge_begin());
    }

    size_type index = node.index();
    nodes.erase(nodes.begin() + index);
    adjacency.erase(adjacency.begin()+index);
    for(auto it = adjacency.begin(); it != adjacency.end(); ++it) {
      for(auto innerIt = it->begin(); innerIt != it->end(); ++innerIt) {
	if(innerIt->first > index)
	  --innerIt->first;
      }
    }
    for(auto it = edges.begin(); it != edges.end(); ++it) {
      if(it->first > index)
	--it->first;
      if(it->second > index)
	--it->second;
    }
    return 1;
  }

  /** 
   * @brief Remove the given node and the edges from the graph
   *
   * @param[in] n_it is a valid interator
   * @return A NodeIterator pointing to node after the erased node
   *
   * @pre @a n_it is valid node iterator and @a n_it != node_end().
   * @pre graph.node(i).index() == i for all 0 <= i <= num_nodes()
   * @post old num_nodes() = new num_nodes() + 1.
   * @post old num_edges() = new num_edges() + old @a node.degree()
   * @post has_node(*n_it) == false;
   * @post All other node objects but the deleted one remain valid.
   * @post All other edge objects but those linked to the node @a node remian valid.
   * @post All outstanding iterator of the graph are invalidated
   * Complexity: O(num_nodes()) given a sparse graph
   */
  node_iterator remove_node(node_iterator n_it) {
    size_type index = (*n_it).index();
    if(remove_node(*n_it)) {
      node_iterator it = node_begin();
      for(; it != node_end(); ++it)
	if(index == (*it).index())
	  return it;
    }
    return n_it;
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
    //using iterator_category = std::forward_iterator_tag;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      return Node(parent_g, index);
      //return parent_g->node(index);
    }
    
    /**
     * @brief NodeIterator is incremented to point to next node
     *
     * @return NodeIterator that points to the next node object
     */
    NodeIterator& operator++() {
      ++index;
      return *this;
    }

    /**
     * @brief Determines if two NodeIterators are the equivalent
     *
     * @param[in] _n_ A NodeIterator object
     * @return bool Evaluates to true if:
     *         1. _n_ and this belong to the same graph
     *         2. n.index == index in the the global order
     */
    bool operator==(const NodeIterator& n) const {
      return (n.parent_g == parent_g) && (n.index == index);
    }

  private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph* parent_g;
    size_type index;

    /**
     * @brief Private constructor for a NodeIterator
     *
     * @param[in] _graph_ptr_ A pointer to an object of type Graph
     * @param[in] _index_arg_ Used to represent the _i_th node
     * @return NodeIterator that belongs to _graph_pter_ and points to the
     *         _i_th Node
     *
     * @pre 0 <= _i_ < num_nodes()
     */
    NodeIterator(const Graph* graph_ptr, size_type index_arg) :
      parent_g(const_cast<Graph*>(graph_ptr)),
      index(index_arg){
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /** Returns the first NodeIterator
   * @return NodeIterator that points to the first Node in the global order if graph is nonempty,
   *                       otherwise returns a nullptr
   */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Returns a NodeIterator that indicates the end of the nodes
   * @return NodeIterator that is the nullptr
   */
  NodeIterator node_end() const {
    return NodeIterator(this, num_nodes());
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
    // Supply definitions AND SPECIFICATIONS for:

    /**
     * @return The Edge that the this points to
     */
    Edge operator*() const{
      return Edge(parent_g, std::get<1>(parent_g->adjacency[node_index][edge_index]), node_index, std::get<0>(parent_g->adjacency[node_index][edge_index]));
    }

    /**
     * @return IncidentIterator such that it points to the next adjacent node
     */
    IncidentIterator& operator++(){
      ++edge_index;
      return *this;
    }

    /**
     * @brief Determines if two IncidentIterators are equivalent
     *
     * @param[in] _n_ An IncidentIterator
     * @return bool value
     *
     * @post bool evaluated to true if:
     *         1. n.parent_g == parent_g
     *         2. n.node_index == node_index
     *         3. n.edge_index == edge_index
     *       are satisfied
     */
    bool operator==(const IncidentIterator& n) const{
      return (n.parent_g == parent_g) && (n.node_index == node_index) && (n.edge_index == edge_index);
    }

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* parent_g;
    size_type node_index;
    size_type edge_index;

    /**
     * @brief Private constructor for Incident iterator
     *
     * @param[in] _graph_ptr_ A Graph object
     * @param[in] _node_index_ A node index value
     * @param[in] _edge_index_  An edge index value
     * @return IncidentIterator
     *
     * @pre 0 <= _node_index_ < num_nodes()
     * @pre 0 <= _edge_index_ < degree()
     */
    IncidentIterator(Graph* graph_ptr, size_type node_index_, size_type edge_index_) :
      parent_g(const_cast<Graph*>(graph_ptr)),
      node_index(node_index_),
      edge_index(edge_index_) {
    }
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
    // Supply definitions AND SPECIFICATIONS for:

    /**
     * @return Edge that this points to
     */
    Edge operator*() const {
      return (*inci_it);
      //return Edge(parent_g, edge_index, parent_g->edges[edge_index].first, parent_g->edges[edge_index].second);
    }

    /**
     * @return NodeIterator that points to the next NodeIterator
     */
    EdgeIterator& operator++() {
      while(node_it != parent_g->node_end()) {
	auto n = *node_it;
	while(inci_it != n.edge_end()) {
	  ++inci_it;
	  if(inci_it == n.edge_end()) {
	    break;
	  }
	  else if(n.node_idx <= (*inci_it).node2().node_idx) {
	    return *this;
	  }
	}
	++node_it;
	if(node_it != parent_g->node_end()) {
	  inci_it = (*node_it).edge_begin();
	}
      }
      return *this;
    };
      /**
      ++edge_index;
      return *this;
    };
      **/

    /**
     * @brief Determines if two EdgeIterators are equivalent
     *
     * @param[in] edge_it An EdgeIterator object
     * @return bool evaluated true only if
     *         1. EdgeIterator and this belong to the same graph
     *         2. the Edge EdgeIterator and this point to have the same edge number in a global order
     */
    bool operator==(const EdgeIterator& edge_it) const {
      return ((edge_it.parent_g == parent_g) && (edge_it.node_it == node_it));
      //return ((edge_it.parent_g == parent_g) && (edge_it.edge_index == edge_index));
    }

  private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* parent_g;
    node_iterator node_it;
    incident_iterator inci_it;
    //size_type edge_index;

    /**
     * @brief Private constructor for EdgeIterator
     *
     * @param[in] _graph_ptr_ A pointer to A Graph object
     * @param[in] _edge_index_ An edge index calue
     * @return An EdgeIterator that points to the Edge with _edge_index_ in global order
     *
     * @pre 0 <= _edge_index_ < num_edges()
     */
    EdgeIterator(const Graph* graph_ptr, node_iterator node_it_, incident_iterator inci_it_) ://, size_type edge_index_) :
      parent_g(const_cast<Graph*>(graph_ptr)),
      node_it(node_it_),
      inci_it(inci_it_) {
      //edge_index(edge_index_) {
    }
  };
  
  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /**
   * @return EdgeIterator which points to the first edge
   */
  EdgeIterator edge_begin() const {
    node_iterator node_it = node_begin();
    auto inci_it = (*node_it).edge_begin();
    return EdgeIterator(this, node_it, inci_it);
    //return EdgeIterator(this, 0);
  }

  /**
   * @return EdgeIterator that indicate sthe end of the edges
   */
  EdgeIterator edge_end() const {
    node_iterator node_it = node_end();
    incident_iterator inci_it = (*node_begin()).edge_begin();
    return EdgeIterator(this, node_it, inci_it);
    //return EdgeIterator(this, this->num_edges());
  }

};

#endif // CME212_GRAPH_HPP

