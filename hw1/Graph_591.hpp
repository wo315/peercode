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

  struct internal_node;
  struct internal_edge;

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

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
  Graph() 
	: nodes_(), num_nodes_(0), edges_(), num_edges_(0){
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph(){
      nodes_.clear();
      edges_.clear();
  }

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered<Node>{
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
      return graph_->nodes_[uid_].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /**
     * @brief Get the reference to the value of the node.
     *
     * @return The reference to the value of the node.
     *
     * This function can also be used as both setter
     * and getter of the value since reference is
     * returned.
     */
    node_value_type& value(){
        return graph_->nodes_[uid_].val;
    }

    /**
     * @brief Get the reference to the value of the node.
     *
     * @return The reference to the value of the node.
     *
     * Since the value here is const, this function can only
     * be used as the getter of the value, not the setter of
     * the value.
     */
    const node_value_type& value() const{
        return graph_->nodes_[uid_].val;
    }

    /**
     * @brief Get the degree of the node.
     *
     * @return The degree of the node.
     *
     * This function uses the adjacency map, which only contains
     * nodes that are not isolated as keys. Hence will directly
     * return 0 if a node's index is not key in the map.
     */
    size_type degree() const{
        auto curr = graph_->adjacency.find(uid_);
        if(curr == graph_->adjacency.end()){
            return 0;
        }else{
            return curr -> second.size();
        }
    }

    /**
     * @brief Returns an Incident Iterator at first edge in the edge list of the node.
     *
     * @return incident_iterator pointing at the first edge.
     */
    incident_iterator edge_begin(){
        return IncidentIterator(graph_, this, 0);
    }

    /**
     * @brief Returns an Incident Iterator at last edge in the edge list of the node.
     *
     * @return incident_iterator pointing at the last edge.
     */
    incident_iterator edge_end(){
        return IncidentIterator(graph_, this, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->graph_ == n.graph_ && uid_ == n.index()){
            return true;
        }else{
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
      // HW0: YOUR CODE HERE
      if(uid_ < n.index()){
            return true;
        }else{
            return false;
        }
    }

   private:
     // Pointer back to the Graph container
     Graph* graph_;
     // This node's unique identification number
     size_type uid_;
     /** Private Constructor */
     Node(const Graph* graph, size_type uid)
         : graph_(const_cast<Graph*>(graph)), uid_(uid) {
     }
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return num_nodes_;
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
  Node add_node(const Point& position, const node_value_type& nodeval = node_value_type()) {
    // HW0: YOUR CODE HERE
    // HW1: modified hw0 code
    ++num_nodes_;
    Node node = Node(this, num_nodes_-1);
    internal_node i_node;
    i_node.p = position;
    i_node.uid = num_nodes_-1;
    i_node.val = nodeval;
    nodes_.push_back(i_node);
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.index() < size()){
          return true;
      }else{
          return false;
      }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_nodes_);
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
  class Edge: private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge(){}

    /** Return a node of this Edge */
    Node node1() const{
      // HW0: YOUR CODE HERE
      return graph_->edges_[uid_].node1;
    }

    /** Return the other node of this Edge */
    Node node2() const{
      // HW0: YOUR CODE HERE
      return graph_->edges_[uid_].node2;
    }

    size_type index() const{
      // HW0: YOUR CODE HERE
      return uid_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const{
      if(uid_ == e.index()){
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
    bool operator<(const Edge& e) const{
      if(uid_ < e.index()){
	return true;
      }else{
	return false;
      }
    }

   private:
     // Pointer back to the Graph container
     Graph* graph_;
     // This edge's unique identification number
     size_type uid_;
     /** Private Constructor */
     Edge(const Graph* graph, size_type uid)
         : graph_(const_cast<Graph*>(graph)), uid_(uid){
     }
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
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
    assert(i < num_edges_);
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const{
    // HW0: YOUR CODE HERE
    auto node_a = adjacency.find(a.index());
    if(node_a == adjacency.end()){
        return false;
    }
    std::map<size_type, size_type> a_adj = node_a -> second; // map of adjacent node id : edge id
    auto node_b = a_adj.find(b.index());
    if (node_b != a_adj.end()){
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
  Edge add_edge(const Node& a, const Node& b){
    // HW0: YOUR CODE HERE
    if(has_edge(a,b)){
        size_type i = adjacency[a.index()][b.index()];
        edges_[i].node1 = a;
        edges_[i].node2 = b;
        return Edge(this, i);
    }
    ++num_edges_;
    Edge edge = Edge(this, num_edges_-1);
    internal_edge i_edge;
    i_edge.node1 = a;
    i_edge.node2 = b;
    i_edge.uid = num_edges_-1;
    edges_.push_back(i_edge);
    adjacency[a.index()][b.index()] = i_edge.uid;
    adjacency[b.index()][a.index()] = i_edge.uid;
    return edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear(){
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
    num_nodes_ = 0;
    num_edges_ = 0;
    adjacency.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable<NodeIterator>{
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
     * @brief Dereference the current node.
     *
     * @return A Node object at the position of the NodeIterator.
     */
    Node operator*() const{
        return Node(graph_, node_id_);
    }

    /**
     * @brief The NodeIterator goes to the next node.
     *
     * @return Reference to the NodeIterator at the next position of the current node.
     */
    NodeIterator& operator++(){
        node_id_ ++;
        return *this;
    }

    /**
     * @brief Compares if two NodeIterators are the same.
     *
     * @param n_i[in] Reference to the other NodeIterator.
     * @return A boolean indicating if this and the other NodeIterator are the same one.
     */
    bool operator==(const NodeIterator& n_i) const{
        return (n_i.graph_ == graph_) && (n_i.node_id_ == node_id_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const graph_type* graph_;
    size_type node_id_;

    /**
     * @brief Constructor of the NodeIterator.
     *
     * @param graph[in] A pointer to the const graph_type.
     * @param node_id[in] Id of the node to be pointed at.
     *
     * Initialize graph_ to graph, node_id_ to node_id.
     */
    NodeIterator(const graph_type* graph, size_type node_id){
        graph_ = graph;
        node_id_ = node_id;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
   * @brief A node iterator to the first node in the graph.
   *
   * @return node_iterator pointing to the first node in the graph.
   */
  node_iterator node_begin() const{
      return NodeIterator(this, 0);
  }

  /**
   * @brief A node iterator to the last node in the graph.
   *
   * @return node_iterator pointing to the last node in the graph.
   */
  node_iterator node_end() const{
      return NodeIterator(this, this->num_nodes_);
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /**
     * @brief Dereference the current edge.
     *
     * @return The edge object that's pointed by the iterator.
     */
    Edge operator*() const{
        size_type node_id_ = node_->index();
        std::map<size_type, size_type> node_adj = graph_->adjacency.find(node_id_)->second;
        auto it = node_adj.begin();
        std::advance(it, incident_id_);
        size_type edge_id_ = it -> second;
        return Edge(graph_, edge_id_);
    }

    /**
     * @brief The iterator goes to the next edge.
     *
     * @return Reference to the IncidentIterator pointing at the next edge.
     */
    IncidentIterator& operator++(){
        incident_id_++;
        return *this;
    }

    /**
     * @brief Check if the two iterators are the same.
     *
     * @param i_i[in] Reference to the other const IncidentIterator.
     * @return An boolean indicating if this and the other IncidentIterator are pointing at the same edge.
     */
    bool operator==(const IncidentIterator& i_i) const{
        return (i_i.graph_ == graph_) && (i_i.node_ == node_) && (i_i.incident_id_ == incident_id_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const graph_type* graph_;
    node_type* node_;
    size_type incident_id_;

    /**
     * @brief Constructor of the IncidentIterator.
     * @param graph[in] A pointer to the graph.
     * @param node[in] A pointer to the node.
     * @param incident_id[in] The index incident edge.
     */
    IncidentIterator(const graph_type* graph, node_type* node, size_type incident_id) {
        graph_ = graph;
        node_ = node;
        incident_id_ = incident_id;
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /**
     * @brief Dereference the edge.
     *
     * @return The Edge pointed by the iterator.
     */
    Edge operator*() const{
        return Edge(graph_, edge_id_);
    }

    /**
     * @brief Iterator goes to the next edge.
     *
     * @return Reference to the iterator pointing at the next edge.
     */
    EdgeIterator& operator++(){
        edge_id_ ++;
        return *this;
    }

    /**
     * @brief Compare if two EdgeIterators are the same.
     *
     * @param e_i[in] Reference to the other const EdgeIterator
     * @return A boolean indicating if this and the other EdgeIterator are the same.
     */
    bool operator==(const EdgeIterator& e_i) const{
        return (e_i.graph_ == graph_) && (e_i.edge_id_ == edge_id_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const graph_type* graph_;
    size_type edge_id_;

    /**
     * @brief Constructor of the EdgeIterator.
     *
     * @param graph[in] Pointer to the const graph.
     * @param edge_id[in] Edge id.
     */
    EdgeIterator(const graph_type* graph, size_type edge_id){
          graph_ = graph;
          edge_id_ = edge_id;
      }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
   * @brief Get the edge_iterator pointing to the first edge.
   *
   * @return An edge iterator pointing to the first edge.
   */
  edge_iterator edge_begin() const{
      return EdgeIterator(this, 0);
  }

  /**
   * @brief Get the edge_iterator pointing to the last edge.
   *
   * @return An edge iterator pointing to the last edge.
   */
  edge_iterator edge_end() const{
      return EdgeIterator(this, this->num_edges_);
  }


 private:
   struct internal_node {
        Point p;
        size_type uid;
        node_value_type val;
      };
    std::vector<internal_node> nodes_;
    size_type num_nodes_;

   struct internal_edge {
        Node node1;
        Node node2;
	    size_type uid;
      };
    std::vector<internal_edge> edges_;
    size_type num_edges_;
    std::map<size_type, std::map<size_type, size_type>> adjacency; // key: node id, value: map <node id, edge id>
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
