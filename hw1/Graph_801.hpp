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

  /** type values contained by nodes (e.g. mass, temperature or color) */
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
      : num_nodes_(0), num_edges_(0) {}

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      return graph_ptr_->nodes_[index()].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    // HW1: YOUR CODE HERE

    /**  Return a reference to this node's value of type node_value_type */
    node_value_type& value() {
      return graph_ptr_->nodes_[index()].value();
    }

    /** Return a reference to this node's value of type const node_value_type */
    const node_value_type& value() const {
      return graph_ptr_->nodes_[index()].value();
    }
  
    /** Return this node's degree (number of adjecent nodes) of type size_type */
    size_type degree() const {
      return graph_ptr_->nodes_[index()].degree_;
    }
    
    /** begin iterator for this node corresponding to its "first" edge */
    incident_iterator edge_begin() const {
      return ( IncidentIterator(*this, 0) );
    }

    /** end iterator for this node corresponding to its "last" edge */
    incident_iterator edge_end() const {
      return ( IncidentIterator(*this, degree()) );
    }

    /** Test whether this node and @a n are equal.
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ptr_ == n.graph_ptr_) && (index() == n.index());
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
      return index() < n.index();
    }

   private:
    friend class Graph;
    graph_type* graph_ptr_; 
    size_type index_;
    Node(const graph_type* graph_ptr, size_type index)
        : graph_ptr_(const_cast<graph_type*>(graph_ptr)), index_(index) {};
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
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
  Node add_node(const Point& position, const node_value_type& value=node_value_type()) {
    nodes_.push_back(InternalNode(this, position, value));
    return Node(this, num_nodes_++);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.index() < size()) && (this == n.graph_ptr_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if ( i >= size() )
      throw std::runtime_error("Node index out of range.");
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
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return graph_ptr_->node(n1_index_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_ptr_->node(n2_index_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_ptr_ == e.graph_ptr_) && (index_) == e.index_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return index_ < e.index_;
    }    

   private:
    friend class Graph;
    graph_type* graph_ptr_;
    size_type index_;
    size_type n1_index_;
    size_type n2_index_;

    Edge(graph_type* graph_ptr, size_type index, size_type n1_index, size_type n2_index)
      : graph_ptr_(graph_ptr), index_(index), n1_index_(n1_index), n2_index_(n2_index) {};
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
    if ( i >= num_edges() )
      throw std::runtime_error("Edge index out of range.");
    return edges_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return ( nodes_[a.index()].is_adjacent(b) );
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

    // return edge of this index if exists
    int edge_idx = nodes_[a.index()].adjacency_index(b); 
    if (edge_idx != -1) 
        return Edge(this, (size_type)edge_idx, a.index(), b.index());

    // otherwise add it
    nodes_[a.index()].add_adjacent(b.index(), num_edges());
    nodes_[b.index()].add_adjacent(a.index(), num_edges());
    Edge e = Edge(this, num_edges_++, a.index(), b.index());
    edges_.push_back(e);
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    num_nodes_ = 0;
    num_edges_ = 0;
  }

  //
  // Node Iterator
  //
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : 
    private totally_ordered<NodeIterator> {

   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    // HW1 #2: YOUR CODE HERE
    /** Dereferencing the iterator returns the corresponding node in O(1).
     *  @post (*result_node_iterator).index() == graph.node(node_iterator.index_)
     */
    Node operator*() const {
      return Node(graph_ptr_, index_);
    }

    /** Increment the iterator by incrementing the corresponding nodes index in O(1).
     *  @post (*result_node_iterator).index() == (*old_node_iterator).index()+1 
    */
    NodeIterator& operator++() {
      ++index_;
      return *this;
    }

    /** Equal comparison for the iterator in O(1).
     * We consider two iterators equal if they're from the same graph and same node index.
     */
    bool operator==(const NodeIterator& node_it) const {
      return (graph_ptr_==node_it.graph_ptr_) && (index_ == node_it.index_);
    }

   private:
    // HW1 #2: YOUR CODE HERE
    friend class Graph;
    friend class EdgeIterator;
    graph_type* graph_ptr_;
    unsigned int index_;
    NodeIterator(const graph_type* graph_ptr, size_type index)
      : graph_ptr_(const_cast<graph_type*>(graph_ptr)), index_(index) {};
  };

  // HW1 #2: YOUR CODE HERE

  /** Starting node iterator for the graph.
   * @post *result_node_iterator == graph.node(0)
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Ending node iterator for the graph.
   * Returns the the iterator corresponding to one passed the nodes
   * @post *result_node_iterator == graph.node(num_nodes) (undefined)
   */
  node_iterator node_end() const {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {

   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    // HW1 #3: YOUR CODE HERE

    /** Dereferencing the iterator returns the corresponding edge.
     * That is, the edge between the given node that we're iterating over and 
     * current adjacenct node. O(1) time.
     * @post restult_edge.node1() == this node
     * @post restult_edge.node2() == current adjacent node
     */
    Edge operator*() const {
      const InternalNode* in_ptr = &(graph_ptr_->nodes_[node_.index()]);
      size_type edge_idx = in_ptr->edge_idxs_[inc_idx_];
      size_type adj_node_idx = in_ptr->adj_nodes_[inc_idx_];
      return Edge(graph_ptr_, edge_idx, node_.index(), adj_node_idx);
    }

    /** Increment the iterator by incrementing the corresponding adjacent nodes index 
    */
    IncidentIterator& operator++() {
      inc_idx_++;
      return *this;
    }

    /** Equal comparison for the iterator. 
     * We consider two iterators equal if they're from the same graph,
     * same node index and same adjacent node index.
     */
    bool operator==(const IncidentIterator& ii) const {
      return (graph_ptr_ == ii.graph_ptr_) && (node_ == ii.node_) && (inc_idx_ == ii.inc_idx_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Node node_;
    size_type inc_idx_;
    graph_type* graph_ptr_;
    
    IncidentIterator(Node node, size_type inc_idx)
      : node_(node), inc_idx_(inc_idx) {
        graph_ptr_ = node.graph_ptr_;
      };
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
    EdgeIterator() {}

    // HW1 #5: YOUR CODE HERE

    /** Dereferencing the edge iterator comes down do dereferencing it's current 
     *  incident iterator which returns an edge. Edge-ordering of the nodes is 
     *  dependent on the ordering of the nodes:
     *  @post restult_edge.node1() < restult_edge.node2()
     *  so restult_edge.node1() is not necessarily node corresponding to node iterator
     */
    Edge operator*() const { 
      return *inc_it_;
    }

    /** Incrementing the edge iterator. 
     * While we still haven't reached the end of the graph,
     * while we still haven't reached the end of the nodes adjecent edges,
     * if the dereferenced edge has node1() < node2() :return that edge
     * else: keep going 
     * 
     * @post there exists no edge e' s.t. ++e' == result_edge ...
     *  ... except for this particular pre_incremented_edge 
    */
    EdgeIterator& operator++() {
      ++inc_it_;
      while ( n_it_ != (graph_ptr_->node_end()) ) {
        while ( inc_it_ != (*n_it_).edge_end() ) {
          if ( (*(inc_it_)).node1() < (*(inc_it_)).node2() ) // avoid double visits
            return *this;
          ++inc_it_;
        }
        inc_it_ = IncidentIterator(*(++n_it_), 0);
      }
      return *this;
    }

    /** Equal comparison for the iterator. 
     * We consider two iterators equal if they're from the same graph and same adjacent node index.
     * Node that node index and adjecent node index are tied together, meaning 
     * inc_it == inc_it' => n_it_ == n_it_'
     */
    bool operator==(const EdgeIterator& edge_it) const {
      return (graph_ptr_ == edge_it.graph_ptr_) && (inc_it_ == edge_it.inc_it_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_ptr_;
    NodeIterator n_it_;
    IncidentIterator inc_it_;
    EdgeIterator(const graph_type* graph_ptr, NodeIterator node_it, IncidentIterator inc_it)
      : graph_ptr_(const_cast<graph_type*>(graph_ptr)), n_it_(node_it), inc_it_(inc_it) {};
  };

  // HW1 #5: YOUR CODE HERE
  edge_iterator edge_begin() const {
    NodeIterator start_node_it = node_begin();
    IncidentIterator start_inc_it = (*start_node_it).edge_begin();
    return EdgeIterator(this, start_node_it, start_inc_it);
  }

  edge_iterator edge_end() const {
    NodeIterator end_node_it = node_end();
    IncidentIterator end_inc_it = (*end_node_it).edge_begin();
    return EdgeIterator(this, end_node_it, end_inc_it);
  }

 private:

  class InternalNode {
  /** @class Graph::InternalNode
  * @brief Helper class for creating an ordering for edge's nodes.
  */
   public:
    InternalNode() {}

    /** Return this node's value of type node_value_type */
    node_value_type& value() {
       return value_;
    }

    /** Finding the adjacency index to n, returns -1 if not found.
     * @pre @a n is from same graph as this node
     * @post -1 <= return_index < degree <= graph.num_nodes
     * Complexity is O(this node's degree)
     */
    int adjacency_index(const Node& n) const {
      for (size_type a=0; a<degree_; a++) {
        if (adj_nodes_[a] == n.index())
          return a;
      }
      return -1;
    }

    /** Determine if n is adjacent to this node
     *  @pre @a n is from same graph as this node
     *  Uses adjacency_index() so complexity is same: O(this node's degree)
     */
    bool is_adjacent(const Node& n) const { 
      return adjacency_index(n) != -1;
    }

   private:
    friend class Graph;

    /** Method for adding node indices to the adjecency list 
     * @pre @a adj_idx < num_nodes 
     * @post adj_nodes_.size() == edge_idxs_.size() == degree_ == old_degree + 1
     * @post result_node.index() = adj_idx
    */
    Node add_adjacent(size_type adj_idx, size_type edge_idx) {
      adj_nodes_.push_back(adj_idx);
      edge_idxs_.push_back(edge_idx);
      ++degree_;
      return graph_ptr_->node(adj_idx);
    }

    // internal node's members
    std::vector<size_type> adj_nodes_ {};
    std::vector<size_type> edge_idxs_ {};
    graph_type* graph_ptr_;
    const Point position_;
    node_value_type value_;
    size_type degree_;

    InternalNode(graph_type* graph_ptr, const Point& position, const node_value_type& value)
    : graph_ptr_(graph_ptr), position_(position), value_(value), degree_(0) {};
  };

  // graph members
  size_type num_nodes_;
  size_type num_edges_;
  std::vector<InternalNode> nodes_ {};
  std::vector<Edge> edges_ {};
};

#endif // CME212_GRAPH_HPP
