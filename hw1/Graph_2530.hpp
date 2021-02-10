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
  using node_value_type = V;

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
    // HW0: YOUR CODE HERE
    // zero or empty initialize all member variables
    : size_(0), num_edges_(0), points_(), nodes_(), edges_(), node_values_(), adjacency() {
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->points_[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** @brief Return this node's value
     */
    node_value_type& value() {
      return graph_->node_values_[index_];
    }

    /** @brief Return this const node's value
     */
    const node_value_type& value () const {
      return graph_->node_values_[index_];
    }

    /** @brief Return this node's degree -- the number
     * of nodes adjacent to it in the graph.
     */
    size_type degree() const {
      auto iter = graph_->adjacency.find(index_);
       return (iter->second).size();
    }

    /** @brief Return an incident iterator for this node that
     * points to the beginning of a collection contaning
     * this node's adjacent nodes.
     */
    incident_iterator edge_begin() const {
      auto iter = graph_->adjacency.find(index_);  
       return IncidentIterator(graph_, (iter->second).begin(), index_);
     }

    /** @brief Return an incident iterator for this node that
     * points to after the end of a collection contaning
     * this node's adjacent nodes.
     */
    incident_iterator edge_end() const {
       auto iter = graph_->adjacency.find(index_); 
       return IncidentIterator(graph_, (iter->second).end(), 0);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((n.graph_ == this->graph_) && (n.index_ == this->index_)) {
          return true;
      }
      return false;
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
      // return if this node's index is less than n
      if (this->index_ < n.index_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // pointer to node's graph
    Graph* graph_;

    // index of node in graph
    size_type index_;


    /** @brief private constructor for a node.
     *
     * @param[in] graph a pointer to the graph.
     * @param[in] index the index of the node.
     *
     * @pre index >= 0.
     * @post Construct a valid node @a n of the graph with n.index_ = index.
     */
    Node(const graph_type* graph, size_type index)
      : graph_(const_cast<graph_type*>(graph)), index_(index) {
      }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    Node n(this, size_);
    nodes_.push_back(n);
    points_.push_back(position);
    node_values_.push_back(value);
    std::vector<size_type> empty;
    adjacency.insert({size_, empty});
    size_ += 1;

    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (nodes_[n.index_] == n) {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < size_);
    return nodes_[i];
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[n1_];
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[n2_];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      /** both edges have the same nodes regardless of order
       * so first check if node 1 and node 2 of each edge match;
       * otherwise check if node 1 of edge e and node 2 of 
       * this edge match and if node 2 of edge e and node 1 of 
       * this edge match.
       */ 
      
      if ((e.n1_ == this->n1_) && (e.n2_ == this->n2_)) {
        return true;
      } else if ((e.n1_ == this->n2_) && (e.n2_ == this->n1_)) {
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

      /** if edges are not equal then 
       * check if the sum of this edge's node indices is
       * less than the sum of e's node indices. If yes,
       * return true; otherwise return false.
       */ 
      size_type sum_e = e.n1_ + e.n2_;
      size_type sum_this = this->n1_ + this->n2_;

      if (*this == e) {
        return false;
      }
      else if (sum_this < sum_e) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // index for first node
    size_type n1_;

    // index for second node
    size_type n2_;

    // pointer to graph
    Graph* graph_;

    /** @brief private constructor for an edge.
     *
     * @param[in] graph a pointer to the graph.
     * @param[in] n1 the index of the first node.
     * @param[in] n2 the index of the second node.
     *
     * @pre index >= 0.
     * @post Construct a valid edge @a e of the graph.
     */
    Edge(const Graph* graph, size_type n1, size_type n2)
      : graph_(const_cast<Graph*>(graph)), n1_(n1), n2_(n2) {
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
    assert(i < num_edges());
    return edges_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for (Edge e : edges_) {
      if ((e.node1() == a) && (e.node2() == b)) {
        return true;
      }
      else if ((e.node1() == b) && (e.node2() == a)) {
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
    // HW0: YOUR CODE HERE
    // check if pre-conditions are fulfilled.
    if ((a.graph_ == b.graph_) && (a.graph_ == this) && (a.index_ != b.index_)) {

      /** if edge exists then return the edge found by edge_between;
       * else, create new edge and return it.
       */
      if (has_edge(a, b)) {
        return Edge(this, a.index_, b.index_);
      } else {
        /** if no existing edge that satisfies post-conditions
        * is found then create a new edge.
        */ 
        Edge e(this, a.index_, b.index_);
        adjacency[a.index_].push_back(b.index_);
        adjacency[b.index_].push_back(a.index_);
        edges_.push_back(e);
        num_edges_ += 1;
        return e;
      }      
    }
    return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // clear all vectors and set other variables to 0.
    edges_.clear();
    nodes_.clear();
    points_.clear();
    node_values_.clear();
    adjacency.clear();
    num_edges_ = 0;
    size_ = 0;
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
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** @brief overrides the usual operator*().
     * 
     *  @return the node this iterator points to.
     */
    Node operator*() const {
      return *iter_;
    }
    
    /** @brief overrides the usual operator++().
     * 
     *  @return the node iterator pointing to the next node.
     */
    NodeIterator& operator++() {
      ++iter_;
      return *this;
    }
    
    /** @brief overrides the usual operator==().
     * 
     *  @param[in] i node_iterator to compare this to.
     * 
     *  @return boolean of whether these two iterators are equal.
     */
    bool operator==(const NodeIterator& i) const {
      return iter_ == i.iter_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    // pointer to graph
    Graph* graph_;

    // iterator to vector of nodes
    typename std::vector<Node>::const_iterator iter_;

    /** @brief private constructor for a node iterator.
     *
     * @param[in] graph a pointer to the graph.
     * @param[in] iter iterator of a vector of nodes
     *
     * @post Construct a valid node_iterator 
     */
    NodeIterator(const Graph* graph, typename std::vector<Node>::const_iterator iter)
    : graph_(const_cast<Graph*>(graph)), iter_(iter) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** @brief access to pointer to this graph's first node
   * 
   *  @return iterator pointing to this graph's first node
   */
  node_iterator node_begin() const {
    return NodeIterator(this, nodes_.begin());
  }
  
  /** @brief access to pointer to after this graph's first node
   * 
   *  @return iterator pointing to after this graph's first node
   */
  node_iterator node_end() const {
    return NodeIterator(this, nodes_.end());
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
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** @brief overrides the usual operator*().
     * 
     *  @return the edge this between root node and node the 
     *          iterator points to with node1 = root node
     *          and node2 = incident node. 
     */
    Edge operator*() const {
      size_type incidence_idx = *iter_;
      return Edge(graph_, root_node, incidence_idx);
    }

    /** @brief overrides the usual operator++().
     * 
     *  @return the incident iterator pointing to the next node.
     */
    IncidentIterator& operator++() {
      ++iter_;
      return *this;
    }
    
    /** @brief overrides the usual operator==().
     * 
     *  @param[in] i incident_iterator to compare this to.
     * 
     *  @return boolean of whether these two iterators are equal.
     */
    bool operator==(const IncidentIterator& i) const {
      return i.iter_ == iter_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    // pointer to graph
    Graph* graph_;

    // iterator over a vector of node indices
    std::vector<size_type>::iterator iter_;

    // root node index for the incidence traversal
    size_type root_node;
    

    /** @brief private constructor for an incident iterator.
     *
     * @param[in] graph a pointer to the graph.
     * @param[in] iter iterator of a vector of nodes
     *
     * @post Construct a valid incident_iterator 
     */
    IncidentIterator(const Graph* graph, std::vector<size_type>::iterator iter, size_type index)
    : graph_(const_cast<Graph*>(graph)), iter_(iter), root_node(index) {}

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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

     /** @brief overrides the usual operator*().
     * 
     *  @return the edge this iterator points to.
     */
    Edge operator*() const {
      return *iter_;
    }

    /** @brief overrides the usual operator++().
     * 
     *  @return the edge iterator pointing to the next edge.
     */
    EdgeIterator& operator++() {
      ++iter_;
      return *this;
    }

    /** @brief overrides the usual operator==().
     * 
     *  @param[in] i edge_iterator to compare this to.
     * 
     *  @return boolean of whether these two iterators are equal.
     */
    bool operator==(const EdgeIterator& i) const {
      return (i.iter_ == iter_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    // pointer to graph
    Graph* graph_;

    // iterator over edge vector
    typename std::vector<Edge>::const_iterator iter_;

    /** @brief private constructor for an edge iterator.
     *
     * @param[in] graph a pointer to the graph.
     * @param[in] iter iterator of a vector of edges
     *
     * @post Construct a valid edge_iterator 
     */
    EdgeIterator(const Graph* graph, typename std::vector<Edge>::const_iterator iter)
    : graph_(const_cast<Graph*>(graph)), iter_(iter) {}


  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** @brief return pointer to first edge of this graph
   * 
   *  @return edge iterator pointing to first edge of this graph
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edges_.begin());
  }
  
  /** @brief return pointer after the first edge of this graph
   * 
   *  @return edge iterator pointing after first edge of this graph
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.end());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // vector for points
  std::vector<Point> points_;

  // corresponding vector for nodes
  std::vector<Node> nodes_;

  // corresponding vector for node values
  std::vector<node_value_type> node_values_;

  // vector of edges
  std::vector<Edge> edges_;

  // size of graph
  size_type size_;

  // number of edges
  size_type num_edges_;

  // adjacency list
  std::map<size_type, std::vector<size_type> > adjacency;




};

#endif // CME212_GRAPH_HPP
