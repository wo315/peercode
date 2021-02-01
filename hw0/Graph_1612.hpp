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
class Graph {
 
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    size_ = 0;
    num_edges_ = 0;

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
  class Node {
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
    Node() { ///Invalid constructor because a user should only construct nodes through the Graph class.
  
    }

    /** Return this node's position. */
    /// Because the position info of a node is stored in the Points and Nodes vector,
    /// we can return the position from either vector. But because the given function returns Point, we'll use the point vector.
    const Point& position() const {

      return graph_->points_[index_]; 
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      /// Because a node is uniquely defined by its index and this index is up to date. We return index_
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      /// Because the instruction says to check graph and index, we'll only check these. 
      /// One could also check for graph, index, and point
      if (n.graph_ == this->graph_) {
        if (n.index_ == this->index_) {
          return true;
        }
      }

      (void) n;          // Quiet compiler warning
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
      /// Currently, we can only check for the index of the nodes.
      /// Note this only evaluates the order the nodes were added which is kind of a mute point for < to do.
      
      if (this->index_ < n.index_) {
        return true;
      }

      (void) n;           // Quiet compiler warning
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    //A proper node constructor will have a pointer to the graph it belongs to, and an index. 
    //Note: it won't contain any points to keep it lightweight
    graph_type* graph_;
    size_type index_;

    Node(const graph_type* graph, size_type index)
      : graph_(const_cast<graph_type*> (graph)), index_(index) {

      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    /// Because the graph class has a size attribute, just return that here.
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
  Node add_node(const Point& position) {
    /// Because our node is tracked by the point and node vectors, 
    /// we will create a node to add to the node vector at the end of the existing Node array and add the given point to the point vector.
    /// Because we are using standard vectors, we can use the standard push_back function

    Node new_node(this, size_); //the new node should be at index size_ for 0 index
    nodes_.push_back(new_node);
    points_.push_back(position);
    size_ += 1; 

    (void) position;      // Quiet compiler warning
    return new_node;      
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    /// Note: all nodes have a pointer to their respective graph.
    /// So, if the Node's graph is the same as this graph we are in, then the Node exists.
    /// This is a much better approach [O(1)] than iterating through each node in a graph and checking for equality

    if (n.graph_ == this) {
      return true;
    }
    (void) n;            // Quiet compiler warning
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    /// Because, there is a Node vector of all nodes, just return the ith element of that vector
    assert(i < size());
    (void) i;             // Quiet compiler warning
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() { //This is the invalid constructor.
    
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->nodes_[idx1_];      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->nodes_[idx2_];      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning

      /// This is checking the reference equality instead of the value equality 
      /// which would have checked the equality of the nodes at the end of the edges
      return ((this->idx1_ == e.idx1_ || this->idx1_ == e.idx2_) && 
        (this->idx2_ == e.idx2_ || this->idx2_ == e.idx1_));   
    
   }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      /// This compares the length of the edges.
      (void) e;           // Quiet compiler warning
      size_type sum_given = e.idx1_ + e.idx2_;
      size_type sum_existing = this->idx1_ + this->idx2_;

      if (sum_existing < sum_given) {
        return true;
      }

      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Constructor
    // An edge connects two nodes. Hence, it is uniquely defined by the indices of those two nodes
    graph_type* graph_;
    size_type idx1_;
    size_type idx2_;


    Edge(const graph_type* graph, size_type idx1, size_type idx2)
      : graph_(const_cast<graph_type*> (graph)), idx1_(idx1), idx2_(idx2) 
      {

      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */


  size_type num_edges() const {
    /// Because we are already tracking the number of edges in Graph, return that here
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    /// because all edges are stored in the edges_ vector, return the ith index.
    assert(i < (num_edges_));
    (void) i;             // Quiet compiler warning
    return edges_[i];        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
   /// Because the edge is undirected, we need to check both points.

    for (size_type i = 0; i < num_edges_; ++i) {
      Edge e = edges_[i];
      Node e1 = e.node1();
      Node e2 = e.node2();

      return ((e1.index_ == a.index_) && (e2.index_ == b.index_)) || 
        ((e1.index_ == b.index_) && (e2.index_ == a.index_));

    }

    (void) a; (void) b;   // Quiet compiler warning
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
    size_type idx_a = a.index_;
    size_type idx_b = b.index_;

    Edge e(this, idx_a, idx_b);

    if(has_edge(a, b)){
      return edges_[idx_b];
    } else {
      edges_.push_back(e);
    }

    num_edges_ += 1;

    (void) a, (void) b;   // Quiet compiler warning
    return e;     
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */

  /// Because all the nodes and edges are stored in a vector, we can just clear the vector

  void clear() {
    nodes_.clear();
    edges_.clear();
    points_.clear();

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

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

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

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
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

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:

  // NODE: represented by points and indices
  // vector of points

  std::vector<Point> points_;

  // vecotr of Nodes
  std::vector<Node> nodes_;

  //EDGE
  //Vector of edges
  std::vector<Edge> edges_;

  //OTHER
  size_type size_; //size of the graph
  size_type num_edges_; //number of edges in a graph
};

#endif // CME212_GRAPH_HPP
