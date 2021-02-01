#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

/////////////////////////////////////////////////////////////
#include <unordered_map>
/////////////////////////////////////////////////////////////

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // HW0: YOUR CODE HERE

  ///////////////////////////////////////////////////////////
  // Predeclare the node struct
  struct internal_node;
  // Predeclare the edge struct
  struct internal_edge;
  // Predeclare the struct needed to have a inverse_edge map
  struct edge_key;
  struct EdgeKeyHasher;
  ///////////////////////////////////////////////////////////

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE

    //////////////////////////////////////////////
    size_ = 0;
    num_edges_ = 0;

    std::unordered_map<size_type, internal_node> nodes_; 
    // default constructor, empty map
    assert(nodes_.empty());

    std::unordered_map<size_type, internal_edge> edges_;
    // default constructor, empty map
    assert(edges_.empty());

    std::unordered_map<edge_key, size_type, EdgeKeyHasher> inverse_edges_;
    // default constructor, empty map
    assert(inverse_edges_.empty());
    //////////////////////////////////////////////
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
    Node() {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      // default constructor for all attributes
      //////////////////////////////////////////////      
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      return fetch().point;
      //////////////////////////////////////////////
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      return fetch().index_node;
      //////////////////////////////////////////////
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
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      if ((graph_ == n.graph_) && (index_node_ == n.index_node_)){  
        //check if both nodes point to the same graph and have same index
        return true;
      }
      //////////////////////////////////////////////
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
      //////////////////////////////////////////////
      assert(graph_ == n.graph_);
      //check if both nodes point to the same graph

      if (index_node_ < n.index_node_){
        return true;
      }
      //////////////////////////////////////////////
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE

    ///////////////////////////////////////////////////
    Graph* graph_; // Associated Graph instance
    size_type index_node_; // Node's unique identification number
    
    /** Private Constructor */
    Node(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_node_(index) {
    }

    /** Helper method to return the appropriate element.
     */
    internal_node& fetch() const {
      auto search = graph_->nodes_.find(index_node_);
      if (search != graph_->nodes_.end()) {
          return search->second;
      }
      assert(false);
    }
    ///////////////////////////////////////////////////

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
    //////////////////////////////////////////////////
    return size_;
    //////////////////////////////////////////////////
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
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////////

    internal_node new_internal_node;
    new_internal_node.point = position;
    new_internal_node.index_node = size_;
    
    nodes_.insert({size_,new_internal_node});

    ;
    // Returns a Node that points to the new element
    return Node(this, size_++); //post increment of size
    /////////////////////////////////////////////////////
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////////
    if (this != n.graph_){
      return false;
    }
    auto search = nodes_.find(n.index_node_); //Done in amortized time O(1)
    if (search != nodes_.end()) {
        return true;
    }
    return false;
    //////////////////////////////////////////////////////
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////////
    assert((0<=(int)i) && (i<=size_));
    auto search = nodes_.find(i); //Done in amortized time O(1)
    if (search != nodes_.end()) {
      return Node(this, i);
    }
    assert(false);
    //////////////////////////////////////////////////////
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
    Edge() {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      // default constructor for all attributes
      //////////////////////////////////////////////   
    }


    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      return fetch().node_1;
      //////////////////////////////////////////////
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      return fetch().node_2;
      //////////////////////////////////////////////
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      if ((graph_ == e.graph_)&& (index_edge_ == e.index_edge_)){
        return true;
      }
      //////////////////////////////////////////////
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      assert(graph_ == e.graph_);
      //check if both edges point to the same graph

      if (index_edge_ < e.index_edge_){
        return true;
      }
      //////////////////////////////////////////////
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE

    ///////////////////////////////////////////////////
    Graph* graph_; // Associated Graph instance
    size_type index_edge_; // Edge's unique identification number

    /** Private Constructor */
    Edge(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_edge_(index) {
    }

    /** Helper method to return the appropriate element.
     */
    internal_edge& fetch() const {
      auto search = graph_->edges_.find(index_edge_);
      if (search != graph_->edges_.end()) {
          return search->second;
      }
      assert(false);
    }

    ///////////////////////////////////////////////////

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
    //////////////////////////////////////////////////
    return num_edges_;
    //////////////////////////////////////////////////
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////////
    assert((0<=(int)i) && (i<=num_edges_));
    auto search = edges_.find(i); //Done in amortized time O(1)
    if (search != edges_.end()) {
      return Edge(this, i);
    }
    assert(false);
    //////////////////////////////////////////////////////
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////////
    if (!(has_node(a)) || !(has_node(b))){ // O(1) complexity
      return false;
    }

    edge_key possible_key;
    
    possible_key.index_node_1=a.index_node_;
    possible_key.index_node_2=b.index_node_;
    auto search = inverse_edges_.find(possible_key); //Done in amortized O(1)
    if (search != inverse_edges_.end()) {
        return true;
    }

    possible_key.index_node_1=b.index_node_;
    possible_key.index_node_2=a.index_node_;
    search = inverse_edges_.find(possible_key); //Done in amortized O(1)
    if (search != inverse_edges_.end()) {
        return true;
    }
    return false;
    //////////////////////////////////////////////////////
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
    ////////////////////////////////////////////////////
    assert((has_node(a)) && (has_node(b))); // O(1) complexity
    // check if a and b are valid nodes for this graph
    
    assert((a.index_node_!=b.index_node_)); 
    // check if a and b are different
    

    edge_key possible_key;
    
    possible_key.index_node_1=a.index_node_;
    possible_key.index_node_2=b.index_node_;
    auto search = inverse_edges_.find(possible_key); //Done in amortized O(1)
    if (search != inverse_edges_.end()) {
        return Edge(this, search->second);
    }

    possible_key.index_node_1=b.index_node_;
    possible_key.index_node_2=a.index_node_;
    search = inverse_edges_.find(possible_key); //Done in amortized O(1)
    if (search != inverse_edges_.end()) {
        return Edge(this, search->second);
    }

    internal_edge new_internal_edge;
    new_internal_edge.node_1 = a;
    new_internal_edge.node_2 = b;
    new_internal_edge.index_egde  = num_edges_;
    edges_.insert({num_edges_,new_internal_edge}); //Done in amortized O(1)
    inverse_edges_.insert({possible_key,num_edges_}); //Done in amortized O(1)
    
    // Returns an Edge that points to the new element
    return Edge(this, num_edges_++); //post-incrementation of num_edges_
    ////////////////////////////////////////////////////
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    ///////////////////////////////////////////////////
    size_=0;
    num_edges_=0;
    nodes_.clear();
    edges_.clear();
    inverse_edges_.clear();
    ///////////////////////////////////////////////////
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

  // HW0: YOUR CODE HERE

  ///////////////////////////////////////////////////////////
  // Internal type for set elements
  struct internal_node {
    Point point;   // The Point held by each node
    size_type index_node;  // The unique indentifier for a node
  };

  struct internal_edge {
    node_type node_1;   // The first node pointed by the edge
    node_type node_2;   // The second node pointed by the edge
    size_type index_egde;  // The unique indentifier for an edge
  };

  struct edge_key { //a hashable struct to be the key of inverse_edges_
    size_type index_node_1;   // The first node pointed by the edge
    size_type index_node_2;   // The second node pointed by the edge
    bool operator==(const edge_key &other) const{
      return (index_node_1 == other.index_node_1 
           && index_node_2 == other.index_node_2);
    }
  };
  struct EdgeKeyHasher{ //Hash function for edge_key. inspired from
  //https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
    std::size_t operator()(const edge_key& k) const{
      size_type hash1 = std::hash<size_type>{}(k.index_node_1); 
      size_type hash2 = std::hash<size_type>{}(k.index_node_2); 
      return hash1 ^ hash2;
    } 
  };

  size_type size_; // number of nodes 
  size_type num_edges_; // number of edges

  std::unordered_map<size_type, internal_node> nodes_;
  
  std::unordered_map<size_type, internal_edge> edges_;
  std::unordered_map<edge_key, size_type, EdgeKeyHasher> inverse_edges_;
  
  

  // Disable copy and assignment of a Graph
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;
  ///////////////////////////////////////////////////////////

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
