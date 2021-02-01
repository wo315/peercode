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
#include <iostream>
#include <assert.h>

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

  // HW0: WIP Check if forward declaration is needed 
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct real_node;
  struct real_edge;


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
  Graph()
    : nodes_storage(), nxt_node_uid_(0), node_cnt(0), node_idx_to_uid()
      , connected_nodes()
      , edges_storage(), nxt_edge_uid_(0), edge_cnt(0), edge_idx_to_uid() {
    // HW0: Initiate empty graph 
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
      // HW0: Public Constructor, create an invalid Node
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: travel to graph, nodes_storage, and use uid as key to retrive
      // the real node, then return it's point. 
      return graph_ -> nodes_storage[uid_].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: travel to graph, nodes_storage, and use uid as key to retrive
      // the real node, then return it's index
      return graph_ -> nodes_storage[uid_].idx;
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
      // HW0: 
      return (this->graph_ == n.graph_) && 
      (this->graph_->nodes_storage[uid_].idx 
        == n.graph_->nodes_storage[uid_].idx);
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
      // HW0: < is overloaded with index comparison. So not possible to have 
      // x == y unless it is the same x (index).
      // Switch to uid if we want this comparison to be invariant. 
      return 
      this->graph_->nodes_storage[uid_].idx < n.graph_->nodes_storage[uid_].idx;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE - In Progress
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Private constructor
    Node(const Graph* graph, size_type uid) : 
      graph_(const_cast<Graph*>(graph)), uid_(uid){

    }
    // graph_ points to the graph that the node belong to
    Graph* graph_;

    // Invariant Unique ID of the node, used to access the map contain
    // in graph that stores all the actual nodes. 
    size_type uid_; 

    // Note that the index of the node is stored in the actual graph

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0:
    return node_cnt;
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
    // HW0:
    
    // determine the new index, uid of the node
    size_type c_idx = size();
    size_type c_uid = nxt_node_uid_;

    // create the real_node structure
    real_node c_node;
    c_node.p = position;
    c_node.idx =  c_idx;

    // insert
    nodes_storage[c_uid] = c_node;
    node_idx_to_uid.push_back(nxt_node_uid_);

    // updateing 
    ++ node_cnt;
    ++ nxt_node_uid_;
    return Node(this, c_uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: WIP, consider if a std::Vector is added to track all uid in graph. 
    // Check if node points to the current graph
    return (this == n.graph_) && (nodes_storage.count(n.uid_) > 0);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: 
    return Node(this, node_idx_to_uid[i]);
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
      // HW0:
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0:
      size_type c_uid;
      if (rev){c_uid = (graph_->edges_storage)[uid_].n1_uid;}
      else{c_uid = (graph_->edges_storage)[uid_].n2_uid;}
      return Node(graph_,c_uid);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0:
      size_type c_uid;
      if (rev){c_uid = (graph_->edges_storage)[uid_].n2_uid;}
      else{c_uid = (graph_->edges_storage)[uid_].n1_uid;}
      return Node(graph_,c_uid);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Comparison 
      return 
      ((graph_ == e.graph_) && 
        (graph_->edges_storage[uid_].conn == 
          e.graph_->edges_storage[e.uid_].conn)
      );

    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: Use UID for comparison 
      return uid_ < e.uid_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge object

    // Construct a new edge object
    Edge(const Graph* graph, size_type uid, bool rev_flag) : 
      graph_(const_cast<Graph*>(graph)), uid_(uid), rev(rev_flag){

    }

    // if rev flag is not give, retieve it from the OG edge in the graph
    Edge(const Graph* graph, size_type uid) : 
      graph_(const_cast<Graph*>(graph)), uid_(uid), rev(){
        rev = ((const_cast<Graph*>(graph))->edges_storage)[uid].input_rev_fl;
    }

    
    // graph_ points to the graph that the node belong to
    Graph* graph_;

    // Invariant Unique ID of the edge
    size_type uid_; 

    // Mode for recording the input nodes uid order
    bool rev;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: 
    return edge_cnt;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: 
    size_type c_uid = edge_idx_to_uid[i];
    return Edge(this, c_uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0:
    std::string conn_name;
    if (a.uid_ < b.uid_) {
      conn_name = std::to_string(a.uid_) + " " + std::to_string(b.uid_);
    } else if (a.uid_ > b.uid_){
      conn_name = std::to_string(b.uid_) + " " + std::to_string(a.uid_);
    } else{
      return false;
    }

    return connected_nodes.find(conn_name) != connected_nodes.end();
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

    size_type c_uid;
    size_type c_n1_uid;
    size_type c_n2_uid;
    bool rev;

    std::string conn_name;
    if (a.uid_ < b.uid_) {
      conn_name = std::to_string(a.uid_) + " " + std::to_string(b.uid_);
      c_n1_uid = a.uid_;
      c_n2_uid = b.uid_;
      rev = true;
    } else if (a.uid_ > b.uid_){
      conn_name = std::to_string(b.uid_) + " " + std::to_string(a.uid_);
      c_n1_uid = b.uid_;
      c_n2_uid = a.uid_;
      rev = false;
    } else{
      assert(false);
    }

    if (connected_nodes.find(conn_name) != connected_nodes.end()){
      // edge already exists, retrive the uid_
      c_uid = connected_nodes[conn_name];
    } else {
      // HW0:
      // create a new edge and insert it 
      // determine the new index, uid
      size_type c_idx = edge_cnt;
      c_uid = nxt_edge_uid_;

      // create the real_edge structure
      real_edge c_edge;
      c_edge.n1_uid = c_n1_uid;
      c_edge.n2_uid = c_n2_uid;
      c_edge.idx =  c_idx;
      c_edge.conn = conn_name;
      c_edge.input_rev_fl = rev;

      // insert into Graph
      edges_storage[c_uid] = c_edge;
      edge_idx_to_uid.push_back(c_uid);
      connected_nodes[conn_name] = c_uid;

      // updating 
      ++ edge_cnt;
      ++ nxt_edge_uid_;
      
    }
    return Edge(this, c_uid, rev);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0:

    // clear nodes
    nodes_storage.clear();
    node_idx_to_uid.clear();
    connected_nodes.clear();

    // clear edges
    edges_storage.clear();
    edge_idx_to_uid.clear();

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

  // HW0: In Progress
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
      
  // Graph related data members


  // Node related data members 
  // real_node is the data structure for the actual nodes
  struct real_node {
    Point p;
    size_type idx; 
  };

  std::unordered_map <size_type, real_node> nodes_storage;
  size_type nxt_node_uid_;
  size_type node_cnt;
  std::vector<size_type> node_idx_to_uid;
  std::unordered_map<std::string, size_type> connected_nodes;

  // Edge related data members
  // real_edge WIP 
  struct real_edge{
    
    size_type n1_uid;
    size_type n2_uid;
    // increasing uid order, with " " in between
    std::string conn;
    size_type idx;
    bool input_rev_fl;
  };

  std::unordered_map<size_type, real_edge> edges_storage;
  size_type nxt_edge_uid_;
  size_type edge_cnt;
  std::vector<size_type> edge_idx_to_uid;
  

  
};

#endif // CME212_GRAPH_HPP
