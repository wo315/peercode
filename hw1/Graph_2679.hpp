#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
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

  struct internal_node;
  struct internal_edge;

 public:
   using node_value_type = V;

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
  Graph(): next_id_(0), nodes_(), edges_(), nodes_pointers_(), edges_pointers_() {}

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
      return fetch().position_;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      for (size_type i = 0; i < graph_->size(); ++i){
      	if (graph_->nodes_[i].uid_ == uid_){
          return i;
        }
      }
      assert(false);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    node_value_type& value(){
      return fetch().value_;
    }

    const node_value_type& value() const {
      return fetch().value_;
    }

    size_type degree() const {
      return graph_->edges_[uid_].size();
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_==n.graph_ && uid_==(n.uid_));
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
      return (graph_<n.graph_) || (graph_==n.graph_ && uid_<(n.uid_));
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_;

    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {}

    // Return the internal node associated to Node
    internal_node& fetch() const {
      return graph_->nodes_[uid_];
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_pointers_.size();
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
    // Define the uid_ of the new node using next_id_ (which is unused)
    // and incrementing this value
    size_type uid = next_id_;
    next_id_++;

    // Create a new internal node
    internal_node internal_node_(uid, position);

    // The internal node is added to the map for easy retrieval from the uid_
    nodes_.emplace(uid, internal_node_);

    // The internal node is added to the list of nodes
    nodes_pointers_.push_back(internal_node_);

    // Intialize the adjancency list of the new node
    std::map<size_type, internal_edge> adjacent_edges;
    edges_.emplace(uid, adjacent_edges);

    return Node(this, uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph_ && nodes_.find(n.uid_)!=nodes_.end());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, nodes_pointers_.at(i).uid_);
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
      return Node(graph_, node1_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_uid_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if(graph_!=e.graph_){
        return false;
      }
      size_type e1n1_uid = fetch().node1_uid_;
      size_type e1n2_uid = fetch().node2_uid_;
      size_type e2n1_uid = fetch().node1_uid_;
      size_type e2n2_uid = fetch().node2_uid_;
      return (e1n1_uid==e2n1_uid&&e1n2_uid==e2n2_uid)||(e1n1_uid==e2n2_uid&&e1n2_uid==e2n1_uid);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return node1()<e.node1() || (node1()==e.node1() && node2()<e.node2());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type node1_uid_;
    size_type node2_uid_;

    Edge(const Graph* graph, size_type node1_uid, size_type node2_uid)
        : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid) {}

    // Returns the internal edge associated to the edge
    internal_edge& fetch() const {
      return graph_->edges_.at(node1_uid_).at(node2_uid_);
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_pointers_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, edges_pointers_.at(i).node1_uid_, edges_pointers_.at(i).node2_uid_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Test in the adjancy lists
    // Firstly test if a is a node of the graph
    // Then tests if b is in the adjancy list of a
    return ((edges_.count(a.uid_)>0) && (edges_.at(a.uid_).count(b.uid_)>0));
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
    if(!has_edge(a, b)){ // Only add an edge if it is not in the graph
      // Create the associated internal edge
      internal_edge internal_edge_(a.uid_, b.uid_);
      // Add the edge to the adjancy lists of a and b
      edges_.at(a.uid_).emplace(b.uid_, internal_edge_);
      edges_.at(b.uid_).emplace(a.uid_, internal_edge_);

      // Add the edge to the list of edges
      edges_pointers_.push_back(internal_edge_);
    }
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    nodes_pointers_.clear();
    edges_pointers_.clear();
    next_id_ = 0;
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

    Node operator*() const {
      return graph_->node(index_);
    }

    NodeIterator& operator++(){
      index_++;
      return *this;
    }

    bool operator==(const NodeIterator& it) const{
      return (graph_==it.graph_)&&(index_==it.index_);
    }

   private:
    friend class Graph;
    // A NodeIterator is a graph and the index in the graph of the current node
    Graph* graph_;
    size_type index_;
    NodeIterator(const Graph* graph, size_type index)
    : graph_(const_cast<Graph*>(graph)), index_(index) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  node_iterator node_end() const {
    return NodeIterator(this, size());
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

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    Edge operator*() const {
      return Edge(graph_, node_uid_, internal_it_->first);
    }
    IncidentIterator& operator++() {
      index_++;
      ++internal_it_;
      return *this;
    }
    bool operator==(const IncidentIterator& it) const {
      return(graph_==it.graph_)&&(node_uid_==it.node_uid_)&&(index_==it.index_);
    }

   private:
    friend class Graph;
    // An incident IncidentIterator is a graph, the uid_ of the indicent node
    // the index of the current edge in the adjancy list and an internal iterator to the current edge
    Graph* graph_;
    size_type node_uid_;
    size_type index_; // from 0 to degree of the node
    typename std::map<size_type, internal_edge>::iterator internal_it_;

    IncidentIterator(const Graph* graph, size_type node_uid, size_type index)
    : graph_(const_cast<Graph*>(graph)), node_uid_(node_uid), index_(index), internal_it_() {
      internal_it_ = graph_->edges_[node_uid_].begin();
      for(size_type i=0; i<index_; i++) {
        internal_it_++;
      }
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

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    Edge operator*() const {
      return graph_->edge(index_);
    }

    EdgeIterator& operator++(){
      index_++;
      return *this;
    }

    bool operator==(const EdgeIterator& it) const{
      return (graph_== it.graph_)&&(index_== it.index_);
    }

   private:
    friend class Graph;
    // An EdgeIterator is a graph and the index in the graph of the current edge
    Graph* graph_;
    size_type index_;
    EdgeIterator(const Graph* graph, size_type index)
    : graph_(const_cast<Graph*>(graph)), index_(index) {}
  };

  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  edge_iterator edge_begin(){
    return EdgeIterator(this, 0);
  }

  edge_iterator edge_end(){
    return EdgeIterator(this, num_edges());
  }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node {
    // Members
    size_type uid_;
    Point position_;
    node_value_type value_;
    // Constructor
    internal_node() {}
    internal_node(size_type uid, Point position):
    uid_(uid), position_(position), value_(0) {}
    internal_node(size_type uid, Point position, node_value_type value):
    uid_(uid), position_(position), value_(value) {}
  };

  struct internal_edge {
    // Members
    size_type node1_uid_;
    size_type node2_uid_;
    // Constructor
    internal_edge(size_type node1_uid, size_type node2_uid):
    node1_uid_(node1_uid), node2_uid_(node2_uid) {}
  };

  size_type next_id_; // Store and increment to ensure each node as a unique uid_
  // nodes and edges are stored using adjancy lists
  // This allows testing has_node/has_edge easily
  // This is useful for incident iteration
  std::map<size_type, internal_node> nodes_;
  std::map<size_type, std::map<size_type, internal_edge>> edges_;
  // nodes_pointers_ and edges_pointers_ store a nodes and edges in a vector
  // They ensure easy access to node(size_type i) and edge(size_type i)
  std::vector<internal_node> nodes_pointers_;
  std::vector<internal_edge> edges_pointers_;


};

#endif // CME212_GRAPH_HPP
