#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>
#include <set>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

// https://piazza.com/class/kjx6t0a2rv4589?cid=175
template <typename V = int>
class Graph {

  public:
    using node_value_type = V;
    using size_type = unsigned;
    class Node;
    class Edge;

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  //later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)


  // three vectors that function like maps. Index links to points to follow
  // a node for point_vector. Index links to internal node copies for 
  // node() calls, and index links to internal edge copies for edge() calls

  std::vector<Point> point_vector;
  std::vector<Graph::Node> node_vector;
  std::vector<Graph::Edge> edge_vector;
  std::vector<std::vector<unsigned>> incident_vector;
  std::vector<V> values_vector;

  // map that stores two node indices as key and edge index as values for
  // has_edge and add_edge

  std::map<std::vector<size_type>, size_type> node_to_edge_map;
  std::set<unsigned> edge_set;

  // store number of nodes and number of edges as graph attributes
  size_type size_;
  size_type edge_size_;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */

  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : size_(0), edge_size_(0) {
    // HW0: YOUR CODE HERE
    // sets the number of nodes and edges of an empty graph to 0 at init
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
    Node()
      : graph_(nullptr) {
      // HW0: YOUR CODE HERE
      // sets graph_ attribute to make clear this is invalid
    }

    /** Return this node's position by accessing the point in internal
        vector
    */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->point_vector[idx_];
    }

    /** Return this node's index, a number in the range [0, graph_size)
        which is stored as unsigned int within node, first checking it
        is in range
    */
    size_type index() const {
      assert(idx_<graph_->size_);
      return idx_;
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
      return idx_ == n.idx_ && graph_ == n.graph_;
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
      // to maintain trichotomy for nodes in different graphs, I use logic
      // that if nodes are of different graphs, the one with a lower memory
      // address will be lesser. If they are the same graph, then it just
      // checks the indices
      if(graph_==n.graph_) {
        return idx_<n.idx_;
      }
      else {
        return graph_<n.graph_;
      }
    }

    /** 
     * @brief Return a reference to a nodes value 
     * @return reference to a node's value that can be modified
     */
    node_value_type& value() {
      return graph_->values_vector[idx_];
    }

    /** 
     * @brief Return const reference to a nodes value 
     * @return reference to a node's value that can't be modified
     */
    const node_value_type& value() const {
      return graph_->values_vector[idx_];
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    //Private pointer to graph and an unsigned int for the node's index
    Graph* graph_;
    size_type idx_;

    // private constructor for a valid node, note it returns a node to the 
    // user and stores an internal copy in the node_vector for use in the 
    // node() method

    Node(const Graph* graph, size_type idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx) {}

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    public:

    /** 
     * @brief Return how many edges a node has 
     * @return Unsigned int with number of edges
     */
    size_type degree() const {
      assert(idx_<graph_->incident_vector.size());
      return graph_->incident_vector[idx_].size();
    }

    /** 
     * @brief Create a begin IncidentIterator pointing to the first edge of
     * a node 
     * @return starting incident_iterator
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_,idx_, 0);
    }

    /** 
     * @brief Create an end IncidentIterator pointing past last edge of
     * a node 
     * @return ending IncidentIterator
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_,idx_, degree());
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
   * @param[in] value for a node, defaults to 0
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */

  Node add_node(const Point& position,\
    const node_value_type& initial_value = node_value_type()) {
    // HW0: YOUR CODE HERE
    // increments number of nodes counter, then adds the point to
    // point vector so that position() can find the point associated with
    // the node, then constructs the node for user and internally node
    // using the private constructor
    ++size_;
    point_vector.push_back(position);
    incident_vector.push_back({});
    Node temp_node = Node(this, size_-1);
    node_vector.push_back(temp_node);
    values_vector.push_back(initial_value);
    return temp_node;        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // checks if n's graph points to this graph
    // Guillermo suggested this in OH to get to O(1)
    return n.graph_==this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // finds the node in node_vector, first making sure a valid index is
    // given
    assert(i<size_); 
    // http://www.cplusplus.com/reference/map/map/at/
    return node_vector[i];        // Invalid node
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

    /* I wrote this for testing but wasn't sure if new public methods were OK
     * so commented it out
     *  size_type index() {
     *    return edge_idx_;
     *  }
     */
    Edge() 
      : graph_(nullptr) {
      // HW0: YOUR CODE HERE
      // initialize to nullptr to ensure invalid
      
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // provides node1's index to node() to get node out of node_vector
      return graph_->node(idx1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // provides node2's index to graph.node() to get node
      return graph_->node(idx2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      //Because constructor ensures edge_idx_ is unique for either direction
      // of an edge, this can just check that idx and graph are same
      return edge_idx_ == e.edge_idx_ && graph_ == e.graph_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // to ensure trichotemy even if edges are for different graphs, it
      // checks index only when graph is same. If graph is different, lower
      // memory address is "less than"
      if(graph_==e.graph_) {
        return edge_idx_ < e.edge_idx_;
      }
      else {
        return graph_ < e.graph_;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    //private attributes of pointer to graph, node1's idx, node2's idx,
    //and an index for the edge itself

    Graph* graph_;
    size_type idx1_;
    size_type idx2_;
    size_type edge_idx_;

    // this ordinary constructor is called when has_edge has returned false
    // so we need a new edge. Note that add_edge first makes sure a<b 
    // so that idx1 for any edge < idx2. This is used so that duplicate
    // nodes cannot be added

    // Constructor adds entry to node_to_edge_map with node1+node2 as key
    // so that has_edge and add_edge can find edge from nodes, since
    // edge_vector only uses edge_idx as the key

    // Similar to node, a copy of the returned edge is stored in
    // edge_vector for the edge() method

    Edge(Graph* graph, size_type idx1, size_type idx2, size_type edge_idx)
      : graph_(const_cast<Graph*>(graph)), idx1_(idx1), idx2_(idx2),\
      edge_idx_(edge_idx) {}

    // this ordinary constructor is called if node1
    // and node2 were provided in the wrong
    // order as an edge already stored.
    // This gives an edge to the user with the same edge_idx as the
    // "flipped" edge but this edge has node1 and node2 in the order the user
    // requested.
    // This constructor
    // does not add anything internally or change any graph attributes
    Edge(Graph* graph, size_type idx1, size_type idx2)
       : graph_(const_cast<Graph*>(graph)), idx1_(idx1), idx2_(idx2),\
         edge_idx_(graph_->node_to_edge_map[{idx2,idx1}]) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // finds the internal edge stored in edge_vector
    return edge_vector[i];        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Because we are given two nodes, we check node_to_edge_map if these
    // nodes match an edge_idx in either order (so will work even if user
    // is using an edge that's flipped relative to node_to_edge_map)
    int edge_count = 0;
    if (a<b) {
      edge_count+=node_to_edge_map.count({a.idx_,b.idx_});
    }
    if (b<a) {
      edge_count+=node_to_edge_map.count({b.idx_,a.idx_});
    }
    // using the int allows this to return false if a==b
    return edge_count;
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

    // if the graph has the edge (in either order), we return it. If they
    // provided a and b in the order that matches node_to_edge_map and
    // edge_vector, we can just grab it. Or we
    // provide them an equivalent edge with node1 and node2 flipped if they
    // asked for the reverse order. Note in second case, we call
    // constructor that does not save this "flipped" edge internally

    // Thanks to Maoguo who kindly helped me think through this logic

    if(has_edge(a,b)) {
      if(a<b) {     
        return edge_vector[node_to_edge_map.at({a.idx_,b.idx_})];
      }
      if(b<a) {
        return Edge(this, a.idx_, b.idx_);
      }
      else {
        return Edge(); // quiet compiler warning but will never be called
      }
    }
    // If the edge does not already exist
    // We call the constructor, then store the new edge
    // new edge internally, checking if we need to flip it first before storing
    // then returning the version they asked for

    else {
      ++edge_size_;
      incident_vector[a.idx_].push_back(edge_size_-1);
      incident_vector[b.idx_].push_back(edge_size_-1);
      Edge temp_edge = Edge(this, a.idx_, b.idx_, edge_size_-1); 
      if (a<b) {
        node_to_edge_map[{a.idx_,b.idx_}]=edge_size_-1;
        edge_vector.push_back(temp_edge);
      }
      if (a>b) {
        Edge flipped_temp = Edge(this, b.idx_, a.idx_, edge_size_-1);
        node_to_edge_map[{b.idx_,a.idx_}]=edge_size_-1;
        edge_vector.push_back(flipped_temp);
      }
      return temp_edge;
    }
  }
 

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // Deletes all internal storage and sets node and edge size to 0, thus
    // invalidating any outside Node and Edges
    // In OH we discussed we don't need to set all client nodes and edge to
    // nullptr, but that this may be implemented in a future HW

    // http://www.cplusplus.com/reference/map/map/clear/
    edge_set.clear();
    values_vector.clear();
    incident_vector.clear();
    point_vector.clear();
    node_vector.clear();
    edge_vector.clear();
    node_to_edge_map.clear();
    size_=0;
    edge_size_=0;
  }


  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: public std::iterator<std::forward_iterator_tag,Node> {
//  class NodeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                    // Element type
    using pointer           = const Node*;             // Pointers to elements
    using reference         = Node&;                   // Reference to elements
    using difference_type   = std::ptrdiff_t;          // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /** Increment the NodeIterator to point to the next node
     * @post ptr_ points to the following node in node_vector
     *
     * @return The new NodeIterator with incremented ptr_
     */
    NodeIterator& operator++() {
      ptr_++;
      return *this;
    }

    /** Dereference the NodeIterator to get the node it points to
     *
     * @return the node the iterator currently points to
     */
    value_type operator*() const {
       return *ptr_;
    }

    /** Check equality against another NodeIterator
     * @param[in] other_node another NodeIterator
     *
     * @return true if both ptr_ attributes point to same Node
     */
    bool operator==(const NodeIterator& other_node) const {
      return ptr_==other_node.ptr_;
    }

    /** Check inequality against another NodeIterator 
     * @param[in] other_node another NodeIterator
     *
     * @return true if ptr_ attributes point to different Nodes
     */
    bool operator!=(const NodeIterator& other_node) const {
      return ptr_!=other_node.ptr_;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    pointer ptr_;

    /** 
     * @brief Constructor for NodeIterator
     * @param[in] ptr to a Node
     * @post Creates a NodeIterator pointing to that node
     */
    NodeIterator(pointer ptr) : ptr_(ptr) {}
  };

  /** 
    * @brief Create a begin NodeIterator that points to the first node in
    * node_vector
    * @return NodeIterator pointing to first node in a graph
   */
  NodeIterator node_begin() const {
    const Node* node_ptr = &node_vector[0];
    return NodeIterator(node_ptr);
  }

  /** 
    * @brief Create an end NodeIterator that points one past the last node in
    * node_vector
    * @return NodeIterator pointing past last node in a graph 
   */
  NodeIterator node_end() const {
    const Node* node_ptr = &node_vector.back()+1;
    return NodeIterator(node_ptr);
  }
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

    /** Dereference the IncidentIterator to get the Edge it points to
     *
     * @return the edge the iterator currently points to, being sure to
     * orient the edge such that node1 == the node that created this iterator
     */
    value_type operator*() const {
      // if node1 of the edge is the node that created this iterator, return
      // the edge in edge_vector referenced by incident_vector
      if(graph_->edge(graph_->incident_vector[node_idx_][edge_number_])\
        .idx1_==node_idx_) {
        return graph_->edge(graph_->incident_vector[node_idx_][edge_number_]);
      }
      // if the edge is flipped, create a new Edge with the same index as the
      // edge stored internally
      else {
        return Edge(graph_, graph_->edge(graph_->incident_vector\
          [node_idx_][edge_number_]).idx2_, graph_->edge(graph_->\
          incident_vector[node_idx_][edge_number_]).idx1_);
      }
    }

    /** Increment the IncidentIterator to point to the next incident edge
     * @post ptr_ points to the following node in incident_vector
     *
     * @return The new IncidentIterator with incremented edge_number
     */
    IncidentIterator& operator++() {
      edge_number_++;
      return *this;
    }

    /** Check equality against another IncidentIterator
     * @param[in] other_iterator another IncidentIterator
     * @pre both incident iterators are valid and point to an edge
     * @return true if both node and edge attributes match, meaning it is the
     * same edge
     */
    bool operator==(const IncidentIterator& other_iterator) const {
      return node_idx_== other_iterator.node_idx_ && edge_number_\
        == other_iterator.edge_number_;
    }

    /** Check inequality against another IncidentIterator
     * @param[in] other_iterator another IncidentIterator
     * @ pre both iterators are valid and point to an edge
     * @return true if both node and edge attributes don't match, meaning
     * it is a differend edge
     */
    bool operator!=(const IncidentIterator& other_iterator) const {
      return node_idx_!= other_iterator.node_idx_ || edge_number_ !=\
        other_iterator.edge_number_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    unsigned node_idx_;
    unsigned edge_number_;
    // HW1 #3: YOUR CODE HERE

    /** 
     * @brief Constructor for IncidentIterator
     * @param[in] Pointer to a graph object
     * @param[in] index of the node that created this iterator
     * @param[in] edge_number to indicate where in incident_vector to point
     * @post Creates an IncidentIterator pointing to that edge
     */
    IncidentIterator(Graph* graph, unsigned idx, unsigned edge_number) :\
      graph_(graph), node_idx_(idx), edge_number_(edge_number) {}

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
    using pointer           = const Edge*;              // Pointers to elements
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

    /** Dereference the EdgeIterator to get the edge it points to
     *
     * @return the edge the iterator currently points to
     */
    value_type operator*() {
      return *ptr_;
    }

    /** Increment the EdgeIterator to point to the next edge
     * @post ptr_ points to the following edge in edge_vector
     *
     * @return The new EdgeIterator with incremented ptr_
     */
    EdgeIterator& operator++() {
      ptr_++;
      return *this;
    }

    /** Check equality against another EdgeIterator
     * @param[in] other_iter another IncidentIterator
     * @pre both incident iterators are valid and point to an edge
     * @return true if pointers point to the same edge
     * same edge
     */
    bool operator==(const EdgeIterator& other_iter) {
      return ptr_==other_iter.ptr_;
    }

    /** Check inequality against another EdgeIterator 
     * @param[in] other_iter another EdgeIterator
     *
     * @return true if ptr_ attributes point to different edges
     */
    bool operator!=(const EdgeIterator& other_iter) {
      return ptr_!=other_iter.ptr_;
    }

   private:
    friend class Graph;
    pointer ptr_;
    
    /** 
     * @brief Constructor for IncidentIterator
     * @param[in] ptr Pointer to an edge
     * @post Creates an EdgeIterator pointing to that edge
     */
    EdgeIterator(pointer ptr) : ptr_(ptr) {}
  };
    
  /** 
    * @brief Create a begin EdgeIterator that points to the first edge in
    * edge_vector
    * @return EdgeIterator pointing to first edge in a graph
   */
  EdgeIterator edge_begin() const {
    const Edge* edge_ptr = &edge_vector[0];
    return EdgeIterator(edge_ptr);
  }

  /** 
    * @brief Create an EdgeIterator that points one past last edge in
    * edge_vector
    * @return EdgeIterator pointing past last edge in a graph
   */
  EdgeIterator edge_end() const {
    const Edge* edge_ptr = &edge_vector.back()+1;
    return EdgeIterator(edge_ptr);
  }


  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
