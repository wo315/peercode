#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <forward_list>
#include <unordered_map>
#include <vector>

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
     *   x = graph.node(0);  this node is valid
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
      // HW0: YOUR CODE HERE

      this->node_parent_graph_ = nullptr; // invalid nodes have no parent graph
    }

    // Public Copy construtor: whenever we copy a node, we keep track of the
    // copy's address with the list node_roster.  This way, we can invalidate
    // each copy when we clear the graph through setting the parent_graph
    // pointer to null
    Node(const Node &n) {
        (n.node_parent_graph_)-> node_roster.emplace_front(this);

        this->node_parent_graph_ = n.node_parent_graph_;
        this->node_index_ = n.node_index_;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      
      return (this->node_parent_graph_)->node_map.at(this->node_index_);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      
      return (this->node_index_);
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
      
      if ((this->node_parent_graph_ == n.node_parent_graph_) && 
           (this->node_index_ == n.node_index_)){
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
      // we will compare nodes by their index, assuming they have the same
      // parent graph
      if ((this->node_index_ < n.node_index_)&&
         (this->node_parent_graph_ == n.node_parent_graph_)){
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

    graph_type* node_parent_graph_; // pointer to the node's parent graph
    size_type node_index_;
    

    //private constructor
   
    Node(const graph_type* parent_graph, size_type idx)
        : node_parent_graph_(const_cast<graph_type*>(parent_graph)),
          node_index_(idx)
          {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE

    // the graph's node map is defined by the injective mapping 
    // node_index -> node_position.  The number of nodes is therefore
    // equal to the size of the map.  The size of an unordered_map
    // can be retrieved in O(1) time.

    return node_map.size();
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

    //node_map is a hash map, so adding a term is O(1) in complexity
    node_map.emplace(num_nodes(), position);
    
    // return a node that holds all the relevant information
    return Node(this, num_nodes() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.node_parent_graph_ == this &&
        (n.node_parent_graph_)->size() > n.node_index_){
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

    assert(i < size());
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      this->edge_parent_graph_ = nullptr; // invalid edges have no parent graph
    }

    // public copy constructor

    Edge(const Edge &e) {
        (e.edge_parent_graph_)->
        edge_roster.emplace_front(this);

        this->edge_parent_graph_ = e.edge_parent_graph_;
        this->edge_index_ = e.edge_index_;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE

      // when an edge connecting nodes n1 and n2 is added to a graph, the pair
      // (n1, n2) as recorded in edge_map always satisfies n1 < n2.
      return Node(this->edge_parent_graph_, 
             (this->edge_parent_graph_)->edge_map.at(this->edge_index_)[0]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->edge_parent_graph_, 
             (this->edge_parent_graph_)->edge_map.at(this->edge_index_)[1]);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if ((this->edge_parent_graph_ == e.edge_parent_graph_) && 
         ((this->node1() == e.node1()) && (this->node2() == e.node2()))){
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
      // we will compare based on edge index
      if (this->edge_parent_graph_ == e.edge_parent_graph_ && 
          this->edge_index_ < e.edge_index_){
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

    // The Edge proxy class will use a pointer to the parent graph along
    // with an unsigned integer representing the edge index in order to 
    // construct an interface with the edge of a graph
    graph_type* edge_parent_graph_;
    size_type edge_index_;

    // private constructor

    Edge(const graph_type* parent_graph, size_type idx)
    : edge_parent_graph_(const_cast<graph_type*>(parent_graph)),
      edge_index_(idx)
      {}
  
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_map.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE

    assert(i < num_edges());
    return Edge(this, i);        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    // the hash map inv_edge_map stores already sorted nodes;
    // complexity is therefore O(1)
  
  if (inv_edge_map.count(std::vector<size_type>{// see if container is nonempty
                         std::min(a,b).node_index_,
                         std::max(a,b).node_index_}) > 0){
        
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
    assert(-(a == b) + 1);                // edges must be distinct
    if (has_edge(a, b)){                  // check if edge is already present
        return Edge(this, inv_edge_map.at(std::vector<size_type>{
                                          std::min(a,b).node_index_,
                                          std::max(a,b).node_index_}));
        }
    else{                            // if the edge is not present, we add it
        edge_map.emplace(num_edges(), std::vector<size_type>{
                                      std::min(a,b).node_index_,
                                      std::max(a,b).node_index_});

        //num_edges() has already increased from using emplace on edge_map
        //the output value of the reversed map is therefore num_edges()-1
        inv_edge_map.emplace(std::vector<size_type>{
                             std::min(a,b).node_index_,
                             std::max(a,b).node_index_}, 
                                              num_edges()-1);
        
        return Edge(this, num_edges() - 1); // return the added edge
        }                         
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

    // we first clear out both node_map, edge_map and inv_edge_map
    node_map.clear();
    edge_map.clear();
    inv_edge_map.clear();

    for (auto const& i : node_roster) {
    (*i).node_parent_graph_ = nullptr;
    }

    node_roster.clear();

    for (auto const& j : edge_roster) {
    (*j).edge_parent_graph_ = nullptr;
    }
  
    edge_roster.clear();
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
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.


  // our map class stores node and edge information with three hash maps
  std::unordered_map< size_type, Point > node_map; 
  // maps node_index -> position

  std::unordered_map< size_type, std::vector< size_type > > edge_map;
  // maps edge_index -> <node1_index, node2_index>, where node1 < node2

  // given an edge e, we will use another hash function so that accessing the
  // nodes spanned by e is O(1) in complexity.  Following Lippman (p. 874), we
  // will define a hash map <node1_index, node2_index> -> edge_index.  Since
  // there is no native hash function for vector<size_type>, we will create
  // our own

  struct inv_edge_hasher {
    size_type operator()(const std::vector<size_type> &v) const {
       return std::hash< size_type>()(v[0])^(std::hash< size_type>()(v[1])<<1);
    }
  };

  std::unordered_map< const std::vector<size_type>, size_type, inv_edge_hasher>
                                                                 inv_edge_map;
  

  // forward lists maintain the addresses of nodes/edges created from
  // our copy constructor, so that their parent_graph pointer attributes
  // can be set to nullptr when the graph is cleared, thereby invalidating
  // these objects.
   
  std::forward_list<Node*> node_roster;
  std::forward_list<Edge*> edge_roster;


};

#endif // CME212_GRAPH_HPP
