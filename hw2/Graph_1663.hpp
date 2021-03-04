#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
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

  /** value type of node */
  using node_value_type = V;
  /** value type of edge */
  using edge_value_type = E;

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
  using size_type = uint32_t;

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
  class Node : private totally_ordered<Node>{
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
    }

    /** Return this node's position. */
    const Point& position() const {
      return  std::get<0>( graph_->nodes.at( index_ ) ) ;
      

    }

    Point& position() {
      return std::get<0>( graph_->nodes.at( index_ ) );
    }



    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodemapper.at(index_);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** @return Return this node's value */
    const node_value_type& value() const {
      return std::get<1>( graph_->nodes.at( index_ ) );
      
    }
    /**@return Returns this node value as non-const */
    node_value_type& value() {
      return std::get<1>( graph_->nodes.at( index_ ) );

    }

    /**
    @param[in] newval  the new node value
    A setter for a node's value in case it wasn't declared as const
    */

    void set_value(node_value_type newval) {
      graph_->nodes.at( index_ ).second = newval;
    }
    /**
    @return the degree of this node
    */
    size_type degree() const {
      if( graph_->edges.count(index_) == 0) return 0;
      return graph_->edges.at(index_).size();
    }

    /**
    @return returns an iterator pointing to the first of the edes incident to this node
    */
    incident_iterator edge_begin() const {
      return IncidentIterator( graph_, index_, graph_->edges[index_].cbegin() );
    }

     /**
    @return returns an iterator pointing to the end of the edes incident to this node
    */
    incident_iterator edge_end() const {
      return IncidentIterator( graph_, index_, graph_->edges[index_].cend() );
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE

      return ( n.index_ == index_ ) && ( graph_ == n.graph_ );
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

      return ( index_ < n.index_ );
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    Graph* graph_ = nullptr;
    size_type index_ = 0;
    Node( const Graph* g, size_type index )
            : graph_( const_cast<Graph*>( g ) ), index_( index ) {
        }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
   
    return nodes.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] val The value to assign the node too
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& val = node_value_type() ) {
    // HW0: YOUR CODE HERE

    nodes[next_index] = std::make_pair( position, val );
    node_uid.push_back(next_index);
    nodemapper[next_index] = node_uid.size() - 1;
    return Node( this, next_index++ );        // Invalid node
  }
  /** Removes node if it exists

  @post has_node( @a n) = false
  @post for all nodes j != @a n, has_edge(n, j) = false
  @post all node and edge iterators outstanding are invalidated
  complexity: O(I), where I is the number of edges incident to @a n
  @return 1 if node was removed, 0 if not ( node does not exist or is invalid)
  */
  size_type remove_node( const Node& n) {
    if( !has_node(n) ) return 0;
    
    for(auto it = n.edge_begin(); it != n.edge_end(); it = n.edge_begin()) {
      remove_edge((*it));
    }
    
    nodes.erase( n.index_ );
    size_type loc = nodemapper.at(n.index_);
    bool flag = loc == node_uid.size() - 1;
    node_uid[loc] = node_uid[node_uid.size() - 1];
    node_uid.pop_back();
    nodemapper.erase(n.index_);
    if(!flag)
      nodemapper.at(node_uid[loc]) = loc;

    return 1;
  }

  /** Removes node if it exists

  @post has_node( @a (*n_it)) = false
  @post for all nodes j != @a *n_it, has_edge(*n_it, j) = false
  @post all node and edge iterators outstanding are invalidated
  complexity: O(I), where I is the number of edges incident to @a *n_it
  @return 1 if node was removed, 0 if not ( node does not exist or is invalid)
  */
  node_iterator remove_node( node_iterator n_it ) {
    remove_node( *n_it );
    return node_begin();
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE

    return ( n.graph_ == this ) && ( nodes.count( n.index_ ) );
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE

    return Node( this, node_uid.at(i) );      
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
  class Edge : private totally_ordered<Node>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HEREN
      return Node( graph_, std::get<0>( graph_->edge_id.at( index_ ) ) );      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node( graph_, std::get<1>( graph_->edge_id.at( index_ ) ) ); ;     
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      //HW0: YOUR CODE HERE
      return e.graph_ == graph_ && index_ == e.index_;
    }
    /**
    Tests inequality between this and another edge
    */

    bool operator!=(const Edge&e ) const {
      return !(*this == e);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {

      //HW0: YOUR CODE HERE
      return (index_ < e.index_ && e.graph_ == graph_) || (e.graph_ != graph_ && index_ <= e.index_);
    }
    bool operator>(const Edge& e) const {
       return (index_ > e.index_ && e.graph_ == graph_) || (e.graph_ != graph_ && index_ > e.index_);
    }

    /**
    @return value associated with this edge as const
    */ 
    const edge_value_type& value() const {
      return graph_->edgedat.at( index_ );
    }

    /**
    @return value associated with this edge as modifiable
    */ 
    edge_value_type& value() {
      return graph_->edgedat.at( index_ );
    }

    void set_value( edge_value_type val ) {
      graph_->edgedat[index_] = val;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_ = nullptr;
    size_type index_ = 0;
    Edge( const Graph* g, size_type index )
            : graph_( const_cast<Graph*>( g ) ), index_( index ) {
        }
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
    return edge_id.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE

    return Edge( this, edge_uid.at(i) );        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    return edges.count( a.index_ ) && edges.at( a.index_ ).count( b.index_ );
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type val = edge_value_type()) {
    // HW0: YOUR CODE HERE

    size_type x, y;
    x = a.index_;
    y = b.index_;
    edgedat[next_edge] = val;
    
    if( has_edge( a, b ) ) {
      size_type eid = edges.at(x).at(y);
      if( std::get<0>(edge_id.at(eid)) == y ) edge_id[eid] = std::make_pair(x, y); 
      return Edge( this, eid );
      
    }
    edge_uid.push_back(next_edge);
    edgemapper[next_edge] = edge_uid.size() - 1;
    edges[x][y] = next_edge;
    edges[y][x] = next_edge;
    edge_id[next_edge] = std::make_pair( x, y );
    return Edge( this, next_edge++ );        // Invalid Edge
  }

  /**
  Removes the edge determined by the two nodes @a a and @a b.  
  @post has_edge(@a a, @a b) = false
  @post all outstanding iterators are invalidated. Any edge object x such that x is equivalent to the edge deleted is invalidated
  Complexity: O(1) expected time
  @post all external ids of edges might have been changed (eg. an edge accessed with edge(x) might not be the same as the one found after an erase operation)
  @return 1 if edge successfully removed, 0 if not (or edge does not exist)
  */
  size_type remove_edge(const Node& a, const Node& b) {

    if( !has_edge( a, b ) ) return 0;
    size_type id = edges.at(a.index_).at(b.index_);
    edges[a.index_].erase(b.index_);
    edges[b.index_].erase(a.index_);
    edgedat.erase(id);
    edge_id.erase(id);
    size_type loc = edgemapper.at(id);
    bool flag = loc == edge_uid.size() - 1;
    edge_uid[loc] = edge_uid[edge_uid.size()-1];
    edge_uid.pop_back();
    edgemapper.erase(id);
    if(!flag)
      edgemapper.at(edge_uid[loc]) = loc;
    return 1;
  }

  /**
  Removes the edge determined by @a e  
  @post has_edge(@a e.node1(), @a e.node2()) = false
  @post all outstanding iterators are invalidated. Any edge object x such that x is equivalent to the edge deleted is invalidated
  Complexity: O(1) expected time
  @post all external ids of edges might have been changed (eg. an edge accessed with edge(x) might not be the same as the one found after an erase operation)
  @return 1 if edge successfully removed, 0 if not (or edge does not exist)
  */ 
  size_type remove_edge(const Edge& e) {
    return remove_edge( e.node1(), e.node2() );
  }

  /**
  Removes the edge pointed to by edge iterator @a e_it  
  @post has_edge(@a (*e_it).node1(), @a (*e_it).node2()) = false
  @post all outstanding iterators are invalidated. Any edge object x such that x is equivalent to the edge deleted is invalidated
  Complexity: O(1) expected time
  @post all external ids of edges might have been changed (eg. an edge accessed with edge(x) might not be the same as the one found after an erase operation)
  @return an edge iterator to the beginning of this graph's edges
  */ 
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge((*e_it).node1(), (*e_it).node2() );
    return edge_begin();
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edges.clear();
    nodes.clear();
    edge_id.clear();
    edge_uid.clear();
    node_uid.clear();
    edgemapper.clear();
    nodemapper.clear();
    edgedat.clear();
    next_index = 0;
    next_edge = 0;
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
    /**
    Return the node pointed too by this iterator
    */
    Node operator*() const {
      return Node( graph_, (*iter).first );
    }
    /**
    Compares equality with this and another NodeIterator
    */

    bool operator==(const NodeIterator& other) const {
      return this->graph_ == other.graph_ && this->iter == other.iter;
    }
    /**
    Increments the iterator to point to the next element
    */
    NodeIterator& operator++() {
      iter++;
      return *this;
    }

    /**
    Checks inequality with this and another NodeIterator
    */
    bool operator!=(const NodeIterator& other ) const {
      return !(*this == other);
    }




   private:
    friend class Graph;
    Graph* graph_ = nullptr;
    typename std::unordered_map<size_type, std::pair<Point, node_value_type>>::const_iterator iter;
    NodeIterator( const Graph* g, typename std::unordered_map<size_type, std::pair<Point, node_value_type>>::const_iterator it)
            : graph_( const_cast<Graph*>( g ) ), iter( it ) {

        }
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  /**
    Returns an iterator to the beginning of the node
    @return the begin iterator
  */
  node_iterator node_begin() const {
    return NodeIterator( this, nodes.cbegin() );

  }

  /**
  creates an iterator pointing to the end of the nodes
  @return the end iterator
  */

  node_iterator node_end() const {
    return NodeIterator( this, nodes.cend() );

  }

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

    /**
    dereferences an Incident Iterator to return the current edge
    */
    Edge operator*() const {
      size_type id = (*iter).second;
      Edge e = Edge( graph_, id );
      if( e.node1().index_ != rootid) graph_->edge_id[id] = std::make_pair( e.node2().index_, e.node1().index_ );

      return e;
    }

    /**
    Checks equality with another IncidentIterator
    */
    bool operator==(const IncidentIterator& other ) const {
      return (this->graph_ == other.graph_ && this->iter == other.iter);
    }

    /**
    Checks inequality with another IncidentIterator
    */

    bool operator!=(const IncidentIterator& other ) const {
      return !(*this == other);
    }

    /**
    increments the IncidentIterator to point to next edge
    */
    IncidentIterator& operator++() {
      iter++;
      return *this;

    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type rootid;
    std::unordered_map<size_type, size_type>::const_iterator iter;
 
    IncidentIterator(const Graph* g, size_type root, typename std::unordered_map<size_type, size_type>::const_iterator it)
        : graph_( const_cast<Graph*>( g )), rootid( root), iter( it ) {

        }
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
    /** 
      @return Returns the edge this iterator currently points to
    */
    Edge operator*() const {
      return Edge( graph_, (*iter).first );
    }
    /**
      Checks to see if this edge iterator is equivalent to another one
      @return true if this is equal to @a other, otherwise false
    */
    bool operator==(const EdgeIterator& other ) const {
      return graph_ == other.graph_ && iter == other.iter;
    }

    EdgeIterator& operator++() {
      iter++;
      return *this;
    }
    /**
    @param other the edgeiterator to compare to
    @return true if this is not equal to the @a other iterator, false otherwise
    */ 
    bool operator!=(const EdgeIterator& other ) const {
      return !(*this == other);
    }



   private:
    friend class Graph;
    Graph* graph_;
    std::unordered_map<size_type, std::pair<size_type, size_type>>::const_iterator iter;
    EdgeIterator( const Graph* g, typename std::unordered_map<size_type, std::pair<size_type, size_type>>::const_iterator it )
        : graph_( const_cast<Graph*>( g ) ), iter( it) {

        }

    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  /**
    @return Returns an beginning iterator over all edges in the graph
  */
  edge_iterator edge_begin() const {
    return EdgeIterator( this, edge_id.cbegin() );
  }

  /**
   @return  Returns an end iterator over all edges in the graph
  */
  edge_iterator edge_end() const {
    return EdgeIterator( this, edge_id.cend() );
  }

 private:
  std::unordered_map<size_type, std::pair<Point, node_value_type>> nodes;
  size_type next_index = 0;
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> edges;
  std::unordered_map<size_type, std::pair<size_type, size_type>> edge_id;
  std::unordered_map<size_type, edge_value_type> edgedat;
  std::vector<size_type> node_uid;
  std::vector<size_type> edge_uid;
  std::unordered_map<size_type, size_type> edgemapper;
  std::unordered_map<size_type, size_type> nodemapper;

  size_type next_edge = 0;

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
