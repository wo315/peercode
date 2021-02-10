#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  using node_value_type = V;

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

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  std::unordered_map<size_type, Point> nodes;
  std::unordered_map<size_type, node_value_type> node_value;
  std::unordered_map<size_type, std::pair<size_type,size_type> > edges;
  std::unordered_map<size_type, std::unordered_map<size_type, size_type> > node_edge;

 public:
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
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_n->nodes[idx_n];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_n;
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
      return (graph_n==n.graph_n && idx_n==n.index());         
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
      return (graph_n==n.graph_n && idx_n<n.index());
    }

    node_value_type& value(){
       return graph_n->node_value[idx_n];
    }

    const node_value_type& value() const{
       return graph_n->node_value[idx_n];
    }

    size_type degree() const{
      IncidentIterator it;
      size_type degree=0;
      for(it=edge_begin(); it!=edge_end(); ++it){
        degree+=1;
      }
      return degree;
    }

    IncidentIterator edge_begin () const{
      return IncidentIterator(graph_n, idx_n, (graph_n->node_edge)[idx_n].begin());
    }

    IncidentIterator edge_end () const{
      return IncidentIterator(graph_n, idx_n, (graph_n->node_edge)[idx_n].end());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class IncidentIterator;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_n;
    size_type idx_n;
    Node(const Graph* graph, size_type idx)
      :graph_n(const_cast<Graph*>(graph)), idx_n(idx) {
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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    nodes.emplace(size(),position);
    node_value.emplace(size(),0);
    return Node(this,size()-1);        
  }

  Node add_node(const Point& position, const node_value_type& value){
    nodes.emplace(size(),position);
    node_value.emplace(size(),value);
    return Node(this,size()-1);   
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.graph_n==this && n.index()<=size()-1);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<size());
    return Node(this,i);        
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
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_e, idx_n1);    
      //return Node();  
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_e, idx_n2);        
      //return Node();
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (graph_e==e.graph_e && ((node1()==e.node1() && node2()==e.node2()) || 
        (node1()==e.node2() && node2()==e.node1())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (graph_e==e.graph_e && idx_e<e.idx_e);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_e;
    size_type idx_e;
    size_type idx_n1; //index of first added node
    size_type idx_n2; //index of second added node
    Edge(const Graph* graph, size_type idx1, size_type idx2, size_type idx3)
      :graph_e(const_cast<Graph*>(graph)), idx_e(idx1), idx_n1(idx2), idx_n2(idx3) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, i, edges.at(i).first, edges.at(i).second);        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    auto n1=node_edge.find(a.idx_n);
    if(n1!=node_edge.end()){
      auto n2=n1->second.find(b.idx_n);
      if(n2!=n1->second.end()){
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
    //(void) a, (void) b;   // Quiet compiler warning
    //return Edge();        // Invalid Edge
    if (has_edge(a, b)) 
      return Edge(this, node_edge[a.idx_n][b.idx_n], a.idx_n, b.idx_n);
    else {
      size_type newidx_e=edges.size();
      std::pair<size_type,size_type> newedge(a.idx_n, b.idx_n);
      edges[newidx_e]=newedge;
      node_edge[a.idx_n][b.idx_n]=newidx_e;
      node_edge[b.idx_n][a.idx_n]=newidx_e;
      return Edge(this, newidx_e, a.idx_n, b.idx_n);
    }
    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes.clear();
    node_value.clear();
    edges.clear();
    node_edge.clear();
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
    NodeIterator(const Graph* graph, size_type ind)
      :graph_ni(const_cast<Graph*>(graph)), ind_n(ind) {
    }

    /**Return the Node of current iterator
    */
    Node operator*() const{
      return Node(graph_ni, ind_n);
    }
    
    /**Increment the iterator to next node
    */
    NodeIterator& operator++(){
      ind_n++;
      return *this;
    }
    
    /**Test whether the two iterators are equal
    */
    bool operator==(const NodeIterator& other) const{
      return (graph_ni==other.graph_ni && ind_n==other.ind_n);
    }

    /**Test whether two iterators are not equal
    */
    bool operator!=(const NodeIterator& other) const{
      return (graph_ni!=other.graph_ni || ind_n!=other.ind_n);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    friend class Node;
    Graph* graph_ni;
    size_type ind_n;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /**Return a node iterator pointing to the beginning of all nodes
  */
  NodeIterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /**Return a node iterator pointing to the end of all nodes
  */
  NodeIterator node_end() const{
    return NodeIterator(this, nodes.size());
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
    Graph* graph_ii;
    size_type ind_n;
    std::unordered_map<size_type,size_type>::const_iterator iter;

    IncidentIterator(const Graph* graph, size_type ind, std::unordered_map<size_type,size_type>::const_iterator it)
      :graph_ii(const_cast<Graph*>(graph)), ind_n(ind), iter(it) {
    }
    
    /**Return the edge of current iterator
    */
    Edge operator*() const{
      size_type n1=(graph_ii->edges)[iter->second].first;
      size_type n2=(graph_ii->edges)[iter->second].second;
      size_type n3;
      if(ind_n==n1) n3=n2;
      else n3=n1;
      return Edge(graph_ii, iter->second, ind_n, n3);
    }
    
    /**Increment the iterator to next edge
    */
    IncidentIterator& operator++(){
      iter=++iter;
      return *this;
    }

    /**Test whether two iterators are equal
    */
    bool operator==(const IncidentIterator& other) const{
      return (graph_ii==other.graph_ii && ind_n==other.ind_n && iter==other.iter);
    }

    /**Test whether two iterators are not equal
    */
    bool operator!=(const IncidentIterator& other) const{
      return (graph_ii!=other.graph_ii || ind_n!=other.ind_n || iter!=other.iter);
    }

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
    EdgeIterator(const Graph* graph, size_type ind)
      :graph_ei(const_cast<Graph*>(graph)), ind_e(ind) {
    }
    
    /**Return the edge of current iterator
    */
    Edge operator*() const{
      return Edge(graph_ei, ind_e, (graph_ei->edges)[ind_e].first, (graph_ei->edges)[ind_e].second);
    }
    
    /**Increment the iterator to next edge
    */
    EdgeIterator& operator++(){
      ind_e++;
      return *this;
    }

    /**Test whether two iterators are equal
    */
    bool operator==(const EdgeIterator& other) const{
      return(graph_ei==other.graph_ei && ind_e==other.ind_e);
    }

    /**Test whether two iterators are not equal
    */
    bool operator!=(const EdgeIterator& other) const{
      return(graph_ei!=other.graph_ei || ind_e!=other.ind_e);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    friend class Edge;
    Graph* graph_ei;
    size_type ind_e;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /**Return an edge iterator pointing to the beginning of all edgess
  */
  EdgeIterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /**Return an edge iterator pointing to the beginning of all edges
  */
  EdgeIterator edge_end() const{
    return EdgeIterator(this, edges.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
