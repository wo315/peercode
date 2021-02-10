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
 * @pre Graph must be initialized with a template type
 * Example: Graph<int> g;
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(): _num_nodes(0), _num_edges(0), 
          _points(), _values(), _e1(),  _e2(), _adj_list(), _edge_num() {}
    // HW0: YOUR CODE HERE

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
      // leaving this blank and the valid node is a private member.
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // more info about private members below.
     return _gr->_points[_node_id];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return _node_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return reference to attribute value
     */
    node_value_type& value() {

      return _gr->_values[_node_id];
    }

    /** Return const reference to attribute value
     */
    const node_value_type& value() const {

      return _gr->_values[_node_id];
    }

    /** Return degree of a node
     */
    size_type degree() const {

      return (_gr->_adj_list[_node_id]).size();
    }

    /** Return Incident Iterator corresponding to
     * beginning of edges connected to this node.
     */
    incident_iterator edge_begin() const {

      return IncidentIterator(0,_node_id,_gr);
    }

    /** Return Incident Iterator corresponding to
     * end of edges connected to this node.
     */
    incident_iterator edge_end() const {

      // get number of edges
      size_type num_inc_edges = degree();
      return IncidentIterator(num_inc_edges,_node_id,_gr);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning

      // we check if the nodes come from the same graph
      if(_gr != n._gr) {return false;}
      // we check if the nodes have the same index
      if(index() == n.index()) {return true;}
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
      (void) n;           // Quiet compiler warning

      // check if the graphs are the same first.
      if(_gr != n._gr) {
        std::cerr<<"The nodes compared don't belong to the same graph!!!";
        return false;
        }
      // checking the indices
      if(index() < n.index()) {return true;}
      return false;
    } 

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // private member tracking this node's ID wrt the graph
    size_type _node_id; 
    // _gr points back to the graph. this way
    // we can communicate with the graph & other objects
    // related to the graph
    graph_type* _gr;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    /** Construct a valid node. */
    Node(size_type _nid, const graph_type* gr) : _node_id(_nid), 
      _gr(const_cast<graph_type*>(gr)) {}
      //if gr comes from a const member function
  

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // defining an invalid behaviour
    if (_num_nodes != _points.size()){
      std::cerr<<"Unexpected behaviour! User should not see this message!!!";
      abort();
    }
    return _num_nodes;
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
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning
    (void) value;

    // add to the _points vector
    size_type new_id = _num_nodes;
    _points.push_back(position);
    _values.push_back(value);
    // since we incremenet nodes by 1, 
    // we dynamically create the adjacency list
    // NOTE: care must be taken when nodes are removed.
    std::vector<size_type> neigh {};
    _adj_list.emplace_back(neigh);
    std::vector<size_type> neigh2 {};
    _edge_num.emplace_back(neigh2);
    // update _num_nodes
    _num_nodes += 1;
    return Node(new_id, this);
  
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning

    // check if the graphs are the same.
    if(this != n._gr) {
        std::cerr<<"The node doesn't belong to the same graph!!!";
        return false;
        }
    // simple check
    if(n._node_id < _num_nodes) {return true;}
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
    (void) i;             // Quiet compiler warning

    // invalid behaviour
    if(i >= _num_nodes) {
      std::cerr<<"invalid node index given to node()!!!\n";
      abort();
    }
    return Node(i, this);

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
      // leave it blank and define a private constructor
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE

      // return node from _e1
      return _n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE

      // return node from _e1
      return _n2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE

      // similar logic as node
      if(_gr != e._gr) {return false;}
      if(_edge_id == e._edge_id) {return true;}
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE

      // similar logic as node
      if(_gr != e._gr) {
        std::cerr<<"The edges compared don't belong to the same graph!!!";
        return false;
        }
      if(_edge_id < e._edge_id) {return true;}
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    friend class Graph::NodeIterator;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // _edge_id indicates edge index in the graph
    size_type _edge_id;
    node_type _n1,_n2;
    // pointer to the graph for access to other objects' info
    graph_type* _gr;

    /** Construct a valid edge. */
    Edge(size_type _eid, node_type n1, node_type n2, const graph_type* gr) 
    : _edge_id(_eid), _n1(n1), _n2(n2),
      _gr(const_cast<graph_type*>(gr)) {}

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return _num_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning

    // invalid behaviour
    if(i >= _num_edges) {
      std::cerr<<"edge index invalid!!!";
      abort();
    }
    return Edge(i, node(_e1[i]), node(_e2[i]), this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning

    // check if nodes are valid
    if(!(has_node(a)) || !(has_node(b))){
      std::cerr<<"invalid node given to has_edge!!!\n";
      abort();
    }

    // adjacency list of node a has all it's neighbors
    std::vector<size_type> vec = _adj_list[a.index()];
    std::vector<const size_type>::iterator it;
    // search through the adj_list[node a]
    for(it = vec.begin(); it != vec.end(); it++) {
      if(*it == b.index()){
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
    (void) a, (void) b;   // Quiet compiler warning

    // important checks - do the nodes exist?
    if(!(has_node(a)) || !(has_node(b))){
      std::cerr<<"invalid node given to add_edge!!!\n";
      abort();
    }
    // avoid self-loops
    if(a==b) {
      std::cerr<<"same node given to add_edge!!!\n";
      abort();
    }

    // logic if edge already exists.
    // similar to has_edge, except use the index information
    if(has_edge(a,b)) {
      std::vector<size_type> vec = _adj_list[a.index()];
      std::vector<const size_type>::iterator it;
      int i=0;
 
      for(it = vec.begin(); it != vec.end(); it++,i++) {
        if(*it == b.index()){
          return Edge(_edge_num[a.index()][i],a,b,this);
        }
      }
    }

    // add nodes to list
    _e1.push_back(a.index());
    _e2.push_back(b.index());
    _num_edges += 1; // update #edges

    // update adjacency list.
    _adj_list[a.index()].push_back(b.index());
    _adj_list[b.index()].push_back(a.index());
    // update edge num which helps with incident iteration.
    _edge_num[a.index()].push_back(_num_edges - 1);
    _edge_num[b.index()].push_back(_num_edges - 1);

    return Edge(_num_edges - 1, a,b, this);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // reset #nodes, #edges
    _num_edges = 0;
    _num_nodes = 0;
    // reset all the vectors
    _points = std::vector<Point>();
    _e1 = std::vector<size_type>();
    _e2 = std::vector<size_type>();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : 
      public std::iterator<std::forward_iterator_tag,NodeIterator>, 
      private totally_ordered<NodeIterator>
  {
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
  
    /** Derefence node iterator
    s*/
    value_type operator*() const {

      return _gr->node(_ptr);
    }

    /** Increment node iterator
     */
    NodeIterator& operator++() {

      _ptr++;
      return *this;
    }

    /** Test whether this node iterator and @a ni are equal
    */
    bool operator==(const NodeIterator& ni) const {
      if(_gr != ni._gr) {return false;}
      if(_ptr == ni._ptr) {return true;}
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type _ptr;
    graph_type* _gr;
    NodeIterator(size_type ptr, const graph_type* gr) 
    : _ptr(ptr),
      _gr(const_cast<graph_type*>(gr)) {}

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return Node Iterator pointing to first node
   */
  node_iterator node_begin() const {

    return node_iterator(0,this);
  }

  /** Return Node Iterator pointing to last node
   */
  node_iterator node_end() const {
    
    return node_iterator(num_nodes(),this);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>,public std::iterator<std::forward_iterator_tag,IncidentIterator>   
  {
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

    /** Derefence incident iterator
    */
    value_type operator*() const {

      // use _edge_num to get this edge id
      // first index is using this node's id
      // ptr is within this node's list.
      size_type edge_idx = _gr->_edge_num[_node_id][_ptr];
      // return _gr->edge(edge_idx);
      size_type node2 = _gr->_adj_list[_node_id][_ptr];
      return Edge(edge_idx, _gr->node(_node_id), _gr->node(node2), _gr);
    }

    /** Increment incident iterator
     */
    IncidentIterator& operator++() {

      _ptr++;
      return *this;
    }

    /** Test whether this incident iterator and @a ii are equal
    */
    bool operator==(const IncidentIterator& ii) const {
      if(_gr != ii._gr) {return false;}
      if(_node_id != ii._node_id) {return false;}
      if(_ptr == ii._ptr) {return true;}
      return false;
    }


   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type _ptr, _node_id;
    graph_type* _gr;
    IncidentIterator(size_type ptr, size_type node_id, 
                     const graph_type* gr) 
    : _ptr(ptr), _node_id(node_id),
      _gr(const_cast<graph_type*>(gr)) {}

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : 
      private totally_ordered<EdgeIterator>,
      public std::iterator<std::forward_iterator_tag,EdgeIterator>
  {
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

      
    /** Derefence edge iterator
    */
    value_type operator*() const {

      return _gr->edge(_ptr);
    }

    /** Increment edge iterator
    */
    EdgeIterator& operator++() {

      _ptr++;
      return *this;
    }

    /** Test whether this edge iterator and @a ei are equal
    */
    bool operator==(const EdgeIterator& ei) const {
      if(_gr != ei._gr) {return false;}
      if(_ptr == ei._ptr) {return true;}
      return false;
    }


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    size_type _ptr;
    graph_type* _gr;
    EdgeIterator(size_type ptr, const graph_type* gr) 
    : _ptr(ptr),_gr(const_cast<graph_type*>(gr)) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return Edge Iterator pointing to first edge
   */
  edge_iterator edge_begin() const {

    return edge_iterator(0,this);
  }

  /** Return Edge Iterator pointing to last edge
   */
  edge_iterator edge_end() const {
    
    return edge_iterator(num_edges(),this);
  }


 private:

    // similar to internal_element in the examples from cme212 github repo.
    // underscore in variable name to indicate private.
    size_type _num_nodes, _num_edges; 
    // _points is a vector of Points
    // nodes are a proxy to these objects
    std::vector<Point> _points;
    std::vector<node_value_type> _values;
    // _e1 and _e2 contain node indices in the order they are added
    // _e1 contains the "smaller" node (as defined in node class).
    // _e2 contains the "larger" node.
    std::vector<size_type> _e1, _e2; // edges as a naive list 
    // adj list for faster neighbor access
    // _edge_num is like adj_list but has edge number
    // the latter helps with incident iterator 
    std::vector<std::vector<size_type>> _adj_list, _edge_num;
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
