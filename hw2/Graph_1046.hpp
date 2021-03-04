#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <climits>
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
template <typename V, typename E>
class Graph {

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  typedef Graph graph_type;
  typedef V node_value_type;
  typedef E edge_value_type;

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

  struct node_data{
    Point p;
    node_value_type val;
    size_type n_idx;
  };

  struct edge_data{
    size_type n1_id;
    size_type n2_id;
    edge_value_type val;
    size_type e_idx;
  };

  //node id to node data
  std::unordered_map<size_type, node_data > nodes;
  //edge id to edge data
  std::unordered_map<size_type, edge_data > edges;
  //node id to (node id to edge id)
  std::unordered_map<size_type, std::unordered_map<size_type, size_type> > node_edge;
  //node index to node id
  std::vector<size_type> idx2id_n;
  //edge index to edge id
  std::vector<size_type> idx2id_e;


 public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph(){}

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

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class IncidentIterator;
    Graph* graph_n;
    size_type n_id;
    Node(const Graph* graph, size_type idx)
      :graph_n(const_cast<Graph*>(graph)), n_id(idx) {
    }

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
    Node() {}

    /*check if this node is valid and not removed from graph
    */
    bool valid() const {
      return (n_id>=0 && n_id<graph_n->nodes.size()
        && graph_n->nodes[n_id].n_idx<graph_n->idx2id_n.size()
        && graph_n->idx2id_n.at(graph_n->nodes[n_id].n_idx)==n_id);
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_n->nodes[n_id].p;
    }

    Point& position (){
      return graph_n->nodes[n_id].p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_n->nodes[n_id].n_idx;
    }

    /*return this node's id*/
    size_type id() const {
      return n_id;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_n==n.graph_n && index()==n.index() && value()==n.value());
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
      return (graph_n==n.graph_n && index()<n.index()) || graph_n!=n.graph_n;
    }

    /*return value of the node
    */
    node_value_type& value(){
      return graph_n->nodes[n_id].val;
    }

    /*return constant value of the node
    */
    const node_value_type& value() const{
      return graph_n->nodes[n_id].val;
    }

    /*return the number of edges connected to the node
    */
    size_type degree() const{
      IncidentIterator it;
      size_type degree=0;
      for(it=edge_begin(); it!=edge_end(); ++it){
        degree+=1;
      }
      return degree;
    }

    /*return an incident iterator pointing to beginning of edges
    */
    IncidentIterator edge_begin () const{
      return IncidentIterator(graph_n, n_id, (graph_n->node_edge)[n_id].begin());
    }

    /*return an incident iterator pointing to end of edges
    */
    IncidentIterator edge_end () const{
      return IncidentIterator(graph_n, n_id, (graph_n->node_edge)[n_id].end());
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return idx2id_n.size();
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
    size_type new_n_idx=idx2id_n.size();
    node_data newnode={position,node_value_type(), new_n_idx};
    size_type new_n_id=nodes.size();
    nodes[new_n_id]=newnode;
    idx2id_n.push_back(new_n_id);
    return Node(this,new_n_id);
  }

  Node add_node(const Point& position, const node_value_type& value){
    size_type new_n_idx=idx2id_n.size();
    node_data newnode={position,value, new_n_idx};
    size_type new_n_id=nodes.size();
    nodes[new_n_id]=newnode;
    idx2id_n.push_back(new_n_id);
    return Node(this,new_n_id);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_n==this && n.valid());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,idx2id_n[i]);
  }

  /*return node's id if it's pinned*/
  size_type pin_id(const Point& pin) const{
    for(auto it = node_begin(); it != node_end(); ++it) {
      auto n = *it;
      if(n.position()==pin) return n.id();
    }
    return UINT_MAX-1;
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

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_e;
    size_type e_id;
    size_type n_id_1; //index of first added node
    size_type n_id_2; //index of second added node
    Edge(const Graph* graph, size_type idx1, size_type idx2, size_type idx3)
      :graph_e(const_cast<Graph*>(graph)), e_id(idx1), n_id_1(idx2), n_id_2(idx3) {
    }

   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /*check if this edge is valid and not removed from graph
    */
    bool valid() const {
      return (e_id>=0 && e_id<graph_e->edges.size()
        && graph_e->edges[e_id].e_idx<graph_e->idx2id_e.size()
        && graph_e->idx2id_e.at(graph_e->edges[e_id].e_idx)==e_id);
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_e, n_id_1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_e, n_id_2);
    }

    /** Return this edge's index, a number in the range [0, num_edges()). */
    size_type index() const {
      return graph_e->edges[e_id].e_idx;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (graph_e==e.graph_e && index()==e.index() && ((node1()==e.node1() && node2()==e.node2()) ||
        (node1()==e.node2() && node2()==e.node1())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (graph_e==e.graph_e && index()<e.index()) || graph_e!=e.graph_e;
    }

    /*return value of the edge
    */
    edge_value_type& value(){
      return graph_e->edges[e_id].val;
    }

    /*return constant value of the edge
    */
    const edge_value_type& value() const{
      return graph_e->edges[e_id].val;
    }

    /*return length of the edge
    */
    double length() const{
      return norm(node1().position()-node2().position());
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return idx2id_e.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, idx2id_e[i], edges.at(idx2id_e[i]).n1_id, edges.at(idx2id_e[i]).n2_id);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    auto n1=node_edge.find(a.n_id);
    if(n1!=node_edge.end()){
      auto n2=n1->second.find(b.n_id);
      if(n2!=n1->second.end()){
        size_type e_id=n2->second;
        return (e_id>=0 && e_id<edges.size()
          && edges.at(e_id).e_idx<idx2id_e.size()
          && idx2id_e.at(edges.at(e_id).e_idx)==e_id);
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
    if (has_edge(a, b))
      return Edge(this, node_edge[a.n_id][b.n_id], a.n_id, b.n_id);
    else {
      size_type new_e_id=edges.size();
      size_type new_e_idx=idx2id_e.size();
      edge_data newedge={a.n_id, b.n_id, edge_value_type(), new_e_idx};
      edges[new_e_id]=newedge;
      node_edge[a.n_id][b.n_id]=new_e_id;
      node_edge[b.n_id][a.n_id]=new_e_id;
      idx2id_e.push_back(new_e_id);
      return Edge(this, new_e_id, a.n_id, b.n_id);
    }
  }

  /*add edge when there is edge value
  */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value){
    if (has_edge(a, b))
      return Edge(this, node_edge[a.n_id][b.n_id], a.n_id, b.n_id);
    else {
      size_type new_e_id=edges.size();
      size_type new_e_idx=idx2id_e.size();
      edge_data newedge={a.n_id, b.n_id, value, new_e_idx};
      edges[new_e_id]=newedge;
      node_edge[a.n_id][b.n_id]=new_e_id;
      node_edge[b.n_id][a.n_id]=new_e_id;
      idx2id_e.push_back(new_e_id);
      return Edge(this, new_e_id, a.n_id, b.n_id);
    }
  }

  /** Remove an edge from the graph given two nodes, or return 0 if it does not exist.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return 1 if edge existed and successfully removed; 0 otherwise
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate the previously largest edge index
   *
   * Complexity: O(1)
   */
  size_type remove_edge(const Node& a, const Node& b){
    if(has_edge(a, b)){
      size_type e_id_rm=node_edge[a.n_id][b.n_id];
      size_type e_idx_rm=edges[e_id_rm].e_idx;
      size_type last_e_id=idx2id_e.back();
      edges[e_id_rm].e_idx=UINT_MAX;
      idx2id_e[e_idx_rm]=last_e_id;
      edges[last_e_id].e_idx=e_idx_rm;
      idx2id_e.pop_back();
      return 1;
    }
    return 0;
  }

  /** Remove an edge from the graph given the edge, or return 0 if it does not exist.
   * @return 1 if edge existed and successfully removed; 0 otherwise
   * @post !(@a e.valid())
   * @post If old (@a e.valid()), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate the previously largest edge index
   *
   * Complexity: O(1)
   */
  size_type remove_edge(const Edge& e){
    if(e.valid()){
      size_type e_id_rm=e.e_id;
      size_type e_idx_rm=edges[e_id_rm].e_idx;
      size_type last_e_id=idx2id_e.back();
      edges[e_id_rm].e_idx=UINT_MAX;
      idx2id_e[e_idx_rm]=last_e_id;
      edges[last_e_id].e_idx=e_idx_rm;
      idx2id_e.pop_back();
      return 1;
    }
    return 0;
  }

  /** Remove an edge from the graph given the edge iterator
   * @return edge iterator pointing to beginning of all edges if edge existed and
   * successfully removed; return given iterator otherwise
   * @post !(*(@a e_it).valid())
   * @post If old (*(@a e_it).valid()), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate the previously largest edge index
   *
   * Complexity: O(1)
   */
  edge_iterator remove_edge(edge_iterator e_it){
    size_type r=remove_edge(*e_it);
    if(r) return edge_begin();
    else return e_it;
  }

  /** Remove a node from the graph given the node
   * @return 1 if node existed and successfully removed; 0 otherwise
   * @post !(@a n.valid())
   * @post If old (@a n.valid()), new num_nodes() == old num_nodes() - 1.
   *       Else,                  new num_nodes() == old num_nodes().
   *
   * Can invalidate the previously largest node index
   *
   * Complexity: O(1)
   */
  size_type remove_node(const Node& n){
    if(n.valid()){
      size_type n_id_rm=n.n_id;
      size_type n_idx_rm=nodes[n_id_rm].n_idx;
      size_type last_n_id=idx2id_n.back();
      nodes[n_id_rm].n_idx=UINT_MAX;
      idx2id_n[n_idx_rm]=last_n_id;
      nodes[last_n_id].n_idx=n_idx_rm;
      idx2id_n.pop_back();
      for(auto it=n.edge_begin(); it!=n.edge_end(); ++it){
        //--style_0
        //--unused variable r
        //--START
        size_type r=remove_edge(*it);
        //--END
      }
      return 1;
    }
    return 0;
  }

  /** Remove a node from the graph given the node iterator
   * @return node iterator pointing to beginning of all nodes if node existed and
   * successfully removed; return given iterator otherwise
   * @post !(*(@a n_it).valid())
   * @post If old (*(@a n_it).valid()), new num_nodes() == old num_nodes() - 1.
   *       Else,                        new num_nodes() == old num_nodes().
   *
   * Can invalidate the previously largest node index
   *
   * Complexity: O(1)
   */
  //--functionality_1
  //--node_begin may already been removed
  //--START
  node_iterator remove_node(node_iterator n_it){
    size_type r=remove_node(*n_it);
    if(r) return node_begin();
    else return n_it;
  }
  //--END

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear();
    node_edge.clear();
    idx2id_e.clear();
    idx2id_n.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {

   private:
    friend class Graph;
    friend class Node;
    Graph* graph_ni;
    size_type n_idx;

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

    NodeIterator(const Graph* graph, size_type ind)
      :graph_ni(const_cast<Graph*>(graph)), n_idx(ind) {
    }

    /**Return the Node of current iterator
    */
    Node operator*() const{
      return Node(graph_ni, graph_ni->idx2id_n[n_idx]);
    }

    /**Increment the iterator to next node
    */
    NodeIterator& operator++(){
      n_idx++;
      return *this;
    }

    /**Test whether the two iterators are equal
    */
    bool operator==(const NodeIterator& other) const{
      return (graph_ni==other.graph_ni && n_idx==other.n_idx);
    }

    /**Test whether two iterators are not equal
    */
    bool operator!=(const NodeIterator& other) const{
      return (graph_ni!=other.graph_ni || n_idx!=other.n_idx);
    }

  };

  /**Return a node iterator pointing to the beginning of all nodes
  */
  NodeIterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /**Return a node iterator pointing to the end of all nodes
  */
  NodeIterator node_end() const{
    return NodeIterator(this, idx2id_n.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {

   private:
    friend class Graph;

   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    Graph* graph_ii;
    size_type n_id;
    std::unordered_map<size_type,size_type>::const_iterator iter;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    IncidentIterator(const Graph* graph, size_type ind, std::unordered_map<size_type,size_type>::const_iterator it)
      :graph_ii(const_cast<Graph*>(graph)), n_id(ind), iter(it) {
    }

    /**Return the edge of current iterator
    */
    Edge operator*() const{
      size_type n1=(graph_ii->edges)[iter->second].n1_id;
      size_type n2=(graph_ii->edges)[iter->second].n2_id;
      size_type n3;
      if(n_id==n1) n3=n2;
      else n3=n1;
      return Edge(graph_ii, iter->second, n_id, n3);
    }

    /**Increment the iterator to next valid edge
    */
    IncidentIterator& operator++(){
      iter=++iter;
      if(iter!=graph_ii->node_edge[n_id].end()){
        size_type e_id=iter->second;
        bool valid=(e_id>=0 && e_id<graph_ii->edges.size()
                    && graph_ii->edges[e_id].e_idx<graph_ii->idx2id_e.size()
                    && graph_ii->idx2id_e.at(graph_ii->edges[e_id].e_idx)==e_id);
        //go to next edge if not valid
        while(!valid){
          iter=++iter;
          if(iter==graph_ii->node_edge[n_id].end()) break;
          else{
            e_id=iter->second;
            valid=(e_id>=0 && e_id<graph_ii->edges.size()
                  && graph_ii->edges[e_id].e_idx<graph_ii->idx2id_e.size()
                  && graph_ii->idx2id_e.at(graph_ii->edges[e_id].e_idx)==e_id);
          }
        }
      }
      return *this;
    }

    /**Test whether two iterators are equal
    */
    bool operator==(const IncidentIterator& other) const{
      return (graph_ii==other.graph_ii && n_id==other.n_id && iter==other.iter);
    }

    /**Test whether two iterators are not equal
    */
    bool operator!=(const IncidentIterator& other) const{
      return (graph_ii!=other.graph_ii || n_id!=other.n_id || iter!=other.iter);
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {

   private:
    friend class Graph;
    friend class Edge;
    Graph* graph_ei;
    size_type e_idx;

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

    EdgeIterator(const Graph* graph, size_type ind)
      :graph_ei(const_cast<Graph*>(graph)), e_idx(ind) {
    }

    /**Return the edge of current iterator
    */
    Edge operator*() const{
      size_type e_id=(graph_ei->idx2id_e)[e_idx];
      return Edge(graph_ei,e_id,(graph_ei->edges)[e_id].n1_id, (graph_ei->edges)[e_id].n2_id);
    }

    /**Increment the iterator to next edge
    */
    EdgeIterator& operator++(){
      e_idx++;
      return *this;
    }

    /**Test whether two iterators are equal
    */
    bool operator==(const EdgeIterator& other) const{
      return(graph_ei==other.graph_ei && e_idx==other.e_idx);
    }

    /**Test whether two iterators are not equal
    */
    bool operator!=(const EdgeIterator& other) const{
      return(graph_ei!=other.graph_ei || e_idx!=other.e_idx);
    }

  };

  /**Return an edge iterator pointing to the beginning of all edgess
  */
  EdgeIterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /**Return an edge iterator pointing to the beginning of all edges
  */
  EdgeIterator edge_end() const{
    return EdgeIterator(this, idx2id_e.size());
  }

};

#endif // CME212_GRAPH_HPP
