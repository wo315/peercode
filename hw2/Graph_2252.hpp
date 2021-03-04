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
template <typename V, typename E>
class Graph {
 private:

  // Forward Declaration
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

  // Bring in template for node type
  //using node_value_type = V;
  typedef V node_value_type;
  // Bring in template for edge value type
  typedef E edge_value_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_storage(), nxt_node_uid_(1), node_cnt(0), node_idx_to_uid()
      , conn(), nodes_to_neighbors()
      , edges_storage(), nxt_edge_uid_(1), edge_cnt(0), edge_idx_to_uid() {
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
  class Node : private totally_ordered<Node> {
   public:

    /** Public Constructor : Construct an invalid node.
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


    /** Return this node's position.
     * @brief Function to obtain the underlying point object of a node
     * @pre The graph contains a valid real_node with the current uid_
     * @return Const Reference to the underlying Point object of the node
    */
    const Point& position() const {
      // HW0: travel to graph, nodes_storage, and use uid as key to retrive
      // the real node, then return it's point.
      return graph_ -> nodes_storage[uid_].p;
    }

    /** Return this node's position.
     * @brief Function to obtain the underlying point object of a node
     * @pre The graph contains a valid real_node with the current uid_
     * @return Reference to the underlying Point object of the node
    */
    Point& position() {
      // HW0: travel to graph, nodes_storage, and use uid as key to retrive
      // the real node, then return it's point.
      return graph_ -> nodes_storage[uid_].p;
    }



    /** Return this node's index
     * @brief Return current node's index in the graph.
     * @pre The graph contains a valid real_node with the current uid_
     * @return The node's index, a number in the range [0, graph_size).
    */
    size_type index() const {
      // HW0: travel to graph, nodes_storage, and use uid as key to retrive
      // the real node, then return it's index
      return graph_ -> nodes_storage[uid_].idx;
    }


    /** Return the node_value_type object stored in the graph.
     * @brief Obtain the val object stored in the node.
     * @pre The graph contains a valid real_node with the current uid_
     * @return A reference to the node_value object associated with the
     *        current node.
    */
    node_value_type& value(){
      return this->graph_->nodes_storage[uid_].val;
    }


    /** Return the node_value_type object stored in the graph.
     * @brief Obtain the val object stored in the node.
     * @pre The graph contains a valid real_node with the current uid_
     * @return A const reference to the node_value object associated with the
     *        current node.
    */
    const node_value_type& value() const {
      return this->graph_->nodes_storage[uid_].val;
    }


    /** Return the number of edges connected to graph
     * @brief Obtain the val object stored in the node.
     * @pre The graph contains a valid real_node with the current uid_
     * @return A const reference to the node_value object associated with the
     *        current node.
    */
    size_type degree() const {
      size_type e_cnt = graph_->nodes_to_neighbors[uid_].size();
      return e_cnt;
    }


    /** Incident Iterator Begin
     * @brief Crete the beginning of the incident iterator
     * @pre The graph contains a valid real_node with the current uid_
     * @return incident_iterator object with 0 as the beginning index
    */
    incident_iterator edge_begin() const {
        return incident_iterator(const_cast<Graph*>(graph_), 0, uid_);
      }


    /** Incident Iterator End
     * @brief Crete the end of the incident iterator
     * @pre The graph contains a valid real_node with the current uid_
     * @return incident_iterator object. Dereference this iterator has
     * undefined behavior.
    */
    incident_iterator edge_end() const {
      size_type max_e_idx = graph_->nodes_to_neighbors[uid_].size();
      return incident_iterator(const_cast<Graph*>(graph_), max_e_idx, uid_);
    }


    /** Node == operator
     * @brief Overloaded == operator
     * @pre both the node uid_ and the input _n_ uid_ are valid.
     * @return boolean value on whether the two nodes have the same graph and
     * the same index.
    */
    bool operator==(const Node& n) const {
      // HW0:
      return (this->graph_ == n.graph_) &&
      (this->graph_->nodes_storage[uid_].idx
        == n.graph_->nodes_storage[n.uid_].idx);
    }


    /** Node < operator
     * @brief Overload < operator
     * @pre The graph contains a valid real_node with the current uid_
     * @return true if the current node's index < input node's index
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
      size_type cidx = this->graph_->nodes_storage[uid_].idx;
      size_type nidx = n.graph_->nodes_storage[n.uid_].idx;
      return cidx < nidx;
    }

   private:
   //public:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;


    /** Private Node Constructor for the valid Node.
     * @return Graph::node object. A proxy to the actual node.
    */
    Node(const Graph* graph, size_type uid) :
      graph_(const_cast<Graph*>(graph)), uid_(uid){}


    // graph_ points to the graph that the node belong to
    Graph* graph_;
    // Invariant Unique ID of the node, used to access the map contain
    // in graph that stores all the actual nodes.
    size_type uid_;
  };


  /** Graph Size function
   * @return return the number of nodes in the graph in O(1).
  */
  size_type size() const {
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
  Node add_node(const Point& position
    , const node_value_type& node_val = node_value_type()) {
    // determine the new index, uid of the node
    size_type c_idx = size();
    size_type c_uid = nxt_node_uid_;

    // create the real_node structure
    real_node c_node;
    c_node.p = position;
    c_node.idx =  c_idx;
    c_node.val = node_val;

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
    return (this == n.graph_) && (nodes_storage.count(n.uid_) > 0);
  }


  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }


    /** Return a node of this Edge
     * @pre The graph contains a valid edge with the current uid_
     * @return Internally n1_uid will always be smaller than n2_uid. If
     * rev (decreasing flag) is true, n1_uid will be returned. Otherwise n2_uid
     * will be returned.
    */
    Node node1() const {
      size_type c_uid;
      if (rev){c_uid = (graph_->edges_storage)[uid_].n1_uid;}
      else{c_uid = (graph_->edges_storage)[uid_].n2_uid;}
      return Node(graph_,c_uid);
    }


    /** Return the other node of this Edge
     * @pre The graph contains a valid edge with the current uid_
     * @return The other node of the edge.
    */
    Node node2() const {
      size_type c_uid;
      if (rev){c_uid = (graph_->edges_storage)[uid_].n2_uid;}
      else{c_uid = (graph_->edges_storage)[uid_].n1_uid;}
      return Node(graph_,c_uid);
    }


    /** Test whether this edge and @a e are equal.
     * @pre The graph contains a valid edge with the current uid_ and @a's uid_
     * @return Whether two edges have the same graph and same nodes (undirected)
    */
    bool operator==(const Edge& e) const {
      size_type e1_n1_uid = graph_->edges_storage[uid_].n1_uid;
      size_type e1_n2_uid = graph_->edges_storage[uid_].n2_uid;
      size_type e2_n1_uid = e.graph_->edges_storage[e.uid_].n1_uid;
      size_type e2_n2_uid = e.graph_->edges_storage[e.uid_].n2_uid;

      return
      ((graph_ == e.graph_) &&
        (
          ((e1_n1_uid == e2_n1_uid ) && (e1_n2_uid == e2_n2_uid))
          ||
          ((e1_n1_uid == e2_n2_uid ) && (e1_n2_uid == e2_n1_uid))
        )
      );
    }


    /** Test whether this edge is less than @a e in a global order.
     * @pre The graph contains a valid edge with the current uid_ and @a's uid_
     * @return uid comparison of two edges.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
    */
    bool operator<(const Edge& e) const {
      //HW0: Use UID for comparison
      if (graph_ != e.graph_){
        // non-sense case, comparing edge between two different graph
        return true;
      }
      return uid_ < e.uid_;
    }


    /** Return the length of the edge */
    double length() const {
      return (graph_->edges_storage)[uid_].len;
    }

    /** Return the edge_value_type object stored in the edge.
     * @brief Obtain the val object stored in the edge.
     * @pre The graph contains a valid reall_edge with the current uid_
     * @return A reference to the edge_value object associated with the
     *        current edge.
    */
    edge_value_type& value(){
      return this->graph_->edges_storage[uid_].val;
    }


    /** Return the node_value_type object stored in the graph.
     *
    */
    const edge_value_type& value() const {
      return this->graph_->edges_storage[uid_].val;
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    /** Private Edge Constructor .
     * @param[in] graph The graph that the edge belong to.
     * @param[in] uid The unique invariant ID of the edge.
     * @param[in] rev_flag Whether node1() method should return n1_uid (rev_falg
     *  = true) or n2_uid (rev_flag = false).
     * @return Graph::edge object. A proxy to the actual edge.
    */
    Edge(const Graph* graph, size_type uid, bool rev_flag) :
      graph_(const_cast<Graph*>(graph)), uid_(uid), rev(rev_flag){

    }


    /** Private Special Edge Constructor for existing edge
     * @param[in] graph The graph that the edge belong to.
     * @param[in] uid The unique invariant ID of the edge.
     * @pre The graph contains edge with the uid _uid_.
     * @return Graph::edge object. A proxy to the actual edge. The rev flag of
     * this edge is the same as the ones stored in the graph.
    */
    Edge(const Graph* graph, size_type uid) :
      graph_(const_cast<Graph*>(graph)), uid_(uid), rev(){
        assert((*graph).edges_storage.find(uid) != (*graph).edges_storage.end());
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
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return edge_cnt;
  }


  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    size_type c_uid = edge_idx_to_uid[i];
    return Edge(this, c_uid);
  }


  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) {
    return conn[a.uid_].find(b.uid_) != conn[a.uid_].end();
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
  Edge add_edge(const Node& a, const Node& b
    , const edge_value_type& edge_val = edge_value_type()) {

    size_type c_uid;
    size_type c_n1_uid;
    size_type c_n2_uid;
    bool rev;

    if (a.uid_ < b.uid_) {
      c_n1_uid = a.uid_;
      c_n2_uid = b.uid_;
      rev = true;
    } else if (a.uid_ > b.uid_){
      c_n1_uid = b.uid_;
      c_n2_uid = a.uid_;
      rev = false;
    } else{
      assert(false);
    }

    if (conn[a.uid_].find(b.uid_) != conn[a.uid_].end()){
      // edge already exists, retrive the uid_
      c_uid = conn[a.uid_][b.uid_];
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
      c_edge.input_rev_fl = rev;
      c_edge.len = norm(a.position() - b.position());
      c_edge.val = edge_val;



      // insert into Graph
      edges_storage[c_uid] = c_edge;
      edge_idx_to_uid.push_back(c_uid);
      conn[a.uid_][b.uid_] = c_uid;
      conn[b.uid_][a.uid_] = c_uid;
      nodes_to_neighbors[a.uid_].push_back(b.uid_);
      nodes_to_neighbors[b.uid_].push_back(a.uid_);

      // updating
      ++ edge_cnt;
      ++ nxt_edge_uid_;

    }
    return Edge(this, c_uid, rev);
  }


  /** Remove Edge by nodes
   * @param[in] n1 node 1
   * @param[in] n2 node 2
   * @pre The input node 1, 2 have to be VALID node of the graph. The edge
   *       itself can exist or not exist in the graph
   * @return 1 if edge is successfully removed, 0 otherwise
   * @post Graph will not contain undirectional edge _n1_-_n2_. All associate
   *        edge objects are invalidated. Using invalidated edge object have
   *        undefined behavior / can cause error. Removing an edge will cause
   *        the size of the edges in graph decrease and will also cause
   *        a iterator pointing to the deleted edge to be pointed to a new edge
   *        so using -- to back out one edge is recommended when using
   *        edge remove inside of the edge iteerator loop.
   *        Note -- is not implemented yet, please check node iterator
   *        -- for similiar idea.
   *
  */
  size_type remove_edge(const Node& n1, const Node& n2){
    size_type edge_uid = conn[n1.uid_][n2.uid_];
    if (edge_uid == 0){
        // edge doesn't exists
        // erase the new key,val pair due to the if condition
        conn[n1.uid_].erase(n2.uid_);
        return 0;
    } else{
        // Update conn
        conn[n1.uid_].erase(n2.uid_);
        conn[n2.uid_].erase(n1.uid_);

        // Update nodes_to_neighbors
        std::vector<size_type>& nvec_n1 = nodes_to_neighbors[n1.uid_];
        std::vector<size_type>& nvec_n2 = nodes_to_neighbors[n2.uid_];


        nvec_n1.erase(
          std::remove(nvec_n1.begin(), nvec_n1.end(), n2.uid_)
          , nvec_n1.end());

        nvec_n2.erase(
          std::remove(nvec_n2.begin(), nvec_n2.end(), n1.uid_)
          , nvec_n2.end());

        // Update edge_idx_to_uid through pop and swap
        size_type e_idx = edges_storage[edge_uid].idx;

        if (edge_idx_to_uid.size() > 1){
          // pop and swap
          std::iter_swap(edge_idx_to_uid.begin() + e_idx
            , edge_idx_to_uid.rbegin());
          edge_idx_to_uid.pop_back();
          // update the index of the new edge that is in edge_idx_to_uid
          edges_storage[edge_idx_to_uid[e_idx]].idx = e_idx;

        } else {
          edge_idx_to_uid.clear();
        }

        // Update edges_storage
        edges_storage.erase(edge_uid);

        // Update edge count
        edge_cnt -= 1;
        return 1;
    }
  }


  /** Remove Edge by edge object
   * @param[in] e Lightweight edge object
   * @pre e has to point to a valid graph.
   * @return 1 if edge is successfully removed, 0 otherwise
   * @post Graph will not contain undirectional edge e. All associate
   *        edge objects are invalidated. Using invalidated edge object have
   *        undefined behavior / can cause error. Removing an edge will cause
   *        the size of the edges in graph decrease and will also cause
   *        a iterator pointing to the deleted edge to be pointed to a new edge
   *        so using -- to back out one edge is recommended when using
   *        edge remove inside of the edge iteerator loop.
   *        Note -- is not implemented yet, please check node iterator
   *        -- for similiar idea.
  */
  size_type remove_edge(const Edge& e){
    size_type edge_uid = e.uid_;

    auto it = edges_storage.find(edge_uid);

    if (it == edges_storage.end()){
      // edge doesn't exists
      return 0;

    } else{
      Node& n1 = e.node1();
      Node& n2 = e.node2();

      // Update conn
      conn[n1.uid_].erase(n2.uid_);
      conn[n2.uid_].erase(n1.uid_);

      // Update nodes_to_neighbors
      std::vector<size_type>& nvec_n1 = nodes_to_neighbors[n1.uid_];
      std::vector<size_type>& nvec_n2 = nodes_to_neighbors[n2.uid_];

      nvec_n1.erase(
        std::remove(nvec_n1.begin(), nvec_n1.end(), n2.uid_)
        , nvec_n1.end());

      nvec_n2.erase(
        std::remove(nvec_n2.begin(), nvec_n2.end(), n1.uid_)
        , nvec_n2.end());

      // Update edge_idx_to_uid through pop and swap
      size_type e_idx = edges_storage[edge_uid].idx;

      if (edge_idx_to_uid.size() > 1){
        // pop and swap
        std::iter_swap(edge_idx_to_uid.begin() + e_idx
          , edge_idx_to_uid.rbegin());
        edge_idx_to_uid.pop_back();
        // update the index of the new edge that is in edge_idx_to_uid
        edges_storage[edge_idx_to_uid[e_idx]].idx = e_idx;

      } else {
        edge_idx_to_uid.clear();
      }

      // Update edges_storage
      edges_storage.erase(edge_uid);

      // Update edge count
      edge_cnt -= 1;
      return 1;
    }
  }


  /** Remove Edge by iterator
   * @param[in] e_it edge iterator
   * @pre e_it has to be a VALID edge iterator
   * @return A valid edge iterator that either points to another random edge or
   *            the end of the edge iterator.
   * @post Graph will not contain undirectional edge e. All associate
   *        edge objects are invalidated. Using invalidated edge object have
   *        undefined behavior / can cause error. Removing an edge will cause
   *        the size of the edges in graph decrease and will also cause
   *        a iterator pointing to the deleted edge to be pointed to a new edge
   *        so using -- to back out one edge is recommended when using
   *        edge remove inside of the edge iteerator loop.
   *        Note -- is not implemented yet, please check node iterator
   *        -- for similiar idea.
   *
  */
  edge_iterator remove_edge(edge_iterator e_it){
    Edge& e = *e_it;
    size_type edge_uid = e.uid_;

    Node& n1 = e.node1();
    Node& n2 = e.node2();

    // Update conn
    conn[n1.uid_].erase(n2.uid_);
    conn[n2.uid_].erase(n1.uid_);

    // Update nodes_to_neighbors
    std::vector<size_type>& nvec_n1 = nodes_to_neighbors[n1.uid_];
    std::vector<size_type>& nvec_n2 = nodes_to_neighbors[n2.uid_];

    nvec_n1.erase(
      std::remove(nvec_n1.begin(), nvec_n1.end(), n2.uid_)
      , nvec_n1.end());

    nvec_n2.erase(
      std::remove(nvec_n2.begin(), nvec_n2.end(), n1.uid_)
      , nvec_n2.end());

    // Update edge_idx_to_uid through pop and swap
    size_type e_idx = edges_storage[edge_uid].idx;

    if (edge_idx_to_uid.size() > 1){
      // pop and swap
      std::iter_swap(edge_idx_to_uid.begin() + e_idx
        , edge_idx_to_uid.rbegin());
      edge_idx_to_uid.pop_back();
      // update the index of the new edge that is in edge_idx_to_uid
      edges_storage[edge_idx_to_uid[e_idx]].idx = e_idx;

    } else {
      edge_idx_to_uid.clear();
    }

    // Update edges_storage
    edges_storage.erase(edge_uid);

    // Update edge count
    edge_cnt -= 1;

    // return the original pointer, which now points to a different edge or
    // the end of the edge iterator
    if (edge_cnt == 0){
      return edge_end();
    }
    else{
      return e_it;
    }
  }


  /** Remove Node by node object
   * @param[in] n Lightweight edge object
   * @pre n has to point to a valid graph.
   * @return 1 if node is successfully removed, 0 otherwise
   * @post Graph will not contain node _n_. All associate
   *        node objects are invalidated. Using invalidated node object have
   *        undefined behavior / can cause error. All edges connected with node
   *        n is deleted from the graph, and their edges objects are
   *        invalidated.  Removing an node will cause
   *        shrink the size of the graph and will also cause a iterator pointing
   *        to the deleted node to be pointed to a new node so using -- to back
   *        out one node is recommended when using node remove inside of the
   *        node iteerator loop. Warning: due to deleting node will delete
   *        associated edges, using node erase during edge iterator loop will
   *        have undefined behavior.
   *
  */
  size_type remove_node(const Node& n){
    size_type node_uid = n.uid_;
    auto it = nodes_storage.find(node_uid);

    if (it == nodes_storage.end()){
      // node doesn't exists
      return 0;
    } else{
      // Edge removal phase
      std::vector<size_type> neighbors = nodes_to_neighbors[node_uid];

      for (auto it = neighbors.begin(); it != neighbors.end(); ++it){
        remove_edge(n, Node(this, *it));
      }
      // Node removal phase

      // update connection
      conn.erase(node_uid);
      // update nodes to neighbors
      nodes_to_neighbors.erase(node_uid);
      // update node_idx_to_uid through pop and swap
      size_type n_idx = nodes_storage[node_uid].idx;

      if (node_idx_to_uid.size() > 1){
        // pop and swap
        std::iter_swap(node_idx_to_uid.begin() + n_idx
          , node_idx_to_uid.rbegin());
        node_idx_to_uid.pop_back();
        // update the index of the new edge that is in edge_idx_to_uid
        nodes_storage[node_idx_to_uid[n_idx]].idx = n_idx;

      } else {
        node_idx_to_uid.clear();
      }
      // update node storage
      nodes_storage.erase(node_uid);

      // update node count
      node_cnt -= 1;

      return 1;
    }
  }



  /** Remove Node by iterator
   * @param[in] n_it node iterator
   * @pre n_it has to be a VALID edge iterator
   * @return A valid node iterator pointing to another random node, or the
   *          end of the node iterator.
   * @post Graph will not contain node _n_it_ points to. All associate
   *        node objects are invalidated. Using invalidated node object have
   *        undefined behavior / can cause error. All edges connected with node
   *        n is deleted from the graph, and their edges objects are
   *        invalidated. Removing an node will cause
   *        shrink the size of the graph and will also cause a iterator pointing
   *        to the deleted node to be pointed to a new node so using -- to back
   *        out one node is recommended when using node remove inside of the
   *        node iteerator loop. Warning: due to deleting node will delete
   *        associated edges, using node erase during edge iterator loop will
   *        have undefined behavior.
   *
  */
  node_iterator remove_node(node_iterator n_it){

    //--functionality_1
    //--cannot bind rvalue to non-const reference
    //--START
    Node& n = *n_it;
    //--END
    size_type node_uid = n.uid_;
    auto it = nodes_storage.find(node_uid);

    // Edge removal phase
    std::vector<size_type> neighbors = nodes_to_neighbors[node_uid];

    for (auto it = neighbors.begin(); it != neighbors.end(); ++it){
      remove_edge(n, Node(this, *it));
    }

    // Node removal phase

    // update connection
    conn.erase(node_uid);
    // update nodes to neighbors
    nodes_to_neighbors.erase(node_uid);
    // update node_idx_to_uid through pop and swap
    size_type n_idx = nodes_storage[node_uid].idx;
    if (node_idx_to_uid.size() > 1){
      // pop and swap
      std::iter_swap(node_idx_to_uid.begin() + n_idx
        , node_idx_to_uid.rbegin());
      node_idx_to_uid.pop_back();
      // update the index of the new edge that is in edge_idx_to_uid
      nodes_storage[node_idx_to_uid[n_idx]].idx = n_idx;

    } else {
      node_idx_to_uid.clear();
    }

    // update node storage
    nodes_storage.erase(node_uid);

    // update node count
    node_cnt -= 1;

    if (node_cnt == 0){
      return node_end();
    }
    else{
      return n_it;
    }
  }




  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // clear nodes
    nodes_storage.clear();
    node_idx_to_uid.clear();
    conn.clear();
    node_cnt = 0;
    nodes_to_neighbors.clear();

    // clear edges
    edges_storage.clear();
    edge_idx_to_uid.clear();
    edge_cnt = 0;

  }





  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
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


    /** Dereference operator overload
    * @return Node object with the current index
    */
    Node operator*() const{
      size_type n_uid = (graph_->node_idx_to_uid[cidx]);
      return Node(graph_, n_uid);
    }


    /** ++ operator overload
    */
    NodeIterator& operator++(){
      cidx ++;
      return *this;
    }

    /** ++ operator overload
    */
    NodeIterator& operator--(){
      cidx --;
      return *this;
    }

    /** == operator overload
    * @return Comparison on wither the iterator points to the same graph
    * and have the same index.
    */
    bool operator==(const NodeIterator& iter) const {
      return (graph_ == iter.graph_ && cidx == iter.cidx);
    }


   private:
    friend class Graph;

    /** Private Constructor for the Iterator
     * @param[in] input_graph The graph iterator operates on
     * @param[in] idx The index of the current iterator's node.
     *
    */
    NodeIterator(Graph* input_graph, size_type idx) :
      graph_{input_graph}, cidx{idx} {}

    // iterator contains a pointer to the graph, and the index
    Graph* graph_;
    size_type cidx;

  };


  /** Create the beginning of the node iterator
  */
  node_iterator node_begin() const {
    return NodeIterator(const_cast<Graph*>(this), 0);
  }


  /** Create the end of the node iterator
  */
  node_iterator node_end() const {
    return NodeIterator(const_cast<Graph*>(this), node_cnt);
  }




  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    /** Dereference operator overload
    * @return Edge object, edge.n1 will be the current node, and n2 will be
    * the other node.
    */
    Edge operator*() const {
      size_type n2_uid = graph_->nodes_to_neighbors[n1_uid][cidx];
      size_type e_uid = graph_->conn[n1_uid][n2_uid];
      bool rev;
      if (n1_uid < n2_uid)
        rev = true;
      else
        rev = false;
      return Edge(graph_, e_uid, rev);
    }


    /** ++ Operator Overload
    */
    IncidentIterator& operator++(){
      cidx++;
      return *this;
    }


    /** == Operator Overload
    */
    bool operator==(const IncidentIterator& iter) const{
      return ((graph_ == iter.graph_) && (cidx == iter.cidx)
        && (n1_uid == iter.n1_uid));
    }

   private:
    friend class Graph;

    /** Private Incident Iterator Construction
    * @param[in] input_graph The graph iterator operats on
    * @param[in] idx The current idx_th neighbors of n1uid.
    * @param[in] n1uid The root node that we wish to find its edges
    */
    IncidentIterator(Graph* input_graph, size_type idx, size_type n1uid) :
    graph_{input_graph}, cidx{idx}, n1_uid{n1uid}{}

    Graph* graph_;
    size_type cidx;
    size_type n1_uid;

  };




  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
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


    /** Dereference operator overload
    */
    Edge operator*() const {
      size_type e_uid = graph_->edge_idx_to_uid[cidx];
      bool rev_flag = graph_->edges_storage[e_uid].input_rev_fl;
      return Edge(graph_, e_uid, rev_flag);
    }


    /** ++ operator overload
    */
    EdgeIterator& operator++() {
      cidx++;
      return *this;
    }


    /** == operator overload
    * @return Node object with the current index
    */
    bool operator==(const EdgeIterator& iter) const {
      return (graph_ == iter.graph_) && (cidx == iter.cidx);
    }


   private:
    friend class Graph;

    /** Private Constructor
     * param[in] input_graph The graph iterator operates on
     * param[in] idx The index of the current edge iterator points to.
    */
    EdgeIterator(Graph* input_graph, size_type idx) :
      graph_{input_graph}, cidx{idx}{}

    Graph* graph_;
    size_type cidx;

  };



  // edge_iterator edge_begin() const
  edge_iterator edge_begin() const {
    return EdgeIterator(const_cast<Graph*>(this), 0);
  }

  // edge_iterator edge_end() const
  edge_iterator edge_end() const {
    return EdgeIterator(const_cast<Graph*>(this), edge_cnt);
  }



 private:

  // Node related data members
  // real_node is the data structure for the actual nodes
  struct real_node {
    Point p;
    size_type idx;
    node_value_type val;
  };

  std::unordered_map <size_type, real_node> nodes_storage;
  size_type nxt_node_uid_;
  size_type node_cnt;
  std::vector<size_type> node_idx_to_uid;
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> conn;
  std::unordered_map<size_type, std::vector<size_type>> nodes_to_neighbors;

  // Edge related data members
  struct real_edge{
    size_type n1_uid;
    size_type n2_uid;
    size_type idx;
    bool input_rev_fl;

    double len;
    edge_value_type val;
  };

  std::unordered_map<size_type, real_edge> edges_storage;
  size_type nxt_edge_uid_;
  size_type edge_cnt;
  std::vector<size_type> edge_idx_to_uid;



};

#endif // CME212_GRAPH_HPP
