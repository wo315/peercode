#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
//--functionality_0
//--missing header <map>
//--START
#include <map>
//-END
#include <tuple>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 * @tparam V represents a value for each node
 * @tparam E represents a value for each edge
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = float>
class Graph {
 private:
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
  /** Synonym for template value type */
  using node_value_type = V;
  /** Synonym for template value type */
  using edge_value_type = E;
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
  using uid_type = int;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    node_uid_counter_ = 0;
    edge_uid_counter_ = 0;
    //Assign each graph a unique id for comparison of objects in different graphs
    g_id_ = CME212::random();

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
    Node() {
      uid_ = -1;
    }

    /** Return this node's position. */
    //Should throw an error if invalid node
    const Point& position() const {
      return (graph_->node_vec_).at(uid_).point_;
    }

    /** Return this node's position */
    //Returns a reference so position can be changed
    Point& position(){
      return (graph_->node_vec_).at(uid_).point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    //Should allow for invalid node to be formed with sentinel value -1
    size_type index() const {
      return (uid_ >= 0) ? (graph_->node_vec_).at(uid_).index_ : -1;
    }

    /** Return this node's value. */
    //Should throw an error if invalid node
    node_value_type& value(){
      return (graph_->node_vec_).at(uid_).node_value_;
    }
    const node_value_type& value() const{
      return (graph_->node_vec_).at(uid_).node_value_;
    }

    /** Return the degree of this node. */
    //Should throw an error if invalid node
    size_type degree() const{
      //Check if valid index
      return graph_->adj_vec_.at(this->index()).size();
    }

    /** IncidentIterator begin(). */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, this->index(), graph_->adj_vec_.at(this->index()).begin());
    }

    /** IncidentIterator end(). */
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, this->index(), graph_->adj_vec_.at(this->index()).end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (this->graph_ == n.graph_) & (this->index() == n.index());
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     * Ordering: Graph_id, then index
     */
    bool operator<(const Node& n) const {
      if (this->graph_ == n.graph_){
        return this->index() < n.index();
        }
      return (this->graph_->g_id_ < n.graph_->g_id_);
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* graph_;
    uid_type uid_;
    // Private Constructor
    Node(const Graph* graph, uid_type uid) : graph_(const_cast<Graph*>(graph)), uid_(uid) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return index_to_uid_.size();
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

  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()){
    node_vec_.push_back(nodeinfo(position, node_value, index_to_uid_.size()));
    // Add new incident_vec_ entry for the new node
    adj_vec_.push_back(std::map<uid_type,uid_type>{});
    index_to_uid_.push_back(node_uid_counter_);
    return Node(this, node_uid_counter_++);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    //Check that it's a member of this graph and graph hasn't been cleared
    if ((this != n.graph_) or (this->node_vec_.size() <= (unsigned)n.uid_)){
      return false;
    }
    //Check if this node's index is valid (may be invalid if n was deleted)
    if (this->index_to_uid_.size() <= n.index()){
      return false;
    }
    //Check if this node's uid matches the expected node id for this node's index
    return (this->index_to_uid_[n.index()] == n.uid_);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, index_to_uid_[i]);
  }

/** Remove a node
 * @param Node n, the node to be removed
 * @return 0 if the node is not removed, 1 else
 * @post: If n is not in this graph, no changes occur
 *        If n is in this graph:
 *           -The data for n is unreachable via standard graph methods for other nodes
 *             or for future node generation,
 *             although the information for n remains in node_vec_ and may be reached via n
 *           -n itself retains it's information, but it's index is no longer valid
 *             and operations on this node are invalid and unretained within the graph
 *           -num_nodes() -= 1
 *           -The last node in the graph (n_l|n_l.index()>n_j.index() for any j != l)
 *             takes n's index
 *           -All other nodes are left unchanged and valid
 *           -Any edges that are incident to n are removed, subject to the
 *             conditions of remove_edge
 *           -num_edges() -= n.degree()
 *           -Any node_iterator that points to a node with index < size(valid_nodes)-1
 *             remainds valid, with the caveat that the iterator pointing to n now points to n_l
 *           -A node_iterator pointing to the end (node_end()) is invalidated
 *           -Any incident iterators for n or n_l are invalidated  (i.e. iterators over nodes adj. to n or n_l)
 *           -Any incident iterator for node n_j pointing to n is invalidated (i.e. iterators over nodes adj. to n_j pointing to edge (n_j, n))
 *           -Any other incident iterator remains valid, since incidents are iteratred as <int, map_it>
 *              (see https://www.cplusplus.com/reference/map/map/erase/)
 * Complexity: O(1) to remove the node itself + O(d*log(d)) where d is n.degree() for incident edges
 *              Assuming d << n, O(d*log(d)) < O(n)
 */
  size_type remove_node(const Node& n){
    if (! has_node(n)){
      return 0;
    }
    size_type node_index = n.index();
  //Find all incident edges
    auto incident_edges_ = adj_vec_.at(node_index);
    for(auto it_ = incident_edges_.begin(); it_!=incident_edges_.end(); ++it_){
      //delete edges from n to the node of this node. it_-><node_uid_, edge_uid_>
      remove_edge(n, Node(this, it_->first));
    }
    //Swap adjacency indices to match the node update
    std::swap(adj_vec_[node_index], adj_vec_[adj_vec_.size()-1]);
    std::swap(index_to_uid_[node_index], index_to_uid_[index_to_uid_.size()-1]);
    //Remove node and it's adj_vec_ vector
    adj_vec_.pop_back();
    index_to_uid_.pop_back();
    //Update newly moved node's index
    node_vec_[index_to_uid_[node_index]].index_ = node_index;
    return 1;
  }
/*Wrapper for remove_node that accepts a node_iterator
 * @param n_it a node_iterator to the node to be deleted
 * @pre: n_it is a valid iterator that can be derefenced
 * @return: a valid node_iterator
 *          if a node is deleted: return an iterator to the same location
 *            as the deleted node, which now stores node n_l if n isn't n_l
 *          if n is n_l, returns an iterator == node_end()
 *            This is valid since iterators just stored indices as ints
 *          if a node is not deleted:
 *            Return an iterator to one past the end of the nodes
 * @post: see remove_node(Node)
 */
  node_iterator remove_node(node_iterator n_it){
    size_type node_it_index = n_it.it_;
    size_type s = remove_node(*n_it);
    if(s){
      //Just an int
      return NodeIterator(this, node_it_index);
    }
    return this->node_end();
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
  class Edge : private totally_ordered<Edge>  {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return the length of an edge*/
    double length() const {
      return norm(this->node1().position()-this->node2().position());
    }

    /** Return a node of this Edge */
    Node node1() const {
      //If edge is given as (a,b) but we want (b,a), return b
      if (flip_){
        return Node(graph_, graph_->edge_vec_[uid_].node_b_uid_);
      }
      return Node(graph_, graph_->edge_vec_[uid_].node_a_uid_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //If edge is given as (a,b) but we want (b,a), return a
      if (flip_){
        return Node(graph_, graph_->edge_vec_[uid_].node_a_uid_);
      }
      return Node(graph_, graph_->edge_vec_[uid_].node_b_uid_);
    }

    /** Return the value associated with this edge*/
    //Throws an error if this isn't a valid edge
    edge_value_type& value(){
      return graph_->edge_vec_.at(uid_).value_;
    }

    const edge_value_type& value() const {
      return graph_->edge_vec_.at(uid_).value_;
    }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     * Checks wether the edges are from the same graph, then whether they are the same edge.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ != e.graph_){
        return false;
      }
      return (this->uid_ == e.uid_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     * Compares graph_id, then edge_id
     */
    bool operator<(const Edge& e) const {
      if (this->graph_ == e.graph_){
       return (this->uid_ < e.uid_);
      }
      return (this->graph_->g_id_ < e.graph_->g_id_);
      }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph_;
    uid_type uid_;
    //Used to determine which node to return
    bool flip_;
    // Private constructor
    Edge(const Graph* graph, uid_type uid, bool flip = false) : graph_(const_cast<Graph*>(graph)), uid_(uid), flip_(flip){};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return valid_edges_.size();
    }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   * Throws an error is i is not a valid edge index
   */
  Edge edge(size_type i) const {
    return Edge(this, valid_edges_.at(i));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   * Actual complexity: O(a.degree())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Verify valid node for first call
    // If a is not a valid node of this graph, throw an error
    return ((adj_vec_.at(a.index())).find(b.uid_) != adj_vec_[a.index()].end());
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
   * Complexity: O(max(b.degree(),a.degree())
   */
  Edge add_edge(const Node& a, const Node& b) {
    //Verify that these are valid nodes of the graph
    if ((! has_node(a)) || (!has_node(b))){
      //Return an invalid edge
      return Edge();
    }
    //Check to see if this edge already exists
    auto edge_uid_it_ = (adj_vec_.at(a.index())).find(b.uid_);
    bool flip = false;
    //If edge already exists
    if (edge_uid_it_ != adj_vec_[a.index()].end()){
      //If b isn't the second node of this edge, flip the nodes
      //Edge_uid_it_-> <node_uid_, edge_uid_>
      //node_uid_ is guaranteed to match, but the edge explicitly stores
      // one node as "a" and one as "b"
      if (b.uid_ != edge_vec_[edge_uid_it_->second].node_b_uid_){
        flip = true;
      }
      return Edge(this, edge_uid_it_->second, flip);
    }
    // Indices are guaranteed by has_edge.
    edge_vec_.push_back(edgeinfo(valid_edges_.size(), edge_uid_counter_, a.uid_, b.uid_));
    valid_edges_.push_back(edge_uid_counter_);

    adj_vec_[a.index()].insert(std::pair<uid_type, uid_type>(b.uid_, edge_uid_counter_));
    adj_vec_[b.index()].insert(std::pair<uid_type, uid_type>(a.uid_, edge_uid_counter_));
    return Edge(this, edge_uid_counter_++);
  }



  /** Remove an edge
   *  @param a the first node of the edge to be removed
   *  @param b the second node of the edge to be removed
   *  @return: 0 if the edge is not removed, else 1
   *  @pre: a and b are valid nodes of this graph
   *        In the case of remove_node, a and b were valid
   *        before remove_node was called
   *  @post: if the edge is not a valid edge(a and b are not connected)
   *          nothing is changed
   *          if a and b define a valid edge e:
   *            -The information for edge e is unreachable via standard graph methods
   *              but may remain accessible via any existing objects that represent e
   *              and stille exists in edge_vec_
   *            -num_edges() -= 1
   *            -The last edge in the graph (e_l|e_l.index()>e_j.index() for any j != l)
   *             takes e's index
   *            -All other edges are left unchanged and valid and retain their index
   *            -Any incident_iterators not pointing to e remain valid
   *              (see https://www.cplusplus.com/reference/map/map/erase/)
   *            -Any edge_iterator that points to a edge with index < size(valid_edges)-1
   *             remainds valid, with the caveat that the iterator pointing to e now points to e_l
   *            -A node_iterator pointing to the end (edge_end()) is invalidated
   * Complexity: O(log(max(a.degree(),b.degree()))
   */
  size_type remove_edge(const Node& a, const Node& b){
    //Returns 0 if nothing is erased
    //O(log(a.degree()))
    size_type s = adj_vec_.at(a.index()).erase(b.uid_);
    if(s>0){
      //get edge_uid from other adjacency
      //O(log(b.degree()))
      auto del_edge_uid_ = (adj_vec_.at(b.index())).find(a.uid_)->second;
      //delete second implementation
      adj_vec_.at(b.index()).erase(a.uid_);
      //Get edge's index
      size_type del_edge_index_ = edge_vec_[del_edge_uid_].index_;
      //swap last element with this element
      std::swap(valid_edges_[del_edge_index_], valid_edges_[valid_edges_.size()-1]);
      //Update the index for our new edge
      edge_vec_[valid_edges_[del_edge_index_]].index_ = del_edge_index_;
      //Delete the last element
      valid_edges_.pop_back();
    }
    return s;
  }

  /** Remove an edge by the edge directly
    * Wraps remove_edge(node, node)
    * @pre: Edge is a member of this graph with two valid nodes a and b
    *       a and b need not be connected, but must be actual nodes
    */
  size_type remove_edge(const Edge& e){
    Node a = e.node1();
    Node b = e.node2();
    return remove_edge(a, b);
  }
  /** Remove an edge by it's iterator
    * Wraps remove_edge(node, node)
    * @return: if successful:
    *            an EdgeIterator pointing to the same location (since EdgeIterator uses ints)
    *             if e was thas last edge, now points to edge_end()
    *             otherwise points to a e_l
    *           else:
    *             an EdgeIterator pointing to the end of the valid_edges
    */
  edge_iterator remove_edge(edge_iterator e_it){
    //Store the index that we are going to swap to
    size_type edge_index_ = e_it.edge_index_;
    size_type s = remove_edge(*e_it);
    if (s){
      //Return to where we were, which now points to a new edge with the same index
      //If this was the last element, we now point one past the end
      return EdgeIterator(this, edge_index_);
    }
    return this->edge_end();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0s
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  //--functionality_1
  //--should also reset edge/node counters
  //--START
  void clear() {
    node_vec_.clear();
    adj_vec_.clear();
    edge_vec_.clear();
    valid_edges_.clear();
    index_to_uid_.clear();
  }
  //--END

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>  {
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

    /** Dereference NodeIterator. */
    Node operator*() const{
      return Node(graph_, graph_->index_to_uid_[it_]);
    }

    /** Increment NodeIterator. */
    NodeIterator& operator++(){
      ++it_;
      return *this;
    }

    /** Check NodeIterator equivalence.
     * Checks that it points to the same graph and the same iterator.
     */
    bool operator==(const NodeIterator& node_it_b) const{
      return ((graph_ == node_it_b.graph_) & (it_ == node_it_b.it_));
    }

   private:
    friend class Graph;
    Graph* graph_;
    /**  We only call nodes by index (not by a reference to an object), so iterate over indices.
     *   The index refers to both the point (point_vec_) and the value (node_values_).
     *   It can also be used to call incident edges to a given node.
     */
    size_type it_;
    // Private constructor
    NodeIterator (const Graph* graph, size_type it = 0) : graph_(const_cast<Graph*>(graph)), it_(it) {}

  };

  /** Return a NodeIterator pointing to the first node. */
  node_iterator node_begin() const{
    return(NodeIterator(this, 0));
  }

  /** Return a NodeIterator pointing past the last node. */
  node_iterator node_end() const{
    return(NodeIterator(this, index_to_uid_.size()));
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

    /** Derefernce IncidentIterator. */
    Edge operator*() const{
      bool flip_ = (graph_->edge_vec_[it_->second].node_a_uid_ == it_->first);
      return Edge(graph_, it_->second, flip_);
    }

    /** Increment IncidentIterator. */
    IncidentIterator& operator++(){
      it_++;
      return *this;
    }

    /** Check for IncidentIterator Equality.
     *  Checks that it points to the same graph and the same indices.
     */
    bool operator==(const IncidentIterator& incident_2) const{
      return ((graph_ == incident_2.graph_) & (node_1_index_ == incident_2.node_1_index_)&(it_ == incident_2.it_));
    }
   private:
    friend class Graph;
    Graph* graph_;
    // The first node is constant and is also the index of node_1.
    const size_type node_1_index_;
    /** Here, we iterate over the indices of <b>a vector of nodes incident to node_1 </b>
     * We must do it this way because we only store incident nodes and therefore the actual index
     * is not ordered.
     * Ex: incident_vec_[node_1_index_][node_2_index_] = index of (node_2_index_)th node incident to node_1_.
     */
    std::map<uid_type, uid_type>::iterator it_;
    // Private constructor
    IncidentIterator(const Graph* graph, const size_type node_1_index, std::map<uid_type, uid_type>::iterator it) : graph_(const_cast<Graph*>(graph)), node_1_index_(node_1_index),it_(it) {}
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
    /** Dereference EdgeIterator. */
    Edge operator*() const{
      size_type tmp_uid = graph_->valid_edges_[edge_index_];
      return Edge(graph_, graph_->edge_vec_[tmp_uid].edge_uid_);
    }

    /** Increment EdgeIterator. */
    EdgeIterator& operator++(){
      edge_index_++;
      return *this;
    }

    /** Check for EdgeIterator Equality.
     *  Checks that it points to the same graph and the same edges.
     */
    bool operator==(const EdgeIterator& edge_iterator_2) const{
      return ((graph_ == edge_iterator_2.graph_) & (edge_index_ == edge_iterator_2.edge_index_));
    }

   private:
    friend class Graph;
    Graph* graph_;
    /**  We can call edges by index (rather than by its nodes), so iterate over indices.
     *   If we iterated over incident edges, we would have to check to not double count.
     *   By calling by index, we guarantee to count each edge once and only once.
     */
    uint edge_index_;
    // Private onstructor
    EdgeIterator(const Graph* graph, uint edge_index = 0) : graph_(const_cast<Graph*>(graph)), edge_index_(edge_index) {}
  };

  /** Return EdgeIterator pointing to first edge */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Return EdgeIterator pointing past last edge */
  edge_iterator edge_end() const {
    return EdgeIterator(this, valid_edges_.size());
  }

 private:
  /**Struct to contain nodes info
   * Holds a point, value, and index
   */
  struct nodeinfo{
    Point point_;
    node_value_type node_value_;
    size_type index_;
    nodeinfo(Point point, node_value_type node_value, size_type index) : point_(point), node_value_(node_value), index_(index) {};
  };

  /**Struct to contain edges info
   * Holds an index, uid_, ordered uid_ for each node,
   * and edge_value
   */
  struct edgeinfo{
    size_type index_;
    uid_type edge_uid_;
    uid_type node_a_uid_;
    uid_type node_b_uid_;
    edge_value_type value_;
    edgeinfo(size_type index, uid_type edge_uid, uid_type node_a_uid, uid_type node_b_uid, edge_value_type value = edge_value_type{}):
      index_(index), edge_uid_(edge_uid), node_a_uid_(node_a_uid), node_b_uid_(node_b_uid), value_(value) {}
  };

  // Vector containing node information - never deleted (except via clear), always holds all (current and past) nodes info
  std::vector<nodeinfo> node_vec_;
  // Assign each node a uid_
  uid_type node_uid_counter_;
  // Assign each edge a uid_
  uid_type edge_uid_counter_;
  //Mapping of node's index -> uid_
  std::vector<uid_type> index_to_uid_;


  // Vector containing edge information - never deleted (except via clear), always holds all (current and past) edges info
  std::vector<edgeinfo> edge_vec_;
  // Vector holding the indices of valid edges; edge variant of index_to_uid_
  std::vector<uid_type> valid_edges_;

  // Incident map stored within a vector
  // Vector contains one map per valid node, index by node.index()
  // Map contains a pair of <node_uid_, edge_uid_> for any edge given by
  // Edge(a.index(), b.uid_) where a indexes the vector and b indexes the map
  std::vector<std::map<uid_type, uid_type>> adj_vec_;
  //g_id to compare nodes, edges between disjoint graphs
  double g_id_;
};

#endif // CME212_GRAPH_HPP
