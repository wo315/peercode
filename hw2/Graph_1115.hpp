#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <algorithm>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 * @pre   The value type _V_ must have a default initializer.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = int>
class Graph {
 private:

  /** Predeclaration of Internal Node type. */
  struct InternalNode_;
  /** Synonym for InternalNode_ (following STL conventions) **/
  using internal_node_type = InternalNode_;

  /** Predeclaration of Internal Edge type. **/
  struct InternalEdge_;
  /** Synonym for InternalEdge_ (following STL conventions) **/
  using internal_edge_type = InternalEdge_;

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

  /** Synonym for V. */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Synonym for E. */
  using edge_value_type = E;

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
  /** Predeclaration of type of container for InteranlNode_, i.e.,
   *  Graph::nodes_.
   *
   *  Nodes are stored in an array primarily to facilitate fast lookup. (And
   *  occasionally fast---amortized---insertion at the end of the array.
   *  Deletion is relatively rare.)
   */
  using node_container_type = std::vector<internal_node_type>;

  /** Predeclaration of type of node_container_type iterators, which are used
   *  internally to implement node iterators.
   */
  using internal_node_iterator = typename node_container_type::const_iterator;

  /** Predeclaration of type of container for user-facing node IDs, i.e.,
   *  Graph::idx_.
   *
   *  The idx_ -> uid_ mapping is stored as a vector where the uid_ of the node
   *  with idx_ == k is equal to idx_[k].
   */
  using id_container_type = std::vector<size_type>;

  /** Predeclaration of type of container for InternalEdge_, i.e.,
   *  Graph::edges_.
   *
   *  Edges are stored as a modified adjacency list, where the edge set
   *  corresponding to a node is stored as a map to balance the cost of
   *  relatively frequent insertions, deletions, and lookups.
   */
  using adjacency_list = std::map<size_type, internal_edge_type>;
  using edge_container_type = std::vector<adjacency_list>;

  /** Predeclaration of type of edge indices, i.e., a pair of size types which
   *  contain the unique identifiers for two nodes.
   */
  struct EdgeIndex_;
  using edge_index = EdgeIndex_;

  /** Predeclaration of type of container for incident edges to node, i.e.,
   *  Graph::InternalNode_::incident_edges.
   */
  using incident_container_type = std::vector<edge_index>;

  /** Predeclaration of convenience aliases of types of iterators used
   * internally in the implementation of edge_iterator and elsewhere.
   */
  using l_it_type = typename edge_container_type::const_iterator;
  using e_it_type = typename adjacency_list::const_iterator;
  using alv_type = typename adjacency_list::value_type;
  using ie_it_type = typename incident_container_type::const_iterator;


 public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
  class Node : public totally_ordered<Node> {
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

    /** Return this node's position in 3-space. */
    const Point& position() const {
      return fetch().position;
    }

    /** Return this node's position in 3-space as a modifiable reference. */
    Point& position() {
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size)
     * recording the order in which it was added to the graph. */
    size_type index() const {
      return idx_;
    }
    
    /** @brief Fetch a reference to the value of a node. */
    node_value_type& value() {
      return fetch().value;
    }

    /** @brief Fetch a constant reference to the value of a node. */
    const node_value_type& value() const {
      return fetch().value;
    }

    /** @brief Return the degree of a node. (I.e., the number of incident
     *         edges.)
     */
    size_type degree() const {
      return fetch().incident_edges.size();
    }

    /** @brief  Yields an iterator pointing to the first edge incident to a
     *          node.
     */
    incident_iterator edge_begin() const {
      InternalNode_& n = fetch();
      return IncidentIterator(graph_, n.uid, n.incident_edges.begin());
    }

    /** @brief  Synonym for edge_begin. */
    incident_iterator begin() const {
      InternalNode_& n = fetch();
      return IncidentIterator(graph_, n.uid, n.incident_edges.begin());
    }

    /** @brief  Yields an iterator pointing to the last edge incident to a
     *          node.
     */
    incident_iterator edge_end() const {
      InternalNode_& n = fetch();
      return IncidentIterator(graph_, n.uid, n.incident_edges.end());
    }

    /** @brief  Synonym for edge_end. */
    incident_iterator end() const {
      InternalNode_& n = fetch();
      return IncidentIterator(graph_, n.uid, n.incident_edges.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && idx_ == n.idx_);
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
      // Dictionary order in terms of graph location and then index.
      return graph_ != n.graph_ ? graph_ < n.graph_ : idx_ < n.idx_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // NOTE: Nodes have two identifiers: an `idx_`, which is used to index valid
    // nodes, and a `uid_`, which is only stored on the underlying
    // InternalNode_, which represents where the node's data are stored in
    // memory.
    
    size_type idx_;             // This element's unique identifier.
    graph_type * graph_;        // Pointer to parent graph.
 
    /** Check that this node points to a node that hasn't been invalidated. */
    // NOTE: Code from Lecture 9.
    bool valid() {
      size_type uid_ = fetch().uid;
                       // uid in  range.
      return uid_ >= 0 && uid_ < graph_->nodes_.size()
                       // idx  in  range.
                       && graph_->nodes_[uid_].idx < graph_->i2u_.size()
                       // uid  in  sync.
                       && graph_->i2u_[graph_->nodes_[uid_].idx] == uid_;
    }

    /** Private Constructor */
    Node(size_type idx, const graph_type* graph)
      : idx_(idx), graph_(const_cast<graph_type*>(graph)) {
    }

    /** Helper method to return the corresponding InternalNode_ */
    internal_node_type& fetch() const {
      return graph_->nodes_[graph_->i2u_[idx_]];
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return i2u_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  node_type add_node(const Point& position, const node_value_type& value = {}) {
    // The uid of the node is equal to the total number of nodes, while the idx
    // (i.e., the external id) is equal to the number of valid nodes.
    size_type uid = nodes_.size();
    size_type idx = i2u_.size();
    // Data for newest node goes to back of nodes_ vector.
    nodes_.push_back(InternalNode_(this, uid, idx, position, value)); 
    // Its idx to uid mapping goes to end of idx_ vector.
    i2u_.push_back(uid);
    // An empty edge map is initialized at the end of the edges container.
    edges_.push_back(std::map<size_type, internal_edge_type> {});
    
    // Return the new node.
    return Node(idx, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // A node is only in the graph if and only if its parent graph is this graph
    // and its index is in the range [0, num_nodes())
    return (this == n.graph_ && n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(const size_type i) const {
    // We can only return a node if the index is valid, i.e., in the range [0,
    // num_nodes()).
    assert(i < num_nodes());
    return Node(i, this);
  }

  /** @brief Method for removing nodes from the graph.
   *  @param[in]  n Node to remove.
   *  @return     1 if the node was removed, or 0 if the node was not removed
   *              (e.g., because the node is no longer valid or belongs to a
   *              different graph.)
   *  @post       If n is a valid node in the graph, pre.size() = post.size() +
   *              1. Otherwise pre.size() = post.size().
   *  @post       If n is a valid node in the graph, then post contains the same
   *              nodes as pre with the sole exception of n. Otherwise, pre and
   *              post contain the same nodes.
   *  @post       If n is a valid node in the graph, pre.num_edges() =
   *              post.num_edges() + n.degree(). Otherwise, pre.num_edges() =
   *              post.num_edges().
   *  
   *  Complexity: O(k * log(k_max)) where k is n.degree() and k_max is the
   *              maximum degree in the graph.
   *
   *  NOTE: Any edge iterators, or incident edge iterators pointing to edges in
   *  this graph that exist prior to calling this method may be invalidated
   *  after its execution. Any node iterators pointing to nodes in this graph
   *  with larger indices may be invalidated after this method's execution. Any
   *  existing Edge objects or Node objects containing nodes with larger indices
   *  may be invalidated after this method's execution.
   */
  size_type remove_node(const Node& n) {
    // Early return if node is not present
    if (! has_node(n)) return 0;

    // Create iterator to node's mapping from idx to uid
    internal_node_iterator it = nodes_.begin();
    std::advance(it, n.fetch().uid);
    remove_node(NodeIterator(this, it, nodes_.end()));

    // Return 1 to indicate success.
    return 1;
  }

  /** @brief Method for removing nodes from the graph.
   *  @param[in]  n_it  Iterator pointing to node to remove.
   *  @pre        n_it  Points to a valid node in the graph.
   *  @return     A valid pointer to a valid node in the graph, or a valid
   *              pointer equal to post.node_end().
   *  @post       *return == *(++n_it) if ++n_it is dereferenceable. Otherwise,
   *              return == g.node_end().
   *  @post       If n_it points to a valid node in the graph, pre.size() =
   *              post.size() + 1. Otherwise pre.size() = post.size().
   *  @post       If n_it points to a valid node in the graph, then post
   *              contains the same nodes as pre with the sole exception of
   *              *n_it. Otherwise, pre and post contain the same nodes.
   *  @post       If n_it points to a valid node in the graph, pre.num_edges() =
   *              post.num_edges() + (*n_it).degree(). Otherwise, pre.num_edges() =
   *              post.num_edges().
   *  @post       If n_it points to a valid node in the graph, pre contains the
   *              same edges as post except exactly the edges containing *n_it.
   *              Otherwise, pre and post contain the same edges.
   *  @post       If n is a node such that n.index() > (*n_it).index(), then
   *              pre.node(n.index()) = post.node(n.index() - 1).
   *
   *  Complexity: O(k * log(k_max)) where k is n.degree() and k_max is the
   *              maximum degree in the graph.
   *
   *  NOTE: Any edge iterators, or incident edge iterators pointing to edges in
   *  this graph that exist prior to calling this method may be invalidated
   *  after its execution. Any node iterators pointing to nodes in this graph
   *  with larger indices may be invalidated after this method's execution. Any
   *  existing Edge objects or Node objects containing nodes with larger indices
   *  may be invalidated after this method's execution.
   */
  node_iterator remove_node(node_iterator n_it) {
    // Remove each edge adjacent to the current node from the adjacency list of
    // nodes. NOTE: We cannot iterate through the adjacency list directly,
    // since each edge removal updates the underlying adjacency list.
    while ((*n_it).degree() > 0) {
      // Dereference the first edge.
      edge_type e = *(*n_it).begin();
      remove_edge(e);
    }

    // Decrement the internal idx of all subsequent nodes.
    auto subs_n_it = i2u_.begin();
    std::advance(subs_n_it, 1 + (*n_it).index());
    std::for_each(subs_n_it, i2u_.end(),
                  [this] (const size_type n) {--(this->nodes_[n].idx);});

    // Remove the current node from idx_.
    auto id_n_it = i2u_.begin();
    std::advance(id_n_it, (*n_it).index());
    i2u_.erase(id_n_it);

    // Pointer to subsequent node is still valid.
    return ++n_it;
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
  class Edge : public totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge. By default, returns the node with the
     * smaller index, unless `reverse_` is set. */
    node_type node1() const {
      return Node(reverse_ ? idx_.second : idx_.first, graph_);
    }

    /** Return the other node of this Edge. By default, returns the node with
     * the larger index, unless `reverse_` is set. */
    node_type node2() const {
      return Node(reverse_ ? idx_.first : idx_.second, graph_);
    }

    /** @brief Fetch a reference to the value of an edge. */
    edge_value_type& value() {
      return fetch().value;
    }

    /** @brief Fetch a constant reference to the value of an edge. */
    const edge_value_type& value() const {
      return fetch().value;
    }

    /** @brief Calculates the length of the edge. */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return idx_ == e.idx_ && graph_ == e.graph_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Delegate to operator<(const InternalEdge_& e).
      return graph_ == e.graph_ ? idx_ < e.idx_ : graph_ < e.graph_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // NOTE: The "orientation" of an edge---i.e., which node is returned by
    // node1() and which by node2()---is not a property of the underlying data,
    // which is symmetric in the two nodes, but rather of this wrapper. This
    // orientation is stored in `reverse_`.

    // NOTE: Edges have two identifiers: an `idx_`, which consists of the
    // corresponding idx_ values of the two nodes the edge connects, and a uid_,
    // which is used interally to store the edges in memory.

    edge_index idx_;            // This element's (external) unique identifier.
    graph_type * graph_;        // Pointer to parent graph.
    bool reverse_;              // Reverse order nodes are output if set.

    /** Private Constructor */
    Edge(edge_index idx, const graph_type* graph, const bool reverse = false)
      : idx_(idx), graph_(const_cast<graph_type*>(graph)), reverse_(reverse) {
    }

    /** Helper method to return the corresponding InternalEdge_ */
    internal_edge_type& fetch() const {
      adjacency_list l = graph_->edges_[graph_->i2u_[idx_.first]];
      return l.at(graph_->i2u_[idx_.second]);
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // Sequentially add size of each adjacency list.
    // NOTE: Since edges of invalid nodes are deleted, this method is correct,
    // even though it iterates over the adjacency lists of invalid nodes.
    auto accumulator = [] (const int& a, const adjacency_list& b) {
      return a + b.size();
    };

    return std::accumulate(edges_.begin(), edges_.end(), (size_type) 0,
                           accumulator);
  }

  /** Return the edge with (size_type) index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  edge_type edge(size_type i) const {
    // NOTE: Edges are optimized for access by edge_index, not by size_type, and
    // do not have a corresponding size_type index that they return. Therefore,
    // simply iterate through the edges until the index is reached.
    edge_iterator e_it = edge_begin();
    std::advance(e_it, i);
    return *e_it;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // If the edge exists, it must be incident to a. Therefore, loop through all
    // of the edges incident to a. If one of them joins to b, return true.
    // Otherwise, return false.
    
    // Both nodes must be in this graph.
    assert(a.graph_ == this && b.graph_ == this);

    // Since Graph is undirected, edges are indexed in a canonical way as a
    // function of the nodes they join.
    const edge_index uid(a.fetch().uid, b.fetch().uid);

    // Return if the uid matches any edge adjacent to a.
    auto same_uid = [&uid] (const edge_index i) {return uid == i;};
    return std::any_of(a.fetch().incident_edges.begin(), a.fetch().incident_edges.end(), same_uid);
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
  edge_type add_edge(const Node& a, const Node& b, const edge_value_type v = {}) {
    // `a` and `b` must be distinct nodes of this graph.
    assert(a.graph_ == this && b.graph_ == this && a != b);

    // Since Graph is undirected, edges are indexed in a canonical way as a
    // function of the nodes they join.
    const edge_index uid(a.fetch().uid, b.fetch().uid);

    // If the index of `a` is greater than the index of `b`, reverse to ensure
    // that a is return.node1() and b is return.node2()
    // NOTE: Since the mapping from idx to uid is monotonic, this computation is
    // correct.
    bool reverse = a.index() > b.index();

    // Check that edge does not already exist; if it does, return early. Only
    // have to loop through edges in one of the nodes, since both nodes have to
    // be incident for the edge to join them.
    auto same_uid = [&uid] (const edge_index i) {return uid == i;};
    bool exists = std::any_of(a.fetch().incident_edges.begin(), a.fetch().incident_edges.end(), same_uid);

    // If a match was found, return it; otherwise push a new edge onto the
    // proper location in the edge container and record that it is incident to both
    // nodes.
    if (! exists) {
      edges_[uid.first].insert({uid.second, InternalEdge_(uid, v)});
      a.fetch().incident_edges.push_back(uid);
      b.fetch().incident_edges.push_back(uid);

      // Construct the idx of the new edge.
      edge_index idx(a.index(), b.index());
      return Edge(idx, this, reverse);
    } else {
      edge_index idx(a.index(), b.index());
      return Edge(idx, this, reverse);
    }
  }

  /** @brief Method for removing edges from the graph.
   *  @param[in]  n1  One of the nodes on the edge to remove.
   *  @param[in]  n2  The other node on the edge to remove.
   *  @pre        n1 and n2 are nodes in this graph.
   *  @return     1 if the edge was removed, or 0 if the edge was not removed
   *              (e.g., because there is no edge joining n1 and n2).
   *  @post       If n1 and n2 make up a valid edge in the graph,
   *              pre.num_edges() = post.num_edges() + 1. Otherwise
   *              pre.num_edges() = post.num_edges().
   *  @post       If n1 and n2 make up a valid edge in the graph, then post
   *              contains the same edges as pre with the sole exception of
   *              edge(n1, n2). Otherwise, pre and post contain the same edges.
   *  @post       If n1 and n2 make up a valid edge in the graph, n1.degree()
   *              and n2.degree() both are reduced by 1 after execution.
   *              A node in any other circumstances has the same degree.
   *  @post       If n1 and n2 make up a valid edge, then the edges incident to
   *              n1 and n2 are the same in post as in pre, except for the edge
   *              edge(n1, n2). Otherwise, a node's incident edges are the same.
   *
   *  Complexity: O(log(n1.degre() + n2.degree()))
   *
   *  NOTE: Any existing edge iterators or incident edge iterators pointing to
   *  edges with smallest node index greater than or equal to the smaller node
   *  index of the edge pointed to by this iterator may be invalidated after
   *  this method's execution. Any existing Edge objects or Node objects
   *  contining a node of index greater than the smaller node index of the edge
   *  pointed to by this iterator may be invalidated after this method's
   *  execution.
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // Early return if edge missing.
    if (! has_edge(n1, n2)) return 0;

    // Construct the edge index.
    edge_index uid(n1.fetch().uid, n2.fetch().uid);

    // Construct an EdgeIterator pointing to this edge. This requires a pointer
    // to an adjacency list, a pointer to the correct element in the adjacency
    // list, and a stop sentinel. Edges are stored in the adjacency list
    // corresponding to the smaller index of the two nodes they join.
    l_it_type l_it = edges_.begin();
    std::advance(l_it, uid.first);
    l_it_type stop_iteration = edges_.end();
    // Get pointer to current element.
    auto e_it_finder = [&uid] (const alv_type& v) {
      return v.first == uid.second;
    };
    e_it_type e_it = std::find_if(l_it->begin(), l_it->end(), e_it_finder);

    // Delegate to iterator remove_edge method.
    remove_edge(EdgeIterator(this, l_it, stop_iteration, e_it));

    // Return 1 to indicate success.
    return 1;
  }

  /** @brief Method for removing edges from the graph.
   *  @param[in]  e  The edge to remove.
   *  @pre        e is an edge in this graph.
   *  @return     1 if the edge was removed, or 0 if the edge was not removed
   *              (e.g., because e has already been removed).
   *  @post       If e is a valid edge in the graph, pre.num_edges() =
   *              post.num_edges() + 1. Otherwise pre.num_edges() =
   *              post.num_edges().
   *  @post       If e is a valid edge in the graph, then post contains the same
   *              edges as pre with the sole exception of e. Otherwise, pre and
   *              post contain the same edges.
   *  @post       If e is a valid edge in the graph, e.node1().degree()
   *              and e.node2().degree() both are reduced by 1 after execution.
   *              A node in any other circumstances has the same degree.
   *  @post       If e is a valid edge, then the edges incident to
   *              e.node1() and e.node2() are the same in post as in pre, except
   *              for the edge e. Otherwise, a node's incident edges
   *              are the same.
   *
   *  Complexity: O(log(e.node1().degre() + e.node2().degree()))
   *
   *  NOTE: Any existing edge iterators or incident edge iterators pointing to
   *  edges with smallest node index greater than or equal to the smaller node
   *  index of the edge pointed to by this iterator may be invalidated after
   *  this method's execution. Any existing Edge objects or Node objects
   *  contining a node of index greater than the smaller node index of the edge
   *  pointed to by this iterator may be invalidated after this method's
   *  execution.
   */
  size_type remove_edge(const Edge& e) {
    // Early return if edge missing.
    if (! has_edge(e.node1(), e.node2())) return 0;

    // Use index of current edge to construct iterator pointing to it.
    edge_index uid = e.fetch().uid;

    // Construct an EdgeIterator pointing to this edge. This requires a pointer
    // to an adjacency list, a pointer to the correct element in the adjacency
    // list, and a stop sentinel. Edges are stored in the adjacency list
    // corresponding to the smaller index of the two nodes they join.
    l_it_type l_it = edges_.begin();
    std::advance(l_it, uid.first);
    l_it_type stop_iteration = edges_.end();
    // Get pointer to current element.
    auto e_it_finder = [&uid] (const alv_type& v) {
      return v.first == uid.second;
    };
    e_it_type e_it = std::find_if(l_it->begin(), l_it->end(), e_it_finder);

    // Delegate to iterator remove_edge method.
    remove_edge(EdgeIterator(this, l_it, stop_iteration, e_it));

    // Return 1 to indicate success.
    return 1;
  }

  /** @brief Method for removing edges from the graph.
   *  @param[in]  e_it  Edge iterator pointing to node to remove.
   *  @pre        e_it points to a valid node in the graph.
   *  @post       *return == *(++e_it) if ++e_it is dereferenceable. Otherwise,
   *              return == post.edge_end().
   *  @post       If e_it points to a valid edge in the graph, pre.num_edges() =
   *              post.num_edges() + 1. Otherwise pre.num_edges() =
   *              post.num_edges().
   *  @post       If e_it points to a valid edge in the graph, then post
   *              contains the same edges as pre with the sole exception of
   *              *e_it. Otherwise, pre and post contain the same edges.
   *  @post       If e_it points to a valid edge in the graph, then
   *              (*e_it).node1().degree and (*e_it).node2().degree are
   *              reduced by 1 after execution. The degrees of any other nodes
   *              remain the same. If e_it does not point to a valid edge, the
   *              degrees of every edge remain the same.
   *  @post       If e_it points to a valid edge in the graph, then
   *              the edges incident to (*e_it).node1() and (*e_it).node2() are
   *              the same in post as in pre, except for (*e_it). Otherwise, a
   *              node's incident edges are unchanged.
   *
   *  Complexity: O(log((*e_it).node1().degree() + (*e_it).node2().degree()))
   *
   *  NOTE: Any existing edge iterators or incident edge iterators pointing to
   *  edges with smallest node index greater than or equal to the smaller node
   *  index of the edge pointed to by this iterator may be invalidated after
   *  this method's execution. Any existing Edge objects or Node objects
   *  contining a node of index greater than the smaller node index of the edge
   *  pointed to by this iterator may be invalidated after this method's
   *  execution.
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    // Store the uid of the edge, along with the uids of the nodes the edge joins.
    edge_type e = *e_it;
    size_type node1_uid = e.node1().fetch().uid;
    size_type node2_uid = e.node2().fetch().uid;
    edge_index uid(node1_uid, node2_uid);

    // Copy pointer to current edge in adjacency list, since it has to be
    // advanced before it is removed.
    e_it_type e_it_del = e_it.e_it_;

    // Delete from adjacency lists of first node.
    ie_it_type it;
    it = std::find(nodes_[node1_uid].incident_edges.begin(),
                   nodes_[node1_uid].incident_edges.end(), uid);
    nodes_[node1_uid].incident_edges.erase(it);

    // Delete from adjacency lists of second node.
    it = std::find(nodes_[node2_uid].incident_edges.begin(),
                   nodes_[node2_uid].incident_edges.end(), uid);
    nodes_[node2_uid].incident_edges.erase(it);

    // The next edge in the adjacency list will still be valid after this edge
    // is deleted, so we get its index. The iterators pointing to the particular
    // adjacency list and stop sentinal remain valid.
    ++e_it;
    l_it_type l_it = e_it.l_it_;
    l_it_type stop_iteration = e_it.stop_iteration_;

    // The pointer to the edge itself may be invalidated, so we store its index
    // and find it in the modified container.
    edge_index n_uid = e_it.e_it_->second.uid;
    auto e_it_finder = [&n_uid] (const alv_type& v) {
      return v.first == n_uid.second;
    };
    e_it_type e_it_internal = std::find_if(l_it->begin(), l_it->end(),
                                           e_it_finder);

    // Delete from edge_list. (NOTE: Edges are always stored by full index at
    // the location of the first index.)
    edges_[uid.first].erase(e_it_del);
    
    // Create new iterator pointing at the next edge.
    return EdgeIterator(this, l_it, stop_iteration, e_it_internal);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    i2u_.clear();
    edges_.clear();
  }

  /** @brief  Sets a default value for all nodes.
   *  @param  v Default value for nodes.
   *  @post   n.value() == v for every node n.
   */
  void default_node_value(const node_value_type& v) {
    for (internal_node_type& n : nodes_) n.value = v;
  }

  /** @brief  Sets a default value for all edges.
   *  @param  v Default value for edges.
   *  @post   e.value() == v for every node n.
   */
  void default_edge_value(const edge_value_type& v) {
    for (adjacency_list& l : edges_)
      for (auto& kv : l)
        kv.second.value = v;
  }

  /** @brief  Sets a default value for all edges using a functor.
   *  @param  f functor with signature edge_value_type f(const edge_type&).
   */
  template <typename F>
  void default_edge_value(const F& f) {
    for (adjacency_list& l : edges_)
      for (auto& kv : l)
        kv.second.value = f(Edge(kv.second.uid, this));
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: public equality_comparable<node_iterator> {
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

    /** Dereference a NodeIterator to get the underlying Node. */
    node_type operator*() const {
      return graph_->node(it_->idx);
    }

    /** Advance the NodeIterator to the next node. */
    node_iterator& operator++() {
      // The next InternalNode_ may be invalid, so advance until a valid node or
      // the end is reached.
      do {
        ++it_;
      } while (it_ != stop_iteration_ && ! it_->valid());
      return *this;
    }

    /** Compare two iterators to see if they're equal. */
    bool operator==(const node_iterator& it) const {
      return graph_ == it.graph_ && it_ == it.it_;
    }

   private:
    friend class Graph;

    // Pointer to parent graph.
    const graph_type * graph_;
    // Iterator pointing to current node.
    internal_node_iterator it_;
    // Sentinel value pointing to end of list.
    internal_node_iterator stop_iteration_;

    /** Private constructor from iterator over Graph::node_ */
    NodeIterator(const graph_type * graph, internal_node_iterator it,
                 internal_node_iterator stop_iteration)
        : graph_(graph), it_(it), stop_iteration_(stop_iteration) {
      // The first InternalNode_ may be invalid, so advance until a valid node
      // or the end is reached.
      while (it_ != stop_iteration_ && ! it_->valid())
        ++it_;
    }
  };

  /** @brief  Yields an iterator pointing to the first node in the graph. */
  node_iterator node_begin() const {
    internal_node_iterator it = nodes_.begin();
    internal_node_iterator stop_iteration = nodes_.end();
    return NodeIterator(this, it, stop_iteration);
  }

  /** @brief  Yields an iterator pointing to the last node in the graph. */
  node_iterator node_end() const {
    internal_node_iterator it = nodes_.end();
    internal_node_iterator stop_iteration = it;
    return NodeIterator(this, it, stop_iteration);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: public equality_comparable<incident_iterator> {
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

    /** Dereference an IncidentIterator to get the underlying Node. */
    edge_type operator*() {
      edge_index uid = *it_;
      bool reverse = node_uid_ != uid.first;
      edge_index idx(graph_->nodes_[uid.first].idx, graph_->nodes_[uid.second].idx);
      return Edge(idx, graph_, reverse);
    }

    /** Advance the NodeIterator to the next node. */
    incident_iterator& operator++() {
      ++it_;
      return *this;
    }

    /** Compare two iterators to see if they're equal. */
    bool operator==(const incident_iterator& it) const {
      return (graph_ == it.graph_
              && node_uid_ == it.node_uid_
              && it_ == it.it_);
    }

   private:
    friend class Graph;

    // Pointer to parent graph.
    graph_type * const graph_;
    // Unique identifier of relevant node.
    size_type node_uid_;
    // Iterator pointing to current edge.
    typename incident_container_type::const_iterator it_;

    /** Private constructor from iterator over
     *  Graph::InternalNode_::incident_edges */
    IncidentIterator(graph_type* graph, size_type node_uid,
                     typename incident_container_type::const_iterator it)
        : graph_(graph), node_uid_(node_uid), it_(it) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: public equality_comparable<edge_iterator> {
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

    /** Dereference a EdgeIterator to get the underlying edge. */
    edge_type operator*() const {
      // Convert the node uids into idxs.
      edge_index uid = e_it_->second.uid;
      edge_index idx(graph_->nodes_[uid.first].idx,
                     graph_->nodes_[uid.second].idx);
      return Edge(idx, graph_);
    }

    /** Advance the EdgeIterator to the next edge.
     *  Note that we can only advance e_it_ if there is a valid iterator to
     *  advance it to. This can fail if (1) the current container is empty, or
     *  (2) if we've reached the end of the final container.
     */
    edge_iterator& operator++() {
      assert(l_it_ != stop_iteration_);
      ++e_it_;
      // Advance to the next dereferenceable edge or the end.
      // NOTE: Since edges containg any invalid nodes are removed, the following
      // does not need to check the validity of the node at l_it_.
      while (e_it_ == l_it_->end() && l_it_ != stop_iteration_) {
        ++l_it_;
        e_it_ = l_it_->begin();
      } 
      
      return *this;
    }

    /** Compare two iterators to see if they're equal. Since this is a nested
     * iterator, after the final outer container is exhausted, there is no
     * meaningful value for the inner iterator, and so all iterators whose outer
     * iterator have reached _stop_iteration_ should be considered equal.
     */
    bool operator==(const edge_iterator& it) const {
      if (graph_ != it.graph_ || stop_iteration_ != it.stop_iteration_) {
        return false;
      } else if (l_it_ == it.l_it_) {
        if (l_it_ == stop_iteration_) {
          // Only compare the adjacency list iterators if at the end, since
          // there is not a valid value for e_it_.
          return true;
        } else {
          return e_it_ == it.e_it_;
        }
      } else {
        return false;
      }
    }

   private:
    friend class Graph;

    // Pointer to parent graph.
    const graph_type * graph_;
    // Iterator pointing to current adjacency list.
    l_it_type l_it_;
    // Sentinel value pointing to end of list.
    l_it_type stop_iteration_;
    // Iterator pointing to current edge.
    e_it_type e_it_;

    /** Private constructor from iterator over Graph::edge_ */
    EdgeIterator(const graph_type * graph, l_it_type l_it, l_it_type
                 stop_iteration, e_it_type e_it)
      : graph_(graph), l_it_(l_it), stop_iteration_(stop_iteration),
      e_it_(e_it) {
        // Begin by advancing to a dereferenceable edge.
        while (e_it_ == l_it_->end() && l_it_ != stop_iteration_) {
          ++l_it_;
          e_it_ = l_it_->begin();
        } 
    }
  };

  /** @brief  Yields an iterator pointing to the first edge in the graph. */
  edge_iterator edge_begin() const {
    l_it_type l_it = edges_.begin();
    l_it_type stop_iteration = edges_.end();
    e_it_type e_it = l_it->begin();
    return EdgeIterator(this, l_it, stop_iteration, e_it);
  }

  /** @brief  Yields an iterator pointing to the last edge in the graph. */
  edge_iterator edge_end() const {
    l_it_type l_it = edges_.end();
    l_it_type stop_iteration = l_it;
    // NOTE: The value of e_it doesn't matter, so it is set it to an arbitrary
    // value.
    e_it_type e_it {};
    return EdgeIterator(this, l_it, stop_iteration, e_it);
  }

 private:

  //
  // Edge Index
  //

  /** @class Graph::EdgeIndex_
   *  @Internal abstraction for representing the index of an edge as a duple of
   *  size_type ints.
   *
   *  NOTE: EdgeIndex_ uniquely assigns an identifier to edges by always putting
   *  the node with the smaller index in `first`. The "orientation" of an edge
   *  is properly a property of the Edge proxy, rather than the underlying data,
   *  and is encoded by the Edge::reverse_ member.
   */
  struct EdgeIndex_ : public std::pair<size_type, size_type> {
    /** Default constructor. Creates invalid EdgeIndex_. */
    EdgeIndex_() {}

    /** Real constructor. */
    EdgeIndex_(const size_type& node1, const size_type& node2)
      : std::pair<size_type, size_type>(std::min(node1, node2),
                                        std::max(node1, node2)) {
    }
  };

  //
  // Internal Node
  //

  /** @class Graph::InternalNode_
   * @brief Class containing node internals. Proxied by Graph::Node.
   *
   * NOTE: Since InternalNode_ does not store a reference to containing graph,
   * there is no notion of equality: two nodes with identical internals should
   * not be equal if they belong to different graphs. Therefore, equality is
   * only implemented at the level of the Node wrapper.
   */
  struct InternalNode_ {
    graph_type * const graph;               // Pointer to containing graph.
                                            // (NOTE: This is necessary to check
                                            // the validity of InternalNode_s,
                                            // which cannot be delegated to
                                            // Nodes.
    const size_type uid;                    // The unique identifier for a node.
    size_type idx;                          // The external identifier for a node.
    Point position;                         // The position of the node.
    node_value_type value;                  // The value contained by the node.
    std::vector<edge_index> incident_edges; // Indices of incident edges.
                                            // NOTE: Incident edges are stored
                                            // by private uid to avoid
                                            // unecessary updates.

    /** Constructor for InternalNode_ */
    InternalNode_(const graph_type * graph, const size_type& uid,
                  const size_type& idx, const Point& position,
                  node_value_type value)
      : graph(const_cast<graph_type*>(graph)), uid(uid), idx(idx),
      position(position), value(value), incident_edges() {
    }

    /** Checks whether InternalNode_ represents a valid node. */
    bool valid() const {
      return uid >= 0 && uid < graph->nodes_.size()
                      && idx < graph->i2u_.size()
                      && graph->i2u_[idx] == uid;
    }
  };

  //
  // Internal Edge
  //

  /** @class Graph::InternalEdge_
   * @brief Class containing edge internals. Proxied by Graph::Edge */
  struct InternalEdge_ : private totally_ordered<internal_edge_type> {
    const edge_index uid;                   // The unique identifier for an edge.
    edge_value_type value;                  // The value contained by the edge.

    /** Constructor for InternalEdge_ */
    InternalEdge_(const edge_index uid, edge_value_type value)
      : uid(uid), value(value) {
    }
  };

  node_container_type nodes_;  // Collection of nodes in graph.
  id_container_type i2u_;      // Collection of user-facing node IDs.
  edge_container_type edges_;  // Collection of edges in graph.
};

#endif // CME212_GRAPH_HPP
