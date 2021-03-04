#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <stdexcept>
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
template <typename V=int, typename E=int>
class Graph {

 public:
  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  /** Type of additional attribute given to nodes (added for HW1) */
  using node_value_type = V;

  /** Type of additional attribute given to edges (added for HW2) */
  using edge_value_type = E;

 private:

  // Edge data type (to help with removal and value storage)
  struct base_edge {
    size_type uid; // the index of the edge in edges_
    edge_value_type value; // the value of the edge
  };

  // graph node with uid and position
  struct base_node {
    size_type idx;  // client-facing index of node in graph
    Point position; // physical location of node
    node_value_type value; // user-specified value of node
  };

  std::vector<base_node> nodes_; // map indices nodes/edges
  std::vector<size_type> i2u_;    // store the current active set of nodes

  // Added HW #1 to optimize edge retrieval by index
  std::vector<std::pair<size_type, size_type>> edges_;

  // Added HW #2 to store edge values
  std::unordered_map<size_type,
    std::unordered_map<size_type, base_edge>> edgemap_;

  // define typenames for adjacency list iterators
  using adj_iterator_ =
    typename std::unordered_map<size_type, base_edge>::const_iterator;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes_(), i2u_(), edges_(), edgemap_() {
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
   private:

     // validity check for node, as suggested in Lecture
     bool valid() const {
       return uid_ >= 0 && uid_ < graph_->nodes_.size()     // uid in bounds
         && graph_->nodes_[uid_].idx < graph_->i2u_.size()  // idx in bounds
         && graph_->i2u_[graph_->nodes_[uid_].idx] == uid_; // idx/uid sync
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
    Node() {
    }

    /**
     * @brief return this node's const position
     * @return const Point reference for position of node
     *
     * Complexity: O(1)
     */
    const Point& position() const {
      assert(valid());
      return const_cast<const Point&>(graph_->nodes_[uid_].position);
    }

    /**
     * @brief return this node's modifiable position
     * @return Point reference for the position of the node.
     *
     * Complexity: O(1)
     */
    Point& position() {
      assert(valid());
      return graph_->nodes_[uid_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(valid());
      return graph_->nodes_[uid_].idx;
    }

    /**
     * @brief retrieve the value held by the node
     * @return node_value_type reference for the value held by the node
     *
     * Complexity: O(1)
     */
    node_value_type& value() {
      assert(valid());
      return graph_->nodes_[uid_].value;
    }

    /**
     * @brief const version of value()
     * @return reference to const node_value_type value held by the node
     *
     * Complexity: O(1)
     */
    const node_value_type& value() const {
      assert(valid());
      return const_cast<const node_value_type&>(graph_->nodes_[uid_].value);
    }

    /**
     * @brief number of connections to other nodes
     * @return unsigned degree value of number of edges incident to the node
     *
     * @post sum of node.degree() for every node in the graph is twice the
     *  number of edges in the graph.
     *
     * Complexity: O(1)
     */
    size_type degree() const {
      assert(valid());
      // size of map from adj list
      return graph_->edgemap_[uid_].size();
    }

    /** Iterator to the first incident edge to the node
     *
     * @return incident_iterator that can can traverse all incident edges
     *  from the first edge
     * @post if the node has no incident edges, edge_begin() == edge_end()
     *
     * By convention, begins at first outbound edge to node of higher index
     * before those of lower index than the current node. If there are no
     * higher index node edges, begins at first edge to a lower index node.
     *
     * Complexity: O(1)
     */
    incident_iterator edge_begin() const {
      assert(valid());
      return IncidentIterator(*this, graph_->edgemap_[uid_].begin());
    }

    /** Iterator to the end of the set of incident edges
     *
     * @return incident_iterator that represents the end of the incident set
     * @post if the node has no incident edges, edge_begin() == edge_end()
     *
     * Complexity: O(1)
     */
    incident_iterator edge_end() const {
      assert(valid());
      return IncidentIterator(*this, graph_->edgemap_[uid_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index
     */
    bool operator==(const Node& n) const {
      assert(valid());
      assert(n.valid());

      return (uid_ == n.uid_) && (graph_ == n.graph_);
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
      assert(valid());
      assert(n.valid());

      if (graph_ == n.graph_) return uid_ < n.uid_;
      else return (graph_ < n.graph_); // arbitrary comparison to break ties
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph<V,E>;

    /** Private Node Constructor */
    Node(const Graph* graph, size_type id )
      : graph_(const_cast<Graph*>(graph))
      , uid_(id)
    { }

    Graph* graph_;  // pointer to the first node in the container
    size_type uid_; // unique id (index of node in graph_->nodes_)
  };

  /** Return the number of (valid) nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] _position_ The new node's position
   * @param[in] _value_ The new node's value (defined by graph template)
   * @return    Node that was added to the graph
   *
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized [std::vector push_back]
   */
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type()) {
    size_type uid = nodes_.size(); // size of the full node vector is uid
    size_type idx = num_nodes();   // size of the valid node vector is idx

    nodes_.push_back(base_node{idx, const_cast<Point&>(position), value});
    i2u_.push_back(uid); // this maps idx to uid

    return Node(this, uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (!n.valid()) return false;
    return this == n.graph_; // check if node points to graph
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i >= this->num_nodes()) {
      throw std::out_of_range("accessing out of bounds node");
    }
    return Node(this, i2u_[i]);
  }

  /** Remove given node from its graph
   * @param[in,out] _n_ const reference to node to be removed
   * @return 1 if remove was successful, else 0
   *
   * @pre  _n_ belongs to this graph and is valid
   * @post num_nodes() [post] == num_nodes() [pre] - 1
   * @post num_edges() [post] == num_edges() [pre] - _n_.degree()
   * @post all member methods of _n_ fail due to validity assertion
   * @post has_node(_n_) == false
   *
   * Note: Node removal does not modify the internal uid of nodes,
   *   but it invalidates idx, i.e. no guarantee that for a given
   *   node n0 != _n_, n0.index() [pre] == n0.index() [post] and
   *   likewise node(i) [pre] == node(i) [post] is not guaranteed.
   *
   *   Any Edge that contains _n_ is no longer valid to use.
   *   All Edge indices are also invalidated, i.e. edge(i) [pre] ==
   *   edge(i) [post] is not guaranteed. As a result, node_iterator,
   *   edge_iterator, and any incident_iterator for a node whose edge
   *   has been deleted as a result of the operation are invalidated.
   *   If using remove_node() in an iterator loop, call iterator.begin()
   *   to reset the iterator, or call iterator overload
   *   remove_node(node_iterator). This can be used to reset the forward
   *   iterator to a valid position in the loop.
   *
   * Complexity: O(num_edges() / num_nodes()) on average
   */
  size_type remove_node(const Node& n) {

    if (!has_node(n)) return 0; // ensure node belongs to graph

    IncidentIterator ii = n.edge_begin();
    IncidentIterator ii_end = n.edge_end();

    // remove each incident edge
    while(ii != ii_end) {
      remove_edge(*ii);
      ii = n.edge_begin(); // reset as iterator is invalidated
    }

    // remove node via swap & pop from i2u_ only
    size_type n1_idx = n.index();
    size_type nn_idx = num_nodes() - 1;

    nodes_[i2u_[nn_idx]].idx = n1_idx; // update outward idx for last node
    i2u_[n1_idx] = i2u_[nn_idx];       // move last uid to this idx
    i2u_.pop_back();                   // pop last item from i2u_

    return 1;
  }

  /** Remove node from graph given iterator to node
   * @param _n_it_ node_iterator to node we wish to remove
   * @return iterator to valid node in the graph - node_begin() by default.
   *
   * Calls remove_node(n) on node at (*_n_it_). See above
   * for complete specification. Returns the same iterator due to the
   * swap and pop method. Failing to assign the iterator to this value
   * can result in invalid iteration (i.e. skipping).
   *
   * Complexity: O(num_edges() / num_nodes()) average
   */
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return n_it;
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

   // checking for validity of edge on all methods proved to hinder
   // performance. including validity as a precondition on all member methods.
   // throw exception in value() retrieval only.
   private:
    bool valid() const {
      return node1().graph_->has_edge(node1(), node2());
    }

   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge
     * @return the first node of the edge (as constructed)
     * @pre this.valid() == true
     */
    Node node1() const {
      return node1_;
    }

    /** Return a node of this Edge
     * @return the second node of the edge (as constructed)
     * @pre this.valid() == true
     */
    Node node2() const {
      return node2_;
    }

    /** Return value held at edge
     * @return reference to edge_value_type defined by Graph template
     *
     * @pre  this.valid() == true
     * @post for e1 = Edge(n1, n2) and e2 = Edge(n2, n1),
     *       e1.value() == e2.value()
     *
     * Note: The return of reference allows assignment of edge_value_type
     *   to value(). Due to the undirected nature of the graph, one subtlety
     *   here is that the value is retrieved from and stored at
     *   edgemap_[n1_uid][n2_uid] where n1_uid < n2_uid by convention,
     *   while the value at edgemap_[n2_uid][n1_uid] is the user-defined
     *   default. The operation fails if at least one node is invalid. An
     *   exception is thrown if the edge is invalid.
     *
     * Complexity: O(1) average [std::unordered_map access]
     */
    edge_value_type& value() {
      if (!valid()) {
        throw std::runtime_error("accessing value of invalid edge");
      }
      auto nodes = std::minmax({node1(), node2()});
      size_type n1_uid = node1().graph_->i2u_[nodes.first.index()];
      size_type n2_uid = node1().graph_->i2u_[nodes.second.index()];
      return node1().graph_->edgemap_[n1_uid][n2_uid].value;
    }

    /** Return const value held at edge
     * @return reference to const edge_value_type defined by Graph template
     *
     * @pre  this.valid() == true
     * @post for e1 = Edge(n1, n2) and e2 = Edge(n2, n1),
     *       e1.value() == e2.value() [see Note above]
     *
     * Complexity: O(1) average
     */
    const edge_value_type& value() const {
      if (!valid()) {
        throw std::runtime_error("accessing value of invalid edge");
      }
      edge_value_type& e_val = value();
      return const_cast<const edge_value_type&>(value());
    }

    /** Test whether this edge and @a e are equal.
     * @param[in] _e_ edge against which to compare this edge
     * @pre this.valid() == true and _e_.valid() == true
     *
     * Equal edges represent the same undirected edge between two nodes.
     *
     * (Note that operation will fail if at least one node is invalid)
     */
    bool operator==(const Edge& e) const {
      // sort nodes to compare undirected edges
      auto nodes1 = std::minmax({this->node1(), this->node2()});
      auto nodes2 = std::minmax({e.node1(), e.node2()});
      return (nodes1.first == nodes2.first) && (nodes1.second == nodes2.second);
    }

    /** Test whether this edge is less than @a e in a global order.
     * @param[in] _e_ edge against which to compare this edge
     * @pre this.valid() == true and _e_.valid() == true
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     *
     * (Note that operation will fail if at least one node is invalid)
     */
    bool operator<(const Edge& e) const {

      // sort nodes to compare undirected edges
      auto nodes1 = std::minmax({this->node1(), this->node2()});
      auto nodes2 = std::minmax({e.node1(), e.node2()});

      // determine inequality via "dictionary" order of nodes
      if (nodes1.first == nodes2.first) {
        return nodes1.second < nodes2.second;
      }
      else {
        return nodes1.first < nodes2.first;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph<V, E>;

    Node node1_; // first node
    Node node2_; // second node

    /** Private Edge Constructor */
    Edge(Node node1, Node node2)
      : node1_(node1)
      , node2_(node2)
    { }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Note: By convention, edge(i).node1() < edge(i).node2()
   *
   * Complexity: O(1) average
   */
  Edge edge(size_type i) const {

    if (i >= num_edges()) {
      throw std::out_of_range("accessing out of bounds edge");
    }

    return Edge(node(edges_[i].first), node(edges_[i].second));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Note: If _a_ or _b_ is invalidated, returns false.
   *
   * Complexity: O(1) average [std::unordered_map find, at]
   */
  bool has_edge(const Node& a, const Node& b) const {
    // make sure nodes belong to the graph in question
    if (a == b) return false;
    if (!this->has_node(a) || !this->has_node(b)) return false;

    size_type a_uid = i2u_[a.index()];
    size_type b_uid = i2u_[b.index()];

    // first ensure that node is involved in _any_ edge
    if (edgemap_.find(a_uid) == edgemap_.end()) return false;
    else return edgemap_.at(a_uid).find(b_uid) != edgemap_.at(a_uid).end();
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
   * Complexity: O(1) average [std::unordered_map insert, std::vector push_back]
   */
  Edge add_edge(const Node& a, const Node& b) {

    // extra checks for precondition
    if (a == b) {
      throw std::runtime_error("nodes of an edge must be distinct");
    }
    if (!this->has_node(a) || !this->has_node(b)) {
      throw std::runtime_error("nodes must exist in the same graph");
    }

    Edge e(a, b);

    // if the edge already exists, then we already have a valid edge
    if (!this->has_edge(a, b)) {
      auto e_nodes = std::minmax({a, b});
      size_type e_id = this->num_edges();
      size_type n1_uid = i2u_[e_nodes.first.index()];
      size_type n2_uid = i2u_[e_nodes.second.index()];
      edges_.push_back(std::make_pair(n1_uid, n2_uid)); // add to edge list

      // add edge to both nodes' adjacency lists
      edgemap_[n1_uid].insert(std::make_pair(n2_uid,
        base_edge{e_id, edge_value_type()}));
      edgemap_[n2_uid].insert(std::make_pair(n1_uid,
        base_edge{e_id, edge_value_type()}));
    }
    return e;
  }

  /** Remove edge defined by two nodes from graph
   * @param[in,out] _a_ const reference to first node in edge
   * @param[in,out] _b_ const reference to second node in edge
   * @return        1 if removal was successful, else 0
   *
   * @pre _a_ and _b_ are both valid nodes in the graph
   * @pre has_edge(_a_, _b_) == true
   *
   * @post num_edges() [post] == num_edges() [pre] - 1
   * @post num_nodes() [post] == num_nodes() [pre]
   * @post has_edge(_a_, _b_) == false
   * @post member operations on Edge(_a_, _b_) will be invalid
   *
   * Note: Edge indices are invalidated, i.e. edge(i) [pre]
   *   == edge(i) [post] is not guaranteed. This means that
   *   edge_iterator is invalidated by the operation.
   *   incident_iterators on _a_ and_b_ are also invalidated.
   *   However, node_iterator and incident_iterator on other
   *   nodes remain valid. If using remove_edge() in an iterator loop,
   *   call iterator.begin() to reset the iterator, or use the iterator
   *   overload remove_edge(edge_iterator). The latter can be used to
   *   set the forward iterator to a valid position in the same loop.
   *
   * Complexity: O(1) average [std::unordered_map access and erase]
   */
  size_type remove_edge(const Node& a, const Node& b) {

    if (!has_edge(a, b)) return 0;
    size_type a_uid = i2u_[a.index()];
    size_type b_uid = i2u_[b.index()];

    base_edge e1_base = edgemap_[a_uid][b_uid];

    // move last edge (en) to current edge position in edges_
    auto en_node_uids = edges_[num_edges() - 1];
    edges_[e1_base.uid] = en_node_uids;

    // update edges_ in edgemap_
    edgemap_[en_node_uids.first][en_node_uids.second].uid = e1_base.uid;
    edgemap_[en_node_uids.second][en_node_uids.first].uid = e1_base.uid;

    // remove edge entries in edgemap_
    edgemap_[a_uid].erase(b_uid);
    edgemap_[b_uid].erase(a_uid);

    // remove last edge in edges_
    edges_.pop_back();

    return 1;
  }

  /** Remove given edge from the graph
   * @param _e_ const reference to edge in the graph
   * @return 1 if removal was successful, else 0
   *
   * Calls remove_edge(n1, n2) on nodes in _e_. See above
   * for complete specification.
   *
   * Complexity: O(1) average
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove edge from graph given iterator to edge
   * @param _e_it_ edge_iterator to edge we wish to remove
   * @return iterator to valid edge in graph;
   *
   * Calls remove_edge(n1, n2) on nodes at (*_e_it_). See above
   * for complete specification. Returns the same iterator due to
   * the swap and pop method (failing to set the iterator to this
   * return value may result in invalid iteration - i.e. skipping)
   *
   * Complexity: O(1) average
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge((*e_it).node1(), (*e_it).node2());
    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   *
   * Complexity: O(nodes_.size() + num_edges())
   */
  void clear() {
    nodes_.clear(); // empty the nodes
    i2u_.clear();   // empty the idx to uid mapping
    edges_.clear(); // empty the edges
    edgemap_.clear(); // emty the edge values
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

    /** Dereference Operator
     * @brief dereference operator for the NodeIterator
     * @return Node at the current location of the NodeIterator
     *
     * Complexity: O(1)
     */
    Node operator*() const {
      return Node(graph_, graph_->i2u_[idx_]);
    }

    /** Increment Operator
     * @brief prefix increment operator for the NodeIterator
     * @return reference to NodeIterator at the next position
     *
     * @post prev_iter != ++iter
     *
     * The iterator progresses in order of Node index by convention. Note
     * that incrementing past the end of the container is allowed but behavior
     * is undefined.
     *
     * Complexity: O(1)
     */
    NodeIterator& operator++() {
      idx_++;
      return *this;
    }

    /** Equality Operator
     * @brief equality operator for the NodeIterator
     * @param _node_iter_ const reference to 2nd node iterator for comparison
     * @return True if the Graph and current Node index match
     *
     * Complexity: O(1)
     */
    bool operator==(const NodeIterator& node_iter) const {
      bool graph_match = this->graph_ == node_iter.graph_;
      bool idx_match = this->idx_ == node_iter.idx_;
      return graph_match && idx_match;
    }

   private:
    friend class Graph; // to generate Nodes

    /** NodeIterator Constructor
     * @param[in] _graph_ pointer to parent graph
     * @param[in] _idx_ index of the node in _graph_.i2u_
     */
    NodeIterator(const Graph* graph, size_type idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx) {}
    Graph* graph_;
    size_type idx_;
  };

  /** Iterator to first Node in the Graph
   * @brief iterator to the beginning of the set of nodes in the graph
   * @return node_iterator to the first node (index 0)
   *
   * Complexity: O(1)
   */
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  /** Iterator to end of Graph's Nodes
   * @brief iterator to the end of the set of nodes in the graph
   * @return node_iterator to the end of the set of nodes
   *
   * Complexity: O(1)
   */
  node_iterator node_end() const {
    return node_iterator(this, this->size());
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

    /** Dereference Operator
     * @brief dereference operator for the IncidentIterator
     * @return Edge at the current location of the IncidentIterator
     *
     * @post (*iter).node1() == iter.node1_ (the "base" node)
     *
     * Complexity: O(1) on average
     */
    Edge operator*() const {
      Node node2 = Node(node1_.graph_, (*e_iter_).first);
      return Edge(node1_, node2);
    }

    /** Increment Operator
     * @brief prefix increment operator for the IncidentIterator
     * @return reference to IncidentIterator at the next position
     *
     * @post prev_iter != ++iter
     *
     * The iterator progresses through the vector of higher-index incident
     * nodes before that of lower index nodes. Boundary behavior depends on
     * on the std::vector iterator.
     *
     * Complexity: O(1)
     */
    IncidentIterator& operator++() {
      ++e_iter_;
      return *this;
    }

    /** Equality Operator
     * @brief equality operator for the IncidentIterator
     * @param _i_iter_ const reference to 2nd incident iterator for comparison
     * @return True if the underlying std::vector iterators match
     *
     * Complexity: O(1)
     */
    bool operator==(const IncidentIterator& i_iter) const {
      return this->e_iter_ == i_iter.e_iter_;
    }

   private:
    friend class Graph;
    /** IncidentIterator Constructor
     * @param[in] _node_ base node from which to traverse incident edges
     * @param[in] _iter_ adjacency list (unordered_map) iterator for the node
     */
    IncidentIterator(const Node& node, adj_iterator_ iter)
      : node1_(node)
      , e_iter_(iter)
    {}
    Node node1_;
    adj_iterator_ e_iter_;
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

    /** Dereference Operator
     * @brief dereference operator for the EdgeIterator
     * @return Edge at the current location of the EdgeIterator
     *
     * @post (*it).node1() < (*it).node2() [specific to this implementation]
     *
     * Complexity: O(1)
     */
    Edge operator*() const {
      std::pair<size_type, size_type> n_idxs = graph_->edges_[idx_];
      Node n1 = Node(graph_, n_idxs.first);
      Node n2 = Node(graph_, n_idxs.second);
      return Edge(n1, n2);
    }

    /** Increment Operator
     * @brief prefix increment operator for the EdgeIterator
     * @return reference to EdgeIterator at the next position
     *
     * @post prev_iter != ++iter
     *
     * Complexity: O(1)
     */
    EdgeIterator& operator++() {
      idx_++;
      return *this;
    }

    /** Equality Operator
     * @brief equality operator for the EdgeIterator
     * @param _e_iter_ const reference to second edge iterator for comparison
     * @return True if the underlying graph and edge index match
     *
     * Complexity: O(1)
     */
    bool operator==(const EdgeIterator& e_iter) const {
      return this->graph_ == e_iter.graph_ && this->idx_ == e_iter.idx_;
    }

   private:
    friend class Graph;
    /** EdgeIterator Constructor
     * @param[in] _graph_ pointer to parent graph for the edges
     * @param[in] _e_idx_ index of the edge in _graph_.edges_
     */
    EdgeIterator(const Graph* graph, size_type e_idx)
      : graph_(const_cast<Graph*>(graph))
      , idx_(e_idx)
    {}
    Graph* graph_;
    size_type idx_;
  };

  /** Iterator to first Edge in the Graph
   * @brief iterator to the beginning of the set of edges in the graph
   * @return edge_iterator to the first edge
   *
   * @post if num_edges() == 0, edge_begin() == edge_end()
   *
   * Complexity: O(1)
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Iterator to end of Graph's Edges
   * @brief iterator to the end of the set of edges in the graph
   * @return edge_iterator to the end of the set of nodes
   *
   * @post if num_edges() == 0, edge_end() == edge_begin()
   *
   * The edge iterator is at its end if the underlying node_iterator
   * is equal to its end. We use a placeholder iterator for the edge.
   *
   * Complexity: O(1)
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->num_edges());
  }

};

#endif // CME212_GRAPH_HPP
