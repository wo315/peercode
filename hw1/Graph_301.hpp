#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <list>  // for incident edges
#include <stdexcept>
#include <unordered_map> // for nodes
#include <vector>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V=int>
class Graph {

 public:
  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  /** Type of additional attribute given to nodes (added for HW1) */
  using node_value_type = V;

 private:

  // graph node with uid and position
  struct base_node {

    Graph* p_graph; // pointer to graph (indicator that node is valid)
    size_type uid;  // unique index of node in graph
    Point position; // physical location of node
    node_value_type value; // user specified value

    std::vector<size_type> adj_list; // list of higher-index nodes
    std::vector<size_type> adj_list_lower; // list of lower-index nodes

    base_node(const Graph* graph, size_type id, Point pos, node_value_type val)
      : p_graph(const_cast<Graph*>(graph))
      , uid(id)
      , position(pos)
      , value(val)
      , adj_list()
      , adj_list_lower()
    { }

  };

  // note: maps used for nodes & edges so that location in memory of members
  // is fixed, and pointer reference and deletion are simpler
  std::unordered_map<size_type, base_node> nodes_; // map indices nodes/edges
  std::unordered_map<size_type, std::pair<size_type, size_type>> edges_;

  // define typenames for adjacency list iterators
  using edgelist_iterator_ =
    typename std::vector<size_type>::const_iterator;

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
  Graph() : nodes_(), edges_() {
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // note: behavior for invalid nodes is undefined
      return node_->position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // note: behavior for invalid nodes is undefined
      return node_->uid;
    }

    /**
     * @brief retrieve the value held by the node
     * @return node_value_type reference for the value held by the node
     *
     * Complexity: O(1)
     */
    node_value_type& value() {
      return node_->value;
    }

    /**
     * @brief const version of value()
     * @return reference to const node_value_type value held by the node
     *
     * Complexity: O(1)
     */
    const node_value_type& value() const {
      return const_cast<const node_value_type&>(node_->value);
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
      return node_->adj_list.size() + node_->adj_list_lower.size();
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
      if (node_->adj_list.empty()) {
        return IncidentIterator(*this, node_->adj_list_lower.begin());
      }
      else {
        return IncidentIterator(*this, node_->adj_list.begin());
      }
    }

    /** Iterator to the end of the set of incident edges
     *
     * @return incident_iterator that represents the end of the incident set
     * @post if the node has no incident edges, edge_begin() == edge_end()
     *
     * Complexity: O(1)
     */
    incident_iterator edge_end() const {
      return IncidentIterator(*this, node_->adj_list_lower.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index
     */
    bool operator==(const Node& n) const {
      return (this->index() == n.index()) && (graph_ == n.graph_);
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
      if (graph_ == n.graph_) return this->index() < n.index();
      else return (graph_ < n.graph_); // arbitrary comparison to break ties
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph<V>;

    /** Private Node Constructor */
    Node(const Graph* graph, const base_node* node )
      : graph_(const_cast<Graph*>(graph))
      , node_(const_cast<base_node*>(node))
    { }

    Graph* graph_;  // pointer to the first node in the container
    base_node* node_;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value (defined by graph template)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) on average
   */
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type()) {
    size_type uid = this->num_nodes();
    nodes_.insert(std::make_pair(uid,
      base_node(this, uid, const_cast<Point&>(position), value)));
    return Node(this, &nodes_.at(uid));
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.node_->p_graph != this) return false; // node must be valid
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
    return Node(this, &nodes_.at(i));
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

    /** Return a node of this Edge */
    Node node1() const {
      return node1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      // sort nodes to compare undirected edges
      auto nodes1 = std::minmax({this->node1(), this->node2()});
      auto nodes2 = std::minmax({e.node1(), e.node2()});
      return (nodes1.first == nodes2.first) && (nodes1.second == nodes2.second);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
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
    friend class Graph;

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
   * Complexity: O(1) average
   */
  Edge edge(size_type i) const {

    if (i >= num_edges()) {
      throw std::out_of_range("accessing out of bounds edge");
    }
    return Edge(node(edges_.at(i).first), node(edges_.at(i).second));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: Average O(num_edges() / num_nodes())
   */
  bool has_edge(const Node& a, const Node& b) const {
    // make sure nodes belong to the graph in question
    if (a == b) return false;
    if (!this->has_node(a) || !this->has_node(b)) return false;

    // sort nodes and search outgoing edge from lower-index node
    auto e_nodes = std::minmax({a, b});
    for (size_type id : nodes_.at(e_nodes.first.index()).adj_list)
      if (id == e_nodes.second.index()) return true;
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
   * Complexity: O(num_edges() / num_nodes()) due to check for edge existence
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
      size_type idx1 = e_nodes.first.index();
      size_type idx2 = e_nodes.second.index();
      edges_.insert(std::make_pair(this->num_edges(),
        std::make_pair(idx1, idx2)));
      nodes_.at(idx1).adj_list.push_back(idx2);
      nodes_.at(idx2).adj_list_lower.push_back(idx1);
    }
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear(); // empty the nodes
    edges_.clear(); // empty the edges
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
     * Complexity: O(1) on average (unordered_map access)
     */
    Node operator*() const {
      return Node(graph_, const_cast<const base_node*>(
        &graph_->nodes_.at(idx_)));
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
      // node keys should be all integers in [0, num_nodes). Should incorporate
      // a defensive check in case of later node removal.
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
     * @param[in] _idx_ index of the node in graph.nodes_
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
      return Edge(node1_, node1_.graph_->node(*e_iter_));
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
      if (++e_iter_ == node1_.node_->adj_list.end()) {
        e_iter_ = node1_.node_->adj_list_lower.begin();
      }
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
     * @param[in] _iter_ adjacency list (vector) iterator for the node
     */
    IncidentIterator(const Node& node, edgelist_iterator_ iter)
      : node1_(node)
      , e_iter_(iter)
    {}
    Node node1_;
    edgelist_iterator_ e_iter_;
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
     * Complexity: O(1) on average
     */
    Edge operator*() const {
      std::pair<size_type, size_type> n_idxs = graph_->edges_.at(idx_);
      return Edge(graph_->node(n_idxs.first), graph_->node(n_idxs.second));
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
     * @param[in] _e_idx_ index of the edge in graph.edges_
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
