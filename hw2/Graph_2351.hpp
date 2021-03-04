#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

//--functionality_0
//--great job!
//--END

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
 */
template <typename V, typename E>
class Graph {
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

  using node_value_type = V;

  using edge_value_type = E;

 private:
  /** Container that stores positions as well as the index, i.e. node_id_ */
  std::vector<std::pair<Point*, node_value_type>> nodes_;

  /** Container that stores the adjacency in the graph:
   * Example: [[node_1, node_3], [node_2], ..., [node_1, node_4]
   * nodes inside each inner list have an edge with the node whose id is the index of that list:
   * --> edges for the 0-th inner list: (0,1) and (0,3). */
  std::vector<std::vector<std::pair<size_type, edge_value_type>>> edges_;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

 public:

  /** Construct an empty graph. */
  Graph() : nodes_(std::vector<std::pair<Point*, node_value_type>>(0)),
            edges_(std::vector<std::vector<std::pair<size_type, edge_value_type>>>(0)) {}

  /** Default destructor */
  ~Graph() {
      this->clear();
  };

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
    Node() {}

    /** Return this node's position. */
    Point& position() {
        return *(graph_->nodes_[node_id_].first);
    }

    const Point& position() const {
      return *(graph_->nodes_[node_id_].first);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_id_;
    }

    node_value_type& value() {
        return graph_->nodes_[node_id_].second;
    }

    const node_value_type& value() const {
        return graph_->nodes_[node_id_].second;
    }

    size_type degree() const {
        return graph_->edges_[node_id_].size();
    }

    incident_iterator edge_begin() const {
        return {graph_, node_id_, 0};
    }

    incident_iterator edge_end() const {
        return {graph_, node_id_, degree()};
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (graph_ == n.graph_ and node_id_ == n.node_id_)
          return true;
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
      if ((graph_ == n.graph_ and node_id_ < n.node_id_) or
              (graph_ < n.graph_))
          return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type node_id_;

    /* Construct valid node */
    Node(const graph_type* graph, size_type node_id)
              : graph_(const_cast<graph_type*>(graph)), node_id_(node_id) {}
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
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    auto* p = new Point;
    *p = position;
    std::pair<Point*, node_value_type> node {p, node_value_type()};  // HW1
    nodes_.push_back(node);
    std::vector<std::pair<size_type, edge_value_type>> empty_edge;
    edges_.push_back(empty_edge);
    return {this, size()-1};
  }

  Node add_node(const Point& position, const node_value_type& value) {
    auto* p = new Point;
    *p = position;
    std::pair<Point*, node_value_type> node {p, value};
    nodes_.push_back(node);
    std::vector<std::pair<size_type, edge_value_type>> empty_edge;
    edges_.push_back(empty_edge);
    return {this, size()-1};
    }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph_ == this and n.node_id_ < size())
        return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return {this, i};
  }

  /** Remove node @a n from this Graph.
    * @return 1 if node removal is successful, else 0.
    * @post has_node(@a n) = false
    * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1
    *       Else,                  new num_nodes() == old num_nodes()
    * @post If old has_node(@a n), new num_edges() == old num_edges() - n.degree()
    *       Else,                  new num_edges() == old num_edges()
    *
    * Complexity: O(degree()^2)
    * Invalidates: Nodes with node_id_ == n.node_id_ or node_id_ == old num_nodes() - 1
                   Edges incident to n
                   node_iterators representing the invalid nodes
                   edge_iterators representing the invalid edges
    */
  size_type remove_node(const Node& n) {
    if (!(has_node(n)))
        return 0;
    size_type n_id = n.node_id_, last_id = size()-1;
    for (size_type i = 0; i < n.degree(); ++i) {
        size_type n2_id = edges_[n_id][i].first;
        for (size_type j = 0; j < edges_[n2_id].size(); ++j)
            if (edges_[n2_id][j].first == n_id) {
                edges_[n2_id][j] = edges_[n2_id].back();
                edges_[n2_id].pop_back();
                break;
            }
    }
    nodes_[n_id] = nodes_.back();
    nodes_.pop_back();
    edges_[n_id] = edges_.back();
    edges_.pop_back();
    if (n_id < last_id) {
        for (size_type i = 0; i < edges_[n_id].size(); ++i) {
            size_type n2_id = edges_[n_id][i].first;
            for (size_type j = 0; j < edges_[n2_id].size(); ++j)
                if (edges_[n2_id][j].first == last_id) {
                    edges_[n2_id][j].first = n_id;
                    break;
                }
        }
    }
    return 1;
  }

  /** Remove node @a n = (*n_it) from this Graph.
    * @return 1 if node removal is successful, else 0.
    * @post has_node(@a n) = false
    * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1
    *       Else,                  new num_nodes() == old num_nodes()
    * @post If old has_node(@a n), new num_edges() == old num_edges() - n.degree()
    *       Else,                  new num_edges() == old num_edges()
    *
    * Complexity: O(degree()^2)
    * Invalidates: Nodes with node_id_ == n.node_id_ or node_id_ == old num_nodes() - 1
                   Edges incident to n
                   node_iterators representing the invalid nodes
                   edge_iterators representing the invalid edges
    */
  node_iterator remove_node(node_iterator n_it) {
      if (n_it != node_end())
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
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return {graph_, node1_id_};
    }

    /** Helper function: returns node2_id given node1_id according to edges in graph */
    size_type GetNode2ID() const {
        return graph_->edges_[node1_id_][node2_loc_].first;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return {graph_, GetNode2ID()};
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ == e.graph_ and
              ((node1_id_ == e.node1_id_ and GetNode2ID() == e.GetNode2ID()) or
                      (GetNode2ID() == e.node1_id_ and node1_id_ == e.GetNode2ID())))
          return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      size_type this_min = std::min(node1_id_, GetNode2ID());
      size_type e_min = std::min(e.node1_id_, e.GetNode2ID());
      size_type this_max = std::max(node1_id_, GetNode2ID());
      size_type e_max = std::max(e.node1_id_, e.GetNode2ID());
      if ((graph_ < e.graph_) or
              ((graph_ == e.graph_) and
              ((this_min < e_min) or ((this_min == e_min) and (this_max < e_max))))) {
          return true;
      }
      return false;
    }

    double length() const {
        return norm(node1().position() - node2().position());
    }

    Edge complement() const {
        for (auto ii = node2().edge_begin(); ii != node2().edge_end(); ++ii) {
            if ((*ii).node2() == this->node1())
                return *ii;
        }
        return Edge();
    }

    edge_value_type& value() {
        return graph_->edges_[node1_id_][node2_loc_].second;
    }

    const edge_value_type& value() const {
        return graph_->edges_[node1_id_][node2_loc_].second;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type node1_id_;
    size_type node2_loc_;

    Edge (const graph_type* graph, size_type node1_id, size_type node2_loc)
        : graph_(const_cast<graph_type*>(graph)), node1_id_(node1_id), node2_loc_(node2_loc) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
      size_type count = 0;
      for (const auto& e : edges_)
          count += e.size();
      return count / 2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
      return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type a_id = a.node_id_, b_id = b.node_id_;
    size_type num_neighbors = edges_[a_id].size();
    for (size_type i = 0; i < num_neighbors; i++)
        if (edges_[a_id][i].first == b_id)
            return true;
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
    size_type a_id = a.node_id_, b_id = b.node_id_;

    size_type num_neighbors = edges_[a_id].size();
    if (this->has_edge(a, b))
        for (size_type i = 0; i < num_neighbors; i++)
            if (edges_[a_id][i].first == b_id)
                return {this, a_id, i};

    std::pair<size_type, edge_value_type> b_pair(b_id, edge_value_type());
    std::pair<size_type, edge_value_type> a_pair(a_id, edge_value_type());
    edges_[a_id].push_back(b_pair);
    edges_[b_id].push_back(a_pair);
    return {this, a_id, static_cast<size_type>(edges_[a_id].size()-1)};
  }

  /** Remove edge from this Graph.
    * @return 1 if edge removal is successful, else 0.
    * @post has_edge(@a a, @a b) == false
    * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1
    *       Else,                        new num_edges() == old num_edges()
    *
    * Complexity: O(a.degree() + b.degree())
    * Invalidates: Edge(graph, a.node_id_, old a.degree() - 1);
                   Edge(graph, b.node_id_, old b.degree() - 1);
                   edge_iterators representing the invalid edges;
    */
  size_type remove_edge(const Node& a, const Node& b) {
    if (!(has_edge(a, b)))
        return 0;
    size_type a_id = a.node_id_, b_id = b.node_id_;
    for (size_type i = 0; i != a.degree(); ++i)
        if (edges_[a_id][i].first == b_id) {
            edges_[a_id][i] = edges_[a_id].back();
            edges_[a_id].pop_back();
            break;
        }
    for (size_type i = 0; i != b.degree(); ++i)
        if (edges_[b_id][i].first == a_id) {
            edges_[b_id][i] = edges_[b_id].back();
            edges_[b_id].pop_back();
            break;
        }
    return 1;
  }

  /** Remove edge from this Graph.
    * @return 1 if edge removal is successful, else 0.
    * @post has_edge(@a e.node1(), @a e.node2()) == false
    * @post If old has_edge(@a e.node1(), @a e.node2()), new num_edges() == old num_edges() - 1
    *       Else,                                        new num_edges() == old num_edges()
    *
    * Complexity: O(e.node1().degree() + e.node2().degree())
    * Invalidates: Edge(graph, a.node_id_, old a.degree() - 1);
                   Edge(graph, b.node_id_, old b.degree() - 1);
                   edge_iterators representing the invalid edges;
    */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove edge from this Graph.
    * @return 1 if edge removal is successful, else 0.
    * @post has_edge(@a (*e_it).node1(), @a (*e_it).node2()) == false
    * @post If old has_edge(@a (*e_it).node1(), @a (*e_it).node2()), new num_edges() == old num_edges() - 1
    *       Else,                                                    new num_edges() == old num_edges()
    *
    * Complexity: O((*e_it).node1().degree() + (*e_it).node2().degree())
    * Invalidates: Edge(graph, a.node_id_, old a.degree() - 1);
                   Edge(graph, b.node_id_, old b.degree() - 1);
                   edge_iterators representing the invalid edges;
    */
  edge_iterator remove_edge(edge_iterator e_it) {
    if (e_it != edge_end())
      remove_edge(*e_it);
    return e_it.validate();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    for (auto node : nodes_)
        delete node.first;
    nodes_.clear();
    edges_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
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

    Node operator*() const {
        return {graph_, cur_node_};
    }

    NodeIterator& operator++() {
        ++cur_node_;
        return *this;
    }

    bool operator==(const NodeIterator& node_iter) const {
        return (this->graph_ == node_iter.graph_ && this->cur_node_ == node_iter.cur_node_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type cur_node_;

    NodeIterator(const graph_type* graph, size_type cur_node) :
                graph_(const_cast<graph_type*>(graph)), cur_node_(cur_node) {}
  };

  node_iterator node_begin() const {
      return {this, 0};
  }

  node_iterator node_end() const {
      return {this, num_nodes()};
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
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

    Edge operator*() const {
        return {graph_, cur_node_, cur_adj_};
    }

    IncidentIterator& operator++() {
        ++cur_adj_;
        return *this;
    }

    bool operator==(const IncidentIterator& incident_iter) const {
        return (this->graph_ == incident_iter.graph_ &&
                this->cur_node_ == incident_iter.cur_node_ &&
                this->cur_adj_ == incident_iter.cur_adj_);
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type cur_node_;
    size_type cur_adj_;

    IncidentIterator(const graph_type* graph, size_type cur_node, size_type cur_adj) :
                    graph_(const_cast<graph_type*>(graph)), cur_node_(cur_node), cur_adj_(cur_adj) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
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

    Edge operator*() const {
        return *ii_;
    }

    EdgeIterator& operator++() {
        ++ii_;
        while (ni_ != graph_->node_end()) {
            while (ii_ != (*ni_).edge_end()) {
                if (ii_.cur_node_ < graph_->edges_[ii_.cur_node_][ii_.cur_adj_].first)
                    return *this;
                ++ii_;
            }
            ++ni_;
            ii_ = {graph_, ni_.cur_node_, 0};
        }
        return *this;
    }

    bool operator==(const EdgeIterator& edge_iter) const {
        return (this->graph_ == edge_iter.graph_ &&
                this->ni_ == edge_iter.ni_ &&
                this->ii_ == edge_iter.ii_);
    }

    EdgeIterator& validate() {
        while (ni_ != graph_->node_end()) {
            while (ii_ != (*ni_).edge_end()) {
                if (ii_.cur_node_ < graph_->edges_[ii_.cur_node_][ii_.cur_adj_].first)
                    return *this;
                ++ii_;
            }
            ++ni_;
            ii_ = {graph_, ni_.cur_node_, 0};
        }
        return *this;
    }

   private:
    friend class Graph;
    graph_type* graph_;
    node_iterator ni_;
    incident_iterator ii_;

    EdgeIterator(const graph_type* graph, node_iterator ni, incident_iterator ii) :
                graph_(const_cast<graph_type*>(graph)), ni_(ni), ii_(ii) {}
  };

  edge_iterator edge_begin() const {
      return {this, this->node_begin(), (*this->node_begin()).edge_begin()};
  }

  edge_iterator edge_end() const {
      return {this, this->node_end(), (*this->node_end()).edge_begin()};
  }
};

#endif // CME212_GRAPH_HPP
