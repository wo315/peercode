#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <vector>
#include <utility>

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

  // List of Point objects
  std::vector < std::pair<Point,V> > point_list_;

  // List of indices into point_list_ that are accessible to the user
  //   The indices of user_point_list_ represent the user-facing node ID's
  std::vector<unsigned int> user_point_list_;

  // List of tuples, where each tuple represents a pair of indices
  //      into point_list_
  std::vector<std::pair<std::tuple<unsigned int, unsigned int>, E>> edges_list_;

  // List of indices into edges_list_ that are accessible to the user
  //    The indices of user_edges_list_ represent the user-facing edge ID's
  std::vector<unsigned int> user_edges_list_;

  // Nested map that represents all the Points each Point is connected to
  //       {node index: {connected node index: edge index}}
  std::map<unsigned int, std::map<unsigned int, unsigned int>> point_map_;


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;
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
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : point_list_(),
            user_point_list_(),
            edges_list_(),
            user_edges_list_(),
            point_map_() {
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

    Node() : indx_(-1), node_user_indx_(-1), graph_ptr_() {
     // initialize an empty Point
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_ptr_->point_list_[indx_].first;
    }

    Point& position() {
      return graph_ptr_->point_list_[indx_].first;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // return *(std::find(graph_ptr_->user_point_list_.begin(), graph_ptr_->user_point_list_.end(), indx_));
      return node_user_indx_;
    }
  
    /** Return this node's value. */
    node_value_type& value() {
      return graph_ptr_->point_list_[indx_].second;
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      return graph_ptr_->point_list_[indx_].second;
    }

    /** Return this node's degree, i.e. the number of nodes connected
     * to this node.
     */
    size_type degree() const {
      std::map<unsigned int, unsigned int> connected_nodes
                                          = (graph_ptr_->point_map_).at(indx_);
      return connected_nodes.size();
    }

    /** Set iterator to the beginning of the map of nodes connected
     * to current node.
     * 
     * @return IncidentIterator pointing to first edge of current node
     */
    incident_iterator edge_begin() const {
      std::map<unsigned int, unsigned int>::const_iterator
                current_node_iter = (graph_ptr_->point_map_).at(node_user_indx_).begin();
      return IncidentIterator(current_node_iter, node_user_indx_, graph_ptr_);
    }

    /** Set iterator to the end of the map of nodes connected
     * to current node.
     * 
     * @return IncidentIterator pointing to one step beyond the 
     *         last node connected to this node
     */
    incident_iterator edge_end() const {
      std::map<unsigned int, unsigned int>::const_iterator
                  current_node_iter = (graph_ptr_->point_map_).at(node_user_indx_).end();
      return IncidentIterator(current_node_iter, node_user_indx_, graph_ptr_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      bool same_index = (indx_ == n.indx_);
      bool same_graph = (graph_ptr_ == n.graph_ptr_);
      return (same_index && same_graph);
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
      if (graph_ptr_ < n.graph_ptr_) {
        return true;
      }
      else {
        return indx_ < n.indx_;
      }
    }

   private:
    friend class Graph;

    // Index into point_list_ (list of Point objects)
    size_type indx_;

    // User-facing node index
    size_type node_user_indx_;

    // Pointer to current graph
    graph_type* graph_ptr_;

    /** Node constructor */
    Node(size_type indx, const graph_type* graph_ptr)
      : indx_(indx), graph_ptr_(const_cast<graph_type*>(graph_ptr)) {
        node_user_indx_ = *(std::find(graph_ptr_->user_point_list_.begin(),
                                      graph_ptr_->user_point_list_.end(), indx));
      }

    /** Node constructor with user_index specified */
    Node(size_type indx, size_type node_user_indx, const graph_type* graph_ptr)
      : indx_(indx),
        node_user_indx_(node_user_indx),
        graph_ptr_(const_cast<graph_type*>(graph_ptr)) {
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return user_point_list_.size();
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
  Node add_node(const Point& position,
                const node_value_type& val = node_value_type()) {
    // 1. add to point_list_
    std::pair<Point,V> new_pt (position, val);
    this->point_list_.emplace_back(new_pt);
    
    // 2. get new index
    size_type new_indx = size();

    // 3. add to list of user-facing nodes and get user index
    this->user_point_list_.push_back(new_indx);
    size_type user_indx = user_point_list_.size()-1;
    
    // 4. make the new Node
    Node new_node = Node(new_indx, user_indx, this);
    
    // 5. add new point to the point_map_
    std::map<size_type, size_type> empty_map;
    point_map_[new_indx] = empty_map;
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.node_user_indx_ < num_nodes()) && n.graph_ptr_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < num_nodes()) {
      Node new_node(user_point_list_[i], i, this);
      return new_node;
    }
    else {
      return Node();
    }
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

    /** Construct an invalid Edge */
    Edge() : indx_(-1), user_indx_(-1), node1_indx_(-1), node2_indx_(-1), edge_graph_ptr_() {
    }

    /** Return the value of this Edge */
    edge_value_type& value() {
      return std::get<1>(edge_graph_ptr_->edges_list_[indx_]);
    }
    const edge_value_type& value() const {
      return std::get<1>(edge_graph_ptr_->edges_list_[indx_]);
    }

    /** Return the length of this Edge */
    double length() const {
      size_type n1_indx = std::get<0>(std::get<0>(edge_graph_ptr_->edges_list_[indx_]));
      size_type n2_indx = std::get<1>(std::get<0>(edge_graph_ptr_->edges_list_[indx_]));
      return norm(std::get<0>(edge_graph_ptr_->point_list_[n1_indx])
                    - std::get<0>(edge_graph_ptr_->point_list_[n2_indx]));
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(node1_indx_, edge_graph_ptr_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(node2_indx_, edge_graph_ptr_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      Node e_node1 = e.node1();
      Node e_node2 = e.node2();
      if (node1()==e_node1) {
        if (node2()==e_node2) {
          return true;
        }
        else {
          return false;
        }
      }
      else if (node1() == e_node2) {
        if (node2() == e_node1) {
          return true;
        }
        else {
          return false;
        }
      }
      else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (edge_graph_ptr_ < e.edge_graph_ptr_) {
        return true;
      }
      else {
        return indx_ < e.indx_;
      }
    }

   private:
    friend class Graph;

    // Index into edges_list_ (list of pairs of indices to point_list_)
    size_type indx_;

    // User-facing edge index
    size_type user_indx_;

    // Index of Node1 in point_list_
    size_type node1_indx_;
    size_type node2_indx_;

    // Pointer to current graph
    graph_type* edge_graph_ptr_;

    /** Edge constructor that takes in 2 nodes indices */
    Edge(size_type node1_indx, size_type node2_indx,
                               const graph_type* edge_graph_ptr)
      : indx_(-1), // In this case the index is a dummy variable
        node1_indx_(node1_indx),
        node2_indx_(node2_indx),
        edge_graph_ptr_(const_cast<graph_type*>(edge_graph_ptr)) {
    }

    /** Edge constructor that takes unqiue index and constructs user index */
    Edge(size_type indx, size_type node1_indx,
                         size_type node2_indx, const graph_type* edge_graph_ptr)
      : indx_(indx),
        node1_indx_(node1_indx),
        node2_indx_(node2_indx),
        edge_graph_ptr_(const_cast<graph_type*>(edge_graph_ptr)) {
          user_indx_ = *(std::find(edge_graph_ptr_->user_edges_list_.begin(), edge_graph_ptr_->user_edges_list_.end(), indx));
      }
  
    /** Edge constructor that takes in all attributes */
    Edge(size_type indx, size_type user_indx, size_type node1_indx,
                         size_type node2_indx, const graph_type* edge_graph_ptr)
      : indx_(indx),
        user_indx_(user_indx),
        node1_indx_(node1_indx),
        node2_indx_(node2_indx),
        edge_graph_ptr_(const_cast<graph_type*>(edge_graph_ptr)) {
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return user_edges_list_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const { // this is the user i
    if (i < num_edges()) {
      size_type unique_indx = user_edges_list_[i];
      size_type node1_indx = std::get<0>(std::get<0>(edges_list_[unique_indx]));
      size_type node2_indx = std::get<1>(std::get<0>(edges_list_[unique_indx]));
      return Edge(unique_indx, i, node1_indx, node2_indx, this);
    }
    else {
      return Edge();
    }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    std::map<unsigned int, unsigned int> nodes_of_a = point_map_.at(a.indx_);
    if (nodes_of_a.find(b.indx_) == nodes_of_a.end()) {
      return false;
    }
    else {
      return true;
    }
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
    if (not has_edge(a, b)) {

      // 1. add pair of point indices to list of edges
      std::tuple<unsigned int,unsigned int> point_pair (a.indx_, b.indx_);
      edges_list_.emplace_back(point_pair, E());

      // 2. get unique index of the new edge
      size_type edge_indx = edges_list_.size()-1;

      // 3. add the new edge index to the list of edges that a user can access
      //    and get the user facing index
      user_edges_list_.push_back(edge_indx);
      size_type user_indx = user_edges_list_.size()-1;
      
      // 4. add to point_map_ (both directions!)
      point_map_[a.indx_][b.indx_] = edge_indx;
      point_map_[b.indx_][a.indx_] = edge_indx;

      // 5. return a new edge with edge index and the current graph
      return Edge(edge_indx, user_indx, a.indx_, b.indx_, this);
    }
    else {
      // get and return the existing edge
      size_type existing_edge_indx = point_map_[a.indx_][b.indx_];
      return Edge(existing_edge_indx, a.indx_, b.indx_, this);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    point_list_.clear();
    user_point_list_.clear();
    edges_list_.clear();
    user_edges_list_.clear();
    point_map_.clear();
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
    NodeIterator() : current_indx_(), node_iter_graph_ptr_() {
    }

    /** Dereference operator for NodeIterator class
     * @return current node object
     */
    Node operator*() const {
      return Node(node_iter_graph_ptr_->user_point_list_[current_indx_],
                  current_indx_, node_iter_graph_ptr_);
    }

    /** Increment operator for NodeIterator class
     * @return reference to current NodeIterator object
     */
    NodeIterator& operator++() {
        current_indx_ += 1;
        return *this;
    }

    /** == Comparison operator for NodeIterator class
     * @return true if other NodeIterator is the same as this NodeIterator
     */
    bool operator==(const NodeIterator& n) const {
      bool same_graph = (n.node_iter_graph_ptr_ == node_iter_graph_ptr_);
      bool same_index = (n.current_indx_ == current_indx_);
      return (same_graph && same_index);
    }

    /** != Comparison operator for NodeIterator class
     * @return true if other NodeIterator is not equal to this NodeIterator
     */
    bool operator!=(const NodeIterator& n) const {
      bool diff_graph = (n.node_iter_graph_ptr_ != node_iter_graph_ptr_);
      bool diff_index = (n.current_indx_ != current_indx_);
      return (diff_graph || diff_index);
    }

   private:
    friend class Graph;

    // Current node index
    size_type current_indx_;

    // Pointer to current graph
    graph_type* node_iter_graph_ptr_;

    /** NodeIterator constructor */
    NodeIterator(size_type current_indx, const graph_type* node_iter_graph_ptr)
      : current_indx_(current_indx),
        node_iter_graph_ptr_(const_cast<graph_type*>(node_iter_graph_ptr)) {
    }

  };

  /** Initialize a NodeIterator at the start of the nodes
   * 
   * @return NodeIterator pointing to the first node of this graph
   */
  node_iterator node_begin() const {
    return NodeIterator((size_type)0, this);
  }

  /** Initialize a NodeIterator at the end of the nodes
   * @return NodeIterator pointing to the index after the last
   *         node in this graph
   */
  node_iterator node_end() const {
    return NodeIterator(num_nodes(), this);
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
    IncidentIterator() : current_node_iter_(), root_indx_(), incident_graph_ptr_() {
    }
    
    /** Dereference operator for IncidentIterator class
     * @return current edge between this node and the current connected node
     */
    Edge operator*() const {
      size_type edge_indx = current_node_iter_->second;
      size_type inc_node1_indx = std::get<0>(std::get<0>(incident_graph_ptr_->edges_list_[edge_indx]));
      size_type inc_node2_indx = std::get<1>(std::get<0>(incident_graph_ptr_->edges_list_[edge_indx]));

      size_type other_indx;
      if (root_indx_ == inc_node1_indx) {
        other_indx = inc_node2_indx;
      }
      else {
        other_indx = inc_node1_indx;
      }
      
      return Edge(edge_indx, root_indx_, other_indx, incident_graph_ptr_);
    }

    /** Increment operator for IncidentIterator class
     * @return current IncidentIterator pointing to the next connected node
     */
    IncidentIterator& operator++() {
      current_node_iter_++;
      return *this;
    }

    /** == Comparison operator for IncidentIterator class
     * @return true if other IncidentIterator is the 
     *         same as this IncidentIterator
     */
    bool operator==(const IncidentIterator& other) const {
      bool same_graph = (other.incident_graph_ptr_ == incident_graph_ptr_);
      bool same_iter = (other.current_node_iter_ == current_node_iter_);
      bool same_root = (other.root_indx_ == root_indx_);
      return (same_graph && same_iter && same_root);
    }

    /** != Comparison operator for IncidentIterator class
     * @return true if other IncidentIterator is not the 
     *         same as this IncidentIterator
     */
    bool operator!=(const IncidentIterator& other) const {
      bool diff_graph = (other.incident_graph_ptr_ != incident_graph_ptr_);
      bool diff_iter = (other.current_node_iter_ != current_node_iter_);
      bool diff_root = (other.root_indx_ != root_indx_);
      return (diff_graph || diff_iter || diff_root);
    }

   private:
    friend class Graph;

    // Iterator that points to the current node (connected to the root node)
    std::map<size_type, size_type>::const_iterator current_node_iter_;

    // Index of root node in point_list_;
    size_type root_indx_;

    // Pointer to current graph
    graph_type* incident_graph_ptr_;

    /** IncidentIterator constructor */
    IncidentIterator(std::map<size_type, size_type>::const_iterator current_node_iter,
                     size_type root_indx,
                     const graph_type* incident_graph_ptr)
      : current_node_iter_(current_node_iter),
        root_indx_(root_indx),
        incident_graph_ptr_(const_cast<graph_type*>(incident_graph_ptr)) {
    }
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
    EdgeIterator() : edge_indx_(), edge_iter_graph_ptr_() {
    }

    /** Dereference operator for EdgeIterator class
     * @return current Edge object
     */
    Edge operator*() const {
      size_type edge_node1_indx = std::get<0>(std::get<0>(edge_iter_graph_ptr_->edges_list_[edge_indx_]));
      size_type edge_node2_indx = std::get<1>(std::get<0>(edge_iter_graph_ptr_->edges_list_[edge_indx_]));
      return Edge(edge_indx_, edge_node1_indx, edge_node2_indx, edge_iter_graph_ptr_);
    }

    /** Increment operator for EdgeIterator class
     * @return current EdgeIterator pointing to next edge
     */ 
    EdgeIterator& operator++() {
      edge_indx_ += 1;
      return *this;
    }

    /** == Comparison operator for EdgeIterator class
     * @return true if other EdgeIterator is the same as this EdgeIterator
     */
    bool operator==(const EdgeIterator& other) const {
      bool same_graph = (other.edge_iter_graph_ptr_ == edge_iter_graph_ptr_);
      bool same_indx = (other.edge_indx_ == edge_indx_);
      return (same_graph && same_indx);
    }

    /** != Comparison operator for EdgeIterator class
     * @return true if other EdgeIterator is not the same as this EdgeIterator
     */
    bool operator!=(const EdgeIterator& other) const {
      bool diff_graph = (other.edge_iter_graph_ptr_ != edge_iter_graph_ptr_);
      bool diff_indx = (other.edge_indx_ != edge_indx_);
      return (diff_graph || diff_indx);
    }

   private:
    friend class Graph;

    // Current edge index
    size_type edge_indx_;

    // Pointer to current graph
    graph_type* edge_iter_graph_ptr_;

    /** EdgeIterator constructor */
    EdgeIterator(size_type edge_indx, const graph_type* edge_iter_graph_ptr) 
      : edge_indx_(edge_indx), edge_iter_graph_ptr_(const_cast<graph_type*>(edge_iter_graph_ptr)) {
    }
  };

  /** Set iterator to the first edge of this graph.
   * 
   * @return EdgeIterator pointing to first edge of this graph
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(0, this);
  }

  /** Set iterator to one beyond the last edge of this graph.
   * 
   * @return IncidentIterator pointing to one index after the last 
   *         edge of this graph
   */
  edge_iterator edge_end() const {
    return EdgeIterator(num_edges(), this);
  }

  //
  // Node and Edge Removal
  //

  /** Function to remove a Node from this graph
   * @pre n is a Node of this graph
   * @param n reference to Node object to remove from the graph
   * @post n is an invalid Node
   * @return false if no errors occurred
   */
  size_type remove_node(const Node& n) {
    size_type node_indx = n.index();

    // Remove all edges connected to the node
    if (n.degree() != 0) {
      for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
        const node_type n1 = (*ii).node1();
        const node_type n2 = (*ii).node2();
        remove_edge(n1, n2);
      }
    }

    // Remove the node from the list of points accessible to the user
    user_point_list_.erase(user_point_list_.begin()+node_indx);

    return 1;
  }

  /** Function to remove a Node from this graph
   * @pre n_it points to a Node of this graph
   * @param n_it node_iterator object point to Node to remove from this graph
   * @post n_it is an invalid node_iterator
   * @return node_iterator pointing to the beginning of this graph
   */
  node_iterator remove_node(node_iterator n_it) {
    node_type n = (*n_it);
    size_type node_indx = n.index();

    // Remove all edges connected to the node
    for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
      const node_type n1 = (*ii).node1();
      const node_type n2 = (*ii).node2();
      remove_edge(n1, n2);
    }

    // Remove the node from the list of points accessible to the user
    user_point_list_.erase(user_point_list_.begin()+node_indx);

    return user_post_list_.begin();
  }

  /** Function to remove an Edge from this graph
   * @pre n1 and n2 are Nodes of this graph
   * @pre n1 and n2 are connected by an Edge
   * @param n1 Node object, connected to n2 by an Edge
   * @param n2 Node object, connected n1 by an Edge
   * @post the number of edges in this Graph has decreased by 1
   * @return false if no errors occurred
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    size_type n1_indx = n1.index();
    size_type n2_indx = n2.index();
    size_type edge_indx = point_map_[n1_indx][n2_indx];

    // Remove the edge_indx from the list of user_facing id's
    std::vector<unsigned int>::iterator i = std::find(user_edges_list_.begin(), user_edges_list_.end(), edge_indx);
    if (i == user_edges_list_.end()) {
      return 0;
    }
    user_edges_list_.erase(i);

    // Remove the n2:n1 edge pair from n1's map
    point_map_[n1_indx].erase(n2_indx);

    // Remove the n1:n2 edge pair from n2's map
    point_map_[n2_indx].erase(n1_indx);

    return 1;
    
  }

  /** Function to remove an Edge from this graph
   * @pre e is an Edge of this graph
   * @param e Edge object to remove from this graph
   * @post e is an invalid Edge
   * @return false if no errors occurred
   */
  size_type remove_edge(const Edge& e) {
    size_type n1_indx = e.node1_indx_;
    size_type n2_indx = e.node2_indx_;
    size_type edge_indx = point_map_[n1_indx][n2_indx];

    // Remove the edge_indx from the list of user_facing id's
    std::vector<unsigned int>::iterator i = std::find(user_edges_list_.begin(), user_edges_list_.end(), edge_indx);
    if (i == user_edges_list_.end()) {
      return 0;
    }
    user_edges_list_.erase(i);

    // Remove the n2:n1 edge pair from n1's map
    point_map_[n1_indx].erase(n2_indx);

    // Remove the n1:n2 edge pair from n2's map
    point_map_[n2_indx].erase(n1_indx);

    return 1;

  }

  /** Function to remove an Edge from this graph
   * @pre e_it points to an Edge in this graph
   * @param e_it edge_iterator object pointing to an Edge to remove from this graph
   * @post e_it is an invalid edge_iterator
   * @return an edge_iterator pointing to the first Edge in this graph
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    size_type n1_indx = (*e_it).node1_indx_;
    size_type n2_indx = (*e_it).node2_indx_;
    size_type edge_indx = point_map_[n1_indx][n2_indx];

    // Remove the edge_indx from the list of user_facing id's
    std::vector<unsigned int>::iterator i = std::find(user_edges_list_.begin(), user_edges_list_.end(), edge_indx);
    if (i == user_edges_list_.end()) {
      return 0;
    }
    user_edges_list_.erase(i);

    // Remove the n2:n1 edge pair from n1's map
    point_map_[n1_indx].erase(n2_indx);

    // Remove the n1:n2 edge pair from n2's map
    point_map_[n2_indx].erase(n1_indx);

    return 1;
  }

};

#endif // CME212_GRAPH_HPP
