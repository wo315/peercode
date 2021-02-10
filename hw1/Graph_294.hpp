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
 */
template  <typename V>
class Graph {
 struct internal_node;
 private:

  struct internal_node {
    const Point& position;
    std::string uid;      // The unique identifcation for an element
    internal_node(const Point& pos, std::string id): position(pos), uid(id) {};
  };

  struct internal_edge {
    unsigned a, b, index;
    internal_edge(unsigned a_, unsigned b_, unsigned index_): a(a_), b(b_), index(index_) {};
  };

  std::vector<Point> node_list;
  std::vector<V> values;
  std::vector<std::vector<internal_edge>> edge_adjacency_list;
  std::vector<internal_edge> edges;
  
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using  node_value_type = V;

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
  Graph() {
    node_list = std::vector<Point>();
    edge_adjacency_list = std::vector<std::vector<internal_edge>>();
    edges = std::vector<internal_edge>();
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
    const Graph* graph_;
    size_type index_;
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
      graph_ = nullptr;
      index_ = -1;
    }
    
    /** Return this node's position. */
    const Point& position() const {
      return fetch();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    node_value_type& value() {
      return const_cast<node_value_type&>(graph_->values.at(index_));
    }
    const node_value_type& value() const {
      return graph_->values.at(index_);
    }
    size_type degree() const {
      return graph_->edge_adjacency_list.at(index_).size();
    }
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, 0, index_);
    }
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, graph_->edge_adjacency_list.at(index_).size(), index_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && index_ == n.index_);
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
      if (graph_ == n.graph_) {
        return index_ < n.index_;
      } else{
        return graph_ < n.graph_;
      }
    }

   private:
    

    Point& fetch() const {
      return const_cast<Point&>(graph_->node_list.at(index_));
    }

    Node(const Graph* graph, size_type index): 
      graph_(graph), index_(index) {}

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_list.size();
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
  Node  add_node(const  Point& position, const  node_value_type& v = node_value_type ()) {
    const Point p = Point(position.x, position.y, position.z);
    size_type index = size();
    std::string node_uid = std::to_string(position.x) + std::to_string(position.y) + std::to_string(position.z);
    internal_node new_node = internal_node(p, node_uid);
    node_list.push_back(p);
    std::vector<internal_edge> node_adj_list = std::vector<internal_edge>();
    edge_adjacency_list.push_back(node_adj_list);
    values.push_back(v);
    return Node(this, index); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    (void) n;
    if (n.index() < num_nodes() && n ==  node(n.index())) return true;
    return false;
    //             // Quiet compiler warning
    // return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    (void) i;             // Quiet compiler warning
    return Node(this, i);        // Invalid node
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
    const Graph* graph_;
    size_type index;
    bool flip;
    /** Construct an invalid Edge. */
    Edge() {
      index = -1;
      graph_ = nullptr;
      flip = false;
    }

    /** Return a node of this Edge */
    Node node1() const {
      if (flip) {
        return Node(graph_, graph_->edges.at(index).b);      // Invalid Node
      }
      return Node(graph_, graph_->edges.at(index).a);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (flip) {
        return Node(graph_, graph_->edges.at(index).a);      // Invalid Node
      }
      return Node(graph_, graph_->edges.at(index).b);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ != e.graph_) return false;
      bool cond_1 = node1() == e.node1() && node2() == e.node2();
      bool cond_2 = node1() == e.node2() && node2() == e.node1();
      (void) e;           // Quiet compiler warning
      return cond_1 || cond_2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_) {
        if (std::min(node1(), node2()) == std::min(e.node1(), e.node2())) {
          return std::max(node1(), node2()) < std::max(e.node1(), e.node2());
        } else {
          return std::min(node1(), node2()) < std::min(e.node1(), e.node2());
        }
      } else{
        return graph_ < e.graph_;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Edge(const Graph* graph, size_type i): graph_(graph), index(i) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i);        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (auto & element : edge_adjacency_list.at(a.index())) {
      if (element.a == b.index() || element.b == b.index()) return true;
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
    (void) a, (void) b;   // Quiet compiler warning
    for (auto & element : edge_adjacency_list.at(a.index())) {
      if (element.a == b.index() || element.b == b.index()) {
        auto e = Edge(this, element.index);
        if (e.node1() != a) e.flip = ! e.flip;
        return e;
      }
    }

    size_type a_ind = std::min(a.index(),b.index());
    size_type b_ind = std::max(b.index(),a.index());
    size_type i = edges.size();

    internal_edge e1 = internal_edge(a_ind, b_ind, i);
    internal_edge e2 = internal_edge(b_ind, a_ind, i);
    
    edge_adjacency_list.at(a_ind).push_back(e1);
    edge_adjacency_list.at(b_ind).push_back(e2);

    edges.push_back(e1);

    auto e = Edge(this, i);
    if (e.node1() != a) e.flip = ! e.flip;
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_list = std::vector<Point>();
    edge_adjacency_list = std::vector<std::vector<internal_edge>>();
    edges = std::vector<internal_edge>();
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

    /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
    Node operator*() const { 
      return graph_->node(index);
    }

    /** check equality*/
    bool  operator==(const NodeIterator& x) const {
      return x.graph_ == graph_ && x.index == index; 
    }

    /** check inequality*/
    bool  operator!=(const NodeIterator& x) const {
      return x.graph_ != graph_ || x.index != index; 
    }
    
    /** increment iterator*/
    NodeIterator&  operator++() {
      index ++;
      return *this;//  We  return  an  lvalue  reference .10}
    }

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    NodeIterator(const Graph* graph, size_type i): graph_(graph), index(i) {}

   private:
    friend class Graph;
    const Graph* graph_;
    size_type index;
  };

  /** get start of iterater through nodes*/
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** get end of iterater through nodes*/
  node_iterator node_end() const {
    return NodeIterator(this, size());
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

    /**get current itterate*/
    Edge operator*() const {
      auto e = graph_->edge(graph_->edge_adjacency_list.at(node_index).at(index).index);
      if (e.node1() != graph_->node(node_index)) e.flip = !e.flip;
      return e;
    }
    /** increment iterator*/
    IncidentIterator& operator++() {
      index++;
      return *this;
    }
    /**check equality*/
    bool  operator==(const IncidentIterator& x) const {
      return x.graph_ == graph_ && x.index == index && node_index == x.node_index; 
    }
    /**check inequality*/
    bool  operator!=(const IncidentIterator& x) const {
      return x.graph_ != graph_ || x.index != index || node_index != x.node_index;
    }
    
    IncidentIterator(const Graph* graph, size_type i, size_type n_i): graph_(graph), index(i), node_index(n_i) {}
   private:
    friend class Graph;
    const Graph* graph_;
    size_type index;
    size_type node_index;

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

    /**get current iterate*/
    Edge operator*() const { 
      return graph_->edge(index);
    }
    /**check equality*/
    bool  operator==(const EdgeIterator& x) const {
      return x.graph_ == graph_ && x.index == index; 
    }
    /**check inequality*/
    bool  operator!=(const EdgeIterator& x) const {
      return x.graph_ != graph_ || x.index != index; 
    }
    /**increment iterator*/
    EdgeIterator&  operator++() {
      index ++;
      return *this;//  We  return  an  lvalue  reference .10}
    }
    
    EdgeIterator(const Graph* graph, size_type i): graph_(graph), index(i) {}

   private:
    friend class Graph;
    const Graph* graph_;
    size_type index;
  };

  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges.size());
  }

 private:


};

#endif // CME212_GRAPH_HPP
