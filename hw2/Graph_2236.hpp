#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template  <typename V, typename E>
class Graph {

  typedef V node_value_type;
  typedef E edge_value_type;

 struct internal_node;
 private:

  struct internal_node {
    Point position;
    node_value_type value;
    unsigned index;
    unsigned id_;
    internal_node(const Point& pos, unsigned id, 
                  node_value_type v, unsigned idx): position(pos), id_(id), value(v), index(idx){};
  };

  struct internal_edge {
    unsigned a, b, index;
    float length;
    edge_value_type value;
    internal_edge(): a(0), b(0), index(0), length(0), value() {}
    internal_edge(unsigned a_, unsigned b_, unsigned index_, float l_,
                  edge_value_type v_): a(a_), b(b_), index(index_), length(l_), value(v_) {};
    internal_edge(unsigned a_, unsigned b_, unsigned index_, float l_): a(a_), b(b_), index(index_), length(l_){};
  };

  std::vector<internal_node> nodes;
  std::vector<std::vector<internal_edge>> edge_adjacency_list;
  std::unordered_map<unsigned, internal_edge> edges;
  int curr_edge_index;
  std::vector<unsigned> indexed_node_locs_; // map an "index" of a node to it's index in the list
  
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    nodes = std::vector<internal_node>();
    edge_adjacency_list = std::vector<std::vector<internal_edge>>();
    edges = std::unordered_map<unsigned, internal_edge>();
    indexed_node_locs_ = std::vector<unsigned>();
    curr_edge_index = 0;
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
    size_type id;
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
      id = -1;
    }
    
    /** Return this node's position. */
    const Point& position() const {
      return fetch().position;
    }

    /** Return this node's position. */
    Point& position() {
      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return fetch().index;
    }

    node_value_type& value() {
      return const_cast<node_value_type&>(fetch().value);
    }
    const node_value_type& value() const {
      return fetch().value;
    }
    size_type degree() const {
      return graph_->edge_adjacency_list.at(id).size();
    }
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, 0, index());
    }
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, graph_->edge_adjacency_list.at(id).size(), index());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ && id == n.id);
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
        return id < n.id;
      } else{
        return graph_ < n.graph_;
      }
    }

   private:
    

    internal_node& fetch() const {
      return const_cast<internal_node&>(graph_->nodes[id]);
    }

    Node(const Graph* graph, size_type index): 
      graph_(graph), id(graph_->indexed_node_locs_[index]) {}

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return indexed_node_locs_.size();
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
    Point p = Point(position.x, position.y, position.z);
    size_type index = size();
    size_type id = nodes.size();

    auto n = internal_node(p, id, v, index);

    indexed_node_locs_.push_back(id);
    nodes.push_back(n);
    std::vector<internal_edge> node_adj_list = std::vector<internal_edge>();
    edge_adjacency_list.push_back(node_adj_list);
    return Node(this, index); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.index() >= 0 && n.index() < num_nodes() && n ==  node(n.index())) return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
    Graph* graph_;
    size_type index_;
    bool flip;
    /** Construct an invalid Edge. */
    Edge() {
      index_ = -1;
      graph_ = nullptr;
      flip = false;
    }



    /** Return a node of this Edge */
    Node node1() const {
      if (flip) {
        return Node(graph_, graph_->nodes.at(graph_->edges[index_].b).index);
      }
      return Node(graph_, graph_->nodes.at(graph_->edges[index_].a).index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (flip) {
        return Node(graph_, graph_->nodes[graph_->edges[index_].a].index);
      }
      return Node(graph_, graph_->nodes[graph_->edges[index_].b].index);
    }

    /*get a reference to the value of the edge, much added value*/
    edge_value_type& value () {
      return const_cast<edge_value_type&>(graph_->edges[index_].value);
    }
    /*get a const reference to the value of the edge, such added value*/
    const edge_value_type& value () const {
      return graph_->edges[index_].value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ != e.graph_) return false;
      bool cond_1 = node1() == e.node1() && node2() == e.node2();
      bool cond_2 = node1() == e.node2() && node2() == e.node1();
      return cond_1 || cond_2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_) {
        std::cout<<node1().id<<","<<node2().id<<std::endl;
        if (std::min(node1(), node2()) == std::min(e.node1(), e.node2())) {
          return std::max(node1(), node2()) < std::max(e.node1(), e.node2());
        } else {
          return std::min(node1(), node2()) < std::min(e.node1(), e.node2());
        }
      } else{
        return graph_ < e.graph_;
      }
    }

    double length() const {
      return graph_->edges[index_].length;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    Edge(const Graph* graph, size_type i): graph_(const_cast<Graph*>(graph)), index_(i) {}

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }
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
    for (auto & element : edge_adjacency_list.at(a.id)) {
      if (element.a == b.id || element.b == b.id) {
        auto e = Edge(this, element.index);
        return e;
      }
    }

    ; // map an "index" of a edge to it's index in the list


    size_type a_id = std::min(a.id, b.id);
    size_type b_id = std::max(b.id, a.id);
    unsigned i = curr_edge_index;
    curr_edge_index++;

    float length = norm(a.position() - b.position());

    internal_edge e1 = internal_edge(a_id, b_id, i, length);
    internal_edge e2 = internal_edge(b_id, a_id, i, length);
    
    edge_adjacency_list.at(a_id).push_back(e1);
    edge_adjacency_list.at(b_id).push_back(e2);

    edges[i] = e1;


    auto e = Edge(this, i);
    if (e.node1() != a) e.flip = ! e.flip;
    return e;
  }

  /**make it like the edge never even happend
   * @post if e was in the graph, it will be removed
   *
   * Invalidates all outstanding edge and incidence iterators
   * invalidates all edge references to this edge
   * O(num_nodes)
   * returns 1 if the edge was removed, zero if not
   */
  size_type  remove_edge(const  Edge& e) {

    node_type n1 = e.node1();
    node_type n2 = e.node2();

    return remove_edge(n1, n2);
  }

  /**make it like the edge never even happend
   * @post if there was an edge between n1 and n2 in the graph, it will be removed
   *
   * Invalidates all outstanding edge and incidence  iterators
   * invalidates all edge references to this edge
   * O(num_nodes)
   * returns 1 if the edge was removed, zero if not
   */
  size_type remove_edge(const  Node& n1, const  Node& n2) {
    if (! has_edge(n1, n2)) {
      return 0;
    }
    auto node_1_loc = std::find_if(edge_adjacency_list.at(n1.id).begin(), 
                                   edge_adjacency_list.at(n1.id).end(), [n2](const internal_edge& e) {
      return e.a == n2.id || e.b == n2.id;
    });

    unsigned index = (*node_1_loc).index;

    edge_adjacency_list.at(n1.id).erase(node_1_loc);

    auto node_2_loc = std::find_if(edge_adjacency_list.at(n2.id).begin(), 
                                   edge_adjacency_list.at(n2.id).end(), [n1](const internal_edge& e) {
      return e.a == n1.id || e.b == n1.id;
    });
    edge_adjacency_list.at(n2.id).erase(node_2_loc);

    edges.erase(index);
    return 1;
  }

  /** make it like the node never even happend
   * @post if n was in the graph it will be removed along with any edges incident to it
   *
   * Invalidates all outstanding edge, node, and incidence  iterators
   * invalidates all node references to this node
   * O(1)
   * returns 1 if the node was removed, zero if not
   */
  size_type  remove_node(const  Node & n) {
    if (!has_node(n)) {
      return 0;
    }
    while (edge_adjacency_list.at(n.id).size() > 0) {
      auto index = edge_adjacency_list.at(n.id).at(0).index;
      remove_edge(Edge(this, index));
    }
    auto node_index = n.index();
    auto last_index = size() - 1;
    if (last_index > 0 && node_index != last_index) {
      auto old_id = indexed_node_locs_.at(node_index);
      auto new_id = indexed_node_locs_.at(last_index);
      indexed_node_locs_.at(node_index) = new_id;
      nodes.at(new_id).index = node_index;
    }
    indexed_node_locs_.pop_back();
    nodes.at(n.id).index = -1;
    return 1;
  }

  /** remove all nodes in n_it from the graph
   * @post fgor each n in n_it if n was in the graph it will be removed along with any edges incident to it
   *
   * Invalidates all outstanding edge, node, and incidence  iterators
   * invalidates all node references to this node
   * O(|n_it|)
   */
  node_iterator  remove_node(node_iterator  n_it) {
    remove_node(*n_it);
    return node_begin();
  }

  /** remove all edges in e_it from the graph
   * @post fgor each edge in e_it if n was in the graph it will be removed
   *
   * Invalidates all outstanding edge, node, and incidence  iterators
   * invalidates all edge references
   * O(|e_it|*num_nodes)
   */
  edge_iterator  remove_edge(edge_iterator  e_it) {
    remove_edge(*e_it);
    return edge_begin();
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes = std::vector<internal_node>();
    edge_adjacency_list = std::vector<std::vector<internal_edge>>();
    edges = std::unordered_map<unsigned, internal_edge>();
    curr_edge_index = 0;
    indexed_node_locs_ = std::vector<unsigned>();
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
      auto e = graph_->edge(graph_->edge_adjacency_list.at(graph_->indexed_node_locs_[node_index]).at(index).index);
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
    using internal_iterator = typename std::unordered_map<unsigned, internal_edge>::const_iterator;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /**get current iterate*/
    Edge operator*() const { 
      return Edge(graph_, it->second.index);
    }
    /**check equality*/
    bool  operator==(const EdgeIterator& x) const {
      return x.graph_ == graph_ && x.it == it; 
    }
    /**check inequality*/
    bool  operator!=(const EdgeIterator& x) const {
      return x.graph_ != graph_ || x.it != it; 
    }
    /**increment iterator*/
    EdgeIterator&  operator++() {
      it ++;
      return *this;//  We  return  an  lvalue  reference .10}
    }
    
    EdgeIterator(const Graph* graph, internal_iterator i): graph_(graph), it(i){}

   private:
    friend class Graph;
    const Graph* graph_;
    internal_iterator it;
  };

  edge_iterator edge_begin() const {
    return EdgeIterator(this, edges.begin());
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges.end());
  }

 private:


};

#endif // CME212_GRAPH_HPP
