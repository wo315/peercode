#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** Synonym for node value type*/
  using node_value_type = V;
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
      size_all_node = 0;
      size_all_edge = 0;
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
      size_type uindex() const {
          return node_id_;
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

    /** Return this node's position. */
    const Point& position() const {
        return g_->nodes_[node_id_].p_;
    }
      
    Point& position() {
        return g_->nodes_[node_id_].p_;
    }
      
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return g_->nodes_[node_id_].idx_;
    }

      node_value_type& value() {
          return g_->nodes_[node_id_].v_;
      }
      const node_value_type& value() const {
          return g_->nodes_[node_id_].v_;
      }
      size_type degree() const {
          if (g_->neighbours.count(node_id_)) {
              return g_->neighbours.at(node_id_).size();
          }
          else {
              return 0;
          }
      }
      incident_iterator edge_begin() const {
         return IncidentIterator(g_, *this, 0);
      }
      incident_iterator edge_end() const {
          size_type size = 0;
          if (g_->neighbours.count(node_id_))
              size = g_->neighbours[node_id_].size();
          return IncidentIterator(g_, *this, size);
      }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (g_ == n.g_ && node_id_ == n.node_id_);
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
        if ((g_ == n.g_ && node_id_ < n.node_id_) || g_ != n.g_) {
          return true;
        }
      return false;
    }
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* g_;
    // this is the unique id of the node
    size_type node_id_;
    Node(const graph_type* graph, size_type node_id)
      : g_(const_cast<graph_type*>(graph)), node_id_(node_id) {}
  };

  /** Return the neighbours of a node
   */
    std::vector<Node> getNeighbours(node_type node) {
        std::vector<Node> v;
        if (neighbours.count(node.uindex()) && neighbours[node.uindex()].size() > 0) {
            for (size_type i: neighbours[node.uindex()]) {
                v.push_back(vectorNodes[i]);
            }
        }
        return v;
    }
  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_i2u_.size();
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    Node node = Node(this, size_all_node);
    vectorNodes.push_back(node);
    nodes_.push_back(nodeinfo{position, value, size()});
    node_i2u_.push_back(size_all_node);
    ++size_all_node;
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.uindex() <  size_all_node && nodes_[n.uindex()].idx_ != (size_type) -1)
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
    assert(i < size());
    return vectorNodes[node_i2u_[i]];
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
  private:
      size_type uindex() const {
          return edge_uidx;
      }
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
        return g_->vectorNodes[node1_uidx];
    }

    /** Return the other node of this Edge */
    Node node2() const {
        return g_->vectorNodes[node2_uidx];
    }
    
    const size_type& index() const {
        return g_->edges_[edge_uidx].idx_;
    }
      
    double length() const  {
        return norm(node1().position() - node2().position());
    }
      
    edge_value_type& value() {
        return g_->edges_[edge_uidx].v_;
    }
    const edge_value_type& value() const {
        return g_->edges_[edge_uidx].v_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (g_ == e.g_ && node1_uidx == e.node1_uidx && node2_uidx == e.node2_uidx);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
      bool operator<(const Edge& e) const {
          if ((g_ == e.g_ && edge_uidx < e.edge_uidx) || g_ != e.g_) {
            return true;
          }
        return false;
      }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
      graph_type* g_;
      size_type edge_uidx;
      size_type node1_uidx;
      size_type node2_uidx;
      Edge(graph_type* graph, size_type edge_idx, size_type node1_idx, size_type node2_idx) : g_(graph), edge_uidx(edge_idx), node1_uidx(node1_idx), node2_uidx(node2_idx) {};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
      return edge_i2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return vectorEdges[edge_i2u_[i]];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) {
      if (!neighbours.count(a.uindex()))
          return false;
      std::vector<size_type> neighbours_a = neighbours[a.uindex()];
      return std::find(neighbours_a.begin(), neighbours_a.end(), b.uindex()) != neighbours_a.end();
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
  Edge add_edge(const Node& a, const Node& b, edge_value_type value = edge_value_type()) {
    size_type a_uindex = a.uindex();
    size_type b_uindex = b.uindex();
    edge_type edge = Edge(this, size_all_edge, a_uindex, b_uindex);
    if (!has_edge(a, b)) {
        neighbours[a_uindex].push_back(b_uindex);
        neighbours[b_uindex].push_back(a_uindex);
        vectorEdges.push_back(edge);
        edges_.push_back(edgeinfo{value, num_edges()});
        edge_i2u_.push_back(size_all_edge);
        ++size_all_edge;
    }
      return edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      vectorNodes.clear();
      nodes_.clear();
      node_i2u_.clear();
      vectorEdges.clear();
      neighbours.clear();
      edges_.clear();
      edge_i2u_.clear();
      size_all_node = 0;
      size_all_edge = 0;
  }

  //
  // Node Iterator
  //
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator{
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    Node operator*() const {
        return g_->vectorNodes[g_->node_i2u_[idx_]];
    }
    NodeIterator& operator++() {
        ++idx_;
        return *this;
    }
    bool operator==(const NodeIterator& node_it) const {
        return (g_ == (node_it.g_)) && (idx_ == node_it.idx_);
    }
    bool operator!=(const NodeIterator& node_it) const {
        return !(*this == node_it);
    }

   private:
    friend class Graph;
    const graph_type* g_;
    size_type idx_;
    NodeIterator(const graph_type* g_, size_type idx) : g_(g_), idx_(idx) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const {
      return NodeIterator(this, 0);
  }
  node_iterator node_end() const {
      return NodeIterator(this, num_nodes());
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

    Edge operator*() const {
        node_type node2 = g_->vectorNodes[neighbours_[idx_]];
        edge_type edge = Edge(g_, 0, node_.uindex(), node2.uindex());
        return edge;
    }
    IncidentIterator& operator++() {
        ++idx_;
        return *this;
    }
    bool operator==(const IncidentIterator& incident_it) const {
        return (node_ == incident_it.node_ && idx_ == incident_it.idx_);
    }
    bool operator!=(const IncidentIterator& incident_it) const {
          return !(*this == incident_it);
    }
     
   private:
    friend class Graph;
    graph_type* g_;
    node_type node_;
    size_type idx_;
    std::vector<size_type> neighbours_;
    IncidentIterator(graph_type* graph, node_type node, size_type idx): g_(graph), node_(node), idx_(idx) {
        neighbours_ = g_->neighbours[node_.uindex()];
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
    EdgeIterator() {
    }

    EdgeIterator& operator++() {
        ++idx_;
        return *this;
    }
    Edge operator*() const {
        return g_->vectorEdges[g_->edge_i2u_[idx_]];
    }
    bool operator==(const EdgeIterator& edge_it) const {
        return idx_ == edge_it.idx_;
    }
    bool operator!=(const EdgeIterator& edge_it) const {
        return !(*this == edge_it);
    }
   private:
    friend class Graph;
    const graph_type* g_;
    size_type idx_;
    EdgeIterator(const graph_type* g_, size_type idx): g_(g_), idx_(idx) {}
  };

    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0);
    }
    edge_iterator edge_end() const {
        return EdgeIterator(this, num_edges());
    }
 // removal related code
public:
    /**
     * @brief removes node from a graph, and its incident edges
     *
     * @param Node n
     * @return 1 if node was removed, 0 if not
     * @pre  n of type Node
     * @post result_graph has not node n, nothing if graph didn't have node n
     *
     * Complexity: O(num_nodes + num_edges).
     */
    size_type remove_node(const Node& node) {
        if(has_node(node)){
            size_type node_uidx = node.uindex();
            size_type node_end_uidx = node_i2u_[size() - 1];
            size_type node_idx = node.index();
            node_i2u_[size() - 1] = node_uidx;
            node_i2u_[node_idx] = node_end_uidx;
            nodes_[node_uidx].idx_ = (size_type) -1;
            nodes_[node_end_uidx].idx_ = node_idx;
            auto result = node_i2u_.erase(node_i2u_.end() - 1);
            if (result == node_i2u_.end()) {
                std::vector<node_type> neighbours = getNeighbours(node);
                for (node_type neighbour : neighbours) {
                     remove_edge(neighbour, node);
                }
                return 1;
            } else {
                return 0;
            }
        }
        return 0;
    }
    /**
     * @brief removes node from a graph, and its incident edges
     *
     * @param Node n
     * @return result_graph.node_begin()
     * @pre  n_it of type NodeIterator
     * @post result_graph has not node n, nothing if graph didn't have node n
     *
     * Complexity: O(num_nodes + num_edges).
     */
    node_iterator remove_node(node_iterator n_it) {
        remove_node(*n_it);
        return node_begin();
    }
    /**
     * @brief removes edge from graph
     *
     * @param Node node1, Node node2
     * @return 1 if edge was removed, 0 if not
     * @pre  node1 and node2 of type Node
     * @post result_graph has not edge, nothing if graph didn't have edge e
     *
     * Complexity: O(num_edges).
     */
    size_type remove_edge(const Node& node1, const Node& node2){
        if (has_edge(node1, node2)) {
            for (auto edge : vectorEdges)  {
                if ( ((edge.node1() == node1) && (edge.node2() == node2))
                    || ((edge.node2() == node1) && (edge.node1() == node2)) ) {
                    return remove_edge(edge);
                }
            }
        }
        return 0;
    }
    /**
     * @brief removes edge from graph
     *
     * @param Edge e
     * @return 1 if edge was removed, 0 if not
     * @pre  e of type Edge
     * @post result_graph has not edge, nothing if graph didn't have edge e
     *
     * Complexity: O(num_edges).
     */
    size_type  remove_edge(const Edge& edge) {
        if(has_edge(edge.node1(), edge.node2())){
            size_type node_1_uidx = edge.node1().uindex();
            size_type node_2_uidx = edge.node2().uindex();
            size_type edge_uidx = find_edge_uidx(edge.node1(), edge.node2());
            size_type edge_end_uidx = edge_i2u_[num_edges() - 1];
            size_type edge_idx = edge.index();
            edge_i2u_[num_edges() - 1] = edge_uidx;
            edge_i2u_[edge_idx] = edge_end_uidx;
            edges_[edge_uidx].idx_ = (size_type) -1;
            edges_[edge_end_uidx].idx_ = edge_idx;
            auto result = edge_i2u_.erase(edge_i2u_.end() - 1);
            remove_el_in_list(neighbours[node_1_uidx], node_2_uidx);
            remove_el_in_list(neighbours[node_2_uidx], node_1_uidx);
            if (neighbours.count(node_1_uidx) && neighbours[node_1_uidx].size() == 0) {
                neighbours.erase(node_1_uidx);
            }
            if (neighbours.count(node_2_uidx) && neighbours[node_2_uidx].size() == 0) {
                neighbours.erase(node_2_uidx);
            }
            if (result == edge_i2u_.end()) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return 1;
        }
        return 0;
    }
    /**
     * @brief removes edge from graph
     *
     * @param EdgeIterator e_it
     * @return result_graph.edge_begin()
     * @pre  e_it of type Edgeterator
     * @post result_graph has not edge e anymore, nothing if graph didn't have edge e
     *
     * Complexity: O(num_edges).
     */
    edge_iterator  remove_edge(edge_iterator e_it){
        remove_edge(*e_it);
        return edge_begin();
    }
    /**
     * @brief removes element from a std::vector
     *
     * @param std::vector<size_type> vec
     * @param size_type el
     * @pre  el is an element of vec
     * @post result_vec has not el anymore, nothing if vec didn't have el
     *
     * Complexity: O(vec.size()).
     */
    void remove_el_in_list(std::vector<size_type>& vec, size_type el) {
        std::remove(vec.begin(), vec.end(), el);
        vec.pop_back();
    }
    /**
     * @brief finds the unique id of an edge
     *
     * @param Node a
     * @param Node b
     * @post returns unique id of edge, or size_all_edge if edge doesn't exist
     *
     * Complexity: O(num_edges).
     */
    size_type find_edge_uidx(const node_type& a, const node_type& b) {
        for(size_type i = 0; i < size_all_edge; ++i) {
            edge_type edge = vectorEdges[i];
            if ( (edge.node1() == a && edge.node2() == b) || (edge.node1() == b && edge.node2() == a) ){
                return i;
            }
        }
        return size_all_edge;
    }
 private:
    /*
     Node attributes
     */
    struct nodeinfo {
        Point p_;
        node_value_type v_;
        size_type idx_;
    };
    std::vector<nodeinfo> nodes_;
    std::vector<size_type> node_i2u_;
    std::vector<Node> vectorNodes;
    size_type size_all_node;
    /*
     Edge attributes
     */
    struct edgeinfo {
        edge_value_type v_;
        size_type idx_;
    };
    std::vector<edge_type> vectorEdges;
    std::vector<edgeinfo> edges_;
    std::vector<size_type> edge_i2u_;
    std::map<size_type, std::vector<size_type>> neighbours;
    size_type size_all_edge;
};

#endif // CME212_GRAPH_HPP
