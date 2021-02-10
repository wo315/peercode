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
template <typename V>
class Graph {
 private:

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** Synonym for node value type*/
  using node_value_type = V;
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
      size_ = 0;
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
        return graph_->vectorPoints[node_id_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
      node_value_type& value() {
          return value_;
      }
      const node_value_type& value() const {
          return value_;
      }
      size_type degree() const {
          if (graph_->edges.count(*this)) {
              return graph_->edges.at(*this).size();
          }
          else {
              return 0;
          }
      }
      incident_iterator edge_begin() const {
          return IncidentIterator(*this, graph_->edges.at(*this), 0);
      }
      incident_iterator edge_end() const {
          return IncidentIterator(*this, graph_->edges.at(*this), graph_->edges.at(*this).size());
      }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_ and node_id_ == n.node_id_);
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
      if (graph_ == n.graph_ and node_id_ < n.node_id_)
          return true;
      return false;
    }
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type node_id_;
    node_value_type value_;
    Node(const graph_type* graph, size_type node_id, node_value_type value)
      : graph_(const_cast<graph_type*>(graph)), node_id_(node_id), value_(value) {}
  };

 /**
    modify a node's value
  */
  void modify_node_value(size_type node_idx, node_value_type val) {
      assert(node_idx < size());
      vectorNodes[node_idx].value() = val;
  }
  /** Return the neighbors of a node
   */
    std::vector<Node> getNeighbours(size_type node_idx) {
        return edges[vectorNodes[node_idx]];
    }
    /** Return the neighbors of a node
     */
    V getValue(size_type node_idx) {
        return vectorNodes[node_idx].value();
    }
  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_;
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
    Node node = Node(this, size_, value);
    vectorNodes.push_back(node);
    vectorPoints.push_back(position);
    ++size_;
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.node_id_ < size())
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
    return vectorNodes[i];
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
      return node_1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return node_2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1() == e.node1() and node2() == e.node2()) or (node2() == e.node1() and node1() == e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (node1() < e.node1() and node2() < e.node2())
          return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
      Node node_1;
      Node node_2;
      Edge(Node node_1, Node node_2) : node_1(node_1), node_2(node_2) {};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
      return vectorEdges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    Edge edge = vectorEdges[i];
    return edge;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if (edges.count(a))
      return (std::find(edges.at(a).begin(), edges.at(a).end(), b) != edges.at(a).end());
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
    Edge edge = Edge(a, b);
    if (!has_edge(a, b)) {
     if (std::find(edges[a].begin(), edges[a].end(), b) == edges[a].end()) {
         edges[a].push_back(b);
     }
     if (std::find(edges[b].begin(), edges[b].end(), a) == edges[b].end()) {
         edges[b].push_back(a);
     }
     vectorEdges.push_back(edge);
    }
    return edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      vectorPoints.clear();
      vectorNodes.clear();
      vectorEdges.clear();
      edges.clear();
      size_ = 0;
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
        assert(idx_ <= (*node_ptr_).size());
        return (*node_ptr_)[idx_];
    }
    NodeIterator& operator++ () {
        if (idx_ < (*node_ptr_).size())
            ++idx_;
        return *this;
    }
    bool operator==(const NodeIterator& node_it) const {
        return (node_ptr_ == node_it.node_ptr_) and (idx_ == node_it.idx_);
    }
    bool operator!=(const NodeIterator& node_it) const {
        return !(*this == node_it);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    std::vector<Node>* node_ptr_;
    size_type idx_;
    NodeIterator(const std::vector<Node>* node_ptr, size_type idx) : node_ptr_(const_cast<std::vector<Node>*>(node_ptr)), idx_(idx) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const {
      return NodeIterator(&vectorNodes, 0);
  }
  node_iterator node_end() const {
      return NodeIterator(&vectorNodes, size());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    Edge operator*() const {
        assert(edge_idx_ < connections_.size());
        return Edge(node_, connections_[edge_idx_]);
    }
    IncidentIterator& operator++() {
        if (edge_idx_ < connections_.size())
            ++edge_idx_;
        return *this;
    }
    bool operator==(const IncidentIterator& incident_it) const {
        return (node_ == incident_it.node_ and edge_idx_ == incident_it.edge_idx_);
    }
    bool operator!=(const IncidentIterator& incident_it) const {
          return !(*this == incident_it);
    }
     
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    node_type node_;
    std::vector<node_type> connections_;
    size_type edge_idx_;
    IncidentIterator(node_type node, std::vector<node_type> connections, size_type edge_idx): node_(node), connections_(connections), edge_idx_(edge_idx) {}
    
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    EdgeIterator& operator++() {
        if (idx_ < (*edges_ptr_).size())
            ++idx_;
        return *this;
    }
    Edge operator*() const {
        assert(idx_ < (*edges_ptr_).size());
        return (*edges_ptr_)[idx_];
    }
    bool operator==(const EdgeIterator& edge_it) const {
        return idx_ == edge_it.idx_;
    }
    bool operator!=(const EdgeIterator& edge_it) const {
        return !(*this == edge_it);
    }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const std::vector<Edge>* edges_ptr_;
    size_type idx_;
    EdgeIterator(const std::vector<Edge>* edges_ptr, size_type idx): edges_ptr_(edges_ptr), idx_(idx) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
    edge_iterator edge_begin() const {
        return EdgeIterator(&vectorEdges, 0);
    }
    edge_iterator edge_end() const {
        return EdgeIterator(&vectorEdges, vectorEdges.size());
    }
 private:
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    /*
     Node attributes
     */
    std::vector<Node> vectorNodes;
    std::vector<Point> vectorPoints;
    size_type size_;
    /*
     Edge attributes
     */
    std::map<Node, std::vector<Node>> edges;
    std::vector<Edge> vectorEdges;
};

#endif // CME212_GRAPH_HPP
