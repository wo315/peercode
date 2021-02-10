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

template <typename V>
class Graph {

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  using node_value_type = V;
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
    // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
      idx = 0;
      g = nullptr;
    }

    void updateValue(node_value_type nv) {
      g->n_Vec[idx].v = nv;
      return;
    }

    node_value_type& value() {
      return g->n_Vec[idx].v;
    }

    const node_value_type& value() const {
      return g->n_Vec[idx].v;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g->n_Vec[idx].pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx;
    }

    size_type degree() const {
        return g->n_Vec[idx].incident_edges.size();
    }

    IncidentIterator edge_begin() const {
        return IncidentIterator(this);
    }

    IncidentIterator edge_end() const {
        return IncidentIterator(nullptr);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->g == n.g && this->idx == n.index()) {
	return true;
      }
      (void) n;          // Quiet compiler warning
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
      // HW0: YOUR CODE HERE
      if (*this == n) {
        return false;
      } else if (n.idx < this->idx) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    unsigned int idx;
    Graph *g;
    Node(Graph *grph, unsigned int i) {g = grph; idx = i;}

  };

  /** Return the number of nodes in the graph.
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return n_Vec.size();
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

/*  Node add_node(const Point& position) {
    size_type n_nodes = this->num_nodes();
    Node n = Node(this, n_nodes);
    n_Vec.push_back(node_(n, position, n_nodes));
    return n;
  }
*/
  Node add_node(const Point& p, const node_value_type& v = node_value_type()) {
    size_type n_nodes = this->num_nodes();
    Node n = Node(this, n_nodes);
    n_Vec.push_back(node_(n, p, n_nodes, v));
    return n;
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE

    if(n.g == this && n_Vec[n.idx].n == n) {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return n_Vec[i].n;
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
      g = nullptr;
      idx = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return g->e_Vec[idx].n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return g->e_Vec[idx].n2;
    }

    size_type index() const {
        return idx;
    }

    /** Test whether this edge and @a e are equal.
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->node1() == e.node1() && this->node2() == e.node2()) {
	return true;
      } else if (this->node1() == e.node2() && this->node2() == e.node1()) {
	return true;
      } else {
	return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (*this == e) {
         return false;
      } else if (this->idx < e.idx) {
         return true;
      }

      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph *g;
    unsigned int idx;
    Edge(Graph *grph, unsigned int i) {g = grph; idx = i;}

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return e_Vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return e_Vec[i].e;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    unsigned int i = 0;
    while (i < this-> num_edges()) {
      if (e_Vec[i].e.node1() == a && e_Vec[i].e.node2() == b) {
        return true;
      } else if (e_Vec[i].e.node1() == b && e_Vec[i].e.node2() == a) {
        return true;
      } else { i++; }
    }

    return false;
  }


  struct node_ {
    const Node n;
    const Point pos;
    unsigned int idx;
    node_value_type v = node_value_type();
    std::vector<unsigned int> incident_edges;
    node_(const Node nt, const Point p, unsigned int i, node_value_type v) :
    	 n(nt), pos(p), idx(i), v(v) {}
    bool operator==(const node_& x) {
        if(n == x.n) { return true; } else { return false; }
    }

  };

  struct edge_ {
    const Edge e;
    Node n1;
    Node n2;
    unsigned int idx;
    void switchNodes() {
      Node temp = n2;
      n2 = n1;
      n1 = temp;
    }
    edge_(const Edge et, const Node n_1, const Node n_2,
		 unsigned int i) : e(et), n1(n_1), n2(n_2), idx(i) {}
  };

  std::vector<node_> n_Vec;
  std::vector<edge_> e_Vec;


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

    unsigned int n_edges = this->num_edges();
    if (!this->has_edge(a, b)) {
      Edge e(this, n_edges);
      e_Vec.push_back(edge_(e, a, b, n_edges));
      n_Vec[a.index()].incident_edges.push_back(n_edges);
      n_Vec[b.index()].incident_edges.push_back(n_edges);
      return e;
    } else {
      unsigned int i = 0;
      bool need_to_flip = false;
      while (i < this->num_edges()) {
        if (e_Vec[i].e.node1() == a && e_Vec[i].e.node2() == b) {
          break;
        } else if (e_Vec[i].e.node1() == b && e_Vec[i].e.node2() == a) {
          need_to_flip = true;
          break;
        } else { i++; }
      }
      if (need_to_flip) {e_Vec[i].e.node1() = a; e_Vec[i].e.node2() = b;}
      return e_Vec[i].e;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    n_Vec.clear();
    e_Vec.clear();
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


    // HW1 #2: YOUR4- CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    const Node* n; // pointer to current node

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
        n = nullptr;
    }

    // Node operator*() const
    Node operator*() const { return *n; };

    // NodeIterator& operator++()
    NodeIterator& operator++() {
        if(n->index()+1 == n->g->num_nodes()) {
	    // current iterator is pointing to last node
            n = nullptr;
        } else {
            n = &n->g->n_Vec[n->index()+1].n;
        }
        return *this;
    }

    bool operator==(const NodeIterator& x) const {
        return (n == x.n);
    }

    bool operator!=(const NodeIterator& x) const {
        return !(n == x.n);
    }

    NodeIterator(const Node* n) : n{n} {}

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

  };

  // HW1 #2: YOUR CODE HERE

  // Supply definitions AND SPECIFICATIONS for:
  NodeIterator node_begin() const {
    return NodeIterator(&n_Vec[0].n);
  //  return NodeIterator(nullptr);

  }

  NodeIterator node_end() const {
    return NodeIterator(nullptr);
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

    const Node* n;
    unsigned int cur_edge_idx;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      n = nullptr;
      cur_edge_idx = 0;
    }


    Edge operator*() const {
      int edge_idx = n->g->n_Vec[n->index()].incident_edges[cur_edge_idx];
      Edge e = n->g->e_Vec[edge_idx].e;
      if (e.node1() != *n) {
        n->g->e_Vec[edge_idx].switchNodes();
      }
      return e;
    }

    IncidentIterator& operator++() {
      if(cur_edge_idx+1 == n->degree()) {
	n = nullptr;
        cur_edge_idx = 0;
      } else {
        cur_edge_idx++;
      }
      return *this;
    }

    bool operator==(const IncidentIterator& x) const {
        return (n == x.n && cur_edge_idx == x.cur_edge_idx);
    }

    bool operator!=(const IncidentIterator& x) const {
        return !(n == x.n && cur_edge_idx == x.cur_edge_idx);
    }

    IncidentIterator(const Node *n) : n{n}, cur_edge_idx{0} {}

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
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

    const Edge *e;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      e = nullptr;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return *e;
    }

    EdgeIterator& operator++() {
        if(e->index()+1 == e->g->num_edges()) {
	    // current iterator is pointing to last node
            e = nullptr;
        } else {
            e = &e->g->e_Vec[e->index()+1].e;
        }
       return *this;

    }

    bool operator==(const EdgeIterator& x) const {
      return (e == x.e);

    }

    bool operator!=(const EdgeIterator& x) const {
      return !(e == x.e);

    }

    EdgeIterator(const Edge* e) : e{e} {}

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  EdgeIterator edge_begin() const {
    return EdgeIterator(&e_Vec[0].e);

  }

  EdgeIterator edge_end() const {
    return EdgeIterator(nullptr);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.



};

#endif // CME212_GRAPH_HPP
