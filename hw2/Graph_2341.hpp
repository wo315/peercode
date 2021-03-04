#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */


#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_set>

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

  struct node_;
  struct edge_;

  std::vector<node_> n_Vec; // contains all nodes
  std::vector<unsigned int> n_i2u; // contains active nodes

  std::vector<edge_> e_Vec; // contains all edges
  std::vector<unsigned int> e_i2u; // contains active nodes

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
      idx = 0;
      g = nullptr;
    }

    /** Update the value associated with the node. */
    void updateValue(node_value_type nv) {
      unsigned int node_uid = g->n_i2u[idx];
      g->n_Vec[node_uid].v = nv;
      return;
    }

    /** Returns reference to value associated with the node. */
    node_value_type& value() {
      unsigned int node_uid = g->n_i2u[idx];
      return g->n_Vec[node_uid].v;
    }

    /** Returns value associated with the node. */
    const node_value_type& value() const {
      unsigned int node_uid = g->n_i2u[idx];
      return g->n_Vec[node_uid].v;
    }

    /** Return this node's position. */
    const Point& position() const {
      unsigned int node_uid = g->n_i2u[idx];
      return g->n_Vec[node_uid].pos;
    }

    /** Returns reference to node's position. */
    Point& position() {
      unsigned int node_uid = g->n_i2u[idx];
      return g->n_Vec[node_uid].pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx;
    }

    /** Update's the node's index to @a i. */
    void updateIndex(unsigned int i) {
      idx = i;
      return;
    }

    /** Returns the number of nodes that this node is connected to. */
    size_type degree() const {
      unsigned int node_uid = g->n_i2u[idx];
      return g->n_Vec[node_uid].incident_edges.size();
    }

    /** Returns a beginning iterator that iterates through edges that are
        connected to the node. */
    IncidentIterator edge_begin() const {
      unsigned int node_uid = g->n_i2u[idx];
      return IncidentIterator(this, g->n_Vec[node_uid].incident_edges.begin());
    }

    /** Returns an ending iterator that iterates through edges that are
        connected to the node. */
    IncidentIterator edge_end() const {
      unsigned int node_uid = g->n_i2u[idx];
      return IncidentIterator(this, g->n_Vec[node_uid].incident_edges.end());
    }

    /** Removes the edge pointed to by @a i_it. Any existing IncidentIterator's
        or EdgeIterators will be made invalid. */
    IncidentIterator remove_edge(IncidentIterator i_it) const {
      g->remove_edge(*i_it);
      return this->edge_begin();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (this->g == n.g & this->g->n_i2u[idx] == n.g->n_i2u[n.index()]) {
	return true;
      }
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
      if (*this == n) {
        return false;
      } else if (g != n.g) {
        return (&g < &n.g);
      } else if (n.idx < this->idx) {
        return true;
      }
      return false;
    }

    /** Returns the number of nodes in the graph. */
    size_type nodes_in_graph() {
      return g->size();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    unsigned int idx;
    Graph *g;

    // Constructor
    Node(Graph *grph, unsigned int i) {g = grph; idx = i;}

  };

  /** Return the number of nodes in the graph.
   * Complexity: O(1).
   */
  size_type size() const {
    return n_i2u.size();
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
  Node add_node(const Point& p, const node_value_type& v = node_value_type()) {
    unsigned int idx = n_i2u.size();
    Node n = Node(this, idx);
    unsigned int uid = n_Vec.size();

    n_i2u.push_back(uid);
    n_Vec.push_back(node_(n, p, v));
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if(n.g == this && n_Vec[n_i2u[n.idx]].n == n) {
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
    return n_Vec[n_i2u[i]].n;
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
      unsigned int edge_uid = g->e_i2u[idx];
      return g->e_Vec[edge_uid].n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      unsigned int edge_uid = g->e_i2u[idx];
      return g->e_Vec[edge_uid].n2;
    }

    /** Return the index of this Edge */
    size_type index() const {
      return idx;
    }

    /** Update the index of this Edge to @a i */
    void updateIndex(unsigned int i) {
      idx = i;
      return;
    }

    /** Return the value of this Edge */
    edge_value_type& value(){
      unsigned int edge_uid = g->e_i2u[idx];
      return g->e_Vec[edge_uid].v;
    }

    /** Return the value of this Edge */
    const edge_value_type& value() const {
      unsigned int edge_uid = g->e_i2u[idx];
      return g->e_Vec[edge_uid].v;
    }

    /** Update the value of this Edge to @a ev */
    void updateValue(edge_value_type ev) {
      unsigned int edge_uid = g->e_i2u[idx];
      g->e_Vec[edge_uid].v = ev;
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
      } else if (g != e.g) {
        return (&g < &e.g);
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

    // Constructor
    Edge(Graph *grph, unsigned int i) {g = grph; idx = i;}

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return e_i2u.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return e_Vec[e_i2u[i]].e;
  }

  /** Test whether two nodes are connected by an edge.
   * If edge exists, the edge's index is returned.
   * Otherwise -1 is returned.
   *
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  int has_edge_idx(const Node& a, const Node& b) const {

    if(a.degree() <= b.degree()) {
      for(auto it = a.edge_begin(); it != a.edge_end(); ++it) {
        if ((*it).node2() == b) { return (*it).index(); }
      }
    } else {
      for(auto it = b.edge_begin(); it != b.edge_end(); ++it) {
        if ((*it).node2() == a) { return (*it).index(); }
      }
    }
    return -1;

  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if(has_edge_idx(a, b) == -1) { return false; } else { return true; }
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

    int idx = has_edge_idx(a, b);

    if (idx == -1) {
      // edge does not exist
      idx = e_i2u.size();
      Edge e(this, idx);
      int uid = e_Vec.size();

      e_i2u.push_back(uid);
      e_Vec.push_back(edge_(e, a, b));

      n_Vec[n_i2u[a.index()]].incident_edges.insert(uid);
      n_Vec[n_i2u[b.index()]].incident_edges.insert(uid);
      return e;

    } else {
      // edge does exist
      Edge e = e_Vec[e_i2u[idx]].e;
      if(e.node1() == b) {
        e_Vec[e_i2u[idx]].switchNodes();
      }
      return e;
    }

  }


  /** Remove an Edge from the graph, if an Edge connecting @a n1 & @a n2 exists.
   * @pre @a n1 and @a n2 are distinct valid nodes of this graph.
   * @return 1 if an Edge is removed, 0 otherwise
   * @post If old has_edge(@a n1, @a n2), new num_edges() == old num_edges() - 1.
   *       Else,                          new num_edges() == old num_edges().
   *
   * Will invalidate any existing EdgeIterator's or IncidentIterator's.
   *
   * Complexity: No more than O(num_nodes()).
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    int eidx = this->has_edge_idx(n1, n2);
    if (eidx != -1) { // edge exists

      // update incident edge set for n1 & n2
      n_Vec[n_i2u[n1.index()]].incident_edges.erase(e_i2u[eidx]);
      n_Vec[n_i2u[n2.index()]].incident_edges.erase(e_i2u[eidx]);

      // remove edge from e_i2u (swap and pop)
      unsigned int uid = e_i2u[e_i2u.size()-1]; // last element in vector
      e_i2u[eidx] = uid;
      e_i2u.pop_back(); // remove last element
      e_Vec[uid].e.updateIndex(eidx); // update index for swapped element

      return 1;
    } else {
      return 0;
    }
  }

  /** Remove an Edge from the graph, if Edge @a e exists in this graph.
   * @pre @a e is a valid Edge.
   * @return 1 if an Edge is removed, 0 otherwise
   * @post If Edge is in graph, new num_edges() == old num_edges() - 1.
   *       Else,                new num_edges() == old num_edges().
   *
   * Will invalidate any existing EdgeIterator's or IncidentIterator's.
   *
   * Complexity: No more than O(num_nodes()).
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove a Node from the graph, if Node @a n exists in this graph.
   * @pre @a n is a valid Node.
   * @return 1 if Node is removed, 0 otherwise
   * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1.
   *       Else,                  new num_nodes() == old num_nodes().
   *
   * Will invalidate any existing NodeIterator's, EdgeIterators
   * or IncidentIterator's.
   *
   * Complexity: No more than O(num_nodes()).
   */
  size_type remove_node (const Node& n) {
    if(this->has_node(n)){ // node exists
      // remove all connected edges
      for(auto it = n.edge_begin(); it != n.edge_end(); ) {
        it = n.remove_edge(it);
      }

      // remove node from n_i2u (swap and pop)
      unsigned int nidx = n.index();
      unsigned int uid = n_i2u[n_i2u.size()-1]; // last element in vector
      n_i2u[nidx] = uid;
      n_i2u.pop_back(); // remove last element
      n_Vec[uid].n.updateIndex(nidx); // update index for swapped element

      return 1;
    } else {
      return 0;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    n_Vec.clear();
    n_i2u.clear();
    e_Vec.clear();
    e_i2u.clear();
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

    /** Returns Node pointed to by NodeIterator */
    Node operator*() const { return *n; };

    /** Increments NodeIterator to point to next valid Node in graph */
    NodeIterator& operator++() {
        if(n->index()+1 == n->g->num_nodes()) {
	    // current iterator is pointing to last node
            n = nullptr;
        } else {
            unsigned int node_uid = n->g->n_i2u[n->index()+1];
            n = &n->g->n_Vec[node_uid].n;
        }
        return *this;
    }

    /** Test whether this NodeIterator and @a x are equal. */
    bool operator==(const NodeIterator& x) const {
        return (n == x.n);
    }

    /** Test whether this NodeIterator and @a x are not equal. */
    bool operator!=(const NodeIterator& x) const {
        return !(n == x.n);
    }

   private:
    friend class Graph;

    /** Construct a valid NodeIterator. */
    NodeIterator(const Node* n) : n{n} {}

  };

  /** Returns a NodeIterator pointing to the first Node in the Graph,
   *  i.e. the Node with an index() == 0
   */
  NodeIterator node_begin() const {
    if(this->num_nodes() == 0) {
      return NodeIterator(nullptr);
    } else {
      return NodeIterator(&n_Vec[n_i2u[0]].n);
    }
  }

  /** Returns a NodeIterator pointing to the past-the-end Node in the Graph.
   *  It does not point to any element, and thus can not be dereferenced.
   */
  NodeIterator node_end() const {
    return NodeIterator(nullptr);
  }


  /** Remove a Node from the graph, if Node pointed to by @a n_it
   * exists in this graph.
   * @pre @a n_it is a valid NodeIterator.
   * @return NodeIterator pointing to valid Node in graph
   * @post If old has_node(@a n), new num_nodes() == old num_nodes() - 1.
   *       Else,                  new num_nodes() == old num_nodes().
   *
   * Will invalidate any existing NodeIterator's, EdgeIterators
   * or IncidentIterator's.
   *
   * Complexity: No more than O(num_nodes()).
   */
  //--functionality_1
  //--node_begin might already been removed
  //--START
  NodeIterator remove_node(NodeIterator n_it) {
    remove_node(*n_it);
    return this->node_begin();
  }
  //--END

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
    std::unordered_set<unsigned int>::iterator iter;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Returns Edge pointed to by IncidentIterator */
    Edge operator*() const {
      unsigned int edge_uid = *iter;
      Edge e = n->g->e_Vec[edge_uid].e;
      if (e.node1() != *n) {
        n->g->e_Vec[edge_uid].switchNodes();
      }
      return e;
    }

    /** Increments IncidentIterator to point to next edge that is
     *  incident to Node @a n.
     */
    IncidentIterator& operator++() {
      ++iter;
      return *this;
    }

    /** Test whether this IncidentIterator and @a x are equal. */
    bool operator==(const IncidentIterator& x) const {
        return (n == x.n && iter == x.iter);
    }

    /** Test whether this IncidentIterator and @a x are not equal. */
    bool operator!=(const IncidentIterator& x) const {
        return !(n == x.n && iter == x.iter);
    }

   private:
    friend class Graph;

    /** Construct a valid IncidentIterator. */
    IncidentIterator(const Node* n_,
      std::unordered_set<unsigned int>::iterator iter_) : n(n_), iter(iter_) { }

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

    /** Returns Edge pointed to by EdgeIterator */
    Edge operator*() const {
      return *e;
    }

    /** Increments EdgeIterator to point to next edge in Graph. */
    EdgeIterator& operator++() {
        if(e->index()+1 == e->g->num_edges()) {
	    // current iterator is pointing to last node
            e = nullptr;
        } else {
            unsigned int edge_uid = e->g->e_i2u[e->index()+1];
            e = &e->g->e_Vec[edge_uid].e;
        }
       return *this;

    }

    /** Test whether this EdgeIterator and @a x are equal. */
    bool operator==(const EdgeIterator& x) const {
      return (e == x.e);

    }

    /** Test whether this EdgeIterator and @a x are not equal. */
    bool operator!=(const EdgeIterator& x) const {
      return !(e == x.e);

    }

   private:
    friend class Graph;

    /** Construct a valid EdgeIterator. */
    EdgeIterator(const Edge* e) : e{e} {}

  };

  /** Returns an EdgeIterator pointing to the first Edge in the Graph,
   *  i.e. the Edge with an index() == 0
   */
  EdgeIterator edge_begin() const {
    if(this->num_edges() == 0){
      return EdgeIterator(nullptr);
    } else {
      return EdgeIterator(&e_Vec[e_i2u[0]].e);
    }
  }

  /** Returns an EdgeIterator pointing to the past-the-end Edge in the Graph.
   *  It does not point to any element, and thus can not be dereferenced.
   */
  EdgeIterator edge_end() const {
    return EdgeIterator(nullptr);
  }

  /** Remove an Edge from the graph, if Edge pointed to by
   * @a e_it exists in this graph.
   *
   * @pre @a e_it is a valid EdgeIterator.
   * @return 1 if an Edge is removed, 0 otherwise
   * @post If Edge is in graph, new num_edges() == old num_edges() - 1.
   *       Else,                new num_edges() == old num_edges().
   *
   * Will invalidate any existing EdgeIterator's or IncidentIterator's.
   *
   * Complexity: No more than O(num_nodes()).
   */
  EdgeIterator remove_edge(EdgeIterator e_it) {
    remove_edge(*e_it);
    return this->edge_begin();
  }


 private:

  struct node_ {
    Node n;
    Point pos;
    node_value_type v = node_value_type();
    std::unordered_set<unsigned int> incident_edges; // contains the edge's uid

    node_(const Node nt, const Point p, node_value_type v) :
    	 n(nt), pos(p), v(v) {}
    bool operator==(const node_& x) {
        if(n == x.n) {
	  return true;
        } else { return false; }
    }

  };

  struct edge_ {
    Edge e;
    Node n1;
    Node n2;
    edge_value_type v = edge_value_type();

    void switchNodes() {
      Node temp = n2;
      n2 = n1;
      n1 = temp;
    }
    edge_(const Edge et, const Node n_1, const Node n_2)
			 : e(et), n1(n_1), n2(n_2) {}
  };


};

#endif // CME212_GRAPH_HPP
