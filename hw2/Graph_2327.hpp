#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cmath>
#include <unordered_set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

// set typename to int by default so if a typename is not
// provided code still works
template <typename V, typename E>
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

  // TODO: typedef vs using??
  typedef E edge_value_type;

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
  class Node: private totally_ordered<Node> {
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
      graph = nullptr;
    }

    Point& position() const {
       //if (idx >= graph->size()) throw "Invalid node";
       return graph->nodes_[id].loc;
    }

    /** Return this node's position. */
    /*const Point& position() const {
      // make sure index is valid and node is still part of the graph
      if (idx >= graph->size()) throw "Invalid node";

      // get the position of the node with given index
      return graph->nodes_[idx].loc;
    }*/

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      //return idx;
      // get current index of the node
      return graph->nodes_[id].idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    /** 
     * @brief return reference to the value of the given node
     *
     * @return reference to node's value
     *
     * @pre Node is a valid node of the current graph
     */
    node_value_type& value() {
      //if (idx >= graph->size()) throw "Invalid node";
      if (!is_valid()) throw "Invalid node";

      return graph->nodes_[id].value;
    }
 
    /** 
     * @brief return constant reference to the value of the given node
     *
     * @return const reference to node's value
     *
     * @pre Node is a valid node of the current graph
     */
    const node_value_type& value() const {
      //if (idx >= graph->size()) throw "Invalid node";
      if (!is_valid()) throw "Invalid node";

      return graph->nodes_[id].value;
    }

    /**
     * @return degree of the given node
     *
     * @pre Node is a valid node of the current graph
     */
    size_type degree() const {
      if (!is_valid()) throw "Invalid node";

      return graph->incident_nodes[id].size();
    }

    /**
     * @return IncidentIterator pointing to the first edge incident to the node
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(this, graph->incident_nodes[id].begin());
    }

    /**
     * @return IncidentIterator indicating the end of the set of incident edges
     */
    incident_iterator edge_end() const {
      return IncidentIterator(this, graph->incident_nodes[id].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      (void) n;          // Quiet compiler warning

      if (n.graph == graph && n.id == id) {
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
      (void) n;           // Quiet compiler warning
     
      if (id < n.id) {
        return true;
      }

      if (graph < n.graph) {
        return true;
      }

      return false;
    }

    bool is_valid() const {
      return (id >= 0 && id < graph->nodes_.size() &&
              index() < graph->size() && graph->i2u[graph->nodes_[id].idx] == id);
   }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    // here 'idx' is a value of the nodes that can be changed. 
    // the indices are assumed to be consecutive, so that if a node is deleted, 
    // re-indexing occurs. 
    Graph* graph;
    size_type id;

    Node(const Graph* g, size_type uid)
        : graph(const_cast<Graph*>(g)), id(uid) {
    }
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return n_nodes;
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
    (void) position;

    node_info n;
    n.idx = i2u.size();
    n.loc = position;
    n.value = value;
    
    // should have that size(i2u) == n_nodes;
    nodes_.push_back(n);
    i2u.push_back(n_encountered_nodes);

    ++n_nodes;    
    ++n_encountered_nodes;

    return Node(this, n_encountered_nodes-1);
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    (void) n;            // Quiet compiler warning
    
    // if index is valid and node's graph pointer points
    // to this graph, graph has the node
    //if (n.index() < n_nodes && n.graph == this) {        
    if (n.graph == this && n.is_valid()) {
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
    (void) i;             // Quiet compiler warning

    if (i >= n_nodes) throw "Invalid node index";

    // get the fixed id of the node
    return Node(this, i2u[i]);

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
      graph = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph, idx1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph, idx2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
        
      if (e.graph == graph && e.idx1 == idx1 && e.idx2 == idx2) {
        return true;
      }
      if (e.graph == graph && e.idx2 == idx1 && e.idx1 == idx2) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      
      // order nodes within an edge and
      // use lexicographic ordering on the node indices

      size_type smaller_edge = idx1 < idx2 ? idx1 : idx2;
      size_type larger_edge = idx1 > idx2 ? idx1 : idx2;

      size_type e_smaller_edge = e.idx1 < e.idx2 ? e.idx1 : e.idx2;
      size_type e_larger_edge = e.idx1 > e.idx2 ? e.idx1 : e.idx2;

      if (smaller_edge < e_smaller_edge) {
        return true;
      }
      if (smaller_edge == e_smaller_edge && larger_edge < e_larger_edge) {
        return true;
      }
      
      if (graph < e.graph) {
        return true;
      }
      return false;
    }

    double length() const {
      return 0.05;
    }

    edge_value_type& value() {
      if (idx1 < idx2) {
        return graph->values[idx1][idx2];
      }
      return graph->values[idx2][idx1];
    }

    const edge_value_type& value() const {
      if (idx1 < idx2) {
        return graph->values[idx1][idx2];
      }
      return graph->values[idx2][idx1];
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    Graph* graph;

    size_type idx1;
    size_type idx2;

    Edge(const Graph* g, size_type uid1, size_type uid2) 
        : graph(const_cast<Graph*>(g)), idx1(uid1), idx2(uid2) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return n_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    (void) i;             // Quiet compiler warning
    
    if (i >= edges_.size()) throw "Invalid edge index";

    return Edge(this, edges_[i].first, edges_[i].second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    (void) a; (void) b;   // Quiet compiler warning

    // get the ids of nodes a and b
    size_type idx1 = a.id;
    size_type idx2 = b.id;

    // order the nodes
    size_type smaller_idx;
    size_type larger_idx;
    
    if (idx1 < idx2) {
      smaller_idx = idx1;
      larger_idx = idx2;
    }
    else {
      smaller_idx = idx2;
      larger_idx = idx1;
    }

    std::pair<size_type, size_type> ed(smaller_idx, larger_idx);

    // check to see if pair is present in `edges_set_`. 
    // O(log(num_edges())) runtime
    return edges_set_.find(ed) != edges_set_.end();

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

    // get fixed indices of nodes a and b
    size_type idx1 = a.id;
    size_type idx2 = b.id;

    if (has_edge(a, b)) {
        return Edge(this, idx1, idx2);
    }

    // order the ids and add the edge to both `edges_` and `edges_set_`
    // so that both containers are up-to-date
    size_type smaller_idx = idx1 < idx2 ? idx1 : idx2;
    size_type larger_idx = idx1 > idx2 ? idx1 : idx2;
    std::pair<size_type, size_type> ed(smaller_idx, larger_idx);
    ++n_edges;
    ++n_encountered_edges;
    
    edges_.push_back(ed);
    edges_set_.insert(ed);
    
    // update values
    if (values.find(smaller_idx) != values.end()) {
        values[smaller_idx].emplace(larger_idx, edge_value_type {});
    }
    else {
        std::unordered_map<size_type, edge_value_type> m;
        m.emplace(larger_idx, edge_value_type {});
        values.emplace(smaller_idx, m);
    }

    // update incident_nodes
    if (incident_nodes.find(idx1) != incident_nodes.end()) {
        // if Node a is already a key in incident_nodes,
        // add idx2 to the corresponding vector
        incident_nodes[idx1].insert(idx2);
    }
    
    else {
        // otherwise, initialize a new vector with idx2
        // and add new (key, value) pair to incident_nodes
        std::unordered_set<size_type> set1{idx2};
        std::pair<size_type, std::unordered_set<size_type>> pair1(idx1, set1);
        incident_nodes.insert(pair1);
        
    }
    
    // same for Node b
    if (incident_nodes.find(idx2) != incident_nodes.end()) {
        incident_nodes[idx2].insert(idx1);
    }
    else {
        std::unordered_set<size_type> set2{idx1};
        std::pair<size_type, std::unordered_set<size_type>> pair2(idx2, set2);
        incident_nodes.insert(pair2);
    }
    return Edge(this, idx1, idx2);
  }

  /*
   * @brief removes given node from the graph
   * 
   * @return 1 if node was removed successfully, 0 if given node is invalid
   * @param n node to remove
   * @post if node is removed successfully,
   *            new size() = old size() - 1
   *            _n_ is invlalidated
   *            new num_edges = old num_edges - n.degree()
   *       else
   *            graph is left unchanged
   */
  size_type remove_node(const Node& n) {
    
    if (!n.is_valid()) {
        return 0;
    }
    
    n_nodes--;

    std::vector<Edge> to_remove {};
    for (incident_iterator it = n.edge_begin(); it != n.edge_end(); ++it) {
        to_remove.push_back(*it);
    }
    for (Edge e : to_remove) {
        remove_edge(e);
    }

    // get index of node
    unsigned n_index = n.index();
    unsigned n_id = n.id;
    unsigned uid_change = i2u[i2u.size()-1];
    std::swap(i2u[n_index], i2u[i2u.size()-1]);

    // since the nodes swapped positions, the new index of the 
    // swapped node will be equal to the old nodes index
    nodes_[uid_change].idx = n_index;

    i2u.pop_back();
    incident_nodes.erase(n_id);

    return 1;
  }

  /*
   * @brief remove a node using NodeIterators
   *
   * @pre n_it is a valid node iterator
   * @param n_it iterator pointing to node to be removed
   * @return iterator pointing to another valid node in the graph
   * @post if n_it is valid,
   *          new size() = old size() - 1
   *          new num_edges = old num_edges - n.degree()
   *          n_it is invalidated
   *       else
   *          graph remains unchanged
   */
  node_iterator remove_node(node_iterator n_it) {
    Node n = *n_it;
    Graph* g = n.graph;

    if (!n.is_valid()) {
         return ++n_it;
    } 
     
    n_nodes--;
     
    std::vector<Edge> to_remove {};
    for (incident_iterator it = n.edge_begin(); it != n.edge_end(); ++it) {
        to_remove.push_back(*it);
    }   
    for (Edge e : to_remove) {
        remove_edge(e);
    }   
     
    // get index of node
    unsigned n_index = n.index();
    unsigned n_id = n.id;
    unsigned uid_change = i2u[i2u.size()-1];
    std::swap(i2u[n_index], i2u[i2u.size()-1]);
     
    // since the nodes swapped positions, the new index of the 
    // swapped node will be equal to the old nodes index
    nodes_[uid_change].idx = n_index;
      
    i2u.pop_back();
    incident_nodes.erase(n_id);

    return NodeIterator(g, n_index);

  }

  /*
   * @brief remove an edge from the graph
   *
   * @param a, b Nodes representing edge to be removed
   * @pre a, b are valid nodes of the graph and are connected by an edge
   * @post if there is an edge in the graph connecting a and b,
   *          new num_edges = old num_edges - 1
   *          Edge connecting a and b is invalidated
   *       else
   *          graph remains unchanged
   * @return 1 if edge was removed, 0 otherwise
   */
  size_type remove_edge(const Node& a, const Node& b) {

    size_type id1 = a.id;
    size_type id2 = b.id;
    
    size_type smaller_id = id1 < id2 ? id1 : id2;
    size_type larger_id = id1 > id2 ? id1 : id2;

    std::pair<size_type, size_type> ed(smaller_id, larger_id);

    auto it = std::find(edges_set_.begin(), edges_set_.end(), ed);
    if (it == edges_set_.end()) {
        // if 'it' isnt in the set, it's not a part of the graph
        return 0;
    }

    n_edges--;

    // O(1)
    incident_nodes[smaller_id].erase(larger_id);
    incident_nodes[larger_id].erase(smaller_id);

    // O(logn)
    edges_set_.erase(ed);

    // linear in number of edges
    auto it2 = std::find(edges_.begin(), edges_.end(), ed);
    std::iter_swap(it2, edges_.end());
    edges_.pop_back();

    values[smaller_id].erase(larger_id);

    return 1;

  }

  /*
   * @brief remove an edge from the graph
   *
   * @param edge to be removed
   * @pre _edge_ is a valid edge of the graph
   * @post if _edge_ is a valid edge of the graph
   *          new num_edges = old num_edges - 1
   *          _edge_ is invalidated
   *       else
   *          graph remains unchanged
   * @return 1 if edge was removed, 0 otherwise
   */
  size_type remove_edge(const Edge& edge) {
    Node a = edge.node1();
    Node b = edge.node2();
    return remove_edge(a, b);
  }

  /*
   * @brief remove an edge from the graph
   *
   * @param e_it edge iterator pointing to edge to be removed
   * @pre _e_it_ is a valid edge iterator
   * @post new num_edges = old num_edges - 1
   * @post _e_it_ and corresponding edge are invalidated
   */
  edge_iterator remove_edge(edge_iterator e_it) {
   
    unsigned index = e_it.index;
    Graph* g = e_it.graph;

    Edge e = *e_it;

    size_type id1 = e.node1();
    size_type id2 = e.node2();
     
    size_type smaller_id = id1 < id2 ? id1 : id2;
    size_type larger_id = id1 > id2 ? id1 : id2;

    std::pair<size_type, size_type> ed(smaller_id, larger_id);

    n_edges--;
 
    // O(1)
    incident_nodes[smaller_id].erase(larger_id);
    incident_nodes[larger_id].erase(smaller_id);
 
    // O(logn)
    edges_set_.erase(ed);
 
    // O(1)
    std::swap(edges_[index], edges_[edges_.size()-1]);
    edges_.pop_back();

    values[smaller_id].erase(larger_id);
    
    return EdgeIterator(g, index);

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    n_nodes = 0;
    n_edges = 0;
    
    nodes_.clear();
    edges_.clear();
    edges_set_.clear();
    incident_nodes.clear();
    values.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private  equality_comparable<NodeIterator> { 
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

    /**
     * @return the Node object being pointed to by iterator
     */
    Node operator*() const {
        // iter loops over elements of i2u
        return Node(graph, graph->i2u[index]);
    }

    /**
     * @brief increments the iterator to point to the next Node
     *
     * @post new Node being pointed to has index one larger than old Node
     */
    NodeIterator& operator++() {
        index++;
        return *this;
    }

    /**
     * @brief checks if this iterator is equal to a given iterator
     *
     * @param[in] the iterator to check equality with
     * @return true if _it_ equals the current iterator, false otherwise
     *
     * @pre _it_ is a valid NodeIterator
     */
    bool operator==(const NodeIterator& it) const {
        return (index == it.index && graph == it.graph);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph;
    size_type index;

    NodeIterator(const Graph* g, size_type i) : graph{const_cast<Graph*>(g)}, index{i} {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
   * @return NodeIterator pointing to the first Node in the graph
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /**
   * @return NodeIterator indicating the end of the collection of Nodes
   */
  node_iterator node_end() const {
    // end index == 1 + index of last node
    return NodeIterator(this, i2u.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  //class IncidentIterator : private totally_ordered<IncidentIterator> {
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    
    /**
     * @return an Edge object e with e.node1() == node being iterated over 
     *         and e.node2() == incident node
     * 
     * complexity: O(1) (unordered map and vector indexing are constant time,
     * Edge constructor is constant time)
     */
    Edge operator*() const {
        // ptr->index() gives the index of the node corresponding to this iterator
        // all nodes incident to this node are stored in g->incident_nodes[ptr->index()]
        /*unsigned node_index = g->incident_nodes[ptr->index()][index];
        return Edge(ptr->graph, ptr->index(), node_index);*/
        
        // both node_id and ptr->id should be the fixed ids
        unsigned node_id = *iter;
        return Edge(ptr->graph, ptr->id, node_id);
    }

    /**
     * @brief increments iterator to point to next indicent edge
     */
    IncidentIterator& operator++() {
        iter++;
        return *this;
    }

    /*
     * @brief checks whether given IncidentIterator is equal to current IncidentIterator
     *
     * @param[in] it IncidentIterator to check equality against
     * @return true if and only if iterators are equal
     *
     * @pre _it_ is a valid IndicentIterator
     */
    bool operator==(const IncidentIterator& it) const {
        return (iter == it.iter && ptr == it.ptr);
    }

   private:
    friend class Graph;
    
    // HW1 #3: YOUR CODE HERE
    Node* ptr;
    std::unordered_set<size_type>::iterator iter;

    IncidentIterator(const Node* n, std::unordered_set<size_type>::iterator i) : ptr{const_cast<Node*>(n)}, iter{i} {}

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  //class EdgeIterator : private totally_ordered<EdgeIterator> {
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /**
     * @return Edge object pointed to by iterator
     */
    Edge operator*() const {
        std::pair<size_type, size_type> ed = graph->edges_[index];
        return Edge(graph, ed.first, ed.second);
    }

    /**
     * @brief increments iterator to point to next Edge
     *
     * @post index of new edge being pointed to is one larger 
             than index of old edge
     */
    EdgeIterator& operator++() {
        index++;
        return *this;
    }

    /**
     * @brief checks equality between current iterator and _it_
     *
     * @param[in] it EdgeIterator to compare with current iterator
     * @return true if and only if _it_ is equal to current iterator
     *
     * @pre _it_ is a valid EdgeIterator
     */
    bool operator==(const EdgeIterator& it) const {
        return (graph == it.graph && index == it.index);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    Graph* graph;
    size_type index;

    EdgeIterator(const Graph* g, size_type i) : graph{const_cast<Graph*>(g)}, index{i} {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
   * @return EdgeIterator pointing to beginning of collection of edges
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /**
   * @return EdgeIterator indicating end of collection of edges
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.size());
  }

 private:

  struct node_info {
    Point loc;
    // fixed id is a feature of the nodes that is hidden to the user.
    // this serves as a way to link the edges and the nodes without 
    // having to worry about reindexing.
    // The node indices are distinct from the fixed_id. 
    //size_type fixed_id;
    size_type idx;
    node_value_type value {};
  };

  // n_encountered_nodes keeps track of how many nodes have been 
  // added to the graph over time.
  // may be much larger than n_nodes, the number of nodes currently in the graph
  size_type n_encountered_nodes = 0;
  size_type n_nodes = 0;
  size_type n_encountered_edges = 0;
  size_type n_edges = 0;
  
  // vector of "node_info"s keeps track of the position of each node, as well as 
  // a unique identifier for each node that stays constant through the lifetime 
  // of the node
  
  // nodes_ has all nodes that have ever been added
  std::vector<node_info> nodes_;
  // current_nodes has just the nodes currently in the graph
  std::vector<unsigned> i2u;

  // "edges_" stores pairs of node indices
  // pairs are always ordered, so that the first element is always smaller
  // than the second element in the pair
  std::vector<std::pair<size_type, size_type>> edges_;

  // edges are also stored in a set. This allows for O(logn) 
  // searching for elements. 
  std::set<std::pair<size_type, size_type>> edges_set_;

  // key: node indices
  // value: vector of nodes connected to key node by an edge
  //std::unordered_map<size_type, std::vector<size_type>> incident_nodes;
  std::unordered_map<size_type, std::unordered_set<size_type>> incident_nodes;

  std::unordered_map<size_type, std::unordered_map<size_type, edge_value_type>> values; 
  
};

#endif // CME212_GRAPH_HPP
