#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <functional>
#include <unordered_set>
#include <unordered_map>
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
    internal_ids = std::unordered_map<size_type,size_type>();
    next_id = 0;
    points = std::vector<Point>();
    values = std::vector<node_value_type>();
    nodes = std::vector<node_type>();
    edge_list = std::vector<edge_type>();
    edge_map = std::unordered_map<size_type,std::vector<edge_type>>();
    edge_set = std::unordered_map<size_type,std::unordered_set<size_type>>();
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
      parent = nullptr;
      id = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      return parent->points[parent->internal_ids.at(id)];
    }

    /** Return this node's value. */
    node_value_type& value() {
      return const_cast<Graph*>(parent)->values[parent->internal_ids.at(id)];
    }

    /** Return this node's value as a const. */
    const node_value_type& value() const {
      return parent->values[parent->internal_ids.at(id)];
    }

    /** Set this node's value. */
    void set_value(node_value_type v) {
        const_cast<Graph*>(parent)->values[parent->internal_ids.at(id)] = v;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return parent->internal_ids.at(id);
    }

    /** Return this node's degree, i.e. the number of edges incident to it. */
    size_type degree() const {
        return parent->edge_map.at(id).size();
    }

    /** Return an IncidentIterator pointing to the first incident edge. */
    incident_iterator edge_begin() const {
        return IncidentIterator(*this, 0, *parent);
    }

    /** Return an IncidentIterator that compares equal to an iterator that has
     *  iterated over every edge incident to this node. */
    incident_iterator edge_end() const {
        return IncidentIterator(*this, parent->edge_map.at(id).size(), *parent);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return parent == n.parent && id == n.id;
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
      if (std::less<Graph*>()(const_cast<Graph*>(n.parent),
                              const_cast<Graph*>(parent))) {
        return false;
      } else if (n.parent == parent) {
        if (n.id <= id) {
          return false;
        }
      }
      return true;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    const Graph* parent;
    size_type id;

    /** Construct a valid node, for use by Graph. */
    Node(const Graph& graph, size_type i) :
         parent(&graph), id(i) {

    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return points.size();
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
  node_type add_node(const Point& position,
                     const node_value_type& v = node_value_type()) {
    internal_ids.insert(std::pair<size_type, size_type>(next_id,points.size()));
    points.push_back(position);
    values.push_back(v);
    node_type node(*this, next_id);
    nodes.push_back(node);
    edge_map.insert(std::pair<size_type,
                    std::vector<edge_type>>(next_id,std::vector<edge_type>()));
    edge_set.insert(std::pair<size_type,
                    std::unordered_set<size_type>>(next_id,
                                                   std::unordered_set<size_type>()));
    next_id++;
    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const node_type& n) const {
    return n.parent == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
    return nodes[i];
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
      node_type n1 = Node();
      node_type n2 = Node();
      node_1 = &n1;
      node_2 = &n2;
    }

    /** Return a node of this Edge */
    node_type node1() const {
      return *node_1;
    }

    /** Return the other node of this Edge */
    node_type node2() const {
      return *node_2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((*node_1 == e.node1() && *node_2 == e.node2()) ||
              (*node_1 == e.node2() && *node_2 == e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      const node_type& n1_other = e.node1();
      if (n1_other < *node_1) {
          return false;
      } else if (*node_1 == n1_other) {
          const node_type& n2_other = e.node2();
          if (n2_other < *node_2 || *node_2 == n2_other) {
              return false;
    	  }
      }
      return true;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    const node_type* node_1;
    const node_type* node_2;

    /** Construct a valid Edge, for use by Graph. */
    Edge(const node_type& n1, const node_type& n2) {
        node_1 = &n1;
        node_2 = &n2;
    }

    /** Returns a copy of this edge, with node_1 and node_2 flipped. */
    Edge flip() {
        return Edge(*node_2, *node_1);
    }
  };

  /** Return the total number of distinct edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_list.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  edge_type edge(size_type i) const {
    return edge_list[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const node_type& a, const node_type& b) const {
    return (edge_set.at(a.id).find(b.id) != edge_set.at(a.id).end());
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
  edge_type add_edge(const node_type& a, const node_type& b) {
    if (!has_edge(a,b)) {
        edge_type new_edge(a,b);
        edge_list.push_back(new_edge);
        edge_set.at(a.id).insert(b.id);
        edge_set.at(b.id).insert(a.id);
        edge_map.at(a.id).insert(std::upper_bound(edge_map.at(a.id).begin(),
                                                  edge_map.at(a.id).end(), new_edge),
                                 new_edge); // keep edge_map in sorted order
        edge_map.at(b.id).insert(std::upper_bound(edge_map.at(b.id).begin(),
                                                  edge_map.at(b.id).end(), new_edge),
                                 new_edge); // keep edge_map in sorted order
        return new_edge;
    } else {
        return edge_type(a,b); // we know it exists, and this will create an equivalent object
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    for (unsigned int i = 0; i < nodes.size(); i++) {
        nodes[i].parent = nullptr;
    }
    internal_ids = std::unordered_map<size_type,size_type>();
    next_id = 0;
    points = std::vector<Point>();
    values = std::vector<node_value_type>();
    nodes = std::vector<node_type>();
    edge_list = std::vector<edge_type>();
    edge_map = std::unordered_map<size_type,std::vector<edge_type>>();
    edge_set = std::unordered_map<size_type,std::unordered_set<size_type>>();
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
      parent = nullptr;
      index = -1;
    }

    /** Dereference the NodeIterator. */
    value_type operator*() const {
        return parent->nodes[index];
    }

    /** Increment the NodeIterator. */
    NodeIterator& operator++() {
        index++;
        return *this;
    }

    /** Test equality between NodeIterators. */
    bool operator==(const NodeIterator& other) const {
        return (parent == other.parent && index == other.index);
    }
   private:
    friend class Graph;

    const Graph* parent;
    int index;

    /** Construct a valid NodeIterator, for use by Graph. */
    NodeIterator(const Graph& graph, int i) {
        parent = &graph;
        index = i;
    }
  };

  /** Returns a NodeIterator pointing to the first Node in the Graph. */
  node_iterator node_begin() const {
      return NodeIterator(*this, 0);
  }

  /** Returns a NodeIterator that is equal to an iterator that has iterated over
   * every Node of the Graph. */
  node_iterator node_end() const {
      return NodeIterator(*this, nodes.size());
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
        parent = nullptr;
        index = -1;
        node = Node();
    }

    /** Dereference the IncidentIterator. */
    value_type operator*() const {
        value_type edge = parent->edge_map.at(node.id)[index];
        if(edge.node1() != node) {
            return edge.flip();
        } else {
            return edge;
        }
    }

    /** Increment the IncidentIterator. */
    IncidentIterator& operator++() {
        index++;
        return *this;
    }

    /** Test equality between IncidentIterators. */
    bool operator==(const IncidentIterator& other) const {
        return (parent == other.parent && node == other.node && index == other.index);
    }

   private:
    friend class Graph;

    const Graph* parent;
    int index;
    const node_type node;

    /** Construct a valid IncidentIterator, for use by Node. */
    IncidentIterator(const node_type& n, int i,const Graph& graph) : node(n) {
        index = i;
        parent = &graph;
    }
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
        parent = nullptr;
        index = -1;
    }

    /** Dereference the EdgeIterator. */
    value_type operator*() const {
        return parent->edge_list[index];

    }

    /** Increment the EdgeIterator. */
    EdgeIterator& operator++() {
        index++;
        return *this;
    }

    /** Test equality between EdgeIterators. */
    bool operator==(const EdgeIterator& other) const {
        return (parent == other.parent && index == other.index);
    }

   private:
    friend class Graph;

    const Graph* parent;
    int index;

    /** Construct a valid EdgeIterator, for use by Graph. */
    EdgeIterator(const Graph& graph, int i) {
        parent = &graph;
        index = i;
    }
  };

  /** Returns a EdgeIterator pointing to the first Edge in the Graph. */
  edge_iterator edge_begin() const {
      return EdgeIterator(*this, 0);
  }

  /** Returns a EdgeIterator that is equal to an iterator that has iterated over
   * every Edge of the Graph. */
  edge_iterator edge_end() const {
      return EdgeIterator(*this, edge_list.size());
  }

 private:

  std::unordered_map<size_type,size_type> internal_ids;
  // map between node ids, and their location in the nodes vector.
  // key: initial id of a node.
  // value: current id of a node.

  size_type next_id;
  // keeps track of the next id to assign to an added node, to ensure uniqueness.

  std::vector<Point> points;
  // one-dimensional vector keeping track of underlying Point classes of nodes,
  // stored at the current id of their node.

  std::vector<node_value_type> values;
  // one-dimensional vector keeping track of underlying values of nodes,
  // stored at the current id of their node.

  std::vector<node_type> nodes;
  // one-dimensional vector keeping track of the already-created node objects.

  std::unordered_map<size_type,std::vector<edge_type>> edge_map;
  // map between node ids, and all the edges incident on that node.
  // key: node id.
  // value: vector of edges incident on that node.

  std::vector<edge_type> edge_list;
  // one-dimensional vector keeping track of edges and their indices.

  std::unordered_map<size_type,std::unordered_set<size_type>> edge_set;
  // map of sets of all unique edges in the graph, for easy existence checks.
  // key: node id.
  // value: set of all node ids that are connected by an edge.

};

#endif // CME212_GRAPH_HPP
