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
template<typename V, typename E>
class Graph {
 private:

  unsigned int num_nodes_, num_edges_;
  struct neighbors;
 
  // These vectors let us go from the id of a node/edge directly to its position in node_vector/edge_vector
  std::vector <unsigned int> node_translate, edge_translate;

  struct pnode;
  std::vector <pnode> node_vector;

  struct pedge;
  std::vector <pedge> edge_vector;

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

  typedef V node_value_type;
  typedef E edge_value_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    num_nodes_ = 0;
    num_edges_ = 0;
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
      return graph_->node_vector[global_index_].point;
    }

    /** Return this node's position. */
    Point& position() {
      return graph_->node_vector[global_index_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->node_vector[global_index_].idx;
    }


    /**
     * @brief Get value of a node.
     *
     * @return node_value.
     */
    node_value_type& value() {
      return graph_->node_vector[global_index_].val;
    }

    /**
     * @brief Get value of a node.
     *
     * @return immutable node_value.
     */
    const node_value_type& value() const {
      return graph_->node_vector[global_index_].val;
    }

    /**
     * @brief Get degree of a node.
     *
     * @return number of neighbors of current node.
     */
    size_type degree() const {
      return graph_->node_vector[global_index_].adj_nodes.size();
    }

    /**
     * @brief Create iterator for start of Edge sequence neighboring current node.
     *
     * @return EdgeIterator object pointing to first element of the sequence.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, this, 0);
    }

    /**
     * @brief Create iterator for end of Edge sequence neighboring current node.
     *
     * @return IncidentIterator object pointing to last element of the sequence.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, this, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((graph_ == n.graph_) and (global_index_ == n.global_index_));
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
      return (global_index_ < n.global_index_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type global_index_;
    Node(const graph_type* graph, size_type id)
      :graph_(const_cast<graph_type*>(graph)), global_index_(id) {}
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_nodes_;
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
  Node add_node(const Point& position, const node_value_type& node_val = node_value_type()) {
    std::vector<neighbors> adj_nodes;
    node_vector.emplace_back(num_nodes_, position, adj_nodes, node_val);
    unsigned int id = node_vector.size() - 1;
    node_translate.push_back(id);
    num_nodes_++;
    return Node(this, id); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (node(n.index()) == n);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i < size()) {
      return Node(this, node_translate[i]);
    }
    // If this statement is reached, the input violated preconditions.
    return Node();
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, nid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, nid2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (((node1() == e.node1()) and (node2() == e.node2())) or
          ((node2() == e.node1()) and (node1() == e.node2())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return ((graph_ != e.graph_) or (eid_ < e.eid_) or (node1() < e.node1()));
    }

    /**
     * @brief Get length of an edge.
     *
     * @return distance from node2 to node1 of a given edge.
     */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /**
     * @brief Get value of an edge.
     *
     * @return value correspopnding to given edge.
     */
    edge_value_type& value() {
      return graph_->edge_vector[eid_].edge_val;
    }

    /**
     * @brief Get value of an edge.
     *
     * @return value correspopnding to given edge.
     */
    const edge_value_type& value() const {
      return graph_->edge_vector[eid_].edge_val;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* graph_;
    size_type eid_, nid1_, nid2_;
    Edge(const graph_type* graph, size_type eid, size_type nid1, size_type nid2)
      :graph_(const_cast<graph_type*>(graph)), eid_(eid),nid1_(nid1), nid2_(nid2) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i < num_edges_) {
      size_type id = edge_translate[i];
      size_type nid1 = edge_vector[id].id1;
      size_type nid2 = edge_vector[id].id2;
      return Edge(this, id, nid1, nid2);
    }
    // If this statement is reached, the input violated preconditions.
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    unsigned int num_adj_nodes = node_vector[a.global_index_].adj_nodes.size();
    for (unsigned int i = 0; i < num_adj_nodes; i++) {
      if (node_vector[a.global_index_].adj_nodes[i].nid == b.global_index_)
        return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& edge_val = edge_value_type()) {
    if (!has_edge(a, b)) {
      edge_vector.emplace_back(num_edges_, a.global_index_, b.global_index_, edge_val);
      unsigned int eid = edge_vector.size() - 1;
      edge_translate.push_back(eid);
      num_edges_++;
      // Add edge from a to b and from b to a
      node_vector[a.global_index_].adj_nodes.emplace_back(eid, b.global_index_);
      node_vector[b.global_index_].adj_nodes.emplace_back(eid, a.global_index_);
      return Edge(this, eid, a.global_index_, b.global_index_);
    }
    return Edge(this, edge_vector.size() - 1, a.global_index_, b.global_index_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Reset all nodes
    num_nodes_ = 0;
    node_translate.clear();
    node_vector.clear();
    // Reset all edges
    num_edges_ = 0;
    edge_translate.clear();    
    edge_vector.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private equality_comparable<NodeIterator>  {
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

    /**
     * @brief Dereference pointer operator for NodeIterator.
     *
     * @return Node object.
     */
    Node operator*() const {
      return graph_->node(global_index_);
    }

    /**
     * @brief Increment operator for NodeIterator.
     *
     * @return NodeIterator object.
     */
    NodeIterator& operator++() {
      global_index_++;
      return *this;
    }

    /**
     * @brief Equality operator for NodeIterator.
     * @param n NodeIterator to compare this object to.
     *
     * @return boolean indicating whether NodeIterator is equal to this object.
     */
    bool operator==(const NodeIterator& n) const {
      return ((graph_ == n.graph_) and (global_index_ == n.global_index_));
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type global_index_;

    /**
     * @brief Constructor for NodeIterator.
     */
    NodeIterator(const graph_type* graph, size_type id)
          :graph_(const_cast<graph_type*>(graph)), global_index_(id) {}

  };

  /**
   * @brief Create iterator for start of Node sequence.
   *
   * @return NodeIterator object pointing to first element of sequence.
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /**
   * @brief Create iterator for end of Node sequence.
   *
   * @return NodeIterator object pointing to last element of sequence.
   */
  node_iterator node_end() const {
    return NodeIterator(this, this->num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private equality_comparable<IncidentIterator> {
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

    /**
     * @brief Dereference pointer operator for IncidentIterator.
     *
     * @return Edge object.
     */    
    Edge operator*() const {
      size_type curr_id = node_->global_index_;
      return Edge(graph_, graph_->node_vector[curr_id].adj_nodes[id_].eid,
                  curr_id, graph_->node_vector[curr_id].adj_nodes[id_].nid);
    }
      
    /**
     * @brief Increment operator for IncidentIterator.
     *
     * @return IncidentIterator object.
     */
    IncidentIterator& operator++() {
      id_++;
      return *this;
    }

    /**
     * @brief Equality operator for IncidentIterator.
     * @param i IncidentIterator to compare this object to.
     *
     * @return boolean indicating whether IncidentIterator is equal to this object.
     */
    bool operator==(const IncidentIterator& i) const {
      return ((graph_ == i.graph_) and (node_ == i.node_) and (id_ == i.id_));
    }

   private:
    friend class Graph;
    graph_type* graph_;
    node_type* node_;
    size_type id_;

    /**
     * @brief Constructor for IncidentIterator.
     */
    IncidentIterator(const graph_type* graph, const node_type* node, size_type id) :
    graph_(const_cast<graph_type*>(graph)), node_(const_cast<node_type*>(node)), id_(id) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private equality_comparable<EdgeIterator> {
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

    /**
     * @brief Dereference pointer operator for EdgeIterator.
     *
     * @return Edge object.
     */    
    Edge operator*() const {
      return graph_->edge(id_);
    }
      
    /**
     * @brief Increment operator for EdgeIterator.
     *
     * @return EdgeIterator object.
     */
    EdgeIterator& operator++() {
      id_++;
      return *this;
    }

    /**
     * @brief Equality operator for EdgeIterator.
     * @param e EdgeIterator to compare this object to.
     *
     * @return boolean indicating whether EdgeIterator is equal to this object.
     */
    bool operator==(const EdgeIterator& e) const {
      return ((graph_ == e.graph_) and (id_ == e.id_));
    }

   private:
    friend class Graph;
    graph_type* graph_;
    size_type id_;

    /**
     * @brief Constructor for EdgeIterator.
     */
    EdgeIterator(const graph_type* graph, size_type id) :
      graph_(const_cast<graph_type*>(graph)), id_(id) {}
  };

  /**
   * @brief Create iterator for start of Edge sequence.
   *
   * @return EdgeIterator object pointing to first element of sequence.
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  };
  
  /**
   * @brief Create iterator for end of Edge sequence.
   *
   * @return EdgeIterator object pointing to last element of sequence.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->num_edges());
  };

  /** Remove node from the graph, if it exists.
  * @pre @a n is a node of this graph
  * @return 1 if n is a valid node, else return 0
  * @post has_node(n) == false
  * @post if the node n was removed, num_nodes() is reduced by 1. Else, it is unchanged.
  *
  * This may invalidate indices of vectors node_translate and adj_nodes.
  * This may not invalidate the indices for the vector node_vector.
  *
  * Complexity is O(num_nodes()).
  */
  size_type remove_node(const Node& n) {
    if (has_node(n)) {
      while (n.degree() > 0) {
        auto e = *(n.edge_begin());
        remove_edge(e);
      }
      node_translate.erase(node_translate.begin() + n.index());
      for (size_type k = n.index(); k < node_translate.size(); ++k) {
        node_vector[node_translate[k]].idx = k;
      }
      num_nodes_ --;
      return 1;
    }
    return 0;
  }

  /** Remove node from the graph, and return a new iterator pointing to the next node.
  * @pre @a n_it is an node_iterator
  * @return a new iterator that points to the next node
  * @post if @a *n_it was removed, num_nodes() is reduced by 1. Else, it is unchanged.
  * @post n_it is not a valid edge iterator
  *
  * This may invalidate indices of vectors node_translate and adj_nodes.
  * This may not invalidate the indices for the vector node_vector.
  *
  * Complexity is O(num_nodes()).
  */
  node_iterator remove_node(node_iterator n_it) {
    auto n = *(n_it);
    remove_node(n);
    node_iterator next_it = n_it;
    return next_it;
  }

  /** Remove edge from the graph, if it exists.
  * @pre @a n1 @a n2 are valid nodes of this graph
  * @return 1 if there exists edge between n1, n2, else return 0
  * @post if the edge between n1 and n2 was removed, num_edges() is reduced by 1. Else, it is unchanged.
  * @post has_edge(n1, n2) == false
  *
  * This may invalidate indices of vectors edge_translate and adj_nodes.
  * This may not invalidate the indices for the vector edge_vector.
  *
  * Complexity is O(num_nodes() + num_edges()).
  */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (has_edge(n1, n2)) {
      for (size_type i = 0; i < n1.degree(); ++i) {
        size_type n1idx = n1.global_index_;
        size_type n2idx = n2.global_index_;
        if (n2idx == node_vector[n1idx].adj_nodes[i].nid) {
          size_type knid = node_vector[n1idx].adj_nodes[i].eid;
          size_type keid = edge_vector[knid].eidx;

          node_vector[n1idx].adj_nodes.erase(
            node_vector[n1idx].adj_nodes.begin() + i);
          
          for (size_type j = 0; j < n2.degree(); ++j) {
            if (node_vector[n2idx].adj_nodes[j].nid == n1idx) {
              node_vector[n2idx].adj_nodes.erase(
                node_vector[n2idx].adj_nodes.begin() + j);
            }
          }

          edge_translate.erase(edge_translate.begin() + keid);
          for (size_type k = keid; k < edge_translate.size(); ++k) {
            edge_vector[edge_translate[k]].eidx = k;
          }
          num_edges_--;
          return 1;
        }
      }
    }
  return 0;
  }

  /** Remove edge from the graph, if it exists.
  * @pre @a e is a valid edge of this graph
  * @return 1 if there exists edge between e.node1(), e.node2(), else return 0
  * @post if the edge between e.node1() and e.node2() was removed, num_edges() is reduced by 1. Else, it is unchanged.
  * @post has_edge(e.node1(), e.node2()) == false
  *
  * This may invalidate indices of vectors edge_translate and adj_nodes.
  * This may not invalidate the indices for the vector edge_vector.
  *
  * Complexity is O(num_nodes() + num_edges()).
  */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from the graph, and return a new iterator pointing to the next edge
  * @pre @a e_it is an edge_iterator
  * @return a new iterator that points to the next node
  * @post e_it is not a valid edge iterator
  * @post if @a *e_it was removed, num_edges() is reduced by 1. Else, it is unchanged.
  *
  * This may invalidate indices of vectors edge_translate and adj_nodes
  * This may not invalidate the indices for the vector edge_vector
  *
  * Complexity is O(num_nodes() + num_edges()).
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    auto e = *(e_it);
    remove_edge(e);
    edge_iterator next_it = e_it;
    return next_it;
  }

 private:

  // Proxy edge
  struct pedge {
    size_type eidx, id1, id2;
    edge_value_type edge_val;
    pedge(size_type eidx, size_type id1, size_type id2, edge_value_type edge_val) 
      :eidx(eidx), id1(id1), id2(id2), edge_val(edge_val) {}
  };

  // Proxy node
  struct pnode {
    size_type idx;
    Point point;
    std::vector<neighbors> adj_nodes;
    node_value_type val;
    pnode(size_type idx, Point point, std::vector<neighbors> adj_nodes, node_value_type val) 
      :idx(idx), point(point), adj_nodes(adj_nodes), val(val) {}
  };

  // Adjacent nodes
  struct neighbors {
     size_type eid, nid;
     neighbors(size_type eid, size_type nid):eid(eid), nid(nid) {}
  };

};

#endif // CME212_GRAPH_HPP
