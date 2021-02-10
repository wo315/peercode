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

  unsigned num_nodes_; // same as size_
  unsigned num_edges_;

  // proxy design
  struct nodeProxy;
  struct edgeProxy;

  // node/edge vector storing its proxy elements
  std::vector <nodeProxy> nodeVector;
  std::vector <edgeProxy> edgeVector;

  // node/edge id to its unique id
  std::vector <unsigned int> node_uids;
  std::vector <unsigned int> edge_uids;

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
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
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
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(this->graph_ != NULL);
      return graph_->nodeVector[node_idx_].P;
      
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(this->graph_ != NULL);
      assert(node_idx_ < graph_->size());
      return graph_->nodeVector[node_idx_].i;
    
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    node_value_type& value() {
      return graph_->nodeVector[node_idx_].value;
    }

    const node_value_type& value() const {
      return graph_->nodeVector[node_idx_].value;
    }

    size_type degree() const {
      return graph_->nodeVector[node_idx_].adjEV.size();
    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, this, 0);
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, this, this->degree());
    } 
       
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (graph_ == n.graph_ and node_idx_ == n.node_idx_){
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
      // HW0: YOUR CODE HERE
      if (graph_ == n.graph_ and node_idx_ < n.node_idx_){
        return true;
      } else if (graph_ < n.graph_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type node_idx_;

    /** private constructor **/
    Node(const graph_type* graph, size_type node_idx)
      : graph_(const_cast<graph_type*>(graph)), node_idx_(node_idx) {   
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, const node_value_type& nval = node_value_type()) {
    // HW0: YOUR CODE HERE
    std::vector<adjEdges> adjEV{};
    nodeVector.emplace_back(num_nodes_, position, adjEV, nval);
    node_uids.push_back(nodeVector.size() - 1);
    ++num_nodes_;
    return Node(this, nodeVector.size() - 1);  
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // (void) n;            // Quiet compiler warning
    if (n.graph_ == this and n.node_idx_ < size()) {
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
    // HW0: YOUR CODE HERE
    if (i >= num_edges_) {
      throw "index invaild!";
    }
    return Node(this, node_uids[i]);
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
  class Edge: private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, n1_idx_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, n2_idx_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // (void) e;           // Quiet compiler warning
      // //HW0: YOUR CODE HERE
      assert(e.graph_ != NULL and this->graph_ != NULL);
      if (node1() == e.node1() and node2() == e.node2()) {
        return true;
      } else if (node2() == e.node1() and node1() == e.node2()) {
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
      // (void) e;           // Quiet compiler warning
      // //HW0: YOUR CODE HERE
      assert(e.graph_ != NULL and this->graph_ != NULL);
      if (graph_ < e.graph_) {
        return true;
      } else if (graph_ == e.graph_ and edge_idx_ < e.edge_idx_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type edge_idx_;
    size_type n1_idx_, n2_idx_;

    /** private constructor **/
    Edge(const graph_type* graph, 
        size_type edge_idx, 
        size_type n1_idx, 
        size_type n2_idx)
      : graph_(const_cast<graph_type*>(graph)), 
              edge_idx_(edge_idx), 
              n1_idx_(n1_idx), n2_idx_(n2_idx) {

    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i >= num_edges_) {
      throw "index invaild!";
    }
    size_type idx = edge_uids[i];
    return Edge(this, idx, edgeVector[idx].id1, edgeVector[idx].id2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // (void) a; (void) b;   // Quiet compiler warning
    assert(this == a.graph_ and this == b.graph_);
    for (size_type i = 0; i < nodeVector[a.node_idx_].adjEV.size(); ++i) {
      if (nodeVector[a.node_idx_].adjEV[i].id2 == b.node_idx_) {
        return true;
      }
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
    // HW0: YOUR CODE HERE
    assert(this == a.graph_ and this == b.graph_);
    assert(a.node_idx_ != b.node_idx_);
    if (!has_edge(a, b)) {
    edgeVector.emplace_back(num_edges_, a.node_idx_, b.node_idx_);
    edge_uids.push_back(edgeVector.size() - 1);
    nodeVector[a.node_idx_].adjEV.emplace_back(edgeVector.size()-1, b.node_idx_);
    nodeVector[b.node_idx_].adjEV.emplace_back(edgeVector.size()-1, a.node_idx_);
    num_edges_++;
    return Edge(this, edgeVector.size() - 1, a.node_idx_, b.node_idx_);       
    }
    return Edge(this, edgeVector.size() - 1, a.node_idx_, b.node_idx_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    num_nodes_ = 0;
    num_edges_ = 0;
    nodeVector.clear();
    node_uids.clear();
    edgeVector.clear();
    edge_uids.clear();
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
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    Node operator*() const {
        return graph_->node(ni_);
    }

    NodeIterator& operator++() {
        ni_++;
        return *this;
    }

    bool operator==(const NodeIterator& n) const {
        return ((graph_ == n.graph_) && (ni_ == n.ni_));
    }

    bool operator!=(const NodeIterator& n) const {
        return not ((graph_ == n.graph_) && (ni_ == n.ni_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;
    size_type ni_;
    NodeIterator(const graph_type* graph, size_type ni1) :
        graph_(const_cast<graph_type*>(graph)), ni_(ni1) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  node_iterator node_end() const {
    return node_iterator(this, this->num_nodes());

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
      size_type id = node_->node_idx_;
      return Edge(graph_, graph_->nodeVector[id].adjEV[iid_].ei, id,\
              graph_->nodeVector[id].adjEV[iid_].id2);
    }

    IncidentIterator& operator++() {
      iid_++;
      return *this;
    }

    bool operator==(const IncidentIterator& ii) const {
      return ((graph_ == ii.graph_) && (node_ == ii.node_) && (iid_ == ii.iid_));
    }

    bool operator!=(const IncidentIterator& ii) const {
      return not ((graph_ == ii.graph_) && (node_ == ii.node_) && (iid_ == ii.iid_));
    }    

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    node_type* node_;
    size_type iid_;

    IncidentIterator(const graph_type* graph, const node_type* node, size_type iid) :
      graph_(const_cast<graph_type*>(graph)), node_(const_cast<node_type*>(node)), iid_(iid) {}
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
    Edge operator*() const {
      return graph_->edge(ei_);
    }

    EdgeIterator& operator++() {
      ei_++;
      return *this;
    }

    bool operator==(const EdgeIterator& n) const {
      return ((graph_ == n.graph_) && (ei_ == n.ei_));
    } 

    bool operator!=(const EdgeIterator& n) const {
      return not ((graph_ == n.graph_) && (ei_ == n.ei_));
    }     

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type ei_;
    EdgeIterator(const graph_type* graph, size_type ei) :
      graph_(const_cast<graph_type*>(graph)), ei_(ei) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const {
    return EdgeIterator(this, this->num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct adjEdges {
    size_type ei;
    size_type id2;

    adjEdges (size_type ei1, size_type id02)
             :ei(ei1), id2(id02) {}
   };

  struct nodeProxy {
    size_type i;
    Point P;
    std::vector<adjEdges> adjEV;
    node_value_type value;

    nodeProxy (size_type ui, Point P1, std::vector<adjEdges> adjEV1, node_value_type v1) 
              :i(ui), P(P1), adjEV(adjEV1), value(v1) {}

   };

  struct edgeProxy {
    size_type ei;
    size_type id1;
    size_type id2;

    edgeProxy (size_type e_ui, size_type id01, size_type id02) 
              :ei(e_ui), id1(id01), id2(id02) {}
   };
 

};

#endif // CME212_GRAPH_HPP
