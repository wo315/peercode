#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <tuple>

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
  using graph_type = Graph<V>;

  /** Type of node value. */
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
  Graph() : points_(), num_nodes_(), num_edges_(), edge_to_nodes_(),
                                            node_to_node_to_edge_(){
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return get_point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return reference to this node's value. */
    node_value_type& value() {
      return get_value();
    }

    /** Return constant reference to this node's value. */
    const node_value_type& value() const {
      return get_value();
    }
    
    /** Return the degree of this node. */
    size_type degree() const {
      if (g_->node_to_node_to_edge_.count(index_)) {
        return g_->node_to_node_to_edge_[index_].size();
      }
      return 0; // if isolated node, return 0
    }
    
    /** Return an incident iterator pointing to the beginning of this node's
     * adjacency list i.e. the first edge. */
    incident_iterator edge_begin() const {
      return IncidentIterator(g_,
                              g_->node_to_node_to_edge_[index_].begin());
    }

    /** Return an incident iterator pointing to the end of this node's
     * adjacency list i.e. past the last edge. */
    incident_iterator edge_end() const {
      return IncidentIterator(g_,
                              g_->node_to_node_to_edge_[index_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return n.g_ == this->g_ and n.index_ == index_;
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
      // For < to be true one of the following must hold:
      // (1) if graph pointers are equal, then n.index < this.index
      // (2) n.graph_pointer < this.graph_pointer (per STD's less)
      return (n.g_ == this->g_ and n.index_ < index()) or
              std::less<Graph*>{}(n.g_, this->g_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* g_;
    size_type index_;

    // Private constructor
    Node(const Graph* g, size_type index)
      : g_(const_cast<Graph*>(g)), index_(index){
    }

    /**
     * @brief Helper method to return position associated with a node.
     * @return Point object associated with a node.
     *
     * @pre The node in question is valid.
     * @post result is the 3D Point corresponding to the node's position.
     */
    Point& get_point() const{
      assert(index_ < g_->points_.size());
      return g_->points_[index_];
    }

   /**
     * @brief Helper method to return value associated with a node.
     * @return Value associated with a node.
     *
     * @pre The node in question is valid.
     * @post result is the value corresponding to the node.
     */
   node_value_type& get_value() {
     assert(index_ < g_->v_.size());
     return g_->v_[index_];
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
  Node add_node(const Point& position, const node_value_type& v = 
                                                 node_value_type()) {
    // HW0: YOUR CODE HERE
    points_.emplace_back(position);
    v_.emplace_back(v);
    num_nodes_ = points_.size();
    return Node(this, num_nodes_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.g_ == this and n.index_ < num_nodes_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_nodes());
    return Node(this, i);
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return gr_->node(nd1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return gr_->node(nd2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return e.gr_ == this->gr_ and e.idx_ == this->idx_ and
             ((e.nd1_ == this->nd1_ and e.nd2_ == this->nd2_) or
             (e.nd1_ == this->nd2_ and e.nd2_ == this->nd1_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // For < to be true one of the following must hold:
      // (1) if graph pointers are equal, then e.index < this.index 
      // (2) e.graph_pointer < this.graph_pointer (per STD's less)
      return (e.gr_ == this->gr_ and e.idx_ < this->idx_) or
              std::less<Graph*>{}(e.gr_, this->gr_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* gr_;
    size_type idx_;
    size_type nd1_;
    size_type nd2_;

    // Private Constructor
    Edge(const Graph* gr, size_type idx, size_type nd1, size_type nd2)
      : gr_(const_cast<Graph*>(gr)), idx_(idx), nd1_(nd1), nd2_(nd2) {
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
    assert(i < num_edges_);
    return Edge(this, i, std::get<0>(edge_to_nodes_[i]),
                         std::get<1>(edge_to_nodes_[i]));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(has_node(a));
    assert(has_node(b));

    return node_to_node_to_edge_.count(a.index_) and
       node_to_node_to_edge_.at(a.index_).count(b.index_);
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
    assert(has_node(a));
    assert(has_node(b));
    assert(a.index_ != b.index_);
    
    size_type edge_ind;
    if (has_edge(a,b)){
      return Edge(this, node_to_node_to_edge_[a.index_][b.index_],
                  a.index_, b.index_);
    }
    else{
      edge_ind = num_edges_;
      edge_to_nodes_.push_back(std::make_tuple(a.index_, b.index_));
      node_to_node_to_edge_[a.index_][b.index_] = edge_ind;
      node_to_node_to_edge_[b.index_][a.index_] = edge_ind;
      num_edges_ += 1;
      return edge(edge_ind);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    points_.clear();
    v_.clear();
    num_nodes_ = 0;
    num_edges_ = 0;
    edge_to_nodes_.clear();
    node_to_node_to_edge_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered <NodeIterator> { 
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
    
    /**
     * @brief * operator for NodeIterator.
     * @return Node object of node to which the NodeIterator is pointing.
     *
     * @pre The node to which the NodeIterator is pointing is valid.
     * @post result is a valid node object.
     */
    Node operator*() const {
      return gra_->node(ind_);
    }
    
    /**
     * @brief ++ operator for NodeIterator.
     * @return Reference to NodeIterator after updation.
     *
     * @pre ++ is called within a valid range of nodes.
     * @post result points to next valid node.
     */
    NodeIterator& operator++() {
      node_it_++;
      ind_++;
      return *this;
    }

    /**
     * @brief == operator for NodeIterator.
     * @param[in] nit  A NodeIterator object.
     * @return boolean true or false.
     *
     * @pre nit must be of type NodeIterator
     * @post result is true iff _nit_ belongs to the same graph as this and
     * points to the same node element in that graph as this; else false.
     */  
    bool operator==(const NodeIterator& nit) const {
      return (nit.gra_ == this->gra_ and nit.node_it_ == this->node_it_ and
                                                    nit.ind_ == this->ind_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph* gra_;
    std::vector<Point>::const_iterator node_it_;
    size_type ind_;

    // Private constructor
    NodeIterator(const Graph* gra, std::vector<Point>::const_iterator node_it,
              size_type ind) : gra_(gra), node_it_(node_it), ind_(ind) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  /**
   * @brief Returns a NodeIterator pointing to the node with index 0 in graph.
   * @return NodeIterator pointing to the node with index 0 in graph.
   */  
  node_iterator node_begin() const {
    return NodeIterator(this, points_.begin(), 0);
  }

  /**
   * @brief Returns a NodeIterator pointing past the last node in graph.
   * @return NodeIterator pointing past the last node in graph.
   */
  node_iterator node_end() const {
    return NodeIterator(this, points_.end(), this->num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered <IncidentIterator> {
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
    /**
     * @brief * operator for IncidentIterator.
     * @return Edge object of edge to which the IncidentIterator is pointing.
     *
     * @pre The edge to which the IncidentIterator is pointing is valid.
     * @post result is a valid edge object.
     * @post the nodes associated with the edge will have the node to which
     * the incident operator is associated as node1 and the other node as node2
     */
    Edge operator*() const {
      size_type node1;
      size_type node2 = incd_it_->first;
      size_type edgeid = incd_it_->second;
      
      // look up the node_ids based on edge_id via graph
      size_type node1id = std::get<0>(grap_->edge_to_nodes_[edgeid]);
      size_type node2id = std::get<1>(grap_->edge_to_nodes_[edgeid]);
      // figure out node1 (the above obtained node<X>id's may be flipped)
      if (node1id != node2) {node1 = node1id;}
      else {node1 = node2id;}
 
      return Edge(grap_, edgeid, node1, node2);
    }

    /**
     * @brief ++ operator for IncidentIterator.
     * @return Reference to IncidentIterator after updation.
     *
     * @pre ++ is called within a valid range of node's adjacency list.
     * @post result points to next valid edge in node's adjacency list.
     */
    IncidentIterator& operator++() {
      incd_it_++;
      return *this;
    }

    /**
     * @brief == operator for IncidentIterator.
     * @param[in] init  A IncidentIterator object.
     * @return boolean true or false.
     *
     * @pre init must be of type IncidentIterator.
     * @post result is true iff _init_ belongs to the same graph as this and
     * points to the same edge element in the same node's adjacency list as
     * this; else false.
     */
    bool operator==(const IncidentIterator& init) const {
      return (init.grap_ == this->grap_ and init.incd_it_ == this->incd_it_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* grap_;
    std::map<size_type, size_type>::const_iterator incd_it_;

    // Private constructor
    IncidentIterator(const Graph* grap,
                     std::map<size_type, size_type>::const_iterator incd_it)
                     : grap_(grap), incd_it_(incd_it) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered <EdgeIterator> {
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
    
    /**
     * @brief * operator for EdgeIterator.
     * @return Edge object of edge to which the EdgeIterator is pointing.
     *
     * @pre The edge to which the EdgeIterator is pointing is valid.
     * @post result is a valid edge object.
     */
    Edge operator*() const { 
      return graph_->edge(eind_);
    } 
 
    /**
     * @brief ++ operator for EdgeIterator.
     * @return Reference to EdgeIterator after updation.
     *
     * @pre ++ is called within a valid range of graph's edges.
     * @post result points to next valid edge in graph, or
     * after all edges have been seen, past the last edge in the graph.
     * @post the edge pointed to is being seen for the first time, i.e. if
     * an edge object corresponding to a->b was previously returned, an object
     * for b->a will not be returned.
     */
    EdgeIterator& operator++() {
      eit_++;
      eind_++;
      return *this;
    }
    
    /**
     * @brief == operator for EdgeIterator.
     * @param[in] eit  A EdgeIterator object.
     * @return boolean true or false.
     *
     * @pre eit must be of type EdgeIterator.
     * @post result is true iff _eit_ belongs to the same graph as this and
     * points to the same edge element as * this; else false.
     */
    bool operator==(const EdgeIterator& eit) const {
      return (eit.graph_ == this->graph_ and eit.eit_ == this->eit_ and
                                                     eit.eind_ == this->eind_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_;
    std::vector<std::tuple<size_type,size_type>>::const_iterator eit_;
    size_type eind_;
    
    // Private constructor
    EdgeIterator(const Graph* graph,
              std::vector<std::tuple<size_type,size_type>>::const_iterator eit,
              size_type eind) : graph_(graph), eit_(eit), eind_(eind) {}

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /** Return an edge iterator pointing to the first edge in the graph. */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, this->edge_to_nodes_.begin(), 0);
  }
  
  /** Return an edge iterator pointing past the last edge in the graph. */
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->edge_to_nodes_.end(), this->num_edges());
  }
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point> points_;
  std::vector<node_value_type> v_;
  size_type num_nodes_;
  size_type num_edges_;
  std::vector<std::tuple<size_type,size_type>> edge_to_nodes_;
  std::map<size_type, std::map<size_type, size_type>> node_to_node_to_edge_;
};

#endif // CME212_GRAPH_HPP
