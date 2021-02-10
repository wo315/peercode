#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <tuple>
#include <map>
#include <vector>

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
  class Node: private totally_ordered<Graph::Node> {
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
    Node() : graph_(nullptr), index_(0) {
    }

    /** Return this node's position. */
    const Point& position() const {
      // TODO: @pre graph exists (not a nullptr)
      return graph_->nodes_.at(index_).position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return index_;
    }

    
    /**
    * @brief Return a non-const reference to the value of a node.
    *
    * @return A non-const reference to the _value__ attribute of the corresponding _internal_node_.
    * 
    * @pre *_this_ is a valid _Node_, i.e. _graph__ != nullptr and 0 <= _index__ < _graph__.size().
    */
    node_value_type& value() {
      return const_cast<Graph*>(graph_)->nodes_.at(index_).value_;
    }

    /**
    * @brief Return a const reference to the value of a node.
    *
    * @return A const reference to the _value__ attribute of the corresponding _internal_node_.
    * 
    * @pre *_this_ is a valid _Node_, i.e. _graph__ != nullptr and 0 <= _index__ < _graph__.size().
    * @post The state of *_this_ was not altered.
    */    
    const node_value_type& value() const {
       return graph_->nodes_.at(index_).value_;   
    }

    /**
    * @brief Get the number of edges incident to a node.
    *
    * @return A size_type value representing the number of edges incident to a node. 
    * 
    * @pre *_this_ is a valid _Node_, i.e. _graph__ != nullptr and 0 <= _index__ < _graph__.size().
    * @post The state of *_this_ was not altered.
    */ 
    size_type degree() const {
      auto it = graph_->node_to_edge_.find(*this);
      if (it == graph_->node_to_edge_.end()) {
        return 0;
      }
      return it->second.size();
    }

    /**
    * @brief Get the begin iterator of edges incident to a node.
    *
    * @return An incident_iterator set to its begin value.
    * 
    * @pre *_this_ is a valid _Node_, i.e. _graph__ != nullptr and 0 <= _index__ < _graph__.size().
    * @post The state of *_this_ was not altered.
    */   
    incident_iterator edge_begin() const {
      return incident_iterator(graph_, *this, 0);
    }

    /**
     * @brief Get the end iterator of edges incident to a node.
     *
     * @return An incident_iterator set to its end value.
     * 
     * @pre *_this_ is a valid _Node_, i.e. _graph__ != nullptr and 0 <= _index__ < _graph__.size().
     * @post The state of *_this_ was not altered.
     */   
    incident_iterator edge_end() const {
       return incident_iterator(graph_, *this, degree());     
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // TODO: @pre same graph
      if (n.index_ == index_) {
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
      // TODO: @pre same graph
      if (index_ < n.index_) {
        return true;
      }
      return false;
    }

    /**
     * @brief Get the vector of edges incident to a node.
     *
     * @param[in] _n_ A const reference to a _Node_ object from which to get the vector of incident edges.
     *
     * @return A std::vector<Edge> that is a value fetched from the _graph__->_node_to_edge__ map
     * using *_this_ as a key.
     * 
     * @pre *_this_ is a valid _Node_, i.e. _graph__ != nullptr and 0 <= _index__ < _graph__.size().
     * @pre _n_ is a valid _Node_, i.e. _n_._graph__ != nullptr and 0 <= _n_._index__ < _n_._graph__.size().
     * @post The state of *_this_ was not altered.
     * @post The state of _n_ was not altered.
     */   
    std::vector<Edge> get_inc_edges(const Node& n) const {
      return graph_->node_to_edge_.at(n);
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the graph
    const graph_type* graph_;
    // Node index
    size_type index_;

    // Private constructor to allow the graph to instantiate valid node objects
    Node(const graph_type* graph, size_type index)
        : graph_(graph),
          index_(index) {
    }

  }; 

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
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
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type()) {
    
    nodes_.emplace_back(position, value);
    return Node(this, size()-1);
  }



  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (this == n.graph_) { 
      if(size() > n.index_) {
        return true;
      }
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
  class Edge: private totally_ordered<Graph::Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() : index_(0),
             node_pair_(std::make_tuple(Node(), Node())) {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return std::get<0>(node_pair_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return std::get<1>(node_pair_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // TODO: @pre: Edge objects that have the same nodees have the same index
      Node my_node1 = std::get<0>(node_pair_);
      Node my_node2 = std::get<1>(node_pair_);
      Node e_node1  = std::get<0>(e.node_pair_);
      Node e_node2  = std::get<1>(e.node_pair_);
      if(((my_node1==e_node1) and (my_node2==e_node2)) or 
        ((my_node1==e_node2) and (my_node2==e_node1))) {
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
      Node my_node1 = std::get<0>(node_pair_);
      Node my_node2 = std::get<1>(node_pair_);
      Node e_node1  = std::get<0>(e.node_pair_);
      Node e_node2  = std::get<1>(e.node_pair_);
      if (std::min(my_node1, my_node2) == std::min(e_node1, e_node2)){
        return (std::max(my_node1, my_node2) < std::max(e_node1, e_node2));
      }
      return (std::min(my_node1, my_node2) < std::min(e_node1, e_node2));
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    size_type index_;
    std::tuple<Node, Node> node_pair_;

    // Private constructor to allow the graph to instantiate valid edge objects
    Edge(size_type index, std::tuple<Node, Node> node_pair): 
                                              index_(index),
                                              node_pair_(node_pair) {
    }
    
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    for(auto e: edges_){
      if (e.index_ == i) {return e;}
    }
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // TODO: this solution is not very elegant since I create a false edge
    // in order to use operator== on edges of the graph
    Edge tmp_edge(0, std::make_tuple(a,b));
    auto it = edges_.find(tmp_edge);
    if (it == edges_.end()) {return false;}
    return true;
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
    // check if edge is in the graph
    Edge edge(num_edges(), std::make_tuple(a,b));
    auto it = edges_.find(edge);
    // if it is in the graph
    if (it != edges_.end()) {return *it;}
    // else add the edge to the graph and return it
    node_to_edge_[a].push_back(edge);
    node_to_edge_[b].push_back(edge);
    edges_.emplace(edge);
    return edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    node_to_edge_.clear();
    return;
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

    /**
     * @brief Get the _Node_ object associated with this node iterator.
     *
     * @return  If old _index__ == graph_.size(): return an invalid _Node_,
     * else: return the _Node_ object associated with this node iterator.
     * 
     * @post The state of *_this_ was not altered.
     */   
    Node operator*() const {
      if (index_ == graph_->size()) {
        return Node();
      }
      return graph_->node(index_);
    }

    /**
     * @brief Increment the node iterator by 1 in place.
     *
     * @return A non-const reference to the *_this_.
     * 
     * @post If old _index__ == graph_.size():  old _index__ == new _index__,
     * else: old _index__ + 1 == new _index__. 
     */
    NodeIterator& operator++() {
      if (index_ != graph_->size()) {
        index_++;
      }
      return *this;
    }


    /**
     * @brief Check equality with another node iterator.
     *
     * @param[in] _n_ const reference to a NodeIterator object.
     *
     * @return Return the result of _index__ == _n_._index__.
     * 
     * @post The state of *_this_ was not altered.
     */
    bool operator==(const NodeIterator& n) const {
      // TODO: @pre same graph
      return n.index_ == index_;
    }

    /**
     * @brief Check the not-equality with another node iterator.
     *
     * @param[in] _n_ const reference to a NodeIterator object.
     * 
     * @return Return the result of _index__ != _n_._index__.
     * 
     * @post The state of *_this_ was not altered.
     */
    bool operator!=(const NodeIterator& n) const {
      // TODO: @pre same graph
      return n.index_ != index_;
    }

   private:
    friend class Graph;

    const graph_type* graph_;
    size_type index_;

    NodeIterator(const graph_type* graph,
                 size_type index):
                   graph_(graph),
                   index_(index) {}
  };

  /**
   * @brief Get the begin node iterator of this graph.
   *
   * @return Return the begin node_iterator of this graph.
   * 
   * @post The state of *_this_ was not altered.
   */
  node_iterator node_begin() const {
    return node_iterator(this, 0);  
  }

  /**
   * @brief Get the end node iterator of this graph.
   *
   * @return Return the end node_iterator of this graph.
   * 
   * @post The state of *_this_ was not altered.
   */
  node_iterator node_end() const {
    return node_iterator(this, size());
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


    /**
     * @brief Get the Edge object associated with a node's incident iterator.
     *
     * @return If _i__ == _node__.degree(): return an invalid Edge,
     * else: return the _i__th incident Edge of _node__ as stored in 
     * _graph__->_node_to_edge__.
     * 
     * @post The state of *_this_ was not altered.
     */
    Edge operator*() const {
      if (i_ == node_.degree()) {
        return Edge();
      }
      Edge edge = graph_->node_to_edge_.at(node_).at(i_);
      Node node1 = std::get<0>(edge.node_pair_);
      if (node_ == node1) {return edge;}
      return Edge(edge.index_, std::make_tuple(node_, node1));
    }

    /**
     * @brief Increment the incident iterator by 1 in place.
     *
     * @return A non-const reference to the *_this_.
     * 
     * @post If old _i__ == _node__.degree():  old _i__ == new _i__,
     * else: old _i__ + 1 == new _i__. 
     */ 
    IncidentIterator& operator++() {
      if (i_ != node_.degree()) {
        i_++;
      }
      return *this;
    }

    /**
     * @brief Check equality with another incident iterator.
     *
     * @param[in] _it_ const reference to a IncidentIterator object.
     *
     * @return Return the result of (_i__ == _it_._i__) and (_node__ == _it_._node__).
     * 
     * @post The state of *_this_ was not altered.
     */
    bool operator==(const IncidentIterator& it) const {
      if (node_ != it.node_) {return false;}
      return (i_==it.i_);
    }

   /**
     * @brief Check not-equality with another incident iterator.
     *
     * @param[in] _it_ const reference to a IncidentIterator object.    
     *
     * @return Return the result of (_i__ != _it_._i__) or (_node__ != _it_._node__).
     * 
     * @post The state of *_this_ was not altered.
     */
    bool operator!=(const IncidentIterator& it) const {
      if (node_ != it.node_) {return true;}
      return (i_ != it.i_);
   }

   private:
    friend class Graph;

    const graph_type* graph_;
    node_type node_;
    size_type i_;


    IncidentIterator(const graph_type* graph,
                     node_type node,
                     size_type i):
                     graph_(graph),
                     node_(node),
                     i_(i) {
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

    /**
     * @brief Get the Edge object associated with a graph's edge iterator.
     *
     * @return The value dereferenced by the wrapped std::set<Edge>::iterator _it__, i.e. *_it__.
     * 
     * @post The state of *_this_ was not altered.
     */ 
    Edge operator*() const {
      return *it_;
    }


    /**
     * @brief Increment this edge iterator by 1 in place.
     * 
     * @post old ++_it__ == new _it__. 
     */
    EdgeIterator& operator++() {
      it_++;
      return *this;
    }

    /**
     * @brief Check equality with another edge iterator.
     *
     * @param[in] _it_ const reference to a EdgeIterator object.
     *
     * @return Return the result of (_it__ == _it_._it__).
     * 
     * @post The state of *_this_ was not altered.
     */
    bool operator==(const EdgeIterator& it) const {
      return (it_ == it.it_);
    }

    /**
     * @brief Check the not-equality with another edge iterator.
     *
     * @param[in] _it_ const reference to a EdgeIterator object.
     *
     * @return Return the result of (_it__ != _it_._it__).
     * 
     * @post The state of *_this_ was not altered.
     */
    bool operator!=(const EdgeIterator& it) const {
      return (it_ != it.it_);
    }


   private:
    friend class Graph;
    
    typename std::set<Edge>::iterator it_;
    
    EdgeIterator(typename std::set<Edge>::iterator it):
                 it_(it) {
    } 
    
  };

  /**
   * @brief Get the begin edge iterator of this graph.
   *
   * @return Return the begin edge_iterator of this graph.
   * 
   * @post The state of *_this_ was not altered.
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(edges_.begin());
  }

   /**
   * @brief Get the end edge iterator of this graph.
   *
   * @return Return the end edge_iterator of this graph.
   * 
   * @post The state of *_this_ was not altered.
   */
  edge_iterator edge_end() const {
    return EdgeIterator(edges_.end());
  }


 private:

  // internal_node struct
  struct internal_node {

    internal_node(Point position): 
                  position_(position) { 
    }

    internal_node(Point position,
                  const node_value_type& value): 
                  position_(position),
                  value_(value) { 
    }

    Point position_;
    node_value_type value_;
  };

  // internal nodes, edges and node-to-edge map
  std::vector<internal_node> nodes_;
  std::set<Edge> edges_;
  std::map<Node,std::vector<Edge>> node_to_edge_;

};

#endif // CME212_GRAPH_HPP
