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
template <typename V, typename E>
class Graph {
 private:
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Node And Edge Value Types */
  typedef V node_value_type;
  typedef E edge_value_type;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_(), next_nid_(0), edges_(), edges_direct_() {
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

    // Copy constructor
    Node(const Node& n)
    	: graph_(n.graph_), nid_(n.nid_) {
    }

    size_type degree() const {
      return (graph_->edges_)[nid_].size(); 
    }

    incident_iterator edge_begin() const {
      return incident_iterator(graph_, (graph_->edges_)[nid_].begin());
    }

    incident_iterator edge_end() const {
      return incident_iterator(graph_, (graph_->edges_)[nid_].end());
    }

    node_value_type& value() {
      return fetch_node().nvalue;
    }
    const node_value_type& value() const {
      return fetch_node().nvalue;
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch_node().location;
    }

    // modifiable reference
    Point& position() {
    	return fetch_node().location;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      for (size_type i = 0; i < graph_->size(); ++i)
    	if (graph_->nodes_[i].nid == nid_)
    	  return i;
      assert(false);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (nid_ == n.nid_);
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
      if (nid_ < n.nid_)
      	return true;
      if ((nid_ == n.nid_) and (graph_ < n.graph_))
      	return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Pointer back to the graph
    Graph* graph_;
    // This node's unique node identifier number
    size_type nid_;
    // Private constructor
    Node(const Graph* graph, size_type nid)
    	: graph_(const_cast<Graph*>(graph)), nid_(nid) {
    }
    // Helper method to return the appropriate node
    internal_node& fetch_node() const {
    	for (size_type i = 0; i < graph_->size(); ++i)
    		if (graph_->nodes_[i].nid == nid_)
    			return graph_->nodes_[i];
    	assert(false);
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
  	internal_node new_node;
  	new_node.location = position;
    new_node.nid = next_nid_++;
    new_node.nvalue = val;
    nodes_.push_back(new_node);
    edges_[new_node.nid] = std::map<size_type, internal_edge>();
    return Node(this, new_node.nid);
  }

  /** Remove a node to the graph. Return 1 if removal was successful.
   * @param[in] n The node to remove.
   * @post The graph size decreases by 1.
   * @post All edges connected to n are removed.
   * @post All node iterators should be considered invalidated.
   * Complexity: O(num_nodes()) amortized operations.
   */
  size_type remove_node(const Node& n) {
  	// Remove all adjacent edges
  	std::vector<Edge> adjacent_edges;
  	for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
  		adjacent_edges.push_back(*it);
  	}
  	for (auto it = adjacent_edges.begin(); it != adjacent_edges.end(); ++it) {
  		remove_edge(*it);
  	}
  	// Remove this node
  	for (auto it = nodes_.begin(); it != nodes_.end(); ++it) {
  		if ((*it).nid == n.nid_) {
  			nodes_.erase(it);
  			break;
  		}
  	}
  	edges_.erase(n.nid_);
  	return 1;
  }

  /** Remove a node to the graph. Return 1 if removal was successful.
   * @param[in] n_it A node iterator pointing to the node of interest.
   * @post The graph size decreases by 1.
   * @post All edges connected to n are removed.
   * @post All node iterators should be considered invalidated.
   * Complexity: O(num_nodes()) amortized operations.
   */
  node_iterator remove_node(node_iterator n_it) {
  	remove_node(*n_it);
  	return n_it;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    for (size_type i = 0; i < size(); ++i)
    	if ((nodes_[i].nid == n.nid_) && (n.graph_ == this))
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
    return Node(this, nodes_[i].nid);
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

    edge_value_type& value() {
      return fetch_edge().evalue;
    }
    const edge_value_type& value() const {
      return fetch_edge().evalue;
    }

    // Return length of this edge as distance between its two nodes
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(fetch_edge().node1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(fetch_edge().node2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((graph_ == e.graph_) and ((nid1_ == e.nid1_ and nid2_ == e.nid2_) or (nid1_ == e.nid2_ and nid2_ == e.nid1_)))
      	return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ < e.graph_)
      	return true;
      if (graph_ > e.graph_)
      	return false;
      size_type e_lowernode, e_highernode, this_lowernode, this_highernode;
      if (e.nid1_ < e.nid2_) {
      	e_lowernode = e.nid1_;
      	e_highernode = e.nid2_;
      } else {
      	e_lowernode = e.nid2_;
      	e_highernode = e.nid1_;
      }
      if (nid1_ < nid2_) {
      	this_lowernode = nid1_;
      	this_highernode = nid2_;
      } else {
      	this_lowernode = nid2_;
      	this_highernode = nid1_;
      }
      if (e_lowernode == this_lowernode) {
      	return (this_highernode < e_highernode);
      }
      return (this_lowernode < e_lowernode);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Pointer back to the graph
    Graph* graph_;
    // This edge's node identifiers
    // NOTE: NID1_ GUARANTEED TO BE < NID2_
    size_type nid1_;
    size_type nid2_;
    // Private constructor
    Edge(const Graph* graph, size_type nid1, size_type nid2)
      : graph_(const_cast<Graph*>(graph)), nid1_(nid1), nid2_(nid2) {
    }
    // Helper method to return the appropriate edge
    internal_edge& fetch_edge() const {
      if (graph_->edges_.count(nid1_) > 0)
      	if (graph_->edges_[nid1_].count(nid2_) > 0)
      		return graph_->edges_[nid1_][nid2_];
      assert(false);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
  	return edges_direct_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
  	assert(i < num_edges());
  	internal_edge ie = *(edges_direct_[i]);
  	return Edge(this, ie.node1.nid_, ie.node2.nid_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   *
   * NOTE: THIS RELIES ON MY IMPLEMENTATION OF THIS GRAPH AS A SYMMETRIC
   *       ADJACENCY MATRIX
   */
  bool has_edge(const Node& a, const Node& b) const {
    if ((edges_.count(a.nid_) > 0) && (edges_.at(a.nid_).count(b.nid_) > 0))
      return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
  	if (has_edge(a, b)) {
  		return Edge(this, a.nid_, b.nid_);
  	}
  	internal_edge new_edge;
  	new_edge.node1 = a;
    new_edge.node2 = b;
    new_edge.evalue = val;
  	edges_[a.nid_][b.nid_] = new_edge; // SYMMETRIC ADJACENCY MATRIX
  	edges_[b.nid_][a.nid_] = new_edge; // SYMMETRIC ADJACENCY MATRIX
  	edges_direct_.push_back(&edges_[a.nid_][b.nid_]);
  	return Edge(this, a.nid_, b.nid_);
  }

  /** Remove an edge in the graph. Return 1 if removal was successful.
   *  0 if unsuccessful.
   * @param[in] n1 One of the nodes in the edge to delete.
   * @param[in] n2 One of the nodes in the edge to delete.
   * @post If edge (n1, n2) or (n2, n1) was in the graph, it will be removed.
   * @post num_edges() will decrease by 1.
   * @post All edge iterators and incident iterators will be considered invalidated.
   * Complexity: O(num_edges()) amortized operations.
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
  	// If this edge does not exist, do nothing
  	if (!has_edge(n1, n2)) {
  		return 0;
  	}
  	// Get pointer to internal_edge
  	internal_edge* this_edge = &edges_[n1.nid_][n2.nid_];
  	internal_edge* this_edge_2 = &edges_[n2.nid_][n1.nid_];
  	// Remove from edges_direct
  	for (auto it = edges_direct_.begin(); it != edges_direct_.end(); ++it) {
  		if (*it == this_edge) {
  			edges_direct_.erase(it);
  			break;
  		}
  		if (*it == this_edge_2) {
  			edges_direct_.erase(it);
  			break;
  		}
  	}
    // Remove from edges_
  	edges_[n1.nid_].erase(n2.nid_);
  	edges_[n2.nid_].erase(n1.nid_);
  	return 1;
  }

  /** Remove an edge in the graph. Returns 1 if successful.
   * @param[in] e The edge to delete.
   * @post Removes the edge from the graph.
   * @post num_edges() will decrease by 1.
   * @post All edge iterators and incident iterators will be considered invalidated.
   * Complexity: O(num_edges()) amortized operations.
   */
  size_type remove_edge(const Edge& e) {
  	// Get pointer to internal_edge
  	internal_edge* this_edge = &edges_[e.nid1_][e.nid2_];
  	internal_edge* this_edge_2 = &edges_[e.nid2_][e.nid1_];
  	// Remove from edges_direct
  	for (auto it = edges_direct_.begin(); it != edges_direct_.end(); ++it) {
  		if (*it == this_edge) {
  			edges_direct_.erase(it);
  			break;
  		}
  		if (*it == this_edge_2) {
  			edges_direct_.erase(it);
  			break;
  		}
  	}
  	// Remove from edges_
  	edges_[e.nid1_].erase(e.nid2_);
  	edges_[e.nid2_].erase(e.nid1_);
  	return 1;
  }
  
  /** Remove an edge in the graph. Returns an edge iterator pointing to the next edge.
   * @param[in] e_it An edge iterator pointing to the edge to remove.
   * @post num_edges() will decrease by 1.
   * @post All edge iterators (other than the one returned by this) and all incident iterators
   *       will be considered invalidated.
   * Complexity: O(num_edges()) amortized operations (because erase is linear in size of
   *             number of elements after erased element).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
  	Edge e = *e_it;
  	// Remove from edges_
  	edges_[e.nid1_].erase(e.nid2_);
  	edges_[e.nid2_].erase(e.nid1_);
  	// Remove from edges_direct and return edge now in its current position
  	edges_direct_.erase(edges_direct_.begin() + e_it.edge_index);
  	return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
  	nodes_.clear();
  	edges_.clear();
  	edges_direct_.clear();
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
    value_type operator*() const {
      return mygraph->node(node_index);
    }

    NodeIterator& operator++() {
      if (node_index + 1 <= mygraph->size())
      	node_index++;
      return *this;
    }

    bool operator==(const NodeIterator& ni) const {
      if ((mygraph == ni.mygraph) && (node_index == ni.node_index))
        return true;
      return false;
    }

    bool operator!=(const NodeIterator& ni) const {
      return !(*(this) == ni);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph* mygraph;
    unsigned int node_index;

    // Private constructor that can be accessed by the Graph class.
    NodeIterator(const Graph* g, unsigned int ni) {
    	mygraph = g;
    	node_index = ni;
    }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const {
  	return NodeIterator(this, 0);
  }

  node_iterator node_end() const {
  	return NodeIterator(this, size());
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
    value_type operator*() const {
      internal_edge currentedge = myiterator->second;
      Graph* g = const_cast<Graph*>(mygraph);
      value_type e = g->add_edge(currentedge.node1, currentedge.node2);
      return e;
    }

    IncidentIterator& operator++() {
      myiterator = ++myiterator;
      return *this;
    }

    bool operator==(const IncidentIterator& ii) const {
      if ((mygraph == ii.mygraph) && (myiterator == ii.myiterator))
      	return true;
      return false;
    }

    bool operator!=(const IncidentIterator& ii) const {
      return!(*(this) == ii);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    const Graph* mygraph;
    typename std::map<Graph::size_type, Graph::internal_edge>::iterator myiterator;

    //Private constructor that can be accessed by the Node class.
    IncidentIterator(const Graph* g, typename std::map<Graph::size_type, Graph::internal_edge>::iterator it) {
      mygraph = g;
      myiterator = it;
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    value_type operator*() const {
      // This makes use of the add_edge behavior that returns existing edge
      Graph::internal_edge ie = *(mygraph->edges_direct_[edge_index]);
      Node n1 = ie.node1;
      Node n2 = ie.node2;
      Graph* g = const_cast<Graph*>(mygraph);
      value_type e = g->add_edge(n1, n2);
      return e;
    }

    EdgeIterator& operator++() {
      if (edge_index +1 <= mygraph->num_edges())
      	edge_index++;
      return *this;
    }

    bool operator==(const EdgeIterator& ei) const {
      if ((mygraph == ei.mygraph) && (edge_index == ei.edge_index))
      	return true;
      return false;
    }

    bool operator!=(const EdgeIterator& ei) const {
      return !(*(this) == ei);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* mygraph;
    unsigned int edge_index;

    // Private constructor that can be accessed by the Graph class
    EdgeIterator(const Graph* g, unsigned int ei) {
      mygraph = g;
      edge_index = ei;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const {
  	return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const {
  	return EdgeIterator(this, num_edges());
  }

 private:
  struct internal_node {
  	Point location; // location of a node
  	size_type nid; // unique identifier of a node
  	node_value_type nvalue; // node value
  };

  std::vector<internal_node> nodes_; // STL container for nodes
  size_type next_nid_; // counter for creating unique nodes

  struct internal_edge {
  	Node node1; // unique identifier of node1
  	Node node2; // unique identifier of node2
  	edge_value_type evalue; // edge value
  };

  std::map<size_type, std::map<size_type, internal_edge>> edges_;
  std::vector<internal_edge*> edges_direct_;

};

#endif // CME212_GRAPH_HPP
