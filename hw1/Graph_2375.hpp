#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <utility>
#include <typeinfo>


#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
// default template type is int
template <typename V = int>
class Graph {
 private:

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

  /** Type of endpoints_ elements corresponding to node indices */
  using PointPair = std::pair<size_type, size_type>;

  using node_value_type = V;

  typedef std::unordered_map<size_type, size_type>::const_iterator incidentitertype;

  typedef std::vector<Point>::const_iterator pitertype;

  typedef std::vector<PointPair>::const_iterator edgeitertype;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    // Declare a graph with no nodes or edges with empty storage
    // structures
    : points_(), attributes_(), num_nodes_(0), \
      num_edges_(0), endpoints_(), edgemap_() {
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
    using size_type = unsigned;
    Node() {
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetchpoint();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return node_id_;
    }

    // reference: https://bit.ly/39Ttzyg    
    /** 
     * @brief Returns a non-const reference to the attribute of a node
     * can be used to modify what the attribute is
     */
    node_value_type& value(){
      return const_cast<node_value_type&>(const_cast<const Node*>(this)\
                                          ->value());
    }

    /**
     * @brief Const member function returning const reference to
     * the attribute associated with a node
     */
    const node_value_type& value() const{
      return graph_->attributes_.at(node_id_);
    }
    
    /**
     * @brief Gets the degree (number of edges) connected to
     * a node
     */
    size_type degree() const {
      return graph_->edgemap_.at(node_id_).size();
    }

    /**
     * @brief Creates a IncidentIterator pointing to the first edge
     * connected to a node. Note since edges are stored as a set for each node,
     * there is no guarantee this will be the same edge as the order they 
     * were inserted
     */
    incident_iterator edge_begin() const {
      incidentitertype it = graph_->edgemap_.at(node_id_).begin();
      return IncidentIterator(graph_, node_id_, it);
    }

    /**
     * @brief Creates a IncidentIterator pointing to the last edge
     * connected to a node. Note since edges are stored as a set for each node,
     * there is no guarantee this will be the same edge as the order they 
     * were inserted
     */
    incident_iterator edge_end() const {
      incidentitertype it = graph_->edgemap_.at(node_id_).end();
      return IncidentIterator(graph_, node_id_, it);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (graph_ == n.graph_ && this->node_id_ == n.node_id_) {
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
      if (this->node_id_ < n.node_id_) {return true;}
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to Graph container
    Graph* graph_;

    // This node's index number
    size_type node_id_;

    /** Private Constructor */
    // why enter as const and use const cast here??
    Node(const Graph* newgraph, size_type new_node_id)
      : graph_(const_cast<Graph*>(newgraph)), node_id_(new_node_id) {
    } 

    // Helper method to return the correct Point
    // associated with the Node
    const Point& fetchpoint() const {
      return graph_->points_.at(node_id_);
    };
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
  Node add_node(const Point& position,
                const node_value_type& v = node_value_type()) {
    // Add to the graph's storage of points
    points_.push_back(position);
    attributes_.push_back(v);
    // create a new empty map to store edges connecting to this node
    std::unordered_map <size_type, size_type> newmap;
    edgemap_.push_back(newmap);

    //update the number of nodes
    this->num_nodes_ += 1;
    return Node(this, num_nodes() - 1);    
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // check if graph the node points to is the same thing
    if (n.graph_ == this) {
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_, node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_, node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // check edges have the same set of nodes
      if ((node1() == e.node1() && node2() == e.node2()) \
        || (node2() == e.node1() && node1() == e.node2())) {
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
      if (this->edge_id_ < e.edge_id_) {
        return true;
      }
      else {
        return false;
      }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Graph this edge is associated with
    Graph* graph_;

    // Index of the Edge
    // This node's index number
    size_type edge_id_;

    // indices of the nodes
    size_type node1_;
    
    size_type node2_;

    // Private Constructor
    Edge(const Graph* newgraph,
         size_type new_edge_id,
         size_type n1,
         size_type n2)
      : graph_(const_cast<Graph*>(newgraph)), edge_id_(new_edge_id), \
      node1_(n1), node2_(n2) {
    } 
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
    size_type n1 = endpoints_.at(i).first;
    size_type n2 = endpoints_.at(i).second;
    return Edge(this, i, n1, n2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    // check if node a contains node b in its map of edges
    if (edgemap_.at(a.node_id_).find(b.node_id_) \
        != edgemap_.at(a.node_id_).end()) {
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
  Edge add_edge(const Node& a, const Node& b) {

    if (has_edge(a,b) == true) {
      // we know there is an edge, so return the edge index
      return Edge(this, edgemap_.at(a.node_id_).find(b.node_id_) -> second, \
                  a.node_id_, b.node_id_);
      }

    else {
      // add the pair of Nodes to our NodePair storage
      endpoints_.push_back(PointPair(a.node_id_, b.node_id_));

      // undirected edge: the edge goes both ways so must add
      // to both a and b's maps of nodes connected to it
      edgemap_.at(a.node_id_)[b.node_id_] = num_edges();
      edgemap_.at(b.node_id_)[a.node_id_] = num_edges();

      //update the number of edges
      this->num_edges_ += 1;
      return Edge(this, num_edges() - 1, a.node_id_, b.node_id_);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

    points_.clear();
    endpoints_.clear();
    edgemap_.clear();
    num_nodes_ = 0;
    num_edges_ = 0;
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
    }

    /** 
     * @brief This dereferences the Node iterator and returns a node
     * @pre The Node iterator does not point to node_end()
     */
    Node operator*() const {
      // find the correct node index by subtracting the beginning
      return Node(graph_, point_iter_ - graph_->points_.begin());
    }

    /**
     * @brief Increment the NodeIterator
     * @pre the NodeIterator doesn't point to node_end() already
     */
    NodeIterator& operator++() {
        ++point_iter_;
        return (*this);
    }

    /**
     * @brief Tests for equality between two NodeIterators
     * @param[in] other_it a valid NodeIterator 
     */
    bool operator==(const NodeIterator& other_it) const {
      return (graph_ == other_it.graph_) && \
       (point_iter_ == other_it.point_iter_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    pitertype point_iter_;

    // Private constructor
    NodeIterator(const Graph* graphpointer, pitertype& pi):
      graph_(const_cast<Graph*> (graphpointer)),
      point_iter_(pi) {
      }
  };

  /**
   * @brief Creates a NodeIterator pointing to the node indexed first
   */
  node_iterator node_begin() const {
    pitertype it = points_.begin();
    return NodeIterator(this, it);
  }

  /**
   * @brief Creates a NodeIterator pointing to one past the last node
   */
  node_iterator node_end() const {
    pitertype it = points_.end();
    return NodeIterator(this, it);
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
    }

    /**
     * @brief Dereferences an IncidentIterator
     * @pre The IncidentIterator is valid and doesn't point past the end
     * of the container
     * @return an incident Edge to the current node,
     * with node1() == current node and node2() == adjacent node
     */
    Edge operator*() const {
      size_type edge_id = (incident_iter_->second);
      int n1idx = incident_iter_->first; // key of the map is the edge

      return Edge(graph_, edge_id, node_index_, n1idx);
    }

    /**
     * @brief Increments the IncidentIterator
     */
    IncidentIterator& operator++() {
      ++incident_iter_;
      return (*this);
    }

    /**
     * @brief Tests for equality between two IncidentIterators
     * @pre _other_it_ is a valid IncidentIterator
     */
    bool operator==(const IncidentIterator& other_it) const{
      return (graph_ == other_it.graph_) && \
          (node_index_ == other_it.node_index_) && \
          (incident_iter_ == other_it.incident_iter_);
    }
    
   private:
    friend class Graph;
    Graph* graph_;
    size_type node_index_;
    incidentitertype incident_iter_;

    // Private constructor
    IncidentIterator(const Graph* graphpointer,
                     size_type id,
                     incidentitertype& ei):
      graph_(const_cast<Graph*> (graphpointer)),
      node_index_(id),
      incident_iter_(ei) {
      }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator>{
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
     * @brief Dereferences an EdgeIterator
     * @pre The EdgeIterator is valid and doesn't point past the end
     * of the container
     * @return An Edge with nodes in no particular order
     */
    Edge operator*() const {
      return Edge(graph_, edge_iter_ - graph_->endpoints_.begin(), \
                  edge_iter_->first, edge_iter_->second);
    }

    /**
     * @brief Increments the EdgeIterator
     */
    EdgeIterator& operator++() {
      ++edge_iter_;
      return (*this);
    }

    /**
     * @brief Tests for equality between EdgeIterators
     * @pre _other_it_ is a valid EdgeIterator
     */
    bool operator==(const EdgeIterator& other_it) const{
      return (graph_ == other_it.graph_) && \
       (edge_iter_ == other_it.edge_iter_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    edgeitertype edge_iter_;

    // Private constructor
    EdgeIterator(const Graph* graphpointer, edgeitertype& ei):
      graph_(const_cast<Graph*> (graphpointer)),
      edge_iter_(ei) {
      }
  };

  /**
   * @brief Creates a EdgeIterator pointing to the edge indexed first
   */
  edge_iterator edge_begin() const {
    edgeitertype it = endpoints_.begin();
    return EdgeIterator(this, it);
  }

  /**
   * @brief Creates a EdgeIterator pointing past the end of the edge 
   * container
   */
  edge_iterator edge_end() const {
    edgeitertype it = endpoints_.end();
    return EdgeIterator(this, it);
  }

 private:
  // Storing the points of the Nodes
  std::vector<Point> points_;
  // Storing the attributes of Nodes
  std::vector<node_value_type> attributes_;
  // Number of Nodes
  size_type num_nodes_;
  // Number of Edges
  size_type num_edges_;
  // Point Pairs corresponding to each edge index
  // ex: index 0 contains the pair of indices in points_
  // corresponding the nodes for edge 0
  std::vector<PointPair> endpoints_;
  // Each element at index i is a map where they keys are the nodes
  // connecting to node i. The values are the corresponding edge index
  // Makes finding edges by index easier at the cost of more storage
  std::vector<std::unordered_map<size_type, size_type> > edgemap_;
};

#endif // CME212_GRAPH_HPP
