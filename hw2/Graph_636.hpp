#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <type_traits>
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
template <typename V = int, typename E = double>
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

  typedef E edge_value_type;

  typedef std::unordered_map<size_type, size_type>::const_iterator \
          incidentitertype;

  typedef std::vector<size_type>::const_iterator pitertype;
  typedef std::vector<size_type>::const_iterator edgeitertype;


  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
    // Declare a graph with no nodes or edges with empty storage
    // structures
    : nodeinfos_(), edgeinfos_(), num_nodes_(0), \
      num_edges_(0), edgemap_() {
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

    /** HW2 modifiable node position */
    Point& position() {
      return const_cast<Point&>(const_cast<const Node*>(this)->position());
    }

    /** Return this node's position. */
    const Point& position() const {
      assert(valid());
      return fetchpoint();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->nodeinfos_.at(node_id_).idx_;
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
      return graph_->nodeinfos_.at(node_id_).v_;
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
      if (this->node_id_ < n.node_id_ || this->graph_ < n.graph_) {return true;}
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
    // INTERNAL ID REPRESENTED HERE
    size_type node_id_;

    /** Private Constructor */
    // why enter as const and use const cast here??
    Node(const Graph* newgraph, size_type new_node_id)
      : graph_(const_cast<Graph*>(newgraph)), node_id_(new_node_id) {
    } 

    // Helper method to return the correct Point
    // associated with the Node
    const Point& fetchpoint() const {
      return graph_->nodeinfos_.at(node_id_).p_;
    };

    // check if valid node in O(1)
    bool valid() const {
      return node_id_>= 0 && node_id_< graph_->nodeinfos_.size() // uid in range. 
            && graph_->nodeinfos_[node_id_].idx_ < graph_->nodei2u_.size() // idx in range. 
            && graph_->nodei2u_[graph_->nodeinfos_[node_id_].idx_] == node_id_; // uid in sync.
    }
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
    // Add to node info
    nodeinfo newnode(position, v, num_nodes());
    nodeinfos_.push_back(newnode);

    // add on the next internal id
    nodei2u_.push_back(nodeinfos_.size() - 1);

    // create a new empty map to store edges connecting to this node
    std::unordered_map <size_type, size_type> newmap;
    edgemap_.push_back(newmap);

    //update the number of nodes
    this->num_nodes_ += 1;

    // return a node with the correct internal id
    return Node(this, nodeinfos_.size() - 1);    
  }

  /**
   * @brief Removes Node and all incident edges connected to the node
   * @pre n is a Node of the graph (gets checked anyway)
   * @post num_nodes_ = num_nodes_ - 1 with a success
   * @post has_node(n) == false
   * @param n Node that is to be removed

   * @return size_type 1 if success, 0 if failure
   * Invalidates node with user facing idx same as _n_
   * (new graph.node(_n_.index()) != old graph.node(_n_.index()))
   * Complexity: O(1) to check node and remove it, num_incident_edges*O(1)
   * to remove all incident edges, so O(num_incident_edges) total.
   */
  size_type remove_node(const Node& n) {
  
    //remove node
    if (has_node(n) == true) {
      // remove all incident edges first
      std::vector<Edge> toremove;
      // must add to another vector since incident iterator is const
      for (IncidentIterator it = n.edge_begin(); it != n.edge_end(); ++it) {
        toremove.push_back(*it);
      }
      // now remove edges here
      for (typename std::vector<Edge>::iterator eit = toremove.begin(); \
          eit!= toremove.end(); ++eit) {
        remove_edge(*eit);
      }

      // find the external id and move the last valid node's info
      // to this external id in nodei2u_
      size_type node_external_id = nodeinfos_.at(n.node_id_).idx_;   
      nodei2u_.at(node_external_id) = std::move(nodei2u_.back());
      //update external id to match
      nodeinfos_.at(nodei2u_.at(node_external_id)).idx_ = node_external_id;
      //delete the last element
      nodei2u_.erase(nodei2u_.end() - 1);
      num_nodes_ -= 1;
      return 1;
    }
    return 0;
  }

  /**
   * @brief Removes Node at current iterator spot
   * @pre *n_it is a Node of the graph (gets checked anyway)
   * @post num_nodes_ = num_nodes_ - 1 with a success
   * @post has_node(*n_it) == false
   * @param n_it Iterator pointing to node that is to be removed

   * @return @a n_it Iterator pointing to same node if failed or
   * pointing to end of @a nodei2u_ container since last element erased
   * Invalidates node with user facing idx same as _n_
   * (new graph.node(_n_.index()) != old graph.node(_n_.index()))
   * Also invalidates node with user facing idx same as _n_
   * (new graph.node(_n_.index()) != old graph.node(_n_.index()))
   * Complexity: O(1) to check node and remove it, num_incident_edges*O()
   * to remove all incident edges, 
   */
  node_iterator remove_node(node_iterator n_it) {
    Node n = *n_it;
    if (has_node(n) == true) {
      // remove all incident edges first
      std::vector<Edge> toremove;
      // must add to another vector since incident iterator is const
      for (IncidentIterator it = n.edge_begin(); it != n.edge_end(); ++it) {
        toremove.push_back(*it);
      }
      // now remove edges here
      for (typename std::vector<Edge>::iterator eit = toremove.begin(); \
          eit!= toremove.end(); ++eit) {
        remove_edge(*eit);
      }

      // find the external id and move the last valid node's info
      // to this external id in nodei2u_
      size_type node_external_id = nodeinfos_.at(n.node_id_).idx_;   
      nodei2u_.at(node_external_id) = std::move(nodei2u_.back());
      //update external id of moved node to match
      nodeinfos_.at(nodei2u_.at(node_external_id)).idx_ = node_external_id;
      //delete the last element from i2u_ and assign to temp
      std::vector<size_type>::const_iterator tmp = \
        nodei2u_.erase(nodei2u_.end() - 1);
      num_nodes_ -= 1;
      return NodeIterator(this, tmp);
    }
    return n_it;
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
    return Node(this, nodei2u_.at(i));
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
      // check node is valid
      assert(Node(graph_, node1_).valid());
      return Node(graph_, node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // check node is valid
      assert(Node(graph_, node2_).valid());
      return Node(graph_, node2_);
    }

    /** HW2 helper func 
    * @brief Returns Euclidean length of the Edge
    */
    double length() const {
      assert(valid());
      return norm_2(node2().position() - node1().position());
    }

    // HW2
    /** 
     * @brief Returns a non-const reference to the attribute of a node
     * can be used to modify what the attribute is
     */
    edge_value_type& value(){
      return const_cast<edge_value_type&>(const_cast<const Edge*>(this)\
                                          ->value());
    }

    /**
     * @brief Const member function returning const reference to
     * the attribute associated with a node
     * Complexity: O(1)
     */
    const edge_value_type& value() const{
      assert(valid());
      return graph_->edgeinfos_.at(edge_id_).v_;
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
      // check if this is a difference graph entirely as well
      if (this->edge_id_ < e.edge_id_ || this->graph_ < e.graph_) {
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

    // Internal edge id
    size_type edge_id_;

    // Internal node idxs
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

    // check if valid edge w/ valid node endpoints
    bool valid() const {
      return edge_id_>= 0 && edge_id_< graph_->edgeinfos_.size() // uid in range. 
            && graph_->edgeinfos_[edge_id_].idx_ < graph_->edgei2u_.size() // idx in range. 
            && graph_->edgei2u_[graph_->edgeinfos_[edge_id_].idx_] == edge_id_ // uid in sync.
            && Node(graph_, node1_).valid() // check node 1 is valid
            && Node(graph_, node2_).valid(); // check node 2 is valid
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
    size_type internal_id = edgei2u_.at(i);
    size_type n1 = (edgeinfos_.at(internal_id).endpoints_).first;
    size_type n2 = (edgeinfos_.at(internal_id).endpoints_).second;
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
  Edge add_edge(const Node& a, const Node& b, \
                const edge_value_type& v = edge_value_type()) {
    if (has_edge(a,b) == true) {
      // we know there is an edge, so return the edge index
      size_type edge_internal_id = edgemap_.at(a.node_id_)\
                                   .find(b.node_id_) -> second;
      return Edge(this, edge_internal_id, \
                  a.node_id_, b.node_id_);
      }

    else {
      // add edge info to storage
      // use internal ids to persist edge information
      edgeinfo newedge(PointPair(a.node_id_, b.node_id_), v, num_edges());
      edgeinfos_.push_back(newedge);
      // map new external ID to latest internal ID
      edgei2u_.push_back(edgeinfos_.size() - 1);
      // undirected edge: the edge goes both ways so must add
      // to both a and b's maps of nodes connected to it
      // mark where the information about the connected points are
      // use internal node ids to persist information
      edgemap_.at(a.node_id_)[b.node_id_] = edgeinfos_.size() - 1;
      edgemap_.at(b.node_id_)[a.node_id_] = edgeinfos_.size() - 1;

      //update the number of edges
      num_edges_ += 1;
      return Edge(this, edgeinfos_.size() - 1, a.node_id_, b.node_id_);
    }
  }

  /**
   * @brief Remove an edge based on its two endpoints
   * @pre There is an edge with endpoints a, b
   * @post num_edges_ = num_edges_ - 1
   * @post has_edge(a, b) == false
   * @param a Endpoint node 1
   * @param b Endpoint node 2
   * @return size_type 1 for a successful removal, 0 for a failed removal.
   * Complexity: O(1) to find the internal and external ids, O(1) to move
   * and erase the last element each and O(1) to locate in map and erase 
   * from edgemap, so O(1) total.
   */
  size_type remove_edge(const Node& a, const Node& b) {
    // find the edge ID
    if (has_edge(a,b) == true) {
      // we know there is an edge, so return the edge index
      size_type edge_internal_id = edgemap_.at(a.node_id_)\
                                   .at(b.node_id_);
      assert(Edge(this, edge_internal_id, a.node_id_, b.node_id_).valid());
      size_type edge_external_id = edgeinfos_.at(edge_internal_id).idx_;

      edgei2u_.at(edge_external_id) = std::move(edgei2u_.back());
      //update external id to match
      edgeinfos_.at(edgei2u_.at(edge_external_id)).idx_ = edge_external_id;

      //delete the last element
      edgei2u_.erase(edgei2u_.end() - 1);
            
      // delete from edgemap_
      edgemap_.at(a.node_id_).erase(b.node_id_);
      edgemap_.at(b.node_id_).erase(a.node_id_);
      num_edges_ -= 1;
      return 1;
    }
    return 0;
  }
  /**
   * @brief Remove an edge given
   * @pre @a e is a valid edge in the graph
   * @post num_edges_ = num_edges_ - 1
   * @post has_edge(e.node1_, e.node2_) == false
   * @param e Valid Edge to delete
   * @return size_type 1 for a successful removal, 0 for a failed removal.
   * Complexity: O(1) to find the internal and external ids, O(1) to move
   * and erase the last element each and O(1) to locate in map and erase 
   * from edgemap, so O(1) total.
   */
  size_type remove_edge(const Edge& e){
    assert(e.valid());
    
    if (has_edge(Node(this, e.node1_), Node(this, e.node2_)) == true) {
      // find the external id and map it to last element
      size_type edge_external_id = edgeinfos_.at(e.edge_id_).idx_;
      edgei2u_.at(edge_external_id) = std::move(edgei2u_.back());
      //update external id of new element to match
      edgeinfos_.at(edgei2u_.at(edge_external_id)).idx_ = edge_external_id;

      //delete the last element
      edgei2u_.erase(edgei2u_.end() - 1);
            
      // delete from edgemap_
      edgemap_.at(e.node1_).erase(e.node2_);
      edgemap_.at(e.node2_).erase(e.node1_);
      num_edges_ -= 1;
      return 1;
    }
    return 0;
  }

  /**
   * @brief Remove the edge pointed to by this edge iterator
   * @pre @a *e_it is a valid edge in the graph
   * @post num_edges_ = num_edges_ - 1
   * @post has_edge((*e_it).node1_, (*e_it).node2_) == false
   * @param e_it iterator pointing to valid Edge to delete
   * @return edge_iterator pointing to either the same edge if failed or
   * edge_iterator pointing to end of edgei2u_ container.
   * Complexity: O(1) to find the internal and external ids, O(1) to move
   * and erase the last element each and O(1) to locate in map and erase 
   * from edgemap, so O(1) total.
   */
  edge_iterator remove_edge(edge_iterator e_it){
    try {
      Edge e = *e_it;
      assert(e.valid());
      if (has_edge(Node(this, e.node1_), Node(this, e.node2_)) == true) {
        size_type edge_external_id = edgeinfos_.at(e.edge_id_).idx_;

        edgei2u_.at(edge_external_id) = std::move(edgei2u_.back());
        //update external id of new element to match
        edgeinfos_.at(edgei2u_.at(edge_external_id)).idx_ = edge_external_id;

        //delete the last element and assign result to tmp
        std::vector<size_type>::const_iterator tmp = \
          edgei2u_.erase(edgei2u_.end() - 1);
              
        // delete from edgemap_
        edgemap_.at(e.node1_).erase(e.node2_);
        edgemap_.at(e.node2_).erase(e.node1_);
        num_edges_ -= 1;
        return EdgeIterator(this, tmp);
      }
      return e_it;
    }
    catch(...){
      throw;
    }
  }
  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

    nodeinfos_.clear();
    edgeinfos_.clear();
    nodei2u_.clear();
    edgei2u_.clear();
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
      return Node(graph_, *point_iter_);
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
    pitertype it = nodei2u_.begin();
    return NodeIterator(this, it);
  }

  /**
   * @brief Creates a NodeIterator pointing to one past the last node
   */
  node_iterator node_end() const {
    pitertype it = nodei2u_.end();
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
      int n2idx = incident_iter_->first; // key of the map is the adj edge

      return Edge(graph_, edge_id, node_index_, n2idx);
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
      size_type edgeid = *edge_iter_;
      PointPair endpoints = graph_->edgeinfos_.at(edgeid).endpoints_; 
      return Edge(graph_, edgeid, \
                  endpoints.first, endpoints.second);
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
    edgeitertype it = edgei2u_.begin();
    return EdgeIterator(this, it);
  }

  /**
   * @brief Creates a EdgeIterator pointing past the end of the edge 
   * container
   */
  edge_iterator edge_end() const {
    edgeitertype it = edgei2u_.end();
    return EdgeIterator(this, it);
  }

  // Store node information in graph
  struct nodeinfo {
    Point p_;
    node_value_type v_;
    // User facing ID
    size_type idx_;    
    
    public:
      //constructor
      nodeinfo(Point p, node_value_type v, size_type idx): p_(p), v_(v), idx_(idx) {}
  };

  // Storage structure for edge info
  struct edgeinfo {
    PointPair endpoints_;
    edge_value_type v_;
    // User facing ID
    size_type idx_;

    public:
      //constructor
      edgeinfo(PointPair endpoints, edge_value_type v, size_type idx): \
               endpoints_(endpoints), v_(v), idx_(idx) {}
  };

 private:
  // Storing node data
  std::vector<nodeinfo> nodeinfos_;
  // Storing edge data
  std::vector<edgeinfo> edgeinfos_;

  // Map from user-facing idx to internal uid for nodes and edges
  std::vector<size_type> nodei2u_;
  std::vector<size_type> edgei2u_;
  // Number of Nodes
  size_type num_nodes_;
  // Number of Edges
  size_type num_edges_;

  // Each element at index i is a map where they keys are the INTERNAL node
  // idxs connecting to node (INTERNAL) i. The values are the corresponding
  // edge index.
  std::vector<std::unordered_map<size_type, size_type> > edgemap_;
};

#endif // CME212_GRAPH_HPP
