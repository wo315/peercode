#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <list> // linked list from STL 

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph
{
private:

  /** Predeclaration of Internal Node type. */
  struct internal_node_;
  /** Synonym for InternalNode_ (following STL conventions) **/
  using internal_node_type = internal_node_;

  /** Predeclaration of Internal edge type */
  struct internal_edge_;
  /** Synonym for internal_edge_ (following STL conventions) **/
  using internal_edge_type = internal_edge_;

public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Synonym for value of a node */ 
  using node_value_type = V;

  /** Synonym for value of an edge */
  using edge_value_type = E;

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

  /** Alias for the type of iterator used to traverse the std::vector container for nodes */ 
  using stl_node_iterator = typename std::vector<unsigned>::const_iterator; 

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Alias for the type of iterator used to traverse the std::vector container for edges */
  using stl_edge_iterator = typename std::vector<unsigned>::const_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Alias for the type of iterator used to traverse the std::vector container for incident edges */
  using stl_incident_iterator = typename std::vector<unsigned>::const_iterator;

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
   * This class is lightweight and thus sizeof(Graph::Node) <= 16 bytes. 
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
      my_graph_ = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().pos_;
    }

    /** Return this node's position (non-const to allow for modification). */
    Point &position()
    {
      return fetch().pos_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return fetch().idx_;
    }

    /** Return this node's degree. */ 
    size_type degree() const{
      return fetch().incident_edges_.size(); 
    }

    /** Return this node's value. */
    node_value_type& value() {
      return fetch().value_; 
    };

    /** Return this node's value (overloaded const version of the above method). */
    const node_value_type& value() const {
      return fetch().value_;
    };

    /** return an iterator pointing at the start of the incident edge list. */
    incident_iterator edge_begin() const
    {
      return IncidentIterator(fetch().incident_edges_.begin(), fetch().my_graph_, index()); 
    }

    /** return an iterator pointing at the end (aka one past the last edge) of the incident edge list. */
    incident_iterator edge_end() const
    {
      return IncidentIterator(fetch().incident_edges_.end(), fetch().my_graph_, index());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((this->my_graph_ == n.my_graph_) & (this->index() == n.index()));
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
      // if from same graph, order by node index. Otherwise, order by 
      // address of each graph. Inspired by HW0 peercode #276.
      return (my_graph_ != n.my_graph_) ? (my_graph_ < n.my_graph_) : (this->index() < n.index());
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the parent graph. 8 bytes of memory. 
    Graph* my_graph_;

    // node's internal ID in the list of all nodes that ever existed for this graph. 
    // 4 bytes of memory. 
    size_type uid_; 

    /** Private Constructor */
    Node(const Graph *graph, size_type uid)
        : my_graph_(const_cast<Graph*>(graph)), uid_(uid)
    {
    }

    // Helper method to return the corresponding internal node. 
    // Inspired by peercode for HW0
    internal_node_type& fetch() const
    {
      return my_graph_->nodes_[uid_];
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return active_nodes_.size(); 
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
  Node add_node(const Point& position, const node_value_type& value=node_value_type()) {

    // compute new node uid and idx 
    size_type num_total_nodes = nodes_.size(); 
    size_type num_active_nodes = num_nodes(); 
    size_type new_node_uid = num_total_nodes;
    size_type new_node_idx = num_active_nodes;

    // update node containers
    nodes_.push_back(internal_node_(this, new_node_idx, new_node_uid, position, value));
    active_nodes_.push_back(new_node_uid); 

    // Returns a proxy to the new node
    return Node(this, new_node_uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return ((this == n.my_graph_) & (n.index() < n.my_graph_->size())); 
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes()); 
    return Node(this, active_nodes_[i]); // return a proxy object for element i, as in proxy_example.cpp 
  }

  /** Remove a node from the graph, returning 1 if the node was removed
   *  succesfully and 0 otherwise (e.g. if the node did not exist).
   * @pre @a n is a valid node of this graph
   * @return int 1 if node removal is succesful, 0 otherwise. 
   * @post has_node( @a n ) == false
   * @post If old has_node( @a n ), new num_nodes() == old num_nodes() - 1.
   *       Else,                        new num_nodes() == old num_nodes().
   *
   * Can invalidate node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). 
   * 
   * Can invalidate node iterators (type node_iterator) -- this is because
   * remove_node() alters the underlying std::vector (active_nodes_) that
   * is backing each node iterator. 
   * 
   * Can invalidate edge indices since remove_edge() is called in this method.
   * Thus, old edge(@a i) might not equal new edge(@a i). 
   * 
   * Can invalidate edge iterators (type edge_iterator) -- this is because
   * remove_edge() is called, which alters the underlying std::vector (active_edges_) 
   * that is backing each edge iterator. 
   *
   * Complexity: O(n.degree() + swapped_node.degree()) in general where @a swapped_node is the
   *             node utilized in the std::iter_swap iteration below. This is thus
   *             O(1) for a sparse graph. 
   */
  size_type remove_node(const Node& n){
    // check if the provided node is present in this graph
    if (has_node(n)) 
    {      

      // remove the edges incident to this node
      while (n.degree()>0){ // O(n.degree()), equal to O(1) for a sparse graph 
        auto ei = n.edge_begin();
        remove_edge(*ei); 
      }

      // note this node's uid and idx
      size_type uid = n.uid_; 
      size_type idx = n.index();

      // remove the node's uid from the vector of active node
      // this can be done in O(1) time by swapping places
      // with the element at the end of the container
      // and popping off the removed-node's uid
      bool no_swap = (active_nodes_.begin() + idx == active_nodes_.end() - 1) ? true : false; 
      std::iter_swap(active_nodes_.begin() + idx, active_nodes_.end() - 1); // O(1) 
      active_nodes_.pop_back(); // O(1) 

      // update the internally stored value of idx_ for the node
      // that was previously at the end of the container
      if (!no_swap) nodes_[active_nodes_[idx]].idx_ = idx;

      // update the internally stored value of node_a and node_b
      // for edges incident to the node previously at the end of the
      // active_nodes_ container (since its idx changed)
      if (!no_swap){
        node_type swapped_node = node(idx);
        for (auto ei = swapped_node.edge_begin(); ei != swapped_node.edge_end(); ++ei){
          // total cost of for loop is O(swapped_node.degree()), equal to O(1) for a sparse graph 
          edge_type e = *ei;

          size_type node_a = e.node1().index(); 
          size_type node_b = e.node2().index(); 
          e.fetch().node_a = std::min(node_a, node_b);
          e.fetch().node_b = std::max(node_a, node_b);
        }
      }

      return 1;
    }
    else
    {
      // no node to remove, return 0
      return 0;
    }
  }

  /** Remove a node from the graph, returning a valid node iterator
   * @pre @a n_it is a valid node iterator of this graph
   * @return a valid node iterator for this graph 
   * @post has_node( @a *n_it ) == false
   * @post If old has_node( @a *n_it ), new num_nodes() == old num_nodes() - 1.
   *       Else,                        new num_nodes() == old num_nodes().
   *
   * Can invalidate node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). 
   * 
   * Can invalidate node iterators (type node_iterator) -- this is because
   * remove_node() alters the underlying std::vector (active_nodes_) that
   * is backing each node iterator. Certainly invalidates n_it by construction. 
   * 
   * Can invalidate edge indices since remove_edge() is called in this method.
   * Thus, old edge(@a i) might not equal new edge(@a i). 
   * 
   * Can invalidate edge iterators (type edge_iterator) -- this is because
   * remove_edge() is called, which alters the underlying std::vector (active_edges_) 
   * that is backing each edge iterator. 
   *
   * Complexity: O(n.degree() + swapped_node.degree()) in general where @a swapped_node is the
   *             node utilized in the std::iter_swap iteration below. This is thus
   *             O(1) for a sparse graph. 
   */
  node_iterator remove_node(node_iterator n_it)
  {
    // assert that this iterator is not graph.node_end()
    // otherwise, there is nothing to remove and the increment
    // operator below will cause a segfault
    assert(n_it != this->node_end() && "invalid node iterator provided");

    // with the way we have implemented remove_node(), the element
    // at this current node's index is one that has not yet
    // been visited. so we return an iterator starting here.
    node_type n = *n_it;
    size_type node_idx = n.index();
    remove_node(n);
    return NodeIterator(this, active_nodes_.begin() + node_idx);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   * This class is lightweight and thus sizeof(Graph::Edge) <= 32 bytes. 
   */
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      my_graph_ = nullptr;
      uid_ = 0;
    }

    /** Return this edge's value. */
    edge_value_type &value()
    {
      return fetch().value_;
    };

    /** Return this edge's value (overloaded const version of the above method). */
    const edge_value_type &value() const
    {
      return fetch().value_;
    };

    /** Return this edge's index, a number in the range [0, num_edges). */
    size_type index() const
    {
      return fetch().idx_;
    }

    /** Compute the length of an edge */
    double length() const
    {
      // length = distance between the nodes of the edge
      return norm(node1().position() - node2().position());
    }

    /** Return node @ a for an edge created via add_edge(Node a, Node b) */
    Node node1() const {
      // determine uids of the two nodes for this edge
      size_type node_a_idx = fetch().node_a;
      size_type node_a_uid = my_graph_->active_nodes_[node_a_idx];
      size_type node_b_idx = fetch().node_b;
      size_type node_b_uid = my_graph_->active_nodes_[node_b_idx];

      // return a proxy node. use reverse_order_ flag to make sure we return the correct node
      return reverse_order_ ? Node(my_graph_, node_b_uid) : Node(my_graph_, node_a_uid);
    }

    /** Return node @ b for an edge created via add_edge(Node a, Node b) */
    Node node2() const {
      // determine uids of the two nodes for this edge
      size_type node_a_idx = fetch().node_a;
      size_type node_a_uid = my_graph_->active_nodes_[node_a_idx];
      size_type node_b_idx = fetch().node_b;
      size_type node_b_uid = my_graph_->active_nodes_[node_b_idx];

      // return a proxy node. use reverse_order_ flag to make sure we return the correct node
      return reverse_order_ ? Node(my_graph_, node_a_uid) : Node(my_graph_, node_b_uid);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((this->my_graph_ == e.my_graph_) & (this->index() == e.index()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // if from same graph, order by edge index. Otherwise, order by
      // address of each graph. Inspired by HW0 peercode #276.
      return (my_graph_ != e.my_graph_) ? (my_graph_ < e.my_graph_) : (this->index() < e.index());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer back to the parent graph (8 bytes of memory) 
    Graph *my_graph_;

    // edge's internal ID in the list of all nodes that ever existed for this graph.
    // 4 bytes of memory.
    size_type uid_;

    // whether to return nodes in ascending or descending order of their index
    // we need this so that we can satisfy e.node1() == @a a and e.node2() == @a b
    // regardless of the relative indices of the input nodes in add_edge() 
    bool reverse_order_; 

    /** Private Constructor */
    Edge(const Graph *graph, size_type uid, bool reverse_order = false)
        : my_graph_(const_cast<Graph*>(graph)), uid_(uid), reverse_order_(reverse_order)
    {
    }

    // Helper method to return the corresponding internal edge. 
    // Inspired by peercode for HW0 
    internal_edge_type& fetch() const
    {
      return my_graph_->edges_[uid_];
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return active_edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, active_edges_[i]); // return a proxy object for element i, as in proxy_example.cpp
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(a.degree()) in general, which is O(1) for a sparse graph. 
   */
  bool has_edge(const Node& a, const Node& b) const {

    // figure out which of a and b has the smaller and larger index
    size_type small_idx = std::min(a.index(), b.index());
    size_type large_idx = std::max(a.index(), b.index());

    // loop through the edges of node a checking to see if the candidate edge
    // connects node a and node b. 
    size_type vec_a_size =  a.fetch().incident_edges_.size(); 
    for (size_type i = 0; i < vec_a_size; i++) // O(a.degree())
    {
      size_type current_edge_uid = a.fetch().incident_edges_[i];
      if (edges_[current_edge_uid].node_a == small_idx && edges_[current_edge_uid].node_b == large_idx)
      {
        return true;
      }
    }

    // If we reach this point, then an appropriate edge was not found 
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
  // Note: revised after HW0 using ideas from peercode to meet runtime complexity 
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& value=edge_value_type()) {

    // check that a and b are present in the graph and that they are not the same node
    assert((this->has_node(a) && this->has_node(b)) && !(a == b));

    // take note of what the new edges index and uid would be
    size_type new_edge_idx = num_edges(); 
    size_type new_edge_uid = active_edges_.size(); 

    // set flag that keeps track of whether or not node a has the smaller index
    // syntax for ternary operator is <condition> ? <true-case> : <false-case> 
    bool reverse_flag = (a.index() < b.index()) ? false : true;

    // loop through the edges of node a checking to see if the candidate edge
    // connects node a and node b. If so, return proxy edge. 
    if (has_edge(a, b)){
      return Edge(this, new_edge_uid, reverse_flag);
    }

    // If we reach this point, then the edge doesn't exist. 
    // In this case, we should push_back the new edge
    // and update the incident_edges of node a and node b
    // we also need to update the number of edges
    edges_.push_back(internal_edge_(this, new_edge_idx, new_edge_uid, a.index(), b.index(),value));
    active_edges_.push_back(new_edge_uid); 
    a.fetch().incident_edges_.push_back(new_edge_uid);
    b.fetch().incident_edges_.push_back(new_edge_uid);

    // Return a proxy to the new edge
    return Edge(this, new_edge_uid, reverse_flag);
  }

  /** Remove an edge from the graph, returning 1 if the edge was removed
   *  succesfully and 0 otherwise (e.g. if the edge did not exist).
   * @pre @a n1 and @a n2 are distinct valid nodes of this graph
   * @return int 1 if edge removal is succesful, 0 otherwise. 
   * @post has_edge(@a n1, @a n2) == false
   * @post If old has_edge(@a n1, @a n2), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). 
   * 
   * Can invalidate edge iterators (type edge_iterator) -- this is because
   * remove_edge() alters the underlying std::vector (active_edges_) that
   * is backing the edge iterator. 
   *
   * Complexity: O(n1.degree() + n2.degree()) in general, which is O(1)
   *             for a sparse graph. 
   */
  size_type remove_edge(const Node& n1, const Node& n2){

    // check if the provided edge is present in this graph 
    // i.e. first make sure both of the nodes belong to this graph
    // and then see if there is an edge connecting them 
    if (has_node(n1) && has_node(n2) && has_edge(n1, n2))
    { // this check is O(1) for a sparse graph, O(n1.degree()) in general 

      // if this condition is met, remove the edge, return 1

      // loop over edges of n1 to find idx of this edge
      size_type edge_idx = 0; // initialize to remove compiler warning
      for (auto ei = n1.edge_begin(); ei != n1.edge_end(); ++ei){ 
        // loop is O(n1.degree()) in general, which is O(1) for a sparse graph 
        
        auto e = *ei; 

        // if n2 is the other node of this edge, store the edge's index
        // note that we are guranteed to enter the below if statement
        // exactly once see we have verified Edge(n1, n2) exists 
        if (e.node2() == n2) edge_idx = e.index(); 
      }
      size_type edge_uid = active_edges_[edge_idx]; 

      // remove the edge's uid from the vector of active edges
      // this can be done in O(1) time by swapping places
      // with the element at the end of the container
      // and popping off the removed-edge's uid
      bool no_swap = (active_edges_.begin() + edge_idx == active_edges_.end() - 1) ? true : false;
      std::iter_swap(active_edges_.begin() + edge_idx, active_edges_.end() - 1); // O(1) 
      active_edges_.pop_back(); // O(1) 

      // update the internally stored value of idx_ for the edge
      // that was previously at the end of the container
      if (!no_swap) edges_[active_edges_[edge_idx]].idx_ = edge_idx;

      // before finishing, we need to update the incident edge container
      // for both node n1 and node n2 (O(n1.degree() + n2.degree()) in general, O(1) if sparse)
      std::vector<size_type>& n1_incident_edges = n1.fetch().incident_edges_;
      size_type pos_to_remove = 0; 
      for (auto ei = n1_incident_edges.begin(); ei != n1_incident_edges.end(); ++ei){
        if (*ei == edge_uid){
          break;
        }else{
          pos_to_remove += 1; 
        }
      }
      std::iter_swap(n1_incident_edges.begin() + pos_to_remove, n1_incident_edges.end() - 1);
      n1_incident_edges.pop_back(); 

      // repeat for n2 
      std::vector<size_type>& n2_incident_edges = n2.fetch().incident_edges_;
      pos_to_remove = 0; 
      for (auto ei = n2_incident_edges.begin(); ei != n2_incident_edges.end(); ++ei){
        if (*ei == edge_uid){
          break;
        }else{
          pos_to_remove += 1; 
        }
      }
      std::iter_swap(n2_incident_edges.begin() + pos_to_remove, n2_incident_edges.end() - 1);
      n2_incident_edges.pop_back(); 

      return 1;
    }
    else
    {

      // no edge to remove, return 0 
      return 0;
    }
  };

  /** Remove an edge from the graph, returning 1 if the edge was removed
   *  succesfully and 0 otherwise (e.g. if the edge did not exist). 
   *  This method works by simply calling the above remove_edge(node_type n1, node_type n2)
   *  method with @a n1 = e.node1() and @a n2 = e.node2().
   * 
   * @pre @a e is a valid edge of this graph
   * @return int 1 if edge removal is succesful, 0 otherwise. 
   * @post has_edge(e.node1(), e.node2()) == false
   * @post If old has_edge(e.node1(), e.node2()), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). 
   * 
   * Can invalidate edge iterators (type edge_iterator) -- this is because
   * remove_edge() alters the underlying std::vector (active_edges_) that
   * is backing the edge iterator. 
   *
   * Complexity: O(n1.degree() + n2.degree()) in general, which is O(1)
   *             for a sparse graph. 
   */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2()); 
  }

  /** Remove an edge from the graph, returning a valid edge iterator. 
   * @pre @a e_it is a valid edge iterator (type @a edge_iterator) of this graph
   * @return a valid edge iterator for this graph
   * @post has_edge(e_it->node1(), e_it->node2()) == false
   * @post If old has_edge(e_it->node1(), e_it->node2()), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). 
   * 
   * Can invalidate edge iterators (type edge_iterator) -- this is because
   * remove_edge() alters the underlying std::vector (active_edges_) that
   * is backing each edge iterator. Certainly invalidates @a e_it by construction. 
   *
   * Complexity: O(n1.degree() + n2.degree()) in general, which is O(1)
   *             for a sparse graph. 
   */
  edge_iterator remove_edge(edge_iterator e_it)
  {
    // assert that this iterator is not graph.edge_end()
    // otherwise, there is nothing to remove and the increment
    // operator below will cause a segfault
    assert(e_it != this->edge_end() && "invalid edge iterator provided");

    // with the way we have implemented remove_edge(), the element
    // at this current edge's index is one that has not yet
    // been visited. so we return an iterator starting here.
    edge_type e = *e_it; 
    size_type edge_idx = e.index();
    remove_edge(e); // O(n1.degree() + n2.degree()) in general, which is O(1) for a sparse graph.
    return EdgeIterator(this, active_edges_.begin() + edge_idx);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

    // Empty out containers
    nodes_.empty(); 
    edges_.empty(); 
    active_nodes_.empty();
    active_edges_.empty();

    // Call clear methods
    nodes_.clear(); 
    edges_.clear(); 
    active_nodes_.clear();
    active_edges_.clear(); 
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

    /** Defines equality between two node iterators as having the same 
     *    underlying STL iterator. */ 
    bool operator==(const node_iterator& node_iter) const
    {
      return node_iter.it_ == it_; 
    }

    /** Defines inequality between two node iterators as having different
     *    underlying STL iterators. */
    bool operator!=(const node_iterator& node_iter) const 
    {
      return node_iter.it_ != it_;
    }

    /* Defines the dereference operator for the node iterator. 
    *   Returns a proxy node having the same graph and index. */ 
    node_type operator*() const 
    {
      return Node(graph_, *it_); // return proxy node
    }

    /** Defines the increment operator for the node iterator.
     *    Increments the underlying STL iterator. */
    node_iterator& operator++()
    {
      ++it_; 
      return *this; // return a reference to this object
    }

   private: 
    // declare Graph as friend so it can access the private
    // member variables of NodeIterator
    friend class Graph;

    // Pointer back to the parent graph.
    const graph_type *graph_;

    /** STL iterator used for the node container */
    stl_node_iterator it_;

    /** Private constructor that can be accessed by the Graph class. */ 
    NodeIterator(const graph_type* graph, stl_node_iterator it) : graph_(graph), it_(it) {}
  };

  // Returns a node iterator pointing at the start of the node list.
  node_iterator node_begin() const
  {
    return node_iterator(this, active_nodes_.begin());
  }

  // Returns a node iterator pointing at the end (aka one past the last node) of the node list.
  node_iterator node_end() const
  {
    return node_iterator(this, active_nodes_.end()); 
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

    /** Defines equality between two incident iterators as having the same 
     *    underlying STL iterator. */
    bool operator==(const incident_iterator& incident_iter) const
    {
      return incident_iter.it_ == it_;
    }

    /** Defines inequality between two incident iterators as having different
     *    underlying STL iterators. */
    bool operator!=(const incident_iterator& incident_iter) const
    {
      return incident_iter.it_ != it_;
    }

    /* Defines the dereference operator for the incident iterator. 
    *   Returns a proxy edge having the same graph and edge index. 
    *   Care is taken to ensure that the edge returned is such that
    *   Edge::node1() will always return the node used to create 
    *   this iterator. */
    edge_type operator*() const
    {
      // To note:
      // - By default node1() returns the node with the smaller index of the edge
      // - This can be flipped by passing the flag "reversed" 
      // Thus, to ensure Edge::node1() returns the spawning node index in this case of 
      // IncidentIterator(), we need to turn on the "reversed" flag if it's the larger index
      size_type current_edge_uid = *it_;
      size_type small_idx = graph_->edges_[current_edge_uid].node_a; // guranteed to be the smaller node index
      size_type large_idx = graph_->edges_[current_edge_uid].node_b; // guranteed to be the larger node index 

      // if the spawning node's index is the smaller index, then nothing needs to be done.
      // otherwise we need to return an edge with the reverse flag turned on. 
      bool reversed = (spawning_node_idx_ == small_idx) ? false : true;

      // return proxy edge
      return Edge(graph_, *it_, reversed); 
    }

    /** Defines the increment operator for the incident iterator.
     *    Increments the underlying STL iterator. */
    incident_iterator operator++()
    {
      ++it_;
      return *this; 
    }

   private:
    friend class Graph;

    // STL iterator for the container for incident edges (which contains uid)
    stl_incident_iterator it_;

    // Pointer to the parent graph. Necessary for dereferencing since the incident edge container only stores indices. 
    const graph_type* graph_;

    // We need to keep track of the index of the spawning node so we can make sure to always return it via Edge::node1()
    size_type spawning_node_idx_;

    // Private constructor that can be accessed by the Graph class.
    IncidentIterator(stl_incident_iterator it, const graph_type *graph, size_type node_idx) 
      : it_(it), graph_(graph), spawning_node_idx_(node_idx) {}
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

    /** Defines equality between two edge iterators as having the same 
     *    underlying STL iterator. */
    bool operator==(const edge_iterator& edge_iter) const
    {
      return edge_iter.it_ == it_;
    }

    /** Defines inequality between two edge iterators as having different
     *    underlying STL iterators. */
    bool operator!=(const edge_iterator& edge_iter) const
    {
      return edge_iter.it_ != it_;
    }

    /* Defines the dereference operator for the edge iterator. 
    *   Returns a proxy edge having the same graph and edge index. */
    edge_type operator*() const
    {
      return Edge(graph_, *it_); // return proxy edge
    }

    /** Defines the increment operator for the edge iterator.
     *    Increments the underlying STL iterator. */
    edge_iterator& operator++()
    {
      ++it_;
      return *this; 
    }

   private:
    friend class Graph;

    // Pointer back to the parent graph
    const graph_type *graph_;

    // STL iterator for the edge container
    stl_edge_iterator it_;

    // Private constructor that can be accessed by the Graph class.
    EdgeIterator(const graph_type* graph, stl_edge_iterator it) : graph_(graph), it_(it) {}

  };

  // Returns an iterator pointing at the start of the edge list.
  edge_iterator edge_begin() const
  {
    return EdgeIterator(this, active_edges_.begin());
  }

  // Returns an iterator pointing at the end (aka one past the last edge) of the edge list.
  edge_iterator edge_end() const
  {
    return EdgeIterator(this, active_edges_.end());
  }

 private:

  // Internal type for nodes 
  // - Inspired by peercode from HW 0 
  struct internal_node_
  {
    // Pointer back to the parent graph
    Graph *my_graph_;

    // node's ID in the global list of nodes for the graph
    // this is the user-facing index 
    size_type idx_;

    // internal index of this node
    size_type uid_;

    // position of the node
    Point pos_;

    // value of this node
    node_value_type value_;

    // keep track of incident edges for this node
    std::vector<size_type> incident_edges_; 

    // constructor for internal_node_
    internal_node_(const Graph *graph, const size_type idx, const size_type uid, const Point pos, const node_value_type value = node_value_type())
        : my_graph_(const_cast<Graph *>(graph)), idx_(idx), uid_(uid), pos_(pos), value_(value), incident_edges_() {}
  };

  // Internal type for edges
  struct internal_edge_
  {
    // Pointer back to the parent graph
    Graph *my_graph_;

    // edge's ID in the global list of edge for the graph
    // this is the user-facing index
    size_type idx_;

    // internal index of this edge
    size_type uid_; 

    size_type node_a; // id (idx) of the first node of the edge. Node with the smaller index. 
    size_type node_b; // id (idx) of the second node of the edge. Node with the larger index.

    // value of this edge
    edge_value_type value_;

    // constructor for internal_edge_. Inspired by peercode for HW0. 
    internal_edge_(const Graph *graph, const size_type idx, const size_type uid, const size_type node_a, const size_type node_b, const edge_value_type value = edge_value_type())
        : my_graph_(const_cast<Graph *>(graph)), idx_(idx), uid_(uid), value_(value) 
        {
          // this ensures internally node_a always refers to the node with the smaller index and node_b always
          // refers to the larger index 
          this->node_a = std::min(node_a, node_b);
          this->node_b = std::max(node_a, node_b);
        }
  };

  // Containers for nodes (their positions) and edges 
  // these containers store info for any nodes/edges ever present 
  std::vector<internal_node_type> nodes_; // indexed by node uid
  std::vector<internal_edge_type> edges_; // idnex by edge uid 

  // Containers to keep track of which nodes are "active"
  // i.e. they store the internal indices (uid) of the active nodes/edges.
  // Thus the size of these containers <= their counterparts above. 
  std::vector<size_type> active_nodes_; // indexed by node idx (user-facing index)
  std::vector<size_type> active_edges_; // indexed by edge idx (user-facing index)
};

#endif // CME212_GRAPH_HPP
