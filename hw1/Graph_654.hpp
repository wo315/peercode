#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <list>
#include <unordered_map>
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
class Graph : private totally_ordered<Graph<V> > {
  public:
    using node_value_type = V;
    using size_type = unsigned;
  private:
    /** Vector of points (which define a node) with the node index being
     * the location in the vector */
    std::vector<Point> graph_points_;
    /*Number of nodes in the graph */
    size_type num_Nodes_;
    
    /** Map of edge indexes as keys, and values being a list of the indexes
     * of the nodes that this edge connects. */
    std::map<size_type, std::pair<size_type, size_type> > graph_edges_;
    
    /** A map of node indexes as keys, with values being a map with node
     * indexes as its keys, with the final values being the edge indexes.*/
    std::map<size_type, std::map<size_type, size_type> > edge_map_;
    /** Number of edges in the graph */
    size_type num_Edges_;
    /** A map with node indexes as keys and node values as values. */
    std::map<size_type, node_value_type> node_values_;
    
  public:
  //
  // PUBLIC TYPE DEFINITIONS
  //
  /** Type of this graph. */
  using graph_type = Graph<node_value_type>;

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
  //using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : graph_points_(), num_Nodes_(0), num_Edges_(0){
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
        return this->fetch();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        return this->my_index_;
    }
    
    /** Return this node's graph pointer*/
    Graph* graph() const {
        return this->my_graph_;
    }

    /** Return the node value type of the graph class in a mutable type. */
    node_value_type& value(){
        return this->my_graph_->node_values_[this->my_index_];
    }
      
    /** Return the node value type of the graph class in a const type. */
    const node_value_type& value() const{
        return this->my_graph_->node_values_[this->my_index_];
    }
      
    /**Return the number of incident edges */
    size_type degree() const {
        return this->my_graph_->edge_map_[my_index_].size();
    }
      
    /**Start of the incident iterator */
    incident_iterator edge_begin() const {
        return IncidentIterator(my_graph_,
                                my_index_,
                                my_graph_->edge_map_[my_index_].end(),
                                my_graph_->edge_map_[my_index_].begin());
    }
      
    /**End of the incident iterato r*/
    incident_iterator edge_end() const{
        return IncidentIterator(my_graph_,
                                my_index_,
                                my_graph_->edge_map_[my_index_].end(),
                                my_graph_->edge_map_[my_index_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.graph() == this->my_graph_ && n.index() == this->my_index_);
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
        if (this->index() < n.index()){
                  return true;
        } else if (this->index() == n.index()){
            return norm(this->position()) < norm(n.position());
        }
        return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* my_graph_;
    size_type my_index_;
      
    /** Private constructor, allows us to pass over the graph and the index value */
    Node(const graph_type* g, size_type i)
    : my_graph_(const_cast<graph_type*>(g)), my_index_(i){
    }
      
    /** Used for proxy implementation. Returns a proxy Point. */
    Point& fetch() const {
        return my_graph_->graph_points_[my_index_];
    }
      
  }; /// End Node class
    
  /** Return the number of nodes in the graph. */
  size_type size() const {
      return this->num_Nodes_;
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
                const node_value_type& vt = node_value_type()) {
      size_type old_num_Nodes = num_Nodes_;
      this->num_Nodes_ += 1;
      this->graph_points_.insert(graph_points_.begin() + old_num_Nodes,
                                 position);
      (void)vt; // Quiet compiler warning. Not used in current implementation.
      return Node(this, old_num_Nodes);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      if (n.my_graph_ == this){
          return true;
      } else {
          return false;
      }
  }
     
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
      // Returning a proxy
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return the first node of this Edge, via a node proxy. */
    Node node1() const {
        return Graph::Node(my_edge_graph_,
                        my_edge_graph_->graph_edges_[my_edge_index_].first);
    }

    /** Return the second node of this Edge, via a node proxy. */
    Node node2() const {
        return Graph::Node(my_edge_graph_,
                        my_edge_graph_->graph_edges_[my_edge_index_].second);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        return (norm(node1().position()) ==
                norm(e.node1().position())) &&
            ((norm(node2().position()) ==
              norm(e.node2().position())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        // Sort by norm of node1, then tiebreak via norm of node2
        if (norm(node1().position()) == norm(e.node1().position())) {
            return norm(node2().position()) < norm(e.node2().position());
        } else {
        return norm(node1().position()) < norm(e.node1().position());
        }
    }

   private:
    friend class Graph;
    Graph* my_edge_graph_;
    size_type my_edge_index_;
    size_type node1_index_;
    size_type node2_index_;
      
    /** Private constructor
     * An edge is defined by a unique index, a graph, and two nods.
     * The nodes are identified by their indexes.
     */
    Edge(const Graph* g, size_type i, const Node& a, const Node& b)
      : my_edge_graph_(const_cast<Graph*>(g)), my_edge_index_(i),
      node1_index_(a.index()), node2_index_(b.index()){}
      
  }; /// End Edge class

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
      return this->num_Edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
      // Return a proxy Edge
      return Edge(this, i, node(graph_edges_.at(i).first),
                  node(graph_edges_.at(i).second));
  }

  /** Return the edge with index @a i, but with nodes in reverse order.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge_swap(size_type i) const {
    // Return a proxy Edge
    return Edge(this, i, node(graph_edges_.at(i).second),
                    node(graph_edges_.at(i).first));
  }
    
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
      size_type a_index = a.index();
      size_type b_index = b.index();
      
      if (edge_map_.find(a_index) != edge_map_.end()){
          if(edge_map_.at(a_index).find(b_index) !=
             edge_map_.at(a_index).end()){
              return true;
          }
          //assert(true);
      } else if (edge_map_.find(b_index) != edge_map_.end()){
          if(edge_map_.at(b_index).find(a_index) !=
                    edge_map_.at(b_index).end()){
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
      // Check to see if the edge exists
      size_type a_index = a.index();
      size_type b_index = b.index();
      if (edge_map_.find(a_index) != edge_map_.end()){
          if(edge_map_.at(a_index).find(b_index) !=
             edge_map_.at(a_index).end()){
              return Edge(this, edge_map_[a_index].at(b_index),a,b);
          }
      } else if (edge_map_.find(b_index) != edge_map_.end()){
          if(edge_map_.at(b_index).find(a_index) !=
             edge_map_.at(b_index).end()){
              return Edge(this, edge_map_[b_index].at(a_index),a,b);
          }
      }
      
      // Create a new edge
      size_type new_edge_index = num_Edges_;
      this->graph_edges_[new_edge_index].first = a.index();
      this->graph_edges_[new_edge_index].second = b.index();
      this->num_Edges_++;
      // Add to our edge_map_ an edge between these nodes
      // Adds edges in both directions to account for undirectness
      this->edge_map_[a.index()][b.index()] = new_edge_index;
      this->edge_map_[b.index()][a.index()] = new_edge_index;
      return Edge(this, new_edge_index,a,b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      num_Nodes_ = 0;
      num_Edges_ = 0;
      graph_points_.clear();
      graph_edges_.clear();
      edge_map_.clear();
  }
    
  //
  // Node Iterator
  //
    
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
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
      
    /** Constructor for a valid NodeIterator. */
    NodeIterator(const graph_type* g, size_type i){
        my_graph_ = const_cast<Graph<node_value_type> *>(g),
        current_iterator_node_index = i;
    }
      
    /** Dereferencing operator. */
    Node operator*() const{
        return my_graph_->node(current_iterator_node_index);
    }
      
    /** Increment operator. */
    NodeIterator& operator++(){
        if (current_iterator_node_index == (my_graph_->num_nodes() - 1)){
            this->current_iterator_node_index = -1;
            return *this;
        }
        ++this->current_iterator_node_index;
        return *this;
    }
    
    /** Equaliity operator. */
    bool operator==(const NodeIterator& x) const {
        if (my_graph_ != x.graph()){
            return false;
        }
        if (this->current_iterator_node_index == x.index()){
            return true;
        }
        return false;
    }
      
    /** Return this NodeIterator's graph pointer.  */
    graph_type* graph() const {
        return this->my_graph_;
    }
      
    /** Return this NodeIterator's current node index. */
    size_type index() const {
        return this->current_iterator_node_index;
    }

   private:
    friend class Graph;
    size_type current_iterator_node_index;
    graph_type* my_graph_;
      
  }; /// End NodeIterator class

  /** Return the beginning NodeIterator*/
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
    
  /** Return the end NodeIterator*/
  node_iterator node_end() const{
    return NodeIterator(this, -1);
  }
    
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
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
    /** Constructor for a valid IncidentIterator. */
    IncidentIterator(const graph_type* g,
                       size_type my_index,
                       std::map<size_type,size_type>::iterator e,
                       std::map<size_type,size_type>::iterator i){
        my_graph_ = const_cast<graph_type*>(g);
        my_node_index_ = my_index;
        my_end_iterator_ = e;
        my_map_iterator_ = i;
    }
      
    /** Dereferencing operator. */
    Edge operator*() {
        Edge edge1 = my_graph_->edge(this->my_map_iterator_->second);
        if (edge1.node1() == node(my_node_index_)){
            return edge1;
        }
        return my_graph_->edge_swap(this->my_map_iterator_->second);
    }
    
    /** Incrementing operator. */
    IncidentIterator& operator++(){
        if (my_map_iterator_ != my_end_iterator_){
            ++my_map_iterator_;
            return *this;
        }
    }
    
    /** Equality operator. */
    bool operator==(const IncidentIterator& x) const {
        return my_map_iterator_ == x.my_map_iterator_();
    }
      
    /** Return the map iterator tied to this IncidentIterator. This map iterator poitns to
     * one of the sub-maps within the edge_map_ map-of-maps. edge_map_ has contains
     * node indexes which point to node indexes which point to edge indexes in a nested map
     * type structure. */
    std::map<size_type,size_type>::iterator get_iterator(){
        return this->my_map_iterator_;
    }
    
   private:
    friend class Graph;
    std::map<size_type,size_type>::iterator my_map_iterator_;
    graph_type* my_graph_;
    size_type my_node_index_;
    std::map<size_type,size_type>::iterator my_end_iterator_;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    using iterator_type = std::map<size_type,
      std::pair<size_type,size_type> >::iterator;
      
    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }
    
    /** Constructor for a valid EdgeIterator. */
    EdgeIterator(const graph_type* g, int i ) {
        my_graph_ = const_cast<graph_type*>(g);
        if (i == -1) {
            my_map_iterator_ = my_graph_->graph_edges_.end();
        } else {
            my_map_iterator_ = my_graph_->graph_edges_.begin();
        }
    }
      
    /** Dereferencing operator. */
    Edge operator*() const {
        return my_graph_->edge(my_map_iterator_->first);
    }
      
    /** Increment operator. */
    EdgeIterator& operator++() {
        if (my_map_iterator_ == my_graph_->graph_edges_.end()){
            return *this;
        }
        ++my_map_iterator_;
        return *this;
    }
     
    /** Equality operator. */
    bool operator==(const EdgeIterator& x) const {
        return my_map_iterator_ == x.get_iterator();
    }
    
    /** Returns an iterator to a map which links edge indexes to node inxeses to node indexes, in
     * a nested map structure. */
    iterator_type get_iterator() const{
          return this->my_map_iterator_;
    }
      
   private:
    friend class Graph;
    iterator_type my_map_iterator_;
    graph_type* my_graph_;
      
  }; /// End of EdgeIterator class

  /** Returns the beginning EdgeIterator. */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }
  
  /** Returns the end EdgeIterator. */
  edge_iterator edge_end() const{
    return EdgeIterator(this, -1);
  }

}; /// End Graph class

#endif // CME212_GRAPH_HPP
