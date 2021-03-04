#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/**@file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <list>
#include <unordered_map>
//--functionality_0
//--missing header <map>
//--START
#include <map>
//--END
#include <vector>
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph : private totally_ordered<Graph<V, E> > {
  public:
    // Declare types
    typedef V node_value_type;
    typedef E edge_value_type;

    /** Type of indexes and sizes. */
    using size_type = unsigned;

  private:
    /** Vector which maps unique node id's to points */
    std::vector<Point> graph_points_;

    /** Vector which maps user-facing node indexes to unique id's*/
    std::vector<size_type> i2u_;

    /** Stores the total number of nodes ever added to this graph,
     *  inclusive of those that might later be deleted. */
    size_type total_num_Nodes_;

    /** Map from unique id's to user-facing indexes*/
    std::map<size_type, size_type> uid2idx_;

    /** Number of active nodes in the graph */
    size_type num_Nodes_;

    /** Map with edge indexes as keys, with corresponding values being
     *  a pair of the unique ids of the nodes that this edge connects. */
    std::map<size_type, std::pair<size_type, size_type> > graph_edges_;

    /** A map of unique node ids as keys, with values being a map with unique node
     * ids as its keys, with the final values being the edge indexes that connect the
     * two nodes together. */
    std::map<size_type, std::map<size_type, size_type> > edge_map_;

    /** Number of active edges in the graph */
    size_type num_Edges_;

    /** A map with unique node indexes as keys and node values as values. */
    std::unordered_map<size_type, node_value_type> node_values_;

    /** A map with edge indexes as keys and edgevalues as values. */
    std::unordered_map<size_type, edge_value_type> edge_values_;

  public:
  //
  // PUBLIC TYPE DEFINITIONS
  //
  /** Type of this graph. */
  using graph_type = Graph<node_value_type, edge_value_type>;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** @function Graph
   *  @brief Graph class constructor
   *  Initializes variables and generates an empty graph.
   */
  Graph()
    : graph_points_(), total_num_Nodes_(0), num_Nodes_(0), num_Edges_(0){
  }

  /** @function ~Graph
   *  @brief Graph class destructor
   *  Destructs graph object
   */
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

    /**@function Node::Node
     * @brief Node constructor
     * Construct an invalid node.
     */
    Node() {
    }

    /**@function Node::position
     * @brief Return a reference to this node's position
     */
    Point& position(){
        return this-> fetch();
    }

    /**@function Node::position
     * @brief Return a constant reference to this node's position
     */
    const Point& position() const {
        return this->fetch();
    }

    /**@function Node::index
     * @brief Return a this node's user-facing index
     * Always in the range [0, graph_size).
     */
    size_type index() const {
        return this->my_index_;
    }

    /**@function Node::graph
     * @brief Return a this node's graph pointer
     */
    Graph* graph() const {
        return this->my_graph_;
    }

    /**@function Node::value
     * @brief Returns a reference to this node's node_value_type
     * @return Returns a reference to a node_value_type object
     * Complexity: Constant average time complexity, since underlying
     * is an unordered_map. Worst case O(num_nodes).
     */
    node_value_type& value(){
        return this->my_graph_->node_values_[my_graph_->i2u_[this->my_index_]];
    }

    /**@function Node::value
     * @brief Returns a constant reference to this node's node_value_type
     * @return Returns a constant reference to a node_value_type object
     * Complexity: Constant average time complexity, since underlying
     * is an unordered_map. Worst case O(num_nodes).
     */
    const node_value_type& value() const{
        return this->my_graph_->node_values_[my_graph_->i2u_[this->my_index_]];
    }

    /**@function Node::degree
     * @brief Return the number of incident edges to this node
     * @return Return the number of incident edges to this node.
     */
    size_type degree() const {
        return this->my_graph_->edge_map_[
                                    my_graph_->i2u_[this->my_index_]].size();
    }

    /**@function Node::edge_begin()
     * @brief Return beginning incident iterator to calling object node.
     * @return Beginning incident iterator for calling object node
     * An incident iterator allows us to iterate over all of the incident
     * edges to a node.
     */
    incident_iterator edge_begin() const {
        return IncidentIterator(my_graph_,
                                my_index_,
                my_graph_->edge_map_[my_graph_->i2u_[this->my_index_]].end(),
                my_graph_->edge_map_[my_graph_->i2u_[this->my_index_]].begin());
    }

    /**@function Node::edge_end()
     * @brief Return end incident iterator to calling object node.
     * @return End incident iterator for calling object node
     */
    incident_iterator edge_end() const{
        return IncidentIterator(my_graph_,
                                my_index_,
                my_graph_->edge_map_[my_graph_->i2u_[this->my_index_]].end(),
                my_graph_->edge_map_[my_graph_->i2u_[this->my_index_]].end());
    }

    /**@function Node::operator==()
     * @brief Equality operator for node objects
     * @param[in]       n   Node
     * @return boolean indicating equality
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n.graph() == this->my_graph_ && n.index() == this->my_index_);
    }

    /**@function Node::operator<()
     * @brief Less than operator for node objects
     * @param[in]       n   Node
     * @return boolean indicating less than status
     * Equal nodes have the same graph and the same index.
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

    /**@function Node::Node
     * @brief Private node constructor
     * @param[in]    g   Pointer to a graph
     * @param[in]    i   User-facing index of node
     * Constructs a valid node
     */
    Node(const graph_type* g, size_type i)
    : my_graph_(const_cast<graph_type*>(g)), my_index_(i){
    }

    /** Returns the point associated with this node by using
     * the node's user-facing idx identifier. */

    /**@function Node::fetch()
     * @brief Grabs the point associated with the calling node.
     * @return Point reference associated with the calling node.
     * Finds point associated with this node in graph datastructure.
     */
    Point& fetch() const {
        return my_graph_->graph_points_[my_graph_->i2u_[my_index_]];
    }
  }; /// End Node class

  /**@function Graph::size()
   * @brief Returns number of active nodes in the graph
   * @return Number of active nodes in the graph
   */
  size_type size() const {
      return this->num_Nodes_;
  }

  /**@function Graph::num_nodes()
   * @brief Returns number of active nodes in the graph
   * @return Number of active nodes in the graph
   */
  size_type num_nodes() const {
      return size();
  }

  /**@function Graph::add_node
   * @brief  Add a node to the graph, returning the added node.
   * @param[in] position   The new node's position
   * @param[in] vt                Value used to initialize this node's node_value_type
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @return New Node object
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& vt = node_value_type()) {
      size_type node_index = num_Nodes_;

      // Add this nodes idx to our vector of unique id's
      this->i2u_.push_back(node_index);
      // Add this node uid to our map of uids to idx's
      this->uid2idx_[total_num_Nodes_] = node_index;

      assert(uid2idx_[i2u_[node_index]] == node_index);

      // Add this node's position to our mapping from unique id's
      // to points.
      this->graph_points_.push_back(position);
      node_values_[total_num_Nodes_] = vt;

      this->num_Nodes_ += 1;
      this->total_num_Nodes_ += 1;
      return Node(this, node_index);
  }

  /**@function Graph::has_node()
   * @brief Determine if a Node belongs to this Graph
   * @param[in]     n   Node who's existence we wish to determine.
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      return n.my_graph_ == this;
  }

  /**@function Graph::node()
   * @brief Return the node with index @a i.
   * @param[in]     i   index of node
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   * @return Node located at user-facing index i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
      // Returning a proxy
      return Node(this, i);
  }

  /**@function Graoh::remove_node()
   * @brief removes a node from the graph
   * @param[in]     n   Node
   * @pre n is a valid (active) node, with a valid index.
   * @post n and all of its incident edges have been
   *  removed from the graph
   * @post node at this index is removed
   * @post all incident edges are removed
   * @post all active node iterators are invalidated.
   * @return 1 if a node was removed, 0 otherwise
   * Removes a node using a reference to the node,
   * including all of its incident edges
   * Complexity: Average O(d), where d is the number if incident edges
   * Iterating through incident edges takes O(d) time, erasing each
   * edge takes on average constant time, and the swap-and-pop
   * operation takes constanttime.
   */
   size_type remove_node(const Node& n){
       // Remove all incident edges
       std::vector<Edge> incident_edges;
       auto first = n.edge_begin();
       auto last = n.edge_end();
       while (first != last){
           incident_edges.push_back(*first);
           ++first;
       }
       // Now we have a vector of all incident edges.
       // Call delete edge on each
       while (incident_edges.size() != 0){
           remove_edge(incident_edges.back());
           incident_edges.pop_back();
       }
       // Remove this node's idx as a key from edge_map_
       this->edge_map_.erase(n.index());

       // Remove node from i2u_
       size_type this_node_idx = n.index();
       size_type last_node_uid = i2u_.back();

       // Swap and Pop
       i2u_[this_node_idx] = last_node_uid;
       i2u_.pop_back();

       // Reindex swapped node
       uid2idx_[last_node_uid] = this_node_idx;

       // Check:
       assert(i2u_[this_node_idx] == last_node_uid);
       num_Nodes_ -= 1;
       return 1;
   }

  //
  // EDGES
  //

  /**@class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /**@function Graph::Edge
     * @brief Edge constructor
     * Construct an invalid Edge.
     */
    Edge() {
    }

    /**@function Edge::value
     * @brief Returns a reference to this edge's edge_value_type
     * @return Reference to a edge_value_type object
     * Complexity: Constant average time complexity, since underlying
     * is an unordered_map. Worst case O(num_edges)
     */
    edge_value_type& value() {
        return this->my_edge_graph_->edge_values_[this->my_edge_index_];
    }

    /**@function Edge::value
     * @brief Returns a constant reference to this edge's edge_value_type
     * @return Constant reference to a edge_value_type object
     * Complexity: Constant average time complexity, since underlying
     * is an unordered_map. Worst case O(num_edges)
     */
    const edge_value_type& value() const {
        return this->my_edge_graph_->edge_values_[this->my_edge_index_];
    }


    /**@function Edge::node1()
     * @brief Returns this edge's first node
     * @return Proxy to this edge's first node.
     * Complexity: O(log(num_edges))
     */
    Node node1() const {
        return Graph::Node(my_edge_graph_, my_edge_graph_->uid2idx_[
                        my_edge_graph_->graph_edges_[my_edge_index_].first]);
    }

    /**@function Edge::node2()
     * @brief Returns this edge's second node
     * @return Proxy to this edge's second node.
     * Complexity: O(log(num_edges))
     */
    Node node2() const {
        return Graph::Node(my_edge_graph_, my_edge_graph_->uid2idx_[
                        my_edge_graph_->graph_edges_[my_edge_index_].second]);
    }

    /**@function Edge::length()
     * @brief Returns this edge's euclidean length
     * @return Edge's euclidean length
     * Complexity: O(log(num_edges))
     */
    double length() const{
        return norm(node1().position() - node2().position());
    }

    /**@function Edge::operator==()
     * @brief Equality operator for edge objects
     * @param[in]      e   Edge to check equality with
     * @return boolean indicating equality
     * Equal edges connect the same nodes.
     * Edge equality is valid even if node1 and node2 are swapped
     */
    bool operator==(const Edge& e) const {
        return ((my_edge_graph_->node(node1_index_) == e.node1() &&
                my_edge_graph_->node(node2_index_) == e.node2()) ||
                (my_edge_graph_->node(node1_index_) == e.node2() &&
                my_edge_graph_->node(node2_index_) == e.node1()));
    }

    /**@function Edge::operator<()
     * @brief Less than operator for edge objects
     * @param[in]       e   Edge to check against
     * @return boolean indicating less than status
     * Smaller edges have shorter euclidean lengths
     */
    //--functionality_1
    //--edge comparion should not be based on position
    //--START
    bool operator<(const Edge& e) const {
        double my_size = norm(my_edge_graph_->node(node1_index_).position() -
                              my_edge_graph_->node(node1_index_).position());
        return my_size < e.length();
    }
    //--END

   private:
    friend class Graph;
    Graph* my_edge_graph_;
    size_type my_edge_index_;
    size_type node1_index_;
    size_type node2_index_;

    /**@function Edge::Edge
     * @brief Private Edge constructor
     * @param[in]    g   Pointer to a graph
     * @param[in]    i   User-facing index of edge
     * @param[in]    a   First node for edge
     * @param[in]    b   Second node for edge
     * @pre graph-pointer g is a valid graph-pointer
     * @pre Nodes a and b are valid nodes.
     * Constructs a valid edge
     */
    Edge(const Graph* g, size_type i, const Node& a, const Node& b)
      : my_edge_graph_(const_cast<Graph*>(g)), my_edge_index_(i),
      node1_index_(a.index()), node2_index_(b.index()){}
  }; /// End Edge class

  /**@fucntion Graph::num_edges()
   * @brief Return the total number of edges in the graph.
   * @return Total number of edges.
   * Complexity: O(1).
   */
  size_type num_edges() const {
      return this->num_Edges_;
  }

  /**@function Graph::edge()
   * @brief Return the edge with index @a i.
   * @param[in]     i   index of edge
   * @pre 0 <= @a i < num_edges()
   * @return Proxy edge equal to the edge with index i.
   * Complexity: O(log(num_nodes) + log(num_edges))
   */
  Edge edge(size_type i) const {
      // Return a proxy Edge
      return Edge(this, i, node(uid2idx_.at(graph_edges_.at(i).first)),
                  node(uid2idx_.at(graph_edges_.at(i).second)));
  }

  /**@function Graph::edge_swap()
   * @brief Return edge with index @a i, but nodes swapped.
   * @param[in]     i   index of edge
   * @pre 0 <= @a i < num_edges()
   * @return Proxy edge equal to the edge with index i, but with nodes swapped.
   * Complexity: O(log(num_nodes) + log(num_edges))
   */
  Edge edge_swap(size_type i) const {
    // Return a proxy Edge
    return Edge(this, i, node(uid2idx_.at(graph_edges_.at(i).second)),
                    node(uid2idx_.at(graph_edges_.at(i).first)));
  }

  /**@function Graph::has_edge()
   * @brief Test whether two nodes are connected by an edge.
   * @param[in]     a   First node of edge we wish to check for
   * @param[in]     b   Second node of edge we wish to check for
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: Average O(log(num_edges) + log(num_nodes))
   */
  bool has_edge(const Node& a, const Node& b) const {
      size_type a_index = i2u_.at(a.index());
      size_type b_index = i2u_.at(b.index());

      if (edge_map_.find(a_index) != edge_map_.end()){
          if(edge_map_.at(a_index).find(b_index) !=
             edge_map_.at(a_index).end()){
              return true;
          }
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
   * @param[in]     a   First node of edge we wish to add
   * @param[in]     b   Second node of edge we wish to add
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: Average O(log(num_edges) + log(num_nodes))
   */
  Edge add_edge(const Node& a, const Node& b) {
      // Check to see if the edge exists
      size_type a_index = i2u_.at(a.index());
      size_type b_index = i2u_.at(b.index());
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
      this->graph_edges_[new_edge_index].first = i2u_.at(a.index());
      this->graph_edges_[new_edge_index].second = i2u_.at(b.index());
      this->num_Edges_+= 1;
      // Add to our edge_map_ an edge between these nodes
      // Adds edges in both directions to account for undirectness
      this->edge_map_[a_index][b_index] = new_edge_index;
      this->edge_map_[b_index][a_index] = new_edge_index;
      return Edge(this, new_edge_index,a,b);
  }

  /**@function Graph::clear()
   * @brief Remove all nodes and edges from this graph.
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

    /**@function Graph::remove_edge()
     * @brief Remove edge with provides nodes
     * @param[in]       n1      First Node of edge
     * @param[in]       n2      Second Node of edge
     * @pre Both n1 and n2 are valid nodes of this graph.
     * @post Edge connecting these nodes is removed
     * @post Any active edge iterator is invalidated.
     * @post Does not guarantee that adding future edges
     *  will not cause problems since this does not keep track
     *  of a unique id. It just removes the edge's index
     *  from the various maps and decrements the num_edges.
     * @return 1 if edge was removed, 0 otherwise
     *
     * Complexity: Average O(log(num_edges) + log(num_nodes))
     */
    size_type remove_edge(const Node& n1, const Node& n2){
        // See if this edge does not exist
        if (!has_edge(n1, n2)){
            // Did not perform a removal
            return 0;
        }

        size_type n1_index = i2u_.at(n1.index());
        size_type n2_index = i2u_.at(n2.index());

        size_type e_index =
                this->edge_map_[n1_index][n2_index];
        // Remove from graph_edges_, a map of edge idx indexes to
        // pairs of node indexes.
        this->graph_edges_.erase(e_index);
        // Remove from datastructure containing data on edge
        //this->edge_values_.erase(e.index);

        // Remove from edge_map_, a map of node indexes with each node index's
        // value being a map of adjacent node indexes, with this child map's
        // value being the edge index connecting the two nodes.
        this->edge_map_[n1_index].erase(n2_index);
        this->edge_map_[n2_index].erase(n1_index);
        this->num_Edges_ -= 1;
        return 1;
    }

    /**@function Graph::remove_edge()
     * @brief Remove edge with provided nodes
     * @param[in]       edge      edge to remove
     * @pre Provided edge is a valid edge of this graph.
     * @post Provied edge is removed
     * @post Any active edge iteratoris invalidated
     * Does not guarantee that adding future edges
     * will not cause problems since this does not keep track
     * of a unique id. It just removes the edge's index
     * from the various maps and decrements the num_edges.
     * @return 1 if preconditions are met
     *
     * Complexity: Average O(log(num_edges) + log(num_nodes))
     */
    size_type remove_edge(const Edge& edge){
        size_type n1_index = i2u_.at(edge.node1().index());
        size_type n2_index = i2u_.at(edge.node2().index());

        size_type e_index =
                this->edge_map_[n1_index][n2_index];
        this->graph_edges_.erase(e_index);
        this->edge_map_[n1_index].erase(n2_index);
        this->edge_map_[n2_index].erase(n1_index);
        this->num_Edges_ -= 1;
        return 1;
    }

  //
  // Node Iterator
  //

  /**@class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /**@function NodeIterator::NodeIterator()
     * @brief Constructs an  invalid NodeIterator
     */
    NodeIterator() {
    }

    /**@function NodeIterator::NodeIterator()
     * @brief Constructs a valid NodeIterator
     * @param[in]   g   pointer to  graph
     * @param[in]   i   Index of underlying node that this iterator points to
     */
    NodeIterator(const graph_type* g, size_type i){
        my_graph_ = const_cast<graph_type *>(g),
        current_iterator_node_index = i;
    }

    /**@function NodeIterator::operator*()
     * @brief Dereferencing operator for calling NodeIterator
     * @return Underlying Node
     */
    Node operator*() const{
        return my_graph_->node(current_iterator_node_index);
    }

    /**@function NodeIterator::operator++()
     * @brief Incrementing operator for calling NodeIterator
     * @return Reference to incremented NodeIterator
     * Returns end iterator if iterator is at last node.
     */
    NodeIterator& operator++(){
        if (current_iterator_node_index == (my_graph_->num_nodes() - 1)){
            this->current_iterator_node_index = -1;
            return *this;
        }
        ++this->current_iterator_node_index;
        return *this;
    }

    /**@function NodeIterator::operator==()
     * @brief Equality operator for calling NodeIterator
     * @param[in]   x       NodeIterator to compare against
     * @return boolean indicating equality
     * Equal node iterators have equal graphs and equal underlying
     * nodes that they point to, verified using user-facing node index.
     */
    bool operator==(const NodeIterator& x) const {
        if (my_graph_ != x.graph()){
            return false;
        }
        if (this->current_iterator_node_index == x.index()){
            return true;
        }
        return false;
    }

    /**@function NodeIterator::graph()
     * @brief Return this NodeIterator's graph pointer.
     * @return A pointer to this NodeIterator's graph.
     */
    graph_type* graph() const {
        return this->my_graph_;
    }

    /**@function NodeIterator::index()
     * @brief Return the index of this iterator's underlying node
     * @return Index of underlying node.
     */
    size_type index() const {
        return this->current_iterator_node_index;
    }

   private:
    friend class Graph;
    size_type current_iterator_node_index;
    graph_type* my_graph_;

  }; /// End NodeIterator class

  /**@function Graph::node_begin()
   * @brief Return the beginning node iterator for this graph
   * @return Beginning node iterator
   */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /**@function Graph::node_end()
   * @brief Return the end Node Iterator for this graph
   * @return End node iterator
   */
  node_iterator node_end() const{
    return NodeIterator(this, -1);
  }
  /**@function Graph::remove_node()
   * @brief Removes node pointed to by provided node iterator.
   * @param     n_it    Iterator which points to the the node to remove
   * @pre Provided node iterator points to a valid (active) node in the graph.
   * @post The node pointed to by this node iterator is removed
   * @post All of this node's incident edges are removed
   * @post All active node iterators are invalidated
   * @return Beginning node iterator for newly editted graph
   * Complexity: Average O(d), where d is the number if incident edges
   * Iterating through incident edges takes O(d) time, erasing each
   * edge takes on average constant time, and the swap-and-pop
   * operation takes constanttime.
   */
  //--functionality_1
  //--node_begin may already been removed
  //--START
  node_iterator remove_node(node_iterator n_it){
    //--functionality_0
    //--typo: remove_edge -> remove_node
    //--START
    remove_node(*n_it);
    //--END
    return node_begin();
  }
  //--END

  //
  // Incident Iterator
  //

  /**@class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /**@function IncidentIterator::IncidentIterator()
     * @brief Constructs an  invalid IncidentIterator
     */
    IncidentIterator() {
    }

    /**@function IncidentIterator::IncidentIterator()
     * @brief Constructs a valid Incidentiterator
     * @param[in]   g   pointer to  graph
     * @param[in]   my_index   Index of underlying node that this iterator points to
     * @param[in]   e   End iterator for child map inside the edge_map_ datastructure.
     * @param[in]   i   Current/beginning iterator for child mpa inside the edge_map_ datastructure.
     */
    IncidentIterator(const graph_type* g,
                       size_type my_index,
                       std::map<size_type,size_type>::iterator e,
                       std::map<size_type,size_type>::iterator i){
        my_graph_ = const_cast<graph_type*>(g);
        my_node_index_ = my_index;
        my_end_iterator_ = e;
        my_map_iterator_ = i;
    }

    /**@function IncidentIterator::operator*()
     * @brief Dereferencing operator for calling IncidentIterator
     * @return Underlying edge that this incident iterator currently points to
     */
    Edge operator*() {
        Edge edge1 = my_graph_->edge(this->my_map_iterator_->second);
        if (edge1.node1() == my_graph_->node(my_node_index_)){
            return edge1;
        }
        // Swap node1 with node 2
        size_type node1index =
        my_graph_->graph_edges_[my_map_iterator_->second].first;
        my_graph_->graph_edges_[my_map_iterator_->second].first =
        my_graph_->graph_edges_[my_map_iterator_->second].second;
        my_graph_->graph_edges_[my_map_iterator_->second].second = node1index;
        // Swap using edge index
        return my_graph_->edge_swap(this->my_map_iterator_->second);
    }

    /**@function IncidentIterator::operator++()
     * @brief Incrementing operator for calling IncidentIterator
     * @return Incremented incident iterator
     * Returns end incident iterator if current iterator points to last element.
     */
    IncidentIterator& operator++(){
        if (my_map_iterator_ != my_end_iterator_){
            ++my_map_iterator_;
            return *this;
        }
        return *this;
    }

    /**@function IncidentIterator::operator==()
     * @brief Equality operator for calling IncidentIterator
     * @param[in]     x       IncidentIterator to compare against
     * @return boolean indicating equality
     * Equal incident iterators point to the same edge.
     */
    bool operator==(const IncidentIterator& x) const {
        return my_map_iterator_ == x.get_iterator();
    }

    /**@function IncidentIterator::get_iterator()
     * @brief Returns underlying map iterator.
     * @return map iterator underlying this incident iterator
     * This function is used for the equality comparision function
     * This map iterator points to one of the sub-maps within the
     * edge_map_ map-of-maps. edge_map_ has contains node
     * indexes which point to node indexes which point to edge
     * indexes in a nested map type structure.
     */
    std::map<size_type,size_type>::iterator get_iterator() const{
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

  /**@class Graph::EdgeIterator
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

    /**@function EdgeIterator::EdgeIterator()
     * @brief Constructs an  invalid EdgeIterator
     */
    EdgeIterator() {
    }

    /**@function EdgeIterator::EdgeIterator()
     * @brief Constructs a valid EdgeIterator
     * @param[in]   g   pointer to  graph
     * @param[in]   i   Index of underlying edge
     * The index is only used to check if the EdgeIterator is
     * at its end.
     */
    EdgeIterator(const graph_type* g, int i ) {
        my_graph_ = const_cast<graph_type*>(g);
        if (i == -1) {
            my_map_iterator_ = my_graph_->graph_edges_.end();
        } else {
            my_map_iterator_ = my_graph_->graph_edges_.begin();
        }
    }

    /**@function EdgeIterator::operator*()
     * @brief Dereferencing operator for calling EdgeIterator
     * @return Underlying Edge that this edge iterator points to
     */
    Edge operator*() const {
        return my_graph_->edge(my_map_iterator_->first);
    }

    /**@function EdgeIterator::operator++()
     * @brief Incrementing operator for calling EdgeIterator
     * @return Incremented edge iterator
     * Returns the end EdgeIterator if iterator points to last edge element.
     */
    EdgeIterator& operator++() {
        if (my_map_iterator_ == my_graph_->graph_edges_.end()){
            return *this;
        }
        ++my_map_iterator_;
        return *this;
    }

    /**@function EdgeIterator::operator==()
     * @brief Equality operator for calling EdgeIterator
     * @param[in]   x   EdgeIterator to compare against
     * @return Boolean indicating equality
     */
    bool operator==(const EdgeIterator& x) const {
        return my_map_iterator_ == x.get_iterator();
    }

    /**@function EdgeIterator::get_iterator()
     * @brief Returns underlying map iterator.
     * @return map iterator underlying this edge iterator.
     * This function is used for the equality comparision function
     */
    iterator_type get_iterator() const{
          return this->my_map_iterator_;
    }

   private:
    friend class Graph;
    iterator_type my_map_iterator_;
    graph_type* my_graph_;

  }; /// End of EdgeIterator class

  /**@function Graph::remove_edge()
   * @brief Remove edge pointed to by provided edge iterator
   * @param[in,out]       e_it      edge iterator
   * @pre Provided edge iterator is valid
   * @post Edge pointed to by provided edge iterator is removed
   * @post All edge iterators active before function call are invalid.
   * Does not guarantee that adding future edges
   * will not cause problems since removal does not keep track
   * of a unique id. It just removes the edge's index
   * from the various maps and decrements num_edges.
   * @return Beginning edge iterator of newly modified graph
   *
   * Complexity: Average O(log(num_edges) + log(num_nodes))
   */
  edge_iterator remove_edge(edge_iterator e_it){
    remove_edge(*e_it);
    return edge_begin();
  }

  /**@function Graph::edge_begin()
   * @brief Return the beginning edge iterator for this graph
   * @return Beginning edge iterator
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /**@function Graph::edge_end()
   * @brief Return the end edge iterator for this graph
   * @return End edge iterator
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this, -1);
  }
}; /// End Graph class

#endif // CME212_GRAPH_HPP
