#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

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
  //internal node that contains all info on node
  struct internal_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  /** type of template node value and edge value*/
  using node_value_type = V;
  using edge_value_type = E;

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
            : size_(0), num_distinct_edges_(0), internal_nodes_(), internal_edges_()  // init internal_nodes and size_
  {
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
  class Node  {
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
          return fetch().position;
      }

      /** Make node's position modifiable by returning reference*/
      Point& position(){
          return fetch().position;
      }

      /** Return this node's index, a number in the range [0, graph_size). */
      size_type index() const {
          return fetch().index;
      }

      /** Return this node's value, in template variable */
      node_value_type& value() {
          return fetch().node_value;
      }

      /** Const return this node's value, in template variable */
      const node_value_type& value() const {
          const node_value_type node_value = fetch().node_value;
          return node_value;
      }

      /** Return the number of incident edges */
      size_type degree() const {
          return fetch().incident_edges_.size();
      }

      /** Start of the incident iterator. */
      incident_iterator edge_begin() const{
          return IncidentIterator(0, this);
      }

      /** End of the incident iterator. */
      incident_iterator edge_end() const{
          return IncidentIterator(this->degree(), this);
      }

      /** Test whether this node and @a n are equal.
       *
       * Equal nodes have the same graph and the same index.
       */
      bool operator==(const Node &n) const {
          if ((fetch().index == n.index()) && (graph_ == n.graph_)) { // n.graph_ i'm a little unsure about
              return true;
          }
          return false;
      }

      bool operator!=(const Node &n) const {
          if ((fetch().index != n.index()) || (graph_ != n.graph_)) { // n.graph_ i'm a little unsure about
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
      bool operator<(const Node &n) const {
          if (fetch().index < n.index()) {
              return true;
          }
          return false;
      }

      bool operator>(const Node &n) const {
          if (fetch().index > n.index()) {
              return true;
          }
          return false;
      }

  private:
      // Allow Graph to access Node's private member data and functions.
      friend class Graph;

      //pointer to graph this node belongs to
      graph_type* graph_;

      //node's unique id
      size_type idx_;

      /** Private Constructor */
      Node(const graph_type *graph, size_type idx)
              : graph_(const_cast<graph_type *>(graph)), idx_(idx) {
      }

      /** Helper method to return the appropriate element.
       * This uses vector indexing
       */
      internal_node& fetch() const {
          return *(graph_->internal_nodes_.at(graph_->i2u_[idx_]));
      }
    };

      /** Return the number of nodes in the graph.
       *
       * Complexity: O(1).
       */
      size_type size() const {
          return size_;
      }

      /** Synonym for size(). */
      size_type num_nodes() const {
          return size();
      }

      /** Add a node to the graph, returning the added node.
       * @param[in] position The new node's position
       * @param[in] node value, the new node's value
       * @post new num_nodes() == old num_nodes() + 1
       * @post result_node.index() == old num_nodes()
       *
       * Complexity: O(1) amortized operations.
       */
      Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
          //add new node to internal nodes
          internal_nodes_.push_back(new internal_node(size_, position, node_value));
          i2u_.push_back(internal_nodes_.size() - 1);
          size_++;
          const Node return_node(this, size_-1);
          return return_node;
      }

      /** Determine if a Node belongs to this Graph
       * @return True if @a n is currently a Node of this Graph
       *
       * Complexity: O(1).
       */
      bool has_node(const Node &n) const {
          if (n.index() < size_) {
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
          return {this, i};        // proxy to node
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
    Edge() {}

    /** Return this edge's index, a number in the range [0, num_edges). */
    size_type index() const {
      return fetch_index();
    }

    /** Return this edge's value */
    edge_value_type& value(){
        return graph_->edge_values_[uid_];
    }

    /** Return this edge's value */
    const edge_value_type& value() const{
        return graph_->edge_values_[uid_];
    }

    /** Return a node of this Edge */
    Node node1() const {
        return Node(graph_, graph_->internal_nodes_[node1_id_]->index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
        return Node(graph_, graph_->internal_nodes_[node2_id_]->index);
    }

    /** Return the euclidian length of this node */
    double length() const{
        return norm(graph_->internal_nodes_[node1_id_]->position - graph_->internal_nodes_[node2_id_]->position);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((this->graph_ == e.graph_)) {
          if ((this->node1_id_ == e.node1_id_) and (this->node2_id_ == e.node2_id_)) {
              return true;
          }

          if ((this->node2_id_ == e.node1_id_) and (this->node1_id_ == e.node2_id_)) {
              return true;
          }
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->graph_ != e.graph_){
          return this->graph_ < e.graph_;
      }
      else{
          if ((this->node1_id_ == e.node1_id_) && (this->node2_id_ == e.node2_id_)){
              return false;
          }
          else if ((this->node1_id_ == e.node2_id_) && (this->node2_id_ == e.node1_id_)){
              return false;
          }

          size_type max1 = std::max(this->node1_id_, this->node2_id_);
          size_type max2 = std::max(e.node1_id_, e.node2_id_);

          if (max1 == max2){
              size_type min1 = std::min(this->node1_id_, this->node2_id_);
              size_type min2 = std::min(e.node1_id_, e.node2_id_);
              return min1 < min2;
          }
          else{
              return max1 < max2;
          }
      }
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    //pointer to graph this node belongs to
    graph_type* graph_;

    //Edge id for easy lookup of value
    size_type uid_;

    //edges's node1 and node2
    size_type  node1_id_;
    size_type node2_id_;

    /** Private Constructor without creating new uid*/
    Edge(const graph_type* graph, size_type node1, size_type node2, size_type uid)
            : graph_(const_cast<graph_type *>(graph)), uid_(uid), node1_id_(node1), node2_id_(node2)  {
    }

    /** Helper method to return the appropriate edge.
    * This uses map key
    */
    size_type fetch_index() const {
      return this->graph_->internal_edges_.at(Edge(graph_, node1_id_, node2_id_, uid_));
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_distinct_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edges_vec_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return internal_edges_.count(Edge(this, this->i2u_[a.index()], this->i2u_[b.index()], 0))>=1 || internal_edges_.count(Edge(this, this->i2u_[b.index()], this->i2u_[a.index()], 0))>=1;
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
    //find the existing index
    if (this->has_edge(a,b)){
        //if we found a,b we can return it's uid
        try {
            return {a.graph_, this->i2u_[a.index()], this->i2u_[b.index()], this->internal_edges_.at(Edge(a.graph_, this->i2u_[a.index()], this->i2u_[b.index()], 0))};
        }
        //if we found b,a we return it's id here
        catch(std::out_of_range& e){
            return {a.graph_, this->i2u_[b.index()], this->i2u_[a.index()], this->internal_edges_.at(Edge(a.graph_, this->i2u_[b.index()], this->i2u_[a.index()], 0))};
        }
    }
    else{
        //add the edge to our map of internal edges and vector of edges
        internal_edges_.insert({Edge(a.graph_, this->i2u_[a.index()], this->i2u_[b.index()], edge_values_.size()), edges_vec_.size()});
        edges_vec_.push_back(Edge(a.graph_, this->i2u_[a.index()], this->i2u_[b.index()], edge_values_.size()));
        edge_values_.push_back(edge_value_type());
        num_distinct_edges_++;

        //add the edge to the incident list of both nodes
        internal_nodes_[this->i2u_[a.index()]]->incident_edges_.push_back(Edge(a.graph_, this->i2u_[a.index()], this->i2u_[b.index()], edge_values_.size()-1));
        internal_nodes_[this->i2u_[b.index()]]->incident_edges_.push_back(Edge(a.graph_, this->i2u_[b.index()], this->i2u_[a.index()], edge_values_.size()-1));
        return Edge(this, this->i2u_[a.index()], this->i2u_[b.index()], edge_values_.size()-1);        // Proxy to edge
    }
  }

    /** Remove a node and it's incident edges from the graph
   * @pre @a n is a valid node of this graph
   * @return 1 if removal was succesful, 0 otherwise
   * @post post_graph.size() = pre_graph.size() -1
   * @post old_graph.num_edges might not equal to new_graph.num_edges
   * @post post_graph.begin() to post_graph.end() will not loop over @n anymore
   * @post g.node(i).index() == i for all i with 0<=0<=g.num_nodes()
   * @post g.node(n.index()) == n
   *
   * Can invalidate outstanding edge, node and iterator objects
   *
   * Complexity: O(len(@a.incident_edges) *  O(remove_edge))
   */
  size_type remove_node(const Node& n){
      size_type idx = n.index();

      if (idx>size_){
          return 0;
      }

      //remove incident edges
      auto it = n.edge_begin();
      while(it!=n.edge_end()){
          remove_edge(*it);
      }

      //remove index from i2u_ by swapping last element
      size_type last_uid = this->i2u_.back();
      this->i2u_[idx] = last_uid;
      this->i2u_.pop_back();

      //change the swapped node's index
      this->internal_nodes_[last_uid]->index = idx;

      --size_;
      return 1;
  }

  /** Remove a node (see above)*/
  node_iterator remove_node(node_iterator n_it){
      remove_node(*n_it);
      return node_begin();
  }

    /** Remove an edge from the graph
   * @pre @a n1 is a valid node of this graph
   * @pre @a n2 is a valid node of this graph
   * @pre edge(@a n1, @a n2) is an existing edge on the graph
   * @return 1 if removal was succesful, 0 otherwise
   * @post post_graph.num_edges() = pre_num_edges() -1
   * @post @a n1.edge_begin() to n1.edge_begin() will not loop over the removed edge anymore
   * @post @a n2.edge_begin() to n2.edge_begin() will not loop over the removed edge anymore
   *
   * Can invalidate outstanding edge, node objects and iterator objects.
   *
   *
   * Complexity: O(max(@a n1.incident_edges, @a n2.incident_edges) + log(num_edges()))
   */
  size_type remove_edge(const Node& n1, const Node& n2){
      //You cannot remove edges that do not exist
      if (has_edge(n1, n2) == false){
          return 0;
      };

      // Remove edge from incident vector of n1 and n2 O(max_num of incident edges)
      for (auto it1 = n1.edge_begin(); it1!=n1.edge_end(); ++it1){
          if ((*it1) == Edge(this, this->i2u_[n1.index()], this->i2u_[n2.index()], 0)){
              internal_nodes_[this->i2u_[n1.index()]]->incident_edges_.erase(internal_nodes_[this->i2u_[n1.index()]]->incident_edges_.begin() + it1.current_edge);
          }
      }

      for (auto it2 = n2.edge_begin(); it2!=n2.edge_end(); ++it2){
          if ((*it2) == Edge(this, this->i2u_[n2.index()], this->i2u_[n1.index()], 0)){
              internal_nodes_[this->i2u_[n2.index()]]->incident_edges_.erase(internal_nodes_[this->i2u_[n2.index()]]->incident_edges_.begin() + it2.current_edge);
          }
      }

      //need to get the index of the edge that corresponds to these two nodes O(log(num_edges)).
      unsigned int idx = this->internal_edges_.at(Edge(this, this->i2u_[n1.index()], this->i2u_[n2.index()], (size_type) 1));

      // Remove edge from the edge vector (by swapping and replacing) (O(1))
      Edge last_edge = edges_vec_.back();
      edges_vec_.back() = edges_vec_[idx];
      edges_vec_[idx] = last_edge;
      edges_vec_.pop_back();

      // Remove edge value from edge value vector (by swapping and replacing) (O(1))
      edge_value_type last_edge_value = edge_values_.back();
      edge_values_.back() = edge_values_[idx];
      edge_values_[idx] = last_edge_value;
      edge_values_.pop_back();

      // Change index of swapped element in vector (O(log(num_edges)))
      internal_edges_[last_edge] = idx;

      // Remove edge from internal_edges_ map O(log(num_edges))
      internal_edges_.erase(Edge(this, this->i2u_[n1.index()], this->i2u_[n2.index()], -1));

      --num_distinct_edges_;
      return 1;
  }

  /** Remove an edge Edge object*/
  size_type remove_edge(const Edge& edge){
      return remove_edge(edge.node1(), edge.node2());
  }

  /** Remove an edge with iterator*/
  edge_iterator remove_edge(edge_iterator e_it){
      remove_edge((*e_it).node1(), (*e_it).node2());
      return edge_begin();
  }


    /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      //clear nodes
      this->internal_nodes_.clear();
      this->size_ = 0;

      //clear edges
      this->num_distinct_edges_ = 0;
      this->internal_edges_.clear();
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

    /** dereference operator*/
    Node operator* () const
    {
        return(Node(this->graph, iterator_index));
    }

    /** increment (++) operator*/
    NodeIterator operator++(){
        iterator_index++;
        return *this;
    }

    bool operator==(const NodeIterator& other_iterator) const{
        return (iterator_index == other_iterator.iterator_index & graph == other_iterator.graph);
    }

   private:
    friend class Graph;

    int iterator_index;
    const graph_type* graph;

    /** private constructor that can be accessed by graph class */
    NodeIterator(int iterator_index_, const graph_type* graph_)
        :iterator_index(iterator_index_), graph(graph_) {}
  };

  /** get's the first iterator object */
  node_iterator node_begin() const{
      return NodeIterator(0, this);
  }

  /** get's the end (e.g. last+1) iterator object */
  node_iterator node_end() const{
      return NodeIterator(this->size_, this);
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

    //dereference operator
    Edge operator*() const{
        return (node->graph_->internal_nodes_[node->graph_->i2u_[node->index()]]->incident_edges_[current_edge]);
    }

    //increment operator that goes 1 beyond incident edges
    IncidentIterator& operator++(){
        if (current_edge < node->graph_->internal_nodes_[node->graph_->i2u_[node->index()]]->incident_edges_.size()){
            ++current_edge;
        }
        return (*this);
    }

    /** decrement (--) operator*/
    IncidentIterator& operator--(){
      --current_edge;
      return *this;
    }

    //Equal if current edge is the same and iterator over the same node
    bool operator==(const IncidentIterator& other) const{
        return ((current_edge == other.current_edge) & (node == other.node));
    }

   private:
    friend class Graph;
    //pointer to edge
    size_type current_edge; //pointer to edge element within internal_node.incident_edges

    //the node of this iterator
    const Node* node; //contains index (so we can use internal_nodes_[node.index()].incident_edges

    //create a private constructor
    IncidentIterator(int current_edge_, const Node* node_)
        : current_edge(current_edge_), node(node_) {};
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

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** dereference operator*/
    Edge operator* () const
    {
      return(this->graph->edges_vec_[edge_index]);
    }

    /** increment (++) operator*/
    EdgeIterator& operator++(){
      edge_index++;
      return *this;
    }

    bool operator==(const EdgeIterator& other_iterator) const{
      return (edge_index == other_iterator.edge_index & graph == other_iterator.graph);
    }

   private:
    friend class Graph;

    // keeps the index of the edge
    int edge_index;

    // keeps the graph
    const graph_type* graph;

    /** private constructor that can be accessed by graph class */
    EdgeIterator(int edge_index_, const graph_type* graph_)
          :edge_index(edge_index_), graph(graph_) {}
  };

  /** get's the first iterator object */
  edge_iterator edge_begin() const{
      return EdgeIterator(0, this);
  }

  /** get's the end (e.g. last+1) iterator object */
  edge_iterator edge_end() const{
      return EdgeIterator(this->edges_vec_.size(), this);
  }

 private:
  //private object that contains all info on node with constructor
  struct internal_node {
      size_type index;   // index of the node
      Point position;    // position of the node
      node_value_type node_value; // node value

      //vector of edges that incident this node
      std::vector<Edge> incident_edges_;

      internal_node(size_type index_, Point position_, const node_value_type& node_value_ = node_value_type())
                :index(index_), position(position_), node_value(const_cast<node_value_type &>(node_value_)) {
      };
  };

  //number of nodes
  size_type size_;

  //number of distinct edges
  size_type num_distinct_edges_;

  //vector that contains pointers to internal node objects
  std::vector<internal_node*> internal_nodes_;

  //vector that translates idx to uid (with only "active" set of nodes)
  std::vector<size_type> i2u_;

  //map that contains edges as keys and index as values
  std::map<Edge, size_type> internal_edges_;

  //vector of edges for faster edge(i)
  std::vector<Edge> edges_vec_;

  //vector of edge values of edge.uid_
  std::vector<edge_value_type> edge_values_;

};

#endif // CME212_GRAPH_HPP
