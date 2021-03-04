#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

/////////////////////////////////////////////////////////////
#include <unordered_map>
#include <vector>
/////////////////////////////////////////////////////////////

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

  // HW0: YOUR CODE HERE

  ///////////////////////////////////////////////////////////
  // Predeclare the node struct
  struct internal_node;
  // Predeclare the edge struct
  struct internal_edge;
  // Predeclare the struct needed to have a inverse_edge map
  struct edge_key;
  struct EdgeKeyHasher;
  ///////////////////////////////////////////////////////////

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

  // HW1: YOUR CODE HERE
  ///////////////////////////////////////////////////////////
  using node_value_type = V;
  ///////////////////////////////////////////////////////////

  // HW2: YOUR CODE HERE
  ///////////////////////////////////////////////////////////
  using edge_value_type = E;
  ///////////////////////////////////////////////////////////

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////
    private_size_ = 0;
    size_ = 0;
    num_edges_ = 0;

    nodes_ = std::unordered_map<size_type, internal_node>(); 
    // default constructor, empty map
    assert(nodes_.empty());

    edges_ = std::unordered_map<size_type, internal_edge>();
    // default constructor, empty map
    assert(edges_.empty());

    inverse_edges_ = std::unordered_map<edge_key, size_type, EdgeKeyHasher>();
    // default constructor, empty map
    assert(inverse_edges_.empty());

    adjacencies_=std::unordered_map<size_type, std::vector<size_type>>();
    // default constructor, empty map
    assert(adjacencies_.empty());
    //////////////////////////////////////////////

    // HW2: YOUR CODE HERE -> public/private node index data structure
    //////////////////////////////////////////////
    valid_nodes_=std::vector<size_type>();
    // default constructor, empty vector
    assert(adjacencies_.empty());
    inverse_valid_nodes_=std::unordered_map<size_type, size_type>();
    // default constructor, empty map
    assert(inverse_valid_nodes_.empty());
    //////////////////////////////////////////////

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
      //////////////////////////////////////////////
      // default constructor for all attributes
      //////////////////////////////////////////////      
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      return fetch().point;
      //////////////////////////////////////////////
    }

    /** Return this node's position, as a reference so that we can modify it */
    Point& position() {
      // HW2: YOUR CODE HERE
      //////////////////////////////////////////////
      return fetch().point;
      //////////////////////////////////////////////
    }  

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // HW2: modified because of private/public index
      //////////////////////////////////////////////
      auto search = graph_->inverse_valid_nodes_.find(index_node_);
      if (search == graph_->inverse_valid_nodes_.end()){
        assert(false);
      }
      return search->second;
      //////////////////////////////////////////////
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    //////////////////////////////////////////////

    /** Return modifiable alias this node's value (type node_value_type) */
    node_value_type& value(){
      return fetch().node_value;
    }
    
    /** Return modifiable alias this node's value (type node_value_type) */
    const node_value_type & value() const {
      return fetch().node_value;
    }

    /** Return this node's degree (number of neigbours) */
    size_type degree() const{
      auto search = graph_->adjacencies_.find(this->index_node_);
      if (search == graph_->adjacencies_.end()){
        assert(false);
      }
      return search->second.size();
    }

    /** Return an iterator pointing to the first neigbour of this node
      * (type incident_iterator) */
    incident_iterator edge_begin() const{
      auto search = graph_->adjacencies_.find(this->index_node_);
      if (search == graph_->adjacencies_.end()){
        assert(false);
      }
      return IncidentIterator(graph_,this->index_node_,search->second.begin());
    }

    /** Return an iterator pointing to one past the last neigbour of this node
      * (type incident_iterator) */
    incident_iterator edge_end() const{
      auto search = graph_->adjacencies_.find(this->index_node_);
      if (search == graph_->adjacencies_.end()){
        assert(false);
      }
      return IncidentIterator(graph_, this->index_node_, search->second.end());
    }    
    //////////////////////////////////////////////

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      if ((graph_ == n.graph_) && (index_node_ == n.index_node_)){  
        //check if both nodes point to the same graph and have same index
        return true;
      }
      //////////////////////////////////////////////
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
      //////////////////////////////////////////////
      //assert(graph_ == n.graph_);
      //check if both nodes point to the same graph
      if (graph_ == n.graph_){
        if (index_node_ < n.index_node_){
          return true;
        }
        else{
          return false;
        }
      }
      return (graph_ < n.graph_);
      //////////////////////////////////////////////
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE

    ///////////////////////////////////////////////////
    Graph* graph_; // Associated Graph instance
    size_type index_node_; // Node's unique id number (private index)
    
    /** Private Constructor */
    Node(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_node_(index) {
    }

    /** Helper method to return the appropriate element.
     */
    internal_node& fetch() const {
      auto search = graph_->nodes_.find(index_node_);
      if (search != graph_->nodes_.end()) {
          return search->second;
      }
      assert(false);
    }
    ///////////////////////////////////////////////////

    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects  
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////
    return size_;
    //////////////////////////////////////////////////
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] position The new node's value (optional, initliazed to default
   * if not provided)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
                const node_value_type& node_value = node_value_type()) {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////////

    internal_node new_internal_node;
    new_internal_node.point = position;
    new_internal_node.index_node = private_size_;
    new_internal_node.node_value = node_value;
    nodes_.insert({private_size_,new_internal_node});
    
    std::vector<size_type> neighbours = std::vector<size_type>();
    assert(neighbours.empty());
    adjacencies_.insert({private_size_,neighbours});
    valid_nodes_.push_back(private_size_);
    inverse_valid_nodes_.insert({private_size_,size_});
    size_++;
    // Returns a Node that points to the new element
    return Node(this, private_size_++); //post increment of private size
    /////////////////////////////////////////////////////
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////////
    if (this != n.graph_){
      return false;
    }
    auto search = inverse_valid_nodes_.find(n.index_node_);
    //Done in amortized time O(1)
    if (search != inverse_valid_nodes_.end()) {
        return true;
    }
    return false;
    //////////////////////////////////////////////////////
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE // Modified HW2 because public/private indexes
    //////////////////////////////////////////////////////
    assert((0<=(int)i) && (i<=size_));
    size_type private_i = valid_nodes_[i];
    auto search = nodes_.find(private_i); //Done in amortized time O(1)
    if (search != nodes_.end()) {
      return Node(this, private_i);
    }
    assert(false);
    //////////////////////////////////////////////////////
  }

  // HW2: YOUR CODE HERE
  //////////////////////////////////////////////

  /** Remove node @a _n_ from the graph, returning 1 if @a _n_ was correctly 
   *  removed, 0 if @a _n_ was not in the graph.
   * 
   * @param[in] _n_ the node to remove, passed as const reference
   * @return 1 if @a _n_ was correctly removed, 0 if @a _n_ was not in the graph
   * 
   * @pre @a _n_ is a valid node for this graph
   * @post If @a _n_ is indeed removed new num_nodes() == old num_nodes() - 1
   *       Else new num_nodes() == old num_nodess()
   * @post has_node(n) == false
   * @post new num_edges() == old (num_edges() - n.degree())
   * @post old n is invalid and should never be used again. Undefined 
   *       behavior otherwise

   * May change public node indexing (ie node indexes given by index() method), 
   * but keeps intact private indexing.
   * Remove (and thus invalidate) all edges whose one endpoint is @a _n_
   *
   * Complexity: O(n.degree()) < O(num_nodes()) on average
   */
  size_type remove_node(const Node& n){
    if (not has_node(n)){
      return 0;
    }
    else{
      // ---- Remove all incident edges ---- //
      std::vector<size_type> neighbours = std::vector<size_type>();
      //we first copy the neigbours indexes as while removing edeges we change  
      //the graph structure and this could provoke indefined behavior to
      //iterate directly on the adjecency list
      auto search = adjacencies_.find(n.index_node_); //private index
      if (search != adjacencies_.end()){
        // we iterate the vector of private indexes
        for (auto it=search->second.begin(); it!=search->second.end(); ++it){
          neighbours.push_back(*it);
        }
      }
      else{
        assert(false);
      }
      Node neigh;
      size_type result;
      for (auto it = neighbours.begin(); it != neighbours.end(); ++it){
        //private constructor with private index
        neigh = Node(this, *it);
        result = remove_edge(n, neigh);
        if (result==0){ // The neigbours was in the adjecency so there must be a
        // removal
          assert(false);
        }
      }
      // ---- Remove all incident edges ---- //



      // ---- Remove node from public/private mapping ---- //
      // ---- The node with highest index take the index of deleted node ---- //
      size_type public_i = n.index();
      size_type private_i = n.index_node_;
      
      if (public_i != size_-1){
        // Swap and pop strategy to have O(1) complexity
        // But this change order of public indexes -> cf specification
        std::iter_swap(valid_nodes_.begin() + public_i,
                      valid_nodes_.end() - 1);
        valid_nodes_.pop_back();
        inverse_valid_nodes_.erase(private_i);
        inverse_valid_nodes_[valid_nodes_[public_i]] = public_i;
      }
      else{ // easy case, the node to delete is already the last vector elem
        valid_nodes_.pop_back();
        inverse_valid_nodes_.erase(private_i);
      }
      // ---- The node with highest index take the index of deleted node ---- //
      // ---- Remove node from public/private mapping ---- //
      


      // ---- Decrease the number of nodes  ---- //
      size_--; //private_size_ remains unchanged
      // ---- Decrease the number of nodes ---- //


      return 1;
      
    }
  }

  


  /** Remove the node n pointed by @a _n_it_ from the graph
   * 
   * @param[in] _n_it_ the iteartor pointing to node to remove
   * @return another node_iterator pointing to another valid node of the graph
   *         if there are still any, otherwise undefined behavior
   * 
   * @pre @a _n_it_ is a valid node_iterator (ie it points to a valid node)
   *      Let's call n the valid node it points to.
   * @post If n is indeed removed new num_nodes() == old num_nodes() - 1
   *       Else new num_nodes() == old num_nodess()
   * @post has_node(n) == false
   * @post new num_edges() == old (num_edges() - n.degree())
   * @post old n is invalid and should never be used again. Undefined 
   *       behavior otherwise
   * May change public node indexing (ie node indexes given by index() method), 
   * but keeps intact private indexing.
   * Remove (and thus invalidate) all edges whose one endpoint is n
   *
   * Complexity: O(n.degree()) < O(num_nodes()) on average
   */
  node_iterator  remove_node(node_iterator  n_it){
    Node pointed_node = *n_it;
    remove_node(pointed_node);
    return this->node_begin();
  }
  //////////////////////////////////////////////




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
      //////////////////////////////////////////////
      // default constructor for all attributes
      //////////////////////////////////////////////   
    }


    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      return Node(graph_, index_node_1_);
      //////////////////////////////////////////////
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      return Node(graph_, index_node_2_);
      //////////////////////////////////////////////
    }

    /** Return this edge's index, a number in the range [0, num_nedges()). */
    size_type index() const {
      // HW1: YOUR CODE HERE
      //////////////////////////////////////////////
      return index_edge_;
      //////////////////////////////////////////////
    }

    /** Return this edge's euclidean length */
    double length() const {
      // HW2: YOUR CODE HERE
      //////////////////////////////////////////////
      return norm(node1().position()-node2().position());
      //////////////////////////////////////////////
    }

    /** Return a mofifiable reference to this edge's value */
    edge_value_type& value() {
      // HW2: YOUR CODE HERE
      //////////////////////////////////////////////
      return fetch().edge_value;
      //////////////////////////////////////////////
    }

    /** Return a const reference to this edge's value */
    const edge_value_type& value() const{
      // HW2: YOUR CODE HERE
      //////////////////////////////////////////////
      return fetch().edge_value;
      //////////////////////////////////////////////
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      if ((graph_ == e.graph_) && (index_edge_ == e.index_edge_)){
        return true;
      }
      //////////////////////////////////////////////
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE
      //////////////////////////////////////////////
      //assert(graph_ == e.graph_);
      //check if both edges point to the same graph 
      // -> taken into account in the if below
      
      if ((graph_ == e.graph_)){
        if (index_edge_ < e.index_edge_){
          return true;
        }
        return false;
      }
      return (graph_ < e.graph_);
      //////////////////////////////////////////////
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE

    ///////////////////////////////////////////////////
    Graph* graph_; // Associated Graph instance
    size_type index_edge_; // Edge's unique identification number
    size_type index_node_1_;
    size_type index_node_2_;

    /** Private Constructor -> uses private node indexes*/
    Edge(const Graph* graph, size_type index,
         size_type index_node_1, size_type index_node_2)
        : graph_(const_cast<Graph*>(graph)), index_edge_(index),
          index_node_1_(index_node_1), index_node_2_(index_node_2) {}

    /** Helper method to return the appropriate element.
     */
    internal_edge& fetch() const {
      auto search = graph_->edges_.find(index_edge_);
      if (search != graph_->edges_.end()) {
          return search->second;
      }
      assert(false);
    }

    ///////////////////////////////////////////////////

    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////
    return num_edges_;
    //////////////////////////////////////////////////
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////////
    assert((0<=(int)i) && (i<=num_edges_));
    auto search = edges_.find(i); //Done in amortized time O(1)
    if (search != edges_.end()) {
      size_type index_n1 = search->second.node_1.index_node_;
      size_type index_n2 = search->second.node_2.index_node_;
      return Edge(this, i, index_n1, index_n2);
    }
    assert(false);
    //////////////////////////////////////////////////////
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //////////////////////////////////////////////////////
    if (!(has_node(a)) || !(has_node(b))){ // O(1) complexity
      return false;
    }

    edge_key possible_key;
    
    possible_key.index_node_1=a.index_node_;
    possible_key.index_node_2=b.index_node_;
    auto search = inverse_edges_.find(possible_key); //Done in amortized O(1)
    if (search != inverse_edges_.end()) {
        return true;
    }

    possible_key.index_node_1=b.index_node_;
    possible_key.index_node_2=a.index_node_;
    search = inverse_edges_.find(possible_key); //Done in amortized O(1)
    if (search != inverse_edges_.end()) {
        return true;
    }
    return false;
    //////////////////////////////////////////////////////
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
    ////////////////////////////////////////////////////
    assert((has_node(a)) && (has_node(b))); // O(1) complexity
    // check if a and b are valid nodes for this graph
    
    assert((a.index_node_!=b.index_node_)); 
    // check if a and b are different
    

    edge_key possible_key;
    size_type index_edge;
    edge_value_type edge_value;
    
    possible_key.index_node_1=a.index_node_; //private index
    possible_key.index_node_2=b.index_node_; //private index
    auto search = inverse_edges_.find(possible_key); //Done in amortized O(1)
    
    if (search != inverse_edges_.end()) {
        index_edge = search->second;
        edge_value = edges_[index_edge].edge_value;
        inverse_edges_.erase(possible_key); //Done in amortized O(1)
        edges_.erase(index_edge);//Done in amortized O(1)
    }
    else{
      possible_key.index_node_1=b.index_node_;
      possible_key.index_node_2=a.index_node_;
      search = inverse_edges_.find(possible_key); //Done in amortized O(1)
      if (search != inverse_edges_.end()) {
          // we must flip node_1 and node_2 in this case
          index_edge = search->second;
          edge_value = edges_[index_edge].edge_value;
          inverse_edges_.erase(possible_key); //Done in amortized O(1)
          edges_.erase(index_edge);//Done in amortized O(1)
      }
      else{ // we add a new edges to edges_ and inverse_edges_
          index_edge = num_edges_;
          num_edges_++;
          // we update adjacencies_
          auto s = adjacencies_.find(a.index_node_);
          assert(s != adjacencies_.end());
          s->second.push_back(b.index_node_);
          
          s = adjacencies_.find(b.index_node_);
          assert(s != adjacencies_.end());
          s->second.push_back(a.index_node_);

          edge_value = edge_value_type(); //default value
      }
    }

    possible_key.index_node_1=a.index_node_;
    possible_key.index_node_2=b.index_node_;

    internal_edge new_internal_edge;
    new_internal_edge.node_1 = a;
    new_internal_edge.node_2 = b;
    new_internal_edge.index_edge  = index_edge;
    new_internal_edge.edge_value  = edge_value; // ok careful. it erases 
    // previous edge value if the edge is already in the graph
    
    edges_.insert({index_edge,new_internal_edge}); //Done in amortized O(1)
    inverse_edges_.insert({possible_key,index_edge}); //Done in amortized O(1)
    
    // Returns an Edge that points to the new element
    return Edge(this, index_edge, a.index_node_, b.index_node_);
    ////////////////////////////////////////////////////
  }

  // HW2: YOUR CODE HERE
  ///////////////////////////////////////////////////

  /** Remove edge formed by @a _a_ and @a _b_ from the graph. Let's call e this
   *  edge. Method returning 1 if e was correctly removed, 0 if e was not in 
   *  the graph.
   * @param[in] @a _a_ and @a _b_ the nodes forming the edge e to remove, 
   *            passed as const reference
   * @return 1 if e was correctly removed, 0 if e was not in the graph
   * 
   * @pre @a _a_ and @a _b_ are valid nodes for this graph
   * @post If e is indeed removed new num_edges() == old num_edges() - 1
   *       Else new num_edges() == old num_edges()
   * @post @a _a_ and @a _b_ are valid nodes for this graph
   * @post has_edge(a,b) == false
   * @post new a.degree() = old (a.degree() - 1) if e was indeed removed
   * @post new b.degree() = old (b.degree() - 1) if e was indeed removed
   * @post old a.index() = new a.index() and old b.index() = new b.index()
   * @post e is now invalid and should never be used again. Undefined 
   *       behavior otherwise
   * 
   * Can invalidate edge indexes -- in other words, old edge(e.index()) might not
   * equal new edge(e.index()). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(1) on average.
   */
  size_type  remove_edge(const  Node& a, const  Node & b){
    if (!(has_edge(a, b))){
      return 0;
    }
    else{

      // ---- Remove edge from our record ---- //
      edge_key possible_key;
      size_type index_edge = num_edges_; // to make sure it fails if we 
      // don't find the right index edge
    
      possible_key.index_node_1=a.index_node_; //private index
      possible_key.index_node_2=b.index_node_; //private index
      auto search1 = inverse_edges_.find(possible_key); //Done in amortized O(1)
      if (search1 != inverse_edges_.end()) {
          index_edge = search1->second;
      }
      else{
        possible_key.index_node_1=b.index_node_;
        possible_key.index_node_2=a.index_node_;
        search1 = inverse_edges_.find(possible_key); //Done in amortized O(1)
        if (search1 != inverse_edges_.end()) {
            index_edge = search1->second;
        }
        else{
          assert(false);
        }
      }
      inverse_edges_.erase(possible_key); //Done in amortized O(1)
      edges_.erase(index_edge);//Done in amortized O(1)
      // ---- Remove edge from our record ---- //



      // ---- Make the highest index edge our new index_edge ---- //
      if (index_edge != num_edges_-1){
        size_type index_highest_edge = num_edges_-1;
        size_type n1;
        size_type n2;
        auto search2 = edges_.find(index_highest_edge); //amortized O(1)
        if (search2 != edges_.end()){ 
            internal_edge int_edge = std::move(search2->second);
            n1 = int_edge.node_1.index_node_; //private index
            n2 = int_edge.node_2.index_node_; //private index
            int_edge.index_edge = index_edge;
            edges_.erase(search2);
            edges_.insert({index_edge, std::move(int_edge)});
        }
        else{
          assert(false);
        }
        possible_key.index_node_1=n1;
        possible_key.index_node_2=n2;
        auto search3 = inverse_edges_.find(possible_key); // amortized O(1)
        if (search3 != inverse_edges_.end()) {
            search3->second = index_edge;
        }
        else{
          possible_key.index_node_1=n2;
          possible_key.index_node_2=n1;
          search3 = inverse_edges_.find(possible_key); //Done in amortized O(1)
          if (search3 != inverse_edges_.end()) {
              search3->second = index_edge;
          }
          else{
            assert(false);
          }
        }
      }
      // ---- Make the highest index edge our new index_edge ---- //

    
      // ---- Remove edge from adjecencies list ---- //
      size_type a_index = a.index_node_; //private index
      size_type b_index = b.index_node_; //private index
      auto search4 = adjacencies_.find(a_index); //amortized O(1)
      if (search4 != adjacencies_.end()) {
        // we use erase-remove idiom in O(degree), ok as the graph is sparse
        search4->second.erase(std::remove(search4->second.begin(), 
                                          search4->second.end(),
                                          b_index),
                              search4->second.end());
      }
      else{
        assert(false);
      }
      search4 = adjacencies_.find(b_index); //amortized O(1)
      if (search4 != adjacencies_.end()) {
        // we use erase-remove idiom in O(degree), ok as the graph is sparse
        search4->second.erase(std::remove(search4->second.begin(), 
                                          search4->second.end(),
                                          a_index),
                              search4->second.end());
      }
      else{
        assert(false);
      }
      // ---- Remove edge from adjecencies list ---- //

      // ---- Decrease by 1 the number of edges ---- //
      num_edges_--;
      // ---- Decrease by 1 the number of edges ---- //

      return 1;
    }
  }



  /** Exact similar specification of 
   *  size_type remove_edge(const Node& a, const  Node & b)
   *  with a = e.node1() and b = e.node2()
   */
  size_type  remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }


  /** Similar specification of 
   *  size_type remove_edge(const Edge& e)
   *  with e = *e_it
   * 
   * Only one change:
   * @return another edge_iterator pointing to another valid edge of the graph
   *         if there are still any, otherwise undefinied behavior
   */
  edge_iterator  remove_edge(edge_iterator  e_it){
    Edge pointed_edge = *e_it;
    remove_edge(pointed_edge);
    return this->edge_begin();
  };

  ///////////////////////////////////////////////////

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    ///////////////////////////////////////////////////
    size_=0;
    num_edges_=0;
    nodes_.clear();
    edges_.clear();
    inverse_edges_.clear();
    adjacencies_.clear();
    valid_nodes_.clear();
    inverse_valid_nodes_.clear();
    ///////////////////////////////////////////////////
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
    ///////////////////////////////////////////////////
    
    
    //Increments to the next node in the Node class.
    /** Increment this NodeIterator to the next element.
     *  new NodeIterator will point to the next Node, 
     *  right after the one pointed by old NodeIterator
     */
    NodeIterator& operator++(){
      node_iter++;
      return *this;
    }
    
    //Defines equality between two iterators
    /** Test whether this NodeIterator and @a node_iter2 are equal.
     *
     * Equal NodeIterator have the same graph and point to the same node.
     */
    bool operator==(const NodeIterator& node_iter2) const{
      return (associated_graph == node_iter2.associated_graph) \
          && (node_iter == node_iter2.node_iter);		
    }
    
    //Defines inequality between two iterators
    /** Test whether this NodeIterator and @a node_iter2 are different.
     *
     * Equal NodeIterator have the same graph and point to the same node.
     */
    bool operator!=(const NodeIterator& node_iter2) const{
      return !(*this == node_iter2);			
    }
    
    //Dereference operator
    /** Dereference this NodeIterator, giving the Node it is pointing to
     *
     */
    Node operator*() const{
      size_type index_node = *node_iter;
      return Node(associated_graph, index_node);		
    }
    ///////////////////////////////////////////////////
    
   private:
    friend class Graph;
   // HW1 #2: YOUR CODE HERE
   ///////////////////////////////////////////////////
    const graph_type* associated_graph;
    typename std::vector<size_type>::const_iterator node_iter;
    
    
    /** Private Constructor */
    NodeIterator(const graph_type* associated_graph_, typename \
                 std::vector<size_type>::const_iterator n_it_){
      associated_graph = associated_graph_;
      node_iter = n_it_;
    }
    ///////////////////////////////////////////////////
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  ///////////////////////////////////////////////////
	
  
  /** Return an iterator pointing to the first node of this graph
    * (type node_iterator) */
	node_iterator node_begin() const{
		return NodeIterator(this, valid_nodes_.begin());
	} 
	
  /** Return an iterator pointing to one past the last node of this graph
    * (type node_iterator) */
	node_iterator node_end() const{
		return NodeIterator(this, valid_nodes_.end());
	}
  ///////////////////////////////////////////////////

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
    ///////////////////////////////////////////////////
    


    /** Increment this IncidentIterator to the next element.
     *  new IncidentIterator will point to the next Edge, 
     *  right after the one pointed by old IncidentIterator
     */
    IncidentIterator& operator++(){
      incident_iter++;
      return *this;
    }
    
    /** Test whether this IncidentIterator and @a ii_iter2 are equal.
     *
     * Equal IncidentIterator have the same associated_graph, same incident_iter
     * and point to the same root_node_index.
     */
    bool operator==(const IncidentIterator& ii_iter2) const{
      return (associated_graph == ii_iter2.associated_graph) \
          && (root_node_index == ii_iter2.root_node_index)
          && (incident_iter == ii_iter2.incident_iter);		
    }
    


    /** Test whether this IncidentIterator and @a ii_iter2 are different.
     *
     * Equal IncidentIterator have the same associated_graph, same incident_iter
     * and point to the same root_node_index.
     */
    bool operator!=(const IncidentIterator& ii_iter2) const{
      return !(*this == ii_iter2);			
    }
    


    //Dereference operator
    /** Dereference this IncidentIterator, giving the Edge it is pointing to
     *
     */
    Edge operator*() const{
      size_type index_neighbour_node = *incident_iter;
      edge_key key;
      size_type index_edge;

      key.index_node_1=root_node_index;
      key.index_node_2=index_neighbour_node;
      auto search = associated_graph->inverse_edges_.find(key);
      if (search != associated_graph->inverse_edges_.end()){
        index_edge = search->second;
      }
      else{
        key.index_node_2=root_node_index;
        key.index_node_1=index_neighbour_node;
        search = associated_graph->inverse_edges_.find(key);        
        if (search != associated_graph->inverse_edges_.end()){
          index_edge = search->second;
        }
        else{
          assert(false); //adjacencies was having an inexistant edge
        }
      }
      return Edge(associated_graph, index_edge, 
                  root_node_index, index_neighbour_node);
    }
    ///////////////////////////////////////////////////

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    ///////////////////////////////////////////////////
    const graph_type* associated_graph;
    size_type root_node_index; //private index
    std::vector<size_type>::const_iterator incident_iter;
    

    /** Private Constructor */
    IncidentIterator(const graph_type* associated_graph_, \
                     size_type root_node_index_, \
                     std::vector<size_type>::const_iterator inc_iter_){
      associated_graph = associated_graph_;
      root_node_index = root_node_index_;
      incident_iter = inc_iter_;
    }
    ///////////////////////////////////////////////////    
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
    ///////////////////////////////////////////////////
    
    /** Increment this EdgeIterator to the next element.
     *  new EdgeIterator will point to the next Edge, 
     *  right after the one pointed by old EdgeIterator
     */
    EdgeIterator& operator++(){
      edge_iter++;
      return *this;
    }
    
    /** Test whether this EdgeIterator and @a edge_iter2 are equal.
     *
     * Equal EdgeIterator have the same associated_graph and same edge_iter.
     */
    bool operator==(const EdgeIterator& edge_iter2) const{
      return (associated_graph == edge_iter2.associated_graph) \
          && (edge_iter == edge_iter2.edge_iter);		
    }
    
    /** Test whether this EdgeIterator and @a edge_iter2 are different.
     *
     * Equal EdgeIterator have the same associated_graph and same edge_iter.
     */
    bool operator!=(const EdgeIterator& edge_iter2) const{
      return !(*this == edge_iter2);			
    }
    
    //Dereference operator
    /** Dereference this EdgeIterator, giving the Edge it is pointing to
     *
     */
    Edge operator*() const{
      size_type index_edge = edge_iter->first;
      return associated_graph->edge(index_edge);		
    }
    ///////////////////////////////////////////////////


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    ///////////////////////////////////////////////////
    const graph_type* associated_graph;
    typename \
        std::unordered_map<size_type, internal_edge>::const_iterator edge_iter;


    /** Private Constructor */
    EdgeIterator(const graph_type* associated_graph_, typename \
           std::unordered_map<size_type, internal_edge>::const_iterator ed_it_){
      associated_graph = associated_graph_;
      edge_iter = ed_it_;
    }
    ///////////////////////////////////////////////////
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  ///////////////////////////////////////////////////
	
  /** Return an iterator pointing to the first edge of this graph
  * (type edge_iterator) */
  edge_iterator edge_begin() const{
		return EdgeIterator(this, edges_.begin());
	} 
	
  /** Return an iterator pointing to one past the last edge of this graph
  * (type edge_iterator) */
	edge_iterator edge_end() const{
		return EdgeIterator(this, edges_.end());
	}  
  ///////////////////////////////////////////////////

 private:

  // HW0: YOUR CODE HERE

  ///////////////////////////////////////////////////////////
  // Internal type for set elements
  struct internal_node {
    Point point;   // The Point held by each node
    size_type index_node;  // The unique indentifier for a node
    node_value_type node_value; // The node value
  };

  struct internal_edge {
    node_type node_1;   // The first node pointed by the edge
    node_type node_2;   // The second node pointed by the edge
    size_type index_edge;  // The unique indentifier for an edge
    edge_value_type edge_value; // The edge value
  };

  struct edge_key { //a hashable struct to be the key of inverse_edges_
    size_type index_node_1;   // The first node pointed by the edge
    size_type index_node_2;   // The second node pointed by the edge
    bool operator==(const edge_key &other) const{
      return (index_node_1 == other.index_node_1 
           && index_node_2 == other.index_node_2);
    }
  };
  struct EdgeKeyHasher{ //Hash function for edge_key. inspired from
  //https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
    std::size_t operator()(const edge_key& k) const{
      size_type hash1 = std::hash<size_type>{}(k.index_node_1); 
      size_type hash2 = std::hash<size_type>{}(k.index_node_2); 
      return hash1 ^ hash2;
    } 
  };

  size_type size_; // number of valid nodes 
  size_type private_size_; // total number of nodes 
  size_type num_edges_; // number of edges

  std::unordered_map<size_type, internal_node> nodes_;
  
  std::unordered_map<size_type, internal_edge> edges_;
  std::unordered_map<edge_key, size_type, EdgeKeyHasher> inverse_edges_;

  std::unordered_map<size_type, std::vector<size_type>> adjacencies_;

  ///////////////////////////////////////////////////////////
  // HW2 -> vector indexed by public index to map to private indexes
  std::vector<size_type> valid_nodes_;
  std::unordered_map<size_type, size_type> inverse_valid_nodes_;
  ///////////////////////////////////////////////////////////
  
  

  // Disable copy and assignment of a Graph
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;
  ///////////////////////////////////////////////////////////

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
