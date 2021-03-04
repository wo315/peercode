#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <forward_list>
#include <list>
#include <tuple>
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
template <typename V, typename E>
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
  using edge_value_type = E;

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
     *   x = graph.node(0);  this node is valid
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() { // a node is invalid if either it has no parent graph, or its
             // parent graph does not recognize the node's id.
      this->node_parent_graph_ = nullptr; 
    }

    // Public Copy construtor: whenever we copy a node, we keep track of the
    // copy's address with the list node_roster.  This way, we can invalidate
    // each copy when we clear the graph through setting the parent_graph
    // pointer to null
    Node(const Node &n) {
        (n.node_parent_graph_)-> node_roster.emplace_front(this);

        this->node_parent_graph_ = n.node_parent_graph_;
        this->node_id_ = n.node_id_;
    }

    /** Return this node's position. */
    const Point& position() const {  
      return std::get<0>((this->node_parent_graph_->node_map).at(node_id_));
    }

    // function to change position of node (automatically updates hash)
    Point& position() {
      return std::get<0>((this->node_parent_graph_->node_map).at(node_id_));
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return std::get<1>((this->node_parent_graph_->node_map).at(node_id_));
    }

    // method used to obtain and set the value of a node
    const node_value_type& value() const{
        return std::get<2>((this->node_parent_graph_->node_map).at(node_id_));
    }

    // we could use const_cast here to avoid repeating ourselves, but this
    // would actually lead to typing more code than repeating the above line
    node_value_type& value(){
        return std::get<2>((this->node_parent_graph_->node_map).at(node_id_));
    }



    size_type degree() const{
        return std::distance(node_parent_graph_->adj_list[index()].begin(), 
               node_parent_graph_->adj_list[index()].end());}

    incident_iterator edge_begin() const{
    return incident_iterator(this, 
           node_parent_graph_->adj_list[index()].begin());}

    incident_iterator edge_end() const{
    return incident_iterator(this, 
           node_parent_graph_->adj_list[index()].end());}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if ((this->node_parent_graph_ == n.node_parent_graph_) && 
          (index() == n.index())){
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
      // we will first compare by the address of their parent graph; if
      // they belong to the same graph, they are then compared by their id

      if (this->node_parent_graph_ < n.node_parent_graph_){
          return true;
      }
      else if(this->node_parent_graph_ == n.node_parent_graph_){
          return (node_id_ < n.node_id_);
      }
      return false;
    }


    bool valid() const {
      if (this->node_parent_graph_ != nullptr){
          return ((this->node_parent_graph_)->node_map.count(node_id_) > 0);
      }
      else{
         return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;


    graph_type* node_parent_graph_; // pointer to the node's parent graph

    size_type node_id_; // unique id of node

    //private constructor
    Node(const graph_type* parent_graph, size_type id)
        : node_parent_graph_(const_cast<graph_type*>(parent_graph)),
          node_id_(id)
          {}
  };



  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // the graph's node map is defined by the injective mapping 
    // node_index -> node_position.  The number of nodes is therefore
    // equal to the size of the map.  The size of an unordered_map
    // can be retrieved in O(1) time.
    return node_map.size();
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
    const node_value_type& val = node_value_type()) {
    //node_map is a hash map, so adding a term is O(1) in complexity
    node_map.emplace(node_count, std::make_tuple(position, num_nodes(), val));

    node_ids.push_back(node_count);
    ++node_count;

    adj_list.push_back(std::list< size_type >{}); // set up adj_list
    //node_vals.push_back(val);
    
    // return a node that holds all the relevant information
    return Node(this, node_count - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.node_parent_graph_ == this &&
        (n.node_parent_graph_)->size() > n.index()){
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
    assert(i < size());
    return Node(this, node_ids[i]);
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
      this->edge_parent_graph_ = nullptr; // invalid edges have no parent graph
    }

    // public copy constructor
    Edge(const Edge &e) {
        (e.edge_parent_graph_)->
        edge_roster.emplace_front(this);

        this->edge_parent_graph_ = e.edge_parent_graph_;
        this->edge_id_ = e.edge_id_;
    }

    double length() const {return norm(node1().position() - 
                           node2().position());}

    // when an edge connecting nodes n1 and n2 is added to a graph, the pair
    // (n1, n2) as recorded in edge_map always satisfies n1 < n2.  However,
    // if an edge is constructed by Edge(a, b), we use the parameter
    // inv_node_indx to determine which nodes to return for node1 and node2.


    // add specpp
    /** Return the index of this Edge */
    size_type index() const {
      return 
      std::get<3>((this->edge_parent_graph_->edge_map.at(this->edge_id_)));
    }

    /** Return a node of this Edge */
    Node node1() const {
      if (this-> inv_node_indx == true){
      return Node(this->edge_parent_graph_, 
      std::get<1>((this->edge_parent_graph_)->edge_map.at(this->edge_id_)));
      }
      else{
      return Node(this->edge_parent_graph_, 
      std::get<0>((this->edge_parent_graph_)->edge_map.at(this->edge_id_)));
      } 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (this-> inv_node_indx == true){
      return Node(this->edge_parent_graph_, 
      std::get<0>((this->edge_parent_graph_)->edge_map.at(this->edge_id_)));
      }else{
      return Node(this->edge_parent_graph_, 
      std::get<1>((this->edge_parent_graph_)->edge_map.at(this->edge_id_)));
    }}

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        if ((this->edge_parent_graph_ == e.edge_parent_graph_) && 
            (this->edge_id_ == e.edge_id_)){
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
      if (this->edge_parent_graph_ < e.edge_parent_graph_){
          return true;
      }
      else if(this->edge_parent_graph_ == e.edge_parent_graph_){
          return (this->edge_id_ < e.edge_id_);
      }
      return false;
    }

    // methods used to obtain and set the value of an edge
    const edge_value_type& value() const{
        return std::get<2>((this->edge_parent_graph_->edge_map).at(edge_id_));
    }

    edge_value_type& value(){
        return std::get<2>((this->edge_parent_graph_->edge_map).at(edge_id_));
    }

    bool valid() const {
      if (this->edge_parent_graph_ != nullptr){
          return ((this->edge_parent_graph_)->edge_map.count(edge_id_) > 0);
      }
      else{
         return false;
      }
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // The Edge proxy class will use a pointer to the parent graph along
    // with an unsigned integer representing the edge index in order to 
    // construct an interface with the edge of a graph
    graph_type* edge_parent_graph_;
    size_type edge_id_;
    bool inv_node_indx;

    //unsigned edge_fingerprint;

    

    // private constructor

    Edge(const graph_type* parent_graph, size_type ei, bool inv = false)
    : edge_parent_graph_(const_cast<graph_type*>(parent_graph)),
      edge_id_(ei), inv_node_indx(inv)
      {}
  
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_map.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return Edge(this, i, false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // the hash map inv_edge_map stores already sorted nodes;
    // complexity is therefore O(1)
  if (inv_edge_map.count(std::vector<size_type>{// see if container is nonempty
                         std::min(a,b).node_id_,
                         std::max(a,b).node_id_}) > 0){
        
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
  Edge add_edge(const Node& a, const Node& b, 
    const edge_value_type& val = edge_value_type()) {
    
    assert(-(a == b) + 1);                // nodes must be distinct
    if (has_edge(a, b)){                  // check if edge is already present
        return Edge(this, inv_edge_map.at(std::vector<size_type>{
                                          std::min(a,b).node_id_,
                                          std::max(a,b).node_id_}), (a>b));
        }
    else{                          // if the edge is not present, we add it
 
        adj_list[a.index()].emplace_front(b.node_id_);
        adj_list[b.index()].emplace_front(a.node_id_);

        edge_ids.push_back(edge_count);
        ++edge_count;

        edge_map.emplace(edge_count - 1, 
            std::make_tuple(std::min(a,b).node_id_, 
                            std::max(a,b).node_id_, 
                            val, num_edges()));

        inv_edge_map.emplace(std::vector<size_type>{
                             std::min(a,b).node_id_,
                             std::max(a,b).node_id_}, 
                                     edge_count - 1);
        
        return Edge(this, edge_count - 1, (a>b)); // return the added edge
        }                         
  }

  //REMOVAL FUNCTIONS//

  /** Remove the given node from the graph, assuming it is valid.
   * @param [in] _n_ is a node object which may belong to the given graph
   * @return "1" if node was valid and removed; "0" if node is not valid
   *
   * @post if _n_ is a valid node, both its incident edged and the node itself
   *  are removed from the graph.
   * 
   * @post INVALIDATED OBJECTS:
   *  Any Node or Edge objects providing interfaces for the removed node
   *    or edges are invalidated, because the graph no longer recognizes the id
   *    associated with each object.
   *  If _n_ is removed, all edge_iterator, incident_iterator, and
   *    node_iterator objects that were pointing to either the removed node
   *    or its incident edges, or are equal to the end() iterator of its
   *    respective type are invalidated.
   *
   * Complexity: O(num_nodes())
   */
  size_type remove_node(const Node& n){

    if (n.valid()){
    // remove all edges from this node //

    // first remove edges from its neighbors
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      const Node& out_node = (*it).node2();
      // search the incident edges of this node for the one connecting to n
      for (auto oit = out_node.edge_begin(); 
                oit != out_node.edge_end(); ++oit){
        if((*oit).node2().node_id_ == n.node_id_){

	  // record index and id of edge to be deleted
          size_type rm_edge_id_ = (*oit).edge_id_; 
          size_type rm_edge_idx_ = Edge(this, rm_edge_id_).index(); 

          // remove from inv_edge_map (O(1) operations)
          inv_edge_map.erase(std::vector<size_type>{
              std::min((*oit).node1().node_id_, (*oit).node2().node_id_),
              std::max((*oit).node1().node_id_, (*oit).node2().node_id_),});

          // remove from adj list entry 
          // (O(1) operations because only one entry in the list is erased)
          adj_list[out_node.index()].erase(oit.neighbor_ptr_);
                        
          // swap and pop from edge_ids vector: O(1); note that if the index
          // is already at the end, this simply removes the last element
          std::swap(edge_ids[rm_edge_idx_], edge_ids[num_edges() - 1]);
          edge_ids.pop_back();

          // rename edge index assuming a nontrivial swap occurred
          if (rm_edge_idx_ != num_edges() - 1){
            std::get<3>(edge_map.at(edge_ids[rm_edge_idx_])) = rm_edge_idx_;
          }
               
          edge_map.erase(rm_edge_id_);

          break; // only one edge to find
        }
      }
    }
    // So far we've searched over at most O(num_edges) edges and conducted at
    // most O(1) operations per edge.  Since the graph is sparse, we have
    // O(num_edges) == O(num_nodes).  Therefore we have conducted O(num_nodes)
    // operations at this point.

    // Now we remove the node from the adj list (O(1) operations in total):
    // first we see whether the node to be removed has index num_nodes() -1
    // if so, there is no reindexing

    if (n.index() == num_nodes() - 1){
        adj_list.pop_back();
        node_ids.pop_back();
        node_map.erase(n.node_id_);
    }
    else{
      // obtain node index (O(1) operations)
      size_type rm_idx = n.index();
    
      // swap-and-pop from adj_list and node_ids (O(degree()) operations to
      // erase the component list in the adj_list vector)
      std::swap(adj_list[rm_idx], adj_list[num_nodes() - 1]);
      adj_list.pop_back();

      // remove entry from node_ids (O(1) operations) 
      // the node whose index is num_nodes() - 1 is changed to the deleted
      // node's index, just like the adj list
      std::swap(node_ids[rm_idx], node_ids[num_nodes() - 1]);
      node_ids.pop_back();

      // rename the swapped node's index in node_map (this changes num_nodes())
      node_map.erase(n.node_id_);
      // update the index
      std::get<1>(node_map.at(node_ids[rm_idx])) = rm_idx; 
    }
      return 1;
  }
  else{ return 0; }  
  }


  /** Remove the referent node from the graph, if it is valid 
   * @param[in] _n_it_ node_iterator object pointing to the node to be removed
   * @pre _n_it_ != graph.node_end()
   *
   * @post == @post remove_node (const Node& ) function defined above
   * @return the node_iterator pointing to the next node after the deleted node
   */
  node_iterator remove_node ( node_iterator n_it ) {
    node_iterator next_iter = node_iterator(this, std::next(n_it.node_ptr_));
    remove_node(*n_it);
    return next_iter;
  }

  //EDGE_REMOVAL

  /** Remove the given edge from the graph, assuming it is valid.
   * @param [in] _n1_, _n2_ are node objects which may belong to the given 
   * graph, and may have an edge between them
   * @return "1" if edge was valid and removed; "0" if edge is not valid
   *
   * @post if _n1_, _n2_ has a valid edge, it is removed from the graph
   * 
   * @post INVALIDATED OBJECTS:
   *  Any Edge objects providing interfaces for the removed edges are
   *    invalidated, because the graph no longer recognizes the id
   *    associated with this edge.
   *  If the edge is removed, all edge_iterator or incident_iterator objects
   *    that were pointing to the removed edge, or are equal to the end()
   *    iterator of its respective type are invalidated.
   *
   * Complexity: O(degree())
   */
  size_type remove_edge(const Node& n1, const Node& n2){
    if (n1.valid() && n2.valid()){
    if (has_edge(n1, n2)) {

      // obtain a reference to the edge id and index
      size_type& rm_edge_id_ = inv_edge_map.at(
                  std::vector<size_type>{std::min(n1.node_id_, n2.node_id_),
                                         std::max(n1.node_id_, n2.node_id_)});
      size_type rm_edge_idx_ = Edge(this, rm_edge_id_).index();

      // remove from edge_ids vector: O(1)
      std::swap(edge_ids[rm_edge_idx_], edge_ids[num_edges() - 1]);
      edge_ids.pop_back();

      // rename swapped edge index
      std::get<3>(edge_map.at(edge_ids[rm_edge_idx_])) = rm_edge_idx_;
              
      // delete from adjacency list: O(degree()) search
      for (auto it = n1.edge_begin(); it != n1.edge_end(); ++it){
        if((*it).node2().node_id_ == n2.node_id_){
          adj_list[n1.index()].erase(it.neighbor_ptr_);
          break;
        }
  
      }
      // see where second node is connected to the first in the adj list
      // O(degree()) search
      for (auto it = n2.edge_begin(); it != n2.edge_end(); ++it){
        if((*it).node2().node_id_ == n1.node_id_){
          adj_list[n2.index()].erase(it.neighbor_ptr_); // O(1) operations
          break;
        }
      }

      // delete from edge_map hashmap
      edge_map.erase(rm_edge_id_);

      // delete from inverse hashmap
      inv_edge_map.erase(std::vector<size_type>{
                         std::min(n1.node_id_,  n2.node_id_),
                         std::max(n1.node_id_, n2.node_id_)});

      // we have iterated over at most O(degree()) edges and conducted O(1)
      // operations for each iteration.  Therefore, removing a valid edge
      // is O(degree()) operations.
      return 1;
    }
    }
    return 0;
  }

  /** Remove the given edge from the graph, assuming it is valid.
   * @param [in] _e_ is an Edge object possibly referring to a valid edge
   * @return "1" if edge was valid and removed; "0" if edge is not valid
   *
   * @post if _e_ is a valid edge, it is removed from the graph
   * @post == @post remove_edge(const Node& , const Node&) defined above
   *
   * Complexity: O(degree())
   */
  size_type remove_edge(const Edge& e){
    return remove_edge(e.node1(), e.node2());
  }

  /** Remove the given valid edge from the graph.
   * @param [in] edge iterator object _e_it_.  Note that _e_it_ points to a valid
   *   edge by construction
   * @pre _e_it_ != graph.edge_end() 
   * @return next iterator from that of the removed edge
   *
   * @post == @post remove_edge(const Node& , const Node&) defined above
   *
   * Complexity: O(degree())
   */
  edge_iterator remove_edge ( edge_iterator e_it ){
    edge_iterator next_iter = edge_iterator(this, std::next(e_it.edge_ptr_));
    remove_edge(*e_it);
    return next_iter; // return the next iterator
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {

    // we first clear out both node_map, edge_map and inv_edge_map
    node_map.clear();
    edge_map.clear();
    inv_edge_map.clear();

    // orphan copy-constructed nodes and edges
    for (auto const& i : node_roster) {
    (*i).node_parent_graph_ = nullptr;
    }

    node_roster.clear();

    for (auto const& j : edge_roster) {
    (*j).edge_parent_graph_ = nullptr;
    }
  
    edge_roster.clear();

    adj_list.clear();
    // note that node_count and 
  }


  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
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


    Node operator*() const{return Node(nodeiter_parent_graph_, 
                          *node_ptr_);}

    NodeIterator& operator++(){
 
        if (node_ptr_ != nodeiter_parent_graph_->node_ids.end()){
            ++node_ptr_;
            }
        return *this;}

    bool operator==(const node_iterator& niter2) const{
                    return (niter2.node_ptr_ == node_ptr_);}

   private:
    friend class Graph;
    
    graph_type* nodeiter_parent_graph_;

    std::vector<size_type >::const_iterator node_ptr_;

    // private constructor



    NodeIterator(const graph_type* parent_graph, 
std::vector<size_type >::const_iterator node_ptr)
          : nodeiter_parent_graph_(const_cast<graph_type*>(parent_graph)),
            node_ptr_(node_ptr)
            {}

  };

  node_iterator node_begin() const {return NodeIterator(this, 
                                                          node_ids.begin());}
  node_iterator node_end() const {return NodeIterator(this, node_ids.end());}

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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

    Edge operator*() const {
    // we just need to check whether the index of the host node is less
    // than its neighbor: if so, we call our Edge constructor with the 
    // 'inv_node_indx' set to true.
      if (parent_node_->node_id_ > *neighbor_ptr_){
            
        return Edge(parent_node_->node_parent_graph_,
               parent_node_->node_parent_graph_->inv_edge_map.at
               (std::vector<size_type>{
               *neighbor_ptr_, parent_node_->node_id_}), true);
      }
      else{
        return Edge(parent_node_->node_parent_graph_,
               parent_node_->node_parent_graph_->inv_edge_map.at
               (std::vector<size_type>{
               parent_node_->node_id_, *neighbor_ptr_}));
      }
    }

    incident_iterator& operator++(){ ++neighbor_ptr_; return *this; }


    bool operator==(const IncidentIterator& iit) const {
      return neighbor_ptr_ == iit.neighbor_ptr_;}

   

   private:
    friend class Graph;
    // private constructor
    IncidentIterator(const Node* n, 
        //std::forward_list<size_type>::const_iterator it) :
        std::list<size_type>::const_iterator it) :
        parent_node_(n), neighbor_ptr_(it) {}

    const Node* parent_node_;

    // pointer to the node's neighbors as listed in the adj_list
    std::list<size_type>::const_iterator neighbor_ptr_;
  };



  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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

    Edge operator*() const{return Edge(edgeiter_parent_graph_,
                           *edge_ptr_);}

    EdgeIterator& operator++(){
        if (edge_ptr_ != edgeiter_parent_graph_->edge_ids.end()){
            ++edge_ptr_;
        }
        return *this;}

    bool operator==(const edge_iterator& eiter2) const{
                    return (eiter2.edge_ptr_ == edge_ptr_);}

   private:

    friend class Graph;

    graph_type* edgeiter_parent_graph_;

    //typename std::map< size_type, 
    //    std::tuple <size_type, size_type, edge_value_type> >::const_iterator
    //    edge_ptr_;

    std::vector<size_type>::const_iterator edge_ptr_;

    // private constructor
    EdgeIterator(const graph_type* parent_graph,
      std::vector<size_type>::const_iterator edge_ptr)
    : edgeiter_parent_graph_(const_cast<graph_type*>(parent_graph)),
      edge_ptr_(edge_ptr)
      {}
  };

  edge_iterator edge_begin() const {return EdgeIterator(this, 
                                                        edge_ids.begin());}
  edge_iterator edge_end() const {return EdgeIterator(this, edge_ids.end());}
 




  // PRIVATE GRAPH ATTRIBUTES //
  private:

  std::vector< unsigned > node_ids;
  unsigned node_count = 0;

  std::vector< unsigned > edge_ids;
  unsigned edge_count = 0;

  std::vector< std::list< size_type > > adj_list;

  // our map class stores node and edge information with three hash maps:

  std::unordered_map< size_type, 
  std::tuple< Point, size_type, node_value_type > > node_map;
  // maps node_id_-> (position, index, value)

  std::unordered_map < size_type, 
  std::tuple < size_type, size_type, edge_value_type, size_type > > edge_map;
  // maps edge_index -> (node1_id, node2_id, val, index) where 
  // node1_id < node2_id

  // given an edge e, we will use another hash function so that accessing the
  // nodes spanned by e is O(1) in complexity.  Following Lippman (p. 874), we
  // will define a hash map <node1_index, node2_index> -> edge_index.  Since
  // there is no native hash function for vector<size_type>, we will create
  // our own

  struct inv_edge_hasher {
    size_type operator()(const std::vector<size_type> &v) const {
       return std::hash< size_type>()(v[0])^(std::hash< size_type>()(v[1])<<1);
    }
  };

  std::unordered_map< const std::vector<size_type>, size_type, inv_edge_hasher>
                                                                 inv_edge_map;

  // forward lists maintain the addresses of nodes/edges created from
  // our copy constructor, so that their parent_graph pointer attributes
  // can be set to nullptr when the graph is cleared, thereby invalidating
  // these objects.
   
  std::forward_list<Node*> node_roster;
  std::forward_list<Edge*> edge_roster;


};

#endif // CME212_GRAPH_HPP
