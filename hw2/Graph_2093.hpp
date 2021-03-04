#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <set>
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
    typedef V node_value_type;
    typedef E edge_value_type;
    std::vector<node_value_type> n_attri;
    std::vector<edge_value_type> e_attri;
    std::vector<unsigned> node_active2uniq;
    std::vector<unsigned> edge_active2uniq;
  //  std::vector<int> node_uniq2active;
//    std::vector<int> edge_uniq2active;
  private:
    struct nodeinfo{
      Point pos;
      unsigned uniq_id;
      unsigned active_id;
    };
    unsigned n_num=0;
    unsigned e_num=0;
    std::map<unsigned, nodeinfo> nodes;//uniq node id
    std::map<unsigned, unsigned> n1; //<uniq idx edge, uniq idx node>
    std::map<unsigned, unsigned> n2;
  //  std::vector<unsigned> node_active2uniq;
   // std::vector<unsigned> edge_active2uniq;
    std::vector<int> node_uniq2active;
    std::vector<int> edge_uniq2active;
    //<node idx, the set of all edge idx connecting to the node>
    //for all nodes(uniq id), all active incident edge with uniq edge id
    std::map<unsigned, std::set<unsigned>> n_neighbor_e;

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
  class Node: private totally_ordered<Node> {
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


    /** RReturn this node's position. The node position is modifiable.*/
    Point& position(){
      assert(graph_!=nullptr);
      unsigned uniq=(*graph_).node_active2uniq.at(idx_);
      return (*graph_).nodes[uniq].pos;
    }


    /** Return this node's position. */
    const Point& position() const {
      //In the map, get the value associated with the key
      assert(graph_ != nullptr);
      unsigned uniq=(*graph_).node_active2uniq.at(idx_);
      return (*graph_).nodes[uniq].pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(graph_ != nullptr);
      assert(idx_ < (*graph_).node_uniq2active.size());
      return (*graph_).node_uniq2active[idx_];  //idx_ must be within the range, by construction
    }


    /** Return the value attribute of this node */
    node_value_type& value(){
      unsigned uniq=(*graph_).node_active2uniq.at(idx_);
      return (*graph_).n_attri[uniq];
    }

    /** Return the value attribute of this node, which cannot be modified */
    const node_value_type& value() const{
      unsigned uniq=(*graph_).node_active2uniq.at(idx_);
      return (*graph_).n_attri[uniq];
    }

    /**
     * @brief Find the degree of this node
     *
     * @return The degree of the node
     * @post It is the number of all incident edges
     */
    size_type degree() const{
      unsigned uniq=(*graph_).node_active2uniq.at(idx_);
      return (*graph_).n_neighbor_e[uniq].size();
    }

    /**
     * @brief Indicate the start of the incident iterator
     * @return The start of the incident iterator of the node
     *
     * For all these edges e, the e.node1() is the current node.
     */
    incident_iterator edge_begin() const{
      unsigned uniq=(*graph_).node_active2uniq.at(idx_);
      std::set<unsigned> incident=(*graph_).n_neighbor_e[uniq];
      for (std::set<unsigned>::iterator it=incident.begin();
           it!=incident.end(); it++){
        unsigned one=(*graph_).n1[*it];
        unsigned other=(*graph_).n2[*it];
        if (one!=uniq && other==uniq){  //other==uniq must hold
          (*graph_).n1[*it]=uniq;
          (*graph_).n2[*it]=one;
        }
      }
      return IncidentIterator(graph_, (*graph_).n_neighbor_e[uniq].begin());
    }


    /**
     * @brief Indicate the end of the incident iterator
     * @return The end of the incident iterator of the node
     *
     * All incident edges are visited once.
     */
    incident_iterator edge_end() const{
      unsigned uniq=(*graph_).node_active2uniq.at(idx_);
      return IncidentIterator(graph_, (*graph_).n_neighbor_e[uniq].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_==n.graph_ && n.index()==idx_);
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
      return (idx_ < n.index()); // compare nodes on the same graph
    }




   private:
    // Allow Graph to access Node's private member data and functions.
     Graph* graph_;
     unsigned idx_; //active node, user id
     Node(const Graph* g, unsigned index):
       graph_(const_cast<Graph*>(g)),idx_(index){
     }
     friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return n_num;//node_active2uniq.size();
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
                const node_value_type& x=node_value_type()) {
    unsigned next_idx=nodes.size();
    unsigned next_active=node_active2uniq.size();
    node_active2uniq.push_back(next_idx);
    node_uniq2active.push_back(next_active);
    nodes[next_idx].pos=position;
    nodes[next_idx].uniq_id=next_idx;
    nodes[next_idx].active_id=next_active;
    n_attri.push_back(x);
    n_num+=1;
    return Node(this, next_idx);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.index()<node_active2uniq.size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(node_uniq2active[i]<(int)node_active2uniq.size());
    return Node(this,node_active2uniq[i]);
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
      unsigned edge_idx;// active edge, user idx
    /** Construct an invalid Edge. */
    Edge() {
      edge_idx=0;
    }

    edge_value_type& value(){
      unsigned uniq=(*graph_).edge_active2uniq.at(edge_idx);
      return (*graph_).e_attri[uniq];
    }


    const edge_value_type& value() const{
      unsigned uniq=(*graph_).edge_active2uniq.at(edge_idx);
      return (*graph_).e_attri[uniq];
    }


    /** Return a node of this Edge */
    Node node1() const {
      unsigned e_uniq=(*graph_).edge_active2uniq[edge_idx];
      unsigned pos=((*graph_).n1.find(e_uniq)->second);
      unsigned n_active=(*graph_).nodes[pos].active_id;
      return Node(graph_, n_active);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      unsigned e_uniq=(*graph_).edge_active2uniq[edge_idx];
      unsigned pos=((*graph_).n2.find(e_uniq)->second);
      unsigned n_active=(*graph_).nodes[pos].active_id;
      return Node(graph_, n_active);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      Node cn1=this->node1();
      Node cn2=this->node2();
      Node en1=e.node1();
      Node en2=e.node2();
      if (cn1==en1 && cn2==en2){
        return true;
      }else if (cn1==en2 && cn2==en1){
        return true;
      }else{
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (edge_idx < e.edge_idx);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
     Graph* graph_;
     Edge(const Graph* g, unsigned i):
       edge_idx(i), graph_(const_cast<Graph*>(g)){
     }
     friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return e_num; //edge_active2uniq.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i<edge_active2uniq.size());
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    unsigned where=find_edge(a,b);
    if (where==edge_active2uniq.size()){
      return false;
    }else{
      return true;
    }
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
    unsigned where=find_edge(a,b);
    if (where!=edge_active2uniq.size()){
      unsigned uniq_e=edge_active2uniq[where];
      n1[uniq_e]=node_active2uniq[a.index()];
      n2[uniq_e]=node_active2uniq[b.index()];
      return Edge(this, where);
    }else{
      e_num+=1;
      unsigned next_uniq=n1.size();
      n1.insert({next_uniq, node_active2uniq[a.index()]});
      n2.insert({next_uniq, node_active2uniq[b.index()]});
      n_neighbor_e[node_active2uniq[a.index()]].insert(next_uniq);
      n_neighbor_e[node_active2uniq[b.index()]].insert(next_uniq);
      edge_active2uniq.push_back(next_uniq);
      edge_uniq2active.push_back(where);
      Point aa=a.position();
      Point bb=b.position();
      e_attri.push_back(norm(aa-bb));
      return Edge(this, where);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    n1.clear();
    n2.clear();
    edge_uniq2active.clear();
    edge_active2uniq.clear();
    node_uniq2active.clear();
    node_active2uniq.clear();
    n_neighbor_e.clear();
    n_attri.clear();
    e_attri.clear();
    n_num=0;
    e_num=0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    typedef typename std::vector<unsigned>::const_iterator iter_;

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    /**
     * @brief Dereference the current pointer
     * @return The node that the node iterator points to.
     *
     * Complexity O(1)
     */
    Node operator*() const{
      unsigned i = *node_ptr_;  //dereference: unique node id
      return Node(graph_ , graph_->node_uniq2active[i]);
    }

    /**
     * @brief Increment the current pointer
     * @return The pointer pointing to the next node.
     *
     * Complexity O(1)
     */
    NodeIterator& operator++(){
      ++node_ptr_;
      return *this;
    }

    /**
     * @brief Check whether two NodeIterator objects are the same
     * @return True if both variables are the same. False ottherwise.
     *
     * Complexity O(1)
     */
    bool operator==(const NodeIterator& other) const{
      return graph_==other.graph_ && node_ptr_==other.node_ptr_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    iter_ node_ptr_;
    NodeIterator(const Graph* g, iter_ nptr):
      graph_(const_cast<Graph*>(g)), node_ptr_(nptr){
    }
  };

  /**
   * @brief Indicate the start of the node iterator
   * @return The start of the node iterator
   */
  node_iterator node_begin() const{
    return NodeIterator(this, node_active2uniq.cbegin());
  }

  /**
   * @brief Indicate the end of the node iterator
   * @return The end of the node iterator
   *
   * All nodes are visited once.
   */
  node_iterator node_end() const{
    return NodeIterator(this, node_active2uniq.cend());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    typedef typename std::set<unsigned>::const_iterator it_;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /**
     * @brief Dereference the current pointer
     * @return The incident edge that the pointer points to
     *
     * Complexity O(1)
     */
    Edge operator*() const{
      unsigned i=*incident_ptr_;
      return Edge(graph_, graph_->edge_uniq2active[i]);
    }


    /**
     * @brief Increment the pointer
     * @return The pointer that points to the next incident edge
     *
     * Complexity O(1)
     */
    IncidentIterator& operator++(){
      ++incident_ptr_;
      return *this;
    }


    /**
     * @brief Check if two IncidentIterator objects are the same
     * @return True if both variables are the same. False otherwise.
     *
     * Complexity O(1)
     */
    bool operator==(const IncidentIterator& other) const{
      return graph_==other.graph_ && incident_ptr_==other.incident_ptr_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    it_ incident_ptr_;
    IncidentIterator(const Graph* g, it_ ptr):
      graph_(const_cast<Graph*>(g)), incident_ptr_(ptr){
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    typedef typename std::vector<unsigned>:: const_iterator iter_;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /**
     * @brief Dereference the pointer
     * @return The edge that the pointer points to
     *
     * Complexity O(1)
     */
    Edge operator*() const{
      unsigned i = *ptr_;
      return Edge(graph_, graph_->edge_uniq2active[i]);
    }

    /**
     * @brief Increment the pointer
     * @return The pointer pointing to the next edge
     *
     * Complexity O(1)
     */
    EdgeIterator& operator++(){
      ++ptr_;
      return *this;
    }

    /**
     * @brief Check if two EdgeIterator objects are the same
     * @return True if all variables are the same. False otherwise.
     *
     * Complexity O(1)
     */
    bool operator==(const EdgeIterator& other) const{
      return graph_==other.graph_
          && ptr_==other.ptr_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    iter_ ptr_;
    EdgeIterator(const Graph* g, iter_ e):
      graph_(const_cast<Graph*>(g)), ptr_(e){
    }
  };

  /**
   * @brief Indicate the start of the edge iterator
   * @return The start of the edge iterator
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, edge_active2uniq.cbegin());
  }

  /**
   * @brief Indicate the end of the edge iterator
   * @return The end of the edge iterator
   *
   * All edges are visited once.
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this, edge_active2uniq.cend());
  }


  /**
   * @brief Remove node n from the graph.
   * @param n Node that will be removed
   * @return unsigned int, 1 means the node is removed successfully
   *         If the node is invalid, return 0.
   * @pre The n.index() should be between 0 and node_active2uniq().size()
   *      node_active2uniq[n.index()] should be between 0 and nodes.size()
   *      node_uniq2active[node_active2uniq[idx]]==idx
   * @post The size of node_active2uniq() decreases by 1
   *       Preconditions still hold.
   *
   * All adjacent edges are removed. Then the node is removed.
   * Complexity O(num_node)
   */
  size_type remove_node(const Node& n){
    if (node_uniq2active[node_active2uniq[n.index()]]==-1 
            || n.index()>=node_active2uniq.size()){
      return 0;
    }
    unsigned i=0;
    for (auto it=n.edge_begin(); it!=n.edge_end(); ++it){
      i=remove_edge(*it);
    }
    (void)i;
    unsigned uniq_idx=node_active2uniq[n.index()];
    unsigned active_idx=n.index();
    node_uniq2active[uniq_idx]=-1;
    unsigned old_size=node_active2uniq.size();
    if(active_idx+1<node_active2uniq.size()){
      node_active2uniq[active_idx]=node_active2uniq.back();
      node_active2uniq.pop_back();
    }else{
      node_active2uniq.pop_back();
    }
    // pop_back() does not work?
    unsigned new_size=node_active2uniq.size();
    assert(old_size==new_size+1);
    nodes[uniq_idx].active_id=-1;
    n_num-=1;
    return 1;
  }


  /**
   * @brief Remove node by calling remove_node(Node&);
   * @param n_it A pointer that points to the node that will be removed
   * @return 1 if the removal is successful. 0 otherwise.
   * 
   * Call remove_node(*n_it) since *n_it is a node.
   * Complexity O(num_node()).
   */
  node_iterator remove_node(node_iterator n_it){
    return remove_node(*n_it);
  }


  /**
   * @brief Remove edge connectted by n1, n2
   * @param n1 One end of the edge
   * @param n2 The other end of the edge
   * @return 1 if the removal is successful. O otherwise.
   * @pre If the edge (n1, n2) does not exist, return 0.
   * @post The number of total edges reduces by 1
   *
   * Complexity O(num_edge())
   */
   size_type remove_edge(const Node& n1, const Node& n2){
     unsigned where=find_edge(n1, n2);
     if(where==edge_active2uniq.size()){
       return 0;
     }
     unsigned uniq_idx=edge_active2uniq[where];
     edge_uniq2active[uniq_idx]=-1;
     if (where+1<edge_active2uniq.size()){
       edge_active2uniq[where]=edge_active2uniq.back();
       edge_active2uniq.pop_back();
     }else{
       edge_active2uniq.pop_back();
     }
     n_neighbor_e[node_active2uniq[n1.index()]].erase(uniq_idx);
     n_neighbor_e[node_active2uniq[n2.index()]].erase(uniq_idx);
     e_num-=1;
     return 1;
   }
   
   
  /**
   * @brief Remove edge
   * @param e The edge object that we want to remove
   * @return 1 is the removal is successful. 0 otherwise.
   *
   * Get the two nodes that connect to this edge. 
   * Then call remove_edge(n1, n2)
   * Complexity O(num_edge())
   */
   size_type remove_edge(const Edge& e){
     Node n1=e.node1();
     Node n2=e.node2();
     return remove_edge(n1, n2);
   }

  /**
   * @brief Remove edge
   * @param e_it An edge_iterator that points to the edge that should be removed
   * @return 1 if the removal is successful. 0 otherwise.
   *
   * Derefernce the iterator we get the Edge object e. Then call remove_edge(e)
   * Complexity O(num_edge())
   */
   edge_iterator remove_edge(edge_iterator e_it){
     Edge e=*e_it;
     return remove_edge(e);
   }





 private:
   unsigned find_edge(const Node& a, const Node& b) const{
     auto ita=n_neighbor_e.find(node_active2uniq[a.index()]);
     auto itb=n_neighbor_e.find(node_active2uniq[b.index()]);
     if (ita==n_neighbor_e.end() || itb==n_neighbor_e.end()){
       return edge_active2uniq.size();
     }
     std::set<unsigned> a_e=ita->second;
     std::set<unsigned> b_e=itb->second;
     std::vector<unsigned> common;
     std::set_intersection(a_e.begin(), a_e.end(),
       b_e.begin(), b_e.end(), std::inserter(common, common.begin()));
     if (common.size()==0){
       return edge_active2uniq.size();
     }
     //assert(common.size()==1);//no repeated edges by construction
     return edge_uniq2active[common[0]];
   }
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
