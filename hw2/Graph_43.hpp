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
 template <typename V,typename E>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  

  //proxy design
  struct NodeProxy;
  struct EdgeProxy;
  struct NeghProxy;//proxy desgin for neghiborhood structure

  //the list of all nodes and all size
  std::vector<NodeProxy> node_list;
  std::vector<EdgeProxy> edge_list;

  //the map from node_index to the position
  //this is used for removing the node
  //I'm not sure will hw1 need us to implement this
  std::vector<unsigned> node_index_list;
  std::vector<unsigned> edge_index_list;

  //store the size of the graph
  unsigned node_num; //node num should be node_index_list.size()
  unsigned edge_num;


  

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
  using node_value_type = V;
  using edge_value_type = E;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    node_num = 0;
    edge_num = 0;
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
  class Node:private totally_ordered<Node> {
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
    }

    /** Return this node's position. */
    Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_pointer->node_list[node_index].P;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_index;
    }

    //HW1: YOUR CODE HERE
    /**implementation of value(), degree(),  **/

    node_value_type& value(){
      return graph_pointer->node_list[node_index].value;
    }

    const node_value_type& value() const{
      return graph_pointer->node_list[node_index].value;
    }

    size_type degree() const{
      return graph_pointer->node_list[node_index].neighborhood_list.size();
    }



    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    incident_iterator edge_begin()const{
      return IncidentIterator(graph_pointer,this,0);
    }

    incident_iterator edge_end()const{
      return IncidentIterator(graph_pointer,this,this->degree());
    }


    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((graph_pointer==n.graph_pointer)&&(node_index==n.node_index));
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
      if(graph_pointer==n.graph_pointer){
        return(node_index<n.node_index);
      }
      return(graph_pointer<n.graph_pointer);//compare grapph
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_pointer;// pointer to the graph which the node belongs to
    size_type node_index;//the id of the node

    //initialize the node
    Node(const graph_type* graph_pointer_val,size_type node_index_val)
      :graph_pointer(const_cast<graph_type*>(graph_pointer_val)),node_index(node_index_val){} 

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_index_list.size();
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
  Node add_node(const Point& position, const node_value_type& node_value=node_value_type()) {
    // HW0: YOUR CODE HERE
    //adding node to node registration list
    node_list.emplace_back(node_list.size(),position,node_value);
    //registrate the node to graph
    node_index_list.emplace_back(node_list.size()-1);
    node_num++;
    return Node(this,node_list.size()-1);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE
    if(n.graph_pointer==this){
      //go through all the registered nodes
      for(auto i = node_index_list.begin();i != node_index_list.end();++i){
        if(*i==n.node_index){
          return true;
        }
      }

    } 
    return false;    
    //return ((n.graph_pointer==this)&&(n.node_index<node_list.size()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if   (i<size()){
      return Node(this,i);
    }
    else{
      return Node();//vaild  node
    }
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
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_pointer,graph_pointer->edge_list[edge_id].node1_id);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_pointer,graph_pointer->edge_list[edge_id].node2_id);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((graph_pointer==e.graph_pointer)&&(edge_id==e.edge_id));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if((graph_pointer==e.graph_pointer)){
        return edge_id < e.edge_id;
      }
      return graph_pointer<e.graph_pointer;
    }

    double length() const{
      Point p1 = node1().position();
      Point p2 = node2().position();
      return norm(p1-p2);
    }

    edge_value_type& value(){
      return graph_pointer->edge_list[edge_id].value;
    }

    const edge_value_type& value() const{
      return graph_pointer->edge_list[edge_id].value;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_pointer;
    size_type edge_id;

    Edge(const graph_type* graph_pointer_val,size_type edge_id_val):graph_pointer(const_cast<graph_type*>(graph_pointer_val)),edge_id(edge_id_val){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_num;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this,i);       
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    for(unsigned iter = 0; iter < node_list[a.node_index].neighborhood_list.size();iter++){
      if(node_list[a.node_index].neighborhood_list[iter].neghibor==b.node_index){
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& edge_val=edge_value_type()) {
    // HW0: YOUR CODE HERE
    if(a==b){
      return Edge();
    }
    if(!has_edge(a,b)){
      if(a<b){//adding node
        edge_list.emplace_back(a.node_index,b.node_index,edge_list.size(),edge_val);
        edge_num++;
        //register to neighborhood list
        node_list[a.node_index].neighborhood_list.emplace_back(edge_list.size()-1,b.node_index);
        node_list[b.node_index].neighborhood_list.emplace_back(edge_list.size()-1,a.node_index);
        edge_index_list.emplace_back(edge_list.size()-1);//register the edge
        return Edge(this,edge_list.size()-1);
      }
      else{
        edge_list.emplace_back(b.node_index,a.node_index,edge_list.size(),edge_val);
        edge_num++;
        node_list[a.node_index].neighborhood_list.emplace_back(edge_list.size()-1,b.node_index);
        node_list[b.node_index].neighborhood_list.emplace_back(edge_list.size()-1,a.node_index);
        edge_index_list.emplace_back(edge_list.size()-1);
        return Edge(this,edge_list.size()-1);
      }

    }
    
    return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_num = 0;
    edge_num = 0;
    //clear the node list and the edge list
    //the swap function can clear the capasity at the same time.
    std::vector<NodeProxy>().swap(node_list);
    std::vector<EdgeProxy>().swap(edge_list);
    std::vector<unsigned>().swap(node_index_list);
    std::vector<unsigned>().swap(edge_index_list);
  }

  size_type remove_edge(const Node& n1,const Node& n2){
    if(has_edge(n1,n2)){
      //adj of node1
      for(size_type i = 0;i < n1.degree();++i){
        if(node_list[n1.node_index].neighborhood_list[i].neghibor==n2.node_index){
          size_type edge_id = node_list[n1.node_index].neighborhood_list[i].edge_id;
          //find out node2 and erase
          node_list[n1.node_index].neighborhood_list.erase(node_list[n1.node_index].neighborhood_list.begin()+i);

          //erase node1 from node2's neigh
          for(size_type j = 0;j < n2.degree();++j){
            if(node_list[n2.node_index].neighborhood_list[j].neghibor==n1.node_index){
              node_list[n2.node_index].neighborhood_list.erase(node_list[n2.node_index].neighborhood_list.begin()+j);
            }
          }


          //erase edge from graph
          for(size_type j = 0; j<n1.graph_pointer->edge_index_list.size();++j){
            if(n1.graph_pointer->edge_index_list[j]==edge_id){
              n1.graph_pointer->edge_index_list.erase(n1.graph_pointer->edge_index_list.begin()+j);
            }
          }
          edge_num -= 1;
          return 1;  

        }
      }
    }
    return 0;
  }

  size_type remove_edge(const Edge& e) {
    Node N1 = e.node1();
    Node N2 = e.node2();
    //using remove_edge(const Node&, const Node&) to deleteedge
    size_type num = remove_edge(N1,N2);
    return num;
  }

  edge_iterator remove_edge(edge_iterator edge_it){
    auto e = *edge_it;//get the edge and delete it
    remove_edge(e);
    edge_iterator edge_it_=edge_it;
    return edge_it_;
  }

  size_type remove_node(const node_type& n){
    if(has_node(n)){
      //erase edge in neigh
      while(n.degree()>0){
        auto e_it = n.edge_begin();
        auto e = *e_it;
        remove_edge(e);
      }

      //erase node
      // we have a node list and node index list
      // we don't change node list --> list of node proxy
      // we only need to change the list that registered to the graph
      // we only delete the index from node index list
      for(size_type j = 0; j < n.graph_pointer->node_index_list.size();++j){
        if(n.graph_pointer->node_index_list[j]==n.node_index){
          n.graph_pointer->node_index_list.erase(n.graph_pointer->node_index_list.begin()+j);
        }
      }

      node_num -= 1;
      return 1;
    }
    return 0;
  }

  node_iterator remove_node(node_iterator n_it){
    node_type n = *n_it;
    remove_node(n);
    node_iterator n_it_ = n_it;
    return n_it;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private equality_comparable<NodeIterator> {
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



    Node operator*()const{
      return Node(graph_pointer,graph_pointer->node_index_list[node_index]);
    }

    NodeIterator& operator++(){
      node_index++;
      return *this;
    }

    bool operator==(const NodeIterator& a)const{
      return ((graph_pointer==a.graph_pointer)&&(node_index==a.node_index));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_pointer;
    size_type node_index;

    NodeIterator(const graph_type* graph_pointer_val,size_type node_index_val):graph_pointer(const_cast<graph_type*>(graph_pointer_val)),node_index(node_index_val){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }
  node_iterator node_end() const{
    return NodeIterator(this,this->node_index_list.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private equality_comparable<IncidentIterator>{
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
    Edge operator*() const{
      int node_index = node_pointer->node_index;
      int iter_index = iter_id;
      
      //std::cout<<"operator*"<<std::endl;
      //std::cout<<graph_pointer->node_list[node_index].id_np<<std::endl;
      std::cout<<iter_id<<","<<graph_pointer->node_list[node_index].neighborhood_list.size()<<std::endl;
      //node id -> node proxy
      //node proxy->neighborhood list->>edge id
      return Edge(graph_pointer,graph_pointer->node_list[node_index].neighborhood_list[iter_index].edge_id);
    }
    IncidentIterator& operator++(){
      iter_id++;
      return *this;
    }
    bool operator==(const IncidentIterator& a)const{
      return ((graph_pointer==a.graph_pointer)&&(node_pointer==a.node_pointer)&&(iter_id==a.iter_id));
    }

   private:
    friend class Graph;
    friend class Node;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_pointer;
    node_type* node_pointer;
    size_type iter_id;

    IncidentIterator(const graph_type* graph_pointer_val,const node_type* node_pointer_val,size_type iter_id_val):graph_pointer(const_cast<graph_type*>(graph_pointer_val)),node_pointer(const_cast<node_type*>(node_pointer_val)),iter_id(iter_id_val){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private equality_comparable<EdgeIterator> {
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
    Edge  operator*() const{
      return Edge(graph_pointer,graph_pointer->edge_index_list[edge_id]);
    }

    EdgeIterator& operator++(){
      edge_id++;
      return *this;
    }

    bool operator==(const EdgeIterator& a)const{
      return ((graph_pointer==a.graph_pointer)&&(edge_id==a.edge_id));
    }




   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_pointer;
    size_type edge_id; 
    EdgeIterator(const graph_type* graph_pointer_val,size_type edge_id_val):graph_pointer(const_cast<graph_type*>(graph_pointer_val)),edge_id(edge_id_val){}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin()const {
    return EdgeIterator(this,0);
  }

  edge_iterator edge_end()const {
    return EdgeIterator(this,this->edge_index_list.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.


  //definition of struct NodeProxy
  struct NodeProxy{
    size_type id_np;//node proxy id
    Point P;
    node_value_type value;

    //neigihborhood list
    std::vector<NeghProxy> neighborhood_list{};

    NodeProxy(size_type id_val, Point P_val,node_value_type value_val):id_np(id_val),P(P_val),value(value_val){}


  };

  /////definition of struuuuct EdgeProxy
  struct EdgeProxy{
    size_type node1_id;
    size_type node2_id;
    size_type edge_id;
    edge_value_type value;

    EdgeProxy(size_type node1_id_val, size_type node2_id_val, size_type edge_id_val,edge_value_type value_val):node1_id(node1_id_val),node2_id(node2_id_val),edge_id(edge_id_val),value(value_val){}
  };


  //proxy design for storing neighborhood
  struct NeghProxy{
    size_type neghibor;
    size_type edge_id;

    NeghProxy(size_type edge_id_val, size_type neghibor_val):neghibor(neghibor_val),edge_id(edge_id_val){}

  };

};

#endif // CME212_GRAPH_HPP
