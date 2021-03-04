#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
//--functionality_0
//--mising header files, e.g. <map>
//--START
#include <map>
//--END

//--functionality_3
//--too many compiler errors, did you tag the wrong commit?
//--END

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
/** In editing my HW1 to be functional for HW2, I got
 *
 * Help from the following files in the PeerCode folder:
 *   Graph_1611.hpp
 *   Graph_1308.hpp
 *   Graph_1275.hpp
 * I found that all three of these files had implemented
 * the Graph class with similar structures to my own,
 * and thus they were very helpful in helping me understand
 * where I went wrong!
 *
 */

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
//typedef std::vector<Node>::const_iterator node_iterator_;
class Graph {
private:

public:
  using node_value_type = V;
  using edge_value_type = E;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //--functionality_1
  //--size_type not defined
  //--START
  std::vector<size_type> public_nodes; // this is the public vector of non-torn nodes.
  std::vector<size_type> public_edges; // this is the public vector of non-torn edges.
  //--END
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

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

  /** Construct an empty graph. */
  Graph() {
    this->nodes = std::vector<Node>(0);
    this->edges = std::vector<Edge>(0);
    this->values_vec = std::vector<node_value_type>(0);
    this->adj_list = std::vector<std::vector<size_type>>(0);
    size_type z = 0;
    this->node_size = z;
    this->edge_size = z;
    this->public_nodes= std::vector<size_type>(0);
    this->public_edges= std::vector<size_type>(0);
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
      //nothing here! define a real one later in private
    }

    /** Return this node's position. */
    const Point& position() const {
      return pos;
    }

    /** HW 2, now need to make node position modifiable*/
    Point& position(){
        return pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx; //idx is a pointer so we need to dereference it
    }

    /** Return reference to node's value */
    node_value_type &value() {
        return g_->values_vec[idx];
    };

    /** Return constant reference to node's value */
    const node_value_type &value() const {
        return g_->values_vec[idx];
    };

    std::vector<Node> neighbours(){
        return g_->adj_list[idx];
    }
    size_type degree() const {
        return g_->adj_list[idx].size();
        //the degree of this node is equal to the length of its entry in the adj list
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {
      //--style_1
      //--delete placeholder code
      //--START
      (void) n;          // Quiet compiler warning
      //--END
      if (idx == n.index()) {
        return true;
      } else { return false; }
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
      (void) n;           // Quiet compiler warning
      if (idx < n.index()) {
        return true;
      } else { return false; }
    }
    /**
     *
     * std::vector<Node> nbrs = n.neighbours();
      size_type i = n.index();
      std::vector<Edge> adj_edge = std::vector<Edge>(0);
      for (auto it = nbrs.begin(); it != nbrs.end(); ++it){
          size_type j = *it.index();
          std::tuple<size_type, size_type> pair1 = {i, j};
          std::tuple<size_type, size_type> pair2 = {j, i};
          Edge e = Edge();
          if (edge_map.find(pair1) == edge_map.end()){
              e = edge_map[pair2];
          }else{
              e = edge_map[pair1];
          }
     */
    Edge find_edge(Node b){
        std::tuple<size_type, size_type> pair1 = {idx, b.index()};
        std::tuple<size_type, size_type> pair2 = {b.index(), idx};
        if (*g_->edge_map.find(pair1) == *g_->edge_map.end()){
            return *g_->edge_map[pair2];
        }else{
            return *g_->edge_map[pair1];
        }
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph* g_; // the graph to which the node belongs
    size_type idx; // the node's index
    Point pos; // a pointer to the node's position/Point. //fix this!!!!
    //node_value_type val; // the value!
    Node(Graph* graph, size_type index, Point position) { //removed & after size_type
      this->idx = index;
      this->g_ = graph;
      this->pos = position; // reference to the address of the Point object
      //this->val = 0;
    }
  };


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_size;
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
  //Node add_node(const Point& position, const node_value_type& = node_value_type());
  Node add_node(const Point &position, const node_value_type& v = node_value_type()) {
    adj_list.push_back(std::vector<size_type>(0)); //add an empty vector to the node_size index of the adjacency list
    const size_type new_ind = node_size;
    const Node new_node = Node(this, new_ind, position); //added const here
    nodes.push_back(new_node); //add a new node to the nodes vector specifying index and position
    values_vec.push_back(v); //add new value to values_vec
    int int_ind = static_cast<int>(node_size);
    //--functionality_0
    //--typo: pushback -> push_back
    //--START
    public_nodes.push_back(int_ind); // this is a NEW vector (for hw2) that keeps index values for "non-torn" nodes.
    //--END
    this->node_size += 1;
    return nodes[int_ind]; // this is the node most recently added to the nodes vector.
  }
  size_type remove_edge(const Node& n1, const Node& n2){
      size_type i = n1.index();
      size_type j = n2.index();
      std::tuple<size_type, size_type> pair1 = {i, j};
      std::tuple<size_type, size_type> pair2 = {j, i};
      Edge e = Edge();
      if (edge_map.find(pair1) == edge_map.end()){
          e = edge_map[pair2];
      }else{
          e = edge_map[pair1];
      }
      auto it = std::find(public_edges.begin(), public_edges.end(), e.index()); // find e's id in the public_nodes vector
      if (it != public_edges.end()){
          std::iter_swap(it , public_edges.end());
          public_edges.pop_back();
          return 1;
      }
      return 0;
  }
  size_type remove_edge(const Edge& e){
      auto it = std::find(public_edges.begin(), public_edges.end(), e.index()); // find e's id in the public_nodes vector
      if (it != public_edges.end()){
          std::iter_swap(it , public_edges.end());
          public_edges.pop_back();
          return 1;
      }
      return 0;
  }
  edge_iterator remove_edge(edge_iterator e_it){
      auto it = std::find(public_edges.begin(), public_edges.end(), *e_it->index()); // find e's id in the public_nodes vector
      std::iter_swap(it , public_edges.end());
      public_edges.pop_back();
      return public_edges.begin();
  }
  size_type remove_node(const Node& n){
      std::vector<Edge> adj_es = adjacent_edges(n);
      for (auto adj_it = adj_es.begin(); adj_it != adj_es.end(); ++adj_es){
          remove_edge(*adj_es); //*adj_es is an edge because adj_es is a pointer to an edge.
      }
      auto it = std::find(public_nodes.begin(), public_nodes.end(), n.index()); // find n's id in the public_nodes vector
      if (it != public_nodes.end()){
          std::iter_swap(it , public_nodes.end());
          public_nodes.pop_back();
          return 1;
      }
      return 0;
  }

  node_iterator remove_node(node_iterator n_it){
      std::vector<Edge> adj_es = adjacent_edges(*n_it);
      for (auto adj_it = adj_es.begin(); adj_it != adj_es.end(); ++adj_es){
          remove_edge(*adj_es); //*adj_es is an edge because adj_es is a pointer to an edge.
      }
      auto it = std::find(public_nodes.begin(), public_nodes.end(), *n_it->index()); // find n's id in the public_nodes vector
      std::iter_swap(it , public_nodes.end());
      public_nodes.pop_back();
      return public_nodes.begin();
  }
  /**
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    adj_list.push_back(std::vector<size_type>(0)); //add an empty vector to the node_size index of the adjacency list
    const size_type new_ind = node_size;
    const Node new_node = Node(new_ind, position); //added const here
    nodes.push_back(new_node); //add a new node to the nodes vector specifying index and position
    //this -> val = v;
    int int_ind = static_cast<int>(node_size);
    this -> node_size+=1;
    return nodes[int_ind]; // this is the node most recently added to the nodes vector.
  }
   */

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    if (n.index() > node_size - 1) {
      return false;
    } else {
      return true;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    int int_ind = static_cast<int>(i);
    return nodes[int_ind];
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
    Edge() { //nothing here, we will do the real constructor later
    }

    /** Return a node of this Edge */
    Node node1() const {
      return left;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return right;
    }

    /** Return the index of this Edge */
    size_type idx() const {
      return index;
    }

    edge_value_type& value() {
        return g_->edge_values[index];
    }
    const edge_value_type& value() const{
        return g_->edge_values[index];
    }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      if (index == e.idx()) {
        return true;
      } else { return false; }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {
      if (index < e.idx()) {
        return true;
      } else { return false; }
    }

    /** HW2, add a length function to Edge */
    double length() const{
        return g_->edge_values[index].v; //edge_values is a vector
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Node left; // pointer to one end of edge
    Node right; //pointer to second end of edge
    size_type index;
    //double len; //length of edge
    Graph* g_; // the graph to which the node belongs

    Edge(const size_type idx, const Node a, const Node b, double l = 1.) {
      this->left = a;
      this->right = b;
      this->index = idx;
      //this ->len = l;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_size;
  }
/** might not use*/
  std::vector<Edge> adjacent_edges(Node n){
      std::vector<Node> nbrs = n.neighbours();
      size_type i = n.index();
      std::vector<Edge> adj_edge = std::vector<Edge>(0);
      for (auto it = nbrs.begin(); it != nbrs.end(); ++it){
          size_type j = *it.index();
          std::tuple<size_type, size_type> pair1 = {i, j};
          std::tuple<size_type, size_type> pair2 = {j, i};
          Edge e = Edge();
          if (edge_map.find(pair1) == edge_map.end()){
              e = edge_map[pair2];
          }else{
              e = edge_map[pair1];
          }
          adj_edge.push_back(e);
      }
      return adj_edge;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edges[i]; //ith edge in edges vector
  }
  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {
    int i = a.index();
    int j = b.index();
    if (adj_list[i].empty()) { //this means that i is not connected to anything, but maybe j is!
      return false; // ith list is empty
    } else if (std::find(adj_list[i].begin(), adj_list[i].end(), j) != adj_list[i].end()) {
      return true; //j is in the ith's adj list
    } else {
      return false;
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


  Edge add_edge(const Node& a, const Node& b, const edge_value_type& e = edge_value_type()) {
    size_type i = a.index();
    size_type j = b.index();
    std::tuple<size_type, size_type> pair1 = {i, j};
    std::tuple<size_type, size_type> pair2 = {j, i};
    if(has_edge(a, b)){
        if (edge_map.find(pair1) == edge_map.end()){ // this means the edge was entered with key (j,i)
            return edges[edge_map[pair2]]; // edge_map[pair2] is the value associated with key (i, j).
        } else{ // this means the edge was entered in the map with key (i,j)
            return edges[edge_map[pair1]]; // edge_map[pair1] is the value associated with key (i, j).
        }
    } else{
        adj_list[i].push_back(j); // connect i to j in the adjacency list
        adj_list[j].push_back(i); //connect j to i in the adjacency list
        Edge new_edge = Edge(edge_size, a, b);
        edges.push_back(new_edge); // add the edge in the edges vector
        edge_values.push_back(e); //add the edge value to the edge_values vec.
        int int_ind = static_cast<int>(edge_size);
        public_edges.pushback(int_ind); // this is a NEW vector (for hw2) that keeps index values for "non-torn" edges.
        edge_map[pair1] = edge_size; //pair1 corresponds to value edge_size.
        edge_size += 1; //update edge_size
        return edges[int_ind]; // return the last added edge.
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    this -> nodes = std::vector< Node>(0);
    this -> edges = std::vector< Edge>(0);
    this -> values_vec = std::vector<node_value_type>(0);
    this -> adj_list = std::vector<std::vector<size_type>>(0);
    size_type z = 0;
    this -> node_size = z;
    this -> edge_size = z;
    this -> edge_map = std::map<std::tuple<size_type, size_type>,  size_type>();
    this->public_nodes= std::vector<size_type>(0);
    this->public_edges= std::vector<size_type>(0);
  }

  //
  // Node Iterator
  //
  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator :private totally_ordered<NodeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    //typename std::vector<Node>::const_iterator node_iterator_;

    /** Construct an invalid NodeIterator. */
    NodeIterator(){
    }

    Node operator*() const{
      return *p;
    }
    NodeIterator& operator++(){
      p++;
      return *this;
    }

    bool operator==(const NodeIterator& x)const{
      return p == x.p;
    }
/**
    bool operator !=(const NodeIterator& x) const{
        return p != x.p;
    }*/
  private:
      friend class Graph;
      pointer p;
      NodeIterator(pointer pointer1){
        this -> p = pointer1;
      }
    // HW1 #2: YOUR CODE HERE
  };
/** syntax help on these iterator functions from peercode
 * Graph_1275.hpp which I found similar to my own implementation */

  node_iterator node_begin() const{
      return NodeIterator(const_cast<Node*>(&nodes[0])); //just use the constant iterator for the vector
  }

  node_iterator node_end() const{
      return NodeIterator(const_cast<Node*>(&nodes[nodes.size()])); //just use the constant iterator for the vector
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered <IncidentIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    //pointer p;
    /** Construct an invalid IncidentIterator. */
    IncidentIterator(){

    };

    Edge operator*() const{
      return *p;
    }
    IncidentIterator& operator++(){
      p++;
      return *this;
    }

    bool operator==(const IncidentIterator& x)const{
      return p == x.p;
    }
/**
    bool operator!=(const IncidentIterator& x)const{
      return p != x.p;
    }*/
    // HW1 #3: YOUR CODE HERe

  private:
    friend class Graph;
    pointer p;
  };
  //typename std::vector<Edge>::const_iterator begin() { return edges.begin();}
  //typename std::vector<Edge>::const_iterator end() { return edges.end(); }
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
    //pointer p;
    /** Construct an invalid EdgeIterator. */
    EdgeIterator(){

    };
    Edge operator*() const{
      return *p;
    }
    EdgeIterator& operator++(){
      p++;
      return *this;
    }

    bool operator==(const EdgeIterator& x)const{
      return p == x.p;
    }
    /**
    bool operator!=(const EdgeIterator& x)const{
      return p != x.p;
    }*/
    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

  private:
    friend class Graph;
    pointer p;
    EdgeIterator(pointer pointer1){
        this -> p = pointer1;
    }

  };

  edge_iterator edge_begin(){
      return EdgeIterator(edges.begin());
  }

  edge_iterator edge_end(){
    return EdgeIterator(edges.end());
  }
  //typename std::vector<Edge>::const_iterator edge_begin() { return edges.begin();}
  //typename std::vector<Edge>::const_iterator edge_end() { return edges.end(); }
  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

private:
  std::vector< Node> nodes; // a vector of references to all the nodes in the graph
  std::vector< Edge> edges; // a vector of references to all the edges in the graph
  std::vector<node_value_type> values_vec;
  std::vector<edge_value_type> edge_values;
  size_type node_size; // an size_type defining the size of the graph. not const because it changes
  size_type edge_size; // an size_type defining the number of edges in the graph. not const because it changes
  std::vector<std::vector<size_type>> adj_list;
  std::map<std::tuple<size_type, size_type>, size_type> edge_map;
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
