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

template <typename V, typename E>
//template <typename V>
class Graph {
private:
  std::vector<Point> points_;
  //int edge_idx_; //move it down
 
 public:
  // PUBLIC TYPE DEFINITIONS

  /** Type of this graph. */
  using graph_type = Graph;

  using node_value_type = V;

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
    size_ = 0;
    num_edges_ = 0;
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
    Node() { ///Invalid constructor because a user should only construct nodes through the Graph class.
    graph_ = nullptr;
    index_ = 0;

    }

    /** Return this node's position. */
    /// Because the position info of a node is stored in the Points and Nodes vector,
    /// we can return the position from either vector. But because the given function returns Point, we'll use the point vector.
    const Point& position() const {
      return graph_->points_[index_]; 
    }

    /* HW 2: make the node position modifiable*/
    Point& position() {
      return graph_->points_[index_]; 
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      /// Because a node is uniquely defined by its index and this index is up to date. Just return index_.
      //return index_;

      std::vector<size_type>::iterator it = std::find(graph_->i2u.begin(), graph_->i2u.end(), index_);
      if (it != graph_->i2u.end()) {
        return std::distance(graph_->i2u.begin(), it);
      }
    }

    node_value_type& value() {
      //std::cout << "node_value: template size: " << graph_->vals_from_template.size() << std::endl;
      //std::cout << "node_value: node index: " << index_ << std::endl;
      return graph_->vals_from_template[index_];
    }

    const node_value_type& value() const {
     return graph_->vals_from_template[index_];
    }

    /* return the number of incident edges */
    size_type degree() const{
      //because nodes to edges maps the nodes to the nodes, the size of the node vector has to be the degree
      //std::cout << "my degree is " << graph_->nodes_to_nodes[index_].size() << std::endl;
      // std::cout << "nodes: " << graph_->nodes_to_nodes[index_].size() << std::endl; 
      // std::cout << "edges: " << graph_->nodes_to_edges[index_].size() << std::endl; 
      return graph_->nodes_to_edges[index_].size(); 
      
      //return graph_->nodes_to_edges[this->index()].size(); 
    }

    /* Point to the start of the incident iterator */
    IncidentIterator edge_begin() const{
      return IncidentIterator(graph_, index_ , 0);
    }

    /* Point to the end of the incident iterator */
    IncidentIterator edge_end() const
    {
      return IncidentIterator(graph_, index_, this->degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      /// Because the instruction says to check graph and index, we'll only check these. 
      /// One could also check for graph, index, and point
      if (n.graph_ == this->graph_ && n.index_ == this->index_) {
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
      /// Currently, we can only check for the index of the nodes.
      /// Note this only evaluates the order the nodes were added which is kind of a mute point for < to do.
      
      if (this->index_ < n.index_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    bool valid() const {
      return index_ >= 0 && index_ < graph_->nodes_.size()
              && graph_->nodes_[index_].index_ < graph_->i2u_.size()
              && graph_->i2u[graph_->nodes_[index_].index_] == index_;
    }
    //A proper node constructor will have a pointer to the graph it belongs to, and an index. 
    //Note: it won't contain any points to keep it lightweight
    graph_type* graph_;
    size_type index_;
    

    Node(const graph_type* graph, size_type index)
      : graph_(const_cast<graph_type*> (graph)), index_(index) {
      }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    /// Because the graph class has a size attribute, just return that here.
    return size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] a set of values passed in as template
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& given_value = node_value_type()) {
    /// Because our node is tracked by the point and node vectors, 
    /// we will create a node to add to the node vector at the end of the existing Node array and add the given point to the point vector.
    /// Because we are using standard vectors, we can use the standard push_back function

    Node new_node(this, size_); //the new node should be at index size_ for 0 index
    nodes_.push_back(new_node);
    //HW2 addition
    i2u.push_back(size_);
    points_.push_back(position);
    nodes_to_nodes.insert(std::pair<int, std::vector<int> > (size_, {})); //to keep track of all the edges of a node
    nodes_to_edges.insert(std::pair<int, std::vector<int> > (size_, {})); //to keep track of all the edges of a node
    size_ += 1; 
    vals_from_template.push_back(given_value);

    return new_node;      
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    /// Note: all nodes have a pointer to their respective graph.
    /// So, if the Node's graph is the same as this graph we are in, then the Node exists.
    /// This is a much better approach [O(1)] than iterating through each node in a graph and checking for equality

    if (n.graph_ == this) {
      return true;
    }
    return false;
  }

  /** Return the node with index @a i of the active nodes.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    /// Because, there is a Node vector of all nodes, just return the ith element of that vector
    // assert(i < size());
    // return nodes_[i];       

    /// HW2 modificaiton
    /// Because i2u tracks the current active nodes, get the unique id of the desired node from 
    /// i2u, then get the node itself from nodes_
    size_type uid = i2u[i];
    return nodes_[uid];
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
    Edge() { //This is the invalid constructor.
    graph_ = nullptr;
    idx1_ = 0;
    idx2_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {

      //return graph_->nodes_[idx1_]; 
      //HW2 modification
      return graph_->node(idx1_);     
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //return graph_->nodes_[idx2_]; 
      //HW2 modificaiton
      return graph_->node(idx2_);     
    }

    /* Return the Length of the edge from template  */

    const edge_value_type& edge_value() const {
      return graph_->vals_from_template_edge[my_edge_idx_];
    }

    edge_value_type& edge_value() {
      return graph_->vals_from_template_edge[my_edge_idx_];
      //return graph_->vals_from_template_edge[my_edge_idx_];
    }

    double length() const {
      return norm(node1().position() - node2().position());
    }

    size_type edge_index() const {
      return my_edge_idx_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      /// This is checking the reference equality instead of the value equality 
      /// which would have checked the equality of the nodes at the end of the edges
      return ((this->idx1_ == e.idx1_ || this->idx1_ == e.idx2_) && 
        (this->idx2_ == e.idx2_ || this->idx2_ == e.idx1_) &&
        this->graph_ == e.graph_);   
   }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (e.graph_ == this->graph_) {
        size_type node1 = e.idx1_;
        size_type node2 = e.idx2_;

        size_type large_node_e = std::max(node1, node2);
        size_type large_node_this = std::max(this->idx1_, this->idx2_);

        // std::cout << "e node idx: " << e.idx1_ << " and " << e.idx2_ << std::endl;
        // std::cout << "large node e: " << large_node_e << std::endl;
        // std::cout << "e this idx: " << this->idx1_ << " and " << this->idx2_ << std::endl;
        // std::cout << "large node this: " << large_node_this << std::endl;
        // std::cout << " " << std::endl;

        if (large_node_this < large_node_e) {
          return true;
        } else if (large_node_this == large_node_e and std::min(this->idx1_, this->idx2_) < std::min(node1, node2)) {
          return true;
        } else {
          return false;
        }
      } else if (this->graph_->size() <= e.graph_->size()){
        return true; //if the two edges aren't in the same graph, return false
      } else {
        return false;
      }
      }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Constructor
    // An edge connects two nodes. Hence, it is uniquely defined by the indices of those two nodes
    graph_type* graph_;
    size_type idx1_;
    size_type idx2_;
    size_type my_edge_idx_;
    


    Edge(const graph_type* graph, size_type idx1, size_type idx2, size_type my_edge_idx)
    //Edge(graph_type* graph, size_type idx1, size_type idx2)
      //: graph_(const_cast<graph_type*> (graph)), idx1_(idx1), idx2_(idx2)
      : graph_(const_cast<graph_type*> (graph)), idx1_(idx1), idx2_(idx2), my_edge_idx_(my_edge_idx)
      {

      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    /// Because we are already tracking the number of edges in Graph, return that here
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    /// because all edges are stored in the edges_ vector, return the ith index.
    assert(i < (num_edges_));
    return edges_[i];        
  }

  // Return the edge that exists between two given nodes -- helper function
  Edge edge_bn_nodes(const Node& a, const Node& b) {
    for (size_type i = 0; i < num_edges_; ++i) {
      Edge e = edges_[i];
      Node e1 = e.node1();
      Node e2 = e.node2();

      if (((e1.index_ == a.index_) && (e2.index_ == b.index_)) || 
        ((e1.index_ == b.index_) && (e2.index_ == a.index_))){
        return e;
      }
    } 
    //return false;
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
   /// Because the edge is undirected, we need to check both points.
    //for (size_type i = 0; i < num_edges_; ++i) {
    for (size_type i = 0; i < edges_.size(); ++i) {
      Edge e = edges_[i];
      Node n1 = e.node1();
      Node n2 = e.node2();
      

      if (((n1.index_ == a.index_) && (n2.index_ == b.index_)) || 
        ((n1.index_ == b.index_) && (n2.index_ == a.index_))){
        //std::cout << "index in edges_: " << i << std::endl;
        //std::cout << "the two nodes: " << e1.index() << " and " << e2.index() << std::endl;
        //std::cout << "has_edge index " << e.edge_index() << std::endl;
        //std::cout << "has edge between: " << e1.index() << " and " << e2.index() << std::endl;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& given_value_edge = edge_value_type()) {
  //Edge add_edge(const node_value_type& a = node_value_type(), const node_value_type& b = node_value_type()) {
    size_type idx_a = a.index_;
    size_type idx_b = b.index_;

    

    if(has_edge(a, b)){
      for (std::vector<int>::iterator it = nodes_to_edges[a.index()].begin(); it != nodes_to_edges[a.index()].end(); ++it) {
        if (edges_[(*it)].node1() == b or edges_[(*it)].node2() == b) {
          return edges_[(*it)];
        }        
      }
    } else {

      /// Using a map of nodes to nodes instead to track each new edge
      // nodes_to_nodes[idx_a].push_back(b);
      // nodes_to_nodes[idx_b].push_back(a);

      nodes_to_nodes[idx_a].push_back(idx_b);
      //std::cout << idx_b << " added to " << idx_a << std::endl;
      nodes_to_nodes[idx_b].push_back(idx_a);
      //std::cout << idx_a << " added to " << idx_b << std::endl;

      nodes_to_edges[idx_a].push_back(num_edges_);
      nodes_to_edges[idx_b].push_back(num_edges_);

      // std::cout << idx_a << " degree from nodes_to_nodes: " << nodes_to_nodes[idx_a].size() << std::endl;
      // std::cout << idx_a << " degree from nodes_to_edges: " << nodes_to_edges[idx_a].size() << std::endl;
      // std::cout << "........: " << std::endl;

      // std::cout << idx_b << " degree from nodes_to_nodes: " << nodes_to_nodes[idx_b].size() << std::endl;
      // std::cout << idx_b << " degree from nodes_to_edges: " << nodes_to_edges[idx_b].size() << std::endl;

      Edge e(this, idx_a, idx_b, num_edges_);
      edges_.push_back(e);

      num_edges_ += 1;
      vals_from_template_edge.push_back(given_value_edge);
      return e;
    }
         
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */ 
  /// Because all the nodes and edges are stored in a vector, we can just clear the vectors
  void clear() {
    nodes_.clear();
    i2u.clear();
    edges_.clear();
    points_.clear();
    size_ = 0;
    num_edges_ = 0;
    vals_from_template.clear();
    vals_from_template_edge.clear();
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
      graph_ = nullptr;
      idx_ = 0;
    }

    public:
      NodeIterator operator++()
      {
        idx_++;
        return *this;
      }

      bool operator==(const NodeIterator& a_node_itr) const
      {
        return a_node_itr.idx_ == idx_ and a_node_itr.graph_ == this->graph_;
      }

      Node operator*() const
      {
        //return graph_->nodes_[idx_];
        //HW2 modification
        return graph_->node(idx_);
      }  

   private:
   friend class Graph;

    graph_type* graph_;
    size_type idx_;

    //constructor
      NodeIterator(const graph_type* graph, int node_idx) 
      : graph_(const_cast<graph_type*> (graph)), idx_(node_idx)
      { 
      } 
  };

   NodeIterator node_begin() const
   {
    return NodeIterator(this, 0);
   }

   NodeIterator node_end() const
   {
    //int node_length = this->nodes_.size();
    //HW2 Modification
    int node_length = this->i2u.size();
    //this returns one past the last element hence why we don't subtract 1 from node_length
    return NodeIterator(this, node_length); 
   }

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

    /** Construct an invalid IncidentIterator. Note: All of these are O(1)*/
    IncidentIterator() {
      graph_ = nullptr;
      node_idx_ = 0;      
    }

    IncidentIterator& operator++(){
      edge_degree_idx_++;
      return *this;
    }

    Edge operator*() { 
      size_type current_node_idx = graph_->nodes_to_nodes[node_idx_][edge_degree_idx_];
      Node current_node = Node(graph_, current_node_idx);
      Edge current_edge = graph_->edge_bn_nodes(Node(graph_, node_idx_), current_node);
      return Edge(graph_, node_idx_, current_node.index(), current_edge.edge_index());
    }

    bool operator==(const IncidentIterator& inc_itr) const {
      return (inc_itr.graph_ == this->graph_ and inc_itr.node_idx_ == node_idx_ and inc_itr.edge_degree_idx_ == edge_degree_idx_);
    }

   private:
    friend class Graph;
    friend Edge Graph::add_edge();
   
    graph_type* graph_;
    size_type node_idx_;
    size_type edge_degree_idx_;

    // Valid Constructor
    IncidentIterator(graph_type* graph, int node_idx, int edge_degree_idx) 
    : graph_(graph), node_idx_(node_idx), edge_degree_idx_(edge_degree_idx)
    {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  //class EdgeIterator {
  class EdgeIterator: private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      graph_ = nullptr;
      edge_idx_ = 0;
    }

    EdgeIterator operator++(){
      edge_idx_++;
      return *this;
    }

    Edge operator*() const{
      return graph_->edges_[edge_idx_];
    }

    bool operator == (EdgeIterator edge_iter) const {
      return (edge_iter.graph_ == this->graph_ and edge_iter.edge_idx_ == edge_idx_);
    }

   private:
    friend class Graph;

    graph_type* graph_;
    size_type edge_idx_;

    EdgeIterator(const graph_type* graph, int edge_idx)
    : graph_(const_cast<graph_type*> (graph)), edge_idx_(edge_idx) {} 
  };

  EdgeIterator edge_begin(){
    return EdgeIterator(this, 0);
  }

  EdgeIterator edge_end(){
    int edge_length = this->edges_.size();
    return EdgeIterator(this, edge_length);
  }


  /* A function that removes a given Edge.
  * param [in] a valid edge
  * param [out] 1 if an edge was removed and 0 else
  *
  * post: the previous indexes and iterators will be invalidated
  *       the edge will be removed from the edge_ container
  *       the edge will be removed from the nodes_to_edges container
  *       the respective nodes will be removed from the respective nodes_to_nodes containers
  *       the number of edges will be reduced by 1
  *
  */
  size_type remove_edge(const Edge& e) {

    size_type starting_edges_size = edges_.size();

    
    if (starting_edges_size > 0) {
      //std::cout << "remove edge started " << std::endl;
      std::cout << edges_.size() << std::endl;

      //std::cout << " " << std::endl;

      // std::cout << "remove_edge entered" << std::endl;
      // std::cout << "starting edge_ size: " << edges_.size() << std::endl;
      /* Swap and pop */

      size_type e_idx = e.edge_index();

      //Part 2: modify the nodes_to_edges map
      Node n1 = e.node1();
      Node n2 = e.node2();

      // std::cout << "pre degree edge: " << n1.position() << " is " << nodes_to_edges[n1.index()].size() << std::endl;
      // std::cout << "pre degree node: " << n1.position() << " is " << nodes_to_nodes[n1.index()].size() << std::endl;
      //std::cout << "pre degree: " << n2.index() << " is " << n2.degree() << std::endl;

      //std::cout << "edge exists? " << has_edge(n2, n1) << std::endl;
      
      //std::cout << "node:degree " << n1.index() << ":" << n1.degree() << " and " << n2.index() << ":" << n2.degree()<< std::endl;
      //retrive the index location of the current edge in n1 and n2
      //if ((nodes_to_edges[n1.index()].size() > 0) and (nodes_to_edges[n2.index()].size() > 0)){

      // if ((nodes_to_edges[n1.index()].size()) == 1) {
      //   nodes_to_edges[n1.index()].clear();
      // }

      // if ((nodes_to_edges[n2.index()].size()) == 1) {
      //   nodes_to_edges[n2.index()].clear();
      // }

      // if ((nodes_to_nodes[n1.index()].size()) == 1) {
      //   nodes_to_nodes[n1.index()].clear();
      // }

      // if ((nodes_to_nodes[n2.index()].size()) == 1) {
      //   nodes_to_nodes[n2.index()].clear();
      // }

      if ((nodes_to_edges[n1.index()].size() > 0) and (nodes_to_edges[n2.index()].size() > 0)){
        std::cout << "> entered first if " << std::endl;

        std::vector<int>::iterator itr_n1 = std::find(nodes_to_edges[n1.index()].begin(), nodes_to_edges[n1.index()].end(), e_idx);
        std::vector<int>::iterator itr_n2 = std::find(nodes_to_edges[n2.index()].begin(), nodes_to_edges[n2.index()].end(), e_idx);
        std::cout << "size[n1] " << nodes_to_edges[n1.index()].size() << std::endl;
        //std::cout << "size[n2] " << nodes_to_edges[n2.index()].size() << std::endl;

        size_type b1 = (itr_n2 != (nodes_to_edges[n2.index()].end()));
        bool b2 =  (itr_n1 != (nodes_to_edges[n1.index()].end()));

        std::cout << "first bool " << b2 << std::endl;
        std::cout << "second bool " << b1 << std::endl;

        if (itr_n1 != (nodes_to_edges[n1.index()].end()) and itr_n2 != (nodes_to_edges[n2.index()].end())) {
          std::cout << "> entered second if " << std::endl;
          size_type edge_in_n1 = std::distance(nodes_to_edges[n1.index()].begin(), itr_n1);
          size_type edge_in_n2 = std::distance(nodes_to_edges[n2.index()].begin(), itr_n2);
          
          // std::cout << "nodes_to_edges["<< n2.index() << "]pre size: " << nodes_to_edges[n2.index()].size() << std::endl;
          // std::cout << "nodes_to_edges[" << n1.index() << "]pre size: " << nodes_to_edges[n1.index()].size() << std::endl;
          //overwrite and pop from the nodes_to_edges map
          nodes_to_edges[n1.index()][edge_in_n1] = nodes_to_edges[n1.index()].back();
          nodes_to_edges[n1.index()].pop_back();

          //std::cout << "size[n2] pre " << nodes_to_edges[n2.index()].size() << std::endl;

          nodes_to_edges[n2.index()][edge_in_n2] = nodes_to_edges[n2.index()].back();
          nodes_to_edges[n2.index()].pop_back();

          //std::cout << "size[n2] post " << nodes_to_edges[n2.index()].size() << std::endl;  
        }
      }
      //std::cout << "nodes_to_edges["<< n2.index() << "]pre size: " << nodes_to_edges[n2.index()].size() << std::endl;
      //std::cout << "nodes_to_edges["<< n1.index() << "]pre size: " << nodes_to_edges[n1.index()].size() << std::endl;

      //Part 3: Modify the nodes to nodes map
      // std::cout << "nodes_to_nodes["<< n2.index() << "]pre size: " << nodes_to_nodes[n2.index()].size() << std::endl;
      // std::cout << "nodes_to_nodes["<< n1.index() << "]pre size: " << nodes_to_nodes[n1.index()].size() << std::endl;

      if (nodes_to_nodes[n1.index()].size() > 0 and nodes_to_nodes[n2.index()].size() > 0) {
        std::cout << "> entered third if " << std::endl;
        std::vector<int>::iterator node_itr_n1 = std::find(nodes_to_nodes[n1.index()].begin(), nodes_to_nodes[n1.index()].end(), n2.index());
        std::vector<int>::iterator node_itr_n2 = std::find(nodes_to_nodes[n2.index()].begin(), nodes_to_nodes[n2.index()].end(), n1.index());

        if (node_itr_n1 != nodes_to_nodes[n1.index()].end() or node_itr_n2 != nodes_to_nodes[n2.index()].end()) {
          std::cout << "> entered fourth if " << std::endl;
          size_type n2_in_n1 = std::distance(nodes_to_nodes[n1.index()].begin(), node_itr_n1);
          size_type n1_in_n2 = std::distance(nodes_to_nodes[n2.index()].begin(), node_itr_n2);

          //Overwrite and pop from the nodes_to_nodes map
          nodes_to_nodes[n1.index()][n2_in_n1] = nodes_to_nodes[n1.index()].back();
          nodes_to_nodes[n1.index()].pop_back();

          nodes_to_nodes[n2.index()][n1_in_n2] = nodes_to_nodes[n2.index()].back();
          nodes_to_nodes[n2.index()].pop_back();
        }
      }

      // std::cout << "post degree edge: " << n1.position() << " is " << nodes_to_edges[n1.index()].size() << std::endl;
      // std::cout << "post degree node: " << n1.position() << " is " << nodes_to_nodes[n1.index()].size() << std::endl;
      //std::cout << "post degree: " << n2.index() << " is " << n2.degree() << std::endl;
            // Part I: modify edges_
      //Find the index of the given nedge

      
      //std::cout << "current edge index " << e_idx <<std::endl;
      //replace the current edge with the last edge in the vector
      //edges_[e_idx] = (*(edges_.end() - 1));
      //std::cout << "edges_ pre size: " << edges_.size() << std::endl;
      edges_[e_idx] = edges_.back();
      //pop out the last edge (which is now the given edge)
      edges_.pop_back();
      //std::cout << "edges_ post size: " << edges_.size() << std::endl; 
      
      num_edges_ -= 1;
      return 1;

      } else {
        //std::cout << "size is zero " << edges_.size() << std::endl;
        return 0;
      }

  }
  
  /* Removes an edge given two nodes
  * param[in] two nodes
  *
  * post: the above remove_edge will be called and all the necessary containers will be modified
          the index, pointer, and other properties of the given index will change
  */
  size_type remove_edge(const Node&a, const Node& b){
    Edge e = edge_bn_nodes(a, b);
    return remove_edge(e);
  }

  /* Removes an edge given an edge iterator
  * param[in] two nodes
  *
  * post: the above remove_edge will be called and all the necessary containers will be modified
          the index, pointer, and other properties of the given index will change
  */
  edge_iterator remove_edge (edge_iterator e_it){
    Edge e = (*e_it);
    remove_edge(e);
    return e_it;
  }

  /* A function that removes a given Node.
  * param [in] a valid node
  * param [out] 1 if a node was removed and 0 else
  *
  * post: the previous indexes and iterators will be invalidated
  *       the node will be removed from the i2u container
  *       the node will be removed from the nodes_to_edges container
  *       the incident nodes will be removed from the respective nodes_to_nodes containers
  *       the number of nodes will be reduced by 1
  *
  */
  size_type remove_node(const Node& n) {
    //start the removal process only if there are active nodes remaining
    std::cout << " " <<std::endl;
    //std::cout << "remove node entered  " << std::endl;
    // std::cout << "graph size pre: " << this->num_nodes() << std::endl;
    // std::cout << "i2u size pre: " << i2u.size() << std::endl;
    if (i2u.size() > 0){
      //First, remove all the incident edges

      //std::cout << "index: degree pre = " << n.index() << " : " << n.degree() <<std::endl;
      
      //size_type result = 1;
      // while(nodes_to_edges[n.index()].size() > 0) {
      //   //std::cout << "degree in while " << n.degree() << std::endl;
      //   //std::cout << "entered while loop " <<std::endl;
      //   auto it = n.edge_begin();
      //   //std::cout << "pre remove edge " <<std::endl;
      //   remove_edge(*it);

      //   //n.degree() = n.degree() - 1;
      // }

      // for(auto idx = 0; idx < nodes_to_edges[n.index()].size(); ++idx){
      //   auto it = n.edge_begin();
      //   size_type result = remove_edge(*it);
      // }

      auto it = n.edge_begin();
      while(it != n.edge_end()){
        remove_edge(*it);
      }

      //std::cout << "index: degree pre = " << n.index() << " : " << n.degree() <<std::endl;

      std::vector<size_type>::iterator itr = std::find(i2u.begin(), i2u.end(), n.index());

      if (itr != i2u.end()) {
        size_type n_i2u_idx = std::distance(i2u.begin(), itr);

        std::map<int, std::vector<int>>::iterator it1 = nodes_to_nodes.find(n_i2u_idx);
        std::map<int, std::vector<int>>::iterator it2 = prev(nodes_to_nodes.end()); //prev(a.end()) find the iterator to the last element
        //swap
        if ((it1 != nodes_to_nodes.end()) && (it2 != nodes_to_nodes.end())) {
          std::swap(it1->second, it2->second);
        }
        //then erase. because this erase is one element, it is constant time
        nodes_to_nodes.erase(prev(nodes_to_nodes.end()));

        //Second, update nodes_to_edges for n
        std::map<int, std::vector<int>>::iterator it3 = nodes_to_edges.find(n_i2u_idx);
        std::map<int, std::vector<int>>::iterator it4 = prev(nodes_to_edges.end()); //prev(a.end()) find the iterator to the last element
        //swap
        if ((it3 != nodes_to_edges.end()) && (it4 != nodes_to_edges.end())) {
          std::swap(it3->second, it4->second);
        }
        //then erase. because this erase is one element, it is constant time
        nodes_to_edges.erase(prev(nodes_to_edges.end()));


        i2u[n_i2u_idx] = i2u.back();
        i2u.pop_back();

      }

      //std::cout << "got to while loop " <<std::endl;


      //std::cout << "edge removal completed." << std::endl;
      size_ -= 1;

      // Second remove the node itself



      //Unlike edges, this n.index() is the value we want in i2u vector which 
      // may not be it's index location. So, let's find the index location

      //std::vector<size_type>::iterator itr = std::find(i2u.begin(), i2u.end(), n.index());

      
        


        // std::cout << "size of i2u: " << i2u.size() << std::endl;
        // std::cout << "node index: " << n.index() << " is located at " << n_i2u_idx << ". Retrived:" << i2u[n_i2u_idx] << std::endl;
        //swap and pop
        //std::cout << "i2u[n_i2u_idx] pre pop " << i2u[n_i2u_idx] << " and back " << i2u.back() << std::endl;
 
        //std::cout << "i2u[n_i2u_idx] post pop " << i2u[n_i2u_idx] << std::endl;

        // Third, adjust graph size
        

        return 1;
      }
      //std::cout << "graph size post: " << this->num_nodes() << std::endl;
      return 0;
    }


  size_type remove_node(node_iterator n_it){
    Node n = (*n_it);
    return remove_node(n);
  }



 private:

  // NODE: represented by points and indices
  // vector of points
  

  // vector of Nodes
  // The comment in my HW0 grading suggested that I wouldn't need this nodes_.
  // I agree that it is not needed since a node is uniquely defined by a point and an int. 
  // However, I have made the concious decision to keep it because it has already been woven into my code
  // and uravelling the current structure would be counter productive.
  
  // HW2:
  // the set of all nodes that have ever been added to the graph
  std::vector<Node> nodes_;

  // a map of currently active nodes and their unique id
  //std::map<size_type, size_type> idx_to_uid;
  std::vector<size_type> i2u;

  //EDGE
  //Vector of edges
  std::vector<Edge> edges_;

  //OTHER
  //size of the graph
  size_type size_; 
  //number of edges in a graph
  size_type num_edges_; 
  //a map containing index_of_node:vector of index of edges of that node
  std::map<int, std::vector<int> > nodes_to_nodes; 
  
  std::map<int, std::vector<int> > nodes_to_edges; 

  // HW2:
  // the set of all nodes that have ever been added to the graph
  std::vector<node_value_type> vals_from_template;

  std::vector<edge_value_type> vals_from_template_edge;
};

#endif // CME212_GRAPH_HPP