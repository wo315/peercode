#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


template <typename V, typename E> 

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

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
  Graph() :nodes_(), edges_(){
    // HW0: YOUR CODE HERE
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
  class Node : private totally_ordered<Node>{
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
      graph_ = NULL; // initialze graph to be null
    }

    Point& position() {
          return (graph_->nodes_)[index_].point;
    }
    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // return Point();
      return graph_->nodes_[index_].point; // position of nodes index
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // return size_type(-1);
      return index_; // retun node index
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    node_value_type& value() {

      return graph_->nodes_[index_].value; 
    }

    const node_value_type& value() const {

      return graph_->nodes_[index_].value; 
    }

    size_type degree() const {
      return graph_->nodes_[index_].adj_list.size(); // degree of connectivity of nodes
    }

    incident_iterator edge_begin() const {

      return IncidentIterator(this, 0); // starting edge of a node
    }

    incident_iterator edge_end() const{
      return IncidentIterator(this,degree());
    }




    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE

      // decide if two nodes are equal
      // 1.same graph
      // 2.same index
      return (this->graph_ == n.graph_) && (this->index() == n.index()); 
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

      return this->index()< n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    friend class LoggingPoint;
    Graph* graph_;
    size_type index_; // node index

    // constructor
    Node (const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  // num of nodes in a graph
  size_type size() const {
    // HW0: YOUR CODE HERE
    // return 0;
    return nodes_.size();
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
  // 
  Node add_node(const Point& position, const node_value_type& myvalue= node_value_type()) {
    // HW0: YOUR CODE HERE

    std::vector<size_type> adj_list;
    nodes_.push_back({position,myvalue,adj_list}); 
    return Node(this,size()-1);

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE

    if(n == node(n.index()))
      return true;
    return false;
  }

  size_type remove_node (const Node& n) {

      // Check that n is a valid node
      assert( n.index() < size());
      assert( Node(this, n.index()) == n);

      // Remove edges connected with node, being careful to not use an iterator
      while(!(nodes_[n.index()].adj_list.empty()))
        remove_edge(Edge(this, nodes_[n.index()].adj_list[0]));

      // Swtich the last node with the one we have
      change_node_id(node(size() - 1), n.index());

      // Remove the end node
      nodes_.pop_back();

      // Return the index of the node that was placed in n's spot
      return n.index();
    }

  size_type remove_edge(const Node& a, const Node& b) {

      for (auto it = a.edge_begin(); it != a.edge_end(); ++it) {

        if (it.node2() == b || it.node1() == b) {

          return remove_edge(*it);
        }
      }

      return 0;
    };

    size_type remove_edge(const Edge& e) {

      // Checks that it could be a valid edge
      assert (e.index() < e.graph_->num_edges());

      // Removes edge e from the adjacency lists of each of its nodes
      // Uses the erase-remove idiom.
      nodes_[e.node1().index()].adj_list.erase(std::remove(
                  nodes_[e.node1().index()].adj_list.begin(),
                  nodes_[e.node1().index()].adj_list.end(), e.index()),
                  nodes_[e.node1().index()].adj_list.end());

      nodes_[e.node2().index()].adj_list.erase(std::remove(
                  nodes_[e.node2().index()].adj_list.begin(),
                  nodes_[e.node2().index()].adj_list.end(), e.index()),
                  nodes_[e.node2().index()].adj_list.end());


      // Changes the uid of the current last endge to the uid of the current edge
      change_edge_uid(Edge(this, num_edges()-1), e.index());

      // Remove the last edge now
      edges_.pop_back();

      // Return the current edge index, now with a new edge in it.
      return e.index();
      }


    void change_edge_uid(Edge e, size_type new_edge_uid) {

    // Checks if there is no work to be done
    if (e.index() == new_edge_uid)
      return;

    // Setting variable for code clarity
    size_type old_edge_uid = e.index();

    // set the new uid
    edges_[new_edge_uid] = edges_[old_edge_uid];


    // Replace the uids in the node adjacency lists
    std::replace(nodes_[e.node1().index()].adj_list.begin(),
                 nodes_[e.node1().index()].adj_list.end(),
                 old_edge_uid, new_edge_uid);

    std::replace(nodes_[e.node2().index()].adj_list.begin(),
                 nodes_[e.node2().index()].adj_list.end(),
                 old_edge_uid, new_edge_uid);

    return;
  }
  
    void change_node_id(Node node, size_type new_node_index) {

    // Check preconditions
    assert (nodes_[new_node_index].adj_list.size() == 0);
    assert (node == Node(this, node.index()));

    // Check if node has the right index already
    size_type old_node_index = node.index();
    if (old_node_index == new_node_index)
      return;

    // Change the edges in the adjacency list to be refering to the new node index
    for (auto i = node.edge_begin(); i != node.edge_end(); ++i) {

      if (edges_[i.index()].index_1 == old_node_index)
        edges_[i.index()].index_1 = new_node_index;

      else if (edges_[i.index()].index_2 == old_node_index)
        edges_[i.index()].index_2 = new_node_index;

    }

    // Switch the nodes
    nodes_[new_node_index] = nodes_[node.index()];

  }
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE

    assert(i<size());
    return(Node(this,i)); // retun the ith node
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
  class Edge :private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      graph_ = NULL; // initialize graph to be null
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // return Node();      // Invalid Node
      size_type index; 

      if (nodes_flipped_)//flip node
        index  = graph_->edges_[edge_index_].index_2;//return node 2

      else
        index  = graph_->edges_[edge_index_].index_1;// return node 1

      return Node(graph_,index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // return Node();      // Invalid Node
      size_type index;
      if (nodes_flipped_)
        index = graph_->edges_[edge_index_].index_1;// return node 1

      else
        index = graph_->edges_[edge_index_].index_2;//  retun node 2

      return Node(graph_, index); 

    }
    double length() const {

      return norm(node1().position() - node2().position());
    }
    edge_value_type &value() {
          return (graph_->edges_)[edge_index_].value;
    }

    const edge_value_type &value() const {
          return (graph_->edges_)[edge_index_].value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      //HW0: YOUR CODE HERE


    
      if (this->graph_ == e.graph_)// if the graphs are the same
      {
        //decide if the two nodes are equal
        if (this->node1() == e.node1() && this->node2() == e.node2()) 
        {
          return true;
        }
        else if (this->node2() == e.node1() && this->node1()== e.node2())
        {
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

      //HW0: YOUR CODE HERE

      return (edge_index_ < e.index() || graph_ < e.get_graph_pointer());

    }

    edge_iterator edge_begin() const {
      return EdgeIterator(this, 0);
    }

  // Returns iterator after last edge in this graph
    edge_iterator edge_end() const {
      return EdgeIterator(this, num_edges());
    }
    Graph *get_graph_pointer() const {
          return graph_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects


    Edge(const Graph* graph, size_type index) // constructor
       : edge_index_(index), graph_(const_cast<Graph*>(graph)), nodes_flipped_(false) {
    }

    size_type index() const {

      return edge_index_; // 
    }

    void flip_nodes() {

      nodes_flipped_ = !nodes_flipped_;
    }



    size_type edge_index_; // 
    Graph* graph_;
    bool nodes_flipped_; // decide if nodes are flipped
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    // return 0;
    return edges_.size();// return num of edges
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE

    assert(i < num_edges());
    (void) i;
    return Edge(this,i); // retun ith edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE


    // For all edges connected by node a, if another edge equals to b, then return true.
    for (auto iter = a.edge_begin(); iter !=a.edge_end(); ++iter)
    {
      if(iter.node2() == b)
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE

   const edge_value_type myvalue= edge_value_type();


   assert (a.index() < size() && b.index() < size());
   assert (a.index() != b.index());


   std::vector<size_type> a_adj_list = nodes_[a.index()].adj_list;// All adjust nodes of node a
    for (std::vector<size_type>::iterator it = a_adj_list.begin(); it != a_adj_list.end();  it++) {

      if (Edge(this, *it).node2() == b)// If node b exists then return edge.
        return Edge(this, *it);

      if (Edge(this, *it).node1() == b) { // return flipped edge
        Edge new_edge = Edge(this, *it);
        // flip edge
        new_edge.flip_nodes();
        return new_edge;
      }
    }

        // idx of new edge
    size_type edge_idx = num_edges(); // 

    // add edge to the adjacent list, as a new edge in the edge vector
    nodes_[b.index()].adj_list.push_back(edge_idx);
    nodes_[a.index()].adj_list.push_back(edge_idx);
    edges_.push_back({a.index(), b.index(), myvalue});

    return Edge(this, edge_idx);



  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    Node operator*() const { // return node
        return graph_->node(curr_node_index_);
     }


    node_iterator& operator++() { // add node
        curr_node_index_++;
        if (curr_node_index_ >= graph_->size()) {
          curr_node_index_ = graph_->size();
        }
        return *this;

     }

    bool operator==(const node_iterator& other_iter) const { 
       if (graph_ == other_iter.graph_ && curr_node_index_ == other_iter.curr_node_index_)
        //if index of nodes are equal
          return true;
       return false;
     }


   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    NodeIterator(const Graph* graph, size_type index)
      : graph_(const_cast<Graph*>(graph)) { // constructor, graph and nodeIterator index

      assert( index <= graph_->size() );

      curr_node_index_ = index;
    }
    Graph* graph_;
    size_type  curr_node_index_; // current node index

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

 node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  // Return the iterator pointing after all nodes
node_iterator node_end() const {
  return NodeIterator(this, size());
}
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator> {
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

    Edge operator*() const { // return current edge
        assert ( curr_edge_index_ < graph_->node(curr_node_index_).degree() );

        Edge curr_edge = Edge(graph_,
                              graph_->nodes_[curr_node_index_].adj_list[curr_edge_index_]);
        if (curr_edge.node1().index() != curr_node_index_)
          curr_edge.flip_nodes();
        return curr_edge;
     }

    incident_iterator& operator++() { // 
       curr_edge_index_++;
       if (curr_edge_index_ >= graph_->node(curr_node_index_).degree())
         curr_edge_index_ = graph_->node(curr_node_index_).degree();

       return *this;
     }


     bool operator==(const incident_iterator& other_iter) const { 
       if (graph_ == other_iter.graph_)
          if(curr_node_index_ == other_iter.curr_node_index_)//if the node index are equal
            if (curr_edge_index_ == other_iter.curr_edge_index_)  //if the edge index are equal
              return true;

       return false;
     }


    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // return Node();      // Invalid Node
      size_type index; 

      if (nodes_flipped_)//flip node
        index  = graph_->edges_[curr_edge_index_].index_2;//return node 2

      else
        index  = graph_->edges_[curr_edge_index_].index_1;// return node 1

      return Node(graph_,index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // return Node();      // Invalid Node
      size_type index;
      if (nodes_flipped_)
        index = graph_->edges_[curr_edge_index_].index_1;// return node 1

      else
        index = graph_->edges_[curr_edge_index_].index_2;//  retun node 2

      return Node(graph_, index); 

    }
    size_type index() const {

      return curr_edge_index_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE


    IncidentIterator(const Node* node, size_type index)
      : curr_node_index_(node->index()) {

      graph_ = const_cast<Graph*>(node->graph_);

      assert( index <= node->degree() );

      curr_edge_index_ = index;

    }
    Graph* graph_;
    size_type curr_node_index_; // current node index
    size_type curr_edge_index_; // current edge index
    bool nodes_flipped_;


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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    Edge operator*() const { 

        assert (curr_edge_index_ < graph_->num_edges());
        return Edge(graph_, curr_edge_index_);// return current edge
    }

    edge_iterator& operator++() {
       if (curr_edge_index_ != graph_->num_edges()) {
         curr_edge_index_++;
       } // 
       return *this;
    }

    bool operator==(const edge_iterator& other_iter) const {
       return (graph_ == other_iter.graph_ && 
                curr_edge_index_ == other_iter.curr_edge_index_);// If the current edge index are equal
    }
    Node node1() const {
      // HW0: YOUR CODE HERE
      // return Node();      // Invalid Node
      size_type index; 

      if (nodes_flipped_)//flip node
        index  = graph_->edges_[curr_edge_index_].index_2;//return node 2

      else
        index  = graph_->edges_[curr_edge_index_].index_1;// return node 1

      return Node(graph_,index);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // return Node();      // Invalid Node
      size_type index;
      if (nodes_flipped_)
        index = graph_->edges_[curr_edge_index_].index_1;// return node 1

      else
        index = graph_->edges_[curr_edge_index_].index_2;//  retun node 2

      return Node(graph_, index); 

    }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    EdgeIterator(const Graph* graph, size_type index) // constructor for edgeIterator
      : graph_(const_cast<Graph*>(graph)), curr_edge_index_(index) {

      // decide if this is a valid edge index, to avoid crossing edges
      assert( index <= graph->num_edges() );
    }

    Graph* graph_;
    size_type curr_edge_index_;
    bool nodes_flipped_;

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

  



 private:

   // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

// 
  struct edge_type_ {
    size_type index_1; // index of node 1
    size_type index_2;  //  index of node 2
    edge_value_type value;  //  edge length

    bool operator<(edge_type_ e2) const {
      return std::tie(index_1, index_2) < std::tie(e2.index_1, e2.index_2);
    }
  };



  struct node_type_ { // 
    Point point;
    node_value_type value;
    std::vector<size_type> adj_list; // adjacent nodes of current node
  };

  // Use vector to store edges and nodes
  std::vector<node_type_> nodes_;
  std::vector<edge_type_> edges_;

  // cannot use copy of graph for constructor
  Graph(const Graph&) = delete;

  // cannot use & operator for graph
  Graph& operator=(const Graph&) = delete;
};

#endif // CME212_GRAPH_HPP
