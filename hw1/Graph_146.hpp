#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

using namespace std;
/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
  
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;
  
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
  Graph(){
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
  class Node:private totally_ordered<Node>
  {
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
	  return graph_->points[node_idx];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
	  return node_idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
	
    /** Return this node's value. */
    node_value_type& value(){
      return graph_->val[node_idx];
	}
	
    /** Return this node's value. */
    const node_value_type& value() const{
      return graph_->val[node_idx];
	}
	
    /** Return the number of incident edges to this node. */
    size_type degree() const{
      return graph_->node2node[node_idx].size();
	}
	
    /** Return an iterator pointing to the first incident edge. */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_,node_idx,graph_->node2node[node_idx].begin());
	}
    
    /** Return an iterator pointing to the last incident edge. */
	incident_iterator edge_end() const{
      return IncidentIterator(graph_,node_idx,graph_->node2node[node_idx].end());
	}

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
	  return (graph_==n.graph_ and node_idx==n.node_idx);
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
	  return (node_idx<n.node_idx);
    }
	
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
	Graph* graph_;
	size_type node_idx;
	Node(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), node_idx(idx){
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const { 
    return points.size();
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
  Node add_node(const Point& position, const node_value_type&  v=node_value_type ()) {
    points.push_back(position);
    val.push_back(v);
    Node n=Node(this,size()-1);
    return n; 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n.graph_==this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this,i);
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
  class Edge:private totally_ordered<Edge>
  {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_,n1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,n2); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
	  return (graph_==e.graph_ and 
	         (node1()==e.node1() and node2()==e.node2()) or 
             (node1()==e.node2() and node2()==e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
	  return (edge_idx<e.edge_idx);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

	Graph* graph_;
	size_type edge_idx;
	size_type n1;
	size_type n2;
	Edge(const Graph* graph, size_type idx, size_type n1idx, size_type n2idx)
        : graph_(const_cast<Graph*>(graph)), edge_idx(idx), n1(n1idx), n2(n2idx){
    }	
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {	  
    return edge2node.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this,i,edge2node.at(i)[0],edge2node.at(i)[1]); 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    auto s1=node2node.find(a.index());
    if (s1!=node2node.end()){
      auto s2=s1->second.find(b.index());
      if(s2!=s1->second.end()){
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
    if (has_edge(a, b)){
      size_type idx=node2node[a.index()][b.index()];
      return Edge(this, idx, a.index(), b.index());
    }

    size_type i=edge2node.size();
    vector<size_type> idx_nodes{a.index(),b.index()};
    edge2node[i]=idx_nodes;
    node2node[a.index()][b.index()]=i;
    node2node[b.index()][a.index()]=i;
    return Edge(this,i,a.index(), b.index()); 
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    points.clear();
    val.clear();
    edge2node.clear();
    node2node.clear();
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
	
    /** Dereference the NodeIterator. */
    Node operator*() const{
      return Node(graph_,node_idx);
    }
    
    /** Pre-increment operator. */
    NodeIterator& operator++(){
      node_idx++;
      return *this;
    }
	
    /** Test NodeIterator equality. */
    bool operator==(const NodeIterator& NodeIter) const{
      return (graph_==NodeIter.graph_ and node_idx==NodeIter.node_idx);
    }
	
    /** Test NodeIterator inequality. */   
    bool operator!=(const NodeIterator& NodeIter) const{
      return !(NodeIter == *this);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_idx;

    NodeIterator(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), node_idx(idx) {
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return an iterator pointing to the first node. */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  /** Return an iterator pointing to the last node. */
  node_iterator node_end() const{
    return NodeIterator(this, size());
  }

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
	
    /** Dereference the IncidentIterator. */
    Edge operator*() const{
      size_type n1=graph_->edge2node[it->second][0];
      size_type n2=graph_->edge2node[it->second][1];
      if (node_idx==n1){
        return Edge(graph_,it->second,n1,n2);
      }
      else{
        return Edge(graph_,it->second,n2,n1);
      }
    }
    
    /** Pre-increment operator. */	
    IncidentIterator& operator++(){		
      it++;
      return *this;
    }

    /** Test IncidentIterator equality. */		 
    bool operator==(const IncidentIterator& IncidentIter) const{
      return (graph_==IncidentIter.graph_ and
              node_idx==IncidentIter.node_idx and 
              it==IncidentIter.it);
    }
	
    /** Test IncidentIterator inequality. */ 
	bool operator!=(const IncidentIterator& IncidentIter) const{
      return !(IncidentIter == *this);
    }

   private:
    friend class Graph;
	Graph* graph_;
    size_type node_idx;
    std::unordered_map<size_type,size_type>::iterator it;
	
    IncidentIterator(const Graph* graph, size_type idx, std::unordered_map<size_type,size_type>::iterator iter)
        : graph_(const_cast<Graph*>(graph)), node_idx(idx), it(iter) {
    }
	
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
    
    /** Dereference the EdgeIterator. */
    Edge operator*() const{
      return Edge(graph_,edge_idx, graph_->edge2node[edge_idx][0],graph_->edge2node[edge_idx][1]);
	}
	
    /**Pre-increment operator.*/
    EdgeIterator& operator++(){
      edge_idx++;
      return *this;
	}
	
    /**Test EdgeIterator equality.*/
    bool operator==(const EdgeIterator& EdgeIter) const{
      return (graph_==EdgeIter.graph_ and edge_idx==EdgeIter.edge_idx);
    }

    /**Test EdgeIterator inequality.*/
    bool operator!=(const EdgeIterator& EdgeIter) const{
      return !(EdgeIter == *this);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
	Graph* graph_;
	size_type edge_idx;
	
	EdgeIterator(const Graph* graph, size_type idx)
        : graph_(const_cast<Graph*>(graph)), edge_idx(idx) {
    }
	
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return an iterator pointing to the first edge. */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  /** Return an iterator pointing to the last edge. */
  edge_iterator edge_end() const{
    return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  // helper functions, data members, and so forth.
  vector<Point> points;
  vector<node_value_type> val;
  unordered_map<size_type, vector<size_type>> edge2node;
  unordered_map<size_type, unordered_map<size_type,size_type>> node2node;
   
};

#endif // CME212_GRAPH_HPP
