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
template <typename V = int>
class Graph
{

  struct internal_node;
  struct internal_edge;

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
  Graph()
  { 
    // Initialize empty data structures
    nodes_ = std::vector<internal_node>(0);
    edges_ = std::vector<internal_edge>(0);
    adj_list_ = std::vector<std::vector<size_type>>(0);
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
  class Node : private totally_ordered<Node>
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
    Node()
    {
    }

    /** Return this node's position. */
    const Point &position() const
    {
      return (graph_->nodes_)[idx_].pt;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {
      return idx_;
    }

    /** Return node value, can be modified */
    node_value_type &value(){
      return (graph_->nodes_)[idx_].node_val;
    };

    /** Return node value, access only */
    const node_value_type &value() const{
      return (graph_->nodes_)[idx_].node_val;
    };

    /** Return the number of the neighbors the current node has */
    size_type degree() const
    {
      return (graph_->adj_list_)[idx_].size();
    };

    /** Returns an incident iterator pointing at the start of the incident */
    incident_iterator edge_begin() const
    {
      return incident_iterator(graph_, idx_, 0);
    };

    /** Returns an incident iterator pointing at the end of the incident */
    incident_iterator edge_end() const
    {
      return incident_iterator(graph_, idx_, degree());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const
    {
      return graph_ == n.graph_ && index() == n.index();
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node &n) const
    {
      // Order by comparing node index
      if (graph_ < n.graph_){
        return true;
      }

      return graph_ == n.graph_ && index() < n.index();
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph container
    Graph *graph_;
    // Unique index
    size_type idx_;

    /** Private Constructor */
    Node(const Graph *graph, size_type idx) : graph_(const_cast<Graph *>(graph)), idx_(idx){}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_val The new node's node value. Defult is set.
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_val = node_value_type()){
    // Update adjacency map
    adj_list_.push_back(std::vector<size_type>(0)); // create an empty vector to store the neighbors

    // Add new node to the node vector
    internal_node new_node = internal_node{position, node_val};
    nodes_.push_back(new_node);
    return Node(this, num_nodes() - 1);
  };

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const
  {
    return n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
    return Node(this, i);
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
  class Edge : private totally_ordered<Edge>
  {
  public:
    /** Construct an invalid Edge. */
    Edge()
    {
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      return Node(graph_, node_idx_1_);
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      return Node(graph_, node_idx_2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const
    {
      return (node_idx_1_ == e.node_idx_1_ && node_idx_2_ == e.node_idx_2_) || (node_idx_1_ == e.node_idx_2_ && node_idx_2_ == e.node_idx_1_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const
    { 

      // returns the order of graph if the edges are from different graphs
      if (graph_ < e.graph_){
        return true;
      }

      if (graph_ > e.graph_){
        return false;
      }

      // Otherwise compare the min and max node indices
      size_type min_idx = std::min(node_idx_1_, node_idx_2_);
      size_type e_min_idx = std::min(e.node_idx_1_, e.node_idx_2_);
      size_type max_idx = std::max(node_idx_1_, node_idx_2_);
      size_type e_max_idx = std::max(e.node_idx_1_, e.node_idx_2_);

      return min_idx < e_min_idx || (min_idx == e_min_idx && max_idx < e_max_idx);
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Pointer back to the Graph container
    Graph *graph_;
    size_type node_idx_1_; // indices for the nodes
    size_type node_idx_2_;

    /** Private Consstructor */
    Edge(const Graph *graph, size_type idx_1, size_type idx_2) : graph_(const_cast<Graph *>(graph)), node_idx_1_(idx_1), node_idx_2_(idx_2){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    assert(i < num_edges());
    return Edge(this,edges_[i].n1_idx,edges_[i].n2_idx);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    // Look up in the adjacency map
    assert(a.index()<size());
    assert(b.index()<size());
    std::vector<size_type> neighbors = adj_list_.at(a.index());
    return std::find(neighbors.begin(), neighbors.end(),b.index()) != neighbors.end();

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
  Edge add_edge(const Node &a, const Node &b)
  {
    // Create new edge if it doesn't exist
    if (!has_edge(a, b))
    {
      edges_.push_back(internal_edge{a.index(), b.index()});                       

      // Add to adjacency list (both directions)
      adj_list_[a.index()].push_back(b.index());
      adj_list_[b.index()].push_back(a.index());
    }

    return Edge(this, a.index(),b.index());
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    // Clear data structures for nodes and edges
    nodes_.clear();
    edges_.clear();
    adj_list_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
    {
    }

    //Dereference operator
    Node operator*() const
    {
      return Node(graph_, node_idx_);
    }

    // Increment
    NodeIterator &operator++()
    {
      ++node_idx_;
      return *this;
    }

    // Equality operator 
    bool operator==(const NodeIterator &other) const
    {
      return graph_ == other.graph_ && node_idx_ == other.node_idx_;
    }

  private:
    friend class Graph;
    Graph *graph_; // pointer to graph
    size_type node_idx_; // current node index

    //Private constructor that can be accessed by the Graph class
    NodeIterator(const Graph *graph, size_t node_idx) : graph_(const_cast<Graph *>(graph)), node_idx_(node_idx){}
  };

  // Return an iterator pointing at the start of the node (in the entire graph)
  node_iterator node_begin() const
  {
    return NodeIterator(this, 0);
  }

  // Return an iterator pointing at the end of the node
  node_iterator node_end() const
  {
    return NodeIterator(this, num_nodes());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator()
    {
    }

    // Dereference
    Edge operator*() const
    {
      // Get the index of the current neighbor node
      size_type neighbor_node_idx = (graph_->adj_list_).at(node_idx_).at(neighbor_idx_);

      return Edge(graph_, node_idx_,neighbor_node_idx);
    }

    // Increment
    IncidentIterator &operator++()
    {
      ++neighbor_idx_;
      return *this;
    };

    // Equality
    bool operator==(const IncidentIterator &other) const
    {
      return graph_ == other.graph_ && node_idx_ == other.node_idx_ && neighbor_idx_ == other.neighbor_idx_;
    };

  private:
    friend class Graph;
    Graph *graph_; // graph pointer
    size_type node_idx_; // current node to iterate on 
    size_type neighbor_idx_; // current neighbor to the node

    // Private constructor
    IncidentIterator(const Graph *graph, size_t node_idx, size_t neighbor_idx) : graph_(const_cast<Graph *>(graph)), node_idx_(node_idx), neighbor_idx_(neighbor_idx){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator()
    {
    }

    // Dereference
    Edge operator*() const
    { 
      size_type node_1 = graph_->edges_[edge_idx_].n1_idx;
      size_type node_2 = graph_->edges_[edge_idx_].n2_idx;
      return Edge(graph_, node_1, node_2);
    }

     // Increment
    EdgeIterator &operator++()
    {
      ++edge_idx_;
      return *this;
    }

    // Equality
    bool operator==(const EdgeIterator &other) const
    {
      return graph_ == other.graph_ && edge_idx_ == other.edge_idx_;
    }

  private:
    friend class Graph;
    Graph *graph_; // pointer to graph
    size_type edge_idx_; // edge index (which we keep track of in edges_)

    // Private constructor
    EdgeIterator(const Graph *graph, size_t edge_idx) : graph_(const_cast<Graph *>(graph)), edge_idx_(edge_idx){}
  };

  // Pointer to the starting edge
  edge_iterator edge_begin() const{
    return edge_iterator(this,0);
  }

  // Pointer to ending edge
  edge_iterator edge_end() const{
    return edge_iterator(this,num_edges());
  }

private:
  // Internal type for nodes (stores point info)
  struct internal_node
  {
    Point pt; // point object of the node
    node_value_type node_val;
  };

  // Internal type for edges (stores node indices)
  struct internal_edge
  {
    size_type n1_idx; // index of node 1
    size_type n2_idx; // index of node 2
  };

  // Use a vector to store internal nodes
  std::vector<internal_node> nodes_;
  std::vector<std::vector<size_type>> adj_list_; // Adjacency lists to keep track of nodes and their neighbors
  std::vector<internal_edge> edges_;            // Maps edge index to node indices

};

#endif // CME212_GRAPH_HPP
