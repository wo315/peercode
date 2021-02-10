

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
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** Predeclaration of Node type. */
  class Node;
    node_value_type& value();
    const node_value_type& value () const;

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
    : nodes_(), size_(0), edge_node1_(), edge_node2_(), size_edge_(0), 
      adj_list_(), adj_list_index_(), nodes_values_() {
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g_->nodes_[uid_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    //* Store value of node by indexing. */
    node_value_type& value() {
       return g_->nodes_values_[uid_];
    }

    //* Store value of node as const by indexing. */ 
    const node_value_type& value() const {
      return g_->nodes_values_[uid_];
    }

    //* Store the degree of a node, i.e., number of incident edges. */
    size_type degree() const {
      return g_->adj_list_[uid_].size();
    }

    //* Iterator to start at edge with index 0 for vertices incident to node. */
    IncidentIterator edge_begin() const {
      return IncidentIterator(this->g_, uid_, g_->adj_list_index_[uid_], g_->adj_list_[uid_], 0);
    }

    //* Iterator to end at edge with index one beyond the last edge for vertices incident to node. */
    IncidentIterator edge_end() const {
      return IncidentIterator(this->g_, uid_, g_->adj_list_index_[uid_], g_->adj_list_[uid_], degree());
    }    

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->g_ == n.g_ && this->index() == n.index()){
        return true;
      }
      else{
        return false;
      }
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
      if (this->g_ == n.g_ && this->index() < n.index()) {
        return true;
      }
      else{
        return false;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    //Pointer back to graph container
    Graph* g_;
    //This Node's unique id
    size_type uid_;
    //Private Constructor
    Node(const Graph* g, size_type uid)
      : g_(const_cast<Graph*>(g)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_;
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
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    // HW0: YOUR CODE HERE
    //Append to set of nodes in Graph
    nodes_.emplace_back(position);
    //Update size and index
    ++size_;
    //Update adjacency list with node andn index of node
    adj_list_.push_back(std::vector<size_type> ());
    adj_list_index_.push_back(std::vector<size_type> ());
    //Update nodes values
    nodes_values_.push_back(node_value_type());
    return Node(this, size_ - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this == n.g_) {
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
    // HW0: YOUR CODE HERE
    assert(i<num_nodes());
    return Node(this, i);    // Quiet compiler warning
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(g__, uid_1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(g__, uid_2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (this->g__ == e.g__ && this->uid_edge_ == e.uid_edge_) {
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
      //HW0: YOUR CODE HERE
      if (this->g__ == e.g__ && this->uid_edge_ < e.uid_edge_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* g__;
    //Pointer to nodes 1 and 2
    size_type uid_1;
    size_type uid_2;
    //Pointer to edge
    size_type uid_edge_;
    //Private Constructor
    Edge(const Graph* g, size_type uid1, size_type uid2, size_type uid_edge)
      : g__(const_cast<Graph*>(g)), uid_1(uid1), uid_2(uid2), uid_edge_(uid_edge) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return size_edge_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<num_edges());
    return Edge(this, edge_node1_[i], edge_node2_[i], i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (has_node(a) && has_node(b)) {
      std::vector<size_type> small_vec = adj_list_[a.index()];
      for (size_type j = 0; j<small_vec.size(); ++j) {
        if (small_vec[j] == b.index()) {
          return true;
        }
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
    // HW0: YOUR CODE HERE
    if (has_node(a) && has_node(b) && !(a==b)) {
      if (!(has_edge(a, b)) && !(has_edge(b, a))) {
        edge_node1_.push_back(a.index());
        edge_node2_.push_back(b.index());
        adj_list_[a.index()].push_back(b.index());
        adj_list_[b.index()].push_back(a.index());
        adj_list_index_[a.index()].push_back(size_edge_);
        adj_list_index_[b.index()].push_back(size_edge_);
        ++size_edge_;
        return Edge(this, a.index(), b.index(), size_edge_-1);
      }
      else {
        std::vector<size_type> small_vec = adj_list_[a.index()];
        for (size_type k = 0; k<small_vec.size(); ++k) {
          if (small_vec[k] == b.index()) {
            return Edge(this, a.index(), b.index(), k);
          }
        }
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
    std::vector<Point> nodes_;
    size_ = 0;
    std::vector<double> edge_node1_; 
    std::vector<double> edge_node2_; 
    size_edge_ = 0;
    std::vector<std::vector<size_type>> adj_list_;
    std::vector<std::vector<size_type>> adj_list_index_;
    std::vector<node_value_type> nodes_values_;
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

    // HW1 #2: YOUR CODE HERE

   /** Increment a node by one. */
   NodeIterator& operator++() {
     ++node_counter;
     return *this;
   }

   /** Check if two nodes are equal by checking their index and graph. */
   bool operator==(const NodeIterator& node_iter) const {
     return g_node_iter == node_iter.g_node_iter && node_iter.node_counter == node_counter;
   }

   /** Dereference the pointer to a node. */
   Node operator*() const {
     return Node(g_node_iter, node_counter);
   }

   private:
    friend class Graph;
    Graph* g_node_iter;
    size_type node_counter;
    NodeIterator(const Graph* g_node_iter_, size_type node_counter_) {
      node_counter = node_counter_;
      g_node_iter = const_cast<Graph*>(g_node_iter_);
    }
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE

  /** Iterator to start at the node with index 0. */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Iterator to end at node with index one beyond the last node. */
  NodeIterator node_end() const {
    return NodeIterator(this, size());
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

   /** Dereference the pointer to an edge. */
   Edge operator*() const {
     return Edge(g_inc_iter, node_id, connect_nodes[edge_node_ptr], edge_indices[edge_node_ptr]);
   }

   /** Increment an incident edge by one. */
   IncidentIterator& operator++() {
      ++edge_node_ptr;
      return *this;
   }

   /** Check if two edges are equal by checking if the pointer, graphs and nodes are the same. */
   bool operator==(const IncidentIterator& inc_iter) const {
      return inc_iter.edge_node_ptr == edge_node_ptr && g_inc_iter == inc_iter.g_inc_iter && node_id == inc_iter.node_id;
   }

   private:
    friend class Graph;
    Graph* g_inc_iter;
    size_type edge_node_ptr;
    std::vector<size_type> edge_indices;
    std::vector<size_type> connect_nodes;
    size_type node_id;
    IncidentIterator(const Graph* g_inc_iter_, size_type node_id_, std::vector<size_type> edge_indices_,
                     std::vector<size_type> connect_nodes_, size_type edge_node_ptr_) {
      edge_node_ptr = edge_node_ptr_;
      g_inc_iter = const_cast<Graph*>(g_inc_iter_);
      edge_indices = edge_indices_;
      connect_nodes = connect_nodes_;
      node_id = node_id_;
    }
    // HW1 #3: YOUR CODE HERE
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

   /** Dereference the pointer to an edge. */
   Edge operator*() const {
     return g_edge_iter->edge(edge_counter);
   }

   /** Increment an edge by one. */
   EdgeIterator& operator++() {
     ++edge_counter;
     return *this;
   }

   /** Check if two edges are equal by comparing the graph and edge index. */
   bool operator==(const EdgeIterator& edge_val) const {
     return edge_counter == edge_val.edge_counter && g_edge_iter == edge_val.g_edge_iter;
   }

   private:
    friend class Graph;
    Graph* g_edge_iter;
    int edge_counter;
    EdgeIterator(const Graph* g_edge_iter_, size_type edge_counter_) {
      edge_counter = edge_counter_;
      g_edge_iter = const_cast<Graph*>(g_edge_iter_);
    }
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE

  /** Iterator to start at the edge with index 0. */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Iterator to end at edge with index one beyond the last edge. */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  std::vector<Point> nodes_; //Store points in graph
  size_type size_; //Store number of nodes
  std::vector<size_type> edge_node1_; //Store first node in edge
  std::vector<size_type> edge_node2_; //Store second node in edge
  size_type size_edge_; //Store number of edges
  std::vector<std::vector<size_type>> adj_list_;
  std::vector<std::vector<size_type>> adj_list_index_;
  std::vector<node_value_type> nodes_values_;
};

#endif // CME212_GRAPH_HPP
